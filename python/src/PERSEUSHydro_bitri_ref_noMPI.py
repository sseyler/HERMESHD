#      DG-Hydro in 3 dimensions: Extension of the finite volume PERSEUS Algorithm (Physics of the Extended-mhd Relaxation System using an Efficient Upwind Scheme,
#                                by Matt Martin and Charles Seyler) to the Discontinuous Galerkin Method specialized to hydrodynamics.
#
#                DG Algorithm by Xuan Zhao, Nat Hamlin and Charles Seyler.
#
#                Solves the 3D compressible Euler equations
#
#                Variables:    There are 5 dependent variables: density (rh), velocity (vx,vy,vz), energy (en)
#
#                Units:        A value of unity for each variable or parameter corresponds to the following dimensional units:
#
#                              Length              L0
#                              Time                t0
#                              Number density      n0
#                              Velocity            v0
#                              Temperature         te0
#
#-------------------------------------------------------------------------------------------------------------------------------------------

import numpy as np

# Indices for dynamical variables
rh  = 0
mx, my, mz  = 1, 2, 3
en  = 4
pxx, pyy, pzz = 5, 6, 7
pxy, pxz, pyz = 8, 9, 10
nQ  = 11  # total number of dynamical variables

################################################################################
# Grid dimensions and DG basis points
#------------------------------------
# The jump in accuracy between the linear basis (nbasis = 4) and quadratic basis (nbasis = 10)
# is much greater than the jump in accuracy between quadratic and cubic (nbasis = 20) bases.
################################################################################
nx, ny, nz = 20, 20, 1

ngu     = 0
nbasis  = 8
nbastot = 27
# nbasis = 4: {1,x,y,z}
# nbasis = 10: {1,x,y,z, P_2(x),P_2(y),P_2(z), yz, zx, xy}
# nbasis = 20: nbasis10 + {xyz,xP2(y),yP2(x),xP2(z),zP2(x),yP2(z),zP2(y),P3(x),P3(y),P3(z)}

# iquad: number of Gaussian quadrature points per direction.
# should not be less than ipoly, where ipoly is the maximum Legendre polynomial
# order used, otherwise unstable. iquad should not be larger than ipoly + 1,
# which is an exact Gaussian quadrature for the Legendre polynomial used. Thus
# there are only 2 cases of iquad per given nbasis.
# Both cases give similar results, but iquad = ipoly + 1 case is formally more
# accurate.
iquad = 2
nedge = iquad

# nface: number of quadrature points per cell face.
# npg: number of internal points per cell.
nface = iquad*iquad
npg   = nface*iquad
nfe   = 2*nface
npge  = 6*nface
nslim = npg + 6*nface

# Boundary condition parameters:
#   if "2" then periodic. MPI does this for you. If "0" then set_bc subroutine
#   used to set BCs
xlbc = 2
xhbc = 2
ylbc = 2
yhbc = 2
zlbc = 2
zhbc = 2

# Set output frequency and temporal integration order
ntout  = 100
iorder = 2

# Choose Riemann solver for computing fluxes.
#   Select a solver by setting its value to 1. If all solvers are set to 0, then
#   LLF is used for fluxes.
#
#   Note: LLF is very diffusive for the hydro problem. Roe and HLLC are much
#   less diffusive than LLF and give very similar results with similar cpu
#   overhead. Only HLLC is set up to handle water EOS (ieos = 2).
ihllc = 1
iroe  = 0
ieos  = 1


# to restart from a checkpoint, set iread to 1 or 2 (when using odd/even scheme)
iread  = 0
iwrite = 0
fpre   = 'Qout'
resuming = False

lx = 100.0
ly = 100.0
lz = 100.0/120.0

# Final time
tf = 1000.0

# Predefined parameters
pi = 4.0*atan(1.0)
aindex = 5.0/3.0  # adiabatic index
mu = 18.0         # atomic mass (water)

aindm1 = aindex - 1.0
cp     = aindex/(aindex - 1.0)
clt    = 2.0
vis    = 1.0e-1
epsi   = 5.0
amplv  = 0.0
amplm  = 0.0
amplen = 10.0

# dimensional (base) units (expressed in MKS)
L0 = 1.0e-9
t0 = 1.0e-12
n0 = 3.32e28

# derived units
v0     = L0/t0
p0     = mu*1.67e-27*n0*v0**2
te0    = p0/n0/1.6e-19
P_1    = 2.15e9/7.2/p0
P_base = 1.01e5/p0  # atmospheric pressure

# rh_min is a minimum density to be used for ideal gas law EOS and rh_min is the minimum density below which the
# pressure becomes negative for the water EOS: P = P_1*(density**7.2 - 1.) + P_base
# The DG based subroutine "limiter" prevents the density from falling below rh_mult*rh_min.
# Note: the EOS for water is likely to be a critical player in getting the fluctuating hydrodynamics correct.
# The EOS used here is a simple one called Tait-Murnaghan.  There are other much more sophisicated EOS's
# some of which take into account ionic solutions. It would be worthwhile to investigate this further and
# experiment with different EOS's.
rh_floor = 1.0e-1
rh_mult  = 1.01
rh_min   = rh_mult*(1.0 - P_base/P_1)**(1.0/7.2)

T_floor = 0.026/te0
P_floor = T_floor*rh_floor

# real, dimension(nx,ny,nz,nQ,nbasis) ::  Q_r0, Q_r1, Q_r2, Q_r3
# real, dimension(nx,ny,nz,nQ,nbasis) ::  glflux_r, source_r
# real, dimension(nx,ny,nz,nQ,nbasis) :: integral_r
Q_r0 = np.zeros((nx,ny,nz,nQ,nbasis))
Q_r1 = np.zeros((nx,ny,nz,nQ,nbasis))
Q_r2 = np.zeros((nx,ny,nz,nQ,nbasis))
Q_r3 = np.zeros((nx,ny,nz,nQ,nbasis))
glflux_r   = np.zeros((nx,ny,nz,nQ,nbasis))
source_r   = np.zeros((nx,ny,nz,nQ,nbasis))
integral_r = np.zeros((nx,ny,nz,nQ,nbasis))

# real eta(nx,ny,nz,npg),den0(nx,ny,nz),Ez0,Zdy(nx,ny,nz,npg)
# real flux_x(nface,1:nx+1,ny,nz,1:nQ), flux_y(nface,nx,1:ny+1,nz,1:nQ), flux_z(nface,nx,ny,1:nz+1,1:nQ)
# real cfrx(nface,nQ),cfry(nface,nQ),cfrz(nface,nQ)
Ez0 = 0.0
den0 = np.zeros((nx,ny,nz))
eta = np.zeros((nx,ny,nz,npg))
Zdy = np.zeros((nx,ny,nz,npg))

flux_x = np.zeros((nface, 1:nx+1, ny,     nz,     1:nQ))
flux_y = np.zeros((nface, nx,     1:ny+1, nz,     1:nQ))
flux_z = np.zeros((nface, nx,     ny,     1:nz+1, 1:nQ))

cfrx = np.zeros((nface,nQ))
cfry = np.zeros((nface,nQ))
cfrz = np.zeros((nface,nQ))

# Deal with the logical thing (try numpy mask arrays)
logical MMask(nx,ny,nz),BMask(nx,ny,nz)
xcell = np.zeros((npg))
ycell = np.zeros((npg))
zcell = np.zeros((npg))
xface = np.zeros((npge))
yface = np.zeros((npge))
zface = np.zeros((npge))
# integer ticks, count_rate, count_max
# real t1, t2, t3, t4, elapsed_time, t_start, t_stop, dtoriginal
# real t,dt,tout,dtout,vf,dxi,dyi,dzi,loc_lxd,loc_lyd,loc_lzd,check_Iz,sl, dz, dy, dx
# real lxd,lxu,lyd,lyu,lzd,lzu,pin_rad,pin_height,foil_thickness,rh_foil,rh_fluid,pin_rad_in,pin_rad_out,rim_rad
# real disk_rad,disk_height,foil_rad,buf_rad,buf_z,dish_height,foil_height
# real gpz_rad,rh_gpz,kappa

Qxhigh_ext = np.zeros((ny,nz,nface,nQ))
Qxlow_int  = np.zeros((ny,nz,nface,nQ))
Qxlow_ext  = np.zeros((ny,nz,nface,nQ))
Qxhigh_int = np.zeros((ny,nz,nface,nQ))

Qyhigh_ext = np.zeros((nx,nz,nface,nQ))
Qylow_int  = np.zeros((nx,nz,nface,nQ))
Qylow_ext  = np.zeros((nx,nz,nface,nQ))
Qyhigh_int = np.zeros((nx,nz,nface,nQ))

Qzhigh_ext = np.zeros((nx,ny,nface,nQ))
Qzlow_int  = np.zeros((nx,ny,nface,nQ))
Qzlow_ext  = np.zeros((nx,ny,nface,nQ))
Qzhigh_int = np.zeros((nx,ny,nface,nQ))

mxa = np.zeros((3))
mya = np.zeros((3))
mza = np.zeros((3))
kroe = np.zeros((nface))
# niter,iseed


# Parameters relating to quadratures and basis functions.
#   wgt1d: quadrature weights for 1-D integration
#   wgt2d: quadrature weights for 2-D integration
#   wgt3d: quadrature weights for 3-D integration
wgt1d = np.zeros((5))
wgt2d = np.zeros((30))
wgt3d = np.zeros((100))
cbasis = np.zeros((nbastot))

# real, dimension(nface,nbastot) :: bfvals_zp, bfvals_zm, bfvals_yp, bfvals_ym, bfvals_xp, bfvals_xm
# real bf_faces(nslim,nbastot), bfvals_int(npg,nbastot),xquad(20)
bfvals_zp = np.zeros((nface,nbastot))
bfvals_zm = np.zeros((nface,nbastot))
bfvals_yp = np.zeros((nface,nbastot))
bfvals_ym = np.zeros((nface,nbastot))
bfvals_xp = np.zeros((nface,nbastot))
bfvals_xm = np.zeros((nface,nbastot))

bf_faces   = np.zeros((nslim,nbastot))
bfvals_int = np.zeros((npg,nbastot))
xquad = np.zeros((20))


# real bval_int_wgt(npg,nbastot),wgtbfvals_xp(nface,nbastot),wgtbfvals_yp(nface,nbastot),wgtbfvals_zp(nface,nbastot)
# real wgtbfvals_xm(nface,nbastot),wgtbfvals_ym(nface,nbastot),wgtbfvals_zm(nface,nbastot)
# real wgtbf_xmp(4,2,nbastot),wgtbf_ymp(4,2,nbastot),wgtbf_zmp(4,2,nbastot)
# real sumx,sumy,sumz
# integer i2f,i01,i2fa
bval_int_wgt = np.zeros((npg,nbastot))
wgtbfvals_xp = np.zeros((nface,nbastot))
wgtbfvals_yp = np.zeros((nface,nbastot))
wgtbfvals_zp = np.zeros((nface,nbastot))
wgtbfvals_xm = np.zeros((nface,nbastot))
wgtbfvals_ym = np.zeros((nface,nbastot))
wgtbfvals_zm = np.zeros((nface,nbastot))
wgtbf_xmp = np.zeros((4,2,nbastot))
wgtbf_ymp = np.zeros((4,2,nbastot))
wgtbf_zmp = np.zeros((4,2,nbastot))

# integer,parameter :: kx=2,ky=3,kz=4,kyz=5,kzx=6,kxy=7,kxyz=8,kxx=9,kyy=10,kzz=11,kyzz=12,kzxx=13,kxyy=14
# integer,parameter :: kyyz=15,kzzx=16,kxxy=17,kyyzz=18,kzzxx=19,kxxyy=20,kyzxx=21,kzxyy=22,kxyzz=23
# integer,parameter :: kxyyzz=24,kyzzxx=25,kzxxyy=26,kxxyyzz=27
kx   = 1
ky   = 2
kz   = 3
kyz  = 4
kzx  = 5
kxy  = 6
kxyz = 7
kxx  = 8
kyy  = 9
kzz  = 10

kyzz  = 11
kzxx  = 12
kxyy  = 13
kyyz  = 14
kzzx  = 15
kxxy  = 16

kyyzz = 17
kzzxx = 18
kxxyy = 19
kyzxx = 20
kzxyy = 21
kxyzz = 22

kxyyzz  = 23
kyzzxx  = 24
kzxxyy  = 25
kxxyyzz = 26


# Parameters for VTK output.
nvtk  = 2
nvtk2 = nvtk*nvtk
nvtk3 = nvtk*nvtk*nvtk
nnx   = nx*nvtk * np.ones((I4P))  # FIXME
nny   = ny*nvtk * np.ones((I4P))  # FIXME
nnz   = nz*nvtk * np.ones((I4P))  # FIXME
bfvtk    = np.zeros((nvtk3,nbastot))
bfvtk_dx = np.zeros((nvtk3,nbastot))
bfvtk_dy = np.zeros((nvtk3,nbastot))
bfvtk_dz = np.zeros((nvtk3,nbastot))
xgrid = np.zeros((20))
# dxvtk,dyvtk,dzvtk

################################################################################
# MPI definitions
#------------------------------------
# mpi_nx = 6
# mpi_ny = 6
# print_mpi = 0

# integer iam,ierr,mpi_nz,numprocs,reorder,cartcomm,mpi_P,mpi_Q,mpi_R
# integer dims(3),coords(3),periods(3),nbrs(6),reqs(4),stats(MPI_STATUS_SIZE,4)
# integer,parameter:: NORTH=1,SOUTH=2,EAST=3,WEST=4,UP=5,DOWN=6,MPI_TT=MPI_REAL4
# dims    = np.zeros((3))
# coords  = np.zeros((3))
# periods = np.zeros((3))
# nbrs    = np.zeros((6))
# reqs    = np.zeros((4))
# stats   = np.zeros((MPI_STATUS_SIZE,4))
#
# NORTH = 1
# SOUTH = 2
# EAST  = 3
# WEST  = 4
# UP    = 5
# DOWN  = 6
# MPI_TT = MPI_REAL4
################################################################################

# real cflm
if nbasis <=  8: cflm = 0.14
if nbasis == 27: cflm = 0.1
if nbasis == 64: cflm = 0.08

# Initialize grid sizes and local lengths
cbasis[0]         = 1.      # coefficient for basis function {1}
cbasis[kx:kz+1]   = 3.      # coefficients for basis functions {x,y,z}
cbasis[kyz:kxy+1] = 9.      # coefficients for basis functions {yz,zx,xy}
cbasis[kxyz]      = 27.     # coefficient for basis function {xyz}

cbasis[kxx:kzz+1]       = 5.    # coefficients for basis functions {P2(x),P2(y),P2(z)}
cbasis[kyzz:kxyy+1]     = 15.   # coefficients for basis functions {yP2(z),zP2(x),xP2(y)}
cbasis[kyyz:kxxy+1]     = 15.   # coefficients for basis functions {P2(y)z,P2(z)y,P2(z)x}
cbasis[kyyzz:kxxyy+1]   = 25.   # coefficients for basis functions {P2(y)P2(z),P2(z)P2(x),P2(x)P2(y)}
cbasis[kyzxx:kxyzz+1]   = 45.   # coefficients for basis functions {yzP_2(x),zxP_2(y),xyP_2(z)}
cbasis[kxyyzz:kzxxyy+1] = 75.   # coefficients for basis functions {xP2(y)P2(z),yP2(z)P2(x),zP2(x)P2(y)}
cbasis[kxxyyzz]         = 125.  # coefficients for basis functions {P2(x)P2(y)P2(z)}

################################################################################
# MPI stuff
#------------------------------------
# call MPI_Init ( ierr )
# call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
#
# mpi_nz = numprocs/(mpi_nx*mpi_ny)
#
# dims[0] = mpi_nx
# dims[1] = mpi_ny
# dims[2] = mpi_nz
################################################################################


periods[:] = 0
if xhbc == 2:
    periods[0] = 1
if yhbc == 2:
    periods[1] = 1
if zhbc == 2:
    periods[2] = 1
reorder = 1

################################################################################
# MPI stuff
#------------------------------------
# call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, reorder,cartcomm, ierr)
# call MPI_COMM_RANK (cartcomm, iam, ierr )
# call MPI_CART_COORDS(cartcomm, iam, 3, coords, ierr)
# mpi_P = coords[0] + 1
# mpi_Q = coords[1] + 1
# mpi_R = coords[2] + 1
# call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(WEST), nbrs(EAST), ierr)
# call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(SOUTH), nbrs(NORTH), ierr)
# call MPI_CART_SHIFT(cartcomm, 2, 1, nbrs(DOWN), nbrs(UP), ierr)
################################################################################

half_length = lx/2.0
lxd = -half_length
lxu = half_length
half_length = ly/2.0
lyd = -half_length
lyu = half_length
half_length = lz/2.0
lzd = -half_length
lzu = half_length

dxi = nx/(lxu-lxd)
dyi = ny/(lyu-lyd)
dzi = nz/(lzu-lzd)

dx = 1./dxi
dy = 1./dyi
dz = 1./dzi
################################################################################
# MPI stuff
#------------------------------------
# loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
# loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
# loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz

rh_fluid = 1.0

# indices used in poynting() for computing Poynting Maxwell flux

mxa[0:3] = np.array([mx, my, mz])
mya[0:3] = np.array([my, mz, mx])
mza[0:3] = np.array([mz, mx, my])

t          = 0.0
dt         = cflm*dx/clt
dtoriginal = dt
nout  = 0
niter = 0
dtout = tf/ntout

# Evaluate local cell values of basis functions on cell interior and faces.
# This is done for 1, 2, or 3 point Gaussian quadrature.
call set_bfvals_3D  # FIXME

for ir in xrange(nbasis):
    for ipg in xrange(npg):
        bval_int_wgt[ipg,ir] = wgt3d[ipg]*bfvals_int[ipg,ir]

for ir in xrange(nbasis):
   wgtbfvals_xp[0:nface,ir] = wgt2d[0:nface]*bfvals_xp[0:nface,ir]
   wgtbfvals_yp[0:nface,ir] = wgt2d[0:nface]*bfvals_yp[0:nface,ir]
   wgtbfvals_zp[0:nface,ir] = wgt2d[0:nface]*bfvals_zp[0:nface,ir]
   wgtbfvals_xm[0:nface,ir] = wgt2d[0:nface]*bfvals_xm[0:nface,ir]
   wgtbfvals_ym[0:nface,ir] = wgt2d[0:nface]*bfvals_ym[0:nface,ir]
   wgtbfvals_zm[0:nface,ir] = wgt2d[0:nface]*bfvals_zm[0:nface,ir]

for ir in xrange(nbasis):
   wgtbf_xmp[0:5,0,ir] = -0.25*cbasis[ir]*dxi*wgtbfvals_xm[0:nface,ir]
   wgtbf_ymp[0:5,0,ir] = -0.25*cbasis[ir]*dyi*wgtbfvals_ym[0:nface,ir]
   wgtbf_zmp[0:5,0,ir] = -0.25*cbasis[ir]*dzi*wgtbfvals_zm[0:nface,ir]
   wgtbf_xmp[0:5,1,ir] =  0.25*cbasis[ir]*dxi*wgtbfvals_xp[0:nface,ir]
   wgtbf_ymp[0:5,1,ir] =  0.25*cbasis[ir]*dyi*wgtbfvals_yp[0:nface,ir]
   wgtbf_zmp[0:5,1,ir] =  0.25*cbasis[ir]*dzi*wgtbfvals_zp[0:nface,ir]

call init_random_seed(123456789)  # FIXME
################################################################################
# MPI stuff
#------------------------------------
# iseed = 1317345*mpi_P + 5438432*mpi_Q + 38472613*mpi_R

# if iam == print_mpi:
#     print 'total dim= ', nx, ny, nz
#     print 'te0 is: ',    te0
#     print 'dx is: ',     ly/(ny)*L0
#     print 'iquad is: ',  iquad
#     print 'nbasis is: ', nbasis
################################################################################

call system_clock(ticks, count_rate, count_max)  # FIXME
t_start = ticks*1.0/count_rate

if iread == 0:
    initial_condition(Q_r0, MMask)
else:
    # This applies only if the initial data are being read from an input file.
    # - If resuming a run, keep the previous clock (i.e., t at nout) running.
    # - If not resuming a run, treat input as initial conditions at t=0, nout=0.
    ################################################################################
    # MPI stuff
    #------------------------------------
    # call readQ(fpre,iam,iread,Q_r0,t_p,dt_p,nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p)
    ################################################################################

    if resuming:
        t    = t_p
        dt   = dt_p
        nout = nout_p
    # Note, nout=1 corresponds to t=dt, but nout=2 corresponds to t~dtout, etc.
    if nout > 1:
        dtout_p = t_p / (nout_p - 1)
    else:  # Automatically pass consistency check
        dtout_p = dtout
    ################################################################################
    # MPI stuff
    #------------------------------------
    # if iam == print_mpi:
    #     print 'resuming = ', resuming
    #     print 't = ', t
    #     print 'dt = ', dt
    #     print 'nout = ', nout
    #     print 'dtout_p = ', dtout_p, ' dtout = ', dtout

    ################################################################################
    # MPI stuff
    #------------------------------------
    # Quit if dtout is incompatible with input t/(nout-1)
    if abs(dtout_p-dtout)/dt_p > 1.01:
        # if iam == print_mpi:
        #     print 'Bad restart, non-matching dtout'
        call exit(-1)  # FIXME
    # if (mpi_nx_p != mpi_nx) or (mpi_ny_p != mpi_ny) or (mpi_nz_p != mpi_nz):
    #     if iam == print_mpi:
    #         print 'Bad restart, non-matching mpi_nx, mpi_ny, or mpi_nz'
    #     call exit(-1)
    ################################################################################


call system_clock( ticks, count_rate, count_max )  # FIXME
t1 = 1.0*ticks / count_rate
call output_vtk(Q_r0,nout,iam)  # FIXME

while t < tf:

    # if (mod(niter,200) == 0 .and. iam == print_mpi) print *,'niter,t,dt = ',niter,t,dt,dtout*nout
    niter = niter + 1
    call get_min_dt(dt)

    if iorder == 2:
        call prep_advance(Q_r0)
        call advance_time_level_gl(Q_r0,Q_r1)

        call prep_advance(Q_r1)
        call advance_time_level_gl(Q_r1,Q_r2)

        Q_r0 = 0.5*(Q_r0 + Q_r2)

    # if iorder == 3:
    #     call prep_advance(Q_r0)
    #     call advance_time_level_gl(Q_r0,Q_r1)
    #
    #     call prep_advance(Q_r1)
    #     call advance_time_level_gl(Q_r1,Q_r2)
    #
    #     Q_r3 = 0.75*Q_r0 + 0.25*Q_r2
    #
    #     call prep_advance(Q_r3)
    #     call advance_time_level_gl(Q_r3,Q_r2)
    #
    #     Q_r0 = (Q_r0 + 2.*Q_r2)/3.

    t = t+dt

    if t > dtout*nout:
        nout = nout + 1
        if iam == print_mpi:
            call system_clock(ticks, count_rate, count_max)  # FIXME
            t2 = 1.0*ticks/count_rate
            print 'Iteration time', (t2-t1), 'seconds'
            t1 = t2
            print 'nout = ', nout
            print "        t= ",t*100.,"         dt= ",dt

        ################################################################################
        # MPI stuff
        #------------------------------------
        # call MPI_BARRIER(cartcomm,ierr)
        call output_vtk(Q_r0,nout,iam)  # FIXME

        # write checkpoint files; assign an odd/even id to ensure last two sets are kept
        if iwrite == 1:
            ioe = 2 - mod(nout,2)  # FIXME
            call writeQ(fpre,iam,ioe,Q_r0,t,dt,nout,mpi_nx,mpi_ny,mpi_nz)  # FIXME

        ################################################################################
        # MPI stuff
        #------------------------------------
        # call MPI_BARRIER(cartcomm,ierr)

        if iam == print_mpi:
            call system_clock(ticks, count_rate, count_max)  # FIXME
            t2 = ticks/count_rate
            print 'Output time', (t2-t1), 'seconds'
            t1 = t2

################################################################################
# MPI stuff
#------------------------------------
# call MPI_Finalize (ierr)

# contains

#------------------------------------------
# FIXME
def init_random_seed(iseed):
    # integer :: i, n, clock,iseed
    # integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    if iseed == 0:
        seed = clock*(iam+1) + 37 * (/ (i - 1, i = 1, n) /)
    else:
        seed = iseed*(iam+1)
    call random_seed(put = seed)

#   print *,seed(1)
    deallocate(seed)

#-------------------------------------------------------------------------------

def prep_advance(Q_ri, bfvals_xm, bfvals_xp,
                       bfvals_ym, bfvals_yp,
                       bfvals_zm, bfvals_zp):
    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    #
    # EXT LIBS:
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * bf_faces (limiter)
    #  * bfvals_xm, bfvals_xp (prepare_exchange)
    #  * bfvals_ym, bfvals_yp (prepare_exchange)
    #  * bfvals_zm, bfvals_zp (prepare_exchange)
    #  * Qxlow_ext, Qxlow_int, Qxhigh_ext, Qxhigh_int (set_bc)
    #  * Qylow_ext, Qylow_int, Qyhigh_ext, Qyhigh_int (set_bc)
    #  * Qzlow_ext, Qzlow_int, Qzhigh_ext, Qzhigh_int (set_bc)
    #  * nx, ny, nz
    #  * mx, my, mz
    #  * xlbc, ylbc, zlbc
    #  * t
    # LOCAL:
    #  * Q_ri (pass-by-ref w/ shape (nx,ny,nz,nQ,nbasis))
    ############################################################################
    # real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri

    if ieos == 1:
        call limiter(Q_ri, bf_faces)
    if ieos == 2:
        call limiter2(Q_ri, bf_faces)

    prepare_exchange(Q_ri, bfvals_xm, bfvals_xp,
                           bfvals_ym, bfvals_yp,
                           bfvals_zm, bfvals_zp)
    set_bc(Qxlow_ext, Qxlow_int, Qxhigh_ext, Qxhigh_int,
           Qylow_ext, Qylow_int, Qyhigh_ext, Qyhigh_int,
           Qzlow_ext, Qzlow_int, Qzhigh_ext, Qzhigh_int)
    call flux_cal(Q_ri, bfvals_xm, bfvals_xp,
                        bfvals_ym, bfvals_yp,
                        bfvals_zm, bfvals_zp))
    call innerintegral(Q_ri)
    call glflux
    call source_calc(Q_ri,t)

#-------------------------------------------------------------------------------

def get_min_dt(dt):
   # real dt,dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0,vex,vey,vez,vem,vem0,dni,dn,vx,vy,vz,Pr,sqdni,vacc,vacc0,cs
   # integer :: i,j,k,main_proc=0,mpi_size=1,loc_reqs(numprocs-1),loc_stats(MPI_STATUS_SIZE,numprocs-1)

   vmag = 0.

    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx):
                dn = Q_r0[i,j,k,rh,0]
                dni = 1./dn
                vx = Q_r0[i,j,k,mx,0]*dni
                vy = Q_r0[i,j,k,my,0]*dni
                vz = Q_r0[i,j,k,mz,0]*dni
                if ieos == 1:
                    cs = sqrt(aindex*(Q_r0[i,j,k,en,0]*dni - 0.5*(vx**2 + vy**2 + vz**2)))
                if ieos == 2:
                    cs = sqrt(7.2*P_1*dn**6.2 + T_floor)

                vmag0 = max(abs(vx)+cs, abs(vy)+cs, abs(vz)+cs)
                if vmag0 > vmag:
                    vmag = vmag0

   vmax = vmag*dxi
   dt_min = 1.0*cflm/vmax  # time step is determined by the maximum flow + sound speed in the domain

    call MPI_BARRIER(cartcomm,ierr)
    if iam == main_proc:
        do i=1,numprocs-1
            call MPI_IRecv(dt_val(i),mpi_size,MPI_TT,i,0,cartcomm,loc_reqs(i),ierr)
        enddo
        call MPI_WaitAll(numprocs-1,loc_reqs,loc_stats,ierr)
        do i=1,numprocs-1
            dt_min=min(dt_min,dt_val(i))
        enddo
        do i=1,numprocs-1
            call MPI_ISend(dt_min,mpi_size,MPI_TT,i,0,cartcomm,loc_reqs(i),ierr)
        enddo
        call MPI_WaitAll(numprocs-1,loc_reqs,loc_stats,ierr)
    else:
        call MPI_ISend(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
        call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
        call MPI_Wait(reqs(1),stats(:,1),ierr)

    dt = dt_min

#----------------------------------------------------------------------------------------------

def initial_condition(Q_r0, MMask):
    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    #
    # EXT LIBS:
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * rh, en
    #  * T_floor, aindex, rh_fluid
    # LOCAL:
    #  * wtev
    ############################################################################

    wtev = T_floor

    Q_r0[:,:,:,:,:] = 0.0
    Q_r0[:,:,:,rh,0] = rh_floor
    Q_r0[:,:,:,en,0] = T_floor*rh_floor/(aindex - 1.)
    MMask[:,:,:] = False

    fill_fluid(Q_r0)

#-------------------------------------------------------

    # subroutine fill_fluid

#-------------------------------------------------------

def fill_fluid(Q_r0):

    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    # * xc(), yc(), zc()
    #
    # EXT LIBS:
    # * np.cosh(), np.uniform()
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * nx, ny, nz
    #  * rh, mx, my, mz, en
    #  * T_floor, aindex, rh_fluid, lyu
    # LOCAL:
    #  * i, j, k (loop vars)
    #  * wtev, rnum, xcc, ycc, zcc
    ############################################################################

    ############################################################################
    # MPI stuff
    #------------------------------------
    # iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

    wtev = T_floor

    # test problem is an unstable flow jet in x with velocity perturbations in y
    for i in xrange(nx):
        for j in xrange(ny):
            for k in xrange(nz):

                rnum = np.uniform(low=-0.5, high=0.5)

                Q_r0[i,j,k,rh,0] = rh_floor
                Q_r0[i,j,k,en,0] = wtev*Q_r0[i,j,k,rh,0]/(aindex-1.0)

                xcc = xc(i)
                ycc = yc(j)
                zcc = zc(k)

                Q_r0[i,j,k,rh,0] = rh_fluid
                Q_r0[i,j,k,my,0] = 0.001*rnum/np.cosh(20*ycc/lyu)**2
                Q_r0[i,j,k,mz,0] = 0.0
                Q_r0[i,j,k,mx,0] = 1.0*Q_r0[i,j,k,rh,0]/np.cosh(20*ycc/lyu)/1.0
                Q_r0[i,j,k,en,0] = wtev*Q_r0[i,j,k,rh,0]/(aindex-1.0)           \
                                 + 0.5*(Q_r0[i,j,k,mx,0]**2                     \
                                      + Q_r0[i,j,k,my,0]**2                     \
                                      + Q_r0[i,j,k,mz,0]**2)/Q_r0[i,j,k,rh,0]
#----------------------------------------------------------------------------------------------

def set_bc(Qxlow_ext, Qxlow_int, Qxhigh_ext, Qxhigh_int,
           Qylow_ext, Qylow_int, Qyhigh_ext, Qyhigh_int,
           Qzlow_ext, Qzlow_int, Qzhigh_ext, Qzhigh_int):
    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    # * r()
    #
    # EXT LIBS:
    # * np.min(), np.max()
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * Qxlow_ext, Qxlow_int, Qxhigh_ext, Qxhigh_int
    #  * Qylow_ext, Qylow_int, Qyhigh_ext, Qyhigh_int
    #  * Qzlow_ext, Qzlow_int, Qzhigh_ext, Qzhigh_int
    #  * nx, ny, nz
    #  * mx, my, mz
    #  * xlbc, ylbc, zlbc
    #  * nface, lxu, dx
    # LOCAL:
    #  * i, j, k, i4 (loop vars)
    #  * disk_rad (optional)
    ############################################################################

    #---------------------------------------------------------
    # Set BCs for x-boundaries
    if xlbc != 2:
        set_bc_x(Qxlow_ext, Qxlow_int, Qxhigh_ext, Qxhigh_int)

    #---------------------------------------------------------
    # Set BCs for x-boundaries
    if ylbc != 2:
        set_bc_y(Qylow_ext, Qylow_int, Qyhigh_ext, Qyhigh_int)

    #---------------------------------------------------------
    # Set BCs for x-boundaries
    if zlbc != 2:
        set_bc_z(Qzlow_ext, Qzlow_int, Qzhigh_ext, Qzhigh_int)


def set_bc_x(Qxlow_ext, Qxlow_int, Qxhigh_ext, Qxhigh_int):
    #--------------------------------------------------------
    # Set B.C.'s at bottom x-boundary
    for k in xrange(nz):
        for j in xrange(ny):

            for i4 in xrange(nface):
                Qxlow_ext[j,k,i4,:] = Qxlow_int[j,k,i4,:]

            if np.max(Qxlow_ext[j,k,0:nface,mx]) > 0.0:
                Qxlow_ext[j,k,0:nface,mx] = 0.

    #--------------------------------------------------------
    # Set B.C.'s at top x-boundary
    for k in xrange(nz):
        for j in xrange(ny):

            for i4 in xrange(nface):
                Qxhigh_ext[j,k,i4,:] = Qxhigh_int[j,k,i4,:]
            if np.min(Qxhigh_ext[j,k,0:nface,mx]) < 0.0:
                Qxhigh_ext[j,k,0:nface,mx] = 0.


def set_bc_y(Qylow_ext, Qylow_int, Qyhigh_ext, Qyhigh_int):
    #--------------------------------------------------------
    # Set B.C.'s at bottom y-boundary
    for k in xrange(nz):
        for i in xrange(nx):

            for i4 in xrange(nface):
                Qylow_ext[i,k,i4,:] = Qylow_int[i,k,i4,:]
            if np.max(Qylow_ext[i,k,0:nface,my]) > 0.0:
                Qylow_ext[i,k,0:nface,my] = 0.

    #--------------------------------------------------------
    # Set B.C.'s at top y-boundary
    for k in xrange(nz):
        for i in xrange(nx):

            for i4 in xrange(nface):
                Qyhigh_ext[i,k,i4,:] = Qyhigh_int[i,k,i4,:]
            if np.min(Qyhigh_ext[i,k,0:nface,my]) < 0.0:
                Qyhigh_ext[i,k,0:nface,my] = 0.


def set_bc_z(Qzlow_ext, Qzlow_int, Qzhigh_ext, Qzhigh_int):
    #--------------------------------------------------------
    # Set B.C.'s at bottom z-boundary
    for j in xrange(ny):
        for i in xrange(nx):

            if r(i,j) > disk_rad and r(i,j) < (lxu - 2.0*dx):
                for i4 in xrange(nface):
                    Qzlow_ext[i,j,i4,:] = Qzlow_int[i,j,i4,:]
            if np.max(Qzlow_ext[i,j,0:nface,mz]) > 0.0:
                Qzlow_ext[i,j,0:nface,mz] = 0.

    #-----------------------------------------------------------
    # Set B.C.'s at top z-boundary
    for j in xrange(ny):
        for i in xrange(nx):

            Qzhigh_ext[i,j,0:nface,:] = Qzhigh_int[i,j,0:nface,:]
            if np.min(Qzhigh_ext[i,j,0:nface,mz]) < 0.0:
                Qzhigh_ext[i,j,0:nface,mz] = 0

#----------------------------------------------------------------------------------------------


def advance_time_level_gl(Q_ri,Q_rp):
    # integer i,j,k,ieq,ir
    # real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri, Q_rp
    # real Q_xtemp, Q_ytemp, Q_ztemp, dti
    # real c1d3, c1d5

    c1d3 = 1./3.
    c1d5 = 1./5.
    dti = 1./dt

    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx):

                for ieq in xrange(nQ):
                    for ir in xrange(nbasis):
                        Q_rp[i,j,k,ieq,ir] = Q_ri[i,j,k,ieq,ir] - dt*glflux_r[i,j,k,ieq,ir] + dt*source_r[i,j,k,ieq,ir]

                for ieq in xrange(nQ):
                    if Q_rp[i,j,k,ieq,0] != Q_rp[i,j,k,ieq,0]:
                        print 'NaN. Bailing out...','  xc  =',xc[i],'  yc  =',yc[j],'  zc  =',zc[k],'  ieq  =',ieq
                        call exit(-1)

#----------------------------------------------------------------------------------------------

def source_calc(Q_ri,t):
    # integer i,j,k,ieq,ipg,ir
    # real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri
    # real t,source(npg,nQ),Qin(nQ),dn,dni,Zin,vx,vy,vz,alpha,temp,dne,eta_a,Teev,Tiev,etaJ2,nuei,TemR,Tev,vmax
    # real Tcoef,fac,en_floor,gyro,Pres,rh_buf

    source_r(:,:,:,:,:) = 0.0
    source(:,:) = 0.0
    en_floor = P_floor/(aindex - 1.)

    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx):

                for ipg in xrange(npg):

                    for ieq in xrange(nQ):
                        Qin[ieq] = sum(bfvals_int[ipg,0:nbasis]*Q_ri[i,j,k,ieq,0:nbasis])

                    # Sources for the fluid variables. Nonzero if randomly forced.

                    source[ipg,rh] = 0  # amplen*(ran(iseed) - 0.5)
                    source[ipg,mx] = 0  # amplm*(ran(iseed) - 0.5)
                    source[ipg,my] = 0  # amplm*(ran(iseed) - 0.5)
                    source[ipg,mz] = 0  # amplm*(ran(iseed) - 0.5)
                    source[ipg,en] = 0  # amplen*(ran(iseed) - 0.5)

                    # Sources for the viscous stress tensor.  The viscous stress is computed using hyperbolic
                    # relaxation to a parabolic problem. This is a generalization of the problem:
                    #
                    # partial_t u = partial_x v ,  eps * partial_t v - v = D * partial_x u
                    #
                    # where eps is a small parameter and D is a diffusion coefficient. In the relaxation limit
                    # eps -> 0, this becomes partial_t u = D * partial_x^2 u
                    # The following are the sources for this problem where epsi = 1/eps

                    source[ipg,pxx] = -epsi*Qin[pxx] #+ amplv*(ran(iseed) - 0.5)
                    source[ipg,pyy] = -epsi*Qin[pyy] #+ amplv*(ran(iseed) - 0.5)
                    source[ipg,pzz] = -epsi*Qin[pzz] #+ amplv*(ran(iseed) - 0.5)
                    source[ipg,pxy] = -epsi*Qin[pxy] #+ amplv*(ran(iseed) - 0.5)
                    source[ipg,pxz] = -epsi*Qin[pxz] #+ amplv*(ran(iseed) - 0.5)
                    source[ipg,pyz] = -epsi*Qin[pyz] #+ amplv*(ran(iseed) - 0.5)


                for ir in xrange(nbasis):
                    for ieq in xrange(nQ):
                        source_r[i,j,k,ieq,ir] = 0.125*cbasis[ir]*sum(bval_int_wgt[0:npg,ir]*source[0:npg,ieq])

end subroutine source_calc

#----------------------------------------------------------------------------------------------

def flux_calc_pnts_r(Qpnts_r,fpnts_r,ixyz,npnts):

    # Calculate the flux "fpnts_r" in direction "ixyz" (x, y, or z) at a set of
    # points corresponding to conserved quantities "Qpnts_r".

    # ixyz=0: x-direction
    # ixyz=1: y-direction
    # ixyz=2: z-direction

    # integer ife, ixyz,npnts
    # real, dimension(npnts,nQ):: Qpnts_r, fpnts_r
    # real dn,dni,vx,vy,vz,P,asqr,fac,Pre,dnei,Psol,dx2,Tem,smsq,nu,c2d3,c4d3
    # real ampx,ampy,ampz,ampd

    nu = epsi*vis
    c2d3 = 2./3.
    c4d3 = 4./3.

    ampx = 0*amplv*(ran(iseed) - 0.5)
    ampy = 0*amplv*(ran(iseed) - 0.5)
    ampz = 0*amplv*(ran(iseed) - 0.5)
    ampd = 0*amplen*(ran(iseed) - 0.5)

    for ife in xrange(npnts):

        dn = Qpnts_r[ife,rh]
        dni = 1./dn
        smsq = Qpnts_r[ife,mx]**2 + Qpnts_r[ife,my]**2 + Qpnts_r[ife,mz]**2

        P = (aindex - 1.)*(Qpnts_r[ife,en] - 0.5*dni*smsq)
        if ieos == 2:
            P = P_1*(dn**7.2 - 1.) + P_base + P
        if P < P_floor:
            P = P_floor

        vx = Qpnts_r[ife,mx]*dni
        vy = Qpnts_r[ife,my]*dni
        vz = Qpnts_r[ife,mz]*dni

        if ixyz == 0:

            fpnts_r[ife,rh] = Qpnts_r[ife,mx]

            fpnts_r[ife,mx] = Qpnts_r[ife,mx]*vx + P + Qpnts_r[ife,pxx] + ampx
            fpnts_r[ife,my] = Qpnts_r[ife,my]*vx + Qpnts_r[ife,pxy]     + ampy
            fpnts_r[ife,mz] = Qpnts_r[ife,mz]*vx + Qpnts_r[ife,pxz]     + ampz

            fpnts_r[ife,en] = (Qpnts_r[ife,en] + P)*vx - (Qpnts_r[ife,pxx]*vx + Qpnts_r[ife,pxy]*vy + Qpnts_r[ife,pxz]*vz)

            fpnts_r[ife,pxx] = c4d3*nu*vx
            fpnts_r[ife,pyy] = -c2d3*nu*vx
            fpnts_r[ife,pzz] = -c2d3*nu*vx

            fpnts_r[ife,pxy] = nu*vy
            fpnts_r[ife,pxz] = nu*vz
            fpnts_r[ife,pyz] = 0

        if ixyz == 1:

            fpnts_r(ife,rh) = Qpnts_r[ife,mxa[ixyz]]

            fpnts_r[ife,mx] = Qpnts_r[ife,mx]*vy + Qpnts_r[ife,pxy]     + ampx
            fpnts_r[ife,my] = Qpnts_r[ife,my]*vy + P + Qpnts_r[ife,pyy] + ampy
            fpnts_r[ife,mz] = Qpnts_r[ife,mz]*vy + Qpnts_r[ife,pyz]     + ampz

            fpnts_r[ife,en] = (Qpnts_r[ife,en] + P)*vy - (Qpnts_r[ife,pyy]*vy + Qpnts_r[ife,pxy]*vx + Qpnts_r[ife,pyz]*vz)

            fpnts_r[ife,pxx] = -c2d3*nu*vy
            fpnts_r[ife,pyy] = c4d3*nu*vy
            fpnts_r[ife,pzz] = -c2d3*nu*vy

            fpnts_r[ife,pxy] = nu*vx
            fpnts_r[ife,pxz] = 0
            fpnts_r[ife,pyz] = nu*vz

        if ixyz == 2:

            fpnts_r[ife,rh] = Qpnts_r[ife,mz]

            fpnts_r[ife,mx] = Qpnts_r[ife,mx]*vz + Qpnts_r[ife,pxz]     + ampx
            fpnts_r[ife,my] = Qpnts_r[ife,my]*vz + Qpnts_r[ife,pyz]     + ampy
            fpnts_r[ife,mz] = Qpnts_r[ife,mz]*vz + P + Qpnts_r[ife,pzz] + ampz

            fpnts_r[ife,en] = (Qpnts_r[ife,en] + P)*vz - (Qpnts_r[ife,pzz]*vz + Qpnts_r[ife,pxz]*vx + Qpnts_r[ife,pyz]*vy)

            fpnts_r[ife,pxx] = -c2d3*nu*vz
            fpnts_r[ife,pyy] = -c2d3*nu*vz
            fpnts_r[ife,pzz] = c4d3*nu*vz

            fpnts_r[ife,pxy] = 0.
            fpnts_r[ife,pxz] = nu*vx
            fpnts_r[ife,pyz] = nu*vy

end subroutine

#----------------------------------------------------------------------------------------------

def flux_calc_pnts_r2(Qpnts_r,fpnts_r,ixyz,npnts):

    # Calculate the flux "fpnts_r" in direction "ixyz" (x, y, or z) at a set of
    # points corresponding to conserved quantities "Qpnts_r".

    # ixyz=1: x-direction
    # ixyz=2: y-direction
    # ixyz=3: z-direction

    # integer ife, ixyz,npnts
    # real, dimension(npnts,nQ):: Qpnts_r, fpnts_r
    # real dn,dni,vr,P,asqr,fac,Pre,dnei,Psol,dx2,Tem,smsq

    for ife in xrange(npnts):

        dn = Qpnts_r[ife,rh]
        dni = 1./dn

        smsq = Qpnts_r[ife,mx]**2 + Qpnts_r[ife,my]**2 + Qpnts_r[ife,mz]**2
        vr = Qpnts_r[ife,mxa[ixyz]]*dni

        P = (aindex - 1.)*(Qpnts_r[ife,en] - 0.5*dni*smsq)
        if P < P_floor:
            P = P_floor

        fpnts_r[ife,rh] = Qpnts_r[ife,mxa[ixyz]]
        fpnts_r[ife,mx:mz+1] = Qpnts_r[ife,mx:mz+1]*vr
        fpnts_r[ife,en] = (Qpnts_r[ife,en] + P)*vr

        fpnts_r[ife,mxa[ixyz]] = fpnts_r[ife,mxa[ixyz]) + P

end subroutine

#----------------------------------------------------------------------------------------------

def flux_cal(Q_ri,
             bfvals_xm, bfvals_xp,
             bfvals_ym, bfvals_yp,
             bfvals_zm, bfvals_zp):
    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    #  * flux_hllc(), flux_roe(), flux_calc_pnts_r()
    #
    # EXT LIBS:
    #  * np.sum()
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * nx, ny, nz, nQ
    #  * rh, mx, my, mz, en
    #  * rh_floor, epsi, npge, nface
    #  * ihllc, iroe
    #  * kroe(nface)
    #  * bf_faces(nslim,nbastot)
    #  * Q_ri (pass-by-ref, shape (nx,ny,nz,nQ,nbasis))
    # LOCAL:
    #  * i, j, k, ieq, ipnt, i4 (loop vars)
    #  * iback, jleft, kdown
    #  * bfvals_zm, bfvals_zp  (shape (nface,nbastot))
    #  * bfvals_ym, bfvals_yp  (shape (nface,nbastot))
    #  * bfvals_xm, bfvals_xp  (shape (nface,nbastot))
    #  * Qface_x, Qface_y, Qface_z  (shape (nfe,nQ))
    #  * fface_x, fface_y, fface_z  (shape (nfe,nQ))
    #  * fhllc_x ,fhllc_y ,fhllc_z  (shape (nface,5))
    #  * cwavex, cwavey, cwavez  (shape (nfe))
    #  * qvin (shape (nQ))
    ############################################################################

    kroe[:] = 1

    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx+1):

                iback = i-1

                if i > 1:
                    for ieq in xrange(nQ):
                        for ipnt in xrange(nface):
                            Qface_x[ipnt,ieq] = sum(bfvals_xp[ipnt,0:nbasis]*Q_ri[iback,j,k,ieq,0:nbasis])

                if i == 1:
                    for ieq in xrange(nQ):
                        Qface_x[0:nface,ieq] = Qxlow_ext[j,k,0:nface,ieq]

                if i < nx+1:
                    for ieq in xrange(nQ):
                        for ipnt in xrange(nface):
                            Qface_x[ipnt+nface,ieq] = sum(bfvals_xm[ipnt,0:nbasis]*Q_ri[i,j,k,ieq,0:nbasis])

                if i == nx+1:
                    for ieq in xrange(nQ):
                        Qface_x[nface:nfe,ieq] = Qxhigh_ext[j,k,0:nface,ieq]

                call flux_calc_pnts_r(Qface_x,fface_x,0,nfe)

                if iroe != 1 and ihllc != 1:

                    for i4 in xrange(nfe):
                        for ieq in xrange(nQ):
                            qvin[ieq] = Qface_x[i4,ieq]
                        cwavex[i4] = cfcal[qvin,0]

                    for i4 in xrange(nface):
                        cfrx[i4,rh:en] = max(cwavex[i4],cwavex[i4+nface])

                for ieq in xrange(nQ):
                    for i4 in xrange(nface):
                        i4p = i4 + nface
                        flux_x[i4,i,j,k,ieq] = 0.5*(fface_x[i4,ieq] + fface_x[i4p,ieq]) - 0.5*cfrx[i4,ieq]*(Qface_x[i4p,ieq] - Qface_x[i4,ieq])

                kroe[0:nface] = 1

                if ihllc == 1:
                    call flux_hllc(Qface_x,fface_x,fhllc_x,0)
                if iroe == 1:
                    call flux_roe(Qface_x,fface_x,fhllc_x,0)

                if ihllc == 1 or iroe == 1:
                    for ieq in xrange(en):
                        for i4 in xrange(nface):
                            if kroe[i4] > 0:
                                flux_x[i4,i,j,k,ieq] = fhllc_x[i4,ieq]

#----------------------------------------------------

    for k in xrange(nz):
        for j in xrange(ny+1):
            jleft = j-1
            for i in xrange(nx):

                if j > 1:
                    for ieq in xrange(nQ):
                        for ipnt in xrange(nface):
                              Qface_y[ipnt,ieq] = sum(bfvals_yp[ipnt,0:nbasis]*Q_ri[i,jleft,k,ieq,0:nbasis])
                if j == 1:
                    for ieq in xrange(nQ):
                       Qface_y[0:nface,ieq] = Qylow_ext[i,k,0:nface,ieq]

                if j < ny+1:
                    for ieq in xrange(nQ):
                        for ipnt in xrange(nface):
                            Qface_y[ipnt+nface,ieq] = sum(bfvals_ym[ipnt,0:nbasis]*Q_ri[i,j,k,ieq,0:nbasis])
                if j == ny+1:
                    for ieq in xrange(nQ):
                       Qface_y[nface:nfe,ieq] = Qyhigh_ext[i,k,0:nface,ieq]

                call flux_calc_pnts_r(Qface_y,fface_y,1,nfe)

                if iroe != 1 and ihllc != 1:

                    for i4 in xrange(nfe):
                        for ieq in xrange(nQ):
                            qvin[ieq] = Qface_y[i4,ieq]
                        cwavey[i4] = cfcal[qvin,1]

                    for i4 in xrange(nface):
                        cfry[i4,rh:en] = max(cwavey[i4],cwavey[i4+nface])

                for ieq in xrange(nQ):
                    for i4 in xrange(nface):
                        i4p = i4 + nface
                        flux_y[i4,i,j,k,ieq] = 0.5*(fface_y[i4,ieq] + fface_y[i4p,ieq]) - 0.5*cfry[i4,ieq]*(Qface_y[i4p,ieq] - Qface_y[i4,ieq])

                kroe[0:nface] = 1

                if ihllc == 1:
                    call flux_hllc(Qface_y,fface_y,fhllc_y,1)
                if iroe == 1:
                    call flux_roe(Qface_y,fface_y,fhllc_y,1)

                if ihllc == 1 or iroe == 1:
                    for ieq in xrange(en):
                        for i4 in xrange(nface):
                            if kroe[i4] > 0:
                                flux_y[i4,i,j,k,ieq] = fhllc_y[i4,ieq]

#-------------------------------------------------------

    for k in xrange(nz+1):
        kdown = k-1
        for j in xrange(ny):
            for i in xrange(nx):

                if k > 1:
                    for ieq in xrange(nQ):
                        for ipnt in xrange(nface):
                            Qface_z[ipnt,ieq] = sum(bfvals_zp[ipnt,0:nbasis]*Q_ri[i,j,kdown,ieq,0:nbasis])
                if k == 1:
                    for ieq in xrange(nQ):
                        Qface_z[0:nface,ieq] = Qzlow_ext[i,j,0:nface,ieq]

                if k < nz+1:
                    for ieq in xrange(nQ):
                        for ipnt in xrange(nface):
                            Qface_z[ipnt+nface,ieq] = sum(bfvals_zm[ipnt,0:nbasis]*Q_ri[i,j,k,ieq,0:nbasis])
                if k == nz+1:
                    for ieq in xrange(nQ):
                        Qface_z[nface:nfe,ieq] = Qzhigh_ext[i,j,0:nface,ieq]

                call flux_calc_pnts_r(Qface_z,fface_z,2,nfe)

                if iroe != 1 and ihllc != 1:
                    for i4 in xrange(nfe):
                        for ieq in xrange(nQ):
                            qvin[ieq] = Qface_z[i4,ieq]
                        cwavez[i4] = cfcal[qvin,2]

                    for i4 in xrange(nface):
                        cfrz[i4,rh:en] = max(cwavez[i4],cwavez[i4+nface])

                for ieq in xrange(nQ):
                    for i4 in xrange(nface):
                        i4p = i4 + nface
                        flux_z[i4,i,j,k,ieq] = 0.5*(fface_z[i4,ieq] + fface_z[i4p,ieq]) - 0.5*cfrz[i4,ieq]*(Qface_z[i4p,ieq] - Qface_z[i4,ieq])

                kroe[0:nface] = 1

                if ihllc == 1:
                    call flux_hllc(Qface_z,fface_z,fhllc_z,2)
                if iroe == 1:
                    call flux_roe(Qface_z,fface_z,fhllc_z,2)

                if ihllc == 1 or iroe == 1:
                    for ieq in xrange(en):
                        for i4 in xrange(nface):
                        if kroe[i4] > 0:
                            flux_z[i4,i,j,k,ieq] = fhllc_z[i4,ieq]

#----------------------------------------------------------------------------------------------

def flux_hllc(Qlr,flr,fhllc,ixyz):

#   Compute ion or electron fluxes (density, momentum, energy) using
#   nonrelativistic HD HLLC approximate Riemann solver developed by Batten, 1997,
#   *****    "On the Choice of Wavespeeds for the HLLC Riemann Solver"

#   This takes into account the left and right propagating shocks along with
#   the contact (tangential) discontinuity.

    # real Qlr(nfe,nQ),flr(nfe,nQ),fhllc(nface,5)
    # real sm_num(nface),sm_den(nface),qtilde(nface,5),rtrho(nfe),rtrho_i(nface),qsq(nfe)
    # real s_lr(nfe),ctilde(nface),hlr(nfe),cslr(nfe),ctsq(nface),Zi, mfact, csfac
    # real aq(nfe),bq(nfe),Qstar(nfe,6),fstar(nfe,6),pstar(nface),s_m(nface)
    # real rhov(nfe),vlr(nfe,3),plr(nfe),slsm_i(nface),rho_i,qslr(nfe),sq_lr(nfe),slrm_i(nfe),B2(nfe),cf(nfe)
    # integer ixyz,i4,i4p1,nr,jie,k,k2,ieq,iparr,iperp1,iperp2,ibatten
    # integer rhj,mxj,myj,mzj,enj,psj,ivar(5),ipassive,nhll,ib1,ib2

    nhll = 5
    rhj = rh
    mxj = mx
    myj = my
    mzj = mz
    enj = en

    iparr  = mxa[ixyz]
    iperp1 = mya[ixyz]
    iperp2 = mza[ixyz]

    ivar[0] = rhj
    ivar[1] = iparr
    ivar[2] = iperp1
    ivar[3] = iperp2
    ivar[4] = enj

    # Indices 1:4 denote the left state, while indices 5:8 denote the right state.

   ibatten = 0

    for k in xrange(nfe):
        rhov[k] = Qlr[k,rhj]
        rho_i = 1.0/rhov[k]
        vlr[k,0] = Qlr[k,iparr]*rho_i     # velocity parallel to direction of flux computation
        vlr[k,1] = Qlr[k,iperp1]*rho_i        # velocity in perpendicular direction 1
        vlr[k,2] = Qlr[k,iperp2]*rho_i        # velocity in perpendicular direction 2
        qsq[k] = vlr[k,0]**2 + vlr[k,1]**2 + vlr[k,2]**2
        plr[k] = (aindex - 1.0)*(Qlr(k,enj) - 0.5*rhov[k]*qsq[k])     # pressure
         if ieos == 2:
             plr[k] = P_1*(rhov[k]**7.2 - 1.) + P_base + plr[k]
        rtrho[k] = sqrt(rhov[k])

    for k in xrange(nface):
        k2 = k + nface
        if ieos == 2:
            cslr[k] = vlr[k,0] - sqrt(7.2*P_1*rhov[k]**6.2 + plr[k]*rho_i)     # lambda_M(Q_l)
            cslr[k2] = vlr[k2,0] + sqrt(7.2*P_1*rhov[k2]**6.2 + plr[k2]*rho_i)     # lambda_P(Q_ri)
        else:
            cslr[k] = vlr[k,0] - sqrt(aindex*plr[k]/rhov[k])       # lambda_M(Q_l)
            cslr[k2] = vlr[k2,0] + sqrt(aindex*plr[k2]/rhov[k2] )      # lambda_P(Q_r)

    if ibatten == 1:     # compute wave speeds using Roe averages following Batten, 1997

        for k in xrange(nface):
            k2 = k + nface
            rtrho_i[k] = 1.0/(rtrho[k] + rtrho[k2])
            qtilde[k,0] = (rtrho[k]*vlr[k,0] + rtrho[k2]*vlr[k2,0])*rtrho_i[k]
            qtilde[k,1] = (rtrho[k]*vlr[k,1] + rtrho[k2]*vlr[k2,1])*rtrho_i[k]
            qtilde[k,2] = (rtrho[k]*vlr[k,2] + rtrho[k2]*vlr[k2,2])*rtrho_i[k]
            qsq[k] = qtilde[k,0]**2 + qtilde[k,1]**2 + qtilde[k,2]**2
            hlr[k] = (Qlr[k,enj] + plr[k])/rhov[k]
            hlr[k2] = (Qlr[k2,enj] + plr[k2])/rhov[k2]
            qtilde[k,3] = (rtrho[k]*hlr[k] + rtrho[k2]*hlr[k2])*rtrho_i[k]
            ctsq(k) = (aindex - 1.0)*(qtilde[k,3] - 0.5*qsq[k])
        if minval(ctsq) > 0.0:
            ctilde = sqrt(ctsq)
            qslr(1:nface) = qtilde(1:nface,1) - ctilde(1:nface)     # lambda_M(Q_Roe)
            qslr(nface+1:nfe) = qtilde(nface+1:nfe,1) + ctilde(nface+1:nfe)    # lambda_P(Q_Roe)
        if minval(ctsq) < 0.0:
            ibatten = 0

    if ibatten == 0:
        for k in xrange(nface):
            k2 = k + nface
            if ieos == 2:
                qslr[k] = vlr[k2,0] - sqrt(7.2*P_1*rhov[k2]**6.2 + plr[k2]*rho_i)      # lambda_M(Q_r)
                qslr[k2] = vlr[k,0] + sqrt(7.2*P_1*rhov[k]**6.2 + plr[k]*rho_i)    # lambda_P(Q_l)
            else:
                qslr[k] = vlr[k2,0] - sqrt(aindex*plr[k2]/rhov[k2])    # lambda_M(Q_r)
                qslr[k2] = vlr[k,0] + sqrt(aindex*plr[k]/rhov[k])      # lambda_P(Q_l)

#     Calculate the slow and fast wavespeeds S_L, S_R of the left, right propagating shocks.

    for k in xrange(nface):
        k2 = k + nface
        s_lr[k] = min(cslr[k],qslr[k])         # S_L = min(lambda_M(Q_l),lambda_M(Q_r or Q_Roe))
        s_lr[k2] = max(cslr[k2],qslr[k2])      # S_R = max(lambda_P(Q_r),lambda_P(Q_l or Q_Roe))
        sm_num[k] = rhov[k2]*vlr[k2,0]*(s_lr[k2] - vlr[k2,0]) - rhov[k]*vlr[k,0]*(s_lr[k] - vlr[k,0])
        sm_num[k] = sm_num[k] + plr[k] - plr[k2]
        sm_den[k] = rhov[k2]*(s_lr[k2] - vlr[k2,0]) - rhov[k]*(s_lr[k] - vlr[k,0])
    where (sm_den .eq. 0.0) sm_den = rh_floor

#      Calculate the wavespeed S_M of the contact discontinuity.

    for k in xrange(nface):
        s_m[k] = sm_num[k]/sm_den[k]          # Eq. (34) of Batten, 1997
        pstar[k] = rhov[k]*(vlr[k,0] - s_lr[k])*(vlr[k,0] - s_m[k]) + plr[k]
                        # Eq. (36) of Batten, 1997

#      Now, calculate Q_l* and Q_r* in order to calculate F_l* and F_r*.

    for k in xrange(nfe):
        if k <= nface:
            i4 = k
        if k > nface:
            i4 = k - nface
        sm_den(1) = s_lr[k] - s_m[i4]     # S_{L,R} - S_M

        if sm_den(1) == 0.0:
            sm_den(1) = rh_floor

        slrm_i[k] = 1.0/sm_den(1)
        sq_lr[k] = s_lr[k] - vlr[k,0]         # S_{L,R} - q_{l,r}
        Qstar[k,0] = rhov[k]*sq_lr[k]*slrm_i[k]       # Eq. (35) of Batten ,1997
        Qstar[k,1] = (sq_lr[k]*Qlr[k,iparr] + pstar[i4] - plr[k])*slrm_i[k]   # Eq. (37-39) of Batten, 1997
        Qstar[k,2] = sq_lr[k]*Qlr[k,iperp1]*slrm_i[k]     # Eq. (37-39) of Batten, 1997
        Qstar[k,3] = sq_lr[k]*Qlr[k,iperp2]*slrm_i[k]     # Eq. (37-39) of Batten, 1997
        Qstar[k,4] = (sq_lr[k]*Qlr[k,enj] - plr[k]*vlr[k,0] + pstar[i4]*s_m[i4])*slrm_i[k]
                            # Eq. (40) of Batten, 1997

        for ieq in xrange(nhll):
            fstar[k,ieq] = flr[k,ivar[ieq]] + s_lr[k]*(Qstar[k,ieq] - Qlr[k,ivar[ieq]])
                                # Eq. (29) of Batten, 1997

    # Finally, calculate the HLLC fluxes from F_l, F_l*, F_r*, and F_r...

    for i4 in xrange(nface):                     # Use Eq. (26) of Batten ,1997
        if s_lr[i4] > 0.0:             # if S_L > 0
            for ieq in xrange(nhll):
                fhllc[i4,ivar[ieq]] = flr[i4,ivar[ieq]]     # F_HLLC = F_l
        if s_lr[i4] <= 0.0 and 0.0 .lt. s_m[i4]:      # if S_L <= 0 < S_M
            for ieq in xrange(nhll):
                fhllc[i4,ivar[ieq]] = fstar[i4,ieq]        # F_HLLC = F_l*
        if s_m[i4] <= 0.0 and 0.0 .le. s_lr[i4+nface]:        # if S_M <= 0 <= S_R
            for ieq in xrange(nhll):
                fhllc[i4,ivar[ieq]] = fstar[i4+nface,ieq]       # F_HLLC = F_r*
        if s_lr[i4+nface] < 0.0:           # if S_R < 0
            for ieq in xrange(nhll):
                fhllc[i4,ivar[ieq]] = flr[i4+nface,ivar[ieq]]   # F_HLLC = F_r

end subroutine flux_hllc

    #----------------------------------------------------------------

def flux_roe(Qlr,flr,froef,ixyz):
    # integer ixyz,i9,j9,ilr,iparr,iperp
    # real Qlr(nfe,nQ),flr(nfe,nQ),froef(nface,5)
    # real dflux(nface,5),dwr(nface,5),dsk(nface),Qri
    # real evec(5,5),swr(2,nface,5),skr(2,nface),asq,lam9(5),a9(5),sk9
    # real a9e,ea2i,vsq,Pri,vels(3),dnii,en_floor
    # real vc(3),hc,dele(nface),delrho(nface),delmom1(nface),delmom2(nface),delmom3(nface),sc9,vcsq,vcdel,csq,gmi,bsq


    # approximate Roe solver for non-relativistic two-fluid, from Eulderink and Mellema, 1994

    # Starting from HD equations for a given fluid (either electrons or ions), this routine
    # computes "dflux_i = (cQ)_i+1 - (cQ)_i" terms by characteristic wave decomposition using
    # a Roe matrix.

    # INPUT is Qr(,1)=conserved density, Qr(,2)=energy, Qr(,3:5)=momentum
    #      Qpr(,1)=density, Qpr(,2)=pressure, Qpr(,3:5)=velocity
    #      n = number of spatial points in given direction (x or y)

    #      ixyz = 0: flux is being computed in x-direction
    #      ixyz = 2: flux is being computed in y-direction


    # OUTPUT is  "dflux_i = (cQ)_i+1 - (cQ)_i"
    #   dflux(,1)=density term, dflux(,2)=energy term, dflux(,3:5)=momentum terms



    # evec(,1) through evec(,5) are the eigenvectors of the Roe matrix.

    en_floor = P_floor*aindm1
    gmi = 1./(aindex - 1.)

    evec(:,:) = 0.
    evec(1,1) = 1.
    evec(1,2) = 1.
    evec(1,3) = 1.
    evec(4,4) = 1.
    evec(5,5) = 1.


    for i9 in xrange(nface):
        kroe[i9] = 1
        do j9=1,2
            if j9 == 1:
                ilr = i9
            if j9 == 2:
                ilr = i9 + nface

            # Make sure density, pressure, and energy are at least at floor values.

            dnii = 1./Qlr[ilr,rh]
            vels[0] = Qlr[ilr,mxa[ixyz]]*dnii
            vels[1] = Qlr[ilr,mya[ixyz]]*dnii
            vels[2] = Qlr[ilr,mza[ixyz]]*dnii
            vsq = vels[0]**2 + vels[1]**2 + vels[2]**2
            Pri = aindm1*(Qlr[ilr,en] - 0.5*Qlr[ilr,rh]*vsq)

            skr[j9,i9] = sqrt(Qlr[ilr,rh])
            swr[j9,i9,0] = skr[j9,i9]*vels[0]  ! sqrt(rho) * v_x
            swr[j9,i9,1] = skr[j9,i9]*vels[1]  ! sqrt(rho) * v_y
            swr[j9,i9,2] = skr[j9,i9]*vels[2]  ! sqrt(rho) * v_z
            swr[j9,i9,3] = 0.5*skr[j9,i9]*(vsq + cp*Pri/Qlr[ilr,rh])
        end do

    for i9 in xrange(nface):
        Qri = 1./Qlr[i9,rh]    # Increments in conserved quantities are normalized w.r.t. density.

        delrho[i9] = Qlr[i9+nface,rh]*Qri - 1. # delrho = increment in conserved density
        dele[i9] = (Qlr[i9+nface,en] - Qlr[i9,en])*Qri    # *one_mime(jie)
        # dele = increment in conserved energy

        delmom1[i9] = (Qlr[i9+nface,mxa[ixyz]] - Qlr[i9,mxa[ixyz]])*Qri    # del1 = increment in x-momentum
        delmom2[i9] = (Qlr[i9+nface,mya[ixyz]] - Qlr[i9,mya[ixyz]])*Qri    # del2 = increment in y-momentum
        delmom3[i9] = (Qlr[i9+nface,mza[ixyz]] - Qlr[i9,mza[ixyz]])*Qri    # del3 = increment in z-momentum

        #    dwr(i,1:3) = 0.5*[sqrt(rho_i)v_{1:3,i} + sqrt(rho_{i+1})v_{1:3,i+1}]
        #    dwr(i,4) = 0.5*[sqrt(rho_i) enthalpy(i)/rho(i) + sqrt(rho_{i+1}) enthalpy(i+1)/rho(i+1)]

        for j9 in xrange(4):
            dwr[i9,j9] = 0.5*(swr[1,i9,j9] + swr[0,i9,j9])

        dsk[i9] = 0.5*(skr[1,i9] + skr[0,i9])

    # The Roe average of a quantity is the arithmetic average
    # between neighboring cells weighted by the square root of density.
    # For example, for "v_x" at position "i" the Roe average "v_{cx}" is

    #    v_{cx} = [sqrt(rho_i)v_{xi} + sqrt(rho_{i+1}) v_{x,i+1}]/[sqrt(rho_i) + sqrt(rho_{i+1})]

    for i9 in xrange(nface):
        vc[0] = dwr(i9,1)/dsk(i9) # component 1 of Roe-averaged velocity (x if jie=1, y if jie=2)
        vc[1] = dwr(i9,2)/dsk(i9) # component 2 of Roe-averaged velocity (y if jie=1, x if jie=2)
        vc[2] = dwr(i9,3)/dsk(i9) # component 3 of Roe-averaged velocity (z-component)
        hc = dwr(i9,4)/dsk(i9)    # Roe-averaged enthalpy/density

        vcsq = vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]
        asq = aindm1*(hc - 0.5*vcsq)  # squared sound speed
        if asq <= 0.0:
            kroe(i9) = 0    # asq = (aindex - 1.0d0)*hc
        sc9 = sqrt(asq)   # sound speed


        # Define the characteristic speeds (eigenvalues of the Roe matrix).
        lam9[0] = abs(vc[0] - sc9)
        lam9[1] = abs(vc[0] + sc9)
        lam9[2] = abs(vc[0])
        lam9[3] = lam9[3]
        lam9[4] = lam9[3]

        # Define the eigenvectors evec(,1)...evec(,5) of the Roe matrix.
        evec[1,0] = hc - sc9*vc[0]
        evec[2,0] = vc[0] - sc9
        evec[3,0] = vc[1]
        evec[4,0] = vc[2]

        evec[1,1] = hc + sc9*vc[0]
        evec[2,1] = vc[0] + sc9
        evec[3,1] = vc[1]
        evec[4,1] = vc[2]

        evec[1,2] = 0.5*vcsq
        evec[2,2] = vc[0]
        evec[3,2] = vc[1]
        evec[4,2] = vc[2]
        evec[1,3] = vc[1]
        evec[1,4] = vc[2]

        # Define a few intermediate variables needed for computing the expansion coefficients.

        ea2i = aindm1/asq
        vcdel = vc[0]*delmom1[i9] + vc[1]*delmom2[i9] + vc[2]*delmom3[i9] - dele[i9]
        a9e = 0.5*ea2i*(0.5*vcsq*delrho[i9] - vcdel)
        sk9 = 0.5*(delmom1[i9] - vc[0]*delrho[i9])/sc9

        # Define the expansion coefficients a9_1...a9_5 such that
        # Delta Q = {delrho,dele,del1,del2,del3} = sum [a9_j evec(,j)]

        a9[0] = a9e - sk9
        a9[1] = a9e + sk9
        a9[2] = ea2i*((hc - vcsq)*delrho[i9] + vcdel)
        a9[3] = delmom2[i9] - vc[1]*delrho[i9]
        a9[4] = delmom3[i9] - vc[2]*delrho[i9]

        # The flux increments "dflux" are now given by   Delta F = sum [a9_j lam_j evec(,j)]

        dflux[i9,0:5] = 0.
        do j9=1,5
            dflux[i9,0] = dflux[i9,0] + a9[j9]*lam9[j9]*evec[0,j9]
            dflux[i9,1] = dflux[i9,1] + a9[j9]*lam9[j9]*evec[1,j9]
            dflux[i9,2] = dflux[i9,2] + a9[j9]*lam9[j9]*evec[2,j9]
            dflux[i9,3] = dflux[i9,3] + a9[j9]*lam9[j9]*evec[3,j9]
            dflux[i9,4] = dflux[i9,4] + a9[j9]*lam9[j9]*evec[4,j9]
        enddo
        dflux[i9,0] = dflux[i9,0]*Qlr[i9,rh]   # flux increment in density
        dflux[i9,1] = dflux[i9,1]*Qlr[i9,rh]   # flux increment in energy
        dflux[i9,2] = dflux[i9,2]*Qlr[i9,rh]   # flux increment in parallel momentum
        dflux[i9,3] = dflux[i9,3]*Qlr[i9,rh]   # flux increment in perpendicular momentum
        dflux[i9,4] = dflux[i9,4]*Qlr[i9,rh]   # flux increment in z-momentum

        if kroe(i9) > 0:
            froef[i9,rh] = 0.5*(flr[i9,rh] + flr[i9+nface,rh] - dflux[i9,0])
            froef[i9,en] = 0.5*(flr[i9,en] + flr[i9+nface,en] - dflux[i9,1])
            froef[i9,mxa[ixyz]] = 0.5*(flr[i9,mxa[ixyz]] + flr[i9+nface,mxa[ixyz]] - dflux[i9,2])
            froef[i9,mya[ixyz]] = 0.5*(flr[i9,mya[ixyz]] + flr[i9+nface,mya[ixyz]] - dflux[i9,3])
            froef[i9,mza[ixyz]] = 0.5*(flr[i9,mza[ixyz]] + flr[i9+nface,mza[ixyz]] - dflux[i9,4])

end subroutine flux_roe

#----------------------------------------------------------------------------------------------

def innerintegral(Q_r):
    # real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
    # integer i,j,k,ieq,ipg,ir
    # real Qinner(npg,nQ),finner_x(npg,nQ), finner_y(npg,nQ), finner_z(npg,nQ), int_r(nbastot,nQ)

    integral_r(:,:,:,:,:) = 0.

    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx):

                for ieq in xrange(nQ):
                    for ipg in xrange(npg):
                        Qinner[ipg,ieq] = sum(bfvals_int[ipg,0:nbasis]*Q_r[i,j,k,ieq,0:nbasis])

                call flux_calc_pnts_r(Qinner,finner_x,1,npg)
                call flux_calc_pnts_r(Qinner,finner_y,2,npg)
                call flux_calc_pnts_r(Qinner,finner_z,3,npg)

                for ieq in xrange(nQ):

                    int_r[kx,ieq] = 0.25*cbasis[kx]*dxi*sum(wgt3d[0:npg]*finner_x[0:npg,ieq])
                    int_r[ky,ieq] = 0.25*cbasis[ky]*dyi*sum(wgt3d[0:npg]*finner_y[0:npg,ieq])
                    int_r[kz,ieq] = 0.25*cbasis[kz]*dzi*sum(wgt3d[0:npg]*finner_z[0:npg,ieq])

                    if nbasis > 4:
                        int_r[kyz,ieq]  = 0.25*cbasis[kyz]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kz]*finner_y[0:npg,ieq])  \
                              + 0.25*cbasis[kyz]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,ky]*finner_z[0:npg,ieq])
                        int_r[kzx,ieq]  = 0.25*cbasis[kzx]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kz]*finner_x[0:npg,ieq])  \
                              + 0.25*cbasis[kzx]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kx]*finner_z[0:npg,ieq])
                        int_r[kxy,ieq]  = 0.25*cbasis[kxy]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,ky]*finner_x[0:npg,ieq])  \
                              + 0.25*cbasis[kxy]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kx]*finner_y[0:npg,ieq])
                        int_r[kxyz,ieq] = 0.25*cbasis[kxyz]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kyz]*finner_x[0:npg,ieq])  \
                              + 0.25*cbasis[kxyz]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kzx]*finner_y[0:npg,ieq])  \
                              + 0.25*cbasis[kxyz]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kxy]*finner_z[0:npg,ieq])

                    if nbasis > 8:
                        int_r[kyzz,ieq] = 0.25*cbasis[kyzz]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kzz]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kyzz]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kyz]*finner_z[0:npg,ieq])
                        int_r[kzxx,ieq] = 0.25*cbasis[kzxx]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kxx]*finner_z[0:npg,ieq])  \
                          + 0.25*cbasis[kzxx]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kzx]*finner_x[0:npg,ieq])
                        int_r[kxyy,ieq] = 0.25*cbasis[kxyy]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kyy]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kxyy]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxy]*finner_y[0:npg,ieq])
                        int_r[kyyz,ieq] = 0.25*cbasis[kyyz]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kyy]*finner_z[0:npg,ieq])  \
                          + 0.25*cbasis[kyyz]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kyz]*finner_y[0:npg,ieq])
                        int_r[kzzx,ieq] = 0.25*cbasis[kzzx]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kzz]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kzzx]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kzx]*finner_z[0:npg,ieq])
                        int_r[kxxy,ieq] = 0.25*cbasis[kxxy]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kxx]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kxxy]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxy]*finner_x[0:npg,ieq])
                        int_r[kyyzz,ieq] = 0.25*cbasis[kyyzz]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kyzz]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kyyzz]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kyyz]*finner_z[0:npg,ieq])
                        int_r[kzzxx,ieq] = 0.25*cbasis[kzzxx]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kzxx]*finner_z[0:npg,ieq])  \
                          + 0.25*cbasis[kzzxx]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kzzx]*finner_x[0:npg,ieq])
                        int_r[kxxyy,ieq] = 0.25*cbasis[kxxyy]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxyy]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kxxyy]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxxy]*finner_y[0:npg,ieq])
                        int_r[kyzxx,ieq] = 0.25*cbasis[kyzxx]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxyz]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kyzxx]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kzxx]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kyzxx]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kxxy]*finner_z[0:npg,ieq])
                        int_r[kzxyy,ieq] = 0.25*cbasis[kzxyy]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxyz]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kzxyy]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kxyy]*finner_z[0:npg,ieq])  \
                          + 0.25*cbasis[kzxyy]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kyyz]*finner_x[0:npg,ieq])
                        int_r[kxyzz,ieq] = 0.25*cbasis[kxyzz]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxyz]*finner_z[0:npg,ieq])  \
                          + 0.25*cbasis[kxyzz]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kyzz]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kxyzz]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kzzx]*finner_y[0:npg,ieq])
                        int_r[kxyyzz,ieq] = 0.25*cbasis[kxyyzz]*dxi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kyyzz]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kxyyzz]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxyzz]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kxyyzz]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kzxyy]*finner_z[0:npg,ieq])
                        int_r[kyzzxx,ieq] = 0.25*cbasis[kyzzxx]*dyi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kzzxx]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kyzzxx]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kyzxx]*finner_z[0:npg,ieq])  \
                          + 0.25*cbasis[kyzzxx]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxyzz]*finner_x[0:npg,ieq])
                        int_r[kzxxyy,ieq] = 0.25*cbasis[kzxxyy]*dzi*sum(wgt3d[0:npg]*bfvals_int[0:npg,kxxyy]*finner_z[0:npg,ieq])  \
                          + 0.25*cbasis[kzxxyy]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kzxyy]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kzxxyy]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kyzxx]*finner_y[0:npg,ieq])
                        int_r[kxxyyzz,ieq] = 0.25*cbasis[kxxyyzz]*dxi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kxyyzz]*finner_x[0:npg,ieq])  \
                          + 0.25*cbasis[kxxyyzz]*dyi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kyzzxx]*finner_y[0:npg,ieq])  \
                          + 0.25*cbasis[kxxyyzz]*dzi*sum(wgt3d[0:npg]*3.*bfvals_int[0:npg,kzxxyy]*finner_z[0:npg,ieq])

                for ieq in xrange(nQ):
                    for ir in xrange(nbasis):
                        integral_r[i,j,k,ieq,ir] = int_r[ir,ieq]

end subroutine innerintegral


#----------------------------------------------------------------------------------------------

def glflux():
    # integer i,j,k,ieq,ir,iqfa

    for ieq in xrange(nQ):
        do k = 1,nz
        do j = 1,ny

            do i = 1,nx
                sumx = 0.
                sumy = 0.
                sumz = 0.

                do iqfa = 1,nface
                    sumx = sumx + 0.25*dxi*wgt2d(iqfa)*(flux_x(iqfa,i+1,j,k,ieq) - flux_x(iqfa,i,j,k,ieq))
                    sumy = sumy + 0.25*dyi*wgt2d(iqfa)*(flux_y(iqfa,i,j+1,k,ieq) - flux_y(iqfa,i,j,k,ieq))
                    sumz = sumz + 0.25*dzi*wgt2d(iqfa)*(flux_z(iqfa,i,j,k+1,ieq) - flux_z(iqfa,i,j,k,ieq))
                end do

                glflux_r(i,j,k,ieq,1) = sumx + sumy + sumz

            end do

            do ir=2,nbasis
                do i = 1,nx
                    sumx = 0.
                    sumy = 0.
                    sumz = 0.

                    do iqfa = 1,nface
                        sumx = sumx + wgtbf_xmp(iqfa,1,ir)*flux_x(iqfa,i,j,k,ieq) + wgtbf_xmp(iqfa,2,ir)*flux_x(iqfa,i+1,j,k,ieq)
                        sumy = sumy + wgtbf_ymp(iqfa,1,ir)*flux_y(iqfa,i,j,k,ieq) + wgtbf_ymp(iqfa,2,ir)*flux_y(iqfa,i,j+1,k,ieq)
                        sumz = sumz + wgtbf_zmp(iqfa,1,ir)*flux_z(iqfa,i,j,k,ieq) + wgtbf_zmp(iqfa,2,ir)*flux_z(iqfa,i,j,k+1,ieq)
                    end do

                    glflux_r(i,j,k,ieq,ir) = sumx + sumy + sumz - integral_r(i,j,k,ieq,ir)

                end do
            end do

        end do
        end do

end subroutine glflux

#----------------------------------------------------------------------------------------------

def glflux2():
    # integer i,j,k,ieq,ir

    for ieq in xrange(nQ):
        do k = 1,nz
        do j = 1,ny

            do i = 1,nx

                glflux_r(i,j,k,ieq,1) =  0.25*(dxi*(wgt2d[0]*(flux_x(1,i+1,j,k,ieq) - flux_x(1,i,j,k,ieq)))  \
                                + dyi*(wgt2d[0]*(flux_y(1,i,j+1,k,ieq) - flux_y(1,i,j,k,ieq)))  \
                                + dzi*(wgt2d[0]*(flux_z(1,i,j,k+1,ieq) - flux_z(1,i,j,k,ieq)))  \
                                                + dxi*(wgt2d[1]*(flux_x(2,i+1,j,k,ieq) - flux_x(2,i,j,k,ieq)))  \
                                + dyi*(wgt2d[1]*(flux_y(2,i,j+1,k,ieq) - flux_y(2,i,j,k,ieq)))  \
                                + dzi*(wgt2d[1]*(flux_z(2,i,j,k+1,ieq) - flux_z(2,i,j,k,ieq)))  \
                                                + dxi*(wgt2d[2]*(flux_x(3,i+1,j,k,ieq) - flux_x(3,i,j,k,ieq)))  \
                                + dyi*(wgt2d[2]*(flux_y(3,i,j+1,k,ieq) - flux_y(3,i,j,k,ieq)))  \
                                + dzi*(wgt2d[2]*(flux_z(3,i,j,k+1,ieq) - flux_z(3,i,j,k,ieq)))  \
                                                + dxi*(wgt2d[3]*(flux_x(4,i+1,j,k,ieq) - flux_x(4,i,j,k,ieq)))  \
                                + dyi*(wgt2d[3]*(flux_y(4,i,j+1,k,ieq) - flux_y(4,i,j,k,ieq)))  \
                                + dzi*(wgt2d[3]*(flux_z(4,i,j,k+1,ieq) - flux_z(4,i,j,k,ieq))))
            end do

            do ir=2,nbasis
                do i = 1,nx

                    glflux_r(i,j,k,ieq,ir) = wgtbf_xmp(1,2,ir)*flux_x(1,i+1,j,k,ieq) + wgtbf_xmp(1,1,ir)*flux_x(1,i,j,k,ieq)  \
                                   + wgtbf_ymp(1,2,ir)*flux_y(1,i,j+1,k,ieq) + wgtbf_ymp(1,1,ir)*flux_y(1,i,j,k,ieq)  \
                                   + wgtbf_zmp(1,2,ir)*flux_z(1,i,j,k+1,ieq) + wgtbf_zmp(1,1,ir)*flux_z(1,i,j,k,ieq)  \
                                   + wgtbf_xmp(2,2,ir)*flux_x(2,i+1,j,k,ieq) + wgtbf_xmp(2,1,ir)*flux_x(2,i,j,k,ieq)  \
                                   + wgtbf_ymp(2,2,ir)*flux_y(2,i,j+1,k,ieq) + wgtbf_ymp(2,1,ir)*flux_y(2,i,j,k,ieq)  \
                                   + wgtbf_zmp(2,2,ir)*flux_z(2,i,j,k+1,ieq) + wgtbf_zmp(2,1,ir)*flux_z(2,i,j,k,ieq)  \
                                   + wgtbf_xmp(3,2,ir)*flux_x(3,i+1,j,k,ieq) + wgtbf_xmp(3,1,ir)*flux_x(3,i,j,k,ieq)  \
                                   + wgtbf_ymp(3,2,ir)*flux_y(3,i,j+1,k,ieq) + wgtbf_ymp(3,1,ir)*flux_y(3,i,j,k,ieq)  \
                                   + wgtbf_zmp(3,2,ir)*flux_z(3,i,j,k+1,ieq) + wgtbf_zmp(3,1,ir)*flux_z(3,i,j,k,ieq)  \
                                   + wgtbf_xmp(4,2,ir)*flux_x(4,i+1,j,k,ieq) + wgtbf_xmp(4,1,ir)*flux_x(4,i,j,k,ieq)  \
                                   + wgtbf_ymp(4,2,ir)*flux_y(4,i,j+1,k,ieq) + wgtbf_ymp(4,1,ir)*flux_y(4,i,j,k,ieq)  \
                                   + wgtbf_zmp(4,2,ir)*flux_z(4,i,j,k+1,ieq) + wgtbf_zmp(4,1,ir)*flux_z(4,i,j,k,ieq)  \
                                       - integral_r(i,j,k,ieq,ir)

                end do
            end do

        end do
        end do

end subroutine glflux2


#-----------------------------------------------------------!
#*******************calculate freezing speeds***************!
#-----------------------------------------------------------!
def cfcal(Qcf,cases):
    # real function cfcal(Qcf,cases)
    # implicit none
    # integer cases
    # real Qcf(nQ)
    # real Pi, Pe, P, B2, ne, cs
    # real dn,dni,vx,vy,vz,hx,hy,hz,dnei,va2,vf1,lil02,va,fac

    dn = Qcf(rh)
    dni = 1./dn
    vx = Qcf(mx)*dni
    vy = Qcf(my)*dni
    vz = Qcf(mz)*dni

    if ieos == 1:
       P = (aindex - 1.)*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
       cs = sqrt(aindex*P*dni)

    if ieos == 2:
        cs = sqrt(7.2*P_1*dn**6.2)

    select case (cases)
        case (1) # freezing speed in x direction for fluid variable
            return abs(vx) + cs

        case (2) # freezing speed in y direction for fluid variable
            return abs(vy) + cs

        case (3) # freezing speed in z direction for fluid variable
            return abs(vz) + cs
    end select

end function cfcal


#----------------------------------------------------------------------------------

def limiter(Q_ri, bf_faces):
    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    #
    # EXT LIBS:
    #  * np.sum(), np.min(), abs()
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * nx, ny, nz, nQ
    #  * rh, mx, my, mz, en
    #  * rh_floor, T_floor, aindex, npge, nbasis
    #  * bf_faces(nslim,nbastot)
    #  * Q_ri (pass-by-ref, shape (nx,ny,nz,nQ,nbasis))
    # LOCAL:
    #  * i, j, k, ieq, ir, ipge (loop vars)
    #  * Qrhmin, epsi, epsiP, thetaj, dn, dni, Pave
    #  * Qedge(npge,nQ), P(npge)
    ############################################################################

    Qedge = np.zeros((npge,nQ))
    P     = np.zeros((npge))
    epsi  = rh_floor
    epsiP = rh_floor*T_floor

    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx):
                if Q_ri[i,j,k,rh,0] < rh_floor:
                    for ir in xrange(1, nbasis):
                        Q_ri[i,j,k,rh:en+1,ir] = 0.0
                    Q_ri[i,j,k,rh,0] = rh_floor

                else:
                    for ipge in xrange(npge):
                        Qedge[ipge,rh] = np.sum(bf_faces[ipge,0:nbasis]*Q_ri[i,j,k,rh,0:nbasis])

                    Qrhmin = np.min(Qedge[:,rh])
                    if Qrhmin < epsi:
                        theta = (epsi-Q_ri[i,j,k,rh,0])/(Qrhmin-Q_ri[i,j,k,rh,0])
                        theta = 1.0 if theta > 1.0 else 0.0

                        for ir in xrange(1, nbasis):
                            Q_ri[i,j,k,rh,ir] = abs(theta)*Q_ri[i,j,k,rh,ir]

                    Pave = (aindex-1.)*(Q_ri[i,j,k,en,0] - 0.5*(Q_ri[i,j,k,mx,0]**2 + Q_ri[i,j,k,my,0]**2 + Q_ri[i,j,k,mz,0]**2)/Q_ri[i,j,k,rh,0])

                    if Pave < epsiP:
                        for ir in xrange(1, nbasis):
                            Q_ri[i,j,k,rh:en+1,ir] = 0.0

                    else:
                        theta = 1.0
                        for ipge in xrange(npge):
                            for ieq in xrange(5):
                                Qedge[ipge,ieq] = np.sum(bf_faces[ipge,0:nbasis]*Q_ri[i,j,k,ieq,0:nbasis])

                            dn = Qedge[ipge,rh]
                            dni = 1./dn
                            P[ipge] = (aindex - 1.)*(Qedge[ipge,en] - 0.5*(Qedge[ipge,mx]**2+Qedge[ipge,my]**2+Qedge[ipge,mz]**2)*dni)

                            if P[ipge] < epsiP:
                                if Pave != P[ipge]:
                                    thetaj = (Pave - epsiP)/(Pave - P[ipge])
                                    theta = np.min(theta, thetaj)

                        theta = 1.0 if theta > 1.0 else 0.0

                        for ir in xrange(1, nbasis):
                            Q_ri[i,j,k,rh:en+1,ir] = theta*Q_ri[i,j,k,rh:en+1,ir]


#----------------------------------------------------------------------------------

def limiter2(Q_ri, bf_faces):
    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    #
    # EXT LIBS:
    #  * np.sum(), np.min(), abs()
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * nx, ny, nz, nQ
    #  * rh, mx, my, mz, en
    #  * rh_floor, epsi, npge, nbasis
    #  * bf_faces(nslim,nbastot)
    #  * Q_ri (pass-by-ref, shape (nx,ny,nz,nQ,nbasis))
    # LOCAL:
    #  * i, j, k, ieq, ir, ipge (loop vars)
    #  * Qrhmin, theta
    #  * Qedge(npge,nQ)
    ############################################################################

    Qedge = np.zeros((npge,nQ))
    epsi = rh_floor

    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx):

                if Q_ri[i,j,k,rh,0] < epsi:
                    for ir in xrange(1, nbasis):
                        Q_ri[i,j,k,rh:en+1,ir] = 0.0

                    Q_ri[i,j,k,rh,0] = epsi

                else:
                    for ipge in xrange(npge):
                        Qedge[ipge,rh] = np.sum(bf_faces[ipge,0:nbasis]*Q_ri[i,j,k,rh,0:nbasis])

                    Qrhmin = np.min(Qedge[:,rh])

                    if Qrhmin < epsi:
                        theta = (epsi - Q_ri[i,j,k,rh,0])/(Qrhmin - Q_ri[i,j,k,rh,0])
                        theta = 1.0 if theta > 1.0 else 0.0

                        for ir in xrange(1, nbasis):
                            Q_ri[i,j,k,rh,ir] = abs(theta)*Q_ri[i,j,k,rh,ir]

#----------------------------------------------------------------------------------------------

def prepare_exchange(Q_ri,
                     bfvals_xm, bfvals_xp,
                     bfvals_ym, bfvals_yp,
                     bfvals_zm, bfvals_zp):
    ############################################################################
    # Necessary functions
    #------------------------------------
    # PERSEUS:
    #
    # EXT LIBS:
    #  * np.sum()
    #
    # Necessary parameters and variables
    #------------------------------------
    # GLOBAL:
    #  * nx, ny, nz, nQ
    #  * rh, mx, my, mz, en
    #  * rh_floor, epsi, npge, nface
    #  * bf_faces(nslim,nbastot)
    #  * Q_ri (pass-by-ref, shape (nx,ny,nz,nQ,nbasis))
    # LOCAL:
    #  * i, j, k, ieq, ipnt (loop vars)
    #  * bfvals_zm, bfvals_zp  (both shape (nface,nbastot))
    #  * bfvals_ym, bfvals_yp  (both shape (nface,nbastot))
    #  * bfvals_xm, bfvals_xp  (both shape (nface,nbastot))
    ############################################################################
    # integer ieq, i, j, k, ipnt
    # real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r

    for ieq in xrange(nQ):
        for j in xrange(ny):
            for i in xrange(nx):
                for ipnt in xrange(nface):
                    Qzlow_int[i,j,ipnt,ieq] = np.sum(bfvals_zm[ipnt,0:nbasis]*Q_ri[i,j,0,ieq,0:nbasis])
                    Qzhigh_int[i,j,ipnt,ieq] = np.sum(bfvals_zp[ipnt,0:nbasis]*Q_ri[i,j,nz,ieq,0:nbasis])

    for ieq in xrange(nQ):
        for k in xrange(nz):
            for i in xrange(nx):
                for ipnt in xrange(nface):
                    Qylow_int[i,k,ipnt,ieq] = np.sum(bfvals_ym[ipnt,0:nbasis]*Q_ri[i,0,k,ieq,0:nbasis])
                    Qyhigh_int[i,k,ipnt,ieq] = np.sum(bfvals_yp[ipnt,0:nbasis]*Q_ri[i,ny,k,ieq,0:nbasis])

    for ieq in xrange(nQ):
        for k in xrange(nz):
            for j in xrange(ny):
                for ipnt in xrange(nface):
                    Qxlow_int[j,k,ipnt,ieq] = np.sum(bfvals_xm[ipnt,0:nbasis]*Q_ri[0,j,k,ieq,0:nbasis])
                    Qxhigh_int[j,k,ipnt,ieq] = np.sum(bfvals_xp[ipnt,0:nbasis]*Q_ri[nx,j,k,ieq,0:nbasis])

    call exchange_flux  # FIXME

end subroutine

#----------------------------------------------------------------------------------------------

# FIXME
def exchange_flux():
    # integer mpi_size

    call MPI_BARRIER(cartcomm,ierr)

    mpi_size=ny*nz*nface*nQ

    if nbrs(EAST) != MPI_PROC_NULL:
        call MPI_ISend(Qxhigh_int,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(1),ierr)

    if nbrs(WEST) != MPI_PROC_NULL:
        call MPI_IRecv(Qxlow_ext,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(2),ierr)

    if nbrs(EAST) != MPI_PROC_NULL:
        call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(Qxhigh_ext,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(3),ierr)

    if nbrs(WEST) != MPI_PROC_NULL:
        call MPI_Wait(reqs(2),stats(:,2),ierr)
        call MPI_ISend(Qxlow_int,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(4),ierr)

    if nbrs(EAST) != MPI_PROC_NULL:
        call MPI_Wait(reqs(3),stats(:,3),ierr)

    if nbrs(WEST) != MPI_PROC_NULL:
        call MPI_Wait(reqs(4),stats(:,4),ierr)

    if mpi_P == 1 and xlbc == 0:
        Qxlow_ext = Qxlow_int

    if mpi_P == mpi_nx and xhbc == 0:
        Qxhigh_ext = Qxhigh_int

    mpi_size=nface*nx*nz*nQ

    if nbrs(NORTH) != MPI_PROC_NULL:
        call MPI_ISend(Qyhigh_int,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(1),ierr)

    if nbrs(SOUTH) != MPI_PROC_NULL:
        call MPI_IRecv(Qylow_ext,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(2),ierr)

    if nbrs(NORTH) != MPI_PROC_NULL:
        call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(Qyhigh_ext,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(3),ierr)

    if nbrs(SOUTH) != MPI_PROC_NULL:
        call MPI_Wait(reqs(2),stats(:,2),ierr)
        call MPI_ISend(Qylow_int,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(4),ierr)

    if nbrs(NORTH) != MPI_PROC_NULL:
        call MPI_Wait(reqs(3),stats(:,3),ierr)

    if nbrs(SOUTH) != MPI_PROC_NULL:
        call MPI_Wait(reqs(4),stats(:,4),ierr)


    if mpi_Q == 1 and ylbc == 0:
        Qylow_ext = Qylow_int

    if mpi_Q == mpi_ny and yhbc == 0:
        Qyhigh_ext = Qyhigh_int


    mpi_size = nface*nx*ny*nQ

    if nbrs(UP) != MPI_PROC_NULL:
        call MPI_ISend(Qzhigh_int,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(1),ierr)
    if nbrs(DOWN) != MPI_PROC_NULL:
        call MPI_IRecv(Qzlow_ext,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(2),ierr)
    if nbrs(UP) != MPI_PROC_NULL:
        call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(Qzhigh_ext,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(3),ierr)
    if nbrs(DOWN) != MPI_PROC_NULL:
        call MPI_Wait(reqs(2),stats(:,2),ierr)
        call MPI_ISend(Qzlow_int,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(4),ierr)
    if nbrs(UP) != MPI_PROC_NULL:
        call MPI_Wait(reqs(3),stats(:,3),ierr)

    if nbrs(DOWN) != MPI_PROC_NULL:
        call MPI_Wait(reqs(4),stats(:,4),ierr)


    if mpi_R == 1 and zlbc == 0:
        Qzlow_ext = Qzlow_int

    if mpi_R == mpi_nz and zhbc == 0:
        Qzhigh_ext = Qzhigh_int

end subroutine


#-----------------------------------------------------------

def xc(i):
    return (i - 0.5)*dx

#-----------------------------------------------------------

def yc(j):
    return (j - 0.5)*dy

#-----------------------------------------------------------

def zc(k):
    return (k - 0.5)*dz

#-----------------------------------------------------------

def rz(i,j):
    return sqrt(yc(j)**2)

#-----------------------------------------------------------

def r(i,j):
    return (xc(i)**2 + yc(j)**2)**0.5

#-----------------------------------------------------------

def theta(i,j):
    return atan2(yc(j),xc(i))

#-----------------------------------------------------------

def xvtk(i):
    return loc_lxd + (i - 0.5)*dxvtk

#-----------------------------------------------------------

def yvtk(j):
    return loc_lyd + (j - 0.5)*dyvtk

#-----------------------------------------------------------

def zvtk(k):
    return loc_lzd + (k - 0.5)*dzvtk

#-----------------------------------------------------------

def rvtk(i,j):
    return sqrt(xvtk(i)**2 + yvtk(j)**2)

#-----------------------------------------------------------

def thetavtk(i,j):
    return atan2(yvtk(j),xvtk(i))

#------------More realistic COBRA  current driver -----------


def E_z(t):
    # real t
    if t <= tr:
        return Ez0*sin(0.5*pi*t/tr)*cos(0.5*pi*t/tr)
    if t >= tr:
        return 0.0

#------------More realistic COBRA  current driver -----------


def Icur(t):
    # real t
    if t <= tr:
        return Ipeak*sin(0.5*pi*t/tr)**2
    if t >= tr:
        return Ipeak

    # if t >= tr and t <= 2*tr:
    #     return Ipeak
    # if t >= 2*tr and t <= 3*tr:
    #     return Ipeak*sin(0.5*pi*(t-tr)/tr)
    # if t >= 3*tr:
    #     return 0.

#-----------------------------------------------------------

def output_vtk0(Qin,nout,iam):
    # real Qin(nx,ny,nz,nQ,nbasis)
    # integer nout
    # integer(I4P) , parameter :: nb=1*ngu,sprd=0*1
###   integer(I4P), parameter:: nnx=nx-2*nb,nny=ny-2*nb,nnz=nz-2*nb
    # real(R4P), dimension(nnx+1):: x_xml_rect
    # real(R4P), dimension(nny+1):: y_xml_rect
    # real(R4P), dimension(nnz+1):: z_xml_rect
    # real(R4P), dimension(nnx*nny*nnz):: var_xml_val_x
    # real(R4P), dimension(nnx*nny*nnz):: var_xml_val_y
    # real(R4P), dimension(nnx*nny*nnz):: var_xml_val_z
    # real P, vx, vy, vz,  dni
    # integer(I4P):: E_IO,i,j,k,l,num,iam
    # character (50) :: out_name
    # character (4) :: tname
    # character (5) :: tname1
    # character (4) :: pname
    # character (5) :: pname1


    num=nout+10000

    # write(tname1,'(i5)')num  # FIX THIS
    tname=tname1(2:5)
    tname = trim(tname)
    tname = adjustr(tname)

    num=iam+10000

    # write(pname1,'(i5)')num  # FIX THIS
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)
#   out_name='/data/data5/perseus_p'//pname//'_t'//tname//'.vtr'
    out_name='data/perseus_p'//pname//'_t'//tname//'.vtr'
#        print *, out_name
    out_name = trim(out_name)
    out_name = adjustr(out_name)

    E_IO = VTK_INI_XML(output_format = 'BINARY',              &
                       filename      = out_name, &
                       mesh_topology = 'RectilinearGrid',     &
                       nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1)


    do i=1+nb,nx-nb+1
        x_xml_rect(i)=(xc(i) - 0.5/dxi)*1e-3
    enddo
    do i=1+nb,ny-nb+1
        y_xml_rect(i)=(yc(i) - 0.5/dyi)*1e-3
    enddo
    do i=1+nb,nz-nb+1
        z_xml_rect(i)=(zc(i) - 0.5/dzi)*1e-3
    enddo


    E_IO = VTK_GEO_XML(nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1, &
                         X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)

    E_IO = VTK_DAT_XML(var_location     = 'cell', &
                       var_block_action = 'OPEN')

    do i=1+nb,nx-nb
        do j=1+nb,ny-nb
            do k=1+nb,nz-nb
                l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
                var_xml_val_x(l)=log(Qin(i,j,k,rh,1)*n0)/log(10.)

            enddo
        enddo
    enddo
    E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                       varname = 'Log Density',                    &
                       var     = var_xml_val_x)

    do i=1+nb,nx-nb
        do j=1+nb,ny-nb
            do k=1+nb,nz-nb
                l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
                var_xml_val_x(l)=Qin(i,j,k,rh,1)*n0
            enddo
        enddo
    enddo
    E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                       varname = 'Density',                    &
                       var     = var_xml_val_x)

    do i=1+nb,nx-nb
        do j=1+nb,ny-nb
            do k=1+nb,nz-nb
                l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
                dni = 1./Qin(i,j,k,rh,1)
                vx = Qin(i,j,k,mx,1)*dni
                vy = Qin(i,j,k,my,1)*dni
                vz = Qin(i,j,k,mz,1)*dni
                P = (aindex - 1)*(Qin(i,j,k,en,1) - 0.5*Qin(i,j,k,rh,1)*(vx**2 + vy**2 + vz**2))*dni
                var_xml_val_x(l)=P*te0
            enddo
        enddo
    enddo
    E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                       varname = 'Temperature',                    &
                       var     = var_xml_val_x)

    do i=1+nb,nx-nb
        do j=1+nb,ny-nb
            do k=1+nb,nz-nb
                l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
                dni = 1./Qin(i,j,k,rh,1)
                vx = Qin(i,j,k,mx,1)*dni
                vy = Qin(i,j,k,my,1)*dni
                vz = Qin(i,j,k,mz,1)*dni
                if ieos == 1:
                    P = (aindex - 1)*(Qin(i,j,k,en,1) - 0.5*Qin(i,j,k,rh,1)*(vx**2 + vy**2 + vz**2))

                var_xml_val_x(l)=P
            enddo
        enddo
    enddo
    E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                       varname = 'pressure',                    &
                       var     = var_xml_val_x)

    do i=1+nb,nx-nb
        do j=1+nb,ny-nb
            do k=1+nb,nz-nb
                l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
                var_xml_val_x(l)=v0*Qin(i,j,k,mx,1)/Qin(i,j,k,rh,1)
                var_xml_val_y(l)=v0*Qin(i,j,k,my,1)/Qin(i,j,k,rh,1)
                var_xml_val_z(l)=v0*Qin(i,j,k,mz,1)/Qin(i,j,k,rh,1)
            enddo
        enddo
    enddo
    E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                       varname = 'Velocity',                    &
                       varX     = var_xml_val_x, varY     = var_xml_val_y, varZ     = var_xml_val_Z )

    do i=1+nb,nx-nb
        do j=1+nb,ny-nb
            do k=1+nb,nz-nb
                l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
                    var_xml_val_x(l)= -2*(Q_r2(i,j,k,my,2)-Q_r2(i,j,k,mx,3))/Q_r2(i,j,k,rh,1)*dxi/t0
            enddo
        enddo
    enddo
    E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                       varname = 'Vorticity',                    &
                       var     = var_xml_val_x)


    E_IO = VTK_DAT_XML(var_location     = 'cell', &
                       var_block_action = 'Close')
    E_IO = VTK_GEO_XML()
    E_IO = VTK_END_XML()

end subroutine output_vtk0

#-----------------------------------------------------------

def output_vtk(Qin,nout,iam):
    # real Qin(nx,ny,nz,nQ,nbasis)
    # integer nout
    # integer(I4P) , parameter :: nb=1*ngu,sprd=0*1
    #
    # real(R4P), dimension(nnx+1) :: x_xml_rect
    # real(R4P), dimension(nny+1) :: y_xml_rect
    # real(R4P), dimension(nnz+1) :: z_xml_rect
    # real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_x
    # real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_y
    # real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_z
    # real(R4P), dimension(nnx,nny,nnz,nQ) :: qvtk
    # real(R4P), dimension(nnx,nny,nnz) :: qvtk_dxvy,qvtk_dyvx
    # real P, vx, vy, vz,  dni, dxrh,dyrh,dxmy,dymx
    # integer(I4P):: E_IO,i,j,k,l,num,iam,igrid,ir,jr,kr,ib,jb,kb,ieq
    # character (50) :: out_name
    # character (4) :: tname
    # character (5) :: tname1
    # character (4) :: pname
    # character (5) :: pname1


    num=nout+10000

    # write(tname1,'(i5)')num  # FIX THIS
    tname=tname1(2:5)
    tname = trim(tname)
    tname = adjustr(tname)

    num=iam+10000

    # write(pname1,'(i5)')num  # FIX THIS
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)
    out_name='/data/data/perseus_p'//pname//'_t'//tname//'.vtr'
    # out_name='data/perseus_p'//pname//'_t'//tname//'.vtr'
    # print *, out_name
    out_name = trim(out_name)
    out_name = adjustr(out_name)

    E_IO = VTK_INI_XML(output_format = 'BINARY',              &
                       filename      = out_name, &
                       mesh_topology = 'RectilinearGrid',     &
                       nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1)




        do i=1+nb,nnx-nb+1
            x_xml_rect(i-nb)=(xvtk(i) - 0.5*dxvtk)*1e-3
        enddo
        do j=1+nb,nny-nb+1
            y_xml_rect(j-nb)=(yvtk(j) - 0.5*dyvtk)*1e-3
        enddo
        do k=1+nb,nnz-nb+1
            z_xml_rect(k-nb)=(zvtk(k) - 0.5*dzvtk)*1e-3
        enddo


        E_IO = VTK_GEO_XML(nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1, &
                             X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)

        E_IO = VTK_DAT_XML(var_location     = 'cell', &
                           var_block_action = 'OPEN')


        do i=1+nb,nnx-nb
            ir = int(i/nvtk) + 1
            ib = i - nvtk*int(i/nvtk)
            if ib == 0:
                ib = nvtk
                ir = ir - 1

            do j=1+nb,nny-nb
                jr = int(j/nvtk) + 1
                jb = j - nvtk*int(j/nvtk)
                if jb == 0:
                    jb = nvtk
                    jr = jr - 1

                do k=1+nb,nnz-nb
                    kr = int(k/nvtk) + 1
                    kb = k - nvtk*int(k/nvtk)
                    if kb == 0:
                        kb = nvtk
                        kr = kr - 1

                    igrid = nvtk2*(ib-1) + nvtk*(jb-1) + kb
                    do ieq=1,nQ
                            qvtk(i,j,k,ieq) = sum(bfvtk(igrid,1:nbasis)*Qin(ir,jr,kr,ieq,1:nbasis))
                    end do

                   dxrh = sum(bfvtk_dx(igrid,1:nbasis)*Qin(ir,jr,kr,rh,1:nbasis))
                   dyrh = sum(bfvtk_dy(igrid,1:nbasis)*Qin(ir,jr,kr,rh,1:nbasis))
                   dxmy = sum(bfvtk_dx(igrid,1:nbasis)*Qin(ir,jr,kr,my,1:nbasis))
                   dymx = sum(bfvtk_dy(igrid,1:nbasis)*Qin(ir,jr,kr,mx,1:nbasis))
                   qvtk_dxvy(i,j,k) = (qvtk(i,j,k,rh)*dxmy - qvtk(i,j,k,my)*dxrh)/qvtk(i,j,k,rh)**2
                   qvtk_dyvx(i,j,k) = (qvtk(i,j,k,rh)*dymx - qvtk(i,j,k,mx)*dyrh)/qvtk(i,j,k,rh)**2
                end do
            end do
        end do

        do i=1+nb,nnx-nb
            do j=1+nb,nny-nb
                do k=1+nb,nnz-nb
                    l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
                    var_xml_val_x(l)=log(qvtk(i,j,k,rh)*n0)/log(10.)
                enddo
            enddo
        enddo
        E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                           varname = 'Log Density',                    &
                           var     = var_xml_val_x)

        do i=1+nb,nnx-nb
            do j=1+nb,nny-nb
                do k=1+nb,nnz-nb
                    l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
                    var_xml_val_x(l)=qvtk(i,j,k,rh)*n0
                enddo
            enddo
        enddo
        E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                           varname = 'Density',                    &
                           var     = var_xml_val_x)


        do i=1+nb,nnx-nb
            do j=1+nb,nny-nb
                do k=1+nb,nnz-nb
                    l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
                    dni = 1./qvtk(i,j,k,rh)
                    vx = qvtk(i,j,k,mx)*dni
                    vy = qvtk(i,j,k,my)*dni
                    vz = qvtk(i,j,k,mz)*dni
                    P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))*dni
                    var_xml_val_x(l)=P*te0
                enddo
            enddo
        enddo
        E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                           varname = 'Ideal EOS Temperature',                    &
                           var     = var_xml_val_x)

        do i=1+nb,nnx-nb
            do j=1+nb,nny-nb
                do k=1+nb,nnz-nb
                    l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
                    dni = 1./qvtk(i,j,k,rh)
                    vx = qvtk(i,j,k,mx)*dni
                    vy = qvtk(i,j,k,my)*dni
                    vz = qvtk(i,j,k,mz)*dni
                    if ieos == 1:
                        P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))
                    if ieos == 2:
                        P = P_1*(qvtk(i,j,k,rh)**7.2 - 1.) + P_base
                        if P < P_floor:
                            P = P_floor
                    var_xml_val_x(l)=P*P0
                enddo
            enddo
        enddo
        E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                           varname = 'pressure',                    &
                           var     = var_xml_val_x)

        do i=1+nb,nnx-nb
            do j=1+nb,nny-nb
                do k=1+nb,nnz-nb
                    l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
                    var_xml_val_x(l)=v0*qvtk(i,j,k,mx)/qvtk(i,j,k,rh)
                    var_xml_val_y(l)=v0*qvtk(i,j,k,my)/qvtk(i,j,k,rh)
                    var_xml_val_z(l)=v0*qvtk(i,j,k,mz)/qvtk(i,j,k,rh)
                enddo
            enddo
        enddo
        E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                           varname = 'Velocity',                    &
                           varX     = var_xml_val_x, varY     = var_xml_val_y, varZ     = var_xml_val_Z )

        do i=1+nb,nnx-nb
            do j=1+nb,nny-nb
                do k=1+nb,nnz-nb
                    l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
                        var_xml_val_x(l)= -2.*(qvtk_dxvy(i,j,k)-qvtk_dyvx(i,j,k))*dxi/t0
                enddo
            enddo
        enddo
        E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                           varname = 'Vorticity',                    &
                           var     = var_xml_val_x)


        E_IO = VTK_DAT_XML(var_location     = 'cell', &
                           var_block_action = 'Close')
        E_IO = VTK_GEO_XML()
        E_IO = VTK_END_XML()

        end subroutine output_vtk

#--------------------------------------------------------------------------------

def writeQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow):
    # real :: Qin(nx,ny,nz,nQ,nbasis),tnow,dtnow
    # integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow,nump,numd,qq,k,j,i,ir
    # character (4) :: fprefix,pname,dname
    # character (5) :: pname1,dname1
    # character (30) :: fname2

    nump = iam + 10000

    # write(pname1,'(i5)')nump  # FIX THIS
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)

    numd = iddump + 10000

    # write(dname1,'(i5)') numd  # FIX THIS
    dname = dname1(2:5)
    dname = trim(dname)
    dname = adjustr(dname)

    fname2 = 'data/'//fprefix//'_p'//pname//'_d'//dname//'.dat'
#   print *,'fname2 ',fname2

    # open(unit=3,file=fname2)  # FIX THIS
    with open(unit=3,file=fname2):
        pass  # FIX THIS

#   open(unit = 10, file = 'data/perseus_t'//dname//'_p'//pname//'.bin',form = 'unformatted',access = 'stream')

    for ir in xrange(nbasis):
        do qq=1,nQ
            do k=1,nz
                do j=1,ny
                    # write(3,*) (Qin(i,j,k,qq,ir),i=1,nx)  # FIX THIS
                enddo
            enddo
        enddo

    # write(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow  # FIX THIS

    # close(3)  # FIX THIS

    end subroutine writeQ

#------------------------------------------------------------------------------

def readQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow):
    # real :: Qin(nx,ny,nz,nQ,nbasis),tnow,dtnow
    # integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow,nump,numd,qq,k,j,i,ir
    # character (4) :: fprefix,pname,dname
    # character (5) :: pname1,dname1
    # character (30) :: fname2

    nump = irank + 10000

    # write(pname1,'(i5)')nump  # FIX THIS
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)

    numd = iddump + 10000

    # write(dname1,'(i5)')numd  # FIX THIS
    dname=dname1(2:5)
    dname = trim(dname)
    dname = adjustr(dname)

    fname2 = 'data/'//fpre//'_p'//pname//'_d'//dname//'.dat'

    # open(unit=3,file=fname2,action='read')  # FIX THIS

    for ir in xrange(nbasis):
        do qq=1,nQ
            do k=1,nz
                do j=1,ny
                    # read(3,*) (Qin(i,j,k,qq,ir),i=1,nx)  # FIX THIS
                enddo
            enddo
        enddo

    # read(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow  # FIX THIS

    # close(3)  # FIX THIS

end subroutine readQ

#--------------------------------------------------------------------------------

def set_bfvals_3D():
        # Defines local basis function values and weights for 1, 2, or 3-point Gaussian quadrature.
        # Basis functions are evaluated in cell interior and on cell faces.

    call set_vtk_vals_3D()    # Define basis function values on a 3D grid INTERNAL to the cell.
    call set_internal_vals_3D()     # Define basis function values at quadrature points INTERNAL to cell.
    call set_face_vals_3D()   # Define local basis function values at quadrature points on a cell face.
    call set_weights_3D()     # Define weights for integral approximation using Gaussian quadrature.

end subroutine set_bfvals_3D

#----------------------------------------------------------
def set_vtk_vals_3D():
    # Define basis function values on a rectangular 3D grid INTERNAL to the cell.
    # For use in output_vtk().

    # integer ixyz,i,igrid

    dxvtk = 1./(nvtk*dxi)
    dyvtk = 1./(nvtk*dyi)
    dzvtk = 1./(nvtk*dzi)

    if nvtk == 1:
        xgrid(1) = 0.
    if nvtk == 2:
        xgrid(1) = -0.5
        xgrid(2) = 0.5
    if nvtk == 3:
        xgrid(1) = -2./3.
        xgrid(2) = 0.
        xgrid(3) = 2./3.
    if nvtk == 4:
        xgrid(1) = -3./4.
        xgrid(2) = -1./4.
        xgrid(3) = 1./4.
        xgrid(4) = 3./4.


        bfvtk(1,1) = 1.     # basis function = 1
        bfvtk(1,kx) = 0.        # basis function = x
        bfvtk(1,ky) = 0.        # basis function = y
        bfvtk(1,kz) = 0.        # basis function = z


        bfvtk(1:nvtk3,1) = 1.       # basis function = 1
        do i=1,nvtk
            bfvtk((i-1)*nvtk+1:i*nvtk,ky) = xgrid(i)     # basis function = y
            bfvtk((i-1)*nvtk+1:i*nvtk,kz) = xgrid(1:nvtk)        # basis function = z
        end do
        do i=1,nvtk
            bfvtk((i-1)*nvtk2+1:i*nvtk2,kx) = xgrid(i)       # basis function = x
            bfvtk((i-1)*nvtk2+1:i*nvtk2,ky) = bfvtk(1:nvtk2,ky)      # basis function = y
            bfvtk((i-1)*nvtk2+1:i*nvtk2,kz) = bfvtk(1:nvtk2,kz)      # basis function = z
        end do

    do i=0,2
        bfvtk(1:nvtk3,kxx+i) = 1.5*bfvtk(1:nvtk3,kx+i)**2 - 0.5    # basis function = P_2(s)
        # bfvtk(1:nvtk3,kxxx+i) = 2.5*bfvtk(1:nvtk3,kx+i)**3 - 1.5*bfvtk(1:nvtk3,kx+i)   # basis function = P3(s)
    end do

    bfvtk(1:nvtk3,kyz) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kz)        # basis function = yz
    bfvtk(1:nvtk3,kzx) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kx)        # basis function = zx
    bfvtk(1:nvtk3,kxy) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,ky)        # basis function = xy
    bfvtk(1:nvtk3,kxyz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kz)     # basis function = xyz

    bfvtk(1:nvtk3,kyyz) = bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kz)     # basis function = P2(y)z
    bfvtk(1:nvtk3,kyzz) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kzz)     # basis function = P2(z)y
    bfvtk(1:nvtk3,kzzx) = bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kx)     # basis function = P2(z)x
    bfvtk(1:nvtk3,kzxx) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kxx)     # basis function = P2(x)z
    bfvtk(1:nvtk3,kxxy) = bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,ky)     # basis function = P2(x)y
    bfvtk(1:nvtk3,kxyy) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)     # basis function = P2(y)x
    bfvtk(1:nvtk3,kyyzz) = bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)     # basis function = P_2(y)P_2(z)
    bfvtk(1:nvtk3,kzzxx) = bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kxx)     # basis function = P_2(z)P_2(x)
    bfvtk(1:nvtk3,kxxyy) = bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,kyy)     # basis function = P_2(x)P_2(y)
    bfvtk(1:nvtk3,kyzxx) = bfvtk(1:nvtk3,kyz)*bfvtk(1:nvtk3,kxx)     # basis function = yz P_2(x)
    bfvtk(1:nvtk3,kzxyy) = bfvtk(1:nvtk3,kzx)*bfvtk(1:nvtk3,kyy)     # basis function = zx P_2(y)
    bfvtk(1:nvtk3,kxyzz) = bfvtk(1:nvtk3,kxy)*bfvtk(1:nvtk3,kzz)     # basis function = xy P_2(z)
    bfvtk(1:nvtk3,kxyyzz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)   # basis function = x P_2(y)P_2(z)
    bfvtk(1:nvtk3,kyzzxx) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kxx)   # basis function = y P_2(z)P_2(x)
    bfvtk(1:nvtk3,kzxxyy) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,kyy)   # basis function = z P_2(x)P_2(y)
    bfvtk(1:nvtk3,kxxyyzz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)  # basis function = P_2(x)P_2(y)P_2(z)

    do igrid=1,nvtk3
        bfvtk_dx(igrid,1) = 0.
        bfvtk_dx(igrid,kx) = 1.
        bfvtk_dx(igrid,ky:kz) = 0.
        bfvtk_dx(igrid,kyz) = 0.
        bfvtk_dx(igrid,kzx) = bfvtk(igrid,kz)
        bfvtk_dx(igrid,kxy) = bfvtk(igrid,ky)
        bfvtk_dx(igrid,kxyz) = bfvtk(igrid,kyz)

        bfvtk_dx(igrid,kxx) = 3.*bfvtk(igrid,kx)
        bfvtk_dx(igrid,kyy:kzz) = 0.
        bfvtk_dx(igrid,kyyz) = 0.
        bfvtk_dx(igrid,kyzz) = 0.
        bfvtk_dx(igrid,kzzx) = bfvtk(igrid,kzz)
        bfvtk_dx(igrid,kzxx) = bfvtk(igrid,kz)*3.*bfvtk(igrid,kx)
        bfvtk_dx(igrid,kxxy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,ky)
        bfvtk_dx(igrid,kxyy) = bfvtk(igrid,kyy)
        bfvtk_dx(igrid,kyyzz) = 0.
        bfvtk_dx(igrid,kzzxx) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kzz)
        bfvtk_dx(igrid,kxxyy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyy)
        bfvtk_dx(igrid,kyzxx) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyz)
        bfvtk_dx(igrid,kzxyy) = bfvtk(igrid,kz)*bfvtk(igrid,kyy)
        bfvtk_dx(igrid,kxyzz) = bfvtk(igrid,ky)*bfvtk(igrid,kzz)
        bfvtk_dx(igrid,kxyyzz) = bfvtk(igrid,kyyzz)
        bfvtk_dx(igrid,kyzzxx) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyzz)
        bfvtk_dx(igrid,kzxxyy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyyz)
        bfvtk_dx(igrid,kxxyyzz) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyy)*bfvtk(igrid,kzz)

        bfvtk_dy(igrid,1) = 0.
        bfvtk_dy(igrid,ky) = 1.
        bfvtk_dy(igrid,kx) = 0.
        bfvtk_dy(igrid,kz) = 0.
        bfvtk_dy(igrid,kyz) = bfvtk(igrid,kz)
        bfvtk_dy(igrid,kzx) = 0.
        bfvtk_dy(igrid,kxy) = bfvtk(igrid,kx)
        bfvtk_dy(igrid,kxyz) = bfvtk(igrid,kzx)

        bfvtk_dy(igrid,kyy) = 3.*bfvtk(igrid,ky)
        bfvtk_dy(igrid,kxx) = 0.
        bfvtk_dy(igrid,kzz) = 0.
        bfvtk_dy(igrid,kyzz) = bfvtk(igrid,kzz)
        bfvtk_dy(igrid,kyyz) = bfvtk(igrid,kz)*3.*bfvtk(igrid,ky)
        bfvtk_dy(igrid,kzzx) = 0.
        bfvtk_dy(igrid,kzxx) = 0.
        bfvtk_dy(igrid,kxxy) = bfvtk(igrid,kxx)
        bfvtk_dy(igrid,kxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kx)
        bfvtk_dy(igrid,kyyzz) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzz)
        bfvtk_dy(igrid,kzzxx) = 0.
        bfvtk_dy(igrid,kxxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kxx)
        bfvtk_dy(igrid,kyzxx) = bfvtk(igrid,kz)*bfvtk(igrid,kxx)
        bfvtk_dy(igrid,kzxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzx)
        bfvtk_dy(igrid,kxyzz) = bfvtk(igrid,kx)*bfvtk(igrid,kzz)
        bfvtk_dy(igrid,kxyyzz) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzzx)
        bfvtk_dy(igrid,kyzzxx) = bfvtk(igrid,kzzxx)
        bfvtk_dy(igrid,kzxxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzxx)
        bfvtk_dy(igrid,kxxyyzz) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzz)*bfvtk(igrid,kxx)


        bfvtk_dz(igrid,1) = 0.
        bfvtk_dz(igrid,kz) = 1.
        bfvtk_dz(igrid,kx) = 0.
        bfvtk_dz(igrid,ky) = 0.
        bfvtk_dz(igrid,kyz) = bfvtk(igrid,ky)
        bfvtk_dz(igrid,kzx) = bfvtk(igrid,kx)
        bfvtk_dz(igrid,kxy) = 0.
        bfvtk_dz(igrid,kxyz) = bfvtk(igrid,kxy)

        bfvtk_dz(igrid,kzz) = 3.*bfvtk(igrid,kz)
        bfvtk_dz(igrid,kxx) = 0.
        bfvtk_dz(igrid,kyy) = 0.

        bfvtk_dz(igrid,kyzz) = bfvtk(igrid,ky)*3.*bfvtk(igrid,kz)
        bfvtk_dz(igrid,kyyz) = bfvtk(igrid,kyy)
        bfvtk_dz(igrid,kzzx) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kx)
        bfvtk_dz(igrid,kzxx) = bfvtk(igrid,kxx)
        bfvtk_dz(igrid,kxxy) = 0.
        bfvtk_dz(igrid,kxyy) = 0.
        bfvtk_dz(igrid,kyyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kyy)
        bfvtk_dz(igrid,kzzxx) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxx)
        bfvtk_dz(igrid,kxxyy) = 0.
        bfvtk_dz(igrid,kyzxx) = bfvtk(igrid,ky)*bfvtk(igrid,kxx)
        bfvtk_dz(igrid,kzxyy) = bfvtk(igrid,kx)*bfvtk(igrid,kyy)
        bfvtk_dz(igrid,kxyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxy)
        bfvtk_dz(igrid,kxyyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxyy)
        bfvtk_dz(igrid,kyzzxx) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxxy)
        bfvtk_dz(igrid,kzxxyy) = bfvtk(igrid,kxxyy)
        bfvtk_dz(igrid,kxxyyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxx)*bfvtk(igrid,kyy)


        # bfvtk_dx(igrid,kxxx) = 7.5*bfvtk(igrid,kx)**2 - 1.5
        # bfvtk_dx(igrid,kyyy:kzzz) = 0.
        # bfvtk_dy(igrid,kyyy) = 7.5*bfvtk(igrid,ky)**2 - 1.5
        # bfvtk_dy(igrid,kxxx) = 0.
        # bfvtk_dy(igrid,kzzz) = 0.
        # bfvtk_dz(igrid,kzzz) = 7.5*bfvtk(igrid,kz)**2 - 1.5
        # bfvtk_dz(igrid,kxxx) = 0.
        # bfvtk_dz(igrid,kyyy) = 0.
    end do


end subroutine set_vtk_vals_3D

#----------------------------------------------------------

#----------------------------------------------------------
def set_internal_vals_3D():
    # Define basis function values at quadrature points INTERNAL to cell.

    # integer ixyz, i
    # real c15d5,c1dsq3,xq4p,xq4m

    c15d5 = sqrt(15.)/5.
    c1dsq3 = 1./sqrt(3.)

    xq4p = sqrt(525. + 70.*sqrt(30.))/35.
    xq4m = sqrt(525. - 70.*sqrt(30.))/35.

    if iquad == 2:
        xquad[0] = -c1dsq3
        xquad[2] = c1dsq3
    if iquad == 3:
        xquad[0] = -c15d5
        xquad[1] = 0.
        xquad[2] = c15d5
    if iquad == 4:
        xquad[0] = -xq4p
        xquad[1] = -xq4m
        xquad[2] = xq4m
        xquad[3] = xq4p

    if iquad == 1:          # 2-point Gaussian quadrature
        bfvals_int[0,0] = 1.        # basis function = 1
        bfvals_int[0,kx] = 0.       # basis function = x
        bfvals_int[0,ky] = 0.       # basis function = y
        bfvals_int[0,kz] = 0.       # basis function = z


    if iquad > 1:
        bfvals_int[0:npg,0] = 1.        # basis function = 1
        for i in xrange(nedge):
            bfvals_int[i*nedge:i*nedge,ky] = xquad[i]      # basis function = y
            bfvals_int[i*nedge:i*nedge,kz] = xquad[0:nedge]        # basis function = z
        for i in xrange(nedge):
            bfvals_int[i*nface:i*nface,kx] = xquad[i]      # basis function = x
            bfvals_int[i*nface:i*nface,ky] = bfvals_int[0:nface,ky]        # basis function = y
            bfvals_int[i*nface:i*nface,kz] = bfvals_int[0:nface,kz]        # basis function = z

    for xrange(3):
        bfvals_int[0:npg,kxx+i] = 1.5*bfvals_int[0:npg,kx+i]**2 - 0.5    # basis function = P_2(s)
        # bfvals_int(1:npg,kxxx+i) = 2.5*bfvals_int(1:npg,kx+i)**3 - 1.5*bfvals_int(1:npg,kx+i)   # basis function = P_3(s)

    bfvals_int[0:npg,kyz] = bfvals_int[0:npg,ky]*bfvals_int[0:npg,kz]       # basis function = yz
    bfvals_int[0:npg,kzx] = bfvals_int[0:npg,kz]*bfvals_int[0:npg,kx]       # basis function = zx
    bfvals_int[0:npg,kxy] = bfvals_int[0:npg,kx]*bfvals_int[0:npg,ky]       # basis function = xy
    bfvals_int[0:npg,kxyz] = bfvals_int[0:npg,kx]*bfvals_int[0:npg,ky]*bfvals_int[0:npg,kz]     # basis function = xyz


    bfvals_int[0:npg,kyyz] = bfvals_int[0:npg,kyy]*bfvals_int[0:npg,kz]     # basis function = P_2(y)z
    bfvals_int[0:npg,kyzz] = bfvals_int[0:npg,ky]*bfvals_int[0:npg,kzz]     # basis function = P_2(z)y
    bfvals_int[0:npg,kzzx] = bfvals_int[0:npg,kzz]*bfvals_int[0:npg,kx]     # basis function = P_2(z)x
    bfvals_int[0:npg,kzxx] = bfvals_int[0:npg,kz]*bfvals_int[0:npg,kxx]     # basis function = P_2(x)z
    bfvals_int[0:npg,kxxy] = bfvals_int[0:npg,kxx]*bfvals_int[0:npg,ky]     # basis function = P_2(x)y
    bfvals_int[0:npg,kxyy] = bfvals_int[0:npg,kx]*bfvals_int[0:npg,kyy]     # basis function = P_2(y)x
    bfvals_int[0:npg,kyyzz] = bfvals_int[0:npg,kyy]*bfvals_int[0:npg,kzz]     # basis function = P_2(y)P_2(z)
    bfvals_int[0:npg,kzzxx] = bfvals_int[0:npg,kzz]*bfvals_int[0:npg,kxx]     # basis function = P_2(z)P_2(x)
    bfvals_int[0:npg,kxxyy] = bfvals_int[0:npg,kxx]*bfvals_int[0:npg,kyy]     # basis function = P_2(x)P_2(y)
    bfvals_int[0:npg,kyzxx] = bfvals_int[0:npg,kyz]*bfvals_int[0:npg,kxx]     # basis function = yz P_2(x)
    bfvals_int[0:npg,kzxyy] = bfvals_int[0:npg,kzx]*bfvals_int[0:npg,kyy]     # basis function = zx P_2(y)
    bfvals_int[0:npg,kxyzz] = bfvals_int[0:npg,kxy]*bfvals_int[0:npg,kzz]     # basis function = xy P_2(z)
    bfvals_int[0:npg,kxyyzz] = bfvals_int[0:npg,kx]*bfvals_int[0:npg,kyy]*bfvals_int[0:npg,kzz]   # basis function = x P_2(y)P_2(z)
    bfvals_int[0:npg,kyzzxx] = bfvals_int[0:npg,ky]*bfvals_int[0:npg,kzz]*bfvals_int[0:npg,kxx]   # basis function = y P_2(z)P_2(x)
    bfvals_int[0:npg,kzxxyy] = bfvals_int[0:npg,kz]*bfvals_int[0:npg,kxx]*bfvals_int[0:npg,kyy]   # basis function = z P_2(x)P_2(y)
    bfvals_int[0:npg,kxxyyzz] = bfvals_int[0:npg,kx]*bfvals_int[0:npg,kyy]*bfvals_int[0:npg,kzz]  # basis function = P_2(x)P_2(y)P_2(z)

end subroutine set_internal_vals_3D

#----------------------------------------------------------
def set_face_vals_3D():
    # Define local basis function values at quadrature points on a cell face.
    # Used in flux_cal() and prepare_exchange() for interpolating from cell center onto cell face.
    # Also used in glflux().

    # integer ixyz,i
    # real c15d5,c1dsq3,xq4p,xq4m

    c15d5 = sqrt(15.)/5.
    c1dsq3 = 1./sqrt(3.)

    xq4p = sqrt(525. + 70.*sqrt(30.))/35.
    xq4m = sqrt(525. - 70.*sqrt(30.))/35.

    if iquad == 2:
        xquad[0] = -c1dsq3
        xquad[1] = c1dsq3
    if iquad == 3:
        xquad[0] = -c15d5
        xquad[1] = 0.
        xquad[2] = c15d5
    if iquad == 4:
        xquad[0] = -xq4p
        xquad[1] = -xq4m
        xquad[2] = xq4m
        xquad[3] = xq4p

    # bfvals_rsp:  positive rs-face
    # bfvals_rsm:  negative rs-face
    # bfvals_rsp(,1):  value=1 on positive rs-face
    # bfvals_rsp(,kx):  x-values on positive rs-face
    # bfvals_rsp(,ky):  y-values on positive rs-face
    # bfvals_rsp(,kz):  z-values on positive rs-face
    # bfvals_rsp(,kxx):  P_2(x)-values on positive rs-face
    # bfvals_rsp(,kyy):  P_2(y)-values on positive rs-face
    # bfvals_rsp(,kzz):  P_2(z)-values on positive rs-face
    # bfvals_rsp(,kyz):  yz-values on positive rs-face
    # bfvals_rsp(,kzx):  zx-values on positive rs-face
    # bfvals_rsp(,kxy):  xy-values on positive rs-face

    bfvals_zp(1:nface,1) = 1.
    if iquad == 1:
        bfvals_zp(1,kx) = 0.
        bfvals_zp(1,ky) = 0.

    if iquad > 1:
        do i=1,nedge
            bfvals_zp((i-1)*nedge+1:i*nedge,kx) = xquad(1:nedge)
            bfvals_zp((i-1)*nedge+1:i*nedge,ky) = xquad(i)
        end do

    bfvals_zp[0:nface,kz] = 1.
    bfvals_zp(1:nface,kyz) = bfvals_zp[0:nface,ky]*bfvals_zp[0:nface,kz]
    bfvals_zp(1:nface,kzx) = bfvals_zp[0:nface,kz]*bfvals_zp[0:nface,kx]
    bfvals_zp(1:nface,kxy) = bfvals_zp[0:nface,kx]*bfvals_zp[0:nface,ky]
    bfvals_zp(1:nface,kxyz) = bfvals_zp[0:nface,kx]*bfvals_zp[0:nface,ky]*bfvals_zp[0:nface,kz]     # basis function = xyz

    bfvals_yp[0:nface,0] = 1.
    bfvals_yp[0:nface,kx] = bfvals_zp[0:nface,ky]
    bfvals_yp[0:nface,ky] = 1.
    bfvals_yp[0:nface,kz] = bfvals_zp[0:nface,kx]
    bfvals_yp[0:nface,kyz] = bfvals_yp[0:nface,ky]*bfvals_yp[0:nface,kz]
    bfvals_yp[0:nface,kzx] = bfvals_yp[0:nface,kz]*bfvals_yp[0:nface,kx]
    bfvals_yp[0:nface,kxy] = bfvals_yp[0:nface,kx]*bfvals_yp[0:nface,ky]
    bfvals_yp[0:nface,kxyz] = bfvals_yp[0:nface,kx]*bfvals_yp[0:nface,ky]*bfvals_yp[0:nface,kz]     # basis function = xyz

    bfvals_xp[0:nface,0] = 1.
    bfvals_xp[0:nface,kx] = 1.
    bfvals_xp[0:nface,ky] = bfvals_zp[0:nface,kx]
    bfvals_xp[0:nface,kz] = bfvals_zp[0:nface,ky]
    bfvals_xp[0:nface,kyz] = bfvals_xp[0:nface,ky]*bfvals_xp[0:nface,kz]
    bfvals_xp[0:nface,kzx] = bfvals_xp[0:nface,kz]*bfvals_xp[0:nface,kx]
    bfvals_xp[0:nface,kxy] = bfvals_xp[0:nface,kx]*bfvals_xp[0:nface,ky]
    bfvals_xp[0:nface,kxyz] = bfvals_xp[0:nface,kx]*bfvals_xp[0:nface,ky]*bfvals_xp[0:nface,kz]     # basis function = xyz

    do i=1,3
        bfvals_xp(1:nface,kxx+i-1) = 1.5*bfvals_xp(1:nface,kx+i-1)**2 - 0.5
        bfvals_yp(1:nface,kxx+i-1) = 1.5*bfvals_yp(1:nface,kx+i-1)**2 - 0.5
        bfvals_zp(1:nface,kxx+i-1) = 1.5*bfvals_zp(1:nface,kx+i-1)**2 - 0.5
        # bfvals_xp(1:nface,kxxx+i-1) = 2.5*bfvals_xp(1:nface,kx+i-1)**3 - 1.5*bfvals_xp(1:nface,kx+i-1)
        # bfvals_yp(1:nface,kxxx+i-1) = 2.5*bfvals_yp(1:nface,kx+i-1)**3 - 1.5*bfvals_yp(1:nface,kx+i-1)
        # bfvals_zp(1:nface,kxxx+i-1) = 2.5*bfvals_zp(1:nface,kx+i-1)**3 - 1.5*bfvals_zp(1:nface,kx+i-1)
    end do


    bfvals_xp[0:nface,kyyz] = bfvals_xp[0:nface,kyy]*bfvals_xp[0:nface,kz]     # basis function = P_2(y)z
    bfvals_xp[0:nface,kyzz] = bfvals_xp[0:nface,ky]*bfvals_xp[0:nface,kzz]     # basis function = P_2(z)y
    bfvals_xp[0:nface,kzzx] = bfvals_xp[0:nface,kzz]*bfvals_xp[0:nface,kx]     # basis function = P_2(z)x
    bfvals_xp[0:nface,kzxx] = bfvals_xp[0:nface,kz]*bfvals_xp[0:nface,kxx]     # basis function = P_2(x)z
    bfvals_xp[0:nface,kxxy] = bfvals_xp[0:nface,kxx]*bfvals_xp[0:nface,ky]     # basis function = P_2(x)y
    bfvals_xp[0:nface,kxyy] = bfvals_xp[0:nface,kx]*bfvals_xp[0:nface,kyy]     # basis function = P_2(y)x
    bfvals_xp[0:nface,kyyzz] = bfvals_xp[0:nface,kyy]*bfvals_xp[0:nface,kzz]     # basis function = P_2(y)P_2(z)
    bfvals_xp[0:nface,kzzxx] = bfvals_xp[0:nface,kzz]*bfvals_xp[0:nface,kxx]     # basis function = P_2(z)P_2(x)
    bfvals_xp[0:nface,kxxyy] = bfvals_xp[0:nface,kxx]*bfvals_xp[0:nface,kyy]     # basis function = P_2(x)P_2(y)
    bfvals_xp[0:nface,kyzxx] = bfvals_xp[0:nface,kyz]*bfvals_xp[0:nface,kxx]     # basis function = yz P_2(x)
    bfvals_xp[0:nface,kzxyy] = bfvals_xp[0:nface,kzx]*bfvals_xp[0:nface,kyy]     # basis function = zx P_2(y)
    bfvals_xp[0:nface,kxyzz] = bfvals_xp[0:nface,kxy]*bfvals_xp[0:nface,kzz]     # basis function = xy P_2(z)
    bfvals_xp[0:nface,kxyyzz] = bfvals_xp[0:nface,kx]*bfvals_xp[0:nface,kyy]*bfvals_xp[0:nface,kzz]   # basis function = x P_2(y)P_2(z)
    bfvals_xp[0:nface,kyzzxx] = bfvals_xp[0:nface,ky]*bfvals_xp[0:nface,kzz]*bfvals_xp[0:nface,kxx]   # basis function = y P_2(z)P_2(x)
    bfvals_xp[0:nface,kzxxyy] = bfvals_xp[0:nface,kz]*bfvals_xp[0:nface,kxx]*bfvals_xp[0:nface,kyy]   # basis function = z P_2(x)P_2(y)
    bfvals_xp[0:nface,kxxyyzz] = bfvals_xp[0:nface,kx]*bfvals_xp[0:nface,kyy]*bfvals_xp[0:nface,kzz]  # basis function = P_2(x)P_2(y)P_2(z)

    bfvals_xm = bfvals_xp
    bfvals_xm[0:nface,kx] = -bfvals_xp[0:nface,kx]
    bfvals_xm[0:nface,kzx] = -bfvals_xp[0:nface,kzx]
    bfvals_xm[0:nface,kxy] = -bfvals_xp[0:nface,kxy]
    bfvals_xm[0:nface,kxyz] = -bfvals_xp[0:nface,kxyz]
    bfvals_xm[0:nface,kzzx] = -bfvals_xp[0:nface,kzzx]
    bfvals_xm[0:nface,kxyy] = -bfvals_xp[0:nface,kxyy]
    # bfvals_xm(1:nface,kxxx) = -bfvals_xp(1:nface,kxxx)
    bfvals_xm[0:nface,kzxyy] = -bfvals_xp[0:nface,kzxyy]
    bfvals_xm[0:nface,kxyzz] = -bfvals_xp[0:nface,kxyzz]
    bfvals_xm[0:nface,kxyyzz] = -bfvals_xp[0:nface,kxyyzz]

    bfvals_yp[0:nface,kyyz] = bfvals_yp[0:nface,kyy]*bfvals_yp[0:nface,kz]     # basis function = P_2(y)z
    bfvals_yp[0:nface,kyzz] = bfvals_yp[0:nface,ky]*bfvals_yp[0:nface,kzz]     # basis function = P_2(z)y
    bfvals_yp[0:nface,kzzx] = bfvals_yp[0:nface,kzz]*bfvals_yp[0:nface,kx]     # basis function = P_2(z)x
    bfvals_yp[0:nface,kzxx] = bfvals_yp[0:nface,kz]*bfvals_yp[0:nface,kxx]     # basis function = P_2(x)z
    bfvals_yp[0:nface,kxxy] = bfvals_yp[0:nface,kxx]*bfvals_yp[0:nface,ky]     # basis function = P_2(x)y
    bfvals_yp[0:nface,kxyy] = bfvals_yp[0:nface,kx]*bfvals_yp[0:nface,kyy]     # basis function = P_2(y)x
    bfvals_yp[0:nface,kyyzz] = bfvals_yp[0:nface,kyy]*bfvals_yp[0:nface,kzz]     # basis function = P_2(y)P_2(z)
    bfvals_yp[0:nface,kzzxx] = bfvals_yp[0:nface,kzz]*bfvals_yp[0:nface,kxx]     # basis function = P_2(z)P_2(x)
    bfvals_yp[0:nface,kxxyy] = bfvals_yp[0:nface,kxx]*bfvals_yp[0:nface,kyy]     # basis function = P_2(x)P_2(y)
    bfvals_yp[0:nface,kyzxx] = bfvals_yp[0:nface,kyz]*bfvals_yp[0:nface,kxx]     # basis function = yz P_2(x)
    bfvals_yp[0:nface,kzxyy] = bfvals_yp[0:nface,kzx]*bfvals_yp[0:nface,kyy]     # basis function = zx P_2(y)
    bfvals_yp[0:nface,kxyzz] = bfvals_yp[0:nface,kxy]*bfvals_yp[0:nface,kzz]     # basis function = xy P_2(z)
    bfvals_yp[0:nface,kxyyzz] = bfvals_yp[0:nface,kx]*bfvals_yp[0:nface,kyy]*bfvals_yp[0:nface,kzz]   # basis function = x P_2(y)P_2(z)
    bfvals_yp[0:nface,kyzzxx] = bfvals_yp[0:nface,ky]*bfvals_yp[0:nface,kzz]*bfvals_yp[0:nface,kxx]   # basis function = y P_2(z)P_2(x)
    bfvals_yp[0:nface,kzxxyy] = bfvals_yp[0:nface,kz]*bfvals_yp[0:nface,kxx]*bfvals_yp[0:nface,kyy]   # basis function = z P_2(x)P_2(y)
    bfvals_yp[0:nface,kxxyyzz] = bfvals_yp[0:nface,kx]*bfvals_yp[0:nface,kyy]*bfvals_yp[0:nface,kzz]  # basis function = P_2(x)P_2(y)P_2(z)

    bfvals_ym = bfvals_yp
    bfvals_ym[0:nface,ky] = -bfvals_yp[0:nface,ky]
    bfvals_ym[0:nface,kyz] = -bfvals_yp[0:nface,kyz]
    bfvals_ym[0:nface,kxy] = -bfvals_yp[0:nface,kxy]
    bfvals_ym[0:nface,kxyz] = -bfvals_yp[0:nface,kxyz]
    bfvals_ym[0:nface,kyzz] = -bfvals_yp[0:nface,kyzz]
    bfvals_ym[0:nface,kxxy] = -bfvals_yp[0:nface,kxxy]
    # bfvals_ym(1:nface,kyyy) = -bfvals_yp(1:nface,kyyy)
    bfvals_ym[0:nface,kyzxx] = -bfvals_yp[0:nface,kyzxx]
    bfvals_ym[0:nface,kxyzz] = -bfvals_yp[0:nface,kxyzz]
    bfvals_ym[0:nface,kyzzxx] = -bfvals_yp[0:nface,kyzzxx]

    bfvals_zp[0:nface,kyyz] = bfvals_zp[0:nface,kyy]*bfvals_zp[0:nface,kz]     # basis function = P_2(y)z
    bfvals_zp[0:nface,kyzz] = bfvals_zp[0:nface,ky]*bfvals_zp[0:nface,kzz]     # basis function = P_2(z)y
    bfvals_zp[0:nface,kzzx] = bfvals_zp[0:nface,kzz]*bfvals_zp[0:nface,kx]     # basis function = P_2(z)x
    bfvals_zp[0:nface,kzxx] = bfvals_zp[0:nface,kz]*bfvals_zp[0:nface,kxx]     # basis function = P_2(x)z
    bfvals_zp[0:nface,kxxy] = bfvals_zp[0:nface,kxx]*bfvals_zp[0:nface,ky]     # basis function = P_2(x)y
    bfvals_zp[0:nface,kxyy] = bfvals_zp[0:nface,kx]*bfvals_zp[0:nface,kyy]     # basis function = P_2(y)x
    bfvals_zp[0:nface,kyyzz] = bfvals_zp[0:nface,kyy]*bfvals_zp[0:nface,kzz]     # basis function = P_2(y)P_2(z)
    bfvals_zp[0:nface,kzzxx] = bfvals_zp[0:nface,kzz]*bfvals_zp[0:nface,kxx]     # basis function = P_2(z)P_2(x)
    bfvals_zp[0:nface,kxxyy] = bfvals_zp[0:nface,kxx]*bfvals_zp[0:nface,kyy]     # basis function = P_2(x)P_2(y)
    bfvals_zp[0:nface,kyzxx] = bfvals_zp[0:nface,kyz]*bfvals_zp[0:nface,kxx]     # basis function = yz P_2(x)
    bfvals_zp[0:nface,kzxyy] = bfvals_zp[0:nface,kzx]*bfvals_zp[0:nface,kyy]     # basis function = zx P_2(y)
    bfvals_zp[0:nface,kxyzz] = bfvals_zp[0:nface,kxy]*bfvals_zp[0:nface,kzz]     # basis function = xy P_2(z)
    bfvals_zp[0:nface,kxyyzz] = bfvals_zp[0:nface,kx]*bfvals_zp[0:nface,kyy]*bfvals_zp[0:nface,kzz]   # basis function = x P_2(y)P_2(z)
    bfvals_zp[0:nface,kyzzxx] = bfvals_zp[0:nface,ky]*bfvals_zp[0:nface,kzz]*bfvals_zp[0:nface,kxx]   # basis function = y P_2(z)P_2(x)
    bfvals_zp[0:nface,kzxxyy] = bfvals_zp[0:nface,kz]*bfvals_zp[0:nface,kxx]*bfvals_zp[0:nface,kyy]   # basis function = z P_2(x)P_2(y)
    bfvals_zp[0:nface,kxxyyzz] = bfvals_zp[0:nface,kx]*bfvals_zp[0:nface,kyy]*bfvals_zp[0:nface,kzz]  # basis function = P_2(x)P_2(y)P_2(z)


    bfvals_zm = bfvals_zp
    bfvals_zm[0:nface,kz] = -bfvals_zp[0:nface,kz]
    bfvals_zm[0:nface,kyz] = -bfvals_zp[0:nface,kyz]
    bfvals_zm[0:nface,kzx] = -bfvals_zp[0:nface,kzx]
    bfvals_zm[0:nface,kxyz] = -bfvals_zp[0:nface,kxyz]
    bfvals_zm[0:nface,kyyz] = -bfvals_zp[0:nface,kyyz]
    bfvals_zm[0:nface,kzxx] = -bfvals_zp[0:nface,kzxx]
    # bfvals_zm(1:nface,kzzz) = -bfvals_zp(1:nface,kzzz)
    bfvals_zm[0:nface,kyzxx] = -bfvals_zp[0:nface,kyzxx]
    bfvals_zm[0:nface,kzxyy] = -bfvals_zp[0:nface,kzxyy]
    bfvals_zm[0:nface,kzxxyy] = -bfvals_zp[0:nface,kzxxyy]

    # Organize local basis values on faces into 1-D vectors.
    # Used in limiter() and max_lim().

    bf_faces[0:nslim,0] = 1.        # basis function = 1

    do ixyz=kx,kz
        bf_faces[0:nface,ixyz] = bfvals_xm[0:nface,ixyz]        # basis function = x,y,z
        bf_faces(nface+1:2*nface,ixyz) = bfvals_xp(1:nface,ixyz)        # basis function = x,y,z
        bf_faces(2*nface+1:3*nface,ixyz) = bfvals_ym[0:nface,ixyz]      # basis function = x,y,z
        bf_faces(3*nface+1:4*nface,ixyz) = bfvals_yp[0:nface,ixyz]      # basis function = x,y,z
        bf_faces(4*nface+1:5*nface,ixyz) = bfvals_zm[0:nface,ixyz]      # basis function = x,y,z
        bf_faces(5*nface+1:6*nface,ixyz) = bfvals_zp[0:nface,ixyz]      # basis function = x,y,z
        bf_faces(6*nface+1:nslim,ixyz) = bfvals_int[0:npg,ixyz]     # basis function = x,y,z
    end do

    bf_faces[0:nslim,kyz] = bf_faces[0:nslim,ky]*bf_faces[0:nslim,kz]     # basis function = yz
    bf_faces[0:nslim,kzx] = bf_faces[0:nslim,kz]*bf_faces[0:nslim,kx]     # basis function = zx
    bf_faces[0:nslim,kxy] = bf_faces[0:nslim,kx]*bf_faces[0:nslim,ky]     # basis function = xy
    bf_faces[0:nslim,kxyz] = bf_faces[0:nslim,kx]*bf_faces[0:nslim,ky]*bf_faces[0:nslim,kz]     # basis function = xyz

    for i in xrange(3):
        bf_faces[0:nslim,kxx+i] = 1.5*bf_faces[0:nslim,kx+i)]*2 - 0.5    # basis function = P_2(s)
        # bf_faces(1:nslim,kxxx+i) = 2.5*bf_faces(1:nslim,kx+i)**3 - 1.5*bf_faces(1:nslim,kx+i)   # basis function = P_3(s)

    bf_faces[0:nslim,kyyz] = bf_faces[0:nslim,kyy]*bf_faces[0:nslim,kz]     # basis function = P_2(y)z
    bf_faces[0:nslim,kyzz] = bf_faces[0:nslim,ky]*bf_faces[0:nslim,kzz]     # basis function = P_2(z)y
    bf_faces[0:nslim,kzzx] = bf_faces[0:nslim,kzz]*bf_faces[0:nslim,kx]     # basis function = P_2(z)x
    bf_faces[0:nslim,kzxx] = bf_faces[0:nslim,kz]*bf_faces[0:nslim,kxx]     # basis function = P_2(x)z
    bf_faces[0:nslim,kxxy] = bf_faces[0:nslim,kxx]*bf_faces[0:nslim,ky]     # basis function = P_2(x)y
    bf_faces[0:nslim,kxyy] = bf_faces[0:nslim,kx]*bf_faces[0:nslim,kyy]     # basis function = P_2(y)x
    bf_faces[0:nslim,kyyzz] = bf_faces[0:nslim,kyy]*bf_faces[0:nslim,kzz]     # basis function = P_2(y)P_2(z)
    bf_faces[0:nslim,kzzxx] = bf_faces[0:nslim,kzz]*bf_faces[0:nslim,kxx]     # basis function = P_2(z)P_2(x)
    bf_faces[0:nslim,kxxyy] = bf_faces[0:nslim,kxx]*bf_faces[0:nslim,kyy]     # basis function = P_2(x)P_2(y)
    bf_faces[0:nslim,kyzxx] = bf_faces[0:nslim,kyz]*bf_faces[0:nslim,kxx]     # basis function = yz P_2(x)
    bf_faces[0:nslim,kzxyy] = bf_faces[0:nslim,kzx]*bf_faces[0:nslim,kyy]     # basis function = zx P_2(y)
    bf_faces[0:nslim,kxyzz] = bf_faces[0:nslim,kxy]*bf_faces[0:nslim,kzz]     # basis function = xy P_2(z)
    bf_faces[0:nslim,kxyyzz] = bf_faces[0:nslim,kx]*bf_faces[0:nslim,kyy]*bf_faces[0:nslim,kzz]   # basis function = x P_2(y)P_2(z)
    bf_faces[0:nslim,kyzzxx] = bf_faces[0:nslim,ky]*bf_faces[0:nslim,kzz]*bf_faces[0:nslim,kxx]   # basis function = y P_2(z)P_2(x)
    bf_faces[0:nslim,kzxxyy] = bf_faces[0:nslim,kz]*bf_faces[0:nslim,kxx]*bf_faces[0:nslim,kyy]   # basis function = z P_2(x)P_2(y)
    bf_faces[0:nslim,kxxyyzz] = bf_faces[0:nslim,kx]*bf_faces[0:nslim,kyy]*bf_faces[0:nslim,kzz]  # basis function = P_2(x)P_2(y)P_2(z)

end subroutine set_face_vals_3D

#----------------------------------------------------------
def set_weights_3D():

    # integer i
    # real wq4p,wq4m

    wq4p = (18. + sqrt(30.))/36.
    wq4m = (18. - sqrt(30.))/36.

    # Define weights for integral approximation using Gaussian quadrature.

    #  Define weights for 1-D integration

    if iquad == 1:     # 1-point quadrature
        wgt1d[0] = 2.
    if iquad == 2:     # 2-point quadrature
        wgt1d[0:2] = 1.
    if iquad == 3:     # 3-point quadrature
        wgt1d[0] = 5./9.
        wgt1d[1] = 8./9.
        wgt1d[2] = 5./9.
    if iquad == 4:     # 4-point quadrature
        wgt1d[0] = wq4m
        wgt1d[1] = wq4p
        wgt1d[2] = wq4p
        wgt1d[3] = wq4m

    #  Define weights for 2-D integration

    if iquad == 1:     # 1-point quadrature
        wgt2d[0] = 4.

    if iquad == 2:     # 2-point quadrature
        wgt2d(1:4) = 1.

    if iquad >= 3:
        do i= 1,nedge
            wgt2d((i-1)*nedge+1:i*nedge) = wgt1d(1:nedge)*wgt1d(i)
        end do

    #  Define weights for 3-D integration

    if iquad == 1:     # 1-point quadrature
        wgt3d(1) = 8.

    if iquad == 2:     # 2-point quadrature
        wgt3d(1:8) = 1.

    if iquad >= 3:
        do i= 1,nedge
            wgt3d((i-1)*nface+1:i*nface) = wgt2d(1:nface)*wgt1d(i)
        end do

end subroutine set_weights_3D

#--------------------------------------------------------------------------------

end program
