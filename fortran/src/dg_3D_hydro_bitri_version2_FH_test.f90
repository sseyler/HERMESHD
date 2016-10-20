!      DG-Hydro in 3 dimensions: Extension of the finite volume PERSEUS Algorithm (Physics of the Extended-mhd Relaxation System using an Efficient Upwind Scheme,
!                                by Matt Martin and Charles Seyler) to the Discontinuous Galerkin Method specialized to hydrodynamics.
!
!                DG Algorithm by Xuan Zhao, Nat Hamlin and Charles Seyler.
!
!                Solves the 3D compressible Euler equations
!
!                Variables:    There are 5 dependent variables: density (rh), velocity (vx,vy,vz), energy (en)
!
!                Units:        A value of unity for each variable or parameter corresponds to the following dimensional units:
!
!                              Length              L0
!                              Time                t0
!                              Number density      n0
!                              Velocity            v0
!                              Temperature         te0
!
!-------------------------------------------------------------------------------------------------------------------------------------------

program main
use lib_vtk_io
include '/nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include/mpif.h'
include '/nfs/packages/opt/Linux_x86_64/intel/17.0/fortran/compilers_and_libraries_2017.0.098/linux/mkl/include'

integer, parameter :: rh=1, mx=2, my=3, mz=4, en=5, pxx=6, pyy=7, pzz=8, pxy=9, pxz=10, pyz=11, nQ=11

integer, parameter :: nx = 20, ny = 20, nz = 1, ngu = 0, nbasis = 8, nbastot = 27
! nbasis = 4: {1,x,y,z}
! nbasis = 10: {1,x,y,z, P_2(x),P_2(y),P_2(z), yz, zx, xy}
! nbasis = 20: nbasis10 + {xyz,xP2(y),yP2(x),xP2(z),zP2(x),yP2(z),zP2(y),P3(x),P3(y),P3(z)}

! The jump in accuracy between the linear basis (nbasis = 4) and quadratic basis (nbasis = 10)
! is much greater than the jump in accuracy between quadratic and cubic (nbasis = 20) bases.

integer, parameter :: iquad=2, nedge = iquad
! iquad: number of Gaussian quadrature points per direction. iquad should not be
! less than ipoly, where ipoly is the maximum Legendre polynomial order used,
! otherwise instability results. iquad should not be larger than ipoly + 1,
! which is an exact Gaussian quadrature for the Legendre polynomial used. Thus
! there are only two cases of iquad for a given nbasis.  Both cases give similar
! results although the iquad = ipoly + 1 case is formally more accurate.

integer, parameter :: nface = iquad*iquad, npg = nface*iquad, nfe = 2*nface, npge = 6*nface, nslim=npg+6*nface
! nface: number of quadrature points per cell face.
! npg: number of internal points per cell.

integer, parameter :: xlbc = 2, xhbc = 2, ylbc = 2, yhbc = 2, zlbc = 2, zhbc = 2
! Boundary condition parameters: if  = 2 then periodic.  MPI does this for you.
! If  = 0, then the set_bc subroutine is used to prescribe BCs

integer, parameter :: ntout = 100, iorder = 2

integer, parameter :: ihllc = 1, iroe = 0, ieos = 1
! Choose Riemann solver for computing fluxes.  Set chosen solver = 1.
! If all of them = 0, then LLF is used for fluxes.
    ! LLF is very diffusive for the hydro problem.  Roe and HLLC are much less
    ! diffusive than LLF and give very similar results with similar cpu overhead
    ! Only HLLC is setup to handle water EOS (ieos = 2)

! to restart from a checkpoint, set iread to 1 or 2 (when using the odd/even scheme)
integer, parameter :: iread = 0, iwrite = 0
character (4), parameter :: fpre = 'Qout'
logical, parameter :: resuming = .false.

real, parameter :: lx = 100., ly = 100., lz = 100./120.

real, parameter :: tf = 1000.

    real, parameter :: pi = 4.0*atan(1.0)

    real, parameter :: aindex = 5./3., mu = 18.

real, parameter :: aindm1 = aindex - 1.0, cp = aindex/(aindex - 1.0), clt = 2., vis = 0.e-1, epsi = 5., amplv = 0.0, amplm = 0.0, amplen = 10.

    ! dimensional units (expressed in MKS)
    real, parameter :: L0=1.0e-9, t0=1.0e-12, n0=3.32e28

    ! derived units
    real, parameter :: v0=L0/t0, p0=mu*1.67e-27*n0*v0**2, te0=p0/n0/1.6e-19, P_1 = 2.15e9/7.2/p0, P_base = 1.01e5/p0 ! atmospheric pressure

! rh_min is a minimum density to be used for ideal gas law EOS and rh_min is the
! minimum density below which the pressure becomes negative for the water EOS:
!  P = P_1*(density**7.2 - 1.) + P_base
! The DG based subroutine "limiter" prevents the density from falling below
! rh_mult*rh_min.
! Note: the EOS for water is likely to be a critical player in getting the
! fluctuating hydrodynamics correct. The EOS used here is a simple one called
! Tait-Murnaghan. There are other much more sophisicated EOS's, some of which
! take into account ionic solutions. It would be worthwhile to investigate this
! further and experiment with different EOS's.

real, parameter :: rh_floor = 1.e-1, rh_mult = 1.01, rh_min = rh_mult*(1. - P_base/P_1)**(1./7.2)

real, parameter :: T_floor = 0.026/te0, P_floor = T_floor*rh_floor

real, dimension(nx,ny,nz,nQ,nbasis) ::  Q_r0, Q_r1, Q_r2, Q_r3
real, dimension(nx,ny,nz,nQ,nbasis) ::  glflux_r, source_r
real, dimension(nx,ny,nz,nQ,nbasis) :: integral_r

    real eta(nx,ny,nz,npg),den0(nx,ny,nz),Ez0,Zdy(nx,ny,nz,npg)
real flux_x(nface,1:nx+1,ny,nz,1:nQ), flux_y(nface,nx,1:ny+1,nz,1:nQ), flux_z(nface,nx,ny,1:nz+1,1:nQ)
    real cfrx(nface,nQ),cfry(nface,nQ),cfrz(nface,nQ)
logical MMask(nx,ny,nz),BMask(nx,ny,nz)
real xcell(npg), ycell(npg), zcell(npg), xface(npge), yface(npge), zface(npge)
integer ticks, count_rate, count_max
real t1, t2, t3, t4, elapsed_time, t_start, t_stop, dtoriginal
    real t,dt,tout,dtout,vf,dxi,dyi,dzi,loc_lxd,loc_lyd,loc_lzd,check_Iz,sl, dz, dy, dx
    real lxd,lxu,lyd,lyu,lzd,lzu,pin_rad,pin_height,foil_thickness,rh_foil,rh_fluid,pin_rad_in,pin_rad_out,rim_rad
    real disk_rad,disk_height,foil_rad,buf_rad,buf_z,dish_height,foil_height
    real gpz_rad,rh_gpz,kappa
real Qxhigh_ext(ny,nz,nface,nQ), Qxlow_int(ny,nz,nface,nQ), Qxlow_ext(ny,nz,nface,nQ), Qxhigh_int(ny,nz,nface,nQ)
real Qyhigh_ext(nx,nz,nface,nQ), Qylow_int(nx,nz,nface,nQ), Qylow_ext(nx,nz,nface,nQ), Qyhigh_int(nx,nz,nface,nQ)
real Qzhigh_ext(nx,ny,nface,nQ), Qzlow_int(nx,ny,nface,nQ), Qzlow_ext(nx,ny,nface,nQ), Qzhigh_int(nx,ny,nface,nQ)
integer mxa(3),mya(3),mza(3),kroe(nface),niter,iseed


! Parameters relating to quadratures and basis functions.

real wgt1d(5), wgt2d(30), wgt3d(100), cbasis(nbastot)
! wgt1d: quadrature weights for 1-D integration
! wgt2d: quadrature weights for 2-D integration
! wgt3d: quadrature weights for 3-D integration

real, dimension(nface,nbastot) :: bfvals_zp, bfvals_zm, bfvals_yp, bfvals_ym, bfvals_xp, bfvals_xm
real bf_faces(nslim,nbastot), bfvals_int(npg,nbastot),xquad(20)
    real bval_int_wgt(npg,nbastot),wgtbfvals_xp(nface,nbastot),wgtbfvals_yp(nface,nbastot),wgtbfvals_zp(nface,nbastot)
    real wgtbfvals_xm(nface,nbastot),wgtbfvals_ym(nface,nbastot),wgtbfvals_zm(nface,nbastot)
    real wgtbf_xmp(4,2,nbastot),wgtbf_ymp(4,2,nbastot),wgtbf_zmp(4,2,nbastot)
    real sumx,sumy,sumz
    integer i2f,i01,i2fa

    integer,parameter :: kx=2,ky=3,kz=4,kyz=5,kzx=6,kxy=7,kxyz=8,kxx=9,kyy=10,kzz=11,kyzz=12,kzxx=13,kxyy=14
    integer,parameter :: kyyz=15,kzzx=16,kxxy=17,kyyzz=18,kzzxx=19,kxxyy=20,kyzxx=21,kzxyy=22,kxyzz=23
    integer,parameter :: kxyyzz=24,kyzzxx=25,kzxxyy=26,kxxyyzz=27

! Parameters for VTK output.

integer, parameter :: nvtk=2
integer, parameter :: nvtk2=nvtk*nvtk, nvtk3=nvtk*nvtk*nvtk
integer(I4P), parameter :: nnx=nx*nvtk, nny=ny*nvtk, nnz=nz*nvtk
real, dimension(nvtk3,nbastot) :: bfvtk, bfvtk_dx, bfvtk_dy, bfvtk_dz
real xgrid(20),dxvtk,dyvtk,dzvtk


! MPI definitions

integer :: mpi_nx=4, mpi_ny=4, print_mpi=0

    integer iam,ierr,mpi_nz,numprocs,reorder,cartcomm,mpi_P,mpi_Q,mpi_R
    integer dims(3),coords(3),periods(3),nbrs(6),reqs(4),stats(MPI_STATUS_SIZE,4)
    integer,parameter:: NORTH=1,SOUTH=2,EAST=3,WEST=4,UP=5,DOWN=6,MPI_TT=MPI_REAL4

    real cflm

if (nbasis .le. 8) cflm = 0.14
if (nbasis .eq. 27) cflm = 0.1
if (nbasis .eq. 64) cflm = 0.08

! Initialize grid sizes and local lengths

cbasis(1) = 1.               ! coefficient for basis function {1}
cbasis(kx:kz) = 3.           ! coefficients for basis functions {x,y,z}
cbasis(kyz:kxy) = 9.         ! coefficients for basis functions {yz,zx,xy}
cbasis(kxyz) = 27.           ! coefficient for basis function {xyz}

cbasis(kxx:kzz) = 5.         ! coefficients for basis functions {P2(x),P2(y),P2(z)}
cbasis(kyzz:kxyy) = 15.      ! coefficients for basis functions {yP2(z),zP2(x),xP2(y)}
cbasis(kyyz:kxxy) = 15.      ! coefficients for basis functions {P2(y)z,P2(z)y,P2(z)x}
cbasis(kyyzz:kxxyy) = 25.    ! coefficients for basis functions {P2(y)P2(z),P2(z)P2(x),P2(x)P2(y)}
cbasis(kyzxx:kxyzz) = 45.    ! coefficients for basis functions {yzP_2(x),zxP_2(y),xyP_2(z)}
cbasis(kxyyzz:kzxxyy) = 75.  ! coefficients for basis functions {xP2(y)P2(z),yP2(z)P2(x),zP2(x)P2(y)}
cbasis(kxxyyzz) = 125.       ! coefficients for basis functions {P2(x)P2(y)P2(z)}

call MPI_Init ( ierr )
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

mpi_nz=numprocs/(mpi_nx*mpi_ny)

dims(1)=mpi_nx
dims(2)=mpi_ny
dims(3)=mpi_nz

periods(:)=0
if (xhbc .eq. 2) then
    periods(1)=1
end if
if (yhbc .eq. 2) then
    periods(2)=1
end if
if (zhbc .eq. 2) then
    periods(3)=1
end if
reorder = 1

call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, reorder,cartcomm, ierr)
call MPI_COMM_RANK (cartcomm, iam, ierr )
call MPI_CART_COORDS(cartcomm, iam, 3, coords, ierr)
mpi_P=coords(1)+1
mpi_Q=coords(2)+1
mpi_R=coords(3)+1
call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(WEST), nbrs(EAST), ierr)
call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(SOUTH), nbrs(NORTH), ierr)
call MPI_CART_SHIFT(cartcomm, 2, 1, nbrs(DOWN), nbrs(UP), ierr)

half_length = lx/2.0
lxd = -half_length
lxu = half_length
half_length = ly/2.0
lyd = -half_length
lyu = half_length
half_length = lz/2.0
lzd = -half_length
lzu = half_length

dxi = (nx*mpi_nx)/(lxu-lxd)
dyi = (ny*mpi_ny)/(lyu-lyd)
dzi = (nz*mpi_nz)/(lzu-lzd)

dx = 1./dxi
dy = 1./dyi
dz = 1./dzi
loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz

rh_fluid = 1.

! indices used in poynting() for computing Poynting Maxwell flux

mxa(1) = mx
mxa(2) = my
mxa(3) = mz
mya(1) = my
mya(2) = mz
mya(3) = mx
mza(1) = mz
mza(2) = mx
mza(3) = my

t = 0.
dt = cflm*dx/clt
dtoriginal = dt
nout = 0
niter = 0
dtout = tf/ntout

! Evaluate local cell values of basis functions on cell interior and faces.
! This is done for 1, 2, or 3 point Gaussian quadrature.
call set_bfvals_3D

do ir=1,nbasis
    do ipg=1,npg
        bval_int_wgt(ipg,ir) = wgt3d(ipg)*bfvals_int(ipg,ir)
    end do
end do

do ir=1,nbasis
    wgtbfvals_xp(1:nface,ir) = wgt2d(1:nface)*bfvals_xp(1:nface,ir)
    wgtbfvals_yp(1:nface,ir) = wgt2d(1:nface)*bfvals_yp(1:nface,ir)
    wgtbfvals_zp(1:nface,ir) = wgt2d(1:nface)*bfvals_zp(1:nface,ir)
    wgtbfvals_xm(1:nface,ir) = wgt2d(1:nface)*bfvals_xm(1:nface,ir)
    wgtbfvals_ym(1:nface,ir) = wgt2d(1:nface)*bfvals_ym(1:nface,ir)
    wgtbfvals_zm(1:nface,ir) = wgt2d(1:nface)*bfvals_zm(1:nface,ir)
end do

do ir=1,nbasis
    wgtbf_xmp(1:4,1,ir) = -0.25*cbasis(ir)*dxi*wgtbfvals_xm(1:nface,ir)
    wgtbf_ymp(1:4,1,ir) = -0.25*cbasis(ir)*dyi*wgtbfvals_ym(1:nface,ir)
    wgtbf_zmp(1:4,1,ir) = -0.25*cbasis(ir)*dzi*wgtbfvals_zm(1:nface,ir)
    wgtbf_xmp(1:4,2,ir) = 0.25*cbasis(ir)*dxi*wgtbfvals_xp(1:nface,ir)
    wgtbf_ymp(1:4,2,ir) = 0.25*cbasis(ir)*dyi*wgtbfvals_yp(1:nface,ir)
    wgtbf_zmp(1:4,2,ir) = 0.25*cbasis(ir)*dzi*wgtbfvals_zp(1:nface,ir)
end do

call init_random_seed(123456789)
iseed = 1317345*mpi_P + 5438432*mpi_Q + 38472613*mpi_R

if (iam .eq. print_mpi) then
    print *,'total dim= ',mpi_nx*nx,mpi_ny*ny,mpi_nz*nz
    print *,'mpi dim= ',mpi_nx,mpi_ny,mpi_nz
print *, 'te0 is: ', te0
print *, 'dx is: ', ly/(ny*mpi_ny)*L0
    print *, 'iquad is: ',  iquad
    print *, 'nbasis is: ', nbasis
end if

call system_clock(ticks, count_rate, count_max)
t_start = ticks*1./count_rate

if (iread .eq. 0) then

    call initial_condition

else

    ! This applies only if the initial data are being read from an input file.
    ! - If resuming a run, keep the previous clock (i.e., t at nout) running.
    ! - If not resuming a run, treat input as initial conditions at t=0, nout=0.

    call readQ(fpre,iam,iread,Q_r0,t_p,dt_p,nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p)

    if (resuming) then
        t = t_p
        dt = dt_p
        nout = nout_p
    end if
    ! Note, nout=1 corresponds to t=dt, but nout=2 corresponds to t~dtout, etc.
    if (nout .gt. 1) then
        dtout_p = t_p/(nout_p-1)
    else  ! Automatically pass consistency check
        dtout_p = dtout
    end if
    if (iam .eq. print_mpi) then
        print *, 'resuming = ', resuming
        print *, 't = ', t
        print *, 'dt = ', dt
        print *, 'nout = ', nout
        print *, 'dtout_p = ', dtout_p, ' dtout = ', dtout
        print *, 'mpi_nx_p = ', mpi_nx_p, ' mpi_nx = ', mpi_nx
        print *, 'mpi_ny_p = ', mpi_ny_p, ' mpi_ny = ', mpi_ny
        print *, 'mpi_nz_p = ', mpi_nz_p, ' mpi_nz = ', mpi_nz
    end if
        ! Quit if dtout is incompatible with input t/(nout-1)
    if (abs(dtout_p-dtout)/dt_p > 1.01) then
        if (iam .eq. print_mpi) then
            print *, 'Bad restart, non-matching dtout'
        end if
        call exit(-1)
    end if
    if ((mpi_nx_p .ne. mpi_nx) .or. (mpi_ny_p .ne. mpi_ny) .or. (mpi_nz_p .ne. mpi_nz)) then
        if (iam .eq. print_mpi) then
            print *, 'Bad restart, non-matching mpi_nx, mpi_ny, or mpi_nz'
        end if
        call exit(-1)
    end if

end if


call system_clock( ticks, count_rate, count_max )
t1 = 1.*ticks / count_rate
call output_vtk(Q_r0,nout,iam)

do while (t<tf)

    ! if (mod(niter,200) .eq. 0 .and. iam .eq. print_mpi) print *,'niter,t,dt = ',niter,t,dt,dtout*nout
    niter = niter + 1
    call get_min_dt(dt)

    if (iorder .eq. 2) then
        call prep_advance(Q_r0)
        call advance_time_level_gl(Q_r0,Q_r1)

        call prep_advance(Q_r1)
        call advance_time_level_gl(Q_r1,Q_r2)

        Q_r0 = 0.5*(Q_r0 + Q_r2)
    end if

    if (iorder .eq. 3) then
        call prep_advance(Q_r0)
        call advance_time_level_gl(Q_r0,Q_r1)

        call prep_advance(Q_r1)
        call advance_time_level_gl(Q_r1,Q_r2)

        Q_r3 = 0.75*Q_r0 + 0.25*Q_r2

        call prep_advance(Q_r3)
        call advance_time_level_gl(Q_r3,Q_r2)

        Q_r0 = (Q_r0 + 2.*Q_r2)/3.
    end if

    t = t+dt

    if (t .gt. dtout*nout) then

        nout = nout+1
        if (iam .eq. print_mpi) then
            call system_clock(ticks, count_rate, count_max)
            t2 = 1.*ticks/count_rate
            print *, 'Iteration time', (t2-t1), 'seconds'
            t1 = t2
            print *, 'nout = ', nout
            print *, "        t= ",t*100.,"         dt= ",dt
        end if

        call MPI_BARRIER(cartcomm,ierr)
        call output_vtk(Q_r0,nout,iam)

        ! write checkpoint files; assign an odd/even id to ensure last two sets are kept
        if (iwrite .eq. 1) then
            ioe = 2 - mod(nout,2)
            call writeQ(fpre,iam,ioe,Q_r0,t,dt,nout,mpi_nx,mpi_ny,mpi_nz)
        end if

        call MPI_BARRIER(cartcomm,ierr)

        if (iam .eq. print_mpi) then
            call system_clock(ticks, count_rate, count_max)
            t2 = ticks/count_rate
            print *, 'Output time', (t2-t1), 'seconds'
            t1 = t2
        end if

    end if

end do

call MPI_Finalize (ierr)

contains

!------------------------------------------
subroutine init_random_seed(iseed)

    implicit none
    integer :: i, n, clock,iseed
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    if (iseed .eq. 0) then
        seed =  clock*(iam+1) + 37 * (/ (i - 1, i = 1, n) /)
    else
        seed =  iseed*(iam+1)
    endif
    call random_seed(put = seed)

    ! print *,seed(1)
    deallocate(seed)

end subroutine

!-------------------------------------------------------------------------------

subroutine prep_advance(Q_ri)

    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri

    if(ieos .eq. 2) call limiter(Q_ri)
    call prepare_exchange(Q_ri)
    call set_bc
    call flux_cal(Q_ri)
    call innerintegral(Q_ri)
    call glflux
    call source_calc(Q_ri,t)

end subroutine prep_advance

!-------------------------------------------------------------------------------

subroutine get_min_dt(dt)

   real dt,dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0,vex,vey,vez,vem,vem0,dni,dn,vx,vy,vz,Pr,sqdni,vacc,vacc0,cs
   integer :: i,j,k,main_proc=0,mpi_size=1,loc_reqs(numprocs-1),loc_stats(MPI_STATUS_SIZE,numprocs-1)

   vmag = 0.

    do k=1,nz
        do j=1,ny
            do i=1,nx

                dn = Q_r0(i,j,k,rh,1)
                dni = 1./dn
                vx = Q_r0(i,j,k,mx,1)*dni
                vy = Q_r0(i,j,k,my,1)*dni
                vz = Q_r0(i,j,k,mz,1)*dni
                if(ieos .eq. 1) cs = sqrt(aindex*(Q_r0(i,j,k,en,1)*dni - 0.5*(vx**2 + vy**2 + vz**2)))
                if(ieos .eq. 2) cs = sqrt(7.2*P_1*dn**6.2 + T_floor)

                vmag0 = max(abs(vx)+cs, abs(vy)+cs, abs(vz)+cs)
                if(vmag0 > vmag) vmag = vmag0

            end do
        end do
    end do

    vmax = vmag*dxi
    dt_min = 1*cflm/vmax  ! time step is determined by the maximum flow + sound speed in the domain

    call MPI_BARRIER(cartcomm,ierr)
    if (iam.eq.main_proc) then
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
    else
        call MPI_ISend(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
    call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
    call MPI_Wait(reqs(1),stats(:,1),ierr)
    endif

    dt = dt_min

end subroutine get_min_dt

!----------------------------------------------------------------------------------------------

subroutine initial_condition
    implicit none
    integer i,j,k
    real coeff
    real den, Pre, wtev

    wtev = T_floor

    Q_r0(:,:,:,:,:) = 0.0
    Q_r0(:,:,:,rh,1) = rh_floor
    Q_r0(:,:,:,en,1) = T_floor*rh_floor/(aindex - 1.)
    MMask(:,:,:) = .false.

    call fill_fluid

end subroutine initial_condition

!-------------------------------------------------------

subroutine fill_fluid

!-------------------------------------------------------

    integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,igrid,ieq,inds(nbasis)
    real wirer,x,y,x0,y0,rhline,wtev,rx,ry,rnum,w0
    real qquad(npg,nQ),xcc,ycc,zcc  ! bfint(npg,nbasis),qquadv(npg)
    iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

    wtev = T_floor
    w0 = 0.3

    ! test problem is an unstable flow jet in x with velocity perturbations in y

    do i = 1,nx
        do j = 1,ny
            do k = 1,nz

                call random_number(rand_number)
                rnum = (rand_number - 0.5)

                qquad(:,:) = 0.

                do igrid=1,npg

                    qquad(igrid,rh) = rh_floor
                    qquad(igrid,en) = wtev*qquad(igrid,rh)/(aindex - 1.)

                    xcc = xc(i) + bfvals_int(igrid,kx)*0.5/dxi
                    ycc = yc(j) + bfvals_int(igrid,ky)*0.5/dyi
                    zcc = zc(k) + bfvals_int(igrid,kz)*0.5/dzi

                    qquad(igrid,rh) = rh_fluid
                    qquad(igrid,my) = 0.001*rnum/1.
                    qquad(igrid,mz) = 0.
                    qquad(igrid,mx) = 0.5*qquad(igrid,rh)/cosh(20*ycc/lyu)/1.
                    qquad(igrid,en) = wtev*qquad(igrid,rh)/(aindex - 1.) + 0.5*(qquad(igrid,mx)**2 + qquad(igrid,my)**2 + qquad(igrid,mz)**2)/qquad(igrid,rh)

                end do

                do ieq=rh,en
                    do ir=1,nbasis
                        Q_r0(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(wgt3d(1:npg)*bfvals_int(1:npg,ir)*qquad(1:npg,ieq))
                    end do
                end do

            end do
        end do
    end do

end subroutine fill_fluid

!----------------------------------------------------------------------------------------------

subroutine set_bc

    implicit none
    real P1, P2, den, vz(nface)
    integer i,j,k,ne,ieq,l,i4

    if (mpi_P .eq. 1 .and. xlbc .ne. 2) then
        ! Set B.C.'s at bottom x-boundary.
        do k = 1,nz
            do j = 1,ny

                do i4=1,nface
                    Qxlow_ext(j,k,i4,:) = Qxlow_int(j,k,i4,:)
                end do
                if (maxval(Qxlow_ext(j,k,1:nface,mx)) .gt. 0.) then
                    Qxlow_ext(j,k,1:nface,mx) = 0.
                end if

            end do
        end do
    end if

!---------------------------------------------------------
    if (mpi_P .eq. mpi_nx .and. xlbc .ne. 2) then
        ! Set B.C.'s at top x-boundary.
        do k = 1,nz
            do j = 1,ny

                do i4=1,nface
                    Qxhigh_ext(j,k,i4,:) = Qxhigh_int(j,k,i4,:)
                end do
                if (minval(Qxhigh_ext(j,k,1:nface,mx)) .lt. 0.) then
                    Qxhigh_ext(j,k,1:nface,mx) = 0.
                end if

            end do
        end do
    end if

!----------------------------------------------------
    if (mpi_Q .eq. 1 .and. ylbc .ne. 2) then
        ! Set B.C.'s at bottom y-boundary.
        do k = 1,nz
            do i = 1,nx

                do i4=1,nface
                    Qylow_ext(i,k,i4,:) = Qylow_int(i,k,i4,:)
                end do
                if (maxval(Qylow_ext(i,k,1:nface,my)) .gt. 0.) then
                    Qylow_ext(i,k,1:nface,my) = 0.
                end if

            end do
        end do
    end if

!------------------------------------------------------
    if (mpi_Q .eq. mpi_ny .and. ylbc .ne. 2) then
        ! Set B.C.'s at top y-boundary.
        do k = 1,nz
            do i = 1,nx

                do i4=1,nface
                    Qyhigh_ext(i,k,i4,:) = Qyhigh_int(i,k,i4,:)
                end do
                if (minval(Qyhigh_ext(i,k,1:nface,my)) .lt. 0.) then
                    Qyhigh_ext(i,k,1:nface,my) = 0.
                end if

            end do
        end do
    end if

!--------------------------------------------------------
    if (mpi_R .eq. 1 .and. zlbc .ne. 2) then
        ! Set B.C.'s at bottom z-boundary.
        do j = 1,ny
            do i = 1,nx

                if (r(i,j) .ge. disk_rad .and. r(i,j) .le. lxu - 2.*dx) then
                    do i4=1,nface
                        Qzlow_ext(i,j,i4,:) = Qzlow_int(i,j,i4,:)
                    end do
                end if
                if (maxval(Qzlow_ext(i,j,1:nface,mz)) .gt. 0.) then
                    Qzlow_ext(i,j,1:nface,mz) = 0.
                end if

            end do
        end do

    end if

!-----------------------------------------------------------
    if (mpi_R .eq. mpi_nz .and. zlbc .ne. 2) then
        ! Set B.C.'s at top z-boundary.
        do j = 1,ny
            do i = 1,nx

                Qzhigh_ext(i,j,1:nface,:) = Qzhigh_int(i,j,1:nface,:)
                if (minval(Qzhigh_ext(i,j,1:nface,mz)) .lt. 0.) then
                    Qzhigh_ext(i,j,1:nface,mz) = 0
                end if

            end do
        end do
    end if

end subroutine set_bc

!----------------------------------------------------------------------------------------------

subroutine advance_time_level_gl(Q_ri,Q_rp)

    implicit none
    integer i,j,k,ieq,ir
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri, Q_rp
    real Q_xtemp, Q_ytemp, Q_ztemp, dti
    real c1d3, c1d5

    c1d3 = 1./3.
    c1d5 = 1./5.
    dti = 1./dt

    do k = 1,nz
        do j = 1,ny
            do i = 1,nx

                do ieq = 1,nQ
                    do ir=1,nbasis
                        Q_rp(i,j,k,ieq,ir) = Q_ri(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir) + dt*source_r(i,j,k,ieq,ir)
                    end do
                end do

                do ieq = 1,nQ
                    if ( Q_rp(i,j,k,ieq,1) .ne. Q_rp(i,j,k,ieq,1)) then
                        print *,'NaN. Bailing out...','  xc  =',xc(i),'  yc  =',yc(j),'  zc  =',zc(k),'  ieq  =',ieq
                        call exit(-1)
                    endif
                end do

            end do
        end do
    end do

end subroutine advance_time_level_gl

!----------------------------------------------------------------------------------------------

subroutine source_calc(Q_ri,t)

    implicit none
    integer i,j,k,ieq,ipg,ir
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri
    real t,source(npg,nQ),Qin(nQ),dn,dni,Zin,vx,vy,vz,alpha,temp,dne,eta_a,Teev,Tiev,etaJ2,nuei,TemR,Tev,vmax
    real Tcoef,fac,en_floor,gyro,Pres,rh_buf

    source_r(:,:,:,:,:) = 0.0
    source(:,:) = 0.0
    en_floor = P_floor/(aindex - 1.)

    do k = 1,nz
        do j = 1,ny
            do i = 1,nx

                do ipg = 1,npg

                    do ieq = 1,nQ
                        Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,k,ieq,1:nbasis))
                    end do

                    ! Sources for the fluid variables. Nonzero if randomly forced.
                    source(ipg,rh) = 0 !amplen*(ran(iseed) - 0.5)
                    source(ipg,mx) = 0 !amplm*(ran(iseed) - 0.5)
                    source(ipg,my) = 0 !amplm*(ran(iseed) - 0.5)
                    source(ipg,mz) = 0 !amplm*(ran(iseed) - 0.5)
                    source(ipg,en) = 0 !amplen*(ran(iseed) - 0.5)

                    ! Sources for the viscous stress tensor.  The viscous stress is computed using hyperbolic
                    ! relaxation to a parabolic problem. This is a generalization of the problem:
                    !
                    !  partial_t u = partial_x v ,  eps * partial_t v - v = D * partial_x u
                    !
                    ! where eps is a small parameter and D is a diffusion coefficient. In the relaxation limit
                    ! eps -> 0, this becomes partial_t u = D * partial_x^2 u
                    ! The following are the sources for this problem where epsi = 1/eps
                    source(ipg,pxx) = -epsi*Qin(pxx) !+ amplv*(ran(iseed) - 0.5)
                    source(ipg,pyy) = -epsi*Qin(pyy) !+ amplv*(ran(iseed) - 0.5)
                    source(ipg,pzz) = -epsi*Qin(pzz) !+ amplv*(ran(iseed) - 0.5)
                    source(ipg,pxy) = -epsi*Qin(pxy) !+ amplv*(ran(iseed) - 0.5)
                    source(ipg,pxz) = -epsi*Qin(pxz) !+ amplv*(ran(iseed) - 0.5)
                    source(ipg,pyz) = -epsi*Qin(pyz) !+ amplv*(ran(iseed) - 0.5)

                end do

                do ir=1,nbasis
                    do ieq = 1,nQ
                        source_r(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(bval_int_wgt(1:npg,ir)*source(1:npg,ieq))
                    end do
                end do

            end do
        end do
    end do

end subroutine source_calc

!----------------------------------------------------------------------------------------------

subroutine flux_calc_pnts_r(Qpnts_r,fpnts_r,ixyz,npnts)

    ! Calculate the flux "fpnts_r" in direction "ixyz" (x, y, or z) at a set of
    ! points corresponding to conserved quantities "Qpnts_r".

    ! ixyz=1: x-direction
    ! ixyz=2: y-direction
    ! ixyz=3: z-direction

    implicit none
    integer ife, ixyz,npnts
    real, dimension(npnts,nQ):: Qpnts_r, fpnts_r
    real dn,dni,vx,vy,vz,P,asqr,fac,Pre,dnei,Psol,dx2,Tem,smsq,nu,c2d3,c4d3
    real ampx,ampy,ampz,ampd

    nu = epsi*vis
    c2d3 = 2./3.
    c4d3 = 4./3.

    ampx = 0*amplv*(ran(iseed) - 0.5)
    ampy = 0*amplv*(ran(iseed) - 0.5)
    ampz = 0*amplv*(ran(iseed) - 0.5)
    ampd = 0*amplen*(ran(iseed) - 0.5)

    do ife = 1,npnts

        dn = Qpnts_r(ife,rh)
        dni = 1./dn
        smsq = Qpnts_r(ife,mx)**2 + Qpnts_r(ife,my)**2 + Qpnts_r(ife,mz)**2

        P = (aindex - 1.)*(Qpnts_r(ife,en) - 0.5*dni*smsq)
        if (ieos .eq. 2) P = P_1*(dn**7.2 - 1.) + P_base + P
        if (P < P_floor) P = P_floor

        vx = Qpnts_r(ife,mx)*dni
        vy = Qpnts_r(ife,my)*dni
        vz = Qpnts_r(ife,mz)*dni

        if(ixyz .eq. 1) then

            fpnts_r(ife,rh) = Qpnts_r(ife,mx)

            fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vx + P + Qpnts_r(ife,pxx) + ampx
            fpnts_r(ife,my) = Qpnts_r(ife,my)*vx + Qpnts_r(ife,pxy)     + ampy
            fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vx + Qpnts_r(ife,pxz)     + ampz

            fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vx - (Qpnts_r(ife,pxx)*vx + Qpnts_r(ife,pxy)*vy + Qpnts_r(ife,pxz)*vz)

            fpnts_r(ife,pxx) = c4d3*nu*vx
            fpnts_r(ife,pyy) = -c2d3*nu*vx
            fpnts_r(ife,pzz) = -c2d3*nu*vx

            fpnts_r(ife,pxy) = nu*vy
            fpnts_r(ife,pxz) = nu*vz
            fpnts_r(ife,pyz) = 0

        end if
        if(ixyz .eq. 2) then

            fpnts_r(ife,rh) = Qpnts_r(ife,mxa(ixyz))

            fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vy + Qpnts_r(ife,pxy)     + ampx
            fpnts_r(ife,my) = Qpnts_r(ife,my)*vy + P + Qpnts_r(ife,pyy) + ampy
            fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vy + Qpnts_r(ife,pyz)     + ampz

            fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vy - (Qpnts_r(ife,pyy)*vy + Qpnts_r(ife,pxy)*vx + Qpnts_r(ife,pyz)*vz)

            fpnts_r(ife,pxx) = -c2d3*nu*vy
            fpnts_r(ife,pyy) = c4d3*nu*vy
            fpnts_r(ife,pzz) = -c2d3*nu*vy

            fpnts_r(ife,pxy) = nu*vx
            fpnts_r(ife,pxz) = 0
            fpnts_r(ife,pyz) = nu*vz

        end if
        if(ixyz .eq. 3) then

            fpnts_r(ife,rh) = Qpnts_r(ife,mz)

            fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vz + Qpnts_r(ife,pxz)     + ampx
            fpnts_r(ife,my) = Qpnts_r(ife,my)*vz + Qpnts_r(ife,pyz)     + ampy
            fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vz + P + Qpnts_r(ife,pzz) + ampz

            fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vz - (Qpnts_r(ife,pzz)*vz + Qpnts_r(ife,pxz)*vx + Qpnts_r(ife,pyz)*vy)

            fpnts_r(ife,pxx) = -c2d3*nu*vz
            fpnts_r(ife,pyy) = -c2d3*nu*vz
            fpnts_r(ife,pzz) = c4d3*nu*vz

            fpnts_r(ife,pxy) = 0.
            fpnts_r(ife,pxz) = nu*vx
            fpnts_r(ife,pyz) = nu*vy

        end if

    end do

end subroutine

!----------------------------------------------------------------------------------------------

subroutine flux_calc_pnts_r2(Qpnts_r,fpnts_r,ixyz,npnts)
    ! Calculate the flux "fpnts_r" in direction "ixyz" (x, y, or z) at a set of
    ! points corresponding to conserved quantities "Qpnts_r".

    ! ixyz=1: x-direction
    ! ixyz=2: y-direction
    ! ixyz=3: z-direction

    implicit none
    integer ife, ixyz,npnts
    real, dimension(npnts,nQ):: Qpnts_r, fpnts_r
    real dn,dni,vr,P,asqr,fac,Pre,dnei,Psol,dx2,Tem,smsq

    do ife = 1,npnts

        dn = Qpnts_r(ife,rh)
        dni = 1./dn

        smsq = Qpnts_r(ife,mx)**2 + Qpnts_r(ife,my)**2 + Qpnts_r(ife,mz)**2
        vr = Qpnts_r(ife,mxa(ixyz))*dni

        P = (aindex - 1.)*(Qpnts_r(ife,en) - 0.5*dni*smsq)
        if (P < P_floor) P = P_floor

        fpnts_r(ife,rh) = Qpnts_r(ife,mxa(ixyz))
        fpnts_r(ife,mx:mz) = Qpnts_r(ife,mx:mz)*vr
        fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vr

        fpnts_r(ife,mxa(ixyz)) = fpnts_r(ife,mxa(ixyz)) + P

    end do

end subroutine

!----------------------------------------------------------------------------------------------

subroutine flux_cal(Q_r)

    implicit none
    integer i,j,k,ieq,iback,jleft,kdown,i4,i4p,ipnt
    real Qface_x(nfe,nQ), Qface_y(nfe,nQ),Qface_z(nfe,nQ),fface_x(nfe,nQ), fface_y(nfe,nQ), fface_z(nfe,nQ)
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
    real cwavex(nfe),cwavey(nfe),cwavez(nfe),Pre,dni,B2, cfbound, cfx(nface,nQ),cfy(nface,nQ),cfz(nface,nQ)
    real fhllc_x(nface,5),fhllc_y(nface,5),fhllc_z(nface,5),fs(nface,nQ),qvin(nQ)

!--------------------------------------------------------------

    kroe(:) = 1

    do k=1,nz
        do j=1,ny
            do i=1,nx+1

                iback = i-1

                if (i .gt. 1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_x(ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_r(iback,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (i .eq. 1) then
                    do ieq = 1,nQ
                        Qface_x(1:nface,ieq) = Qxlow_ext(j,k,1:nface,ieq)
                    end do
                end if

                if (i .lt. nx+1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_x(ipnt+nface,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if

                if (i .eq. nx+1) then
                    do ieq = 1,nQ
                        Qface_x(nface+1:nfe,ieq) = Qxhigh_ext(j,k,1:nface,ieq)
                    end do
                end if

                call flux_calc_pnts_r(Qface_x,fface_x,1,nfe)

                if(iroe .ne. 1 .and. ihllc .ne. 1) then
                    do i4=1,nfe
                        do ieq=1,nQ
                            qvin(ieq) = Qface_x(i4,ieq)
                        end do
                        cwavex(i4) = cfcal(qvin,1)
                    end do

                    do i4=1,nface
                        cfrx(i4,rh:en) = max(cwavex(i4),cwavex(i4+nface))
                    end do
                end if

                do ieq = 1,nQ
                    do i4=1,nface
                        i4p = i4 + nface
                        flux_x(i4,i,j,k,ieq) = 0.5*(fface_x(i4,ieq) + fface_x(i4p,ieq)) - 0.5*cfrx(i4,ieq)*(Qface_x(i4p,ieq) - Qface_x(i4,ieq))
                    end do
                end do

                kroe(1:nface) = 1

                if (ihllc .eq. 1) call flux_hllc(Qface_x,fface_x,fhllc_x,1)
                if (iroe .eq. 1) call flux_roe(Qface_x,fface_x,fhllc_x,1)

                if (ihllc .eq. 1 .or. iroe .eq. 1) then
                    do ieq = 1,en
                        do i4=1,nface
                            if (kroe(i4) .gt. 0) then
                                flux_x(i4,i,j,k,ieq) = fhllc_x(i4,ieq)
                            end if
                        end do
                    end do
                end if

            end do
        end do
    end do

!----------------------------------------------------

    do k=1,nz
        do j=1,ny+1
            jleft = j-1
            do i=1,nx

                if (j .gt. 1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_y(ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_r(i,jleft,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (j .eq. 1) then
                    do ieq = 1,nQ
                        Qface_y(1:nface,ieq) = Qylow_ext(i,k,1:nface,ieq)
                    end do
                end if

                if (j .lt. ny+1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_y(ipnt+nface,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (j .eq. ny+1) then
                    do ieq = 1,nQ
                        Qface_y(nface+1:nfe,ieq) = Qyhigh_ext(i,k,1:nface,ieq)
                    end do
                end if

                call flux_calc_pnts_r(Qface_y,fface_y,2,nfe)

                if(iroe .ne. 1 .and. ihllc .ne. 1) then
                    do i4=1,nfe
                        do ieq=1,nQ
                            qvin(ieq) = Qface_y(i4,ieq)
                        end do
                        cwavey(i4) = cfcal(qvin,2)
                    end do

                    do i4=1,nface
                        cfry(i4,rh:en) = max(cwavey(i4),cwavey(i4+nface))
                    end do
                end if

                do ieq = 1,nQ
                    do i4=1,nface
                        i4p = i4 + nface
                        flux_y(i4,i,j,k,ieq) = 0.5*(fface_y(i4,ieq) + fface_y(i4p,ieq)) - 0.5*cfry(i4,ieq)*(Qface_y(i4p,ieq) - Qface_y(i4,ieq))
                    end do
                end do

                kroe(1:nface) = 1

                if (ihllc .eq. 1) call flux_hllc(Qface_y,fface_y,fhllc_y,2)
                if (iroe .eq. 1) call flux_roe(Qface_y,fface_y,fhllc_y,2)

                if (ihllc .eq. 1 .or. iroe .eq. 1) then
                    do ieq = 1,en
                        do i4=1,nface
                            if (kroe(i4) .gt. 0) then
                                flux_y(i4,i,j,k,ieq) = fhllc_y(i4,ieq)
                            end if
                        end do
                    end do
                end if

            end do
        end do
    end do

!-------------------------------------------------------

    do k=1,nz+1
        kdown = k-1
        do j=1,ny
            do i=1,nx

                if (k .gt. 1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_z(ipnt,ieq) = sum(bfvals_zp(ipnt,1:nbasis)*Q_r(i,j,kdown,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (k .eq. 1) then
                    do ieq = 1,nQ
                        Qface_z(1:nface,ieq) = Qzlow_ext(i,j,1:nface,ieq)
                    end do
                end if

                if (k .lt. nz+1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_z(ipnt+nface,ieq) = sum(bfvals_zm(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (k .eq. nz+1) then
                    do ieq = 1,nQ
                        Qface_z(nface+1:nfe,ieq) = Qzhigh_ext(i,j,1:nface,ieq)
                    end do
                end if

                call flux_calc_pnts_r(Qface_z,fface_z,3,nfe)

                if(iroe .ne. 1 .and. ihllc .ne. 1) then
                    do i4=1,nfe
                        do ieq=1,nQ
                            qvin(ieq) = Qface_z(i4,ieq)
                        end do
                        cwavez(i4) = cfcal(qvin,3)
                    end do

                    do i4=1,nface
                        cfrz(i4,rh:en) = max(cwavez(i4),cwavez(i4+nface))
                    end do
                end if

                do ieq = 1,nQ
                    do i4=1,nface
                        i4p = i4 + nface
                        flux_z(i4,i,j,k,ieq) = 0.5*(fface_z(i4,ieq) + fface_z(i4p,ieq)) - 0.5*cfrz(i4,ieq)*(Qface_z(i4p,ieq) - Qface_z(i4,ieq))
                    end do
                end do

                kroe(1:nface) = 1

                if (ihllc .eq. 1) call flux_hllc(Qface_z,fface_z,fhllc_z,3)
                if (iroe .eq. 1) call flux_roe(Qface_z,fface_z,fhllc_z,3)

                if (ihllc .eq. 1 .or. iroe .eq. 1) then
                    do ieq = 1,en
                        do i4=1,nface
                            if (kroe(i4) .gt. 0) then
                                flux_z(i4,i,j,k,ieq) = fhllc_z(i4,ieq)
                            end if
                        end do
                    end do
                end if

            end do
        end do
    end do

end subroutine flux_cal

!----------------------------------------------------------------------------------------------

subroutine flux_hllc(Qlr,flr,fhllc,ixyz)

!    Compute ion or electron fluxes (density, momentum, energy) using
!    nonrelativistic HD HLLC approximate Riemann solver developed by Batten, 1997,
!    *****    "On the Choice of Wavespeeds for the HLLC Riemann Solver"

!    This takes into account the left and right propagating shocks along with
!    the contact (tangential) discontinuity.

    implicit none
    real Qlr(nfe,nQ),flr(nfe,nQ),fhllc(nface,5)
    real sm_num(nface),sm_den(nface),qtilde(nface,5),rtrho(nfe),rtrho_i(nface),qsq(nfe)
    real s_lr(nfe),ctilde(nface),hlr(nfe),cslr(nfe),ctsq(nface),Zi, mfact, csfac
    real aq(nfe),bq(nfe),Qstar(nfe,6),fstar(nfe,6),pstar(nface),s_m(nface)
    real rhov(nfe),vlr(nfe,3),plr(nfe),slsm_i(nface),rho_i,qslr(nfe),sq_lr(nfe),slrm_i(nfe),B2(nfe),cf(nfe)
    integer ixyz,i4,i4p1,nr,jie,k,k2,ieq,iparr,iperp1,iperp2,ibatten
    integer rhj,mxj,myj,mzj,enj,psj,ivar(5),ipassive,nhll,ib1,ib2

    nhll = 5
    rhj = rh
    mxj = mx
    myj = my
    mzj = mz
    enj = en

    iparr  = mxa(ixyz)
    iperp1 = mya(ixyz)
    iperp2 = mza(ixyz)

    ivar(1) = rhj
    ivar(2) = iparr
    ivar(3) = iperp1
    ivar(4) = iperp2
    ivar(5) = enj

    ! Indices 1:4 denote the left state, while indices 5:8 denote the right state.
    ibatten = 0

    do k=1,nfe
        rhov(k) = Qlr(k,rhj)
        rho_i = 1.0/rhov(k)
        vlr(k,1) = Qlr(k,iparr)*rho_i        ! velocity parallel to direction of flux computation
        vlr(k,2) = Qlr(k,iperp1)*rho_i        ! velocity in perpendicular direction 1
        vlr(k,3) = Qlr(k,iperp2)*rho_i        ! velocity in perpendicular direction 2
        qsq(k) = vlr(k,1)**2 + vlr(k,2)**2 + vlr(k,3)**2
        plr(k) = (aindex - 1.0)*(Qlr(k,enj) - 0.5*rhov(k)*qsq(k))        ! pressure
        if(ieos .eq. 2) plr(k) = P_1*(rhov(k)**7.2 - 1.) + P_base + plr(k)
        rtrho(k) = sqrt(rhov(k))
    end do

    do k=1,nface
        k2 = k + nface
        if(ieos .eq. 2)then
            cslr(k) = vlr(k,1) - sqrt(7.2*P_1*rhov(k)**6.2 + plr(k)*rho_i)       ! lambda_M(Q_l)
            cslr(k2) = vlr(k2,1) + sqrt(7.2*P_1*rhov(k2)**6.2 + plr(k2)*rho_i)       ! lambda_P(Q_r)
        else
            cslr(k) = vlr(k,1) - sqrt(aindex*plr(k)/rhov(k))       ! lambda_M(Q_l)
            cslr(k2) = vlr(k2,1) + sqrt(aindex*plr(k2)/rhov(k2) )       ! lambda_P(Q_r)
        end if
    end do

    if (ibatten .eq. 1) then  ! compute wave speeds using Roe averages following Batten, 1997

        do k=1,nface
            k2 = k + nface
            rtrho_i(k) = 1.0/(rtrho(k) + rtrho(k2))
            qtilde(k,1) = (rtrho(k)*vlr(k,1) + rtrho(k2)*vlr(k2,1))*rtrho_i(k)
            qtilde(k,2) = (rtrho(k)*vlr(k,2) + rtrho(k2)*vlr(k2,2))*rtrho_i(k)
            qtilde(k,3) = (rtrho(k)*vlr(k,3) + rtrho(k2)*vlr(k2,3))*rtrho_i(k)
            qsq(k) = qtilde(k,1)**2 + qtilde(k,2)**2 + qtilde(k,3)**2
            hlr(k) = (Qlr(k,enj) + plr(k))/rhov(k)
            hlr(k2) = (Qlr(k2,enj) + plr(k2))/rhov(k2)
            qtilde(k,4) = (rtrho(k)*hlr(k) + rtrho(k2)*hlr(k2))*rtrho_i(k)
            ctsq(k) = (aindex - 1.0)*(qtilde(k,4) - 0.5*qsq(k))
        end do
        if (minval(ctsq) .ge. 0.0) then
            ctilde = sqrt(ctsq)
            qslr(1:nface) = qtilde(1:nface,1) - ctilde(1:nface)       ! lambda_M(Q_Roe)
            qslr(nface+1:nfe) = qtilde(nface+1:nfe,1) + ctilde(nface+1:nfe)    ! lambda_P(Q_Roe)
        end if
        if (minval(ctsq) .lt. 0.0) then
            ibatten = 0
        end if

    end if

    if (ibatten .eq. 0) then
        do k=1,nface
            k2 = k + nface
            if(ieos .eq. 2)then
                qslr(k) = vlr(k2,1) - sqrt(7.2*P_1*rhov(k2)**6.2 + plr(k2)*rho_i)       ! lambda_M(Q_r)
                qslr(k2) = vlr(k,1) + sqrt(7.2*P_1*rhov(k)**6.2 + plr(k)*rho_i)       ! lambda_P(Q_l)
            else
                qslr(k) = vlr(k2,1) - sqrt(aindex*plr(k2)/rhov(k2))       ! lambda_M(Q_r)
                qslr(k2) = vlr(k,1) + sqrt(aindex*plr(k)/rhov(k))       ! lambda_P(Q_l)
                end if
        end do
    end if

    ! Calculate the slow and fast wavespeeds S_L, S_R of the left, right propagating shocks.
    do k=1,nface
        k2 = k + nface
        s_lr(k) = min(cslr(k),qslr(k))         ! S_L = min(lambda_M(Q_l),lambda_M(Q_r or Q_Roe))
        s_lr(k2) = max(cslr(k2),qslr(k2))         ! S_R = max(lambda_P(Q_r),lambda_P(Q_l or Q_Roe))
        sm_num(k) = rhov(k2)*vlr(k2,1)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*vlr(k,1)*(s_lr(k) - vlr(k,1))
        sm_num(k) = sm_num(k) + plr(k) - plr(k2)
        sm_den(k) = rhov(k2)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*(s_lr(k) - vlr(k,1))
    end do
    where (sm_den .eq. 0.0) sm_den = rh_floor

    ! Calculate the wavespeed S_M of the contact discontinuity.

    do k=1,nface
        s_m(k) = sm_num(k)/sm_den(k)                                          ! Eq. (34) of Batten, 1997
        pstar(k) = rhov(k)*(vlr(k,1) - s_lr(k))*(vlr(k,1) - s_m(k)) + plr(k)  ! Eq. (36) of Batten, 1997
    end do


    ! Now, calculate Q_l* and Q_r* in order to calculate F_l* and F_r*.
    do k=1,nfe
        if (k .le. nface) i4 = k
        if (k .gt. nface) i4 = k - nface
        sm_den(1) = s_lr(k) - s_m(i4)        ! S_{L,R} - S_M

        if (sm_den(1) .eq. 0.0) then
            sm_den(1) = rh_floor
        end if

        slrm_i(k) = 1.0/sm_den(1)
        sq_lr(k) = s_lr(k) - vlr(k,1)            ! S_{L,R} - q_{l,r}
        Qstar(k,1) = rhov(k)*sq_lr(k)*slrm_i(k)                              ! Eq. (35) of Batten ,1997
        Qstar(k,2) = (sq_lr(k)*Qlr(k,iparr) + pstar(i4) - plr(k))*slrm_i(k)  ! Eq. (37-39) of Batten, 1997
        Qstar(k,3) = sq_lr(k)*Qlr(k,iperp1)*slrm_i(k)                        ! Eq. (37-39) of Batten, 1997
        Qstar(k,4) = sq_lr(k)*Qlr(k,iperp2)*slrm_i(k)                        ! Eq. (37-39) of Batten, 1997
        Qstar(k,5) = (sq_lr(k)*Qlr(k,enj) - plr(k)*vlr(k,1) + pstar(i4)*s_m(i4))*slrm_i(k)
        ! Eq. (40) of Batten, 1997

        do ieq=1,nhll
            fstar(k,ieq) = flr(k,ivar(ieq)) + s_lr(k)*(Qstar(k,ieq) - Qlr(k,ivar(ieq)))  ! Eq. (29) of Batten, 1997
        end do
    end do

    ! Finally, calculate the HLLC fluxes from F_l, F_l*, F_r*, and F_r...
    do i4 = 1,nface                        ! Use Eq. (26) of Batten ,1997
        if (s_lr(i4) .gt. 0.0) then                               ! if S_L > 0
            do ieq=1,nhll
                fhllc(i4,ivar(ieq)) = flr(i4,ivar(ieq))        ! F_HLLC = F_l
            end do
        end if
        if (s_lr(i4) .le. 0.0 .and. 0.0 .lt. s_m(i4)) then        ! if S_L <= 0 < S_M
            do ieq=1,nhll
                fhllc(i4,ivar(ieq)) = fstar(i4,ieq)            ! F_HLLC = F_l*
            end do
        end if
        if (s_m(i4) .le. 0.0 .and. 0.0 .le. s_lr(i4+nface)) then  ! if S_M <= 0 <= S_R
            do ieq=1,nhll
                fhllc(i4,ivar(ieq)) = fstar(i4+nface,ieq)      ! F_HLLC = F_r*
            end do
        end if
        if (s_lr(i4+nface) .lt. 0.0) then                         ! if S_R < 0
            do ieq=1,nhll
                fhllc(i4,ivar(ieq)) = flr(i4+nface,ivar(ieq))  ! F_HLLC = F_r
            end do
        end if

    end do

end subroutine flux_hllc

    !----------------------------------------------------------------

subroutine flux_roe(Qlr,flr,froef,ixyz)

    implicit none
    integer ixyz,i9,j9,ilr,iparr,iperp
    real Qlr(nfe,nQ),flr(nfe,nQ),froef(nface,5)
    real dflux(nface,5),dwr(nface,5),dsk(nface),Qri
    real evec(5,5),swr(2,nface,5),skr(2,nface),asq,lam9(5),a9(5),sk9
    real a9e,ea2i,vsq,Pri,vels(3),dnii,en_floor
    real vc(3),hc,dele(nface),delrho(nface),delmom1(nface),delmom2(nface),delmom3(nface),sc9,vcsq,vcdel,csq,gmi,bsq


    ! approximate Roe solver for non-relativistic two-fluid, from Eulderink and Mellema, 1994

    ! Starting from HD equations for a given fluid (either electrons or ions), this routine
    ! computes "dflux_i = (cQ)_i+1 - (cQ)_i" terms by characteristic wave decomposition using
    ! a Roe matrix.

    ! INPUT is Qr(,1)=conserved density, Qr(,2)=energy, Qr(,3:5)=momentum
    !       Qpr(,1)=density, Qpr(,2)=pressure, Qpr(,3:5)=velocity
    !       n = number of spatial points in given direction (x or y)

    !        ixyz = 1: flux is being computed in x-direction
    !        ixyz = 2: flux is being computed in y-direction


    ! OUTPUT is  "dflux_i = (cQ)_i+1 - (cQ)_i"
    !    dflux(,1)=density term, dflux(,2)=energy term, dflux(,3:5)=momentum terms



    ! evec(,1) through evec(,5) are the eigenvectors of the Roe matrix.

    en_floor = P_floor*aindm1
    gmi = 1./(aindex - 1.)

    evec(:,:) = 0.
    evec(1,1) = 1.
    evec(1,2) = 1.
    evec(1,3) = 1.
    evec(4,4) = 1.
    evec(5,5) = 1.


    do i9=1,nface
        kroe(i9) = 1
        do j9=1,2
            if (j9 .eq. 1) ilr = i9
            if (j9 .eq. 2) ilr = i9 + nface

            ! Make sure density, pressure, and energy are at least at floor values.

            dnii = 1./Qlr(ilr,rh)
            vels(1) = Qlr(ilr,mxa(ixyz))*dnii
            vels(2) = Qlr(ilr,mya(ixyz))*dnii
            vels(3) = Qlr(ilr,mza(ixyz))*dnii
            vsq = vels(1)**2 + vels(2)**2 + vels(3)**2
            Pri = aindm1*(Qlr(ilr,en) - 0.5*Qlr(ilr,rh)*vsq)

            skr(j9,i9) = sqrt(Qlr(ilr,rh))
            swr(j9,i9,1) = skr(j9,i9)*vels(1)  ! sqrt(rho) * v_x
            swr(j9,i9,2) = skr(j9,i9)*vels(2)  ! sqrt(rho) * v_y
            swr(j9,i9,3) = skr(j9,i9)*vels(3)  ! sqrt(rho) * v_z
            swr(j9,i9,4) = 0.5*skr(j9,i9)*(vsq + cp*Pri/Qlr(ilr,rh))
        end do
    end do

    do i9=1,nface
        Qri = 1./Qlr(i9,rh)    ! Increments in conserved quantities are normalized w.r.t. density.

        delrho(i9) = Qlr(i9+nface,rh)*Qri - 1.    ! delrho = increment in conserved density
        dele(i9) = (Qlr(i9+nface,en) - Qlr(i9,en))*Qri    ! *one_mime(jie)
        ! dele = increment in conserved energy

        delmom1(i9) = (Qlr(i9+nface,mxa(ixyz)) - Qlr(i9,mxa(ixyz)))*Qri    ! del1 = increment in x-momentum
        delmom2(i9) = (Qlr(i9+nface,mya(ixyz)) - Qlr(i9,mya(ixyz)))*Qri    ! del2 = increment in y-momentum
        delmom3(i9) = (Qlr(i9+nface,mza(ixyz)) - Qlr(i9,mza(ixyz)))*Qri    ! del3 = increment in z-momentum
        ! dwr(i,1:3) = 0.5*[sqrt(rho_i)v_{1:3,i} + sqrt(rho_{i+1})v_{1:3,i+1}]
        ! dwr(i,4) = 0.5*[sqrt(rho_i) enthalpy(i)/rho(i) + sqrt(rho_{i+1}) enthalpy(i+1)/rho(i+1)]

        do j9=1,4
            dwr(i9,j9) = 0.5*(swr(2,i9,j9) + swr(1,i9,j9))
        enddo
        dsk(i9) = 0.5*(skr(2,i9) + skr(1,i9))

    enddo

    ! The Roe average of a quantity is the arithmetic average
    ! between neighboring cells weighted by the square root of density.
    ! For example, for "v_x" at position "i" the Roe average "v_{cx}" is

    !    v_{cx} = [sqrt(rho_i)v_{xi} + sqrt(rho_{i+1}) v_{x,i+1}]/[sqrt(rho_i) + sqrt(rho_{i+1})]

    do i9=1,nface
        vc(1) = dwr(i9,1)/dsk(i9)    ! component 1 of Roe-averaged velocity (x if jie=1, y if jie=2)
        vc(2) = dwr(i9,2)/dsk(i9)    ! component 2 of Roe-averaged velocity (y if jie=1, x if jie=2)
        vc(3) = dwr(i9,3)/dsk(i9)    ! component 3 of Roe-averaged velocity (z-component)
        hc = dwr(i9,4)/dsk(i9)    ! Roe-averaged enthalpy/density

        vcsq = vc(1)*vc(1) + vc(2)*vc(2) + vc(3)*vc(3)
        asq = aindm1*(hc - 0.5*vcsq)    ! squared sound speed
        if (asq .le. 0.0) then
            kroe(i9) = 0    ! asq = (aindex - 1.0d0)*hc
        end if
        sc9 = sqrt(asq)    ! sound speed


        ! Define the characteristic speeds (eigenvalues of the Roe matrix).
        lam9(1) = abs(vc(1) - sc9)
        lam9(2) = abs(vc(1) + sc9)
        lam9(3) = abs(vc(1))
        lam9(4) = lam9(3)
        lam9(5) = lam9(3)

        ! Define the eigenvectors evec(,1)...evec(,5) of the Roe matrix.
        evec(2,1) = hc - sc9*vc(1)
        evec(3,1) = vc(1) - sc9
        evec(4,1) = vc(2)
        evec(5,1) = vc(3)

        evec(2,2) = hc + sc9*vc(1)
        evec(3,2) = vc(1) + sc9
        evec(4,2) = vc(2)
        evec(5,2) = vc(3)

        evec(2,3) = 0.5*vcsq
        evec(3,3) = vc(1)
        evec(4,3) = vc(2)
        evec(5,3) = vc(3)
        evec(2,4) = vc(2)
        evec(2,5) = vc(3)

        ! Define a few intermediate variables needed for computing the expansion coefficients.

        ea2i = aindm1/asq
        vcdel = vc(1)*delmom1(i9) + vc(2)*delmom2(i9) + vc(3)*delmom3(i9) - dele(i9)
        a9e = 0.5*ea2i*(0.5*vcsq*delrho(i9) - vcdel)
        sk9 = 0.5*(delmom1(i9) - vc(1)*delrho(i9))/sc9

        ! Define the expansion coefficients a9_1...a9_5 such that
        !     Delta Q = {delrho,dele,del1,del2,del3} = sum [a9_j evec(,j)]

        a9(1) = a9e - sk9
        a9(2) = a9e + sk9
        a9(3) = ea2i*((hc - vcsq)*delrho(i9) + vcdel)
        a9(4) = delmom2(i9) - vc(2)*delrho(i9)
        a9(5) = delmom3(i9) - vc(3)*delrho(i9)

        ! The flux increments "dflux" are now given by   Delta F = sum [a9_j lam_j evec(,j)]

        dflux(i9,1:5) = 0.
        do j9=1,5
            dflux(i9,1) = dflux(i9,1) + a9(j9)*lam9(j9)*evec(1,j9)
            dflux(i9,2) = dflux(i9,2) + a9(j9)*lam9(j9)*evec(2,j9)
            dflux(i9,3) = dflux(i9,3) + a9(j9)*lam9(j9)*evec(3,j9)
            dflux(i9,4) = dflux(i9,4) + a9(j9)*lam9(j9)*evec(4,j9)
            dflux(i9,5) = dflux(i9,5) + a9(j9)*lam9(j9)*evec(5,j9)
        enddo
        dflux(i9,1) = dflux(i9,1)*Qlr(i9,rh)     ! flux increment in density
        dflux(i9,2) = dflux(i9,2)*Qlr(i9,rh)     ! flux increment in energy
        dflux(i9,3) = dflux(i9,3)*Qlr(i9,rh)     ! flux increment in parallel momentum
        dflux(i9,4) = dflux(i9,4)*Qlr(i9,rh)     ! flux increment in perpendicular momentum
        dflux(i9,5) = dflux(i9,5)*Qlr(i9,rh)     ! flux increment in z-momentum

        if (kroe(i9) .gt. 0) then
            froef(i9,rh) = 0.5*(flr(i9,rh) + flr(i9+nface,rh) - dflux(i9,1))
            froef(i9,en) = 0.5*(flr(i9,en) + flr(i9+nface,en) - dflux(i9,2))
            froef(i9,mxa(ixyz)) = 0.5*(flr(i9,mxa(ixyz)) + flr(i9+nface,mxa(ixyz)) - dflux(i9,3))
            froef(i9,mya(ixyz)) = 0.5*(flr(i9,mya(ixyz)) + flr(i9+nface,mya(ixyz)) - dflux(i9,4))
            froef(i9,mza(ixyz)) = 0.5*(flr(i9,mza(ixyz)) + flr(i9+nface,mza(ixyz)) - dflux(i9,5))
        end if

    end do

end subroutine flux_roe

!----------------------------------------------------------------------------------------------

subroutine innerintegral(Q_r)

    implicit none
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
    integer i,j,k,ieq,ipg,ir
    real Qinner(npg,nQ),finner_x(npg,nQ), finner_y(npg,nQ), finner_z(npg,nQ), int_r(nbastot,nQ)

    integral_r(:,:,:,:,:) = 0.

    do k = 1,nz
    do j = 1,ny
    do i = 1,nx

        do ieq = 1,nQ
            do ipg = 1,npg
                Qinner(ipg,ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
            end do
        end do

        call flux_calc_pnts_r(Qinner,finner_x,1,npg)
        call flux_calc_pnts_r(Qinner,finner_y,2,npg)
        call flux_calc_pnts_r(Qinner,finner_z,3,npg)

        do ieq = 1,nQ

            int_r(kx,ieq) = 0.25*cbasis(kx)*dxi*sum(wgt3d(1:npg)*finner_x(1:npg,ieq))
            int_r(ky,ieq) = 0.25*cbasis(ky)*dyi*sum(wgt3d(1:npg)*finner_y(1:npg,ieq))
            int_r(kz,ieq) = 0.25*cbasis(kz)*dzi*sum(wgt3d(1:npg)*finner_z(1:npg,ieq))

            if (nbasis .gt. 4) then
                int_r(kyz,ieq)  = 0.25*cbasis(kyz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kz)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kyz)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,ky)*finner_z(1:npg,ieq))
                int_r(kzx,ieq)  = 0.25*cbasis(kzx)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kz)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kzx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kx)*finner_z(1:npg,ieq))
                int_r(kxy,ieq)  = 0.25*cbasis(kxy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,ky)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kxy)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kx)*finner_y(1:npg,ieq))
                int_r(kxyz,ieq) = 0.25*cbasis(kxyz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyz)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kxyz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzx)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kxyz)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxy)*finner_z(1:npg,ieq))
            end if

            if (nbasis .gt. 8) then
                int_r(kyzz,ieq) = 0.25*cbasis(kyzz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzz)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyz)*finner_z(1:npg,ieq))
                int_r(kzxx,ieq) = 0.25*cbasis(kzxx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxx)*finner_z(1:npg,ieq)) &
                  + 0.25*cbasis(kzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzx)*finner_x(1:npg,ieq))
                int_r(kxyy,ieq) = 0.25*cbasis(kxyy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyy)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxy)*finner_y(1:npg,ieq))
                int_r(kyyz,ieq) = 0.25*cbasis(kyyz)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyy)*finner_z(1:npg,ieq)) &
                  + 0.25*cbasis(kyyz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyz)*finner_y(1:npg,ieq))
                int_r(kzzx,ieq) = 0.25*cbasis(kzzx)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzz)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kzzx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzx)*finner_z(1:npg,ieq))
                int_r(kxxy,ieq) = 0.25*cbasis(kxxy)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxx)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kxxy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxy)*finner_x(1:npg,ieq))
                int_r(kyyzz,ieq) = 0.25*cbasis(kyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzz)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyyz)*finner_z(1:npg,ieq))
                int_r(kzzxx,ieq) = 0.25*cbasis(kzzxx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxx)*finner_z(1:npg,ieq)) &
                  + 0.25*cbasis(kzzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzzx)*finner_x(1:npg,ieq))
                int_r(kxxyy,ieq) = 0.25*cbasis(kxxyy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyy)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kxxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxxy)*finner_y(1:npg,ieq))
                int_r(kyzxx,ieq) = 0.25*cbasis(kyzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kyzxx)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzxx)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kyzxx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxxy)*finner_z(1:npg,ieq))
                int_r(kzxyy,ieq) = 0.25*cbasis(kzxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kzxyy)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxyy)*finner_z(1:npg,ieq)) &
                  + 0.25*cbasis(kzxyy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyyz)*finner_x(1:npg,ieq))
                int_r(kxyzz,ieq) = 0.25*cbasis(kxyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_z(1:npg,ieq)) &
                  + 0.25*cbasis(kxyzz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyzz)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kxyzz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzzx)*finner_y(1:npg,ieq))
                int_r(kxyyzz,ieq) = 0.25*cbasis(kxyyzz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyyzz)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kxyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyzz)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kxyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxyy)*finner_z(1:npg,ieq))
                int_r(kyzzxx,ieq) = 0.25*cbasis(kyzzxx)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzzxx)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kyzzxx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzxx)*finner_z(1:npg,ieq)) &
                  + 0.25*cbasis(kyzzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyzz)*finner_x(1:npg,ieq))
                int_r(kzxxyy,ieq) = 0.25*cbasis(kzxxyy)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxxyy)*finner_z(1:npg,ieq)) &
                  + 0.25*cbasis(kzxxyy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxyy)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kzxxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzxx)*finner_y(1:npg,ieq))
                int_r(kxxyyzz,ieq) = 0.25*cbasis(kxxyyzz)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyyzz)*finner_x(1:npg,ieq)) &
                  + 0.25*cbasis(kxxyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzzxx)*finner_y(1:npg,ieq)) &
                  + 0.25*cbasis(kxxyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxxyy)*finner_z(1:npg,ieq))
            end if
        end do

        do ieq = 1,nQ
            do ir=1,nbasis
                integral_r(i,j,k,ieq,ir) = int_r(ir,ieq)
            end do
        end do

    end do
    end do
    end do

end subroutine innerintegral


!----------------------------------------------------------------------------------------------

subroutine glflux

    implicit none
    integer i,j,k,ieq,ir,iqfa

    do ieq = 1,nQ
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
    end do

end subroutine glflux

!----------------------------------------------------------------------------------------------

subroutine glflux2

    implicit none
    integer i,j,k,ieq,ir

    do ieq = 1,nQ
    do k = 1,nz
    do j = 1,ny

        do i = 1,nx

            glflux_r(i,j,k,ieq,1) = 0.25*(dxi*(wgt2d(1)*(flux_x(1,i+1,j,k,ieq) - flux_x(1,i,j,k,ieq)))  &
                                        + dyi*(wgt2d(1)*(flux_y(1,i,j+1,k,ieq) - flux_y(1,i,j,k,ieq)))  &
                                        + dzi*(wgt2d(1)*(flux_z(1,i,j,k+1,ieq) - flux_z(1,i,j,k,ieq)))  &
                                        + dxi*(wgt2d(2)*(flux_x(2,i+1,j,k,ieq) - flux_x(2,i,j,k,ieq)))  &
                                        + dyi*(wgt2d(2)*(flux_y(2,i,j+1,k,ieq) - flux_y(2,i,j,k,ieq)))  &
                                        + dzi*(wgt2d(2)*(flux_z(2,i,j,k+1,ieq) - flux_z(2,i,j,k,ieq)))  &
                                        + dxi*(wgt2d(3)*(flux_x(3,i+1,j,k,ieq) - flux_x(3,i,j,k,ieq)))  &
                                        + dyi*(wgt2d(3)*(flux_y(3,i,j+1,k,ieq) - flux_y(3,i,j,k,ieq)))  &
                                        + dzi*(wgt2d(3)*(flux_z(3,i,j,k+1,ieq) - flux_z(3,i,j,k,ieq)))  &
                                        + dxi*(wgt2d(4)*(flux_x(4,i+1,j,k,ieq) - flux_x(4,i,j,k,ieq)))  &
                                        + dyi*(wgt2d(4)*(flux_y(4,i,j+1,k,ieq) - flux_y(4,i,j,k,ieq)))  &
                                        + dzi*(wgt2d(4)*(flux_z(4,i,j,k+1,ieq) - flux_z(4,i,j,k,ieq))))
        end do

        do ir=2,nbasis
            do i = 1,nx

                glflux_r(i,j,k,ieq,ir) = wgtbf_xmp(1,2,ir)*flux_x(1,i+1,j,k,ieq) + wgtbf_xmp(1,1,ir)*flux_x(1,i,j,k,ieq)  &
                               + wgtbf_ymp(1,2,ir)*flux_y(1,i,j+1,k,ieq) + wgtbf_ymp(1,1,ir)*flux_y(1,i,j,k,ieq)  &
                               + wgtbf_zmp(1,2,ir)*flux_z(1,i,j,k+1,ieq) + wgtbf_zmp(1,1,ir)*flux_z(1,i,j,k,ieq)  &
                               + wgtbf_xmp(2,2,ir)*flux_x(2,i+1,j,k,ieq) + wgtbf_xmp(2,1,ir)*flux_x(2,i,j,k,ieq)  &
                               + wgtbf_ymp(2,2,ir)*flux_y(2,i,j+1,k,ieq) + wgtbf_ymp(2,1,ir)*flux_y(2,i,j,k,ieq)  &
                               + wgtbf_zmp(2,2,ir)*flux_z(2,i,j,k+1,ieq) + wgtbf_zmp(2,1,ir)*flux_z(2,i,j,k,ieq)  &
                               + wgtbf_xmp(3,2,ir)*flux_x(3,i+1,j,k,ieq) + wgtbf_xmp(3,1,ir)*flux_x(3,i,j,k,ieq)  &
                               + wgtbf_ymp(3,2,ir)*flux_y(3,i,j+1,k,ieq) + wgtbf_ymp(3,1,ir)*flux_y(3,i,j,k,ieq)  &
                               + wgtbf_zmp(3,2,ir)*flux_z(3,i,j,k+1,ieq) + wgtbf_zmp(3,1,ir)*flux_z(3,i,j,k,ieq)  &
                               + wgtbf_xmp(4,2,ir)*flux_x(4,i+1,j,k,ieq) + wgtbf_xmp(4,1,ir)*flux_x(4,i,j,k,ieq)  &
                               + wgtbf_ymp(4,2,ir)*flux_y(4,i,j+1,k,ieq) + wgtbf_ymp(4,1,ir)*flux_y(4,i,j,k,ieq)  &
                               + wgtbf_zmp(4,2,ir)*flux_z(4,i,j,k+1,ieq) + wgtbf_zmp(4,1,ir)*flux_z(4,i,j,k,ieq)  &
                                   - integral_r(i,j,k,ieq,ir)

            end do
        end do

    end do
    end do
    end do

end subroutine glflux2

!-----------------------------------------------------------!
!*******************calculate freezing speeds***************!
!-----------------------------------------------------------!
real function cfcal(Qcf,cases)

    implicit none
    integer cases
    real Qcf(nQ)
    real Pi, Pe, P, B2, ne, cs
    real dn,dni,vx,vy,vz,hx,hy,hz,dnei,va2,vf1,lil02,va,fac

    dn = Qcf(rh)
    dni = 1./dn
    vx = Qcf(mx)*dni
    vy = Qcf(my)*dni
    vz = Qcf(mz)*dni

    if (ieos .eq. 1) then
        P = (aindex - 1.)*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
        cs = sqrt(aindex*P*dni)
    end if

    if (ieos .eq. 2) then
        cs = sqrt(7.2*P_1*dn**6.2)
    end if

    select case (cases)
        case (1) !freezing speed in x direction for fluid variable
            cfcal = abs(vx) + cs

        case (2) !freezing speed in y direction for fluid variable
            cfcal = abs(vy) + cs

        case (3) !freezing speed in z direction for fluid variable
            cfcal = abs(vz) + cs
    end select

end function cfcal

!----------------------------------------------------------------------------------

subroutine limiter(Q_r)

    implicit none
    integer i, j, k, ieq, ipge, minindex, ir
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
    real Qedge(npge,nQ),theta,Qmin(nQ), deltaQ(nQ)
    real eps, Qrhmin, QPmin, P(npge), Pave, dn, dni, epsiP, thetaj
    real*8 a, b, c

    eps = rh_min

    do k = 1,nz
    do j = 1,ny
    do i = 1,nx

        if (Q_r(i,j,k,rh,1) < eps) then

            do ir=2,nbasis
                Q_r(i,j,k,rh:en,ir) = 0.0
            end do
            Q_r(i,j,k,rh,1) = eps

        else

            do ipge = 1,npge
                Qedge(ipge,rh) = sum(bf_faces(ipge,1:nbasis)*Q_r(i,j,k,rh,1:nbasis))
            end do

            Qrhmin = minval(Qedge(:,rh))

            if (Qrhmin < eps) then
                theta = (eps - Q_r(i,j,k,rh,1))/(Qrhmin - Q_r(i,j,k,rh,1))
                if (theta .gt. 1.) then
                    theta = 1.
                end if

                if (theta .lt. 0) then
                    theta = 0.
                end if
                do ir=2,nbasis
                    Q_r(i,j,k,rh,ir) = abs(theta)*Q_r(i,j,k,rh,ir)
                end do
            end if

        end if

    end do
    end do
    end do

end subroutine limiter

!----------------------------------------------------------------------------------------------

subroutine prepare_exchange(Q_r)

    integer ieq, i, j, k, ipnt
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r

    do ieq = 1,nQ
    do j = 1,ny
    do i = 1,nx
        do ipnt=1,nface
            Qzlow_int(i,j,ipnt,ieq) = sum(bfvals_zm(ipnt,1:nbasis)*Q_r(i,j,1,ieq,1:nbasis))
            Qzhigh_int(i,j,ipnt,ieq) = sum(bfvals_zp(ipnt,1:nbasis)*Q_r(i,j,nz,ieq,1:nbasis))
        end do
    end do
    end do
    end do

    do ieq = 1,nQ
    do k = 1,nz
    do i = 1,nx
        do ipnt=1,nface
            Qylow_int(i,k,ipnt,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_r(i,1,k,ieq,1:nbasis))
            Qyhigh_int(i,k,ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_r(i,ny,k,ieq,1:nbasis))
        end do
    end do
    end do
    end do

    do ieq = 1,nQ
    do k = 1,nz
    do j = 1,ny
        do ipnt=1,nface
            Qxlow_int(j,k,ipnt,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_r(1,j,k,ieq,1:nbasis))
            Qxhigh_int(j,k,ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_r(nx,j,k,ieq,1:nbasis))
        end do
    end do
    end do
    end do
    call exchange_flux

end subroutine

!----------------------------------------------------------------------------------------------


subroutine exchange_flux

    integer mpi_size

    call MPI_BARRIER(cartcomm,ierr)

     mpi_size=ny*nz*nface*nQ

    if (nbrs(EAST) .ne. MPI_PROC_NULL) then
        call MPI_ISend(Qxhigh_int,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(1),ierr)
    endif

    if (nbrs(WEST) .ne. MPI_PROC_NULL) then
        call MPI_IRecv(Qxlow_ext,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(2),ierr)
    endif

    if (nbrs(EAST) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(Qxhigh_ext,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(3),ierr)
    endif

    if (nbrs(WEST) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(2),stats(:,2),ierr)
        call MPI_ISend(Qxlow_int,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(4),ierr)
    endif

    if (nbrs(EAST) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(3),stats(:,3),ierr)
    endif

    if (nbrs(WEST) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(4),stats(:,4),ierr)
    endif

    if (mpi_P .eq. 1 .and. xlbc .eq. 0) then
        Qxlow_ext = Qxlow_int
    end if

    if (mpi_P .eq. mpi_nx .and. xhbc .eq. 0) then
        Qxhigh_ext = Qxhigh_int
    end if

    mpi_size=nface*nx*nz*nQ

    if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
        call MPI_ISend(Qyhigh_int,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(1),ierr)
    endif

    if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
        call MPI_IRecv(Qylow_ext,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(2),ierr)
    endif

    if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(Qyhigh_ext,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(3),ierr)
    endif

    if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(2),stats(:,2),ierr)
        call MPI_ISend(Qylow_int,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(4),ierr)
    endif

    if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(3),stats(:,3),ierr)
    endif

    if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(4),stats(:,4),ierr)
    endif

    if (mpi_Q .eq. 1 .and. ylbc .eq. 0) then
        Qylow_ext = Qylow_int
    end if

    if (mpi_Q .eq. mpi_ny .and. yhbc .eq. 0) then
        Qyhigh_ext = Qyhigh_int
    end if

    mpi_size=nface*nx*ny*nQ

    if (nbrs(UP) .ne. MPI_PROC_NULL) then
        call MPI_ISend(Qzhigh_int,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(1),ierr)
    endif
    if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
        call MPI_IRecv(Qzlow_ext,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(2),ierr)
    endif
    if (nbrs(UP) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(1),stats(:,1),ierr)
        call MPI_IRecv(Qzhigh_ext,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(3),ierr)
    endif
    if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(2),stats(:,2),ierr)
        call MPI_ISend(Qzlow_int,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(4),ierr)
    endif
    if (nbrs(UP) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(3),stats(:,3),ierr)
    endif

    if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
        call MPI_Wait(reqs(4),stats(:,4),ierr)
    endif


    if (mpi_R .eq. 1 .and. zlbc .eq. 0) then
        Qzlow_ext = Qzlow_int
    end if

    if (mpi_R .eq. mpi_nz .and. zhbc .eq. 0) then
        Qzhigh_ext = Qzhigh_int
    end if

end subroutine


!-----------------------------------------------------------

    real function xc(i)
        integer i
        xc = loc_lxd + (i - 0.5)*dx
    end function xc

!-----------------------------------------------------------

    real function yc(j)
        integer j
        yc = loc_lyd + (j - 0.5)*dy
    end function yc

!-----------------------------------------------------------

    real function zc(k)
        integer k
        zc = loc_lzd + (k - 0.5)*dz
    end function zc

!-----------------------------------------------------------

    real function rz(i,j)
        integer i,j
        rz = sqrt(yc(j)**2)
    end function rz

!-----------------------------------------------------------

    real function r(i,j)
        integer i,j,k
        r = sqrt(xc(i)**2 + yc(j)**2)
    end function r

!-----------------------------------------------------------

    real function theta(i,j)
        integer i,j
        theta = atan2(yc(j),xc(i))
    end function theta

!-----------------------------------------------------------

    real function xvtk(i)
        integer i
        xvtk = loc_lxd + (i - 0.5)*dxvtk
    end function xvtk

!-----------------------------------------------------------

    real function yvtk(j)
        integer j
        yvtk = loc_lyd + (j - 0.5)*dyvtk
    end function yvtk

!-----------------------------------------------------------

    real function zvtk(k)
        integer k
        zvtk = loc_lzd + (k - 0.5)*dzvtk
    end function zvtk

!-----------------------------------------------------------

    real function rvtk(i,j)
        integer i,j,k
        rvtk = sqrt(xvtk(i)**2 + yvtk(j)**2)
    end function rvtk

!-----------------------------------------------------------

    real function thetavtk(i,j)
        integer i,j
        thetavtk = atan2(yvtk(j),xvtk(i))
    end function thetavtk


!------------More realistic COBRA  current driver -----------

    real function E_z(t)
        real t
        if(t.le.tr) E_z = Ez0*sin(0.5*pi*t/tr)*cos(0.5*pi*t/tr)
        if(t.ge.tr) E_z = 0.
    end function E_z

!------------More realistic COBRA  current driver -----------


    real function Icur(t)
        real t
        if(t.le.tr) Icur = Ipeak*sin(0.5*pi*t/tr)**2
        if(t.ge.tr) Icur = Ipeak
!        if(t.ge.tr.and.t.le.2*tr)Icur = Ipeak
!        if(t.ge.2*tr.and.t.le.3*tr)Icur = Ipeak*sin(0.5*pi*(t-tr)/tr)
!        if(t.ge.3*tr)Icur = 0.
    end function Icur

!-----------------------------------------------------------

    subroutine output_vtk0(Qin,nout,iam)

    implicit none
    real Qin(nx,ny,nz,nQ,nbasis)
    integer nout
    integer(I4P) , parameter :: nb=1*ngu,sprd=0*1
!    integer(I4P), parameter:: nnx=nx-2*nb,nny=ny-2*nb,nnz=nz-2*nb
    real(R4P), dimension(nnx+1):: x_xml_rect
    real(R4P), dimension(nny+1):: y_xml_rect
    real(R4P), dimension(nnz+1):: z_xml_rect
    real(R4P), dimension(nnx*nny*nnz):: var_xml_val_x
    real(R4P), dimension(nnx*nny*nnz):: var_xml_val_y
    real(R4P), dimension(nnx*nny*nnz):: var_xml_val_z
    real P, vx, vy, vz,  dni
    integer(I4P):: E_IO,i,j,k,l,num,iam
    character (50) :: out_name
    character (4) :: tname
    character (5) :: tname1
    character (4) :: pname
    character (5) :: pname1

    num=nout+10000

    write(tname1,'(i5)')num
    tname=tname1(2:5)
    tname = trim(tname)
    tname = adjustr(tname)

    num=iam+10000

    write(pname1,'(i5)')num
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)
    ! out_name='/data/data5/perseus_p'//pname//'_t'//tname//'.vtr'
    out_name='data/perseus_p'//pname//'_t'//tname//'.vtr'
    ! print *, out_name
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
                if(ieos .eq. 1)P = (aindex - 1)*(Qin(i,j,k,en,1) - 0.5*Qin(i,j,k,rh,1)*(vx**2 + vy**2 + vz**2))

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

!-----------------------------------------------------------

subroutine output_vtk(Qin,nout,iam)

    implicit none
    real Qin(nx,ny,nz,nQ,nbasis)
    integer nout
    integer(I4P) , parameter :: nb=1*ngu,sprd=0*1

    real(R4P), dimension(nnx+1) :: x_xml_rect
    real(R4P), dimension(nny+1) :: y_xml_rect
    real(R4P), dimension(nnz+1) :: z_xml_rect
    real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_x
    real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_y
    real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_z
    real(R4P), dimension(nnx,nny,nnz,nQ) :: qvtk
    real(R4P), dimension(nnx,nny,nnz) :: qvtk_dxvy,qvtk_dyvx
    real P, vx, vy, vz,  dni, dxrh,dyrh,dxmy,dymx
    integer(I4P):: E_IO,i,j,k,l,num,iam,igrid,ir,jr,kr,ib,jb,kb,ieq
    character (50) :: out_name
    character (4) :: tname
    character (5) :: tname1
    character (4) :: pname
    character (5) :: pname1

    num=nout+10000

    write(tname1,'(i5)')num
    tname=tname1(2:5)
    tname = trim(tname)
    tname = adjustr(tname)

    num=iam+10000

    write(pname1,'(i5)')num
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)
    ! out_name='/data/data/perseus_p'//pname//'_t'//tname//'.vtr'
    out_name='data/perseus_p'//pname//'_t'//tname//'.vtr'
    ! print *, out_name
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
        if (ib .eq. 0) then
            ib = nvtk
            ir = ir - 1
        end if

        do j=1+nb,nny-nb

            jr = int(j/nvtk) + 1
            jb = j - nvtk*int(j/nvtk)
            if (jb .eq. 0) then
                jb = nvtk
                jr = jr - 1
            end if

            do k=1+nb,nnz-nb

                kr = int(k/nvtk) + 1
                kb = k - nvtk*int(k/nvtk)
                if (kb .eq. 0) then
                    kb = nvtk
                    kr = kr - 1
                end if
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
                l = (i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
                dni = 1./qvtk(i,j,k,rh)
                vx = qvtk(i,j,k,mx)*dni
                vy = qvtk(i,j,k,my)*dni
                vz = qvtk(i,j,k,mz)*dni
                if(ieos .eq. 1) P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))
                if(ieos .eq. 2)then
                    P = P_1*(qvtk(i,j,k,rh)**7.2 - 1.) + P_base
                    if(P < P_floor) P= P_floor
                end if
                var_xml_val_x(l) = P*P0
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

!--------------------------------------------------------------------------------

subroutine writeQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow)

    implicit none
    real :: Qin(nx,ny,nz,nQ,nbasis),tnow,dtnow
    integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow,nump,numd,qq,k,j,i,ir
    character (4) :: fprefix,pname,dname
    character (5) :: pname1,dname1
    character (30) :: fname2

    nump = iam + 10000

    write(pname1,'(i5)')nump
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)

    numd = iddump + 10000

    write(dname1,'(i5)')numd
    dname=dname1(2:5)
    dname = trim(dname)
    dname = adjustr(dname)

    fname2 = 'data/'//fprefix//'_p'//pname//'_d'//dname//'.dat'
    ! print *,'fname2 ',fname2

    open(unit=3,file=fname2)

    ! open(unit = 10, file = 'data/perseus_t'//dname//'_p'//pname//'.bin',form = 'unformatted',access = 'stream')

    do ir=1,nbasis
        do qq=1,nQ
            do k=1,nz
                do j=1,ny
                    write(3,*) (Qin(i,j,k,qq,ir),i=1,nx)
                enddo
            enddo
        enddo
    enddo

    write(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow
    close(3)

end subroutine writeQ

!------------------------------------------------------------------------------

subroutine readQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow)

    implicit none
    real :: Qin(nx,ny,nz,nQ,nbasis),tnow,dtnow
    integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow,nump,numd,qq,k,j,i,ir
    character (4) :: fprefix,pname,dname
    character (5) :: pname1,dname1
    character (30) :: fname2

    nump = irank + 10000

    write(pname1,'(i5)')nump
    pname=pname1(2:5)
    pname = trim(pname)
    pname = adjustr(pname)

    numd = iddump + 10000

    write(dname1,'(i5)')numd
    dname=dname1(2:5)
    dname = trim(dname)
    dname = adjustr(dname)

    fname2 = 'data/'//fpre//'_p'//pname//'_d'//dname//'.dat'

    open(unit=3,file=fname2,action='read')

    do ir=1,nbasis
    do qq=1,nQ
        do k=1,nz
            do j=1,ny
                read(3,*) (Qin(i,j,k,qq,ir),i=1,nx)
            enddo
        enddo
    enddo
    enddo

    read(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow
    close(3)

end subroutine readQ

!--------------------------------------------------------------------------------

subroutine set_bfvals_3D
    ! Defines local basis function values and weights for 1, 2, or 3-point Gaussian quadrature.
    ! Basis functions are evaluated in cell interior and on cell faces.

    implicit none

    call set_vtk_vals_3D()    ! Define basis function values on a 3D grid INTERNAL to the cell.
    call set_internal_vals_3D()     ! Define basis function values at quadrature points INTERNAL to cell.
    call set_face_vals_3D()   ! Define local basis function values at quadrature points on a cell face.
    call set_weights_3D()     ! Define weights for integral approximation using Gaussian quadrature.

end subroutine set_bfvals_3D

!----------------------------------------------------------
subroutine set_vtk_vals_3D
    ! Define basis function values on a rectangular 3D grid INTERNAL to the cell.
    ! For use in output_vtk().

    implicit none
    integer ixyz,i,igrid

    dxvtk = 1./(nvtk*dxi)
    dyvtk = 1./(nvtk*dyi)
    dzvtk = 1./(nvtk*dzi)

    if (nvtk .eq. 1) xgrid(1) = 0.

    if (nvtk .eq. 2) then
        xgrid(1) = -0.5
        xgrid(2) = 0.5
    end if
    if (nvtk .eq. 3) then
        xgrid(1) = -2./3.
        xgrid(2) = 0.
        xgrid(3) = 2./3.
    end if
    if (nvtk .eq. 4) then
        xgrid(1) = -3./4.
        xgrid(2) = -1./4.
        xgrid(3) = 1./4.
        xgrid(4) = 3./4.
    end if


    bfvtk(1,1) = 1.        ! basis function = 1
    bfvtk(1,kx) = 0.        ! basis function = x
    bfvtk(1,ky) = 0.        ! basis function = y
    bfvtk(1,kz) = 0.        ! basis function = z

    bfvtk(1:nvtk3,1) = 1.        ! basis function = 1
    do i=1,nvtk
        bfvtk((i-1)*nvtk+1:i*nvtk,ky) = xgrid(i)        ! basis function = y
        bfvtk((i-1)*nvtk+1:i*nvtk,kz) = xgrid(1:nvtk)        ! basis function = z
    end do
    do i=1,nvtk
        bfvtk((i-1)*nvtk2+1:i*nvtk2,kx) = xgrid(i)        ! basis function = x
        bfvtk((i-1)*nvtk2+1:i*nvtk2,ky) = bfvtk(1:nvtk2,ky)        ! basis function = y
        bfvtk((i-1)*nvtk2+1:i*nvtk2,kz) = bfvtk(1:nvtk2,kz)        ! basis function = z
    end do

    do i=0,2
        bfvtk(1:nvtk3,kxx+i) = 1.5*bfvtk(1:nvtk3,kx+i)**2 - 0.5    ! basis function = P_2(s)
        ! bfvtk(1:nvtk3,kxxx+i) = 2.5*bfvtk(1:nvtk3,kx+i)**3 - 1.5*bfvtk(1:nvtk3,kx+i)   ! basis function = P3(s)
    end do

    bfvtk(1:nvtk3,kyz) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kz)        ! basis function = yz
    bfvtk(1:nvtk3,kzx) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kx)        ! basis function = zx
    bfvtk(1:nvtk3,kxy) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,ky)        ! basis function = xy
    bfvtk(1:nvtk3,kxyz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kz)     ! basis function = xyz

    bfvtk(1:nvtk3,kyyz) = bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kz)     ! basis function = P2(y)z
    bfvtk(1:nvtk3,kyzz) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kzz)     ! basis function = P2(z)y
    bfvtk(1:nvtk3,kzzx) = bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kx)     ! basis function = P2(z)x
    bfvtk(1:nvtk3,kzxx) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kxx)     ! basis function = P2(x)z
    bfvtk(1:nvtk3,kxxy) = bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,ky)     ! basis function = P2(x)y
    bfvtk(1:nvtk3,kxyy) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)     ! basis function = P2(y)x
    bfvtk(1:nvtk3,kyyzz) = bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)     ! basis function = P_2(y)P_2(z)
    bfvtk(1:nvtk3,kzzxx) = bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kxx)     ! basis function = P_2(z)P_2(x)
    bfvtk(1:nvtk3,kxxyy) = bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,kyy)     ! basis function = P_2(x)P_2(y)
    bfvtk(1:nvtk3,kyzxx) = bfvtk(1:nvtk3,kyz)*bfvtk(1:nvtk3,kxx)     ! basis function = yz P_2(x)
    bfvtk(1:nvtk3,kzxyy) = bfvtk(1:nvtk3,kzx)*bfvtk(1:nvtk3,kyy)     ! basis function = zx P_2(y)
    bfvtk(1:nvtk3,kxyzz) = bfvtk(1:nvtk3,kxy)*bfvtk(1:nvtk3,kzz)     ! basis function = xy P_2(z)
    bfvtk(1:nvtk3,kxyyzz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)   ! basis function = x P_2(y)P_2(z)
    bfvtk(1:nvtk3,kyzzxx) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kxx)   ! basis function = y P_2(z)P_2(x)
    bfvtk(1:nvtk3,kzxxyy) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,kyy)   ! basis function = z P_2(x)P_2(y)
    bfvtk(1:nvtk3,kxxyyzz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

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

        ! bfvtk_dx(igrid,kxxx) = 7.5*bfvtk(igrid,kx)**2 - 1.5
        ! bfvtk_dx(igrid,kyyy:kzzz) = 0.
        ! bfvtk_dy(igrid,kyyy) = 7.5*bfvtk(igrid,ky)**2 - 1.5
        ! bfvtk_dy(igrid,kxxx) = 0.
        ! bfvtk_dy(igrid,kzzz) = 0.
        ! bfvtk_dz(igrid,kzzz) = 7.5*bfvtk(igrid,kz)**2 - 1.5
        ! bfvtk_dz(igrid,kxxx) = 0.
        ! bfvtk_dz(igrid,kyyy) = 0.
    end do

end subroutine set_vtk_vals_3D

!----------------------------------------------------------

!----------------------------------------------------------
subroutine set_internal_vals_3D
    ! Define basis function values at quadrature points INTERNAL to cell.

    implicit none
    integer ixyz, i
    real c15d5,c1dsq3,xq4p,xq4m

    c15d5 = sqrt(15.)/5.
    c1dsq3 = 1./sqrt(3.)

    xq4p = sqrt(525. + 70.*sqrt(30.))/35.
    xq4m = sqrt(525. - 70.*sqrt(30.))/35.


    if (iquad .eq. 2) then
        xquad(1) = -c1dsq3
        xquad(2) = c1dsq3
    end if
    if (iquad .eq. 3) then
        xquad(1) = -c15d5
        xquad(2) = 0.
        xquad(3) = c15d5
    end if
    if (iquad .eq. 4) then
        xquad(1) = -xq4p
        xquad(2) = -xq4m
        xquad(3) = xq4m
        xquad(4) = xq4p
    end if

    if (iquad .eq. 1) then            ! 2-point Gaussian quadrature
        bfvals_int(1,1) = 1.        ! basis function = 1
        bfvals_int(1,kx) = 0.        ! basis function = x
        bfvals_int(1,ky) = 0.        ! basis function = y
        bfvals_int(1,kz) = 0.        ! basis function = z
    end if


    if (iquad .gt. 1) then
        bfvals_int(1:npg,1) = 1.        ! basis function = 1
        do i=1,nedge
            bfvals_int((i-1)*nedge+1:i*nedge,ky) = xquad(i)        ! basis function = y
            bfvals_int((i-1)*nedge+1:i*nedge,kz) = xquad(1:nedge)        ! basis function = z
        end do
        do i=1,nedge
            bfvals_int((i-1)*nface+1:i*nface,kx) = xquad(i)        ! basis function = x
            bfvals_int((i-1)*nface+1:i*nface,ky) = bfvals_int(1:nface,ky)        ! basis function = y
            bfvals_int((i-1)*nface+1:i*nface,kz) = bfvals_int(1:nface,kz)        ! basis function = z
        end do
    end if

    do i=0,2
        bfvals_int(1:npg,kxx+i) = 1.5*bfvals_int(1:npg,kx+i)**2 - 0.5    ! basis function = P_2(s)
        ! bfvals_int(1:npg,kxxx+i) = 2.5*bfvals_int(1:npg,kx+i)**3 - 1.5*bfvals_int(1:npg,kx+i)   ! basis function = P_3(s)
    end do

    bfvals_int(1:npg,kyz) = bfvals_int(1:npg,ky)*bfvals_int(1:npg,kz)        ! basis function = yz
    bfvals_int(1:npg,kzx) = bfvals_int(1:npg,kz)*bfvals_int(1:npg,kx)        ! basis function = zx
    bfvals_int(1:npg,kxy) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,ky)        ! basis function = xy
    bfvals_int(1:npg,kxyz) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,ky)*bfvals_int(1:npg,kz)     ! basis function = xyz

    bfvals_int(1:npg,kyyz) = bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kz)     ! basis function = P_2(y)z
    bfvals_int(1:npg,kyzz) = bfvals_int(1:npg,ky)*bfvals_int(1:npg,kzz)     ! basis function = P_2(z)y
    bfvals_int(1:npg,kzzx) = bfvals_int(1:npg,kzz)*bfvals_int(1:npg,kx)     ! basis function = P_2(z)x
    bfvals_int(1:npg,kzxx) = bfvals_int(1:npg,kz)*bfvals_int(1:npg,kxx)     ! basis function = P_2(x)z
    bfvals_int(1:npg,kxxy) = bfvals_int(1:npg,kxx)*bfvals_int(1:npg,ky)     ! basis function = P_2(x)y
    bfvals_int(1:npg,kxyy) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,kyy)     ! basis function = P_2(y)x
    bfvals_int(1:npg,kyyzz) = bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kzz)     ! basis function = P_2(y)P_2(z)
    bfvals_int(1:npg,kzzxx) = bfvals_int(1:npg,kzz)*bfvals_int(1:npg,kxx)     ! basis function = P_2(z)P_2(x)
    bfvals_int(1:npg,kxxyy) = bfvals_int(1:npg,kxx)*bfvals_int(1:npg,kyy)     ! basis function = P_2(x)P_2(y)
    bfvals_int(1:npg,kyzxx) = bfvals_int(1:npg,kyz)*bfvals_int(1:npg,kxx)     ! basis function = yz P_2(x)
    bfvals_int(1:npg,kzxyy) = bfvals_int(1:npg,kzx)*bfvals_int(1:npg,kyy)     ! basis function = zx P_2(y)
    bfvals_int(1:npg,kxyzz) = bfvals_int(1:npg,kxy)*bfvals_int(1:npg,kzz)     ! basis function = xy P_2(z)
    bfvals_int(1:npg,kxyyzz) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kzz)   ! basis function = x P_2(y)P_2(z)
    bfvals_int(1:npg,kyzzxx) = bfvals_int(1:npg,ky)*bfvals_int(1:npg,kzz)*bfvals_int(1:npg,kxx)   ! basis function = y P_2(z)P_2(x)
    bfvals_int(1:npg,kzxxyy) = bfvals_int(1:npg,kz)*bfvals_int(1:npg,kxx)*bfvals_int(1:npg,kyy)   ! basis function = z P_2(x)P_2(y)
    bfvals_int(1:npg,kxxyyzz) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

end subroutine set_internal_vals_3D

!----------------------------------------------------------
subroutine set_face_vals_3D
    ! Define local basis function values at quadrature points on a cell face.
    ! Used in flux_cal() and prepare_exchange() for interpolating from cell center onto cell face.
    ! Also used in glflux().

    implicit none
    integer ixyz,i
    real c15d5,c1dsq3,xq4p,xq4m

    c15d5 = sqrt(15.)/5.
    c1dsq3 = 1./sqrt(3.)

    xq4p = sqrt(525. + 70.*sqrt(30.))/35.
    xq4m = sqrt(525. - 70.*sqrt(30.))/35.

    if (iquad .eq. 2) then
        xquad(1) = -c1dsq3
        xquad(2) = c1dsq3
    end if
    if (iquad .eq. 3) then
        xquad(1) = -c15d5
        xquad(2) = 0.
        xquad(3) = c15d5
    end if
    if (iquad .eq. 4) then
        xquad(1) = -xq4p
        xquad(2) = -xq4m
        xquad(3) = xq4m
        xquad(4) = xq4p
    end if

    ! bfvals_rsp:  positive rs-face
    ! bfvals_rsm:  negative rs-face
    ! bfvals_rsp(,1):  value=1 on positive rs-face
    ! bfvals_rsp(,kx):  x-values on positive rs-face
    ! bfvals_rsp(,ky):  y-values on positive rs-face
    ! bfvals_rsp(,kz):  z-values on positive rs-face
    ! bfvals_rsp(,kxx):  P_2(x)-values on positive rs-face
    ! bfvals_rsp(,kyy):  P_2(y)-values on positive rs-face
    ! bfvals_rsp(,kzz):  P_2(z)-values on positive rs-face
    ! bfvals_rsp(,kyz):  yz-values on positive rs-face
    ! bfvals_rsp(,kzx):  zx-values on positive rs-face
    ! bfvals_rsp(,kxy):  xy-values on positive rs-face

    bfvals_zp(1:nface,1) = 1.
    if (iquad .eq. 1) then
        bfvals_zp(1,kx) = 0.
        bfvals_zp(1,ky) = 0.
    end if

    if (iquad .gt. 1) then
        do i=1,nedge
            bfvals_zp((i-1)*nedge+1:i*nedge,kx) = xquad(1:nedge)
            bfvals_zp((i-1)*nedge+1:i*nedge,ky) = xquad(i)
        end do
    end if

    bfvals_zp(1:nface,kz) = 1.
    bfvals_zp(1:nface,kyz) = bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kz)
    bfvals_zp(1:nface,kzx) = bfvals_zp(1:nface,kz)*bfvals_zp(1:nface,kx)
    bfvals_zp(1:nface,kxy) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,ky)
    bfvals_zp(1:nface,kxyz) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kz)     ! basis function = xyz

    bfvals_yp(1:nface,1) = 1.
    bfvals_yp(1:nface,kx) = bfvals_zp(1:nface,ky)
    bfvals_yp(1:nface,ky) = 1.
    bfvals_yp(1:nface,kz) = bfvals_zp(1:nface,kx)
    bfvals_yp(1:nface,kyz) = bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kz)
    bfvals_yp(1:nface,kzx) = bfvals_yp(1:nface,kz)*bfvals_yp(1:nface,kx)
    bfvals_yp(1:nface,kxy) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,ky)
    bfvals_yp(1:nface,kxyz) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kz)     ! basis function = xyz

    bfvals_xp(1:nface,1) = 1.
    bfvals_xp(1:nface,kx) = 1.
    bfvals_xp(1:nface,ky) = bfvals_zp(1:nface,kx)
    bfvals_xp(1:nface,kz) = bfvals_zp(1:nface,ky)
    bfvals_xp(1:nface,kyz) = bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kz)
    bfvals_xp(1:nface,kzx) = bfvals_xp(1:nface,kz)*bfvals_xp(1:nface,kx)
    bfvals_xp(1:nface,kxy) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,ky)
    bfvals_xp(1:nface,kxyz) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kz)     ! basis function = xyz

    do i=1,3
        bfvals_xp(1:nface,kxx+i-1) = 1.5*bfvals_xp(1:nface,kx+i-1)**2 - 0.5
        bfvals_yp(1:nface,kxx+i-1) = 1.5*bfvals_yp(1:nface,kx+i-1)**2 - 0.5
        bfvals_zp(1:nface,kxx+i-1) = 1.5*bfvals_zp(1:nface,kx+i-1)**2 - 0.5
        ! bfvals_xp(1:nface,kxxx+i-1) = 2.5*bfvals_xp(1:nface,kx+i-1)**3 - 1.5*bfvals_xp(1:nface,kx+i-1)
        ! bfvals_yp(1:nface,kxxx+i-1) = 2.5*bfvals_yp(1:nface,kx+i-1)**3 - 1.5*bfvals_yp(1:nface,kx+i-1)
        ! bfvals_zp(1:nface,kxxx+i-1) = 2.5*bfvals_zp(1:nface,kx+i-1)**3 - 1.5*bfvals_zp(1:nface,kx+i-1)
    end do

    bfvals_xp(1:nface,kyyz) = bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kz)     ! basis function = P_2(y)z
    bfvals_xp(1:nface,kyzz) = bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kzz)     ! basis function = P_2(z)y
    bfvals_xp(1:nface,kzzx) = bfvals_xp(1:nface,kzz)*bfvals_xp(1:nface,kx)     ! basis function = P_2(z)x
    bfvals_xp(1:nface,kzxx) = bfvals_xp(1:nface,kz)*bfvals_xp(1:nface,kxx)     ! basis function = P_2(x)z
    bfvals_xp(1:nface,kxxy) = bfvals_xp(1:nface,kxx)*bfvals_xp(1:nface,ky)     ! basis function = P_2(x)y
    bfvals_xp(1:nface,kxyy) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,kyy)     ! basis function = P_2(y)x
    bfvals_xp(1:nface,kyyzz) = bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kzz)     ! basis function = P_2(y)P_2(z)
    bfvals_xp(1:nface,kzzxx) = bfvals_xp(1:nface,kzz)*bfvals_xp(1:nface,kxx)     ! basis function = P_2(z)P_2(x)
    bfvals_xp(1:nface,kxxyy) = bfvals_xp(1:nface,kxx)*bfvals_xp(1:nface,kyy)     ! basis function = P_2(x)P_2(y)
    bfvals_xp(1:nface,kyzxx) = bfvals_xp(1:nface,kyz)*bfvals_xp(1:nface,kxx)     ! basis function = yz P_2(x)
    bfvals_xp(1:nface,kzxyy) = bfvals_xp(1:nface,kzx)*bfvals_xp(1:nface,kyy)     ! basis function = zx P_2(y)
    bfvals_xp(1:nface,kxyzz) = bfvals_xp(1:nface,kxy)*bfvals_xp(1:nface,kzz)     ! basis function = xy P_2(z)
    bfvals_xp(1:nface,kxyyzz) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kzz)   ! basis function = x P_2(y)P_2(z)
    bfvals_xp(1:nface,kyzzxx) = bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kzz)*bfvals_xp(1:nface,kxx)   ! basis function = y P_2(z)P_2(x)
    bfvals_xp(1:nface,kzxxyy) = bfvals_xp(1:nface,kz)*bfvals_xp(1:nface,kxx)*bfvals_xp(1:nface,kyy)   ! basis function = z P_2(x)P_2(y)
    bfvals_xp(1:nface,kxxyyzz) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

    bfvals_xm = bfvals_xp
    bfvals_xm(1:nface,kx) = -bfvals_xp(1:nface,kx)
    bfvals_xm(1:nface,kzx) = -bfvals_xp(1:nface,kzx)
    bfvals_xm(1:nface,kxy) = -bfvals_xp(1:nface,kxy)
    bfvals_xm(1:nface,kxyz) = -bfvals_xp(1:nface,kxyz)
    bfvals_xm(1:nface,kzzx) = -bfvals_xp(1:nface,kzzx)
    bfvals_xm(1:nface,kxyy) = -bfvals_xp(1:nface,kxyy)
    ! bfvals_xm(1:nface,kxxx) = -bfvals_xp(1:nface,kxxx)
    bfvals_xm(1:nface,kzxyy) = -bfvals_xp(1:nface,kzxyy)
    bfvals_xm(1:nface,kxyzz) = -bfvals_xp(1:nface,kxyzz)
    bfvals_xm(1:nface,kxyyzz) = -bfvals_xp(1:nface,kxyyzz)

    bfvals_yp(1:nface,kyyz) = bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kz)     ! basis function = P_2(y)z
    bfvals_yp(1:nface,kyzz) = bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kzz)     ! basis function = P_2(z)y
    bfvals_yp(1:nface,kzzx) = bfvals_yp(1:nface,kzz)*bfvals_yp(1:nface,kx)     ! basis function = P_2(z)x
    bfvals_yp(1:nface,kzxx) = bfvals_yp(1:nface,kz)*bfvals_yp(1:nface,kxx)     ! basis function = P_2(x)z
    bfvals_yp(1:nface,kxxy) = bfvals_yp(1:nface,kxx)*bfvals_yp(1:nface,ky)     ! basis function = P_2(x)y
    bfvals_yp(1:nface,kxyy) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,kyy)     ! basis function = P_2(y)x
    bfvals_yp(1:nface,kyyzz) = bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kzz)     ! basis function = P_2(y)P_2(z)
    bfvals_yp(1:nface,kzzxx) = bfvals_yp(1:nface,kzz)*bfvals_yp(1:nface,kxx)     ! basis function = P_2(z)P_2(x)
    bfvals_yp(1:nface,kxxyy) = bfvals_yp(1:nface,kxx)*bfvals_yp(1:nface,kyy)     ! basis function = P_2(x)P_2(y)
    bfvals_yp(1:nface,kyzxx) = bfvals_yp(1:nface,kyz)*bfvals_yp(1:nface,kxx)     ! basis function = yz P_2(x)
    bfvals_yp(1:nface,kzxyy) = bfvals_yp(1:nface,kzx)*bfvals_yp(1:nface,kyy)     ! basis function = zx P_2(y)
    bfvals_yp(1:nface,kxyzz) = bfvals_yp(1:nface,kxy)*bfvals_yp(1:nface,kzz)     ! basis function = xy P_2(z)
    bfvals_yp(1:nface,kxyyzz) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kzz)   ! basis function = x P_2(y)P_2(z)
    bfvals_yp(1:nface,kyzzxx) = bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kzz)*bfvals_yp(1:nface,kxx)   ! basis function = y P_2(z)P_2(x)
    bfvals_yp(1:nface,kzxxyy) = bfvals_yp(1:nface,kz)*bfvals_yp(1:nface,kxx)*bfvals_yp(1:nface,kyy)   ! basis function = z P_2(x)P_2(y)
    bfvals_yp(1:nface,kxxyyzz) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

    bfvals_ym = bfvals_yp
    bfvals_ym(1:nface,ky) = -bfvals_yp(1:nface,ky)
    bfvals_ym(1:nface,kyz) = -bfvals_yp(1:nface,kyz)
    bfvals_ym(1:nface,kxy) = -bfvals_yp(1:nface,kxy)
    bfvals_ym(1:nface,kxyz) = -bfvals_yp(1:nface,kxyz)
    bfvals_ym(1:nface,kyzz) = -bfvals_yp(1:nface,kyzz)
    bfvals_ym(1:nface,kxxy) = -bfvals_yp(1:nface,kxxy)
    ! bfvals_ym(1:nface,kyyy) = -bfvals_yp(1:nface,kyyy)
    bfvals_ym(1:nface,kyzxx) = -bfvals_yp(1:nface,kyzxx)
    bfvals_ym(1:nface,kxyzz) = -bfvals_yp(1:nface,kxyzz)
    bfvals_ym(1:nface,kyzzxx) = -bfvals_yp(1:nface,kyzzxx)

    bfvals_zp(1:nface,kyyz) = bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kz)     ! basis function = P_2(y)z
    bfvals_zp(1:nface,kyzz) = bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kzz)     ! basis function = P_2(z)y
    bfvals_zp(1:nface,kzzx) = bfvals_zp(1:nface,kzz)*bfvals_zp(1:nface,kx)     ! basis function = P_2(z)x
    bfvals_zp(1:nface,kzxx) = bfvals_zp(1:nface,kz)*bfvals_zp(1:nface,kxx)     ! basis function = P_2(x)z
    bfvals_zp(1:nface,kxxy) = bfvals_zp(1:nface,kxx)*bfvals_zp(1:nface,ky)     ! basis function = P_2(x)y
    bfvals_zp(1:nface,kxyy) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,kyy)     ! basis function = P_2(y)x
    bfvals_zp(1:nface,kyyzz) = bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kzz)     ! basis function = P_2(y)P_2(z)
    bfvals_zp(1:nface,kzzxx) = bfvals_zp(1:nface,kzz)*bfvals_zp(1:nface,kxx)     ! basis function = P_2(z)P_2(x)
    bfvals_zp(1:nface,kxxyy) = bfvals_zp(1:nface,kxx)*bfvals_zp(1:nface,kyy)     ! basis function = P_2(x)P_2(y)
    bfvals_zp(1:nface,kyzxx) = bfvals_zp(1:nface,kyz)*bfvals_zp(1:nface,kxx)     ! basis function = yz P_2(x)
    bfvals_zp(1:nface,kzxyy) = bfvals_zp(1:nface,kzx)*bfvals_zp(1:nface,kyy)     ! basis function = zx P_2(y)
    bfvals_zp(1:nface,kxyzz) = bfvals_zp(1:nface,kxy)*bfvals_zp(1:nface,kzz)     ! basis function = xy P_2(z)
    bfvals_zp(1:nface,kxyyzz) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kzz)   ! basis function = x P_2(y)P_2(z)
    bfvals_zp(1:nface,kyzzxx) = bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kzz)*bfvals_zp(1:nface,kxx)   ! basis function = y P_2(z)P_2(x)
    bfvals_zp(1:nface,kzxxyy) = bfvals_zp(1:nface,kz)*bfvals_zp(1:nface,kxx)*bfvals_zp(1:nface,kyy)   ! basis function = z P_2(x)P_2(y)
    bfvals_zp(1:nface,kxxyyzz) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)


    bfvals_zm = bfvals_zp
    bfvals_zm(1:nface,kz) = -bfvals_zp(1:nface,kz)
    bfvals_zm(1:nface,kyz) = -bfvals_zp(1:nface,kyz)
    bfvals_zm(1:nface,kzx) = -bfvals_zp(1:nface,kzx)
    bfvals_zm(1:nface,kxyz) = -bfvals_zp(1:nface,kxyz)
    bfvals_zm(1:nface,kyyz) = -bfvals_zp(1:nface,kyyz)
    bfvals_zm(1:nface,kzxx) = -bfvals_zp(1:nface,kzxx)
    ! bfvals_zm(1:nface,kzzz) = -bfvals_zp(1:nface,kzzz)
    bfvals_zm(1:nface,kyzxx) = -bfvals_zp(1:nface,kyzxx)
    bfvals_zm(1:nface,kzxyy) = -bfvals_zp(1:nface,kzxyy)
    bfvals_zm(1:nface,kzxxyy) = -bfvals_zp(1:nface,kzxxyy)

    ! Organize local basis values on faces into 1-D vectors.
    ! Used in limiter() and max_lim().

    bf_faces(1:nslim,1) = 1.        ! basis function = 1

    do ixyz=kx,kz
        bf_faces(1:nface,ixyz) = bfvals_xm(1:nface,ixyz)        ! basis function = x,y,z
        bf_faces(nface+1:2*nface,ixyz) = bfvals_xp(1:nface,ixyz)        ! basis function = x,y,z
        bf_faces(2*nface+1:3*nface,ixyz) = bfvals_ym(1:nface,ixyz)        ! basis function = x,y,z
        bf_faces(3*nface+1:4*nface,ixyz) = bfvals_yp(1:nface,ixyz)        ! basis function = x,y,z
        bf_faces(4*nface+1:5*nface,ixyz) = bfvals_zm(1:nface,ixyz)        ! basis function = x,y,z
        bf_faces(5*nface+1:6*nface,ixyz) = bfvals_zp(1:nface,ixyz)        ! basis function = x,y,z
        bf_faces(6*nface+1:nslim,ixyz) = bfvals_int(1:npg,ixyz)        ! basis function = x,y,z
    end do

    bf_faces(1:nslim,kyz) = bf_faces(1:nslim,ky)*bf_faces(1:nslim,kz)     ! basis function = yz
    bf_faces(1:nslim,kzx) = bf_faces(1:nslim,kz)*bf_faces(1:nslim,kx)     ! basis function = zx
    bf_faces(1:nslim,kxy) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,ky)     ! basis function = xy
    bf_faces(1:nslim,kxyz) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,ky)*bf_faces(1:nslim,kz)     ! basis function = xyz

    do i=0,2
        bf_faces(1:nslim,kxx+i) = 1.5*bf_faces(1:nslim,kx+i)**2 - 0.5    ! basis function = P_2(s)
        ! bf_faces(1:nslim,kxxx+i) = 2.5*bf_faces(1:nslim,kx+i)**3 - 1.5*bf_faces(1:nslim,kx+i)   ! basis function = P_3(s)
    end do


    bf_faces(1:nslim,kyyz) = bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kz)     ! basis function = P_2(y)z
    bf_faces(1:nslim,kyzz) = bf_faces(1:nslim,ky)*bf_faces(1:nslim,kzz)     ! basis function = P_2(z)y
    bf_faces(1:nslim,kzzx) = bf_faces(1:nslim,kzz)*bf_faces(1:nslim,kx)     ! basis function = P_2(z)x
    bf_faces(1:nslim,kzxx) = bf_faces(1:nslim,kz)*bf_faces(1:nslim,kxx)     ! basis function = P_2(x)z
    bf_faces(1:nslim,kxxy) = bf_faces(1:nslim,kxx)*bf_faces(1:nslim,ky)     ! basis function = P_2(x)y
    bf_faces(1:nslim,kxyy) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,kyy)     ! basis function = P_2(y)x
    bf_faces(1:nslim,kyyzz) = bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kzz)     ! basis function = P_2(y)P_2(z)
    bf_faces(1:nslim,kzzxx) = bf_faces(1:nslim,kzz)*bf_faces(1:nslim,kxx)     ! basis function = P_2(z)P_2(x)
    bf_faces(1:nslim,kxxyy) = bf_faces(1:nslim,kxx)*bf_faces(1:nslim,kyy)     ! basis function = P_2(x)P_2(y)
    bf_faces(1:nslim,kyzxx) = bf_faces(1:nslim,kyz)*bf_faces(1:nslim,kxx)     ! basis function = yz P_2(x)
    bf_faces(1:nslim,kzxyy) = bf_faces(1:nslim,kzx)*bf_faces(1:nslim,kyy)     ! basis function = zx P_2(y)
    bf_faces(1:nslim,kxyzz) = bf_faces(1:nslim,kxy)*bf_faces(1:nslim,kzz)     ! basis function = xy P_2(z)
    bf_faces(1:nslim,kxyyzz) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kzz)   ! basis function = x P_2(y)P_2(z)
    bf_faces(1:nslim,kyzzxx) = bf_faces(1:nslim,ky)*bf_faces(1:nslim,kzz)*bf_faces(1:nslim,kxx)   ! basis function = y P_2(z)P_2(x)
    bf_faces(1:nslim,kzxxyy) = bf_faces(1:nslim,kz)*bf_faces(1:nslim,kxx)*bf_faces(1:nslim,kyy)   ! basis function = z P_2(x)P_2(y)
    bf_faces(1:nslim,kxxyyzz) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

end subroutine set_face_vals_3D

!----------------------------------------------------------
subroutine set_weights_3D

    integer i
    real wq4p,wq4m

    wq4p = (18. + sqrt(30.))/36.
    wq4m = (18. - sqrt(30.))/36.

    ! Define weights for integral approximation using Gaussian quadrature.

    ! Define weights for 1-D integration
    if (iquad .eq. 1) then        ! 1-point quadrature
        wgt1d(1) = 2.
    end if
    if (iquad .eq. 2) then        ! 2-point quadrature
        wgt1d(1:2) = 1.
    end if
    if (iquad .eq. 3) then        ! 3-point quadrature
        wgt1d(1) = 5./9.
        wgt1d(2) = 8./9.
        wgt1d(3) = 5./9.
    end if
    if (iquad .eq. 4) then        ! 4-point quadrature
        wgt1d(1) = wq4m
        wgt1d(2) = wq4p
        wgt1d(3) = wq4p
        wgt1d(4) = wq4m
    end if

    ! Define weights for 2-D integration
    if (iquad .eq. 1) then        ! 1-point quadrature
        wgt2d(1) = 4.
    end if
    if (iquad .eq. 2) then        ! 2-point quadrature
        wgt2d(1:4) = 1.
    end if

    if (iquad .ge. 3) then
        do i= 1,nedge
            wgt2d((i-1)*nedge+1:i*nedge) = wgt1d(1:nedge)*wgt1d(i)
        end do
    end if

    ! Define weights for 3-D integration
    if (iquad .eq. 1) then        ! 1-point quadrature
        wgt3d(1) = 8.
    end if
    if (iquad .eq. 2) then        ! 2-point quadrature
        wgt3d(1:8) = 1.
    end if

    if (iquad .ge. 3) then
        do i= 1,nedge
            wgt3d((i-1)*nface+1:i*nface) = wgt2d(1:nface)*wgt1d(i)
        end do
    end if

end subroutine set_weights_3D

!--------------------------------------------------------------------------------

end program
