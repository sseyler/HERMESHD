module parameters

    use input

    use lib_vtk_io
    use MKL_VSL_TYPE
    use MKL_VSL

    include '/nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include/mpif.h'

    integer, parameter :: rh  = 1                     ! density
    integer, parameter :: mx  = 2, my  = 3,  mz  = 4  ! momenta
    integer, parameter :: en  = 5                     ! scalar energy
    integer, parameter :: pxx = 6, pyy = 7,  pzz = 8  ! isotropic stress
    integer, parameter :: pxy = 9, pxz = 10, pyz = 11 ! deviatoric stress
    integer, parameter :: nQ  = 11                    ! number of field variables

    !===========================================================================
    ! Spatial resolution -- # grid cells and DG basis order
    !------------------------------------------------------------
    ! The jump in accuracy b/w the linear basis (nbasis=4) and quadratic basis
    ! (nbasis=10) is much greater than jump b/w quadratic and cubic (nbasis=20).
    !   nbasis = 4:  {1,x,y,z}
    !   nbasis = 10: nbasis4  + {P_2(x),P_2(y),P_2(z), yz, zx, xy}
    !   nbasis = 20: nbasis10 + {xyz,xP2(y),yP2(x),xP2(z),
    !                                zP2(x),yP2(z),zP2(y),P3(x),P3(y),P3(z)}
    integer, parameter :: nbasis  = 8
    integer, parameter :: nbastot = 27  ! TODO: only used in setup + innerintegral()

    integer, parameter :: ngu = 0  ! TODO: only used in output_vtk()

    ! iquad: # of Gaussian quadrature points per direction. iquad should not be:
    !   < ipoly (max Legendre polynomial order used) --> unstable
    !   > ipoly+1 --> an exact Gaussian quadrature for Legendre poly used
    ! Thus there are only two cases of iquad for a given nbasis. Both give similar
    ! results although iquad = ipoly + 1 is formally more accurate.
    integer, parameter :: nedge = iquad
    integer, parameter :: nface = iquad*iquad    ! number of quadrature points per cell face
    integer, parameter :: npg   = nface*iquad    ! number of internal points per cell
    integer, parameter :: nfe   = 2*nface        ! number of cell face quad points per direction
    integer, parameter :: npge  = 6*nface        ! total number of cell face quad points
    integer, parameter :: nslim = npg + 6*nface  ! total number of quad points per cell
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Constants, and physical and numerical parameters
    !------------------------------------------------------------
    ! Useful constants
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: sqrt2 = 2.0**0.5, sqrt2i = 1.0/sqrt2
    real, parameter :: c1d5 = 1./5.
    real, parameter :: c1d3 = 1./3., c2d3 = 2./3., c4d3 = 4./3.
    real, parameter :: eV_per_K = 8.61728e-5

    ! Useful derived parameters
    real, parameter :: aindm1 = aindex - 1. ! gamma - 1
    real, parameter :: cp = aindex/aindm1   ! specific heat at constant pressure

    ! Dimensional units -- expressed in MKS. NOTE: temperature (te0) in eV!
    real, parameter :: L0 = 1.0e0  ! 1.0e-9                 ! length
    real, parameter :: t0 = 1.0e0  ! 1.0e-12                ! time
    real, parameter :: n0 = 1.0e18 ! 3.32e28                ! number density

    ! Derived units
    real, parameter :: v0  = L0/t0                 ! velocity
    real, parameter :: p0  = mu*1.67e-27*n0*v0**2  ! pressure
    real, parameter :: te0 = p0/n0/1.6e-19         ! temperature (eV, not K!)

    ! rh_min is a min density to be used for ideal gas EOS, rh_min is min density
    ! below which the pressure becomes negative for the MT water EOS.
    ! The DG-based subroutine "limiter" keeps density above rh_mult*rh_min.
    real, parameter :: rh_floor = 1.0e-5     ! 1.0e-1 was old value
    real, parameter :: T_floor  = (te*eV_per_K)/te0  ! 0.02585 eV ~ 300 K
    real, parameter :: P_floor  = T_floor*rh_floor

    ! Murnaghan-Tait EOS
    !   P = P_1*(density**7.2 - 1.) + P_base
    ! Note: the EOS for water is likely to be a critical player in getting the
    ! fluctuating hydrodynamics correct. There are much more sophisicated EOS's,
    ! some of which account for ionic solutions. Would be worthwhile to
    ! further investigate and experiment with different EOS's.
    real, parameter :: n_tm = 7.2  ! 7.2 (or 7.15) for water
    real, parameter :: P_1 = 2.15e9/n_tm/p0, P_base = 1.01e5/p0 ! atmo pressure
    real, parameter :: rh_mult = 1.01, rh_min = rh_mult*(1.0-P_base/P_1)**(1./n_tm)
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Miscellaneous stuff for viscosity, heat transport, GRM generation + more
    !------------------------------------------------------------
    ! NOTE: New way of handling relaxation system:
    !   coll -- plays same role as epsi (from old approach)
    !   colvis = coll*vis = density*temp  (vis is dynamic viscosity, e.g. eta)
    real :: coll, colvis    ! values set in set_ic
    real :: c2d3cv, c4d3cv  ! probably want to set in set_ic (not currently used)

    ! NOTE: Old way of handling relaxation system
    real, parameter :: nu = epsi*vis
    real, parameter :: c2d3nu=c2d3*nu, c4d3nu=c4d3*nu

    real, parameter :: T_base     = 300.0**eV_per_K/te0 ! te*eV_per_K/te0  ! temp (isothermal)
    real, parameter :: eta_base   = vis    ! dynamic viscosity
    real, parameter :: zeta_base  = 0.  ! bulk viscosity---will need to adjust this!
    real, parameter :: kappa_base = 1.e-1

    real, parameter :: eta_sd   = (2.*eta_base*T_base)**0.5  ! stdev of fluctuations for shear viscosity terms
    real, parameter :: zeta_sd  = (zeta_base*T_base/3.)**0.5  ! stdev of fluctuations for bulk viscosity term
    real, parameter :: kappa_sd = (2.*kappa_base*T_base**2)**0.5
    !===========================================================================


    !===========================================================================
    ! MKL VSL parameters
    !------------------------------------------------------------
    real vsl_errcode
    TYPE (VSL_STREAM_STATE) :: vsl_stream

    integer, parameter :: vsl_brng   = VSL_BRNG_MCG31
    integer, parameter :: vsl_method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
    real, parameter :: vsl_mean  = 0.0
    real, parameter :: vsl_sigma = 1.0
    !===========================================================================


    !===========================================================================
    ! Parameters relating to basis functions & VTK output
    !------------------------------------------------------------
    integer, dimension(3) :: mxa ,mya ,mza  ! TODO: currently only used in flux_hllc

    ! Basis function flags
    ! TODO: these variables are in:
    !   * initialize.f90 (setup)
    !   * innerintegral
    !   * set_vtk_vals_3D, set_internal_vals_3D, set_face_vals_3D
    ! integer, parameter :: kx     = 2, ky    = 3, kz    = 4
    ! integer, parameter :: kyz    = 5, kzx   = 6, kxy   = 7
    ! integer, parameter :: kxyz   = 8
    ! integer, parameter :: kxx    = 9, kyy   =10, kzz   =11
    ! integer, parameter :: kyzz   =12, kzxx  =13, kxyy  =14
    ! integer, parameter :: kyyz   =15, kzzx  =16, kxxy  =17
    ! integer, parameter :: kyyzz  =18, kzzxx =19, kxxyy =20
    ! integer, parameter :: kyzxx  =21, kzxyy =22, kxyzz =23
    ! integer, parameter :: kxyyzz =24, kyzzxx=25, kzxxyy=26
    ! integer, parameter :: kxxyyzz=27
    integer kx,ky,kz, kyz,kzx,kxy, kxyz, kxx,kyy,kzz, kyzz,kzxx,kxyy
    integer kyyz,kzzx,kxxy, kyyzz,kzzxx,kxxyy, kyzxx,kzxyy,kxyzz
    integer kxyyzz,kyzzxx,kzxxyy, kxxyyzz, kxxx,kyyy,kzzz

    ! TODO: only in output & basis_funcs
    integer, parameter :: nvtk  = 1    ! (was 2)
    integer, parameter :: nvtk2 = nvtk*nvtk
    integer, parameter :: nvtk3 = nvtk*nvtk*nvtk

    ! TODO: only in set_vtk_vals_3D, output_vtk
    real, dimension(nvtk3,nbastot) :: bfvtk, bfvtk_dx, bfvtk_dy, bfvtk_dz

    ! TODO: only in set_vtk_vals_3D, xvtk/yvtk/zvtk (parameters), output_vtk
    real dxvtk,dyvtk,dzvtk
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Arrays for field variables, fluxes/inner-integrals, and sources, and time(s)
    !------------------------------------------------------------
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r0, Q_r1, Q_r2, Q_r3
    real, dimension(nx,ny,nz,nQ,nbasis) :: glflux_r, source_r, integral_r
    !===========================================================================


    !===========================================================================
    ! Time(s)
    !------------------------------------------------------------
    real t, dt, dtout, sqrt_dVdt_i ! Inv sq-root of (dV*dt), dV = grid cell volume
    !===========================================================================


    !===========================================================================
    ! Helper variables (initialized here)
    !------------------------------------------------------------
    real t1,t2,t_start,t_stop,dtoriginal  ! used for timing (dtoriginal optional)
    real lxd,lxu,lyd,lyu,lzd,lzu  ! used in init + indirectly used by the grid coord functions
    real loc_lxd,loc_lyd,loc_lzd  ! used directly by the grid coord functions
    real dz, dy, dx, dxi, dyi, dzi, dVi  ! used throughout + directly by grid coord functions
    !===========================================================================


    !===========================================================================
    ! MPI definitions
    !------------------------------------------------------------
    integer, parameter :: print_mpi = 0  ! sets the MPI rank that will do printing
    integer :: mpi_nz  ! used in apply_boundaries + init (set implicitly via mpi_nx/mpi_ny)
    integer :: mpi_P,mpi_Q,mpi_R  ! used in exchange_flux and apply_BC + init
    integer :: numprocs  ! used in get_min_dt + init
    integer :: iam,ierr  ! used all over (wherever MPI stuff seems to be)
    integer :: cartcomm  ! used all over (wherever MPI stuff seems to be)

    integer, parameter :: NORTH = 1  ! used in exchange_flux + init
    integer, parameter :: SOUTH = 2  ! used in exchange_flux + init
    integer, parameter :: EAST  = 3  ! used in exchange_flux + init
    integer, parameter :: WEST  = 4  ! used in exchange_flux + init
    integer, parameter :: UP    = 5  ! used in exchange_flux + init
    integer, parameter :: DOWN  = 6  ! used in exchange_flux + init
    integer, parameter :: MPI_TT = MPI_REAL4  ! used in exchange_flux and get_min_dt + init

    integer nbrs(6)  ! only used in init + exchange_flux for MPI things
    integer reqs(4),stats(MPI_STATUS_SIZE,4)  ! only used in get_min_dt + exchange_flux for MPI things
    !===========================================================================

contains

    !-----------------------------------------------------------
    !   Return the x coordinate of (the center of) cell i
    !     Note: based on the location of this MPI domain (loc_lxd)
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

end module parameters
