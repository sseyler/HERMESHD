!****** PARAMS.F90 *******************************************************************

! TEST RESULTS SUMMARY USING VARIOUS #s OF BASIS FUNCTIONS & QUADRATURE POINTS
!
!   The nbasis = 8 and iquad = 2 case gives extremely close results to nbasis = 10
!   and iquad = 3. The (8,2) case takes 21 seconds per iteration and (10,3) takes 35
!   seconds. Close examination of the output indicates that the (10,3) case has
!   slightly better resolution/sharpness/definition as seen in the transitions from
!   high to low or low to high vorticity.  The (27,3) case has sightly better
!   resolution/sharpness/definition than the (10,3) case in about the same
!   proportion as the (10,3) case does compared to the (8,2) case.  However, the
!   time for an iteration is about 190 sec, which is about 4x higher than the (10,3)
!   case.  It would appear for this test there is no significant advantage of using
!   27 elements over the 10 elements. There of ourse could be applications in which
!   the 27 elements does provide significant improvement of accuracy. Finally in
!   comparing the cases (3,27) to (4,20) there does not appear to be any real
!   differences, only that the (4,20) case was about 15% faster.
!
! Use iquad = 2
!   nbasis = 4:  {1, x, y, z}
!   nbasis = 8:  {1, x, y, z, yz, zx, xy, xyz}
!
! Use iquad = 3
!   nbasis = 10: {1, x, y, z, yz, zx, xy, P_2(x), P_2(y), P_2(z)}
!   nbasis = 27: {1, x, y, z, yz, zx, xy, xyz, P2(x), P2(y), P2(z),
!                 yP2(z), zP2(x), xP2(y), P2(y)z, P2(z)x, P2(x)y,
!                 P2(y)P2(z), P2(z)P2(x), P2(x)P2(y), yzP2(x), zxP2(y), xyP2(z),
!                 xP2(y)P2(z), yP2(z)P2(x), zP2(x)P2(y), P2(x)P2(y)P2(z)}
!
! Use iquad = 4
!   nbasis = 20: nbasis10 + {xyz, xP2(y), yP2(x), xP2(z),
!                                 zP2(x), yP2(z), zP2(y), P3(x), P3(y), P3(z)}
!*******************************************************************************
module params

    use input

    ! NOTE: "mpif.h" is in /nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include
    ! include '/nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include/mpif.h'
    include 'mpif.h'

    integer, parameter :: rh = 1                    ! density
    integer, parameter :: mx = 2, my = 3,  mz = 4   ! vector momentum
    integer, parameter :: en = 5                    ! scalar energy
    integer, parameter :: exx= 6, eyy= 7,  ezz= 8   ! isotropic stress
    integer, parameter :: exy= 9, exz= 10, eyz= 11  ! deviatoric stress
    integer, parameter :: nQ = 11                   ! number of field variables

    integer, parameter :: nbastot = 35
    integer, parameter :: ngu = 0  ! TODO: only used in output_vtk()

    ! iquad: # of Gaussian quadrature points per direction. iquad should not be:
    !   < ipoly (max Legendre polynomial order used) --> unstable
    !   > ipoly+1 --> an exact Gaussian quadrature for Legendre poly used
    ! Thus there are only two cases of iquad for a given nbasis. Both give similar
    ! results although iquad = ipoly + 1 is formally more accurate.
    integer, parameter :: nedge = iquad
    integer, parameter :: nface = iquad*iquad    ! # of quadrature points per cell face
    integer, parameter :: npg   = nface*iquad    ! # of internal points per cell
    integer, parameter :: nfe   = 2*nface        ! # of cell face quad points per direction
    integer, parameter :: npge  = 6*nface        ! total # of cell face quad points
    integer, parameter :: nslim = npg + 6*nface  ! total # of quad points per cell

    integer, parameter :: nface1 = nface + 1     ! Compact shortcut for nface + 1
    ! Grid cell dimensions parameters (mostly for Sflux_xyz)
    integer, parameter :: nx1 = nx + 1
    integer, parameter :: ny1 = ny + 1
    integer, parameter :: nz1 = nz + 1
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Constants, and physical and numerical parameters
    !------------------------------------------------------------
    ! Physical system dimensions
    real, parameter :: lxd = -lxu
    real, parameter :: lyd = -lyu
    real, parameter :: lzd = -lzu
    real, parameter :: lx = lxu-lxd
    real, parameter :: ly = lyu-lyd
    real, parameter :: lz = lzu-lzd

    ! Useful constants
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: sqrt2 = 2**0.5, sqrt2i = 1.0/sqrt2
    real, parameter :: c1d5 = 1./5.
    real, parameter :: c1d3 = 1./3., c2d3 = 2./3., c4d3 = 4./3.
    real, parameter :: eV_per_K = 8.61728e-5

    ! Useful derived parameters
    real, parameter :: aindm1 = aindex - 1. ! gamma - 1
    real, parameter :: cp = aindex/aindm1   ! specific heat at constant pressure

    ! Derived units
    real, parameter :: v0  = L0/t0                 ! velocity
    real, parameter :: p0  = mu*1.67e-27*n0*v0**2  ! pressure
    real, parameter :: te0 = p0/n0/1.6e-19         ! temperature (eV, not K!)

    ! rh_min is a min density to be used for ideal gas EOS, rh_min is min density
    ! below which the pressure becomes negative for the MT water EOS.
    ! The DG-based subroutine "limiter" keeps density above rh_mult*rh_min.
    real, parameter :: rh_floor = 1.0e-3     ! 1.0e-1 was old value
    real, parameter :: T_floor  = TK*eV_per_K/te0  ! 0.02585 eV ~ 300 K
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
    !   coleta = coll*eta = density*temp  (eta is dynamic viscosity, or sometimes mu)
    real :: nu, nueta, c2d3ne, c4d3ne
    real :: coll, coleta   ! values set in set_ic
    real :: c2d3cv, c4d3cv ! probably want to set in set_ic (not currently used)

    ! NOTE: Old way of handling relaxation system
    ! real, parameter :: nu = epsi*vis
    ! real, parameter :: c2d3nu=c2d3*nu, c4d3nu=c4d3*nu
    real, parameter :: gr         = 0 !9.807e3*t0**2/L0  ! DEBUG: gravitational accel
    real, parameter :: T_base     = TK*eV_per_K/te0 ! te*eV_per_K/te0  ! temp (isothermal)
    real, parameter :: eta_base   = eta/(p0*t0)     ! dynamic viscosity
    real, parameter :: zeta_base  = zeta/(p0*t0)    ! bulk viscosity---will need to adjust this!
    real, parameter :: kappa_base = 1.e-1

    ! NOTE: factors of sqrt(2) below are due to temporal(?) averaging; isn't this
    !       just due to averaging two random samples on both sides of the face of
    !       a grid cell? (Decreases std by sqrt(2), so multiply to correct it.)
    real, parameter :: eta_sd   = (2  *  eta_base*T_base)**0.5    ! sqrt2 fluc std for shear visc
    real, parameter :: zeta_sd  = (c1d3*zeta_base*T_base)**0.5    ! sqrt2 fluc std for bulk visc
    real, parameter :: kappa_sd = (2 * kappa_base*T_base**2)**0.5 ! sqrt2    !===========================================================================


    !===========================================================================
    ! Time(s)
    !------------------------------------------------------------
    real t, dt, dtout ! Inv sq-root of (dV*dt), dV = grid cell volume
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

end module params
