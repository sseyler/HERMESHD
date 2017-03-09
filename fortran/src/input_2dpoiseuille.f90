!****** INPUT.F90 ************************************************************************
module input

    ! Physical system dimensions
    real, parameter :: lx = 2.2e0 !1.0e6 !4.1e2
    real, parameter :: ly = 4.1e-1 !1.0e6 !4.1e2
    real, parameter :: lz = ly/120. !1.0e6/120. !4.1e2/120.

    ! Number of Gaussian quadrature points per spatial dimension
    integer, parameter :: iquad   = 2
    integer, parameter :: nbasis  = 8 ! 8
    integer, parameter :: nbastot = 27 !30  ! TODO: only used in setup + innerintegral()

    ! Grid cell dimensions per MPI domain
    integer, parameter :: nx = 24  ! 55 (mpi_nx = 8)
    integer, parameter :: ny = 18  ! 41 (mpi_ny = 2)
    integer, parameter :: nz = 1

    ! Set number of MPI domains per spatial dimension
    integer :: mpi_nx = 8
    integer :: mpi_ny = 2

    ! Temporal integration order
    !   * 2 or 'heun' for 2nd-order RK
    !   * 3 or 'shu-osher' for 3rd-order RK
    integer, parameter :: iorder = 2
    character(*), parameter :: iname = 'heun'

    ! Fluctuating hydrodynamics
    logical, parameter :: llns = .false.

    ! Initial conditions
    ! character(*), parameter :: icname = ''
    integer, parameter :: icid = 5

    ! Boundary conditions
    character(*), parameter :: xlobc = 'wall' !'periodic'
    character(*), parameter :: xhibc = 'outflow' !'periodic'
    character(*), parameter :: ylobc = 'noslip'  !'periodic'
    character(*), parameter :: yhibc = 'noslip'  !'periodic'
    character(*), parameter :: zlobc = 'periodic'
    character(*), parameter :: zhibc = 'periodic'

    ! Simulation time
    real, parameter :: tf = 2.0e-2 !8.5e4

    ! Console output frequency
    integer, parameter :: ntout = 200

    ! Riemann solver
    ! If all of them = 0, then LLF is used for fluxes.
    !   LLF is very diffusive for the hydro problem. Roe and HLLC are much less
    !   diffusive than LLF and give very similar results with similar cpu overhead
    !   Only HLLC is setup to handle water EOS (ieos = 2)
    logical, parameter :: ihllc = .true.

    ! Equation of state:
    !   * 1 for ideal gas
    !   * 2 for Murnaghan-Tait for water
    integer, parameter :: ieos = 1

    ! Thermodynamic, constitutive, and transport parameters
    real, parameter :: TK     = 300.0     ! in Kelvin
    real, parameter :: mu     = 18.0     ! AMU per molecule
    real, parameter :: aindex = 5./3.    ! adiabatic index (gamma)
    real, parameter :: clt    = 2.0      ! numerical speed of sound

    ! Viscosity control
    !   * 0 for explicit integration
    !   * 1 for semi-implicit integration of stress terms
    !   * 2 for full 10-moment formulation (NOTE: not finished!)
    integer, parameter :: ivis = 1
    real, parameter    :: vis  = 1.0e-3  !1.79e-6  ! dynamic viscosity (air @ 15 C)
    real, parameter    :: epsi = 5.0      ! inverse relaxation coefficient

    ! Output control: location/naming and VTK output
    character (*), parameter :: datadir = "data"
    character (*), parameter :: outname = "test_poi2d_cyl_temp4"
    character (*), parameter :: outdir  = trim(datadir//"/"//outname)

    integer, parameter :: nstdout  = ntout ! density
    integer, parameter :: nstldout = 0     ! log density
    integer, parameter :: nstvout  = ntout ! velocity
    integer, parameter :: nsttout  = ntout ! temperature
    integer, parameter :: nsteiout = 0     ! internal energy density
    integer, parameter :: nstenout = 0     ! total energy density
    integer, parameter :: nstesout = 0     ! entropy density
    integer, parameter :: nstpout  = ntout ! pressure
    integer, parameter :: nststout = 0     ! stress components
    integer, parameter :: nstvrout = 0     ! vorticity

    ! Checkpointing
    !   set iread to 1 or 2 (when using the odd/even scheme)
    integer, parameter :: iread  = 0
    integer, parameter :: iwrite = 0
    logical, parameter :: resuming = .false.
    character (4), parameter :: fpre = 'Qout'

    ! Basis function testing
    logical, parameter :: test_basis = .false.

end module input
