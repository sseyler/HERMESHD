!****** SOD_INPUT.F90 ************************************************************************
module input

    ! Physical system dimensions
    real, parameter :: lx = 1.0e1 !1.0e6 !4.1e2
    real, parameter :: ly = lx/12. !1.0e6 !4.1e2
    real, parameter :: lz = lx/120. !1.0e6/120. !4.1e2/120.

    ! Number of Gaussian quadrature points per spatial dimension
    integer, parameter :: iquad = 2
    integer, parameter :: nbasis  = 8 ! 8
    integer, parameter :: nbastot = 27 !30  ! TODO: only used in setup + innerintegral()

    ! Grid cell dimensions per MPI domain
    integer, parameter :: nx = 50
    integer, parameter :: ny = 1
    integer, parameter :: nz = 1

    ! Set number of MPI domains per spatial dimension
    integer :: mpi_nx = 16
    integer :: mpi_ny = 1

    ! Temporal integration order
    !   * 2 or 'heun' for 2nd-order RK
    !   * 3 or 'shu-osher' for 3rd-order RK
    integer, parameter :: iorder = 2
    character(*), parameter :: iname = 'heun'

    ! Fluctuating hydrodynamics
    logical, parameter :: llns = .false.

    ! Initial conditions
    ! character(*), parameter :: icname = ''
    integer, parameter :: icid = 3

    ! Boundary conditions
    character(*), parameter :: xlobc = 'wall'
    character(*), parameter :: xhibc = 'wall'
    character(*), parameter :: ylobc = 'periodic'
    character(*), parameter :: yhibc = 'periodic'
    character(*), parameter :: zlobc = 'periodic'
    character(*), parameter :: zhibc = 'periodic'

    ! Simulation time
    real, parameter :: tf = 8.5e-3 !8.5e4

    ! Console output frequency
    integer, parameter :: ntout = 100

    ! Riemann solver
    logical, parameter :: ihllc = .true.

    ! Equation of state:
    !   * 1 for ideal gas
    !   * 2 for Murnaghan-Tait for water
    integer, parameter :: ieos = 1

    ! Thermodynamic, constitutive, and transport parameters
    real, parameter :: te     = 300.0    ! in Kelvin
    real, parameter :: mu     = 18.0     ! AMU per molecule
    real, parameter :: aindex = 5./3.    ! adiabatic index (gamma)
    real, parameter :: clt    = 2.0      ! numerical speed of sound

    ! Viscosity control
    !   * 0 for explicit integration
    !   * 1 for semi-implicit integration of stress terms
    !   * 2 for full 10-moment formulation (NOTE: not finished!)
    integer, parameter :: ivis = 1
    real, parameter    :: vis  = 1.0e-3  ! dynamic viscosity
    real, parameter    :: epsi = 5.0     ! inverse relaxation coefficient

    ! Output control: location/naming and VTK output
    character (*), parameter :: datadir = "data"
    character (*), parameter :: outname = "test_sod_0"
    character (*), parameter :: outdir  = trim(datadir//"/"//outname)

    logical, parameter :: o_density     = .true.
    logical, parameter :: o_logdensity  = .false.
    logical, parameter :: o_velocity    = .true.
    logical, parameter :: o_temperature = .true.
    logical, parameter :: o_entropy     = .false.
    logical, parameter :: o_pressure    = .true.
    logical, parameter :: o_stress      = .false.
    logical, parameter :: o_vorticity   = .false.

    ! Checkpointing
    !   set iread to 1 or 2 (when using the odd/even scheme)
    integer, parameter :: iread  = 0
    integer, parameter :: iwrite = 0
    logical, parameter :: resuming = .false.
    character (4), parameter :: fpre = 'Qout'

    ! Basis function testing
    logical, parameter :: test_basis = .false.

end module input
