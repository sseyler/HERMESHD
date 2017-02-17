module input

    ! Physical system dimensions
    real, parameter :: lx = 2.2e0 !1.0e6 !4.1e2
    real, parameter :: ly = 8.2e-1 !1.0e6 !4.1e2
    real, parameter :: lz = ly/120. !1.0e6/120. !4.1e2/120.

    ! Number of Gaussian quadrature points per spatial dimension
    integer, parameter :: iquad  = 2

    ! Grid cell dimensions per MPI domain
    integer, parameter :: nx = 20
    integer, parameter :: ny = 30
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
    integer, parameter :: icid = 5
    ! character(*), parameter :: icname = ''

    ! Boundary conditions:
    !   * 0 for periodic (MPI does this for you)
    !   * 1 for outflow (MPI does this for you)
    !   * 2 for wall (vanishing normal velocities)
    !   * 3 for noslip (vanishing normal/tangential velocities)
    character(*), parameter :: xlobc = 'noslip' !'periodic'
    character(*), parameter :: xhibc = 'outflow' !'periodic'
    character(*), parameter :: ylobc = 'noslip'  !'periodic'
    character(*), parameter :: yhibc = 'noslip'  !'periodic'
    character(*), parameter :: zlobc = 'periodic'
    character(*), parameter :: zhibc = 'periodic'

    ! Simulation time
    real, parameter :: tf = 1.0e0 !8.5e4

    ! Console output frequency
    integer, parameter :: ntout = 1000

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
    real, parameter :: te     = 300.0    ! in Kelvin
    real, parameter :: mu     = 18.0     ! AMU per molecule
    real, parameter :: aindex = 5./3.    ! adiabatic index (gamma)
    real, parameter :: clt    = 2.0      ! numerical speed of sound

    ! Viscosity control
    integer, parameter :: ivis = 0       ! 0 for explicit, 1 for semi-implicit
    real, parameter    :: vis  = 1.0e-2  ! dynamic viscosity
    real, parameter    :: epsi = 5.0     ! inverse relaxation coefficient

    ! Output location and naming
    character (*), parameter :: datadir = "data"
    character (*), parameter :: outname = "test_modbc_pipe_2"
    character (*), parameter :: outdir  = trim(datadir//"/"//outname)

    ! Checkpointing
    !   set iread to 1 or 2 (when using the odd/even scheme)
    integer, parameter :: iread  = 0
    integer, parameter :: iwrite = 0
    logical, parameter :: resuming = .false.
    character (4), parameter :: fpre = 'Qout'

end module input
