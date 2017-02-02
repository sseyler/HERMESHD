module input

    ! Physical system dimensions
    real, parameter :: lx = 4.1e2 !4.1e2
    real, parameter :: ly = 4.1e2 !4.1e2
    real, parameter :: lz = 4.1e2/120. !4.1e2/120.

    ! Number of Gaussian quadrature points per spatial dimension
    integer, parameter :: iquad = 2

    ! Grid cell dimensions per MPI domain
    integer, parameter :: nx = 40
    integer, parameter :: ny = 40
    integer, parameter :: nz = 1

    ! Set number of MPI domains per spatial dimension
    integer :: mpi_nx = 4
    integer :: mpi_ny = 4

    ! Temporal integration order
    !   * 2 or 'heun' for 2nd-order RK
    !   * 3 or 'shu-osher' for 3rd-order RK
    integer, parameter :: iorder = 2
    character(*), parameter :: iname = 'heun'

    ! Fluctuating hydrodynamics
    logical, parameter :: llns = .false.

    ! Initial conditions
    integer, parameter :: icid = 1
    ! character(*), parameter :: icname = ''

    ! Boundary conditions:
    !   * 0 for periodic (MPI does this for you)
    !   * 1 for outflow (MPI does this for you)
    !   * 2 for wall (vanishing normal velocities)
    !   * 3 for no-slip (vanishing normal/tangential velocities)
    character(*), parameter :: xlobc = 'periodic' !'periodic'
    character(*), parameter :: xhibc = 'periodic' !'periodic'
    character(*), parameter :: ylobc = 'periodic'  !'periodic'
    character(*), parameter :: yhibc = 'periodic'  !'periodic'
    character(*), parameter :: zlobc = 'periodic'
    character(*), parameter :: zhibc = 'periodic'

    ! Simulation time
    real, parameter :: tf = 5.0e3

    ! Console output frequency
    integer, parameter :: ntout = 200

    ! Riemann solver
    ! If all of them = 0, then LLF is used for fluxes.
    !   LLF is very diffusive for the hydro problem. Roe and HLLC are much less
    !   diffusive than LLF and give very similar results with similar cpu overhead
    !   Only HLLC is setup to handle water EOS (ieos = 2)
    logical, parameter :: ihllc = .true.

    ! Thermodynamic and transport parameters
    real, parameter :: ieos = 1
    real, parameter :: mu = 2.0
    real, parameter :: aindex = 5./3.
    real, parameter :: aindm1 = aindex - 1.0
    real, parameter :: cp = aindex/aindm1

    ! Equation of state and constitutive parameters
    real, parameter :: vis = 1.0e-2
    real, parameter :: epsi = 5.0
    real, parameter :: clt = 2.0

    ! Output location and naming
    character (*), parameter :: datadir="data"
    character (*), parameter :: outname="test_modbc_pipe_0"
    character (*), parameter :: outdir = trim(datadir//"/"//outname)

    ! Checkpointing
    !   set iread to 1 or 2 (when using the odd/even scheme)
    integer, parameter :: iread  = 0
    integer, parameter :: iwrite = 0
    logical, parameter :: resuming = .false.
    character (4), parameter :: fpre = 'Qout'

end module input
