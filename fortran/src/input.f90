module input

    ! Grid cell dimensions per MPI domain
    integer, parameter :: nx = 40
    integer, parameter :: ny = 1
    integer, parameter :: nz = 1

    ! Set number of MPI domains per spatial dimension
    integer :: mpi_nx = 16
    integer :: mpi_ny = 1

    ! Physical system dimensions
    real, parameter :: lx = 1.0e6
    real, parameter :: ly = 1.0e6
    real, parameter :: lz = 1.0e6/120.

    ! Number of Gaussian quadrature points per spatial dimension
    integer, parameter :: iquad = 2

    ! Temporal integration order
    !   * 2 for
    !   * 3 for 3rd-order Runga-Kutta (Shu-Osher)
    integer, parameter :: iorder = 2

    ! Fluctuating hydrodynamics
    logical, parameter :: llns = .false.

    ! Initial conditions
    integer, parameter :: icid = 3

    ! Boundary conditions:
    !   * 0 for set_bc subroutine used to prescribe BCs
    !   * 1 for wall (vanishing normal velocities).
    !   * 2 for periodic (MPI does this for you).
    integer, parameter :: xlbc = 1, xhbc = 1
    integer, parameter :: ylbc = 2, yhbc = 2
    integer, parameter :: zlbc = 2, zhbc = 2

    ! Simulation time
    real, parameter :: tf = 3.3e4

    ! Riemann solver
    integer, parameter :: ihllc = 1, iroe = 0, ieos = 1

    ! Flow parameters
    real, parameter :: vis = 0.0
    real, parameter :: epsi = 5.0
    real, parameter :: clt = 2.0

    ! Thermodynamic and transport parameters
    real, parameter :: mu = 2.0
    real, parameter :: aindex = 5./3.
    real, parameter :: aindm1 = aindex-1.0
    real, parameter :: cp = aindex/aindm1

    ! Output frequency and directory
    integer, parameter :: ntout = 100
    character (*), parameter :: datadir="data"
    character (*), parameter :: outname="sod_1d_1"
    character (*), parameter :: outdir = trim(datadir//"/"//outname)

    ! Checkpointing
    !   set iread to 1 or 2 (when using the odd/even scheme)
    integer, parameter :: iread = 0, iwrite = 0
    character (4), parameter :: fpre = 'Qout'
    logical, parameter :: resuming = .false.

end module input
