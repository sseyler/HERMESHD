!****** SOD_INPUT.F90 ************************************************************************
module input

    ! Number of Gaussian quadrature points per spatial dimension
    integer, parameter :: iquad  = 2
    integer, parameter :: nbasis = 8

    ! Grid cell dimensions per MPI domain
    integer, parameter :: nx = 24
    integer, parameter :: ny = 30
    integer, parameter :: nz = 1

    ! Set number of MPI domains per spatial dimension
    integer, parameter :: mpi_nx = 5
    integer, parameter :: mpi_ny = 4

    ! Physical system dimensions
    real, parameter :: Lbox = 2.0e3
    real, parameter :: lxu =  Lbox/2.
    real, parameter :: lyu =  Lbox/2.
    real, parameter :: lzu =  Lbox/2./(mpi_nx*nx)

    ! Temporal integration order
    !   * 2 or 'heun' for 2nd-order RK
    !   * 3 or 'shu-osher' for 3rd-order RK
    integer, parameter :: iorder = 2
    character(*), parameter :: iname = 'heun'

    ! Initial conditions
    ! character(*), parameter :: icname = ''
    integer, parameter :: icid = 0

    ! Boundary conditions
    character(*), parameter :: xlobc = 'periodic'
    character(*), parameter :: xhibc = 'periodic'
    character(*), parameter :: ylobc = 'periodic'
    character(*), parameter :: yhibc = 'periodic'
    character(*), parameter :: zlobc = 'periodic'
    character(*), parameter :: zhibc = 'periodic'

    ! Simulation time and console output frequency
    real, parameter :: tf = 2.5e2
    integer, parameter :: ntout = 200

    ! Fluid model: equation of state and fluctuations:
    !   * 1 for ideal gas
    !   * 2 for Murnaghan-Tait for water
    integer, parameter :: ieos = 1
    logical, parameter :: llns = .true.

    ! Viscosity control
    !   * 'linear_ex' for explicit integration of linearized 10-moment eqns
    !   * 'linear' for semi-implicit integration of linearized 10-moment eqns
    !   * 'full' for semi-implicit integration of full, nonlinear 10-moment eqns
    character(*), parameter :: ivis = 'full'
    real, parameter    :: eta  = 2.7e-4   ! dynamic viscosity
    real, parameter    :: zeta = 9.0e-5   ! bulk viscosity
    ! real, parameter    :: epsi = 5.0     ! inverse relaxation coefficient

    ! Thermodynamic, constitutive, and transport parameters
    real, parameter :: TK     = 80.0     ! in Kelvin
    real, parameter :: mu     = 40.0     ! AMU per molecule
    real, parameter :: aindex = 5./3.    ! adiabatic index (gamma)
    real, parameter :: clt    = 3.0      ! numerical speed of sound

    ! Dimensional units -- expressed in MKS. NOTE: temperature (te0) in eV!
    real, parameter :: L0 = 2.5e-9                 ! length
    real, parameter :: t0 = 1.0e-9                 ! time
    real, parameter :: n0 = 2.3e28                 ! number density

    ! Riemann solver
    logical, parameter :: ihllc = .true.

    ! Output control: location/naming and VTK output
    character(*), parameter :: basedir = "/scratch/sseyler/WORK/BWAnnualReport2018"
    character(*), parameter :: datadir = "hydrojet"
    character(*), parameter :: outname = trim(ivis//"_"//"fhd_005")
    character(*), parameter :: outdir  = trim(basedir//"/"//datadir//"/"//outname)

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
    character(4), parameter :: fpre = 'Qout'

    ! Basis function testing
    logical, parameter :: test_basis = .false.

    ! Setting initial conditions
    real, parameter :: vx_max = 1.0e0, vy_max = 1.0e0
    real, parameter :: vx_min = 0.0,   vy_min = 0.0

end module input
