module spatial

    implicit none
    private
    public :: grid3D, init_grid3D

    type grid3D
        real :: lx, ly, ly, dx, dy, dz
        integer :: xlo, xhi, ylo, yhi, zlo, zhi
        integer :: nx, ny, nz, nb

        ! dimension(nx,ny,nbasis)
        real, allocatable, dimension(:,:,:,:) :: Qrh
        real, allocatable, dimension(:,:,:,:) :: Qmx
        real, allocatable, dimension(:,:,:,:) :: Qmy
        real, allocatable, dimension(:,:,:,:) :: Qmz
        real, allocatable, dimension(:,:,:,:) :: Qen

        real, allocatable, dimension(:,:,:,:) :: Qpxx
        real, allocatable, dimension(:,:,:,:) :: Qpyy
        real, allocatable, dimension(:,:,:,:) :: Qpzz
        real, allocatable, dimension(:,:,:,:) :: Qpxy
        real, allocatable, dimension(:,:,:,:) :: Qpxz
        real, allocatable, dimension(:,:,:,:) :: Qpyz

    end type grid3D

contains

    subroutine init_grid3D(grid, nx, ny, nz, nb, vis)
        type(grid3D), intent(inout) :: grid
        integer, intent(in) :: nx, ny, nz, nb
        boolean, intent(in) :: vis

        grid%nx = nx
        grid%ny = ny
        grid%nz = nz
        grid%nb = nb

        allocate(grid%Qrh(nx,ny,nz,nb))
        allocate(grid%Qmx(nx,ny,nz,nb))
        allocate(grid%Qmy(nx,ny,nz,nb))
        allocate(grid%Qmz(nx,ny,nz,nb))
        allocate(grid%Qen(nx,ny,nz,nb))

        if (vis) then
            allocate(grid%Qpxx(nx,ny,nz,nb))
            allocate(grid%Qpyy(nx,ny,nz,nb))
            allocate(grid%Qpzz(nx,ny,nz,nb))
            allocate(grid%Qpxy(nx,ny,nz,nb))
            allocate(grid%Qpxz(nx,ny,nz,nb))
            allocate(grid%Qpyz(nx,ny,nz,nb))
        endif

        lxd = -(lx/2.0)
        lxu =  (lx/2.0)
        lyd = -(ly/2.0)
        lyu =  (ly/2.0)
        lzd = -(lz/2.0)
        lzu =  (lz/2.0)

        dxi = (nx*mpi_nx)/(lxu-lxd)
        dyi = (ny*mpi_ny)/(lyu-lyd)
        dzi = (nz*mpi_nz)/(lzu-lzd)
        dx = 1./dxi
        dy = 1./dyi
        dz = 1./dzi

        loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
        loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
        loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz
    end subroutine init_grid3D


    real function xc(i)
        integer i
        xc = loc_lxd + (i - 0.5)*dx
    end function xc

    real function yc(j)
        integer j
        yc = loc_lyd + (j - 0.5)*dy
    end function yc

    real function zc(k)
        integer k
        zc = loc_lzd + (k - 0.5)*dz
    end function zc


end module spatial
