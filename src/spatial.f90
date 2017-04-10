module spatial

    implicit none

    type grid2D
        real :: lx, ly, dx, dy
        integer :: xlo, xhi, ylo, yhi

        ! dimension(nx,ny,nbasis)
        real, allocatable, dimension(:,:,:) :: Qrh
        real, allocatable, dimension(:,:,:) :: Qmx
        real, allocatable, dimension(:,:,:) :: Qmy
        real, allocatable, dimension(:,:,:) :: Qmz
        real, allocatable, dimension(:,:,:) :: Qen

        character(len=32) :: xlbc, xhbc, ylbc, yhbc
    end type grid2D

    public :: grid2D

contains

    subroutine init_grid(grid, nx, ny, nbasis)
        type(grid2D), intent(inout) :: grid
        integer, intent(in) :: nx, ny, nbasis
    end subroutine init_grid


end module spatial
