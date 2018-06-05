module spatial

use input
use params

implicit none

! private
! public :: grid3D, init_grid3D

!===========================================================================
! Arrays for field variables
!------------------------------------------------------------
! real, allocatable, dimension(:,:,:,:,:) :: Q_r0, Q_r1, Q_r2, Q_r3
real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r0, Q_r1, Q_r2, Q_r3
!===========================================================================

real :: dz,dy,dx, dxi,dyi,dzi, dvi  ! used throughout + directly by grid coord functions
! real lxd,lxu,lyd,lyu,lzd,lzu  ! used in init + indirectly used by the grid coord functions
real :: loc_lxd,loc_lyd,loc_lzd  ! used directly by the grid coord functions

real :: dxid4,dyid4,dzid4

! type grid3D
!     real :: lx, ly, ly, dx, dy, dz
!     integer :: xlo, xhi, ylo, yhi, zlo, zhi
!     integer :: nx, ny, nz, nb
!
!     ! dimension(nx,ny,nbasis)
!     real, allocatable, dimension(:,:,:,:) :: Qrh
!     real, allocatable, dimension(:,:,:,:) :: Qmx
!     real, allocatable, dimension(:,:,:,:) :: Qmy
!     real, allocatable, dimension(:,:,:,:) :: Qmz
!     real, allocatable, dimension(:,:,:,:) :: Qen
!
!     real, allocatable, dimension(:,:,:,:) :: Qexx
!     real, allocatable, dimension(:,:,:,:) :: Qeyy
!     real, allocatable, dimension(:,:,:,:) :: Qezz
!     real, allocatable, dimension(:,:,:,:) :: Qexy
!     real, allocatable, dimension(:,:,:,:) :: Qexz
!     real, allocatable, dimension(:,:,:,:) :: Qeyz
!
! end type grid3D

contains

    ! subroutine init_grid3D(grid, nx, ny, nz, nb, vis)
    !     type(grid3D), intent(inout) :: grid
    !     integer, intent(in) :: nx, ny, nz, nb
    !     boolean, intent(in) :: vis
    !
    !     grid%nx = nx
    !     grid%ny = ny
    !     grid%nz = nz
    !     grid%nb = nb
    !
    !     allocate(grid%Qrh(nx,ny,nz,nb))
    !     allocate(grid%Qmx(nx,ny,nz,nb))
    !     allocate(grid%Qmy(nx,ny,nz,nb))
    !     allocate(grid%Qmz(nx,ny,nz,nb))
    !     allocate(grid%Qen(nx,ny,nz,nb))
    !
    !     if (vis) then
    !         allocate(grid%Qexx(nx,ny,nz,nb))
    !         allocate(grid%Qeyy(nx,ny,nz,nb))
    !         allocate(grid%Qezz(nx,ny,nz,nb))
    !         allocate(grid%Qexy(nx,ny,nz,nb))
    !         allocate(grid%Qexz(nx,ny,nz,nb))
    !         allocate(grid%Qeyz(nx,ny,nz,nb))
    !     endif
    !
    !     lxd = -(lx/2.0)
    !     lxu =  (lx/2.0)
    !     lyd = -(ly/2.0)
    !     lyu =  (ly/2.0)
    !     lzd = -(lz/2.0)
    !     lzu =  (lz/2.0)
    !
    !     dxi = (nx*mpi_nx)/(lxu-lxd)
    !     dyi = (ny*mpi_ny)/(lyu-lyd)
    !     dzi = (nz*mpi_nz)/(lzu-lzd)
    !     dx = 1./dxi
    !     dy = 1./dyi
    !     dz = 1./dzi
    !
    !     loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
    !     loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
    !     loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz
    ! end subroutine init_grid3D


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

    real function rz(i,j)
        integer i,j
        rz = sqrt(yc(j)**2)
    end function rz

    real function r(i,j)
        integer i,j,k
        r = sqrt(xc(i)**2 + yc(j)**2)
    end function r

    real function theta(i,j)
        integer i,j
        theta = atan2(yc(j),xc(i))
    end function theta


end module spatial
