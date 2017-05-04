!***** RESERVOIR.F90 **********************************************************************

!##### reservoir ###########################################################
module reservoir

use params
use spatial

implicit none

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutines for setting BCs
    !------------------------------------------------------------
    abstract interface
        subroutine xres_fcn_ptr(Qxbc_i, Qxbc_e)
            use input, only : ny,nz
            use params, only : mx,my,mz,nface,nQ
            real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_i
            real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_e
        end subroutine xres_fcn_ptr

        ! subroutine ybc_fcn_ptr(Qybc_i, Qybc_e)
        !     use input, only : nx,nz
        !     use params, only : mx,my,mz,nface,nQ
        !     real, dimension(nx,nz,nface,nQ), intent(in) :: Qybc_i
        !     real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_e
        ! end subroutine ybc_fcn_ptr
        !
        ! subroutine zbc_fcn_ptr(Qzbc_i, Qzbc_e)
        !     use input, only : nx,ny
        !     use params, only : mx,my,mz,nface,nQ
        !     real, dimension(nx,ny,nface,nQ), intent(in) :: Qzbc_i
        !     real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_e
        ! end subroutine zbc_fcn_ptr
    end interface
    !---------------------------------------------------------------------------

contains

    !===========================================================================
    ! PERIODIC BC subroutines
    !------------------------------------------------------------
    subroutine xlobc_periodic(Qxlo_i, Qxlo_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxlo_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e
    end subroutine xlobc_periodic

    subroutine xhibc_periodic(Qxhi_i, Qxhi_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxhi_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e
    end subroutine xhibc_periodic
    !--------------------------------------------------------
    subroutine ylobc_periodic(Qylo_i, Qylo_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qylo_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e
    end subroutine ylobc_periodic

    subroutine yhibc_periodic(Qyhi_i, Qyhi_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qyhi_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e
    end subroutine yhibc_periodic
    !--------------------------------------------------------
    subroutine zlobc_periodic(Qzlo_i, Qzlo_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzlo_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e
    end subroutine zlobc_periodic

    subroutine zhibc_periodic(Qzhi_i, Qzhi_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzhi_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e
    end subroutine zhibc_periodic
    !--------------------------------------------------------


    !===========================================================================
    ! WALL BC subroutine
    !------------------------------------------------------------
    subroutine xlobc_wall(Qxlo_i, Qxlo_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxlo_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e
        integer j,k
        do k = 1,nz
        do j = 1,ny
            Qxlo_e(j,k,1:nface,:) = Qxlo_i(j,k,1:nface,:)
            Qxlo_e(j,k,1:nface,mx) = 0.0
        end do
        end do
    end subroutine xlobc_wall
    !--------------------------------------------------------
    subroutine xhibc_wall(Qxhi_i, Qxhi_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxhi_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e
        integer j,k
        do k = 1,nz
        do j = 1,ny
            Qxhi_e(j,k,1:nface,:) = Qxhi_i(j,k,1:nface,:)
            Qxhi_e(j,k,1:nface,mx) = 0.0
        end do
        end do
    end subroutine xhibc_wall
    !--------------------------------------------------------
    subroutine ylobc_wall(Qylo_i, Qylo_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qylo_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e
        integer i,k
        do k = 1,nz
        do i = 1,nx
            Qylo_e(i,k,1:nface,:) = Qylo_i(i,k,1:nface,:)
            Qylo_e(i,k,1:nface,my) = 0.0
        end do
        end do
    end subroutine ylobc_wall

    subroutine yhibc_wall(Qyhi_i, Qyhi_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qyhi_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e
        integer i,k
        do k = 1,nz
        do i = 1,nx
            Qyhi_e(i,k,1:nface,:) = Qyhi_i(i,k,1:nface,:)
            Qyhi_e(i,k,1:nface,my) = 0.0
        end do
        end do
    end subroutine yhibc_wall
    !--------------------------------------------------------
    subroutine zlobc_wall(Qzlo_i, Qzlo_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzlo_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e
        integer i,j
        do j = 1,ny
        do i = 1,nx
            Qzlo_e(i,j,1:nface,:) = Qzlo_i(i,j,1:nface,:)
            Qzlo_e(i,j,1:nface,mz) = 0.0
        end do
        end do
    end subroutine zlobc_wall
    !--------------------------------------------------------
    subroutine zhibc_wall(Qzhi_i, Qzhi_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzhi_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e
        integer i,j
        do j = 1,ny
        do i = 1,nx
            Qzhi_e(i,j,1:nface,:) = Qzhi_i(i,j,1:nface,:)
            Qzhi_e(i,j,1:nface,mz) = 0.0
        end do
        end do
    end subroutine zhibc_wall
    !---------------------------------------------------------------------------

end module reservoir



!##### boundary ################################################################
module boundary

    use params
    use helpers, only : mpi_print
    use boundary_defs
    use boundary_custom

    !===========================================================================
    ! Initialize pointers to subroutines for setting BCs
    !------------------------------------------------------------
    procedure(xbc_fcn_ptr), pointer :: apply_xlobc => null ()
    procedure(xbc_fcn_ptr), pointer :: apply_xhibc => null ()
    procedure(ybc_fcn_ptr), pointer :: apply_ylobc => null ()
    procedure(ybc_fcn_ptr), pointer :: apply_yhibc => null ()
    procedure(zbc_fcn_ptr), pointer :: apply_zlobc => null ()
    procedure(zbc_fcn_ptr), pointer :: apply_zhibc => null ()
    !---------------------------------------------------------------------------

    real :: rnx = rx


    !===========================================================================
    ! Initialize arrays to store reservoir state variables
    !------------------------------------------------------------
    real, dimension(rnx,ny,nz,nface,nQ) :: Qxrlo, Qxrhi

    !---------------------------------------------------------------------------

contains

    !===========================================================================
    ! Apply boundary conditions (specified by user at runtime)
    !------------------------------------------------------------
    subroutine apply_boundaries
        real ux_scale
        if ( mpi_P == 1 ) then
            call apply_xlobc(Qxlo_int, Qxlo_ext)
        end if
        if ( mpi_P == mpi_nx ) then
            call apply_xhibc(Qxhi_int, Qxhi_ext)
        end if
        !----------------------------------------------------
        if ( mpi_Q == 1 ) then
            call apply_ylobc(Qylo_int, Qylo_ext)
        end if
        if ( mpi_Q == mpi_ny ) then
            call apply_yhibc(Qyhi_int, Qyhi_ext)
        end if
        !--------------------------------------------------------
        if ( mpi_R == 1 ) then
            call apply_zlobc(Qzlo_int, Qzlo_ext)
        end if
        if ( mpi_R == mpi_nz ) then
            call apply_zhibc(Qzhi_int, Qzhi_ext)
        end if
    end subroutine apply_boundaries
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Select user-specified BC (subroutines) at runtime
    !------------------------------------------------------------
    subroutine select_x_boundaries(xlobc_flag, xhibc_flag, xlobc_ptr, xhibc_ptr)
        character(*), intent(in) :: xlobc_flag, xhibc_flag
        procedure(xbc_fcn_ptr), pointer :: xlobc_ptr, xhibc_ptr

        select case (xlobc_flag)
            case ('periodic')
                xlobc_ptr => xlobc_periodic
            case ('outflow')
                xlobc_ptr => xlobc_outflow
            case ('wall')
                xlobc_ptr => xlobc_wall
            case ('noslip')
                xlobc_ptr => xlobc_noslip
            case default
                xlobc_ptr => xlobc_periodic
        end select
        select case (xhibc_flag)
            case ('periodic')
                xhibc_ptr => xhibc_periodic
            case ('outflow')
                xhibc_ptr => xhibc_outflow
            case ('wall')
                xhibc_ptr => xhibc_wall
            case ('noslip')
                xhibc_ptr => xhibc_noslip
            case default
                xhibc_ptr => xhibc_periodic
        end select

    end subroutine select_x_boundaries
    !--------------------------------------------------------


end module boundary
