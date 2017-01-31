module boundaries

    use parameters
    use helpers, only : mpi_print

    abstract interface
        subroutine xbc_fcn_ptr(Qxbc_e, Qxbc_i)
            use input, only : ny,nz
            use parameters, only : mx,nface,nQ
            real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_e, Qxbc_i
        end subroutine xbc_fcn_ptr
    end interface

    abstract interface
        subroutine ybc_fcn_ptr(Qybc_e, Qybc_i)
            use input, only : nx,nz
            use parameters, only : my,nface,nQ
            real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_e, Qybc_i
        end subroutine ybc_fcn_ptr
    end interface

    abstract interface
        subroutine zbc_fcn_ptr(Qzbc_e, Qzbc_i)
            use input, only : nx,ny
            use parameters, only : mz,nface,nQ
            real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_e, Qzbc_i
        end subroutine zbc_fcn_ptr
    end interface

    procedure(xbc_fcn_ptr), pointer :: set_xlbc => null ()
    procedure(xbc_fcn_ptr), pointer :: set_xhbc => null ()
    procedure(ybc_fcn_ptr), pointer :: set_ylbc => null ()
    procedure(ybc_fcn_ptr), pointer :: set_yhbc => null ()
    procedure(zbc_fcn_ptr), pointer :: set_zlbc => null ()
    procedure(zbc_fcn_ptr), pointer :: set_zhbc => null ()

    ! Initialize arrays to store boundary conditions
    real Qxhigh_ext(ny,nz,nface,nQ), Qxlow_int(ny,nz,nface,nQ)
    real Qxlow_ext(ny,nz,nface,nQ), Qxhigh_int(ny,nz,nface,nQ)
    real Qyhigh_ext(nx,nz,nface,nQ), Qylow_int(nx,nz,nface,nQ)
    real Qylow_ext(nx,nz,nface,nQ), Qyhigh_int(nx,nz,nface,nQ)
    real Qzhigh_ext(nx,ny,nface,nQ), Qzlow_int(nx,ny,nface,nQ)
    real Qzlow_ext(nx,ny,nface,nQ), Qzhigh_int(nx,ny,nface,nQ)

contains

    subroutine apply_BC
        if ( mpi_P .eq. 1 ) then
            call set_xlbc(Qxlow_int, Qxlow_ext)
        end if
        if ( mpi_P .eq. mpi_nx ) then
            call set_xhbc(Qxhigh_int, Qxhigh_ext)
        end if
        !----------------------------------------------------
        if ( mpi_Q .eq. 1 ) then
            call set_ylbc(Qylow_int, Qylow_ext)
        end if
        if ( mpi_Q .eq. mpi_ny ) then
            call set_yhbc(Qyhigh_int, Qyhigh_ext)
        end if
        !--------------------------------------------------------
        if ( mpi_R .eq. 1 ) then
            call set_zlbc(Qzlow_int, Qzlow_ext)
        end if
        if ( mpi_R .eq. mpi_nz ) then
            call set_zhbc(Qzhigh_int, Qzhigh_ext)
        end if
    end subroutine apply_BC


    !---------------------------------------------------------------------------
    ! subroutine apply_BC
    !     implicit none
    !     integer i,j,k
    !
    !     if ( mpi_P .eq. 1 ) then
    !         do k = 1,nz
    !             do j = 1,ny
    !                 call set_xlbc(j,k)
    !             end do
    !         end do
    !     end if
    !     if ( mpi_P .eq. mpi_nx ) then
    !         do k = 1,nz
    !             do j = 1,ny
    !                 call set_xhbc(j,k)
    !             end do
    !         end do
    !     end if
    !     !----------------------------------------------------
    !     if ( mpi_Q .eq. 1 ) then
    !         do k = 1,nz
    !             do i = 1,nx
    !                 call set_ylbc(i,k)
    !             end do
    !         end do
    !     end if
    !     if ( mpi_Q .eq. mpi_ny ) then
    !         do k = 1,nz
    !             do i = 1,nx
    !                 call set_yhbc(i,k)
    !             end do
    !         end do
    !     end if
    !     !--------------------------------------------------------
    !     if ( mpi_R .eq. 1 ) then
    !         do j = 1,ny
    !             do i = 1,nx
    !                 call set_zlbc(i,j)
    !             end do
    !         end do
    !     end if
    !     if ( mpi_R .eq. mpi_nz ) then
    !         do j = 1,ny
    !             do i = 1,nx
    !                 call set_zhbc(i,j)
    !             end do
    !         end do
    !     end if
    !
    ! end subroutine apply_BC
    !---------------------------------------------------------------------------


    !--------------------------------------------------------
    subroutine set_xlbc_outflow(Qxlo_e, Qxlo_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e, Qxlo_i
        ! integer j,k,i4
        do k = 1,nz
            do j = 1,ny
                do i4=1,nface
                    Qxlo_e(j,k,i4,:) = Qxlo_i(j,k,i4,:)
                end do
                if (maxval(Qxlo_e(j,k,1:nface,mx)) > 0.) then
                    Qxlo_e(j,k,1:nface,mx) = 0.
                end if
            enddo
        enddo
    end subroutine set_xlbc_outflow

    subroutine set_xhbc_outflow(Qxhi_e, Qxhi_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e, Qxhi_i
        ! integer j,k,i4
        do k = 1,nz
            do j = 1,ny
                do i4=1,nface
                    Qxhi_e(j,k,i4,:) = Qxhi_i(j,k,i4,:)
                end do
                if (minval(Qxhi_e(j,k,1:nface,mx)) < 0.) then
                    Qxhi_e(j,k,1:nface,mx) = 0.
                end if
            enddo
        enddo
    end subroutine set_xhbc_outflow
    !--------------------------------------------------------

    !--------------------------------------------------------
    subroutine set_ylbc_outflow(Qylo_e, Qylo_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e, Qylo_i
        ! integer i,k,i4
        do k = 1,nz
            do i = 1,nx
                do i4=1,nface
                    Qylo_e(i,k,i4,:) = Qylo_i(i,k,i4,:)
                end do
                if (maxval(Qylo_e(i,k,1:nface,my)) > 0.) then
                    Qylo_e(i,k,1:nface,my) = 0.
                end if
            enddo
        enddo
    end subroutine set_ylbc_outflow

    subroutine set_yhbc_outflow(Qyhi_e, Qyhi_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e, Qyhi_i
        ! integer i,k,i4
        do k = 1,nz
            do i = 1,nx
                do i4=1,nface
                    Qyhi_e(i,k,i4,:) = Qyhi_i(i,k,i4,:)
                end do
                if (minval(Qyhi_e(i,k,1:nface,my)) < 0.) then
                    Qyhi_e(i,k,1:nface,my) = 0.
                end if
            enddo
        enddo
    end subroutine set_yhbc_outflow
    !--------------------------------------------------------

    !--------------------------------------------------------
    subroutine set_zlbc_outflow(Qzlo_e, Qzlo_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e, Qzlo_i
        ! integer i,j,i4
        do j = 1,ny
            do i = 1,nx
                do i4=1,nface
                    Qzlo_e(i,j,i4,:) = Qzlo_i(i,j,i4,:)
                end do
                if (maxval(Qzlo_e(i,j,1:nface,mz)) > 0.) then
                    Qzlo_e(i,j,1:nface,mz) = 0.
                end if
            enddo
        enddo
    end subroutine set_zlbc_outflow

    subroutine set_zhbc_outflow(Qzhi_e, Qzhi_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e, Qzhi_i
        ! integer i,j,i4
        do j = 1,ny
            do i = 1,nx
                do i4=1,nface
                    Qzhi_e(i,j,i4,:) = Qzhi_i(i,j,i4,:)
                end do
                if (minval(Qzhi_e(i,j,1:nface,mz)) < 0.) then
                    Qzhi_e(i,j,1:nface,mz) = 0.
                end if
            enddo
        enddo
    end subroutine set_zhbc_outflow
    !--------------------------------------------------------


    !--------------------------------------------------------
    subroutine set_xlbc_wall(Qxlo_e, Qxlo_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e, Qxlo_i
        ! integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxlo_e(j,k,1:nface,mx) = 0.
            end do
        end do

    end subroutine set_xlbc_wall

    subroutine set_xhbc_wall(Qxhi_e, Qxhi_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e, Qxhi_i
        ! integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxhi_e(j,k,1:nface,mx) = 0.
            end do
        end do

    end subroutine set_xhbc_wall
    !--------------------------------------------------------

    !--------------------------------------------------------
    subroutine set_ylbc_wall(Qylo_e, Qylo_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e, Qylo_i
        ! integer i,k
        do k = 1,nz
            do i = 1,nx
                Qylo_e(i,k,1:nface,my) = 0.
            end do
        end do

    end subroutine set_ylbc_wall

    subroutine set_yhbc_wall(Qyhi_e, Qyhi_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e, Qyhi_i
        ! integer i,k
        do k = 1,nz
            do i = 1,nx
                Qyhi_e(i,k,1:nface,my) = 0.
            end do
        end do

    end subroutine set_yhbc_wall
    !--------------------------------------------------------

    !--------------------------------------------------------
    subroutine set_zlbc_wall(Qzlo_e, Qzlo_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e, Qzlo_i
        ! integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzlo_e(i,j,1:nface,mz) = 0.
            end do
        end do

    end subroutine set_zlbc_wall

    subroutine set_zhbc_wall(Qzhi_e, Qzhi_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e, Qzhi_i
        ! integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzhi_e(i,j,1:nface,mz) = 0
            end do
        end do

    end subroutine set_zhbc_wall
    !--------------------------------------------------------


    !--------------------------------------------------------
    subroutine set_xlbc_periodic(Qxlo_e, Qxlo_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e, Qxlo_i
    end subroutine set_xlbc_periodic

    subroutine set_xhbc_periodic(Qxhi_e, Qxhi_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e, Qxhi_i
    end subroutine set_xhbc_periodic
    !--------------------------------------------------------

    !--------------------------------------------------------
    subroutine set_ylbc_periodic(Qylo_e, Qylo_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e, Qylo_i
    end subroutine set_ylbc_periodic

    subroutine set_yhbc_periodic(Qyhi_e, Qyhi_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e, Qyhi_i
    end subroutine set_yhbc_periodic
    !--------------------------------------------------------

    !--------------------------------------------------------
    subroutine set_zlbc_periodic(Qzlo_e, Qzlo_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e, Qzlo_i
    end subroutine set_zlbc_periodic

    subroutine set_zhbc_periodic(Qzhi_e, Qzhi_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e, Qzhi_i
    end subroutine set_zhbc_periodic
    !--------------------------------------------------------


    subroutine select_x_BC(xlbc_flag, xhbc_flag, xlbc_ptr, xhbc_ptr)
        character(*), intent(in) :: xlbc_flag, xhbc_flag
        procedure(xbc_fcn_ptr), pointer :: xlbc_ptr, xhbc_ptr

        select case(xlbc_flag)
            case('outflow')
                call mpi_print(mpi_P, 'Selected outflow BCs for lower x boundary')
                xlbc_ptr => set_xlbc_outflow
            case('wall')
                call mpi_print(mpi_P, 'Selected wall BCs for lower x boundary')
                xlbc_ptr => set_xlbc_wall
            case('periodic')
                call mpi_print(mpi_P, 'Selected periodic BCs for lower x boundary')
                xlbc_ptr => set_xlbc_periodic
        end select

        select case(xhbc_flag)
            case('outflow')
                call mpi_print(mpi_P, 'Selected outflow BCs for upper x boundary')
                xhbc_ptr => set_xhbc_outflow
            case('wall')
                call mpi_print(mpi_P, 'Selected wall BCs for upper x boundary')
                xhbc_ptr => set_xhbc_wall
            case('periodic')
                call mpi_print(mpi_P, 'Selected periodic BCs for upper x boundary')
                xhbc_ptr => set_xhbc_periodic
        end select

    end subroutine select_x_BC


    subroutine select_y_BC(ylbc_flag, yhbc_flag, ylbc_ptr, yhbc_ptr)
        character(*), intent(in) :: ylbc_flag, yhbc_flag
        procedure(ybc_fcn_ptr), pointer :: ylbc_ptr, yhbc_ptr

        select case(ylbc_flag)
            case('outflow')
                call mpi_print(mpi_Q, 'Selected outflow BCs for lower y boundary')
                ylbc_ptr => set_ylbc_outflow
            case('wall')
                call mpi_print(mpi_Q, 'Selected wall BCs for lower y boundary')
                ylbc_ptr => set_ylbc_wall
            case('periodic')
                call mpi_print(mpi_Q, 'Selected periodic BCs for lower y boundary')
                ylbc_ptr => set_ylbc_periodic
        end select

        select case(yhbc_flag)
            case('outflow')
                call mpi_print(mpi_Q, 'Selected outflow BCs for upper y boundary')
                yhbc_ptr => set_yhbc_outflow
            case('wall')
                call mpi_print(mpi_Q, 'Selected wall BCs for upper y boundary')
                yhbc_ptr => set_yhbc_wall
            case('periodic')
                call mpi_print(mpi_Q, 'Selected periodic BCs for upper y boundary')
                yhbc_ptr => set_yhbc_periodic
        end select

    end subroutine select_y_BC


    subroutine select_z_BC(zlbc_flag, zhbc_flag, zlbc_ptr, zhbc_ptr)
        character(*), intent(in) :: zlbc_flag, zhbc_flag
        procedure(zbc_fcn_ptr), pointer :: zlbc_ptr, zhbc_ptr

        select case(zlbc_flag)
            case('outflow')
                call mpi_print(mpi_R, 'Selected outflow BCs for lower z boundary')
                zlbc_ptr => set_zlbc_outflow
            case('wall')
                call mpi_print(mpi_R, 'Selected wall BCs for lower z boundary')
                zlbc_ptr => set_zlbc_wall
            case('periodic')
                call mpi_print(mpi_R, 'Selected periodic BCs for lower z boundary')
                zlbc_ptr => set_zlbc_periodic
        end select

        select case(zhbc_flag)
            case('outflow')
                call mpi_print(mpi_R, 'Selected outflow BCs for upper z boundary')
                zhbc_ptr => set_zhbc_outflow
            case('wall')
                call mpi_print(mpi_R, 'Selected wall BCs for upper z boundary')
                zhbc_ptr => set_zhbc_wall
            case('periodic')
                call mpi_print(mpi_R, 'Selected periodic BCs for upper z boundary')
                zhbc_ptr => set_zhbc_periodic
        end select

    end subroutine select_z_BC



end module boundaries
