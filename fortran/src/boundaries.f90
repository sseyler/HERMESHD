module custom_boundary

    use parameters
    use helpers

    !===========================================================================
    ! Masking parameters (for advanced or internal initial/boundary conditions)
    !------------------------------------------------------------
    real, dimension(ny,nz,nface,nQ) :: Qxlow_ext_custom
    real, dimension(nx,ny,nface,nQ) :: Qcyl_ext_c, Qcyl_ext
    logical QMask(nx,ny,nz), MMask(nx,ny,nz)
    !---------------------------------------------------------------------------

contains

    function cyl_in_2d_pipe_mask(xctr, yctr, rad)
        real :: xctr, yctr, rad
        real, dimension(nx,ny,nz) :: cyl_in_2d_pipe_mask

        real :: mask(nx,ny,nz)
        real xp0, yp0, cyl2, rad2
        integer i,j

        xp0 = xc(0)        ! zero-value of x-coordinate (left)
        yp0 = yc(0)        ! zero-value of y-coordinate (bottom)
        xctr = xctr + xp0  ! center of cylinder relative to origin
        yctr = yctr + yp0  ! center of cylinder relative to origin

        mask(:,:,:) = .false.
        rad2 = rad**2

        do j = 1,ny
        do i = 1,nx
            xp = xc(i) - xctr
            yp = yc(j) - yctr
            cyl2 = xp**2 + yp**2

            if ( cyl2 <= rad2 ) then
                mask(i,j,:) = .true.
            end if
        end do
        end do

        cyl_in_2d_pipe_mask(:,:,:) = mask(:,:,:)
        return
    end function cyl_in_2d_pipe_mask


    subroutine set_cyl_in_2d_pipe_boundaries(Qxlo_e, Qcyl_e, Qmask_i, ux_amb)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qcyl_e
        logical, dimension(nx,ny,nz) :: Qmask_i
        real yp,yp0,vx,ux_amb
        integer :: j,k,ieq,i4

        yp0 = yc(0)  ! zero-value of y-coordinate (bottom)

        ! apply no-slip BCs for cylinder (set all momenta to 0)
        do ieq = mx,mz
        do i4 = 1,nface
            where ( Qmask_i(:,:,1) )
                Qcyl_e(:,:,i4,ieq) = 0.0
            end where
        end do
        end do

        ! apply parabolic inflow BCs on lower x wall
        do i4 = 1,nface
            do k = 1,nz
            do j = 1,ny
                yp = yc(j) - yp0
                vx = 4*ux_amb*yp*(ly-yp)/ly**2
                Qxlo_e(j,k,i4,mx) = dn*vx
            end do
            end do
        end do
    end subroutine set_cyl_in_2d_pipe_boundaries


    subroutine apply_xlo_in_2d_pipe_boundaries(Qxlo_i, Qxlo_e)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_i, Qxlo_e
        integer j,k,i4

        do i4 = 1,nface
            do k = 1,nz
            do j = 1,ny
                Qxlo_e(j,k,i4,mx) = Qxlow_ext_custom(j,k,i4,mx)
                Qxlo_e(j,k,i4,my:mz) = 0.0
            end do
            end do
        end do
    end subroutine apply_xlo_in_2d_pipe_boundaries


    subroutine apply_cyl_in_2d_pipe_boundaries(Qcyl_ext_c, Qcyl_e)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qcyl_ext_c, Qcyl_e
        integer j,k

        ! NOTE: This loops matches construction of Qcyl_e in set_cyl_in_2d_pipe_boundaries
        do j = 1,ny
        do i = 1,nx
            Qcyl_e(i,j,1:nface,mx:mz) = Qcyl_ext_c(i,j,1:nface,mx:mz)
        end do
        end do

        ! apply no-slip BCs for cylinder (set all momenta to 0)
        ! do ieq = mx,mz
        ! do i4 = 1,nface
        !     where ( Qmask(:,:,1) )
        !         Qcyl_e(:,:,i4,ieq) = 0.0
        !     end where
        ! end do
        ! end do
    end subroutine apply_cyl_in_2d_pipe_boundaries

end module custom_boundary




module boundaries

    use parameters
    use helpers, only : mpi_print
    use custom_boundary

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutines for setting BCs
    !------------------------------------------------------------
    abstract interface
        subroutine xbc_fcn_ptr(Qxbc_e, Qxbc_i)
            use input, only : ny,nz
            use parameters, only : mx,my,mz,nface,nQ
            real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_e, Qxbc_i
        end subroutine xbc_fcn_ptr

        subroutine ybc_fcn_ptr(Qybc_e, Qybc_i)
            use input, only : nx,nz
            use parameters, only : mx,my,mz,nface,nQ
            real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_e, Qybc_i
        end subroutine ybc_fcn_ptr

        subroutine zbc_fcn_ptr(Qzbc_e, Qzbc_i)
            use input, only : nx,ny
            use parameters, only : mx,my,mz,nface,nQ
            real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_e, Qzbc_i
        end subroutine zbc_fcn_ptr
    end interface
    !---------------------------------------------------------------------------


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


    !===========================================================================
    ! Initialize arrays to store boundary conditions
    !------------------------------------------------------------
    real Qxhigh_ext(ny,nz,nface,nQ), Qxhigh_int(ny,nz,nface,nQ)
    real Qxlow_ext(ny,nz,nface,nQ),   Qxlow_int(ny,nz,nface,nQ)
    real Qyhigh_ext(nx,nz,nface,nQ), Qyhigh_int(nx,nz,nface,nQ)
    real Qylow_ext(nx,nz,nface,nQ),   Qylow_int(nx,nz,nface,nQ)
    real Qzhigh_ext(nx,ny,nface,nQ), Qzhigh_int(nx,ny,nface,nQ)
    real Qzlow_ext(nx,ny,nface,nQ),   Qzlow_int(nx,ny,nface,nQ)
    !---------------------------------------------------------------------------

contains

    !===========================================================================
    ! Apply boundary conditions (specified by user at runtime)
    !------------------------------------------------------------
    subroutine apply_boundaries
        if ( mpi_P .eq. 1 ) then
            call apply_xlobc(Qxlow_int, Qxlow_ext)
        end if
        if ( mpi_P .eq. mpi_nx ) then
            call apply_xhibc(Qxhigh_int, Qxhigh_ext)
        end if
        !----------------------------------------------------
        if ( mpi_Q .eq. 1 ) then
            call apply_ylobc(Qylow_int, Qylow_ext)
        end if
        if ( mpi_Q .eq. mpi_ny ) then
            call apply_yhibc(Qyhigh_int, Qyhigh_ext)
        end if
        !--------------------------------------------------------
        if ( mpi_R .eq. 1 ) then
            call apply_zlobc(Qzlow_int, Qzlow_ext)
        end if
        if ( mpi_R .eq. mpi_nz ) then
            call apply_zhibc(Qzhigh_int, Qzhigh_ext)
        end if

        ! NOTE: Inefficient b/c overwrites 'default' BC application above
        ! NOTE: This is NOT a general implementation for applying custom BCs
        call apply_cyl_in_2d_pipe_boundaries(Qcyl_ext_c, Qcyl_ext)
        if ( mpi_P .eq. 1 ) then
            call apply_xlo_in_2d_pipe_boundaries(Qxlow_int, Qxlow_ext)
        end if
    end subroutine apply_boundaries
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Select user-specified BC (subroutines) at runtime
    !------------------------------------------------------------
    subroutine select_x_boundaries(xlobc_flag, xhibc_flag, xlobc_ptr, xhibc_ptr)
        character(*), intent(in) :: xlobc_flag, xhibc_flag
        procedure(xbc_fcn_ptr), pointer :: xlobc_ptr, xhibc_ptr

        select case(xlobc_flag)
            case('periodic')
                call mpi_print(mpi_P, 'Selected periodic BCs for lower x boundary')
                xlobc_ptr => xlobc_periodic
            case('outflow')
                call mpi_print(mpi_P, 'Selected outflow BCs for lower x boundary')
                xlobc_ptr => xlobc_outflow
            case('wall')
                call mpi_print(mpi_P, 'Selected wall BCs for lower x boundary')
                xlobc_ptr => xlobc_wall
            case('noslip')
                call mpi_print(mpi_P, 'Selected no-slip BCs for lower x boundary')
                xlobc_ptr => xlobc_noslip
        end select
        select case(xhibc_flag)
            case('periodic')
                call mpi_print(mpi_P, 'Selected periodic BCs for upper x boundary')
                xhibc_ptr => xhibc_periodic
            case('outflow')
                call mpi_print(mpi_P, 'Selected outflow BCs for upper x boundary')
                xhibc_ptr => xhibc_outflow
            case('wall')
                call mpi_print(mpi_P, 'Selected wall BCs for upper x boundary')
                xhibc_ptr => xhibc_wall
            case('noslip')
                call mpi_print(mpi_P, 'Selected no-slip BCs for upper x boundary')
                xhibc_ptr => xhibc_noslip
        end select

    end subroutine select_x_boundaries
    !--------------------------------------------------------
    subroutine select_y_boundaries(ylobc_flag, yhibc_flag, ylobc_ptr, yhibc_ptr)
        character(*), intent(in) :: ylobc_flag, yhibc_flag
        procedure(ybc_fcn_ptr), pointer :: ylobc_ptr, yhibc_ptr

        select case(ylobc_flag)
            case('periodic')
                call mpi_print(mpi_Q, 'Selected periodic BCs for lower y boundary')
                ylobc_ptr => ylobc_periodic
            case('outflow')
                call mpi_print(mpi_Q, 'Selected outflow BCs for lower y boundary')
                ylobc_ptr => ylobc_outflow
            case('wall')
                call mpi_print(mpi_Q, 'Selected wall BCs for lower y boundary')
                ylobc_ptr => ylobc_wall
            case('noslip')
                call mpi_print(mpi_Q, 'Selected no-slip BCs for lower y boundary')
                ylobc_ptr => ylobc_noslip
        end select
        select case(yhibc_flag)
            case('periodic')
                call mpi_print(mpi_Q, 'Selected periodic BCs for upper y boundary')
                yhibc_ptr => yhibc_periodic
            case('outflow')
                call mpi_print(mpi_Q, 'Selected outflow BCs for upper y boundary')
                yhibc_ptr => yhibc_outflow
            case('wall')
                call mpi_print(mpi_Q, 'Selected wall BCs for upper y boundary')
                yhibc_ptr => yhibc_wall
            case('noslip')
                call mpi_print(mpi_Q, 'Selected no-slip BCs for upper y boundary')
                yhibc_ptr => yhibc_noslip
        end select

    end subroutine select_y_boundaries
    !--------------------------------------------------------
    subroutine select_z_boundaries(zlobc_flag, zhibc_flag, zlobc_ptr, zhibc_ptr)
        character(*), intent(in) :: zlobc_flag, zhibc_flag
        procedure(zbc_fcn_ptr), pointer :: zlobc_ptr, zhibc_ptr

        select case(zlobc_flag)
            case('periodic')
                call mpi_print(mpi_R, 'Selected periodic BCs for lower z boundary')
                zlobc_ptr => zlobc_periodic
            case('outflow')
                call mpi_print(mpi_R, 'Selected outflow BCs for lower z boundary')
                zlobc_ptr => zlobc_outflow
            case('wall')
                call mpi_print(mpi_R, 'Selected wall BCs for lower z boundary')
                zlobc_ptr => zlobc_wall
            case('noslip')
                call mpi_print(mpi_R, 'Selected no-slip BCs for lower z boundary')
                zlobc_ptr => zlobc_noslip
        end select
        select case(zhibc_flag)
            case('periodic')
                call mpi_print(mpi_R, 'Selected periodic BCs for upper z boundary')
                zhibc_ptr => zhibc_periodic
            case('outflow')
                call mpi_print(mpi_R, 'Selected outflow BCs for upper z boundary')
                zhibc_ptr => zhibc_outflow
            case('wall')
                call mpi_print(mpi_R, 'Selected wall BCs for upper z boundary')
                zhibc_ptr => zhibc_wall
            case('noslip')
                call mpi_print(mpi_R, 'Selected no-slip BCs for upper z boundary')
                zhibc_ptr => zhibc_noslip
        end select

    end subroutine select_z_boundaries
    !---------------------------------------------------------------------------


    !===========================================================================
    ! PERIODIC BC subroutines
    !------------------------------------------------------------
    subroutine xlobc_periodic(Qxlo_e, Qxlo_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e, Qxlo_i
    end subroutine xlobc_periodic

    subroutine xhibc_periodic(Qxhi_e, Qxhi_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e, Qxhi_i
    end subroutine xhibc_periodic
    !--------------------------------------------------------
    subroutine ylobc_periodic(Qylo_e, Qylo_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e, Qylo_i
    end subroutine ylobc_periodic

    subroutine yhibc_periodic(Qyhi_e, Qyhi_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e, Qyhi_i
    end subroutine yhibc_periodic
    !--------------------------------------------------------
    subroutine zlobc_periodic(Qzlo_e, Qzlo_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e, Qzlo_i
    end subroutine zlobc_periodic

    subroutine zhibc_periodic(Qzhi_e, Qzhi_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e, Qzhi_i
    end subroutine zhibc_periodic
    !--------------------------------------------------------


    !===========================================================================
    ! OUTFLOW BC subroutine
    !------------------------------------------------------------
    subroutine xlobc_outflow(Qxlo_e, Qxlo_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e, Qxlo_i
        integer j,k
        do k = 1,nz
            do j = 1,ny
                    Qxlo_e(j,k,1:nface,:) = Qxlo_i(j,k,1:nface,:)
                if (maxval(Qxlo_e(j,k,1:nface,mx)) > 0.) then
                    Qxlo_e(j,k,1:nface,mx) = 0.0
                end if
            enddo
        enddo
    end subroutine xlobc_outflow
    !--------------------------------------------------------
    subroutine xhibc_outflow(Qxhi_e, Qxhi_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e, Qxhi_i
        integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxhi_e(j,k,1:nface,:) = Qxhi_i(j,k,1:nface,:)
                if (minval(Qxhi_e(j,k,1:nface,mx)) < 0.) then
                    Qxhi_e(j,k,1:nface,mx) = 0.0
                end if
            enddo
        enddo
    end subroutine xhibc_outflow
    !--------------------------------------------------------
    subroutine ylobc_outflow(Qylo_e, Qylo_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e, Qylo_i
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qylo_e(i,k,1:nface,:) = Qylo_i(i,k,1:nface,:)
                if (maxval(Qylo_e(i,k,1:nface,my)) > 0.) then
                    Qylo_e(i,k,1:nface,my) = 0.0
                end if
            enddo
        enddo
    end subroutine ylobc_outflow
    !--------------------------------------------------------
    subroutine yhibc_outflow(Qyhi_e, Qyhi_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e, Qyhi_i
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qyhi_e(i,k,1:nface,:) = Qyhi_i(i,k,1:nface,:)
                if (minval(Qyhi_e(i,k,1:nface,my)) < 0.) then
                    Qyhi_e(i,k,1:nface,my) = 0.0
                end if
            enddo
        enddo
    end subroutine yhibc_outflow
    !--------------------------------------------------------
    subroutine zlobc_outflow(Qzlo_e, Qzlo_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e, Qzlo_i
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzlo_e(i,j,1:nface,:) = Qzlo_i(i,j,1:nface,:)
                if (maxval(Qzlo_e(i,j,1:nface,mz)) > 0.) then
                    Qzlo_e(i,j,1:nface,mz) = 0.0
                end if
            enddo
        enddo
    end subroutine zlobc_outflow
    !--------------------------------------------------------
    subroutine zhibc_outflow(Qzhi_e, Qzhi_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e, Qzhi_i
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzhi_e(i,j,1:nface,:) = Qzhi_i(i,j,1:nface,:)
                if (minval(Qzhi_e(i,j,1:nface,mz)) < 0.) then
                    Qzhi_e(i,j,1:nface,mz) = 0.0
                end if
            enddo
        enddo
    end subroutine zhibc_outflow
    !---------------------------------------------------------------------------


    !===========================================================================
    ! WALL BC subroutine
    !------------------------------------------------------------
    subroutine xlobc_wall(Qxlo_e, Qxlo_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e, Qxlo_i
        integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxlo_e(j,k,1:nface,mx) = 0.0
            end do
        end do

    end subroutine xlobc_wall
    !--------------------------------------------------------
    subroutine xhibc_wall(Qxhi_e, Qxhi_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e, Qxhi_i
        integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxhi_e(j,k,1:nface,mx) = 0.0
            end do
        end do

    end subroutine xhibc_wall
    !--------------------------------------------------------
    subroutine ylobc_wall(Qylo_e, Qylo_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e, Qylo_i
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qylo_e(i,k,1:nface,my) = 0.0
            end do
        end do

    end subroutine ylobc_wall

    subroutine yhibc_wall(Qyhi_e, Qyhi_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e, Qyhi_i
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qyhi_e(i,k,1:nface,my) = 0.0
            end do
        end do

    end subroutine yhibc_wall
    !--------------------------------------------------------
    subroutine zlobc_wall(Qzlo_e, Qzlo_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e, Qzlo_i
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzlo_e(i,j,1:nface,mz) = 0.0
            end do
        end do

    end subroutine zlobc_wall
    !--------------------------------------------------------
    subroutine zhibc_wall(Qzhi_e, Qzhi_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e, Qzhi_i
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzhi_e(i,j,1:nface,mz) = 0.0
            end do
        end do

    end subroutine zhibc_wall
    !---------------------------------------------------------------------------


    !===========================================================================
    ! NO-SLIP BC subroutine
    !------------------------------------------------------------
    subroutine xlobc_noslip(Qxlo_e, Qxlo_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e, Qxlo_i
        integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxlo_e(j,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine xlobc_noslip
    !--------------------------------------------------------
    subroutine xhibc_noslip(Qxhi_e, Qxhi_i)
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e, Qxhi_i
        integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxhi_e(j,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine xhibc_noslip
    !--------------------------------------------------------
    subroutine ylobc_noslip(Qylo_e, Qylo_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e, Qylo_i
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qylo_e(i,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine ylobc_noslip

    subroutine yhibc_noslip(Qyhi_e, Qyhi_i)
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e, Qyhi_i
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qyhi_e(i,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine yhibc_noslip
    !--------------------------------------------------------
    subroutine zlobc_noslip(Qzlo_e, Qzlo_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e, Qzlo_i
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzlo_e(i,j,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine zlobc_noslip
    !--------------------------------------------------------
    subroutine zhibc_noslip(Qzhi_e, Qzhi_i)
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e, Qzhi_i
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzhi_e(i,j,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine zhibc_noslip
    !---------------------------------------------------------------------------


end module boundaries
