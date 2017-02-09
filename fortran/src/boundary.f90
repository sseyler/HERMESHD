module boundary_custom

    use parameters
    use helpers

    !===========================================================================
    ! Masking parameters (for advanced or internal initial/boundary conditions)
    !------------------------------------------------------------
    real, dimension(ny,nz,nface,nQ) :: Qxlow_ext_c
    real, dimension(nx,ny,nface,nQ) :: Qcyl_ext_c, Qcyl_ext
    logical QMask(nx,ny,nz), MMask(nx,ny,nz)
    !---------------------------------------------------------------------------

contains

    function cyl_in_2d_pipe_mask(xctr, yctr, rad)
        real, intent(in) :: xctr, yctr, rad
        real, dimension(nx,ny,nz) :: cyl_in_2d_pipe_mask
        real :: mask(nx,ny,nz), cyl2, rad2
        integer i,j

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


    subroutine set_cyl_in_2d_pipe_boundaries(Qmask_i, ux_amb, den, Qxlo_e_c, Qcyl_e_c)
        logical, dimension(nx,ny,nz), intent(in) :: Qmask_i
        real, intent(in) :: ux_amb, den
        real, dimension(ny,nz,nface,nQ), intent(out) :: Qxlo_e_c
        real, dimension(nx,ny,nface,nQ), intent(out) :: Qcyl_e_c

        real yp,yp0,vx
        integer :: j,k,ieq,i4

        yp0 = lyd  ! set zero-value of y-coordinate to domain bottom

        ! apply no-slip BCs for cylinder (set all momenta to 0)
        do ieq = rh,en
        do i4 = 1,nface
            Qcyl_e_c(:,:,i4,ieq) = Q_r0(:,:,1,i4,ieq)  ! Initialize to grid values
            where ( Qmask_i(:,:,1) )
                Qcyl_e_c(:,:,i4,ieq) = 0.0             ! Set density/momenta to zero
            end where
        end do
        end do

        ! apply parabolic inflow BCs on lower x wall
        do i4 = 1,nface
            do k = 1,nz
            do j = 1,ny
                yp = yc(j) - yp0
                vx = 4*ux_amb*yp*(ly-yp)/ly**2
                Qxlo_e_c(j,k,i4,mx) = den*vx
            end do
            end do
        end do
    end subroutine set_cyl_in_2d_pipe_boundaries


    subroutine apply_xlo_in_2d_pipe_boundaries(Qxlo_e_c, Qxlo_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxlo_e_c
        real, dimension(ny,nz,nface,nQ), intent(out) :: Qxlo_e
        integer j,k,i4

        do i4 = 1,nface
            do k = 1,nz
            do j = 1,ny
                Qxlo_e(j,k,i4,mx) = Qxlo_e_c(j,k,i4,mx)
                Qxlo_e(j,k,i4,my:mz) = 0.0
            end do
            end do
        end do
    end subroutine apply_xlo_in_2d_pipe_boundaries


    subroutine apply_cyl_in_2d_pipe_boundaries(Qcyl_e_c, Qcyl_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qcyl_e_c
        real, dimension(nx,ny,nface,nQ), intent(out) :: Qcyl_e
        integer j,k

        ! NOTE: This loops matches construction of Qcyl_e in set_cyl_in_2d_pipe_boundaries
        ! do j = 1,ny
        ! do i = 1,nx
        !     Qcyl_e(i,j,1:nface,rh:mz) = Qcyl_e_c(i,j,1:nface,rh:mz)
        ! end do
        ! end do

        ! apply no-slip BCs for cylinder (set all momenta to 0)
        ! do ieq = rh,mz
        ! do i4 = 1,nface
        !     where ( Qmask(:,:,1) )
        !         Qcyl_e(:,:,i4,ieq) = 0.0
        !     end where
        ! end do
        ! end do

        do i4=1,nface
            do ieq=rh,en
                where (Qmask(:,:,:))
                    flux_x(i4,1:nx,:,:,ieq) = 0.0
                    flux_y(i4,:,1:ny,:,ieq) = 0.0
                    flux_z(i4,:,:,1:nz,ieq) = 0.0
                end where
            end do
        end do
        where (Qmask(:,:,:))
            Q_r0(:,:,:,rh,1) = 2.0
            Q_r0(:,:,:,mx,1) = 0.0
            Q_r0(:,:,:,my,1) = 0.0
            Q_r0(:,:,:,mz,1) = 0.0
        end where
    end subroutine apply_cyl_in_2d_pipe_boundaries

end module boundary_custom


module boundary

    use parameters
    use helpers, only : mpi_print
    use boundary_custom

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutines for setting BCs
    !------------------------------------------------------------
    abstract interface
        ! subroutine bc_fcn_ptr(nc1, nc2, Qbc_i, Qbc_e)
        !     use parameters, only : mx,my,mz,nface,nQ
        !     integer, intent(in) :: nc1, nc2
        !     real, dimension(nc1,nc2,nface,nQ), intent(in) :: Qxbc_i
        !     real, dimension(nc1,nc2,nface,nQ), intent(inout) :: Qxbc_e
        ! end subroutine bc_fcn_ptr

        subroutine xbc_fcn_ptr(Qxbc_i, Qxbc_e)
            use input, only : ny,nz
            use parameters, only : mx,my,mz,nface,nQ
            real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_i
            real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_e
        end subroutine xbc_fcn_ptr

        subroutine ybc_fcn_ptr(Qybc_i, Qybc_e)
            use input, only : nx,nz
            use parameters, only : mx,my,mz,nface,nQ
            real, dimension(nx,nz,nface,nQ), intent(in) :: Qybc_i
            real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_e
        end subroutine ybc_fcn_ptr

        subroutine zbc_fcn_ptr(Qzbc_i, Qzbc_e)
            use input, only : nx,ny
            use parameters, only : mx,my,mz,nface,nQ
            real, dimension(nx,ny,nface,nQ), intent(in) :: Qzbc_i
            real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_e
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
        ! call apply_cyl_in_2d_pipe_boundaries(Qcyl_ext_c, Qcyl_ext)
        ! print *,QMask(:,:,1)
        ! if ( mpi_P .eq. 1 ) then
        !     call apply_xlo_in_2d_pipe_boundaries(Qxlow_ext_c, Qxlow_ext)
            ! print *,Qxlow_ext_c(:,:,1,rh:mx)
            ! print *,''
        ! end if
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
    subroutine select_y_boundaries(ylobc_flag, yhibc_flag, ylobc_ptr, yhibc_ptr)
        character(*), intent(in) :: ylobc_flag, yhibc_flag
        procedure(ybc_fcn_ptr), pointer :: ylobc_ptr, yhibc_ptr

        select case (ylobc_flag)
            case ('periodic')
                ylobc_ptr => ylobc_periodic
            case ('outflow')
                ylobc_ptr => ylobc_outflow
            case ('wall')
                ylobc_ptr => ylobc_wall
            case ('noslip')
                ylobc_ptr => ylobc_noslip
            case default
                ylobc_ptr => ylobc_periodic
        end select
        select case (yhibc_flag)
            case ('periodic')
                yhibc_ptr => yhibc_periodic
            case ('outflow')
                yhibc_ptr => yhibc_outflow
            case ('wall')
                yhibc_ptr => yhibc_wall
            case ('noslip')
                yhibc_ptr => yhibc_noslip
            case default
                yhibc_ptr => yhibc_periodic
        end select

    end subroutine select_y_boundaries
    !--------------------------------------------------------
    subroutine select_z_boundaries(zlobc_flag, zhibc_flag, zlobc_ptr, zhibc_ptr)
        character(*), intent(in) :: zlobc_flag, zhibc_flag
        procedure(zbc_fcn_ptr), pointer :: zlobc_ptr, zhibc_ptr

        select case (zlobc_flag)
            case ('periodic')
                zlobc_ptr => zlobc_periodic
            case ('outflow')
                zlobc_ptr => zlobc_outflow
            case ('wall')
                zlobc_ptr => zlobc_wall
            case ('noslip')
                zlobc_ptr => zlobc_noslip
            case default
                zlobc_ptr => zlobc_periodic
        end select
        select case (zhibc_flag)
            case ('periodic')
                zhibc_ptr => zhibc_periodic
            case ('outflow')
                zhibc_ptr => zhibc_outflow
            case ('wall')
                zhibc_ptr => zhibc_wall
            case ('noslip')
                zhibc_ptr => zhibc_noslip
            case default
                zhibc_ptr => zhibc_periodic
        end select

    end subroutine select_z_boundaries
    !---------------------------------------------------------------------------


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
    ! OUTFLOW BC subroutine
    !------------------------------------------------------------
    subroutine xlobc_outflow(Qxlo_i, Qxlo_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxlo_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e
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
    subroutine xhibc_outflow(Qxhi_i, Qxhi_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxhi_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e
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
    subroutine ylobc_outflow(Qylo_i, Qylo_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qylo_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e
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
    subroutine yhibc_outflow(Qyhi_i, Qyhi_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qyhi_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e
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
    subroutine zlobc_outflow(Qzlo_i, Qzlo_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzlo_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e
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
    subroutine zhibc_outflow(Qzhi_i, Qzhi_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzhi_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e
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


    !===========================================================================
    ! NO-SLIP BC subroutine
    !------------------------------------------------------------
    subroutine xlobc_noslip(Qxlo_i, Qxlo_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxlo_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e
        integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxlo_e(j,k,1:nface,:) = Qxlo_i(j,k,1:nface,:)
                Qxlo_e(j,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine xlobc_noslip
    !--------------------------------------------------------
    subroutine xhibc_noslip(Qxhi_i, Qxhi_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxhi_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxhi_e
        integer j,k
        do k = 1,nz
            do j = 1,ny
                Qxhi_e(j,k,1:nface,:) = Qxhi_i(j,k,1:nface,:)
                Qxhi_e(j,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine xhibc_noslip
    !--------------------------------------------------------
    subroutine ylobc_noslip(Qylo_i, Qylo_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qylo_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qylo_e(i,k,1:nface,:) = Qylo_i(i,k,1:nface,:)
                Qylo_e(i,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine ylobc_noslip

    subroutine yhibc_noslip(Qyhi_i, Qyhi_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qyhi_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qyhi_e
        integer i,k
        do k = 1,nz
            do i = 1,nx
                Qyhi_e(i,k,1:nface,:) = Qyhi_i(i,k,1:nface,:)
                Qyhi_e(i,k,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine yhibc_noslip
    !--------------------------------------------------------
    subroutine zlobc_noslip(Qzlo_i, Qzlo_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzlo_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzlo_e
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzlo_e(i,j,1:nface,:) = Qzlo_i(i,j,1:nface,:)
                Qzlo_e(i,j,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine zlobc_noslip
    !--------------------------------------------------------
    subroutine zhibc_noslip(Qzhi_i, Qzhi_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qzhi_i
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzhi_e
        integer i,j
        do j = 1,ny
            do i = 1,nx
                Qzhi_e(i,j,1:nface,:) = Qzhi_i(i,j,1:nface,:)
                Qzhi_e(i,j,1:nface,mx:mz) = 0.0
            end do
        end do

    end subroutine zhibc_noslip
    !---------------------------------------------------------------------------


end module boundary
