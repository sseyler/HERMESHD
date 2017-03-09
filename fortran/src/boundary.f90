!***** BOUNDARY.F90 **********************************************************************

!##### boundary_defs ###########################################################
module boundary_defs

use parameters

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutines for setting BCs
    !------------------------------------------------------------
    abstract interface
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
    ! OUTFLOW BC subroutine
    !------------------------------------------------------------
    subroutine xlobc_outflow(Qxlo_i, Qxlo_e)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxlo_i
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxlo_e
        integer j,k
        do k = 1,nz
        do j = 1,ny
            Qxlo_e(j,k,1:nface,:) = Qxlo_i(j,k,1:nface,:)
            if (maxval(Qxlo_e(j,k,1:nface,mx)) > 0.0) then
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
        real, dimension(nface) :: dn_i,dni_i, mx_i,my_i,mz_i,vx_i, Ener_i,Ener_e,P_i
        real mflowrate_out, mflowrate_net

        ! Qxlo_e(i,k,1:nface,:) = Qxlo_i(i,k,1:nface,:)
        ! if (minval(Qxlo_e(i,k,1:nface,my)) < 0.0) then
        !     Qxlo_e(i,k,1:nface,my) = 0.0
        ! end if

        do k = 1,nz
        do j = 1,ny
            where(Qxlow_ext(j,k,:,mx) > 0.0) Qxlow_ext(j,k,:,mx) = -Qxlow_int(j,k,:,mx)
        enddo
        enddo

        ! mflowrate_out = 0.0
        ! do k = 1,nz
        ! do j = 1,ny
        !     Qxhi_e(j,k,1:nface,:) = Qxhi_i(j,k,1:nface,:)
        !
        !     mx_i = Qxhi_i(j,k,1:nface,mx)
        !     mflowrate_out = mflowrate_out + sum(mx_i)
        ! enddo
        ! enddo
        ! ! print *,'mflowrate_out', mflowrate_out,'mflowrate', mflowrate
        ! mflowrate_net = mflowrate_out - mflowrate
        ! do k = 1,nz
        ! do j = 1,ny
        !     Qxhi_e(j,k,1:nface,mx) = Qxhi_e(j,k,1:nface,mx) - mflowrate_net/ny/nface
        ! enddo
        ! enddo

        ! print *,Qxhi_e(:,1,1,en)
    end subroutine xhibc_outflow
    !--------------------------------------------------------
    subroutine ylobc_outflow(Qylo_i, Qylo_e)
        real, dimension(nx,nz,nface,nQ), intent(in) :: Qylo_i
        real, dimension(nx,nz,nface,nQ), intent(inout) :: Qylo_e
        integer i,k
        do k = 1,nz
        do i = 1,nx
            Qylo_e(i,k,1:nface,:) = Qylo_i(i,k,1:nface,:)
            if (maxval(Qylo_e(i,k,1:nface,my)) > 0.0) then
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
            if (minval(Qyhi_e(i,k,1:nface,my)) < 0.0) then
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
            if (maxval(Qzlo_e(i,j,1:nface,mz)) > 0.0) then
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
            if (minval(Qzhi_e(i,j,1:nface,mz)) < 0.0) then
                Qzhi_e(i,j,1:nface,mz) = 0.0
            end if
        enddo
        enddo
    end subroutine zhibc_outflow
    !---------------------------------------------------------------------------

end module boundary_defs


!##### boundary_defs2 ##########################################################
! NOTE: This module is currently for testing purposes only. The boundary_defs
!       module should be used to ensure correctness.
! The purpose of this module is to reduce the amount of code needed to specify
! boundary conditions by re-using procedures when setting boundary conditions
! for a given direction. For instance, if we specify no-slip boundary conditions
! for both the high and low y-boundary, then the subroutine to set the no-slip
! BC should be the same for the high and low y-boundary. Ideally, we would think
! only a single subroutine is needed for a particular BC in a given direction,
! so only THREE subroutines would be needed in total for each type of BC for the
! x, y, and z boundaries, rather than SIX (two for high/low in each direction).
!
! However, this doesn't seems to work in practice since the code crashes (after
! a short number of time steps): since the BC subroutines are selected at
! runtime using procedure pointers, there is probably some kind of interference
! between the BC subroutine calls since they are really pointers to the
! underlying selected BC procedure, and/OR the problem could be that MPI doesn't
! play nicely with this approach.
!
! For the time being, there must be a separate BC subroutine for a given BC type
! for each direction (i.e., for, say, no-slip, SIX subroutines are needed to)
! cover high/low BCs for each direction [x,y,z], rather than THREE).
!
! module boundary_defs2
!
! use parameters
!
!     !===========================================================================
!     ! ABSTRACT INTERFACE to subroutines for setting BCs
!     !------------------------------------------------------------
!     abstract interface
!         subroutine xbc_fcn_ptr(Qxbc_int, Qxbc_ext)
!             use input, only : ny,nz
!             use parameters
!             real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_int
!             real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_ext
!         end subroutine xbc_fcn_ptr
!
!         subroutine ybc_fcn_ptr(Qybc_int, Qybc_ext)
!             use input, only : nx,nz
!             use parameters
!             real, dimension(nx,nz,nface,nQ), intent(in) :: Qybc_int
!             real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_ext
!         end subroutine ybc_fcn_ptr
!
!         subroutine zbc_fcn_ptr(Qzbc_int, Qzbc_ext)
!             use input, only : nx,ny
!             use parameters
!             real, dimension(nx,ny,nface,nQ), intent(in) :: Qzbc_int
!             real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_ext
!         end subroutine zbc_fcn_ptr
!     end interface
!     !---------------------------------------------------------------------------
!
! contains
!
!     !===========================================================================
!     ! PERIODIC BC subroutines
!     !------------------------------------------------------------
!     subroutine x_bc_periodic(Qxbc_int, Qxbc_ext)
!         real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_int
!         real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_ext
!     end subroutine x_bc_periodic
!     !--------------------------------------------------------
!     subroutine y_bc_periodic(Qybc_int, Qybc_ext)
!         real, dimension(nx,nz,nface,nQ), intent(in) :: Qybc_int
!         real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_ext
!     end subroutine y_bc_periodic
!     !--------------------------------------------------------
!     subroutine z_bc_periodic(Qzbc_int, Qzbc_ext)
!         real, dimension(nx,ny,nface,nQ), intent(in) :: Qzbc_int
!         real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_ext
!     end subroutine z_bc_periodic
!     !---------------------------------------------------------------------------
!
!
!     !===========================================================================
!     ! NO-SLIP BC subroutine
!     !------------------------------------------------------------
!     subroutine x_bc_noslip(Qxbc_int, Qxbc_ext)
!         real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_int
!         real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_ext
!         integer j,k
!         do k = 1,nz
!         do j = 1,ny
!             Qxbc_ext(j,k,1:nface,:) = Qxbc_int(j,k,1:nface,:)
!             Qxbc_ext(j,k,1:nface,mx:mz) = 0.0
!         end do
!         end do
!     end subroutine x_bc_noslip
!     !--------------------------------------------------------
!     subroutine y_bc_noslip(Qybc_int, Qybc_ext)
!         real, dimension(nx,nz,nface,nQ), intent(in) :: Qybc_int
!         real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_ext
!         integer i,k
!         do k = 1,nz
!         do i = 1,nx
!             Qybc_ext(i,k,1:nface,:) = Qybc_int(i,k,1:nface,:)
!             Qybc_ext(i,k,1:nface,mx:mz) = 0.0
!         end do
!         end do
!     end subroutine y_bc_noslip
!     !--------------------------------------------------------
!     subroutine z_bc_noslip(Qzbc_int, Qzbc_ext)
!         real, dimension(nx,ny,nface,nQ), intent(in) :: Qzbc_int
!         real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_ext
!         integer i,j
!         do j = 1,ny
!         do i = 1,nx
!             Qzbc_ext(i,j,1:nface,:) = Qzbc_int(i,j,1:nface,:)
!             Qzbc_ext(i,j,1:nface,mx:mz) = 0.0
!         end do
!         end do
!     end subroutine z_bc_noslip
!     !---------------------------------------------------------------------------
!
!
!     !===========================================================================
!     ! WALL BC subroutine
!     !------------------------------------------------------------
!     subroutine x_bc_wall(Qxbc_int, Qxbc_ext)
!         real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_int
!         real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_ext
!         integer j,k
!         do k = 1,nz
!         do j = 1,ny
!             Qxbc_ext(j,k,1:nface,:) = Qxbc_int(j,k,1:nface,:)
!             Qxbc_ext(j,k,1:nface,mx) = 0.0
!         end do
!         end do
!     end subroutine x_bc_wall
!     !--------------------------------------------------------
!     subroutine y_bc_wall(Qybc_int, Qybc_ext)
!         real, dimension(nx,nz,nface,nQ), intent(in) :: Qybc_int
!         real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_ext
!         integer i,k
!         do k = 1,nz
!         do i = 1,nx
!             Qybc_ext(i,k,1:nface,:) = Qybc_int(i,k,1:nface,:)
!             Qybc_ext(i,k,1:nface,my) = 0.0
!         end do
!         end do
!     end subroutine y_bc_wall
!     !--------------------------------------------------------
!     subroutine z_bc_wall(Qzbc_int, Qzbc_ext)
!         real, dimension(nx,ny,nface,nQ), intent(in) :: Qzbc_int
!         real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_ext
!         integer i,j
!         do j = 1,ny
!         do i = 1,nx
!             Qzbc_ext(i,j,1:nface,:) = Qzbc_int(i,j,1:nface,:)
!             Qzbc_ext(i,j,1:nface,mz) = 0.0
!         end do
!         end do
!     end subroutine z_bc_wall
!     !---------------------------------------------------------------------------
!
!
!     !===========================================================================
!     ! OUTFLOW BC subroutine
!     !------------------------------------------------------------
!     subroutine x_bc_outflow(Qxbc_int, Qxbc_ext)
!         real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_int
!         real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_ext
!         integer j,k
!         do k = 1,nz
!         do j = 1,ny
!             Qxbc_ext(j,k,1:nface,:) = Qxbc_int(j,k,1:nface,:)
!             if (maxval(Qxbc_ext(j,k,1:nface,mx)) > 0.0) then
!                 Qxbc_ext(j,k,1:nface,mx) = 0.0
!             end if
!         enddo
!         enddo
!     end subroutine x_bc_outflow
!     !--------------------------------------------------------
!     subroutine y_bc_outflow(Qybc_int, Qybc_ext)
!         real, dimension(nx,nz,nface,nQ), intent(in) :: Qybc_int
!         real, dimension(nx,nz,nface,nQ), intent(inout) :: Qybc_ext
!         integer i,k
!         do k = 1,nz
!         do i = 1,nx
!             Qybc_ext(i,k,1:nface,:) = Qybc_int(i,k,1:nface,:)
!             if (maxval(Qybc_ext(i,k,1:nface,my)) > 0.0) then
!                 Qybc_ext(i,k,1:nface,my) = 0.0
!             end if
!         enddo
!         enddo
!     end subroutine y_bc_outflow
!     !--------------------------------------------------------
!     subroutine z_bc_outflow(Qzbc_int, Qzbc_ext)
!         real, dimension(nx,ny,nface,nQ), intent(in) :: Qzbc_int
!         real, dimension(nx,ny,nface,nQ), intent(inout) :: Qzbc_ext
!         integer i,j
!         do j = 1,ny
!         do i = 1,nx
!             Qzbc_ext(i,j,1:nface,:) = Qzbc_int(i,j,1:nface,:)
!             if (maxval(Qzbc_ext(i,j,1:nface,mz)) > 0.0) then
!                 Qzbc_ext(i,j,1:nface,mz) = 0.0
!             end if
!         enddo
!         enddo
!     end subroutine z_bc_outflow
!     !---------------------------------------------------------------------------
!
! end module boundary_defs2



!##### boundary_custom #########################################################
module boundary_custom

    use parameters
    use helpers

    !===========================================================================
    ! Masking parameters (for advanced or internal initial/boundary conditions)
    !------------------------------------------------------------
    real, dimension(ny,nz,nface,nQ) :: Qxlo_ext_def, Qxlo_ext_scale
    real, dimension(nx,ny,nface,nQ) :: Qcyl_ext_c, Qcyl_ext
    logical QMask(nx,ny,nz), MMask(nx,ny,nz)
    !---------------------------------------------------------------------------

contains

    subroutine x_bc_inflow(Qxbc_int, Qxbc_def, Qxbc_ext)
        real, dimension(ny,nz,nface,nQ), intent(in) :: Qxbc_def
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_int
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_ext
        integer j,k,i4

        do i4=1,nface
            do k = 1,nz
            do j = 1,ny
                ! Initially, set all values to internal face values
                Qxbc_ext(j,k,i4,:)     = Qxbc_int(j,k,i4,:)

                ! Set external face values to those specified by Qxbc_def
                Qxbc_ext(j,k,i4,rh)    = Qxbc_def(j,k,i4,rh)

                Qxbc_ext(j,k,i4,mx)    = Qxbc_def(j,k,i4,mx)
                Qxbc_ext(j,k,i4,my:mz) = 0.0
                Qxbc_ext(j,k,i4,en)    = Qxbc_def(j,k,i4,en)
            end do
            end do
        end do
    end subroutine x_bc_inflow


    subroutine set_xlobc_inflow(ux_amb, den, pres, Qxbc_d)
        real, intent(in) :: ux_amb, den, pres
        real, dimension(ny,nz,nface,nQ), intent(inout) :: Qxbc_d

        real yp,yp0,vx
        integer j,k,i4

        yp0 = lyd  ! set zero-value of y-coordinate to domain bottom

        mflowrate = 0.0
        ! apply parabolic inflow BCs on lower x wall
        ! Qxbc_d(:,:,:,:) = 0.0
        do i4 = 1,nface
            do k = 1,nz
            do j = 1,ny
                yp = yc(j) - yp0
                vx = 4*ux_amb*yp*(ly-yp)/ly**2
                mflowrate = mflowrate + den*vx

                Qxbc_d(j,k,i4,rh) = den
                Qxbc_d(j,k,i4,mx) = den*vx
                Qxbc_d(j,k,i4,my) = 0
                Qxbc_d(j,k,i4,mz) = 0
                Qxbc_d(j,k,i4,en) = pres/aindm1 + 0.5*den*vx**2
            end do
            end do
        end do
        print *,'den',den,'pres',pres,'mfr',mflowrate

    end subroutine set_xlobc_inflow


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


    subroutine set_cyl_in_2d_pipe_boundaries(Qmask_i, ux_amb, den, pres, Qxlo_e_c, Qcyl_e_c)
        logical, dimension(nx,ny,nz), intent(in) :: Qmask_i
        real, intent(in) :: ux_amb, den, pres
        real, dimension(ny,nz,nface,nQ), intent(out) :: Qxlo_e_c
        real, dimension(nx,ny,nface,nQ), intent(out) :: Qcyl_e_c

        real yp,yp0,vx
        integer :: j,k,ieq,i4

        yp0 = lyd  ! set zero-value of y-coordinate to domain bottom

        ! apply no-slip BCs for cylinder (set all momenta to 0)
        Qcyl_e_c(:,:,:,:) = 0.0
        do ieq = rh,en
            do i4 = 1,nface
                Qcyl_e_c(:,:,i4,ieq) = Q_r0(:,:,1,i4,ieq)  ! Init to grid values
                where ( Qmask_i(:,:,1) )
                    Qcyl_e_c(:,:,i4,ieq) = 0.0             ! Set dens/momenta to 0
                end where
            end do
        end do

        ! apply parabolic inflow BCs on lower x wall
        Qxlo_e_c(:,:,:,:) = 0.0
        do i4 = 1,nface
            do k = 1,nz
            do j = 1,ny
                yp = yc(j) - yp0
                vx = 4*ux_amb*yp*(ly-yp)/ly**2

                Qxlo_e_c(j,k,i4,rh) = den
                Qxlo_e_c(j,k,i4,mx) = den*vx
                Qxlo_e_c(j,k,i4,my) = 0
                Qxlo_e_c(j,k,i4,mz) = 0
                Qxlo_e_c(j,k,i4,en) = pres/aindm1 + 0.5*den*vx**2
            end do
            end do
        end do
    end subroutine set_cyl_in_2d_pipe_boundaries


    subroutine apply_cyl_in_2d_pipe_boundaries(Qcyl_e_c, Qcyl_e)
        real, dimension(nx,ny,nface,nQ), intent(in) :: Qcyl_e_c
        real, dimension(nx,ny,nface,nQ), intent(inout) :: Qcyl_e  ! TODO: might wanna change

        do ieq=1,nQ
            where (Qmask(:,:,:))
                Q_r0(:,:,:,rh,ieq) = 0.0
                Q_r0(:,:,:,rh,1)   = 1.25
                Q_r0(:,:,:,mx,ieq) = 0.0
                Q_r0(:,:,:,my,ieq) = 0.0
                Q_r0(:,:,:,mz,ieq) = 0.0
            end where
        end do
    end subroutine apply_cyl_in_2d_pipe_boundaries

end module boundary_custom




!##### boundary ################################################################
module boundary

    use parameters
    use helpers, only : mpi_print
    use boundary_defs
    ! use boundary_defs2
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


    !===========================================================================
    ! Initialize arrays to store boundary conditions
    !------------------------------------------------------------
    real, dimension(ny,nz,nface,nQ) :: Qxlo_ext, Qxhi_ext
    real, dimension(ny,nz,nface,nQ) :: Qxlo_int, Qxhi_int
    real, dimension(nx,nz,nface,nQ) :: Qylo_ext, Qyhi_ext
    real, dimension(nx,nz,nface,nQ) :: Qylo_int, Qyhi_int
    real, dimension(nx,ny,nface,nQ) :: Qzlo_ext, Qzhi_ext
    real, dimension(nx,ny,nface,nQ) :: Qzlo_int, Qzhi_int
    !---------------------------------------------------------------------------

contains

    !===========================================================================
    ! Apply boundary conditions (specified by user at runtime)
    !------------------------------------------------------------
    subroutine apply_boundaries
        real ux_scale
        if ( mpi_P == 1 ) then
            call apply_xlobc(Qxlo_int, Qxlo_ext)
            ! if ( t < 0.1*tf ) then
            !     ux_scale = 10.0*t/tf + 0.001
            ! else
            !     ux_scale = 1.0 + 0.001
            ! end if
            ! Qxlo_ext_scale(:,:,:,:) = ux_scale*Qxlo_ext_def(:,:,:,:)
            ! call x_bc_inflow(Qxlo_int, Qxlo_ext_scale, Qxlo_ext)
            call x_bc_inflow(Qxlo_int, Qxlo_ext_def, Qxlo_ext)
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

        ! if (mpi_P == 1) print *,'X velocity:'
        ! if (mpi_P == 1) print *,Qxlo_ext(:,1,1,mx)/Qxlo_ext(:,1,1,rh)
        ! if (mpi_P == 1) print *,''
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

end module boundary
