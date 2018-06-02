!***** PREP_STEP.F90 ******************************************************************
module prep_step

use params
use helpers
use basis

use boundary
use random

! write(*,'(A11,I1,A2,2ES9.1,A3,2ES9.1)') 'Qylo_ext (',iam,'):',Qylo_ext(1:2,1,1,exx),' | ',Qylo_ext(nx-1:nx,1,1,exx)
! write(*,'(A12,I1,A2,2ES9.1,A3,2ES9.1)') 'Qyhi_ext (',iam,'):',Qyhi_ext(1:2,1,1,exx),' | ',Qyhi_ext(nx-1:nx,1,1,exx)
! write(*,'(A11,I1,A2,2ES9.1,A3,2ES9.1)') 'Qylo_int (',iam,'):',Qylo_int(1:2,1,1,exx),' | ',Qylo_int(nx-1:nx,1,1,exx)
! write(*,'(A12,I1,A2,2ES9.1,A3,2ES9.1)') 'Qyhi_int (',iam,'):',Qyhi_int(1:2,1,1,exx),' | ',Qyhi_int(nx-1:nx,1,1,exx)

contains

    !===========================================================================
    ! exchange_flux : set internal BCs from field variables
    !------------------------------------------------------------
    subroutine exchange_flux(Q_in)
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_in
        integer :: ieq, i, j, k, ipnt

        !#########################################################
        ! Step 1a: Get fluxes at the boundaries of this MPI domain
        !---------------------------------------------------------
        do ieq = 1,nQ
        do k = 1,nz
        do j = 1,ny
          do ipnt=1,nface
            Qxlo_int(j,k,ipnt,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_in(1,j,k,ieq,1:nbasis))
            Qxhi_int(j,k,ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_in(nx,j,k,ieq,1:nbasis))
          end do
        end do
        end do
        end do

        do ieq = 1,nQ
        do k = 1,nz
        do i = 1,nx
          do ipnt=1,nface
            Qylo_int(i,k,ipnt,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_in(i,1,k,ieq,1:nbasis))
            Qyhi_int(i,k,ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_in(i,ny,k,ieq,1:nbasis))
          end do
        end do
        end do
        end do

        do ieq = 1,nQ
        do j = 1,ny
        do i = 1,nx
          do ipnt=1,nface
            Qzlo_int(i,j,ipnt,ieq) = sum(bfvals_zm(ipnt,1:nbasis)*Q_in(i,j,1,ieq,1:nbasis))
            Qzhi_int(i,j,ipnt,ieq) = sum(bfvals_zp(ipnt,1:nbasis)*Q_in(i,j,nz,ieq,1:nbasis))
          end do
        end do
        end do
        end do
        !#########################################################

        !#########################################################
        ! Step 1b: Get stochastic fluxes at the boundaries of this MPI domain
        !---------------------------------------------------------
        ! WARNING: indexing might be slow since last two (slow) indices are implicit
        do k = 1,nz
        do j = 1,ny
          do ipnt=1,nface
            Sfxlo(ipnt,j,k,1:3,1:3) = Sflux(ipnt,1,j,k,1:3,1:3)
            ! Sfxhi(ipnt,j,k,1:3,1:3) = Sflux(ipnt,nx1,j,k,1:3,1:3)
          end do
        end do
        end do

        do k = 1,nz
        do i = 1,nx
          do ipnt=1,nface
            Sfylo(ipnt,i,k,1:3,1:3) = Sflux(ipnt,i,1,k,1:3,1:3)
            ! Sfyhi(ipnt,i,k,1:3,1:3) = Sflux(ipnt,i,ny1,k,1:3,1:3)
          end do
        end do
        end do

        do j = 1,ny
        do i = 1,nx
          do ipnt=1,nface
            Sfzlo(ipnt,i,j,1:3,1:3) = Sflux(ipnt,i,j,1,1:3,1:3)
            ! Sfzhi(ipnt,i,j,1:3,1:3) = Sflux(ipnt,i,j,nz1,1:3,1:3)
          end do
        end do
        end do
        !#########################################################

        !#########################################################
        ! Step 2a: Exchange fluxes with (neighboring) MPI domains
        !---------------------------------------------------------
        call perform_flux_exchange

        ! Step 2b: Exchange fluxes with (neighboring) MPI domains
        !---------------------------------------------------------
        call perform_Sflux_exchange
        !#########################################################
    end subroutine exchange_flux
    !---------------------------------------------------------------------------


    !===========================================================================
    ! perform_flux_exchange : exchange fluxes between MPI domains
    !------------------------------------------------------------
    subroutine perform_flux_exchange
        integer :: mpi_size

        call MPI_BARRIER(cartcomm,ierr)

        !---------------------------------------------
        mpi_size = nface*ny*nz*nQ

        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Qxhi_int,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(1),ierr)
        endif
        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Qxlo_ext,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(2),ierr)
        endif

        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(1),stats(:,1),ierr)
            call MPI_IRecv(Qxhi_ext,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(3),ierr)
        endif
        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(2),stats(:,2),ierr)
            call MPI_ISend(Qxlo_int,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(4),ierr)
        endif

        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(3),stats(:,3),ierr)
        endif
        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(4),stats(:,4),ierr)
        endif

        !---------------------------------------------
        mpi_size = nface*nx*nz*nQ

        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Qyhi_int,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(1),ierr)
        endif
        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Qylo_ext,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(2),ierr)
        endif

        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(1),stats(:,1),ierr)
            call MPI_IRecv(Qyhi_ext,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(3),ierr)
        endif
        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(2),stats(:,2),ierr)
            call MPI_ISend(Qylo_int,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(4),ierr)
        endif

        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(3),stats(:,3),ierr)
        endif
        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(4),stats(:,4),ierr)
        endif

        !---------------------------------------------
        mpi_size = nface*nx*ny*nQ

        if (nbrs(UP) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Qzhi_int,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(1),ierr)
        endif
        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Qzlo_ext,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(2),ierr)
        endif

        if (nbrs(UP) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(1),stats(:,1),ierr)
            call MPI_IRecv(Qzhi_ext,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(3),ierr)
        endif
        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(2),stats(:,2),ierr)
            call MPI_ISend(Qzlo_int,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(4),ierr)
        endif

        if (nbrs(UP) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(3),stats(:,3),ierr)
        endif
        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(4),stats(:,4),ierr)
        endif

    end subroutine perform_flux_exchange
    !---------------------------------------------------------------------------

    !===========================================================================
    ! perform_Sflux_exchange : exchange fluxes between MPI domains
    !------------------------------------------------------------
    subroutine perform_Sflux_exchange
        integer :: mpi_size

        call MPI_BARRIER(cartcomm,ierr)

        !---------------------------------------------
        mpi_size = nface*ny*nz*9

        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Sfxlo,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(1),ierr)
        endif
        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Sfxhi,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(2),ierr)
        endif

        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(1),stats(:,1),ierr)
        endif
        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(2),stats(:,2),ierr)
        endif

        !---------------------------------------------
        mpi_size = nface*nx*nz*9

        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Sfylo,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(1),ierr)
        endif
        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Sfyhi,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(2),ierr)
        endif

        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(1),stats(:,1),ierr)
        endif
        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(2),stats(:,2),ierr)
        endif

        !---------------------------------------------
        mpi_size = nface*nx*ny*9

        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Sfzlo,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(1),ierr)
        endif
        if (nbrs(UP) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Sfzhi,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(2),ierr)
        endif

        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(1),stats(:,1),ierr)
        endif
        if (nbrs(UP) .ne. MPI_PROC_NULL) then
            call MPI_Wait(reqs(2),stats(:,2),ierr)
        endif

        ! Set "hi" (upper) values of Sflux for THIS MPI domain based on the
        ! "lo" (lower) values sent by the WEST/SOUTH/DOWN domain:
        do k = 1,nz
        do j = 1,ny
          do ipnt=1,nface
            Sflux(ipnt,nx1,j,k,1:3,1:3) = Sfxhi(ipnt,j,k,1:3,1:3)
          end do
        end do
        end do

        do k = 1,nz
        do i = 1,nx
          do ipnt=1,nface
            Sflux(ipnt,i,ny1,k,1:3,1:3) = Sfyhi(ipnt,i,k,1:3,1:3)
          end do
        end do
        end do

        do j = 1,ny
        do i = 1,nx
          do ipnt=1,nface
            Sflux(ipnt,i,j,nz1,1:3,1:3) = Sfzhi(ipnt,i,j,1:3,1:3)
          end do
        end do
        end do

    end subroutine perform_Sflux_exchange
    !---------------------------------------------------------------------------

end module prep_step
