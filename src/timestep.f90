module timestep

use params
use spatial

implicit none

real :: cflm  ! Initialized in init.f90

contains

    !===========================================================================
    ! Get the smallest time step required by CFL condition
    !------------------------------------------------------------
    real function get_min_dt(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r

        real dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0
        real vex,vey,vez,vem,vem0,dni,dn,vx,vy,vz,Pr,sqdni,vacc,vacc0,cs
        integer :: i,j,k,main_proc=0,mpi_size=1
        integer :: loc_reqs(numprocs-1),loc_stats(MPI_STATUS_SIZE,numprocs-1)
        ! NOTE: uses the global variable cflm (set in initialize.f90)

        vmag = 0.

        do k=1,nz
        do j=1,ny
        do i=1,nx
            dn = Q_r(i,j,k,rh,1)
            dni = 1./dn
            vx = Q_r(i,j,k,mx,1)*dni
            vy = Q_r(i,j,k,my,1)*dni
            vz = Q_r(i,j,k,mz,1)*dni
            if (ieos == 1) cs = sqrt(aindex*(Q_r(i,j,k,en,1)*dni - 0.5*(vx**2 + vy**2 + vz**2)))
            if (ieos == 2) cs = sqrt(7.2*P_1*dn**6.2 + T_floor)

            vmag0 = max( abs(vx)+cs, abs(vy)+cs, abs(vz)+cs )
            if (vmag0 > vmag .and. dn > rh_mult*rh_floor) vmag = vmag0  ! NOTE: from newCES (excluded dn thing)
        end do
        end do
        end do

        vmax = (vmag + cs)*dxi  ! NOTE: from newCES  (was vmag*dxi)
        dt_min = cflm/vmax  ! time step determined by maximum flow + sound speed in the domain

        call MPI_BARRIER(cartcomm,ierr)
        if (iam == main_proc) then
            do i=1,numprocs-1
                call MPI_IRecv(dt_val(i),mpi_size,MPI_TT,i,0,cartcomm,loc_reqs(i),ierr)
            enddo
            call MPI_WaitAll(numprocs-1,loc_reqs,loc_stats,ierr)
            do i=1,numprocs-1
                dt_min = min(dt_min,dt_val(i))
            enddo
            do i=1,numprocs-1
                call MPI_ISend(dt_min,mpi_size,MPI_TT,i,0,cartcomm,loc_reqs(i),ierr)
            enddo
            call MPI_WaitAll(numprocs-1,loc_reqs,loc_stats,ierr)
        else
            call MPI_ISend(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
            call MPI_Wait(reqs(1),stats(:,1),ierr)
            call MPI_IRecv(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
            call MPI_Wait(reqs(1),stats(:,1),ierr)
        end if

        get_min_dt = dt_min
        return
    end function get_min_dt
    !---------------------------------------------------------------------------

end module timestep
