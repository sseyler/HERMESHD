program main

use input!, only : nx,ny,nz
use parameters!, only : nQ,nbasis,Q_r0,Q_r1,Q_r2,Q_r3
use helpers!, only : xc,yc,zc,get_clock_time

use integrator
use boundary_custom
use boundary
use initialcon
use basis_funcs

use initialize

use prepare_step
use sources
use random
use flux
use output


!###############################################################################
! I. SETUP
!----------------------------------------------------------------

!-------------------------------------------------
! 1. Initialize general simulation variables
!-------------------------------------------------
call perform_setup

t_start = get_clock_time()  ! start timer for wall time

!-------------------------------------------------
! 2. Select and set initial conditions
!-------------------------------------------------
if (iread == 0) then
    call set_ic(Q_r0, icid)
else
    call initialize_from_file(Q_r0)
endif

!-------------------------------------------------
! 3. Select integration method
!-------------------------------------------------
call select_integrator(iname, update)

!-------------------------------------------------
! 4. Select boundary conditions
!-------------------------------------------------
call select_x_boundaries(xlobc, xhibc, apply_xlobc, apply_xhibc)
call select_y_boundaries(ylobc, yhibc, apply_ylobc, apply_yhibc)
call select_z_boundaries(zlobc, zhibc, apply_zlobc, apply_zhibc)

t1 = get_clock_time()

!-------------------------------------------------
! 5. Generate initial output
!-------------------------------------------------
call output_vtk(Q_r0, nout, iam)

!-------------------------------------------------------------------------------


!###############################################################################
! II. SIMULATION
!----------------------------------------------------------------
do while( t < tf )

    dt = get_min_dt(Q_r0)
    call update(Q_r0, Q_r1, Q_r2)
    t = t + dt

    if (t > dtout*nout) then
        call generate_output(Q_r0, nout)
    end if

end do
!-------------------------------------------------------------------------------


!###############################################################################
! III. CLEANUP
!----------------------------------------------------------------

!-------------------------------------------------
! 1. De-allocate system resources for RNG
!-------------------------------------------------
vsl_errcode = vsldeletestream( vsl_stream )

!-------------------------------------------------
! 2. MPI cleanup
!-------------------------------------------------
call MPI_Finalize(ierr)

!-------------------------------------------------
! 3. Report simulation wall time
!-------------------------------------------------
t_stop = get_clock_time()  ! start timer for wall time
call report_wall_time(iam, t_stop-t_start)

!-------------------------------------------------------------------------------


contains

    !===========================================================================
    ! Get the smallest time step required by CFL condition
    !------------------------------------------------------------
    real function get_min_dt(Q_r)
        ! uses the global variable cflm (set in initialize.f90)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r

        real dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0
        real vex,vey,vez,vem,vem0,dni,dn,vx,vy,vz,Pr,sqdni,vacc,vacc0,cs
        integer :: i,j,k,main_proc=0,mpi_size=1
        integer :: loc_reqs(numprocs-1),loc_stats(MPI_STATUS_SIZE,numprocs-1)

        vmag = 0.

        do k=1,nz
        do j=1,ny
        do i=1,nx
            dn = Q_r(i,j,k,rh,1)
            dni = 1./dn
            vx = Q_r(i,j,k,mx,1)*dni
            vy = Q_r(i,j,k,my,1)*dni
            vz = Q_r(i,j,k,mz,1)*dni
            if(ieos == 1) cs = sqrt(aindex*(Q_r(i,j,k,en,1)*dni - 0.5*(vx**2 + vy**2 + vz**2)))
            if(ieos == 2) cs = sqrt(7.2*P_1*dn**6.2 + T_floor)

            vmag0 = max( abs(vx)+cs, abs(vy)+cs, abs(vz)+cs )
            if(vmag0 > vmag .and. dn > rh_mult*rh_floor) vmag = vmag0  ! NOTE: from newCES (excluded dn thing)
        end do
        end do
        end do

        vmax = (vmag + cs)*dxi  ! NOTE: from newCES  (was vmag*dxi)
        dt_min = cflm/vmax  ! time step determined by maximum flow + sound speed in the domain

        call MPI_BARRIER(cartcomm,ierr)
        if (iam.eq.main_proc) then
            do i=1,numprocs-1
                call MPI_IRecv(dt_val(i),mpi_size,MPI_TT,i,0,cartcomm,loc_reqs(i),ierr)
            enddo
            call MPI_WaitAll(numprocs-1,loc_reqs,loc_stats,ierr)
            do i=1,numprocs-1
                dt_min=min(dt_min,dt_val(i))
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
        endif

        get_min_dt = dt_min
        return
    end function get_min_dt
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Generate console and VTK output
    !------------------------------------------------------------
    subroutine generate_output(Q_r, nout)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r
        integer, intent(inout) :: nout
        integer ioe

        nout = nout + 1
        if (iam .eq. print_mpi) then
            print *, 'nout = ', nout
            print *, '   t = ',t*100.,'         dt= ',dt
            t2 = get_clock_time()
            print *, '  >> Iteration time', (t2-t1), 'seconds'
            t1 = t2
        end if

        call MPI_BARRIER(cartcomm,ierr)
        call output_vtk(Q_r,nout,iam)

        ! write checkpoint files; assign an odd/even id to ensure last two sets are kept
        if (iwrite) then
            ioe = 2 - mod(nout,2)
            call writeQ(fpre,iam,ioe,Q_r,t,dt,nout,mpi_nx,mpi_ny,mpi_nz)
        end if

        call MPI_BARRIER(cartcomm,ierr)

        if (iam .eq. print_mpi) then
            t2 = get_clock_time()
            print *, '  >> Output time', (t2-t1), 'seconds'
            print *, ''
            t1 = t2
        end if
    end subroutine generate_output
    !---------------------------------------------------------------------------

end program main
