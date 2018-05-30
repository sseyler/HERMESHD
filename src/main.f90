!****** MAIN.F90 *************************************************************************
program main

use input!, only : nx,ny,nz
use params!, only : nQ,nbasis,Q_r0,Q_r1,Q_r2,Q_r3
use helpers!, only : xc,yc,zc,get_clock_time
use spatial
use timestep

use integrator
use boundary_custom
use boundary
use initialcon
use basis

use initialize

use prep_step
use sources
use random
use flux
use output

integer :: nout, comm

!###############################################################################
! I. SETUP
!----------------------------------------------------------------

!-------------------------------------------------
! 1. Initialize general simulation variables
!-------------------------------------------------
call MPI_Init ( ierr )
call initializer(t, dt, nout, MPI_COMM_WORLD)

t_start = get_clock_time()  ! start timer for wall time
dtout = tf/ntout  ! TODO: move this to a more sensible place once output scheme is improved!

!-------------------------------------------------
! 2. Select hydrodynamic model (equations)
!-------------------------------------------------
call select_hydro_model(ivis, Fpt_x, Fpt_y, Fpt_z, Fpts_x, Fpts_y, Fpts_z)  ! DEBUG

!-------------------------------------------------
! 3. Select and set initial conditions
!-------------------------------------------------
if (iread == 0) then
    call set_ic(Q_r0, icid)
else
    call set_ic_from_file(Q_r0, t, dt, dtout, nout)
endif

!-------------------------------------------------
! 4. Select integration method
!-------------------------------------------------
call select_integrator(iname, update)

!-------------------------------------------------
! 5. Select boundary conditions
!-------------------------------------------------
call select_x_boundaries(xlobc, xhibc, apply_xlobc, apply_xhibc)
call select_y_boundaries(ylobc, yhibc, apply_ylobc, apply_yhibc)
call select_z_boundaries(zlobc, zhibc, apply_zlobc, apply_zhibc)

t1 = get_clock_time()

!-------------------------------------------------
! 6. Generate initial output
!-------------------------------------------------
call output_vtk(Q_r0, nout, iam)

!-------------------------------------------------------------------------------


!###############################################################################
! II. SIMULATION
!----------------------------------------------------------------
do while( t < tf )

    dt = get_min_dt(Q_r0)  ! WARNING, this is a HACK
    call update(Q_r0, Q_r1, Q_r2, dt)
    t = t + dt

    call generate_output(Q_r0, t, nout)  ! determines when output should be generated

end do
!-------------------------------------------------------------------------------


!###############################################################################
! III. CLEANUP
!----------------------------------------------------------------

!-------------------------------------------------
! 1. De-allocate system resources for RNG
!-------------------------------------------------
if (llns) call random_cleanup()

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
    ! Generate console and VTK output
    !------------------------------------------------------------
    subroutine generate_output(Q_r, t, nout)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r
        real, intent(in) :: t
        integer, intent(inout) :: nout

        integer ioe

        ! TODO: dtout may be deprecated once improved output scheme is used
        ! TODO: consider using init_temporal_params() to initialize the
        !       dtout-type parameters for each dynamical quantity

        ! dtdout  = tf/nstdout
        ! dtldout = tf/nstldout
        ! dtvout  = tf/nstvout
        ! dttout  = tf/nsttout
        ! dteiout = tf/nsteiout
        ! dtenout = tf/nstenout
        ! dtesout = tf/nstesout
        ! dtpout  = tf/nstpout
        ! dtstout = tf/nststout
        ! dtvrout = tf/nstvrout
        !
        ! if (t > dtdout*ndout) then
        !     call generate_output(Q_r0, nout)
        ! end if

        if (t > dtout*nout) then

            nout = nout + 1
            if (iam == print_mpi) then
                print *, 'nout = ', nout
                print *, '   t = ',t,'         dt= ',dt
                t2 = get_clock_time()
                print *, '  >> Iteration time', (t2-t1), 'seconds'
                t1 = t2
            end if

            call MPI_BARRIER(cartcomm,ierr)
            call output_vtk(Q_r,nout,iam)

            ! write checkpoint files; assign an odd/even id to ensure last two sets are kept
            if (iwrite == 1) then
                ioe = 2 - mod(nout,2)
                call writeQ(fpre,iam,ioe,Q_r,t,dt,nout,mpi_nx,mpi_ny,mpi_nz)
            end if

            call MPI_BARRIER(cartcomm,ierr)

            if (iam == print_mpi) then
                t2 = get_clock_time()
                print *, '  >> Output time', (t2-t1), 'seconds'
                print *, ''
                t1 = t2
            end if

            ! if (mpi_P == 1) print *,Qxlo_ext_def(:,1,1,mx)

        end if
    end subroutine generate_output
    !---------------------------------------------------------------------------

end program main
