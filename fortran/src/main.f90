program main

! use input
! use parameters
! use helpers
use initialize

use ic_mod
use prepare_time_advance_mod
use io_mod
use basis_funcs_mod

! implicit none

call perform_setup()

t1 = get_clock_time()
call output_vtk(Q_r0,nout,iam)

do while (t<tf)

    dt = get_min_dt()
    dti = 1./dt
    sqrt_dVdt_i = (dVi*dti)**0.5  ! can move dVi into eta_sd to remove a multiplication

    select case(iorder)
        case(2)
            call update_2nd_order(Q_r0, Q_r1, Q_r2)
        case(3)
            call update_3rd_order(Q_r0, Q_r1, Q_r2, Q_r3)
    end select

    t = t + dt

    if (t .gt. dtout*nout) then
        call generate_output()
    end if

end do

!-------------------------------------------------------------------------------
! RNG is completed with de-allocation of system resources:
vsl_errcode = vsldeletestream( vsl_stream )

call MPI_Finalize(ierr)

if (iam .eq. print_mpi) then
    print *, 'Wall time', (get_clock_time() - t_start), 'seconds'
end if
!-------------------------------------------------------------------------------

contains


    !-------------------------------------------------------------------------------
    subroutine update_2nd_order(Q_io, Q_1, Q_2)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2

        call euler_step(Q_io, Q_1)
        call euler_step(Q_1, Q_2)
        Q_io = 0.5 * ( Q_io + Q_2 )

    end subroutine update_2nd_order


    subroutine update_3rd_order(Q_io, Q_1, Q_2, Q_3)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2, Q_3

        call euler_step(Q_io, Q_1)
        call euler_step(Q_1, Q_2)
        Q_3 = 0.75*Q_io + 0.25*Q_2

        call euler_step(Q_3, Q_2)  ! re-use the second array
        Q_io = c1d3 * ( Q_io + 2.0*Q_2 )

    end subroutine update_3rd_order
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    subroutine euler_step(Q_in, Q_out)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out

        call prep_advance(Q_in)
        call calc_rhs(Q_in)
        call advance_time_level_gl(Q_in, Q_out)

    end subroutine euler_step


    subroutine prep_advance(Q_io)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        if(ieos .eq. 1) call limiter(Q_io)  ! added in (from "viscosity" version)
        if(ieos .eq. 2) call limiter(Q_io)
        call prepare_exchange(Q_io)
        call set_bc2

    end subroutine prep_advance

    subroutine calc_rhs(Q_io)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        call flux_cal(Q_io)
        call innerintegral(Q_io)
        call glflux  ! glflux currently breaks after "bug fix"
        call source_calc(Q_io)

    end subroutine calc_rhs


    subroutine advance_time_level_gl(Q_in, Q_out)

        implicit none
        integer i,j,k,ieq,ir
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

            do ieq = 1,nQ
            do ir=1,nbasis
                Q_out(i,j,k,ieq,ir) =                                            &
                        Q_in(i,j,k,ieq,ir) - dt*( glflux_r(i,j,k,ieq,ir)        &
                                                - source_r(i,j,k,ieq,ir) )
            end do
            end do

            do ieq = 1,nQ
                if ( Q_out(i,j,k,ieq,1) .ne. Q_out(i,j,k,ieq,1)) then
                    print *,'NaN. Bailing out...'
                    print *,'  xc  =',xc(i),'  yc  =',yc(j),'  zc  =',zc(k),'  ieq  =',ieq
                    call exit(-1)
                endif
            end do

        end do
        end do
        end do

    end subroutine advance_time_level_gl
    !-------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------

    real function get_min_dt()
        implicit none

        real dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0
        real vex,vey,vez,vem,vem0,dni,dn,vx,vy,vz,Pr,sqdni,vacc,vacc0,cs
        integer :: i,j,k,main_proc=0,mpi_size=1
        integer :: loc_reqs(numprocs-1),loc_stats(MPI_STATUS_SIZE,numprocs-1)

        vmag = 0.

        do k=1,nz
            do j=1,ny
                do i=1,nx

                    dn = Q_r0(i,j,k,rh,1)
                    dni = 1./dn
                    vx = Q_r0(i,j,k,mx,1)*dni
                    vy = Q_r0(i,j,k,my,1)*dni
                    vz = Q_r0(i,j,k,mz,1)*dni
                    if(ieos .eq. 1) cs = sqrt(aindex*(Q_r0(i,j,k,en,1)*dni - 0.5*(vx**2 + vy**2 + vz**2)))
                    if(ieos .eq. 2) cs = sqrt(7.2*P_1*dn**6.2 + T_floor)

                    vmag0 = max(abs(vx)+cs, abs(vy)+cs, abs(vz)+cs)
                    if(vmag0 > vmag) vmag = vmag0

                end do
            end do
        end do

        vmax = vmag*dxi
        dt_min = 1*cflm/vmax  ! time step is determined by the maximum flow + sound speed in the domain

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

    end function get_min_dt

!----------------------------------------------------------------------------------------------

    subroutine generate_output()

        implicit none
        integer ioe

        nout = nout+1
        if (iam .eq. print_mpi) then
            print *, 'nout = ', nout
            print *, '   t = ',t*100.,'         dt= ',dt
            t2 = get_clock_time()
            print *, '  >> Iteration time', (t2-t1), 'seconds'
            t1 = t2
        end if

        call MPI_BARRIER(cartcomm,ierr)
        call output_vtk(Q_r0,nout,iam)

        ! write checkpoint files; assign an odd/even id to ensure last two sets are kept
        if (iwrite .eq. 1) then
            ioe = 2 - mod(nout,2)
            call writeQ(fpre,iam,ioe,Q_r0,t,dt,nout,mpi_nx,mpi_ny,mpi_nz)
        end if

        call MPI_BARRIER(cartcomm,ierr)

        if (iam .eq. print_mpi) then
            t2 = get_clock_time()
            print *, '  >> Output time', (t2-t1), 'seconds'
            print *, ''
            t1 = t2
        end if

    end subroutine generate_output

end program main
