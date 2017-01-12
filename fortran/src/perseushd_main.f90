program main

use parameters_mod
use ic_mod
use prepare_time_advance_mod
use io_mod
use basis_funcs_mod

! include '/nfs/packages/opt/Linux_x86_64/openmpi/1.6.3/intel13.0/include/mpif.h'


    if (nbasis .le. 8) cflm = 0.14
    if (nbasis .eq. 27) cflm = 0.1
    if (nbasis .eq. 64) cflm = 0.08

    ! Initialize grid sizes and local lengths

    cbasis(1) = 1.              ! coeff for basis func {1}
    cbasis(kx:kz) = 3.          ! coeff for basis funcs {x,y,z}
    cbasis(kyz:kxy) = 9.        ! coeff for basis funcs {yz,zx,xy}
    cbasis(kxyz) = 27.          ! coeff for basis func {xyz}

    cbasis(kxx:kzz) = 5.        ! coeff for basis funcs {P2(x),P2(y),P2(z)}
    cbasis(kyzz:kxyy) = 15.     ! coeff for basis funcs {yP2(z),zP2(x),xP2(y)}
    cbasis(kyyz:kxxy) = 15.     ! coeff for basis funcs {P2(y)z,P2(z)y,P2(z)x}
    cbasis(kyyzz:kxxyy) = 25.   ! coeff for basis funcs {P2(y)P2(z),P2(z)P2(x),P2(x)P2(y)}
    cbasis(kyzxx:kxyzz) = 45.   ! coeff for basis funcs {yzP_2(x),zxP_2(y),xyP_2(z)}
    cbasis(kxyyzz:kzxxyy) = 75. ! coeff for basis funcs {xP2(y)P2(z),yP2(z)P2(x),zP2(x)P2(y)}
    cbasis(kxxyyzz) = 125.      ! coeff for basis funcs {P2(x)P2(y)P2(z)}

    call MPI_Init ( ierr )
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    mpi_nz=numprocs/(mpi_nx*mpi_ny)

    dims(1)=mpi_nx
    dims(2)=mpi_ny
    dims(3)=mpi_nz

    periods(:)=0
    if (xhbc .eq. 2) then
        periods(1)=1
    end if
    if (yhbc .eq. 2) then
        periods(2)=1
    end if
    if (zhbc .eq. 2) then
        periods(3)=1
    end if
    reorder = 1

    call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, reorder,cartcomm, ierr)
    call MPI_COMM_RANK (cartcomm, iam, ierr )
    call MPI_CART_COORDS(cartcomm, iam, 3, coords, ierr)
    mpi_P=coords(1)+1
    mpi_Q=coords(2)+1
    mpi_R=coords(3)+1
    call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(WEST), nbrs(EAST), ierr)
    call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(SOUTH), nbrs(NORTH), ierr)
    call MPI_CART_SHIFT(cartcomm, 2, 1, nbrs(DOWN), nbrs(UP), ierr)

    half_length = lx/2.0
    lxd = -half_length
    lxu = half_length
    half_length = ly/2.0
    lyd = -half_length
    lyu = half_length
    half_length = lz/2.0
    lzd = -half_length
    lzu = half_length

    dxi = (nx*mpi_nx)/(lxu-lxd)
    dyi = (ny*mpi_ny)/(lyu-lyd)
    dzi = (nz*mpi_nz)/(lzu-lzd)

    dx = 1./dxi
    dy = 1./dyi
    dz = 1./dzi
    loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
    loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
    loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz

    rh_fluid = 1.

    ! indices used in poynting() for computing Poynting Maxwell flux

    mxa(1) = mx
    mxa(2) = my
    mxa(3) = mz
    mya(1) = my
    mya(2) = mz
    mya(3) = mx
    mza(1) = mz
    mza(2) = mx
    mza(3) = my

    t = 0.
    dt = cflm*dx/clt
    dtoriginal = dt
    nout = 0
    niter = 0
    dtout = tf/ntout


    !===============================================================================
    ! NOTE: this is new stuff!
    ! Volume of grid cells (for stochastic forcing terms)
    !   --> this should be chosen to be consistent with the spatial discretization
    !       i.e. stochastic terms are 4 points (iquad == 2) at each cell face, which
    !       are computed using Gaussian quadrature from the internal points

    ! A spatial discretization that gives the correct equilibrium discrete
    ! covariance is said to satisfy the discrete fluctuation-dissipation balance
    ! (DFDB) condition [1, 2].
    ! The condition guarantees that for sufficiently small time steps the statistics
    ! of the discrete fluctuations are consistent with the continuum formulation.
    ! For larger time steps, the difference between the discrete and continuum
    ! covariances will depend on the order of weak accuracy of the temporal
    ! discretization [3].
    !
    ! [1] Donev, et al. On the accurracy... (2010)
    ! [2] Delong, et al. Temporal integrators for... (2013)
    ! [3] Mattingly,, et al. Convergence of numerical... (2010)
    !---------------------------------------------------------------------------
    dV = dx*dy*dz
    dVi = dxi*dyi*dzi
    dti = 1./dt
    sqrt_dVdt_i = (dVi*dti)**0.5  ! Inv sq-root of (dV*dt), dV = grid cell volume
    !===============================================================================


    ! Evaluate local cell values of basis functions on cell interior and faces.
    ! This is done for 1, 2, or 3 point Gaussian quadrature.
    call set_bfvals_3D

    do ir=1,nbasis
        do ipg=1,npg
            bval_int_wgt(ipg,ir) = wgt3d(ipg)*bfvals_int(ipg,ir)
        end do
    end do

    do ir=1,nbasis
        wgtbfvals_xp(1:nface,ir) = wgt2d(1:nface)*bfvals_xp(1:nface,ir)
        wgtbfvals_yp(1:nface,ir) = wgt2d(1:nface)*bfvals_yp(1:nface,ir)
        wgtbfvals_zp(1:nface,ir) = wgt2d(1:nface)*bfvals_zp(1:nface,ir)
        wgtbfvals_xm(1:nface,ir) = wgt2d(1:nface)*bfvals_xm(1:nface,ir)
        wgtbfvals_ym(1:nface,ir) = wgt2d(1:nface)*bfvals_ym(1:nface,ir)
        wgtbfvals_zm(1:nface,ir) = wgt2d(1:nface)*bfvals_zm(1:nface,ir)
    end do

    do ir=1,nbasis
        wgtbf_xmp(1:nface,1,ir) = -0.25*cbasis(ir)*dxi*wgtbfvals_xm(1:nface,ir)
        wgtbf_ymp(1:nface,1,ir) = -0.25*cbasis(ir)*dyi*wgtbfvals_ym(1:nface,ir)
        wgtbf_zmp(1:nface,1,ir) = -0.25*cbasis(ir)*dzi*wgtbfvals_zm(1:nface,ir)
        wgtbf_xmp(1:nface,2,ir) = 0.25*cbasis(ir)*dxi*wgtbfvals_xp(1:nface,ir)
        wgtbf_ymp(1:nface,2,ir) = 0.25*cbasis(ir)*dyi*wgtbfvals_yp(1:nface,ir)
        wgtbf_zmp(1:nface,2,ir) = 0.25*cbasis(ir)*dzi*wgtbfvals_zp(1:nface,ir)
    end do

    call init_random_seed(123456789)
    iseed = 1317345*mpi_P + 5438432*mpi_Q + 38472613*mpi_R


    !===============================================================================
    ! Initialize MKL random number generator
    !---------------------------------------------------------------------------
    ! NOTE: might want to move random number generator initialzation into subroutine
    vsl_errcode = vslnewstream(vsl_stream, vsl_brng, iseed)
    !===============================================================================


    if (iam .eq. print_mpi) then
        print *, ''
        print *, '----------------------------------------------'
        print *, 'Starting simulation...'
        print *, '----------------------------------------------'
        print *, 'total dim= ',mpi_nx*nx,mpi_ny*ny,mpi_nz*nz
        print *, 'mpi dim= ',mpi_nx,mpi_ny,mpi_nz
        print *, 'te0 is: ', te0
        print *, 'dx is: ', ly/(ny*mpi_ny)*L0
        print *, 'iquad is: ',  iquad
        print *, 'nbasis is: ', nbasis
        print *, '----------------------------------------------'
        print *, ''
    end if

    call system_clock(ticks, count_rate, count_max)
    t_start = ticks*1./count_rate

    if (iread .eq. 0) then

        call set_ic(Q_r0, icid)

    else

        ! This applies only if the initial data are being read from an input file.
        ! - If resuming a run, keep the previous clock (i.e., t at nout) running.
        ! - If not resuming a run, treat input as initial conditions at t=0, nout=0.

        call readQ(fpre,iam,iread,Q_r0,t_p,dt_p,nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p)

        if (resuming) then
            t = t_p
            dt = dt_p
            nout = nout_p
        end if
        ! Note, nout=1 corresponds to t=dt, but nout=2 corresponds to t~dtout, etc.
        if (nout .gt. 1) then
            dtout_p = t_p/(nout_p-1)
        else  ! Automatically pass consistency check
            dtout_p = dtout
        end if
        if (iam .eq. print_mpi) then
            print *, 'resuming = ', resuming
            print *, 't = ', t
            print *, 'dt = ', dt
            print *, 'nout = ', nout
            print *, 'dtout_p = ', dtout_p, ' dtout = ', dtout
            print *, 'mpi_nx_p = ', mpi_nx_p, ' mpi_nx = ', mpi_nx
            print *, 'mpi_ny_p = ', mpi_ny_p, ' mpi_ny = ', mpi_ny
            print *, 'mpi_nz_p = ', mpi_nz_p, ' mpi_nz = ', mpi_nz
        end if
            ! Quit if dtout is incompatible with input t/(nout-1)
        if (abs(dtout_p-dtout)/dt_p > 1.01) then
            if (iam .eq. print_mpi) then
                print *, 'Bad restart, non-matching dtout'
            end if
            call exit(-1)
        end if
        if ((mpi_nx_p .ne. mpi_nx) .or. (mpi_ny_p .ne. mpi_ny) .or. (mpi_nz_p .ne. mpi_nz)) then
            if (iam .eq. print_mpi) then
                print *, 'Bad restart, non-matching mpi_nx, mpi_ny, or mpi_nz'
            end if
            call exit(-1)
        end if

    end if


    call system_clock( ticks, count_rate, count_max )
    t1 = 1.*ticks / count_rate
    call output_vtk(Q_r0,nout,iam)

    do while (t<tf)

        ! if (mod(niter,200) .eq. 0 .and. iam .eq. print_mpi) then
            ! print *,'niter,t,dt = ',niter,t,dt,dtout*nout
        ! end if
        niter = niter + 1
        call get_min_dt(dt)
        dti = 1./dt
        sqrt_dVdt_i = (dVi*dti)**0.5  ! can move dVi into eta_sd to remove a multiplication

        if (iorder .eq. 2) then
            call prep_advance(Q_r0)
            call advance_time_level_gl(Q_r0,Q_r1)

            call prep_advance(Q_r1)
            call advance_time_level_gl(Q_r1,Q_r2)

            Q_r0 = 0.5*(Q_r0 + Q_r2)
        end if

        if (iorder .eq. 3) then
            call prep_advance(Q_r0)
            call advance_time_level_gl(Q_r0,Q_r1)

            call prep_advance(Q_r1)
            call advance_time_level_gl(Q_r1,Q_r2)

            Q_r3 = 0.75*Q_r0 + 0.25*Q_r2

            call prep_advance(Q_r3)
            call advance_time_level_gl(Q_r3,Q_r2)

            Q_r0 = (Q_r0 + 2.*Q_r2)/3.
        end if

        t = t+dt

        if (t .gt. dtout*nout) then

            nout = nout+1
            if (iam .eq. print_mpi) then
                print *, 'nout = ', nout
                print *, '   t = ',t*100.,'         dt= ',dt

                call system_clock(ticks, count_rate, count_max)
                t2 = 1.*ticks/count_rate
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
                call system_clock(ticks, count_rate, count_max)
                t2 = 1.*ticks/count_rate
                print *, '  >> Output time', (t2-t1), 'seconds'
                print *, ''
                t1 = t2
            end if

        end if

    end do

    !-------------------------------------------------------------------------------
    ! RNG is completed with de-allocation of system resources:
    vsl_errcode = vsldeletestream( vsl_stream )
    !-------------------------------------------------------------------------------

    call MPI_Finalize (ierr)

    if (iam .eq. print_mpi) then
        call system_clock(ticks, count_rate, count_max)
        print *, 'Wall time', (1.*ticks/count_rate - t_start), 'seconds'
    end if


contains


    !------------------------------------------
    subroutine init_random_seed(iseed)

        implicit none
        integer :: i, n, clock,iseed
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        if (iseed .eq. 0) then
            seed =  clock*(iam+1) + 37 * (/ (i - 1, i = 1, n) /)
        else
            seed =  iseed*(iam+1)
        endif
        call random_seed(put = seed)

        ! print *,seed(1)
        deallocate(seed)

    end subroutine


!-------------------------------------------------------------------------------

    subroutine get_min_dt(dt)

       real dt,dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0
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

        dt = dt_min

    end subroutine get_min_dt

!----------------------------------------------------------------------------------------------

    subroutine prep_advance(Q_ri)

        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_ri

        if(ieos .eq. 1) call limiter(Q_ri)  ! added in (from "viscosity" version)
        if(ieos .eq. 2) call limiter(Q_ri)
        call prepare_exchange(Q_ri)
        call set_bc
        call flux_cal(Q_ri)
        call innerintegral(Q_ri)
        call glflux  ! glflux currently breaks after "bug fix"
        call source_calc(Q_ri)

    end subroutine prep_advance


!----------------------------------------------------------------------------------------------

    subroutine advance_time_level_gl(Q_ri, Q_rp)

        implicit none
        integer i,j,k,ieq,ir
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_ri
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_rp
        real Q_xtemp, Q_ytemp, Q_ztemp

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

            do ieq = 1,nQ
                do ir=1,nbasis
                    Q_rp(i,j,k,ieq,ir) =                                            &
                            Q_ri(i,j,k,ieq,ir) - dt*( glflux_r(i,j,k,ieq,ir)        &
                                                    - source_r(i,j,k,ieq,ir) )
                end do
            end do

            do ieq = 1,nQ
                if ( Q_rp(i,j,k,ieq,1) .ne. Q_rp(i,j,k,ieq,1)) then
                    print *,'NaN. Bailing out...','  xc  =',xc(i),'  yc  =',yc(j),'  zc  =',zc(k),'  ieq  =',ieq
                    call exit(-1)
                endif
            end do

        end do
        end do
        end do

    end subroutine advance_time_level_gl

!----------------------------------------------------------------------------------------------


end program main
