!***** INTEGRATOR.F90 ********************************************************************
module integrator

    use input!, only : nx,ny,nz
    use params!, only : nQ,nbasis,Q_r0,Q_r1,Q_r2,Q_r3
    use helpers

    use prepare_step
    use sources
    use flux

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutine for temporal integration
    !-----------------------------------------------------------------
    abstract interface
        subroutine update_ptr (Q_io, Q_1, Q_2, dt)
            use input, only : nx,ny,nz
            use params, only : nQ,nbasis

            real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
            real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2
            real, intent(inout) :: dt
        end subroutine update_ptr
    end interface
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Initialize pointer to temporal integration subroutine
    !-----------------------------------------------------------------
    procedure (update_ptr), pointer :: update => null ()
    !---------------------------------------------------------------------------

contains

    !===========================================================================
    ! Select user-specified integration method at runtime
    !-----------------------------------------------------------------
    subroutine select_integrator(name, integrator)
        implicit none
        character(*), intent(in) :: name
        procedure(update_ptr), pointer :: integrator

        select case (name)
            case ('heun')
                call mpi_print(iam, 'Selected 2nd-order Runga-Kutta (Heun) integration')
                integrator => RK2
            case ('shu-osher')
                call mpi_print(iam, 'Selected 3rd-order Runga-Kutta (Shu-Osher) integration')
                integrator => RK3
            case default
                call mpi_print(iam, 'Defaulting to 2nd-order Runga-Kutta (Heun) integration')
                integrator => RK2
        end select
    end subroutine
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Temporal integration subroutines (subject to change!)
    !-----------------------------------------------------------------
    subroutine RK2(Q_io, Q_1, Q_2, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2
        real, intent(inout) :: dt

        call euler_step(Q_io, Q_1, dt)
        call euler_step(Q_1, Q_2, dt)
        Q_io = 0.5 * ( Q_io + Q_2 )
    end subroutine RK2

    !----------------------------------------------------
    subroutine RK3(Q_io, Q_1, Q_2, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2
        real, intent(inout) :: dt

        call euler_step(Q_io, Q_1, dt)
        call euler_step(Q_1, Q_2, dt)
        Q_1 = 0.75*Q_io + 0.25*Q_2

        call euler_step(Q_1, Q_2, dt)  ! re-use the second array
        Q_io = c1d3*Q_io + c2d3*Q_2
    end subroutine RK3
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Explicit Euler integration step
    !-----------------------------------------------------------------
    subroutine euler_step(Q_in, Q_out, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real, intent(inout) :: dt

        call prep_advance(Q_in)
        call calc_rhs(Q_in)
        call advance_time_level(Q_in, Q_out, dt)
    end subroutine euler_step


    !----------------------------------------------------
    subroutine prep_advance(Q_io)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        if (ieos == 1 .or. ieos == 2) call limiter(Q_io)
        call exchange_flux(Q_io)
        call apply_boundaries
    end subroutine prep_advance

    !----------------------------------------------------
    subroutine calc_rhs(Q_io)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        call glflux(Q_io)
        call source_calc(Q_io)
    end subroutine calc_rhs

    !----------------------------------------------------
    subroutine advance_time_level(Q_in, Q_out, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real, intent(inout) :: dt

        real P_xx,P_yy,P_zz,P_xy,P_xz,P_yz
        real glf_pxx,glf_pyy,glf_pzz,glf_pxy,glf_pxz,glf_pyz
        real src_pxx,src_pyy,src_pzz,src_pxy,src_pxz,src_pyz
        real Q_sum,glf_sum,src_sum,faci
        integer i,j,k,ieq,ir

        faci = 1./(1. + dt*coll)  ! need to ensure a value is given to coll!

        select case(ivis)
        !-----------------------------------------------------------------------
        case(0) ! Fully explicit integration: forward-Euler for all fields
            do ir=1,nbasis
            do ieq = 1,nQ
                do k = 1,nz
                do j = 1,ny
                do i = 1,nx
                    Q_out(i,j,k,ieq,ir) =                                       &
                        Q_in(i,j,k,ieq,ir) - dt*( glflux_r(i,j,k,ieq,ir)        &
                                                - source_r(i,j,k,ieq,ir) )
                end do
                end do
                end do
            end do
            end do

        !-----------------------------------------------------------------------
        case(1) ! Semi-implicit integration: backward-Euler for stress fields
            do ir=1,nbasis
                do k = 1,nz
                do j = 1,ny
                do i = 1,nx
                    !--- Separate ----------
                    ! Q_out(i,j,k,rh:en,ir) =     Q_in(i,j,k,rh:en,ir)            &
                    !                     - ( glflux_r(i,j,k,rh:en,ir)            &
                    !                       - source_r(i,j,k,rh:en,ir) ) * dt
                    ! Q_out(i,j,k,pxx:pyz,ir) = ( Q_in(i,j,k,pxx:pyz,ir)          &
                    !                    - dt*glflux_r(i,j,k,pxx:pyz,ir)          &
                    !                    + dt*source_r(i,j,k,pxx:pyz,ir) ) * faci

                    !--- Combo -------------
                    Q_out(i,j,k,1:nQ,ir) =                                      &
                        Q_in(i,j,k,1:nQ,ir) - dt*( glflux_r(i,j,k,1:nQ,ir)      &
                                                 - source_r(i,j,k,1:nQ,ir) )
                    Q_out(i,j,k,pxx:pyz,ir) = faci * Q_out(i,j,k,pxx:pyz,ir)
                end do
                end do
                end do
            end do

        !-----------------------------------------------------------------------
        case(2) ! Semi-implicit integration: backward-Euler for stress fields
            do ir=1,nbasis
                do k = 1,nz
                do j = 1,ny
                do i = 1,nx
                    P_xx = Q_in(i,j,k,pxx,ir)
                    P_yy = Q_in(i,j,k,pyy,ir)
                    P_zz = Q_in(i,j,k,pzz,ir)
                    P_xy = Q_in(i,j,k,pxy,ir)
                    P_xz = Q_in(i,j,k,pxz,ir)
                    P_yz = Q_in(i,j,k,pyz,ir)

                    glf_pxx = glflux_r(i,j,k,pxx,ir)
                    glf_pyy = glflux_r(i,j,k,pyy,ir)
                    glf_pzz = glflux_r(i,j,k,pzz,ir)
                    glf_pxy = glflux_r(i,j,k,pxy,ir)
                    glf_pxz = glflux_r(i,j,k,pxz,ir)
                    glf_pyz = glflux_r(i,j,k,pyz,ir)

                    src_pxx = source_r(i,j,k,pxx,ir)
                    src_pyy = source_r(i,j,k,pyy,ir)
                    src_pzz = source_r(i,j,k,pzz,ir)
                    src_pxy = source_r(i,j,k,pxy,ir)
                    src_pxz = source_r(i,j,k,pxz,ir)
                    src_pyz = source_r(i,j,k,pyz,ir)

                    Q_sum   = coll*dt*   ( P_xx    + P_yy    + P_zz    )*c1d3
                    glf_sum = coll*dt**2*( glf_pxx + glf_pyy + glf_pzz )*c1d3
                    src_sum = coll*dt**2*( src_pxx + src_pyy + src_pzz )*c1d3

                    Q_out(i,j,k,pxx,ir) = faci * ( P_xx - dt*glf_pxx + dt*src_pxx &
                                                 + Q_sum -   glf_sum +    src_sum )

                    Q_out(i,j,k,pyy,ir) = faci * ( P_yy - dt*glf_pyy + dt*src_pyy &
                                                 + Q_sum -   glf_sum +    src_sum )

                    Q_out(i,j,k,pzz,ir) = faci * ( P_zz - dt*glf_pzz + dt*src_pzz &
                                                 + Q_sum -   glf_sum +    src_sum )

                    do ieq = pxy,nQ
                        Q_out(i,j,k,ieq,ir) = ( Q_in(i,j,k,ieq,ir)              &
                                           - dt*glflux_r(i,j,k,ieq,ir)          &
                                           + dt*source_r(i,j,k,ieq,ir) ) * faci
                    end do
                end do
                end do
                end do
            end do

        end select

        !-----------------------------------------------------------------------
        ! Check for NaNs; bail out with info.
        do ieq = 1,nQ
        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          if ( Q_out(i,j,k,ieq,1) /= Q_out(i,j,k,ieq,1) ) then
            print *,'------------------------------------------------'
            print *,'NaN. Bailing out...'
            write(*,'(A7,I9,A7,I9,A7,I9)')          '   i = ',   i, '   j = ',    j, '   k = ',k
            write(*,'(A7,ES9.2,A7,ES9.2,A7,ES9.2)') '  xc = ',xc(i),'  yc = ',yc(j), '  zc = ', zc(k)
            write(*,'(A14,I2,A7,I2)') '    >>> iam = ', iam, ' ieq = ', ieq
            print *,''
            call exit(-1)
          endif
        end do
        end do
        end do
        end do

    end subroutine advance_time_level
    !---------------------------------------------------------------------------

end module integrator
