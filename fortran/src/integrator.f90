module integrator

    use input!, only : nx,ny,nz
    use parameters!, only : nQ,nbasis,Q_r0,Q_r1,Q_r2,Q_r3
    use helpers

    use prepare_step
    use sources
    use flux

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutine for temporal integration
    !-----------------------------------------------------------------
    abstract interface
        subroutine update_ptr (Q_io, Q_1, Q_2)
            use input, only : nx,ny,nz
            use parameters, only : nQ,nbasis

            real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
            real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2
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
    subroutine RK2(Q_io, Q_1, Q_2)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2

        call euler_step(Q_io, Q_1)
        call euler_step(Q_1, Q_2)
        Q_io = 0.5 * ( Q_io + Q_2 )
    end subroutine RK2

    !----------------------------------------------------
    subroutine RK3(Q_io, Q_1, Q_2)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2

        call euler_step(Q_io, Q_1)
        call euler_step(Q_1, Q_2)
        Q_1 = 0.75*Q_io + 0.25*Q_2

        call euler_step(Q_1, Q_2)  ! re-use the second array
        Q_io = c1d3 * ( Q_io + 2.0*Q_2 )
    end subroutine RK3
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Explicit Euler integration step
    !-----------------------------------------------------------------
    subroutine euler_step(Q_in, Q_out)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out

        call prep_advance(Q_in)
        call calc_rhs(Q_in)
        call advance_time_level(Q_in, Q_out)
    end subroutine euler_step

    !----------------------------------------------------
    subroutine prep_advance(Q_io)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        if (ieos .eq. 1 .or. ieos .eq. 2) call limiter(Q_io)
        call exchange_flux(Q_io)
        call apply_boundaries
    end subroutine prep_advance

    !----------------------------------------------------
    subroutine calc_rhs(Q_io)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        ! call flux_calc(Q_io)
        ! call innerintegral(Q_io)
        call glflux(Q_io)
        call source_calc(Q_io)
    end subroutine calc_rhs

    !----------------------------------------------------
    subroutine advance_time_level(Q_in, Q_out)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real faci
        integer i,j,k,ieq,ir

        ! faci = 1./(1. + dt*coll)  ! need to ensure a value is given to coll!

        do ir=1,nbasis
        do ieq = 1,nQ
            do k = 1,nz
            do j = 1,ny
            do i = 1,nx
                Q_out(i,j,k,ieq,ir) = Q_in(i,j,k,ieq,ir)                        &
                                    - dt*( glflux_r(i,j,k,ieq,ir)               &
                                         - source_r(i,j,k,ieq,ir) )
            end do
            end do
            end do
        end do
        end do

        ! Implicit solve of stress variables?
        ! do k = 1,nz
        ! do j = 1,ny
        ! do i = 1,nx
        !     do ir=1,nbasis
        !     do ieq = pxx,nQ
        !         Q_out(i,j,k,ieq,ir) = ( Q_in(i,j,k,ieq,ir)                      &
        !                              - dt*glflux_r(i,j,k,ieq,ir)                &
        !                              + dt*source_r(i,j,k,ieq,ir) ) * faci
        !     end do
        !     end do
        ! end do
        ! end do
        ! end do

        do ieq = 1,nQ
            do k = 1,nz
            do j = 1,ny
            do i = 1,nx
                if ( Q_out(i,j,k,ieq,1) /= Q_out(i,j,k,ieq,1)) then
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
