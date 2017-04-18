!***** INTEGRATOR.F90 ********************************************************************
module integrator

use input!, only : nx,ny,nz
use params!, only : nQ,nbasis,Q_r0,Q_r1,Q_r2,Q_r3
use helpers
use spatial
! use timestep

use prepare_step
use sources
use flux

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutine for temporal integration
    !-----------------------------------------------------------------
    abstract interface
        subroutine update_ptr(Q_io, Q_1, Q_2, dt)
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
    end subroutine select_integrator
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
        call check_for_NaNs(Q_out)
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


    !---------------------------------------------------------------------------
    subroutine advance_time_level(Q_in, Q_out, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real, intent(inout) :: dt

        integer i,j,k,ieq,ir
        real fac

        fac = 1.
        if (ivis == 1) fac = 1./(1. + coll*dt)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

            do ieq = 1,en
              do ir=1,nbasis
                Q_out(i,j,k,ieq,ir) = Q_in(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir)
              end do
            end do

            do ieq = exx,nQ
              do ir=1,nbasis
                Q_out(i,j,k,ieq,ir) = (Q_in(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir))*fac
              end do
            end do

        end do
        end do
        end do

        if (ivis == 2) call impl_advance_source(Q_out, dt)

    end subroutine advance_time_level
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! NOTE: may be able to make functions from Qin loop and final loop w/
    !       0.125*cbasis(...) for re-use by impl_advance_source & source_calc
    subroutine impl_advance_source(Q_io, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, intent(inout) :: dt

        real, dimension(npg,nQ) :: source, Qout
        real, dimension(nQ) :: Qin
        integer i,j,k,ieq,ipg,ir
        real dn,dni,vx,vy,vz,P, colldt,fac

        colldt = coll*dt
        fac = 1./(1. + colldt)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            do ipg = 1,npg
                do ieq = 1,nQ
                    Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_io(i,j,k,ieq,1:nbasis))
                end do

                dn = Qin(rh)
                dni = 1./Qin(rh)
                vx = Qin(mx)*dni
                vy = Qin(my)*dni
                vz = Qin(mz)*dni
                P  = aindm1*(Qin(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))

                Qout(ipg,exx) = ( Qin(exx) + colldt*(P + dn*vx**2) ) * fac
                Qout(ipg,eyy) = ( Qin(eyy) + colldt*(P + dn*vy**2) ) * fac
                Qout(ipg,ezz) = ( Qin(ezz) + colldt*(P + dn*vz**2) ) * fac
                Qout(ipg,exy) = ( Qin(exy) + colldt*(    dn*vx*vy) ) * fac
                Qout(ipg,exz) = ( Qin(exz) + colldt*(    dn*vx*vz) ) * fac
                Qout(ipg,eyz) = ( Qin(eyz) + colldt*(    dn*vy*vz) ) * fac
            end do

            do ieq=exx,nQ
                do ir=1,nbasis
                    Q_io(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(wgt3d(1:npg)*bfvals_int(1:npg,ir)*Qout(1:npg,ieq))
                end do
            end do
        end do
        end do
        end do

    end subroutine impl_advance_source
    !---------------------------------------------------------------------------


    !----------------------------------------------------
    subroutine advance_time_level_v0(Q_in, Q_out, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real, intent(inout) :: dt

        real E_xx,E_yy,E_zz,E_xy,E_xz,E_yz
        real glf_exx,glf_eyy,glf_ezz,glf_exy,glf_exz,glf_eyz
        real src_exx,src_eyy,src_ezz,src_exy,src_exz,src_eyz
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
                    Q_out(i,j,k,exx:eyz,ir) = faci * Q_out(i,j,k,exx:eyz,ir)
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
                    E_xx = Q_in(i,j,k,exx,ir)
                    E_yy = Q_in(i,j,k,eyy,ir)
                    E_zz = Q_in(i,j,k,ezz,ir)
                    E_xy = Q_in(i,j,k,exy,ir)
                    E_xz = Q_in(i,j,k,exz,ir)
                    E_yz = Q_in(i,j,k,eyz,ir)

                    glf_exx = glflux_r(i,j,k,exx,ir)
                    glf_eyy = glflux_r(i,j,k,eyy,ir)
                    glf_ezz = glflux_r(i,j,k,ezz,ir)
                    glf_exy = glflux_r(i,j,k,exy,ir)
                    glf_exz = glflux_r(i,j,k,exz,ir)
                    glf_eyz = glflux_r(i,j,k,eyz,ir)

                    src_exx = source_r(i,j,k,exx,ir)
                    src_eyy = source_r(i,j,k,eyy,ir)
                    src_ezz = source_r(i,j,k,ezz,ir)
                    src_exy = source_r(i,j,k,exy,ir)
                    src_exz = source_r(i,j,k,exz,ir)
                    src_eyz = source_r(i,j,k,eyz,ir)

                    Q_sum   = coll*dt*   ( E_xx    + E_yy    + E_zz    )*c1d3
                    glf_sum = coll*dt**2*( glf_exx + glf_eyy + glf_ezz )*c1d3
                    src_sum = coll*dt**2*( src_exx + src_eyy + src_ezz )*c1d3

                    Q_out(i,j,k,exx,ir) = faci * ( E_xx - dt*glf_exx + dt*src_exx &
                                                 + Q_sum -   glf_sum +    src_sum )

                    Q_out(i,j,k,eyy,ir) = faci * ( E_yy - dt*glf_eyy + dt*src_eyy &
                                                 + Q_sum -   glf_sum +    src_sum )

                    Q_out(i,j,k,ezz,ir) = faci * ( E_zz - dt*glf_ezz + dt*src_ezz &
                                                 + Q_sum -   glf_sum +    src_sum )

                    do ieq = exy,nQ
                        Q_out(i,j,k,ieq,ir) = ( Q_in(i,j,k,ieq,ir)              &
                                           - dt*glflux_r(i,j,k,ieq,ir)          &
                                           + dt*source_r(i,j,k,ieq,ir) ) * faci
                    end do
                end do
                end do
                end do
            end do

        end select

    end subroutine advance_time_level_v0
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! Check for NaNs; bail out with info.
    !---------------------------------------------------------------------------
    subroutine check_for_NaNs(Q_io)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        integer i,j,k,ieq

        do ieq = 1,nQ
        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            if ( Q_io(i,j,k,ieq,1) /= Q_io(i,j,k,ieq,1) ) then
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
    end subroutine check_for_NaNs
    !---------------------------------------------------------------------------

end module integrator
