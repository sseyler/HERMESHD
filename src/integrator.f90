!***** INTEGRATOR.F90 ********************************************************************
module integrator

use input!, only : nx,ny,nz
use params!, only : nQ,nbasis,Q_r0,Q_r1,Q_r2,Q_r3
use helpers
use spatial
! use timestep

use prep_step
use sources
use flux

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutine for temporal integration
    !-----------------------------------------------------------------
    abstract interface
        subroutine step_ptr(Q_io, Q_1, Q_2, dt)
            use input, only : nx,ny,nz
            use params, only : nQ,nbasis

            real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
            real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_1, Q_2
            real, intent(inout) :: dt
        end subroutine step_ptr
    end interface
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Initialize pointer to temporal integration subroutine
    !-----------------------------------------------------------------
    procedure (step_ptr), pointer :: step => null ()
    !---------------------------------------------------------------------------

contains

    !===========================================================================
    ! Select user-specified integration method at runtime
    !-----------------------------------------------------------------
    subroutine select_integrator(name, integrator)
        implicit none
        character(*), intent(in) :: name
        procedure(step_ptr), pointer :: integrator

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
        call calc_rhs(Q_in, dt)
        call adv_time_lvl(Q_in, Q_out, dt)
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
    subroutine calc_rhs(Q_io, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, intent(inout) :: dt

        call glflux(Q_io, dt)
        call source_calc(Q_io, dt)
    end subroutine calc_rhs


    !---------------------------------------------------------------------------
    subroutine adv_time_lvl(Q_in, Q_out, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in)  :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real, intent(inout) :: dt
        real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
        integer :: i,j,k,ieq,ir
        real :: fac

        select case(ivis)
        case(0)  ! LINEARNIZED (explicit)
            fac = 1.
        case(1)  ! LINEARNIZED (IMEX)
            fac = 1./(1. + nu*dt)
        case(2)  ! NONLINEAR (IMEX)
            fac = 1
        end select

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          !--------------------------------
          do ieq = 1,en
            do ir=1,nbasis
                ! NOTE: from GFV1
                ! Q_out(i,j,k,ieq,ir) = Q_in(i,j,k,ieq,ir) - dt*(glflux_r(i,j,k,ieq,ir) - source_r(i,j,k,ieq,ir))
                Q_out(i,j,k,ieq,ir) = Q_in(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir)
            end do
          end do
          !--------------------------------
          do ieq = exx,nQ
            do ir=1,nbasis
                ! NOTE: from GFV1
                ! Q_out(i,j,k,ieq,ir) = (Q_in(i,j,k,ieq,ir) - dt*(glflux_r(i,j,k,ieq,ir) - source_r(i,j,k,ieq,ir)))*fac
                Q_out(i,j,k,ieq,ir) = (Q_in(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir))*fac
            end do
          end do
          !--------------------------------
        end do
        end do
        end do

        if (ivis == 2) then ! NOTE: only for FULL NONLINEAR IMEX code
            call adv_time_lvl_src_v1(Q_out, Q_out, dt)
        end if
    end subroutine adv_time_lvl
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    subroutine adv_time_lvl_src_v1(Q_in, Q_out, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in)  :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real, intent(inout) :: dt
        real, dimension(npg,nQ) :: source, Qout
        real, dimension(nQ) :: Qin
        integer i,j,k,ieq,ipg,ir
        real dn,dni,vx,vy,vz,P, vx2,vy2,vz2,vsq,fac

        fac = 1./(1. + nu*dt)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
          do ir=1,nbasis
            Q_out(i,j,k,exx,ir) = ( Q_in(i,j,k,exx,ir) + dt*source_r(i,j,k,exx,ir) ) * fac
            Q_out(i,j,k,eyy,ir) = ( Q_in(i,j,k,eyy,ir) + dt*source_r(i,j,k,eyy,ir) ) * fac
            Q_out(i,j,k,ezz,ir) = ( Q_in(i,j,k,ezz,ir) + dt*source_r(i,j,k,ezz,ir) ) * fac
            Q_out(i,j,k,exy,ir) = ( Q_in(i,j,k,exy,ir) + dt*source_r(i,j,k,exy,ir) ) * fac
            Q_out(i,j,k,exz,ir) = ( Q_in(i,j,k,exz,ir) + dt*source_r(i,j,k,exz,ir) ) * fac
            Q_out(i,j,k,eyz,ir) = ( Q_in(i,j,k,eyz,ir) + dt*source_r(i,j,k,eyz,ir) ) * fac
          end do
        end do
        end do
        end do
    end subroutine adv_time_lvl_src_v1
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! NOTE: may be able to make functions from Qin loop and final loop w/
    !       0.125*cbasis(...) for re-use by adv_time_lvl_src & source_calc
    subroutine adv_time_lvl_src_v0(Q_in, Q_out, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in)  :: Q_in
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_out
        real, intent(inout) :: dt

        real, dimension(npg,nQ) :: source, Qout
        real, dimension(nQ) :: Qin
        integer i,j,k,ieq,ipg,ir
        real dn,dni,vx,vy,vz,P, vx2,vy2,vz2,vsq,nudt,fac

        nudt = nu*dt
        fac = 1./(1. + nudt)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            do ipg = 1,npg
                do ieq = 1,nQ
                    Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_in(i,j,k,ieq,1:nbasis))
                end do
                !--------------------------------
                dn  = Qin(rh)
                dni = 1./Qin(rh)
                vx  = Qin(mx)*dni
                vy  = Qin(my)*dni
                vz  = Qin(mz)*dni
                vx2 = vx*vx
                vy2 = vy*vy
                vz2 = vz*vz
                vsq = vx2 + vy2 + vz2
                P   = aindm1*(Qin(en) - 0.5*dn*vsq)
                if (P < P_floor) P = P_floor  ! NOTE: is this necessary?
                !--------------------------------
                Qout(ipg,exx) = ( Qin(exx) + nudt*(P + dn*vx2  ) ) * fac
                Qout(ipg,eyy) = ( Qin(eyy) + nudt*(P + dn*vy2  ) ) * fac
                Qout(ipg,ezz) = ( Qin(ezz) + nudt*(P + dn*vz2  ) ) * fac
                Qout(ipg,exy) = ( Qin(exy) + nudt*(    dn*vx*vy) ) * fac
                Qout(ipg,exz) = ( Qin(exz) + nudt*(    dn*vx*vz) ) * fac
                Qout(ipg,eyz) = ( Qin(eyz) + nudt*(    dn*vy*vz) ) * fac
            end do

            do ieq=exx,nQ
                do ir=1,nbasis
                    Q_out(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(wgt3d(1:npg)*bfvals_int(1:npg,ir)*Qout(1:npg,ieq))
                end do
            end do
        end do
        end do
        end do
    end subroutine adv_time_lvl_src_v0
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
