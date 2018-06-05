!***** SOURCES.F90 ***********************************************************************
module sources

use input
use params
use helpers
use spatial
use basis
use random

! use random  ! TODO: commented to get working w/o MKL

real, dimension(nx,ny,nz,nQ,nbasis) :: source_r

contains

    subroutine source_calc(Q_io, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, intent(inout) :: dt
        real, dimension(npg,nQ) :: source
        real, dimension(nQ) :: Qin
        real :: dn,dni, vx,vy,vz,vx2,vy2,vz2,vsqd3, pr,te, nudn,nudnd3
        real :: S_xx,S_xy,S_xz,S_yy,S_yz,S_zz
        integer :: i,j,k,ieq,ipg,ir

        source_r(:,:,:,:,:) = 0.0
        source(:,:) = 0.0

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            if (ivis == 'linear_src') then
                S_xx = Sflux(npts_llns,i,j,k,1,1)
                S_yy = Sflux(npts_llns,i,j,k,2,2)
                S_zz = Sflux(npts_llns,i,j,k,3,3)
                S_xy = Sflux(npts_llns,i,j,k,1,2)
                S_xz = Sflux(npts_llns,i,j,k,1,3)
                S_yz = Sflux(npts_llns,i,j,k,2,3)
            end if

            do ipg = 1,npg
                do ieq = 1,nQ
                    Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_io(i,j,k,ieq,1:nbasis))
                end do
                !--------------------------------
                source(ipg,rh) = 0
                source(ipg,mx) = 0
                source(ipg,my) = 0 !-Qin(rh)*gr  ! DEBUG
                source(ipg,mz) = 0
                source(ipg,en) = 0
                !--------------------------------
                select case (ivis)
                case('linear_ex')  ! NOTE: LINEARIZED 10-moment-eqns (explicit)
                    source(ipg,exx) = -nu*Qin(exx)
                    source(ipg,eyy) = -nu*Qin(eyy)
                    source(ipg,ezz) = -nu*Qin(ezz)
                    source(ipg,exy) = -nu*Qin(exy)
                    source(ipg,exz) = -nu*Qin(exz)
                    source(ipg,eyz) = -nu*Qin(eyz)
                case('linear')  ! NOTE: LINEARIZED 10-moment-eqns (IMEX)
                    source(ipg,exx) = 0 ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,eyy) = 0 ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,ezz) = 0 ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,exy) = 0 ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,exz) = 0 ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,eyz) = 0 ! impl solv w/ sources next step in adv_time_lvl
                case('linear_src')
                    source(ipg,exx) = -nu*S_xx ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,eyy) = -nu*S_yy ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,ezz) = -nu*S_zz ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,exy) = -nu*S_xy ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,exz) = -nu*S_xz ! impl solv w/ sources next step in adv_time_lvl
                    source(ipg,eyz) = -nu*S_yz ! impl solv w/ sources next step in adv_time_lvl
                case('full') ! NOTE: FULL NONLINEAR 10-moment-eqns
                    dn  = Qin(rh)
                    dni = 1./dn
                    vx  = Qin(mx)*dni
                    vy  = Qin(my)*dni
                    vz  = Qin(mz)*dni
                    vx2 = vx*vx
                    vy2 = vy*vy
                    vz2 = vz*vz
                    vsqd3 = c1d3*(vx2 + vy2 + vz2)
                    source(ipg,exx) = dn*(vx2 - vsqd3)
                    source(ipg,eyy) = dn*(vy2 - vsqd3)
                    source(ipg,ezz) = dn*(vz2 - vsqd3)
                    source(ipg,exy) = dn*vx*vy
                    source(ipg,exz) = dn*vx*vz
                    source(ipg,eyz) = dn*vy*vz
                end select
            end do
            !--------------------------------
            do ir=1,nbasis
                do ieq = 1,nQ
                    source_r(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(bval_int_wgt(1:npg,ir)*source(1:npg,ieq))
                end do
            end do
        end do
        end do
        end do

    end subroutine source_calc
    !---------------------------------------------------------------------------

end module sources
