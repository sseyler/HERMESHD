!***** SOURCES.F90 ***********************************************************************
module sources

use params
use helpers
use spatial
use basis

! use random  ! TODO: commented to get working w/o MKL

real, dimension(nx,ny,nz,nQ,nbasis) :: source_r

contains

    subroutine source_calc(Q_io)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        real, dimension(npg,nQ) :: source
        real, dimension(nQ) :: Qin
        real dn,dni, vx,vy,vz,vx2,vy2,vz2, pr,te, colldn
        integer i,j,k,ieq,ipg,ir

        source_r(:,:,:,:,:) = 0.0
        source(:,:) = 0.0

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

            do ipg = 1,npg

                do ieq = 1,nQ
                    Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_io(i,j,k,ieq,1:nbasis))
                end do

                source(ipg,rh:en) = 0

                select case (ivis)
                case(-1)  ! NOTE: explicitly solving the linearized 10-moment equations
                    source(ipg,pxx) = -coll*Qin(pxx)
                    source(ipg,pyy) = -coll*Qin(pyy)
                    source(ipg,pzz) = -coll*Qin(pzz)
                    source(ipg,pxy) = -coll*Qin(pxy)
                    source(ipg,pxz) = -coll*Qin(pxz)
                    source(ipg,pyz) = -coll*Qin(pyz)
                case(1)  ! NOTE: solving the linearized 10-moment equations
                    source(ipg,pxx) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,pyy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,pzz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,pxy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,pxz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,pyz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                case(2) ! NOTE: solving the full (nonlinear) 10-moment equations
                    source(ipg,pxx) = 0
                    source(ipg,pyy) = 0
                    source(ipg,pzz) = 0
                    source(ipg,pxy) = 0
                    source(ipg,pxz) = 0
                    source(ipg,pyz) = 0
                case(3) ! NOTE: solving the full (nonlinear) 10-moment equations the OLD way
                    dn = Qin(rh)
                    dni = 1./dn
                    vx = Qin(mx)*dni
                    vy = Qin(my)*dni
                    vz = Qin(mz)*dni
                    vx2 = vx**2
                    vy2 = vy**2
                    vz2 = vz**2
                    colldn = coll*dn

                    source(ipg,pxx) = colldn*(2*vx2 - vy2 - vz2)*c1d3
                    source(ipg,pyy) = colldn*(2*vy2 - vz2 - vx2)*c1d3
                    source(ipg,pzz) = colldn*(2*vz2 - vx2 - vy2)*c1d3
                    source(ipg,pxy) = colldn*vx*vy
                    source(ipg,pxz) = colldn*vx*vz
                    source(ipg,pyz) = colldn*vy*vz
                end select

            end do

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


    subroutine impl_source_calc(Q_i, Q_o)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_i
        real, dimension(nx,ny,nz,nQ,nbasis), intent(out) :: Q_o
        real, dimension(npg,nQ) :: source, Qout
        real, dimension(nQ) :: Qin
        integer i,j,k,ieq,ipg,ir
        real dn,dni,vx,vy,vz,P
        real colldt,fac

        ! if (ivis == 2) call impl_source_calc(Q_r1,Q_r1)
        colldt = coll*dt
        fac = 1./(1. + colldt)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            do ipg = 1,npg
                do ieq = 1,nQ
                    Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_i(i,j,k,ieq,1:nbasis))
                end do

                dn = Qin(rh)
                dni = 1./Qin(rh)
                vx = Qin(mx)*dni
                vy = Qin(my)*dni
                vz = Qin(mz)*dni
                P  = aindm1*(Qin(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))

                Qout(ipg,pxx) = ( Qin(exx) + colldt*(P + dn*vx**2) ) * fac
                Qout(ipg,pyy) = ( Qin(eyy) + colldt*(P + dn*vy**2) ) * fac
                Qout(ipg,pzz) = ( Qin(ezz) + colldt*(P + dn*vz**2) ) * fac
                Qout(ipg,pxy) = ( Qin(exy) + colldt*(    dn*vx*vy) ) * fac
                Qout(ipg,pxz) = ( Qin(exz) + colldt*(    dn*vx*vz) ) * fac
                Qout(ipg,pyz) = ( Qin(eyz) + colldt*(    dn*vy*vz) ) * fac
            end do

            do ieq=pxx,nQ
                do ir=1,nbasis
                    Q_o(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(wgt3d(1:npg)*bfvals_int(1:npg,ir)*Qout(1:npg,ieq))
                end do
            end do
        end do
        end do
        end do

    end subroutine impl_source_calc

end module sources
