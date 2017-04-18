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
                    source(ipg,exx) = -coll*Qin(exx)
                    source(ipg,eyy) = -coll*Qin(eyy)
                    source(ipg,ezz) = -coll*Qin(ezz)
                    source(ipg,exy) = -coll*Qin(exy)
                    source(ipg,exz) = -coll*Qin(exz)
                    source(ipg,eyz) = -coll*Qin(eyz)
                case(1)  ! NOTE: solving the linearized 10-moment equations
                    source(ipg,exx) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,eyy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,ezz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,exy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,exz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    source(ipg,eyz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                case(2) ! NOTE: solving the full (nonlinear) 10-moment equations
                    source(ipg,exx) = 0
                    source(ipg,eyy) = 0
                    source(ipg,ezz) = 0
                    source(ipg,exy) = 0
                    source(ipg,exz) = 0
                    source(ipg,eyz) = 0
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

                    source(ipg,exx) = colldn*(2*vx2 - vy2 - vz2)*c1d3
                    source(ipg,eyy) = colldn*(2*vy2 - vz2 - vx2)*c1d3
                    source(ipg,ezz) = colldn*(2*vz2 - vx2 - vy2)*c1d3
                    source(ipg,exy) = colldn*vx*vy
                    source(ipg,exz) = colldn*vx*vz
                    source(ipg,eyz) = colldn*vy*vz
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

end module sources
