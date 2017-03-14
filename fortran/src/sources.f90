!***** SOURCES.F90 ***********************************************************************
module sources

use parameters
use helpers
use basis_funcs

contains

    subroutine source_calc(Q_ri)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_ri

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
                    Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,k,ieq,1:nbasis))
                end do

                source(ipg,rh:en) = 0

                select case (ivis)
                case(0)  ! NOTE: FOR INVISCID FLOW, MUST ENSURE coll = colvis = 0
                        source(ipg,pxx) = -coll*Qin(pxx)
                        source(ipg,pyy) = -coll*Qin(pyy)
                        source(ipg,pzz) = -coll*Qin(pzz)
                        source(ipg,pxy) = -coll*Qin(pxy)
                        source(ipg,pxz) = -coll*Qin(pxz)
                        source(ipg,pyz) = -coll*Qin(pyz)
                    case(1)
                        source(ipg,pxx) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pyy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pzz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pxy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pxz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pyz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                    case(2) ! NOTE: for doing full 10-moment equations explicitly
                        dn = Qin(rh)
                        dni = 1./dn
                        vx = Qin(mx)*dni
                        vy = Qin(my)*dni
                        vz = Qin(mz)*dni
                        vx2 = vx**2
                        vy2 = vy**2
                        vz2 = vz**2
                        ! pr = (aindex - 1.)*(Qin(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
                        ! te = Pres*dni
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


end module sources
