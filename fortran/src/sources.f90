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
        real dn,dni,vx,vy,vz,en_floor
        integer i,j,k,ieq,ipg,ir

        source_r(:,:,:,:,:) = 0.0
        source(:,:) = 0.0
        en_floor = P_floor/(aindex - 1.)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

            do ipg = 1,npg

                ! NOTE: Needs to be done for each var/eqn with nonzero sources,
                ! i.e., when Qin() is used for something
                ! do ieq = 1,nQ
                !     Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,k,ieq,1:nbasis))
                ! end do
                source(ipg,rh:en) = 0

                select case (ivis)
                case (0)  ! NOTE: THIS DOESN'T MAKE ANY SENSE IF THERE'S NO VISCOSITY
                        do ieq = pxx,pyz
                            Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,k,ieq,1:nbasis))
                        end do
                        source(ipg,pxx) = -coll*Qin(pxx)
                        source(ipg,pyy) = -coll*Qin(pyy)
                        source(ipg,pzz) = -coll*Qin(pzz)
                        source(ipg,pxy) = -coll*Qin(pxy)
                        source(ipg,pxz) = -coll*Qin(pxz)
                        source(ipg,pyz) = -coll*Qin(pyz)
                    case (1)
                        source(ipg,pxx) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pyy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pzz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pxy) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pxz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                        source(ipg,pyz) = 0  ! impl. solv w/ sources at next tstep in advance_time_level_gl
                end select

            end do

            do ir=1,nbasis
                do ieq = 1,nQ
                    source_r(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(bval_int_wgt(1:npg,ir)*source(1:npg,ieq))
                end do
            end do

            ! Copied from 10-moment code... not sure what it's for, yet.
            ! source_r(i,j,k,rh:en,kxx:kyy) = -1.e-1*Q_ri(i,j,k,rh:en,kxx:kyy)

        end do
        end do
        end do

    end subroutine source_calc

!----------------------------------------------------------------------------------------------


end module sources
