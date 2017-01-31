module sources

use parameters
use helpers
use boundaries


contains

    subroutine source_calc(Q_ri)

        implicit none
        integer i,j,k,ieq,ipg,ir
        real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri
        real source(npg,nQ),Qin(nQ),dn,dni,vx,vy,vz,en_floor

        source_r(:,:,:,:,:) = 0.0
        source(:,:) = 0.0
        en_floor = P_floor/(aindex - 1.)

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

            do ipg = 1,npg
                do ieq = pxx,nQ
                    Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,k,ieq,1:nbasis))
                end do
                ! Sources for the fluid variables. Nonzero if randomly forced.
                source(ipg,rh:en) = 0.0 !amplen*(ran(iseed) - 0.5)

                ! Sources for the viscous stress tensor. The viscous stress
                ! is computed using hyperbolic relaxation to a parabolic problem:
                !     partial_t u = partial_x v , eps*partial_t v - v = D*partial_x u
                ! where eps is a small parameter and D is a diffusion coefficient.
                ! In the relaxation limit eps -> 0, this becomes:
                !     partial_t u = D*partial_x^2 u
                ! The following are the sources for this problem where epsi = 1/eps
                source(ipg,pxx) = -epsi*Qin(pxx) !+ amplv*(ran(iseed) - 0.5)
                source(ipg,pyy) = -epsi*Qin(pyy) !+ amplv*(ran(iseed) - 0.5)
                source(ipg,pzz) = -epsi*Qin(pzz) !+ amplv*(ran(iseed) - 0.5)
                source(ipg,pxy) = -epsi*Qin(pxy) !+ amplv*(ran(iseed) - 0.5)
                source(ipg,pxz) = -epsi*Qin(pxz) !+ amplv*(ran(iseed) - 0.5)
                source(ipg,pyz) = -epsi*Qin(pyz) !+ amplv*(ran(iseed) - 0.5)
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

!----------------------------------------------------------------------------------------------


end module sources
