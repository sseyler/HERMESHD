module initial_condition_mod

use parameters_mod
use auxiliary_mod

contains

    subroutine initial_condition(Q_r, id)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,id
        real wtev

        wtev = T_floor

        Q_r(:,:,:,:,:) = 0.0
        Q_r(:,:,:,rh,1) = rh_floor
        Q_r(:,:,:,en,1) = T_floor*rh_floor/(aindex - 1.)
        MMask(:,:,:) = .false.

        if ( id .eq. 1 ) then
            call fill_fluid(Q_r)
        end if
        if ( id .eq. 2 ) then
            call fill_fluid2(Q_r)
        end if
    end subroutine initial_condition

!-------------------------------------------------------

    subroutine fill_fluid(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,ieq
        real x,y,wtev,rnum,w0,jet_strength
        real qquad(npg,nQ),xcc,ycc,zcc  ! bfint(npg,nbasis),qquadv(npg)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        wtev = T_floor
        w0 = 0.3
        jet_strength = 0.5
        ! test problem is an unstable flow jet in x with velocity perturbations in y

        do i = 1,nx
        do j = 1,ny
        do k = 1,nz

            call random_number(rand_number)
            rnum = (rand_number - 0.5)

            !-------------------------------------------------------
            ! from "viscosity" version
            Q_r0(i,j,k,rh,1) = rh_floor
            Q_r0(i,j,k,en,1) = wtev*Q_r0(i,j,k,rh,1)/(aindex - 1.)

            xcc = xc(i)
            ycc = yc(j)
            zcc = zc(k)

            Q_r0(i,j,k,rh,1) = rh_fluid
            Q_r0(i,j,k,my,1) = 0.001*rnum/cosh(20*yc(j)/lyu)**2
            Q_r0(i,j,k,mz,1) = 0.
            Q_r0(i,j,k,mx,1) = jet_strength*Q_r0(i,j,k,rh,1)/cosh(20*ycc/lyu)/1.
            Q_r0(i,j,k,en,1) = wtev*Q_r0(i,j,k,rh,1)/(aindex - 1.)                  &
                             + 0.5*( Q_r0(i,j,k,mx,1)**2                            &
                                  +  Q_r0(i,j,k,my,1)**2                            &
                                  +  Q_r0(i,j,k,mz,1)**2 ) / Q_r0(i,j,k,rh,1)

        end do
        end do
        end do

    end subroutine fill_fluid


!-------------------------------------------------------

    subroutine fill_fluid2(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,igrid,ieq
        real x,y,wtev,rnum,w0
        real qquad(npg,nQ),xcc,ycc,zcc  ! bfint(npg,nbasis),qquadv(npg)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        wtev = T_floor
        w0 = 0.3

        ! test problem is an unstable flow jet in x with velocity perturbations in y

        do i = 1,nx
            do j = 1,ny
                do k = 1,nz

                    call random_number(rand_number)
                    rnum = (rand_number - 0.5)

                    qquad(:,:) = 0.

                    do igrid=1,npg

                        qquad(igrid,rh) = rh_floor
                        qquad(igrid,en) = wtev*qquad(igrid,rh)/(aindex - 1.)

                        xcc = xc(i) + bfvals_int(igrid,kx)*0.5/dxi
                        ycc = yc(j) + bfvals_int(igrid,ky)*0.5/dyi
                        zcc = zc(k) + bfvals_int(igrid,kz)*0.5/dzi

                        qquad(igrid,rh) = rh_fluid
                        qquad(igrid,my) = 0.001*rnum/cosh(20*yc(j)/lyu)**2 !0.001*rnum/1.
                        qquad(igrid,mz) = 0.
                        qquad(igrid,mx) = 1.0*qquad(igrid,rh)/cosh(20*ycc/lyu)/1.
                        qquad(igrid,en) = wtev*qquad(igrid,rh)/(aindex - 1.)        &
                                        + 0.5*(qquad(igrid,mx)**2                   &
                                             + qquad(igrid,my)**2                   &
                                             + qquad(igrid,mz)**2)/qquad(igrid,rh)

                    end do

                    do ieq=rh,en
                        do ir=1,nbasis
                            Q_r0(i,j,k,ieq,ir) = 0.125 * cbasis(ir) * sum(wgt3d(1:npg) &
                                                       * bfvals_int(1:npg,ir)*qquad(1:npg,ieq))
                        end do
                    end do

                end do
            end do
        end do

    end subroutine fill_fluid2


end module initial_condition_mod
