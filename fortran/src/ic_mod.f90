module ic_mod

use parameters_mod

contains

    subroutine set_ic(Q_r, id)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,id
        real wtev

        wtev = T_floor

        Q_r(:,:,:,:,:) = 0.0
        Q_r(:,:,:,rh,1) = rh_floor
        Q_r(:,:,:,en,1) = T_floor*rh_floor/(aindex - 1.)
        MMask(:,:,:) = .false.

        if ( id .eq. 0 ) then
            call fill_fluid2(Q_r)
        end if
        if ( id .eq. 1 ) then
            call hydro_jet(Q_r)
        end if
        if ( id .eq. 2 ) then
            call isentropic_vortex(Q_r)
        end if

    end subroutine set_ic

!-------------------------------------------------------

    subroutine hydro_jet(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,ieq
        real x,y,wtev,rnum,jet_strength,rand_num
        real qquad(npg,nQ),xcc,ycc,zcc  ! bfint(npg,nbasis),qquadv(npg)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        wtev = T_floor
        jet_strength = 1.0

        do i = 1,nx
        do j = 1,ny
        do k = 1,nz

            call random_number(rand_num)
            rnum = (rand_num - 0.5)

            Q_r(i,j,k,rh,1) = rh_floor
            Q_r(i,j,k,en,1) = wtev*Q_r(i,j,k,rh,1)/(aindex - 1.)

            xcc = xc(i)
            ycc = yc(j)
            zcc = zc(k)

            Q_r(i,j,k,rh,1) = rh_fluid
            Q_r(i,j,k,my,1) = 0.005*rnum/cosh(20*ycc/lyu)**2
            Q_r(i,j,k,mz,1) = 0.
            Q_r(i,j,k,mx,1) = jet_strength*Q_r(i,j,k,rh,1)/cosh(20*ycc/lyu)/1.
            Q_r(i,j,k,en,1) = wtev*Q_r(i,j,k,rh,1)/(aindex - 1.)                &
                             + 0.5*( Q_r(i,j,k,mx,1)**2                         &
                                  +  Q_r(i,j,k,my,1)**2                         &
                                  +  Q_r(i,j,k,mz,1)**2 ) / Q_r(i,j,k,rh,1)

        end do
        end do
        end do

    end subroutine hydro_jet

    !-------------------------------------------------------

    subroutine isentropic_vortex(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k
        real rh_amb,vx_amb,vy_amb,vz_amb,p_amb,T_amb,beta
        real xctr,yctr,zctr,xp,yp,r2,delta_vx,delta_vy,delta_T
        real dn,vx,vy,vz,temp

        beta = 5.0             ! vortex strength
        rh_amb = rh_fluid      ! ambient density
        vx_amb = 1.0           ! ambient x-velocity
        vy_amb = 0.0           ! ambient y-velocity
        vz_amb = 0.0           ! ambient z-velocity
        T_amb  = 1.0           ! ambient temperature
        p_amb  = T_amb*rh_amb  ! ambient pressure

        xctr = 0               ! vortex center in x-direction
        yctr = 0               ! vortex center in y-direction
        zctr = 0               ! vortex center in z-direction

        do i = 1,nx
        do j = 1,ny
        do k = 1,nz

            xp = xc(i) - xctr   ! x-value from vortex center
            yp = yc(j) - yctr   ! y-value from vortex center
            r2 = xp**2 + yp**2  ! radial distance from vortex center

            delta_vx = -yp*beta/(2*pi) * exp( (1 - r2)/2 )
            delta_vy =  xp*beta/(2*pi) * exp( (1 - r2)/2 )
            delta_T  = -aindm1*beta**2/(8*aindex*pi**2) * exp(1 - r2)

            vx = vx_amb + delta_vx
            vy = vy_amb + delta_vy
            vz = vz_amb
            temp = T_amb + delta_T
            dn = rh_amb * temp**(1./aindm1)

            Q_r(i,j,k,rh,1) = dn
            Q_r(i,j,k,mx,1) = dn*vx
            Q_r(i,j,k,my,1) = dn*vy
            Q_r(i,j,k,mz,1) = dn*vz
            Q_r(i,j,k,en,1) = temp*dn/aindm1 + 0.5*dn*(vx**2 + vy**2 + vz**2)

        end do
        end do
        end do

    end subroutine isentropic_vortex

    !-------------------------------------------------------

    subroutine sod_shock_tube(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k
        real rh_amb,vx_amb,vy_amb,vz_amb,p_amb,T_amb,beta
        real xctr,yctr,zctr,xp,yp,r2,delta_vx,delta_vy,delta_T
        real dn,vx,vy,vz,temp

        beta = 5.0  ! vortex strength
        rh_left  = 1.0
        rh_right = 0.125
        p_left  =  1.0
        p_right =  0.1
        vx_amb = 0.0
        vy_amb = 0.0
        vz_amb = 0.0
        T_amb  = T_base
        p_amb  = T_amb*rh_amb

        xctr = 0

        do i = 1,nx
        do j = 1,ny
        do k = 1,nz

            xp = xc(i) - xctr
            yp = yc(j) - yctr
            r2 = xp**2 + yp**2

            delta_vx = -yp*beta/(2*pi) * exp( (1 - r2)/2 )
            delta_vy =  xp*beta/(2*pi) * exp( (1 - r2)/2 )
            delta_T  = -aindm1*beta**2/(8*aindex*pi**2) * exp(1 - r2)

            vx = vx_amb + delta_vx
            vy = vy_amb + delta_vy
            vz = vz_amb

            temp = T_amb + delta_T
            dn = rh_amb * temp**(1./aindm1)

            Q_r(i,j,k,rh,1) = dn
            Q_r(i,j,k,mx,1) = dn*vx
            Q_r(i,j,k,my,1) = dn*vy
            Q_r(i,j,k,mz,1) = dn*vz
            Q_r(i,j,k,en,1) = temp*dn/aindm1 + 0.5*dn*(vx**2 + vy**2 + vz**2)

        end do
        end do
        end do

    end subroutine sod_shock_tube

    !-------------------------------------------------------

    subroutine fill_fluid2(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,igrid,ieq
        real x,y,wtev,rnum,rand_num
        real qquad(npg,nQ),xcc,ycc,zcc  ! bfint(npg,nbasis),qquadv(npg)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        wtev = T_floor

        ! test problem is an unstable flow jet in x with velocity perturbations in y

        do i = 1,nx
            do j = 1,ny
                do k = 1,nz

                    call random_number(rand_num)
                    rnum = (rand_num - 0.5)

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
                            Q_r(i,j,k,ieq,ir) = 0.125 * cbasis(ir) * sum(wgt3d(1:npg) &
                                                      * bfvals_int(1:npg,ir)*qquad(1:npg,ieq))
                        end do
                    end do

                end do
            end do
        end do

    end subroutine fill_fluid2


end module ic_mod
