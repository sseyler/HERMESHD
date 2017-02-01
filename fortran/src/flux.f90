module flux

use parameters
use helpers
use boundaries


! NOTE: It might make sense to use global arrays for stochastic
!   stuff at some point, especially if in separate module
! real GRM_x(nface, 1:nx+1, ny,     nz,     3,3)
! real GRM_y(nface, nx,     1:ny+1, nz,     3,3)
! real GRM_z(nface, nx,     ny,     1:nz+1, 3,3)
! real Sflux_x(nface, 1:nx+1, ny,     nz,     3,3)
! real Sflux_y(nface, nx,     1:ny+1, nz,     3,3)
! real Sflux_z(nface, nx,     ny,     1:nz+1, 3,3)


contains

!----------------------------------------------------------------------------------------------
    subroutine get_GRM(GRMpnts_r, npnts)

        implicit none
        integer ife,npnts,b,e,vsl_ndim
        real, dimension(npnts,3,3) :: GRMpnts_r
        real, allocatable :: grn(:)

        vsl_ndim = npnts*3*3
        allocate(grn(vsl_ndim))

        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)

        ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
        do ife = 1,npnts
            e = 9*ife
            b = e - 8
            GRMpnts_r(ife,1,:) = grn(b:b+2)
            GRMpnts_r(ife,2,:) = grn(b+3:b+5)
            GRMpnts_r(ife,3,:) = grn(b+6:e)
        end do

    end subroutine get_GRM
!----------------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------------
    subroutine random_stresses_pnts_r(Spnts_r, npnts)

        implicit none
        integer ife,npnts
        real, dimension(npnts,3,3) :: GRMpnts_r, Spnts_r
        real Gxx,Gyy,Gzz,Gxy,Gxz,Gyz
        real eta_d,zeta_d
        real trG,trGd3,trG_zeta

        eta_d = eta_sd * sqrt_dVdt_i
        zeta_d = zeta_sd * sqrt_dVdt_i
        call get_GRM(GRMpnts_r, npnts)

        ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
        do ife = 1,npnts
            Gxx = GRMpnts_r(ife,1,1)
            Gyy = GRMpnts_r(ife,2,2)
            Gzz = GRMpnts_r(ife,3,3)

            Gxy = sqrt2i*( GRMpnts_r(ife,1,2) + GRMpnts_r(ife,2,1) )
            Gxz = sqrt2i*( GRMpnts_r(ife,1,3) + GRMpnts_r(ife,3,1) )
            Gyz = sqrt2i*( GRMpnts_r(ife,2,3) + GRMpnts_r(ife,3,2) )

            Spnts_r(ife,1,2) = eta_sd*Gxy
            Spnts_r(ife,1,3) = eta_sd*Gxz
            Spnts_r(ife,2,3) = eta_sd*Gyz
            Spnts_r(ife,2,1) = Spnts_r(ife,1,2)
            Spnts_r(ife,3,1) = Spnts_r(ife,1,3)
            Spnts_r(ife,3,2) = Spnts_r(ife,2,3)

            trG = (Gxx + Gyy + Gzz)
            trGd3 = trG/3.0
            trG_zeta = zeta_sd*trG

            Spnts_r(ife,1,1) = eta_sd*(Gxx - trGd3) + trG_zeta
            Spnts_r(ife,2,2) = eta_sd*(Gyy - trGd3) + trG_zeta
            Spnts_r(ife,3,3) = eta_sd*(Gzz - trGd3) + trG_zeta
        end do

    end subroutine random_stresses_pnts_r
!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
    subroutine random_heatflux_pnts_r(Hpnts_r, npnts)

        implicit none
        integer ife,npnts,b,e,vsl_ndim
        real, dimension(npnts,3) :: GRVpnts_r, Hpnts_r
        real Gx,Gy,Gz,kappa_d
        real, allocatable :: grn(:)

        vsl_ndim = npnts*3
        allocate(grn(vsl_ndim))

        kappa_d = kappa_sd * sqrt_dVdt_i

        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)

        ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
        do ife = 1,npnts
            e = 3*ife
            b = e - 2
            GRVpnts_r(ife,1:3) = grn(b:e)

            Gx = GRVpnts_r(ife,1)
            Gy = GRVpnts_r(ife,2)
            Gz = GRVpnts_r(ife,3)

            Hpnts_r(ife,1) = kappa_sd*Gx
            Hpnts_r(ife,2) = kappa_sd*Gy
            Hpnts_r(ife,3) = kappa_sd*Gz
        end do

    end subroutine random_heatflux_pnts_r
!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
    subroutine flux_calc_pnts_r(Qpnts_r,fpnts_r,ixyz,npnts)

        ! Calculate the flux "fpnts_r" in direction "ixyz" (x, y, or z) at a set of
        ! points corresponding to conserved quantities "Qpnts_r":
        !   ixyz=1: x-direction
        !   ixyz=2: y-direction
        !   ixyz=3: z-direction
        implicit none
        integer ife,ixyz,npnts
        real, dimension(npnts,nQ) :: Qpnts_r, fpnts_r
        real dn,dni,vx,vy,vz,P,asqr,fac,Pre,Psol,dx2,Tem,smsq

        real Spnts_r(npnts,3,3), Hpnts_r(npnts,3)
        real Sxx,Syy,Szz,Sxy,Sxz,Syz
        real Qx,Qy,Qz

        Spnts_r(:,:,:) = 0.0
        Hpnts_r(:,:) = 0.0
        if (llns .eq. 1) then
            call random_stresses_pnts_r(Spnts_r, npnts)
        end if

        do ife = 1,npnts

            dn = Qpnts_r(ife,rh)
            dni = 1./dn
            smsq = Qpnts_r(ife,mx)**2 + Qpnts_r(ife,my)**2 + Qpnts_r(ife,mz)**2

            P = (aindex - 1.)*(Qpnts_r(ife,en) - 0.5*dni*smsq)
            if (ieos .eq. 2) P = P_1*(dn**7.2 - 1.) + P_base + P
            if (P < P_floor) P = P_floor

            vx = Qpnts_r(ife,mx)*dni
            vy = Qpnts_r(ife,my)*dni
            vz = Qpnts_r(ife,mz)*dni

            ! NOTE: may not need all values since ixyz choose flux direction
            Sxx = Spnts_r(ife,1,1)
            Syy = Spnts_r(ife,2,2)
            Szz = Spnts_r(ife,3,3)
            Sxy = Spnts_r(ife,1,2)
            Sxz = Spnts_r(ife,1,3)
            Syz = Spnts_r(ife,2,3)

            Qx = Hpnts_r(ife,1)
            Qy = Hpnts_r(ife,2)
            Qz = Hpnts_r(ife,3)

            if(ixyz .eq. 1) then

                fpnts_r(ife,rh) = Qpnts_r(ife,mx)

                ! NOTE: the stress terms may need plus sign in the energy flux!!!
                fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vx + P + Qpnts_r(ife,pxx) - Sxx
                fpnts_r(ife,my) = Qpnts_r(ife,my)*vx     + Qpnts_r(ife,pxy) - Sxy
                fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vx     + Qpnts_r(ife,pxz) - Sxz

                fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vx                          &
                                + ( (Qpnts_r(ife,pxx) - Sxx)*vx                     &
                                +   (Qpnts_r(ife,pxy) - Sxy)*vy                     &
                                +   (Qpnts_r(ife,pxz) - Sxz)*vz )

                fpnts_r(ife,pxx) =  c4d3nu*vx
                fpnts_r(ife,pyy) = -c2d3nu*vx
                fpnts_r(ife,pzz) = -c2d3nu*vx

                fpnts_r(ife,pxy) = nu*vy
                fpnts_r(ife,pxz) = nu*vz
                fpnts_r(ife,pyz) = 0

            end if
            if(ixyz .eq. 2) then

                fpnts_r(ife,rh) = Qpnts_r(ife,mxa(ixyz))

                fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vy     + Qpnts_r(ife,pxy) - Sxy
                fpnts_r(ife,my) = Qpnts_r(ife,my)*vy + P + Qpnts_r(ife,pyy) - Syy
                fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vy     + Qpnts_r(ife,pyz) - Syz

                fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vy                          &
                                + ( (Qpnts_r(ife,pyy) - Syy)*vy                     &
                                +   (Qpnts_r(ife,pxy) - Sxy)*vx                     &
                                +   (Qpnts_r(ife,pyz) - Syz)*vz )

                fpnts_r(ife,pxx) = -c2d3nu*vy
                fpnts_r(ife,pyy) =  c4d3nu*vy
                fpnts_r(ife,pzz) = -c2d3nu*vy

                fpnts_r(ife,pxy) = nu*vx
                fpnts_r(ife,pxz) = 0
                fpnts_r(ife,pyz) = nu*vz

            end if
            if(ixyz .eq. 3) then

                fpnts_r(ife,rh) = Qpnts_r(ife,mz)

                fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vz     + Qpnts_r(ife,pxz) - Sxz
                fpnts_r(ife,my) = Qpnts_r(ife,my)*vz     + Qpnts_r(ife,pyz) - Syz
                fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vz + P + Qpnts_r(ife,pzz) - Szz

                fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vz                          &
                                + ( (Qpnts_r(ife,pzz) - Szz)*vz                     &
                                +   (Qpnts_r(ife,pxz) - Sxz)*vx                     &
                                +   (Qpnts_r(ife,pyz) - Syz)*vy )

                fpnts_r(ife,pxx) = -c2d3nu*vz
                fpnts_r(ife,pyy) = -c2d3nu*vz
                fpnts_r(ife,pzz) =  c4d3nu*vz

                fpnts_r(ife,pxy) = 0.
                fpnts_r(ife,pxz) = nu*vx
                fpnts_r(ife,pyz) = nu*vy

            end if

        end do

    end subroutine

!----------------------------------------------------------------------------------------------

    subroutine flux_cal(Q_r)

        implicit none
        integer i,j,k,ieq,iback,jleft,kdown,i4,i4p,ipnt
        real Qface_x(nfe,nQ), Qface_y(nfe,nQ),Qface_z(nfe,nQ),fface_x(nfe,nQ), fface_y(nfe,nQ), fface_z(nfe,nQ)
        real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
        real cwavex(nfe),cwavey(nfe),cwavez(nfe),Pre,dni,B2, cfbound, cfx(nface,nQ),cfy(nface,nQ),cfz(nface,nQ)
        real fhllc_x(nface,5),fhllc_y(nface,5),fhllc_z(nface,5),fs(nface,nQ),qvin(nQ)

        kroe(:) = 1
        sqrt_dVdt_i = (dVi/dt)**0.5  ! used in calculating random stresses/heat fluxes

        do k=1,nz
            do j=1,ny
            do i=1,nx+1
                iback = i-1

                if (i .gt. 1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_x(ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_r(iback,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (i .eq. 1) then
                    do ieq = 1,nQ
                        Qface_x(1:nface,ieq) = Qxlow_ext(j,k,1:nface,ieq)
                    end do
                end if

                if (i .lt. nx+1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_x(ipnt+nface,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if

                if (i .eq. nx+1) then
                    do ieq = 1,nQ
                        Qface_x(nface+1:nfe,ieq) = Qxhigh_ext(j,k,1:nface,ieq)
                    end do
                end if

                call flux_calc_pnts_r(Qface_x,fface_x,1,nfe)

                if(.not. ihllc) then
                    do i4=1,nfe
                        do ieq=1,nQ
                            qvin(ieq) = Qface_x(i4,ieq)
                        end do
                        cwavex(i4) = cfcal(qvin,1)
                    end do

                    do i4=1,nface
                        cfrx(i4,rh:en) = max(cwavex(i4),cwavex(i4+nface))
                    end do
                end if

                do ieq = 1,nQ
                    do i4=1,nface
                        i4p = i4 + nface
                        flux_x(i4,i,j,k,ieq) = 0.5*(fface_x(i4,ieq) + fface_x(i4p,ieq)) &
                                             - 0.5*cfrx(i4,ieq)*(Qface_x(i4p,ieq) - Qface_x(i4,ieq))
                    end do
                end do

                kroe(1:nface) = 1

                if (ihllc) call flux_hllc(Qface_x,fface_x,fhllc_x,1)

                ! Needs to be done for HLLC and Roe
                if (ihllc) then
                    do ieq = 1,en
                        do i4=1,nface
                            if (kroe(i4) .gt. 0) then
                                flux_x(i4,i,j,k,ieq) = fhllc_x(i4,ieq)
                            end if
                        end do
                    end do
                end if

            end do
            end do
        end do

    !----------------------------------------------------

        do k=1,nz
            do j=1,ny+1
            jleft = j-1
            do i=1,nx

                if (j .gt. 1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_y(ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_r(i,jleft,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (j .eq. 1) then
                    do ieq = 1,nQ
                        Qface_y(1:nface,ieq) = Qylow_ext(i,k,1:nface,ieq)
                    end do
                end if

                if (j .lt. ny+1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_y(ipnt+nface,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (j .eq. ny+1) then
                    do ieq = 1,nQ
                        Qface_y(nface+1:nfe,ieq) = Qyhigh_ext(i,k,1:nface,ieq)
                    end do
                end if

                call flux_calc_pnts_r(Qface_y,fface_y,2,nfe)

                if (.not. ihllc) then
                    do i4=1,nfe
                        do ieq=1,nQ
                            qvin(ieq) = Qface_y(i4,ieq)
                        end do
                        cwavey(i4) = cfcal(qvin,2)
                    end do

                    do i4=1,nface
                        cfry(i4,rh:en) = max(cwavey(i4),cwavey(i4+nface))
                    end do
                end if

                do ieq = 1,nQ
                    do i4=1,nface
                        i4p = i4 + nface
                        flux_y(i4,i,j,k,ieq) = 0.5*(fface_y(i4,ieq) + fface_y(i4p,ieq)) &
                                             - 0.5*cfry(i4,ieq)*(Qface_y(i4p,ieq) - Qface_y(i4,ieq))
                    end do
                end do

                kroe(1:nface) = 1

                if (ihllc) call flux_hllc(Qface_y,fface_y,fhllc_y,2)

                ! Needs to be done for HLLC and Roe
                if (ihllc) then
                    do ieq = 1,en
                        do i4=1,nface
                            if (kroe(i4) .gt. 0) then
                                flux_y(i4,i,j,k,ieq) = fhllc_y(i4,ieq)
                            end if
                        end do
                    end do
                end if

            end do
            end do
        end do

    !-------------------------------------------------------

        do k=1,nz+1
            kdown = k-1
            do j=1,ny
            do i=1,nx

                if (k .gt. 1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_z(ipnt,ieq) = sum( bfvals_zp(ipnt,1:nbasis)   &
                                                    *Q_r(i,j,kdown,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (k .eq. 1) then
                    do ieq = 1,nQ
                        Qface_z(1:nface,ieq) = Qzlow_ext(i,j,1:nface,ieq)
                    end do
                end if

                if (k .lt. nz+1) then
                    do ieq = 1,nQ
                        do ipnt=1,nface
                            Qface_z(ipnt+nface,ieq) = sum( bfvals_zm(ipnt,1:nbasis) &
                                                          *Q_r(i,j,k,ieq,1:nbasis))
                        end do
                    end do
                end if
                if (k .eq. nz+1) then
                    do ieq = 1,nQ
                        Qface_z(nface+1:nfe,ieq) = Qzhigh_ext(i,j,1:nface,ieq)
                    end do
                end if

                call flux_calc_pnts_r(Qface_z,fface_z,3,nfe)

                if (.not. ihllc) then
                    do i4=1,nfe
                        do ieq=1,nQ
                            qvin(ieq) = Qface_z(i4,ieq)
                        end do
                        cwavez(i4) = cfcal(qvin,3)
                    end do

                    do i4=1,nface
                        cfrz(i4,rh:en) = max(cwavez(i4),cwavez(i4+nface))
                    end do
                end if

                do ieq = 1,nQ
                    do i4=1,nface
                        i4p = i4 + nface
                        flux_z(i4,i,j,k,ieq) = 0.5*(fface_z(i4,ieq) + fface_z(i4p,ieq)) &
                                             - 0.5*cfrz(i4,ieq)*(Qface_z(i4p,ieq) - Qface_z(i4,ieq))
                    end do
                end do

                kroe(1:nface) = 1

                if (ihllc) call flux_hllc(Qface_z,fface_z,fhllc_z,3)

                ! Needs to be done for HLLC and Roe
                if (ihllc) then
                    do ieq = 1,en
                        do i4=1,nface
                            if (kroe(i4) .gt. 0) then
                                flux_z(i4,i,j,k,ieq) = fhllc_z(i4,ieq)
                            end if
                        end do
                    end do
                end if

            end do
            end do
        end do

    end subroutine flux_cal

!----------------------------------------------------------------------------------------------

    subroutine flux_hllc(Qlr,flr,fhllc,ixyz)

    !    Compute ion or electron fluxes (density, momentum, energy) using
    !    nonrelativistic HD HLLC approximate Riemann solver developed by Batten, 1997,
    !    *****    "On the Choice of Wavespeeds for the HLLC Riemann Solver"

    !    This takes into account the left and right propagating shocks along with
    !    the contact (tangential) discontinuity.

        implicit none
        real Qlr(nfe,nQ),flr(nfe,nQ),fhllc(nface,5)
        real sm_num(nface),sm_den(nface),qtilde(nface,5),rtrho(nfe),rtrho_i(nface),qsq(nfe)
        real s_lr(nfe),ctilde(nface),hlr(nfe),cslr(nfe),ctsq(nface),Zi, mfact, csfac
        real aq(nfe),bq(nfe),Qstar(nfe,6),fstar(nfe,6),pstar(nface),s_m(nface)
        real rhov(nfe),vlr(nfe,3),plr(nfe),rho_i
        real qslr(nfe),sq_lr(nfe),slrm_i(nfe),B2(nfe),cf(nfe)
        integer ixyz,i4,i4p1,nr,jie,k,k2,ieq,iparr,iperp1,iperp2,ibatten
        integer rhj,mxj,myj,mzj,enj,psj,ivar(5),ipassive,nhll,ib1,ib2

        nhll = 5
        rhj = rh
        mxj = mx
        myj = my
        mzj = mz
        enj = en

        iparr  = mxa(ixyz)
        iperp1 = mya(ixyz)
        iperp2 = mza(ixyz)

        ivar(1) = rhj
        ivar(2) = iparr
        ivar(3) = iperp1
        ivar(4) = iperp2
        ivar(5) = enj

        ! Indices 1:4 denote the left state, while indices 5:8 denote the right state.
        ibatten = 0

        do k=1,nfe
            rhov(k) = Qlr(k,rhj)
            rho_i = 1.0/rhov(k)
            vlr(k,1) = Qlr(k,iparr)*rho_i        ! velocity parallel to direction of flux computation
            vlr(k,2) = Qlr(k,iperp1)*rho_i        ! velocity in perpendicular direction 1
            vlr(k,3) = Qlr(k,iperp2)*rho_i        ! velocity in perpendicular direction 2
            qsq(k) = vlr(k,1)**2 + vlr(k,2)**2 + vlr(k,3)**2
            plr(k) = (aindex - 1.0)*(Qlr(k,enj) - 0.5*rhov(k)*qsq(k))        ! pressure
            if(ieos .eq. 2) plr(k) = P_1*(rhov(k)**7.2 - 1.) + P_base + plr(k)
            rtrho(k) = sqrt(rhov(k))
        end do

        do k=1,nface
            k2 = k + nface
            if(ieos .eq. 2)then
                cslr(k) = vlr(k,1) - sqrt(7.2*P_1*rhov(k)**6.2 + plr(k)*rho_i)       ! lambda_M(Q_l)
                cslr(k2) = vlr(k2,1) + sqrt(7.2*P_1*rhov(k2)**6.2 + plr(k2)*rho_i)       ! lambda_P(Q_r)
            else
                cslr(k) = vlr(k,1) - sqrt(aindex*plr(k)/rhov(k))       ! lambda_M(Q_l)
                cslr(k2) = vlr(k2,1) + sqrt(aindex*plr(k2)/rhov(k2) )       ! lambda_P(Q_r)
            end if
        end do

        if (ibatten .eq. 1) then  ! compute wave speeds using Roe averages following Batten, 1997

            do k=1,nface
                k2 = k + nface
                rtrho_i(k) = 1.0/(rtrho(k) + rtrho(k2))
                qtilde(k,1) = (rtrho(k)*vlr(k,1) + rtrho(k2)*vlr(k2,1))*rtrho_i(k)
                qtilde(k,2) = (rtrho(k)*vlr(k,2) + rtrho(k2)*vlr(k2,2))*rtrho_i(k)
                qtilde(k,3) = (rtrho(k)*vlr(k,3) + rtrho(k2)*vlr(k2,3))*rtrho_i(k)
                qsq(k) = qtilde(k,1)**2 + qtilde(k,2)**2 + qtilde(k,3)**2
                hlr(k) = (Qlr(k,enj) + plr(k))/rhov(k)
                hlr(k2) = (Qlr(k2,enj) + plr(k2))/rhov(k2)
                qtilde(k,4) = (rtrho(k)*hlr(k) + rtrho(k2)*hlr(k2))*rtrho_i(k)
                ctsq(k) = (aindex - 1.0)*(qtilde(k,4) - 0.5*qsq(k))
            end do
            if (minval(ctsq) .ge. 0.0) then
                ctilde = sqrt(ctsq)
                qslr(1:nface) = qtilde(1:nface,1) - ctilde(1:nface)       ! lambda_M(Q_Roe)
                qslr(nface+1:nfe) = qtilde(nface+1:nfe,1) + ctilde(nface+1:nfe)    ! lambda_P(Q_Roe)
            end if
            if (minval(ctsq) .lt. 0.0) then
                ibatten = 0
            end if

        end if

        if (ibatten .eq. 0) then
            do k=1,nface
                k2 = k + nface
                if(ieos .eq. 2)then
                    qslr(k) = vlr(k2,1) - sqrt(7.2*P_1*rhov(k2)**6.2 + plr(k2)*rho_i)       ! lambda_M(Q_r)
                    qslr(k2) = vlr(k,1) + sqrt(7.2*P_1*rhov(k)**6.2 + plr(k)*rho_i)       ! lambda_P(Q_l)
                else
                    qslr(k) = vlr(k2,1) - sqrt(aindex*plr(k2)/rhov(k2))       ! lambda_M(Q_r)
                    qslr(k2) = vlr(k,1) + sqrt(aindex*plr(k)/rhov(k))       ! lambda_P(Q_l)
                    end if
            end do
        end if

        ! Calculate the slow and fast wavespeeds S_L, S_R of the left, right propagating shocks.
        do k=1,nface
            k2 = k + nface
            s_lr(k) = min(cslr(k),qslr(k))         ! S_L = min(lambda_M(Q_l),lambda_M(Q_r or Q_Roe))
            s_lr(k2) = max(cslr(k2),qslr(k2))         ! S_R = max(lambda_P(Q_r),lambda_P(Q_l or Q_Roe))
            sm_num(k) = rhov(k2)*vlr(k2,1)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*vlr(k,1)*(s_lr(k) - vlr(k,1))
            sm_num(k) = sm_num(k) + plr(k) - plr(k2)
            sm_den(k) = rhov(k2)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*(s_lr(k) - vlr(k,1))
        end do
        where (sm_den .eq. 0.0) sm_den = rh_floor

        ! Calculate the wavespeed S_M of the contact discontinuity.

        do k=1,nface
            s_m(k) = sm_num(k)/sm_den(k)                                          ! Eq. (34) of Batten, 1997
            pstar(k) = rhov(k)*(vlr(k,1) - s_lr(k))*(vlr(k,1) - s_m(k)) + plr(k)  ! Eq. (36) of Batten, 1997
        end do


        ! Now, calculate Q_l* and Q_r* in order to calculate F_l* and F_r*.
        do k=1,nfe
            if (k .le. nface) i4 = k
            if (k .gt. nface) i4 = k - nface
            sm_den(1) = s_lr(k) - s_m(i4)        ! S_{L,R} - S_M

            if (sm_den(1) .eq. 0.0) then
                sm_den(1) = rh_floor
            end if

            slrm_i(k) = 1.0/sm_den(1)
            sq_lr(k) = s_lr(k) - vlr(k,1)            ! S_{L,R} - q_{l,r}
            Qstar(k,1) = rhov(k)*sq_lr(k)*slrm_i(k)                              ! Eq. (35) of Batten ,1997
            Qstar(k,2) = (sq_lr(k)*Qlr(k,iparr) + pstar(i4) - plr(k))*slrm_i(k)  ! Eq. (37-39) of Batten, 1997
            Qstar(k,3) = sq_lr(k)*Qlr(k,iperp1)*slrm_i(k)                        ! Eq. (37-39) of Batten, 1997
            Qstar(k,4) = sq_lr(k)*Qlr(k,iperp2)*slrm_i(k)                        ! Eq. (37-39) of Batten, 1997
            Qstar(k,5) = (sq_lr(k)*Qlr(k,enj) - plr(k)*vlr(k,1) + pstar(i4)*s_m(i4))*slrm_i(k)
            ! Eq. (40) of Batten, 1997

            do ieq=1,nhll
                fstar(k,ieq) = flr(k,ivar(ieq)) + s_lr(k)*(Qstar(k,ieq) - Qlr(k,ivar(ieq)))  ! Eq. (29) of Batten, 1997
            end do
        end do

        ! Finally, calculate the HLLC fluxes from F_l, F_l*, F_r*, and F_r...
        do i4 = 1,nface                        ! Use Eq. (26) of Batten ,1997
            if (s_lr(i4) .gt. 0.0) then                               ! if S_L > 0
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = flr(i4,ivar(ieq))        ! F_HLLC = F_l
                end do
            end if
            if (s_lr(i4) .le. 0.0 .and. 0.0 .lt. s_m(i4)) then        ! if S_L <= 0 < S_M
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = fstar(i4,ieq)            ! F_HLLC = F_l*
                end do
            end if
            if (s_m(i4) .le. 0.0 .and. 0.0 .le. s_lr(i4+nface)) then  ! if S_M <= 0 <= S_R
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = fstar(i4+nface,ieq)      ! F_HLLC = F_r*
                end do
            end if
            if (s_lr(i4+nface) .lt. 0.0) then                         ! if S_R < 0
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = flr(i4+nface,ivar(ieq))  ! F_HLLC = F_r
                end do
            end if

        end do

    end subroutine flux_hllc


!----------------------------------------------------------------------------------------------

    subroutine innerintegral(Q_r)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
        integer i,j,k,ieq,ipg,ir
        real Qinner(npg,nQ),finner_x(npg,nQ), finner_y(npg,nQ), finner_z(npg,nQ), int_r(nbastot,nQ)

        integral_r(:,:,:,:,:) = 0.

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

            do ieq = 1,nQ
                do ipg = 1,npg
                    Qinner(ipg,ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                end do
            end do

            call flux_calc_pnts_r(Qinner,finner_x,1,npg)
            call flux_calc_pnts_r(Qinner,finner_y,2,npg)
            call flux_calc_pnts_r(Qinner,finner_z,3,npg)

            do ieq = 1,nQ

                int_r(kx,ieq) = 0.25*cbasis(kx)*dxi*sum(wgt3d(1:npg)*finner_x(1:npg,ieq))
                int_r(ky,ieq) = 0.25*cbasis(ky)*dyi*sum(wgt3d(1:npg)*finner_y(1:npg,ieq))
                int_r(kz,ieq) = 0.25*cbasis(kz)*dzi*sum(wgt3d(1:npg)*finner_z(1:npg,ieq))

                if (nbasis .gt. 4) then
                    int_r(kyz,ieq)  = 0.25*cbasis(kyz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kz)*finner_y(1:npg,ieq)) &
                          + 0.25*cbasis(kyz)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,ky)*finner_z(1:npg,ieq))
                    int_r(kzx,ieq)  = 0.25*cbasis(kzx)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kz)*finner_x(1:npg,ieq)) &
                          + 0.25*cbasis(kzx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kx)*finner_z(1:npg,ieq))
                    int_r(kxy,ieq)  = 0.25*cbasis(kxy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,ky)*finner_x(1:npg,ieq)) &
                          + 0.25*cbasis(kxy)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kx)*finner_y(1:npg,ieq))
                    int_r(kxyz,ieq) = 0.25*cbasis(kxyz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyz)*finner_x(1:npg,ieq)) &
                          + 0.25*cbasis(kxyz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzx)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kxyz)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxy)*finner_z(1:npg,ieq))
                end if

                if (nbasis .gt. 8) then
                    int_r(kyzz,ieq) = 0.25*cbasis(kyzz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzz)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyz)*finner_z(1:npg,ieq))
                    int_r(kzxx,ieq) = 0.25*cbasis(kzxx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxx)*finner_z(1:npg,ieq)) &
                      + 0.25*cbasis(kzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzx)*finner_x(1:npg,ieq))
                    int_r(kxyy,ieq) = 0.25*cbasis(kxyy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyy)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxy)*finner_y(1:npg,ieq))
                    int_r(kyyz,ieq) = 0.25*cbasis(kyyz)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyy)*finner_z(1:npg,ieq)) &
                      + 0.25*cbasis(kyyz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyz)*finner_y(1:npg,ieq))
                    int_r(kzzx,ieq) = 0.25*cbasis(kzzx)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzz)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kzzx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzx)*finner_z(1:npg,ieq))
                    int_r(kxxy,ieq) = 0.25*cbasis(kxxy)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxx)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kxxy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxy)*finner_x(1:npg,ieq))
                    int_r(kyyzz,ieq) = 0.25*cbasis(kyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzz)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyyz)*finner_z(1:npg,ieq))
                    int_r(kzzxx,ieq) = 0.25*cbasis(kzzxx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxx)*finner_z(1:npg,ieq)) &
                      + 0.25*cbasis(kzzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzzx)*finner_x(1:npg,ieq))
                    int_r(kxxyy,ieq) = 0.25*cbasis(kxxyy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyy)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kxxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxxy)*finner_y(1:npg,ieq))
                    int_r(kyzxx,ieq) = 0.25*cbasis(kyzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kyzxx)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzxx)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kyzxx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxxy)*finner_z(1:npg,ieq))
                    int_r(kzxyy,ieq) = 0.25*cbasis(kzxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kzxyy)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxyy)*finner_z(1:npg,ieq)) &
                      + 0.25*cbasis(kzxyy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyyz)*finner_x(1:npg,ieq))
                    int_r(kxyzz,ieq) = 0.25*cbasis(kxyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_z(1:npg,ieq)) &
                      + 0.25*cbasis(kxyzz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyzz)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kxyzz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzzx)*finner_y(1:npg,ieq))
                    int_r(kxyyzz,ieq) = 0.25*cbasis(kxyyzz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyyzz)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kxyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyzz)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kxyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxyy)*finner_z(1:npg,ieq))
                    int_r(kyzzxx,ieq) = 0.25*cbasis(kyzzxx)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzzxx)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kyzzxx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzxx)*finner_z(1:npg,ieq)) &
                      + 0.25*cbasis(kyzzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyzz)*finner_x(1:npg,ieq))
                    int_r(kzxxyy,ieq) = 0.25*cbasis(kzxxyy)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxxyy)*finner_z(1:npg,ieq)) &
                      + 0.25*cbasis(kzxxyy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxyy)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kzxxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzxx)*finner_y(1:npg,ieq))
                    int_r(kxxyyzz,ieq) = 0.25*cbasis(kxxyyzz)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyyzz)*finner_x(1:npg,ieq)) &
                      + 0.25*cbasis(kxxyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzzxx)*finner_y(1:npg,ieq)) &
                      + 0.25*cbasis(kxxyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxxyy)*finner_z(1:npg,ieq))
                end if
            end do

            do ieq = 1,nQ
                do ir=1,nbasis
                    integral_r(i,j,k,ieq,ir) = int_r(ir,ieq)
                end do
            end do

        end do
        end do
        end do

    end subroutine innerintegral


!----------------------------------------------------------------------------------------------

    subroutine glflux

        implicit none
        integer i,j,k,ieq,ir

        do ieq = 1,nQ
        do k = 1,nz
        do j = 1,ny

            do i = 1,nx

                glflux_r(i,j,k,ieq,1) = 0.25*(dxi*(wgt2d(1)*(flux_x(1,i+1,j,k,ieq) - flux_x(1,i,j,k,ieq)))  &
                                            + dyi*(wgt2d(1)*(flux_y(1,i,j+1,k,ieq) - flux_y(1,i,j,k,ieq)))  &
                                            + dzi*(wgt2d(1)*(flux_z(1,i,j,k+1,ieq) - flux_z(1,i,j,k,ieq)))  &
                                            + dxi*(wgt2d(2)*(flux_x(2,i+1,j,k,ieq) - flux_x(2,i,j,k,ieq)))  &
                                            + dyi*(wgt2d(2)*(flux_y(2,i,j+1,k,ieq) - flux_y(2,i,j,k,ieq)))  &
                                            + dzi*(wgt2d(2)*(flux_z(2,i,j,k+1,ieq) - flux_z(2,i,j,k,ieq)))  &
                                            + dxi*(wgt2d(3)*(flux_x(3,i+1,j,k,ieq) - flux_x(3,i,j,k,ieq)))  &
                                            + dyi*(wgt2d(3)*(flux_y(3,i,j+1,k,ieq) - flux_y(3,i,j,k,ieq)))  &
                                            + dzi*(wgt2d(3)*(flux_z(3,i,j,k+1,ieq) - flux_z(3,i,j,k,ieq)))  &
                                            + dxi*(wgt2d(4)*(flux_x(4,i+1,j,k,ieq) - flux_x(4,i,j,k,ieq)))  &
                                            + dyi*(wgt2d(4)*(flux_y(4,i,j+1,k,ieq) - flux_y(4,i,j,k,ieq)))  &
                                            + dzi*(wgt2d(4)*(flux_z(4,i,j,k+1,ieq) - flux_z(4,i,j,k,ieq))))
            end do

            do ir=2,nbasis
                do i = 1,nx

                    glflux_r(i,j,k,ieq,ir) =                                        &
                         wgtbf_xmp(1,2,ir)*flux_x(1,i+1,j,k,ieq) + wgtbf_xmp(1,1,ir)*flux_x(1,i,j,k,ieq)  &
                       + wgtbf_ymp(1,2,ir)*flux_y(1,i,j+1,k,ieq) + wgtbf_ymp(1,1,ir)*flux_y(1,i,j,k,ieq)  &
                       + wgtbf_zmp(1,2,ir)*flux_z(1,i,j,k+1,ieq) + wgtbf_zmp(1,1,ir)*flux_z(1,i,j,k,ieq)  &
                       + wgtbf_xmp(2,2,ir)*flux_x(2,i+1,j,k,ieq) + wgtbf_xmp(2,1,ir)*flux_x(2,i,j,k,ieq)  &
                       + wgtbf_ymp(2,2,ir)*flux_y(2,i,j+1,k,ieq) + wgtbf_ymp(2,1,ir)*flux_y(2,i,j,k,ieq)  &
                       + wgtbf_zmp(2,2,ir)*flux_z(2,i,j,k+1,ieq) + wgtbf_zmp(2,1,ir)*flux_z(2,i,j,k,ieq)  &
                       + wgtbf_xmp(3,2,ir)*flux_x(3,i+1,j,k,ieq) + wgtbf_xmp(3,1,ir)*flux_x(3,i,j,k,ieq)  &
                       + wgtbf_ymp(3,2,ir)*flux_y(3,i,j+1,k,ieq) + wgtbf_ymp(3,1,ir)*flux_y(3,i,j,k,ieq)  &
                       + wgtbf_zmp(3,2,ir)*flux_z(3,i,j,k+1,ieq) + wgtbf_zmp(3,1,ir)*flux_z(3,i,j,k,ieq)  &
                       + wgtbf_xmp(4,2,ir)*flux_x(4,i+1,j,k,ieq) + wgtbf_xmp(4,1,ir)*flux_x(4,i,j,k,ieq)  &
                       + wgtbf_ymp(4,2,ir)*flux_y(4,i,j+1,k,ieq) + wgtbf_ymp(4,1,ir)*flux_y(4,i,j,k,ieq)  &
                       + wgtbf_zmp(4,2,ir)*flux_z(4,i,j,k+1,ieq) + wgtbf_zmp(4,1,ir)*flux_z(4,i,j,k,ieq)  &
                       - integral_r(i,j,k,ieq,ir)

                end do
            end do

        end do
        end do
        end do

    end subroutine glflux

!-----------------------------------------------------------!
!*******************calculate freezing speeds***************!
!-----------------------------------------------------------!
    real function cfcal(Qcf,cases)

        implicit none
        integer cases
        real Qcf(nQ)
        real Pi, Pe, P, B2, ne, cs
        real dn,dni,vx,vy,vz,hx,hy,hz,va2,vf1,lil02,va,fac

        dn = Qcf(rh)
        dni = 1./dn
        vx = Qcf(mx)*dni
        vy = Qcf(my)*dni
        vz = Qcf(mz)*dni

        if (ieos .eq. 1) then
            P = (aindex - 1.)*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
            cs = sqrt(aindex*P*dni)
        end if

        if (ieos .eq. 2) then
            cs = sqrt(7.2*P_1*dn**6.2)
        end if

        select case (cases)
            case (1) !freezing speed in x direction for fluid variable
                cfcal = abs(vx) + cs

            case (2) !freezing speed in y direction for fluid variable
                cfcal = abs(vy) + cs

            case (3) !freezing speed in z direction for fluid variable
                cfcal = abs(vz) + cs
        end select

    end function cfcal


!----------------------------------------------------------------------------------

    subroutine limiter(Q_r)

        implicit none
        integer i, j, k, ieq, ipge, minindex, ir
        real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
        real Qedge(npge,nQ),theta,Qmin(nQ), deltaQ(nQ)
        real epsi, Qrhmin, QPmin, P(npge), Pave, dn, dni, epsiP, thetaj
        real*8 a, b, c

        epsi = rh_floor
        epsiP = rh_floor*T_floor

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            if (Q_r(i,j,k,rh,1) < rh_floor) then
                do ir=2,nbasis
                    Q_r(i,j,k,rh:en,ir) = 0.0
                end do
                Q_r(i,j,k,rh,1) = rh_floor
            else
                do ipge = 1,npge
                    Qedge(ipge,rh) = sum(bf_faces(ipge,1:nbasis)*Q_r(i,j,k,rh,1:nbasis))
                end do

                Qrhmin = minval(Qedge(:,rh))
                if (Qrhmin < epsi) then
                    theta = (epsi-Q_r(i,j,k,rh,1))/(Qrhmin-Q_r(i,j,k,rh,1))
                    if (theta .gt. 1.) then
                        theta = 1.
                    end if

                    if (theta .lt. 0) then
                        theta = 0.
                    end if
                    do ir=2,nbasis
                        Q_r(i,j,k,rh,ir) = abs(theta)*Q_r(i,j,k,rh,ir)
                    end do

                end if

                Pave = (aindex-1.)*(Q_r(i,j,k,en,1) - 0.5*(Q_r(i,j,k,mx,1)**2 + Q_r(i,j,k,my,1)**2 + Q_r(i,j,k,mz,1)**2)/Q_r(i,j,k,rh,1))

                if (Pave < epsiP) then
                    do ir=2,nbasis
                        Q_r(i,j,k,rh:en,ir) = 0.0
                    end do
                else
                    theta = 1.
                    do ipge = 1,npge
                        do ieq = 1,5
                            Qedge(ipge,ieq) = sum(bf_faces(ipge,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                        end do

                        dn = Qedge(ipge,rh)
                        dni = 1./dn
                        P(ipge) = (aindex - 1.)*(Qedge(ipge,en) - 0.5*(Qedge(ipge,mx)**2+Qedge(ipge,my)**2+Qedge(ipge,mz)**2)*dni)

                        if (P(ipge) < epsiP) then
                            if (Pave .ne. P(ipge)) then
                                thetaj = (Pave - epsiP)/(Pave - P(ipge))
                                theta = min(theta,thetaj)
                            end if
                        end if
                    end do
                    if (theta .gt. 1.) then
                        theta = 1.
                    end if

                    if (theta .lt. 0.) then
                        theta = 0.
                    end if
                    do ir=2,nbasis
                        Q_r(i,j,k,rh:en,ir) = theta*Q_r(i,j,k,rh:en,ir)
                    end do
                end if
            end if
        end do
        end do
        end do


    end subroutine limiter


end module flux
