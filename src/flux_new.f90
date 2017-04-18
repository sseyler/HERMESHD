!***** FLUX.F90 **************************************************************************
module flux

use params
use helpers
use spatial
use basis

use boundary
use random

    !===========================================================================
    ! ABSTRACT INTERFACE to subroutine for temporal integration
    !-----------------------------------------------------------------
    abstract interface
        function Fpt_ptr(Qpt)
            use params, only : nQ

            real, dimension(nQ), intent(in)  :: Qpt
            real, dimension(nQ) :: Fpt_ptr
        end function Fpt_ptr
    end interface
    !---------------------------------------------------------------------------

    !===========================================================================
    ! Initialize pointer to temporal integration subroutine
    !-----------------------------------------------------------------
    procedure (Fpt_ptr), pointer :: Fpt_x => null ()
    procedure (Fpt_ptr), pointer :: Fpt_y => null ()
    procedure (Fpt_ptr), pointer :: Fpt_z => null ()
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Global variables
    !-----------------------------------------------------------------
    real, dimension(nx,ny,nz,nQ,nbasis) :: glflux_r, integral_r

    ! Only used by flux_calc (flux_cal) and glflux
    real, dimension(nface,1:nx+1,ny,nz,1:nQ) :: flux_x
    real, dimension(nface,nx,1:ny+1,nz,1:nQ) :: flux_y
    real, dimension(nface,nx,ny,1:nz+1,1:nQ) :: flux_z
    !---------------------------------------------------------------------------


    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------
    ! real, dimension(nface,nQ) :: cfrx,cfry,cfrz         ! in CES's code, these are defined globally
    ! real, dimension(nfe,nQ) :: Qface_x,Qface_y,Qface_z  ! in CES's code, these are inside flux_cal
    ! real, dimension(nfe,nQ) :: fface_x,fface_y,fface_z  ! in CES's code, these are inside flux_cal

contains

    !===========================================================================
    ! Select user-specified integration method at runtime
    !-----------------------------------------------------------------
    subroutine select_hydro_model(ivis, flux_x_pt, flux_y_pt, flux_z_pt)
        implicit none
        integer, intent(in) :: ivis
        procedure(Fpt_ptr), pointer :: flux_x_pt,flux_y_pt,flux_z_pt

        select case (ivis)
            case (0)
                call mpi_print(iam, 'ivis = 0: Selected Euler equations (inviscid fluid)')
                flux_x_pt => Fpt_x_0
                flux_y_pt => Fpt_y_0
                flux_z_pt => Fpt_z_0
            case (1)
                call mpi_print(iam, 'ivis = 1: Selected linearized 10-moment equations')
                flux_x_pt => Fpt_x_1
                flux_y_pt => Fpt_y_1
                flux_z_pt => Fpt_z_1
            case (2)
                ! WARNING/TODO
                call mpi_print(iam, 'Error (ivis = 2): Full 10-moment equations not implemented')
                ! flux_x_pt => Fpt_x_2
                ! flux_y_pt => Fpt_y_2
                ! flux_z_pt => Fpt_z_2
                call exit(-1)
            case default
                call mpi_print(iam, 'ivis not set: Defaulting to Euler equations (inviscid fluid)')
                flux_x_pt => Fpt_x_0
                flux_y_pt => Fpt_y_0
                flux_z_pt => Fpt_z_0
        end select
    end subroutine select_hydro_model
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! subroutine Fpts_x(Qpnts, Fpnts, npnts)
    !     implicit none
    !     real, dimension(npnts,nQ), intent(in)  :: Qpnts
    !     real, dimension(npnts,nQ), intent(out) :: Fpnts
    !     real, dimension(nQ) :: Qpt, Fpt
    !     integer ife
    !
    !     do ife = 1,npnts
    !         Fpnts(ife,1:nQ) = Fpt_x(Qpnts(ife,1:nQ))
    !     end do
    ! end subroutine Fpts_x
    !---------------------------------------------------------------------------

    !===========================================================================
    ! Euler equation (inviscid) fluxes
    !---------------------------------------------------------------------------
    function Fpt_x_0(Qpt)
        ! WARNING: the stochastic terms have been excluded here!!!
        real, dimension(nQ), intent(in)  :: Qpt
        real, dimension(nQ) :: Fpt_x_0
        real :: Fpt(nQ), dn,dni,vx,vy,vz,P

        dn = Qpt(rh)
        dni = 1./dn
        vx = Qpt(mx)*dni
        vy = Qpt(my)*dni
        vz = Qpt(mz)*dni
        P  = aindm1 * ( Qpt(en) - 0.5*dn*(vx**2 + vy**2 + vz**2) )

        Fpt(rh) = Qpt(mx)
        Fpt(mx) = Qpt(mx)*vx + P
        Fpt(my) = Qpt(my)*vx
        Fpt(mz) = Qpt(mz)*vx
        Fpt(en) = (Qpt(en) + P)*vx

        Fpt_x_0 = Fpt(:)
    end function Fpt_x_0
    !---------------------------------------------------------------------------
    function Fpt_y_0(Qpt)
        ! WARNING: the stochastic terms have been excluded here!!!
        real, dimension(nQ), intent(in)  :: Qpt
        real, dimension(nQ) :: Fpt_y_0
        real :: Fpt(nQ), dn,dni,vx,vy,vz,P

        dn = Qpt(rh)
        dni = 1./dn
        vx = Qpt(mx)*dni
        vy = Qpt(my)*dni
        vz = Qpt(mz)*dni
        P  = aindm1 * ( Qpt(en) - 0.5*dn*(vx**2 + vy**2 + vz**2) )

        Fpt(rh) = Qpt(my)
        Fpt(mx) = Qpt(mx)*vy
        Fpt(my) = Qpt(my)*vy + P
        Fpt(mz) = Qpt(mz)*vy
        Fpt(en) = (Qpt(en) + P)*vy

        Fpt_y_0 = Fpt(:)
    end function Fpt_y_0
    !---------------------------------------------------------------------------
    function Fpt_z_0(Qpt)
        ! WARNING: the stochastic terms have been excluded here!!!
        real, dimension(nQ), intent(in)  :: Qpt
        real, dimension(nQ) :: Fpt_z_0
        real :: Fpt(nQ), dn,dni,vx,vy,vz,P

        dn = Qpt(rh)
        dni = 1./dn
        vx = Qpt(mx)*dni
        vy = Qpt(my)*dni
        vz = Qpt(mz)*dni
        P  = aindm1 * ( Qpt(en) - 0.5*dn*(vx**2 + vy**2 + vz**2) )

        Fpt(rh) = Qpt(mz)
        Fpt(mx) = Qpt(mx)*vz
        Fpt(my) = Qpt(my)*vz
        Fpt(mz) = Qpt(mz)*vz + P
        Fpt(en) = (Qpt(en) + P)*vz

        Fpt_z_0 = Fpt(:)
    end function Fpt_z_0
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Linearized 10-moment fluxes
    !---------------------------------------------------------------------------
    function Fpt_x_1(Qpt)
        ! WARNING: the stochastic terms have been excluded here!!!
        real, dimension(nQ), intent(in)  :: Qpt
        real, dimension(nQ) :: Fpt_x_1
        real :: Fpt(nQ), dn,dni,vx,vy,vz,P,P_xx,P_xy,P_xz

        dn = Qpt(rh)
        dni = 1./dn
        vx = Qpt(mx)*dni
        vy = Qpt(my)*dni
        vz = Qpt(mz)*dni
        P    = aindm1 * ( Qpt(en) - 0.5*dn*(vx**2 + vy**2 + vz**2) )
        P_xx = Qpt(exx) - dn*vx**2
        P_xy = Qpt(exy) - dn*vx*vy
        P_xz = Qpt(exz) - dn*vx*vz

        Fpt(rh)  = Qpt(mx)
        Fpt(mx)  = Qpt(mx)*vx + P + Qpt(exx)
        Fpt(my)  = Qpt(my)*vx     + Qpt(exy)
        Fpt(mz)  = Qpt(mz)*vx     + Qpt(exz)
        Fpt(en)  = (Qpt(en) + P_xx)*vx + P_xy*vy + P_xz*vz

        Fpt(exx) =  c4d3cv*vx
        Fpt(eyy) = -c2d3cv*vx
        Fpt(ezz) = -c2d3cv*vx
        Fpt(exy) =  colvis*vy
        Fpt(exz) =  colvis*vz
        Fpt(eyz) =  0

        Fpt_x_1 = Fpt(:)
    end function Fpt_x_1
    !---------------------------------------------------------------------------
    function Fpt_y_1(Qpt)
        ! WARNING: the stochastic terms have been excluded here!!!
        real, dimension(nQ), intent(in)  :: Qpt
        real, dimension(nQ) :: Fpt_y_1
        real :: Fpt(nQ), dn,dni,vx,vy,vz,P,P_xy,P_yy,P_yz

        dn = Qpt(rh)
        dni = 1./dn
        vx = Qpt(mx)*dni
        vy = Qpt(my)*dni
        vz = Qpt(mz)*dni
        P    = aindm1 * ( Qpt(en) - 0.5*dn*(vx**2 + vy**2 + vz**2) )
        P_xy = Qpt(exy) - dn*vx*vy
        P_yy = Qpt(eyy) - dn*vy**2
        P_yz = Qpt(eyz) - dn*vy*vz

        Fpt(rh)  = Qpt(my)
        Fpt(mx)  = Qpt(mx)*vy     + Qpt(exy)
        Fpt(my)  = Qpt(my)*vy + P + Qpt(eyy)
        Fpt(mz)  = Qpt(mz)*vy     + Qpt(eyz)
        Fpt(en)  = (Qpt(en) + P_yy)*vy + P_xy*vx + P_yz*vz

        Fpt(exx) = -c2d3cv*vy
        Fpt(eyy) =  c4d3cv*vy
        Fpt(ezz) = -c2d3cv*vy
        Fpt(exy) =  colvis*vx
        Fpt(exz) =  0
        Fpt(eyz) =  colvis*vz

        Fpt_y_1 = Fpt(:)
    end function Fpt_y_1
    !---------------------------------------------------------------------------
    function Fpt_z_1(Qpt)
        ! WARNING: the stochastic terms have been excluded here!!!
        real, dimension(nQ), intent(in)  :: Qpt
        real, dimension(nQ) :: Fpt_z_1
        real :: Fpt(nQ), dn,dni,vx,vy,vz,P,P_xz,P_yz,P_zz

        dn = Qpt(rh)
        dni = 1./dn
        vx = Qpt(mx)*dni
        vy = Qpt(my)*dni
        vz = Qpt(mz)*dni
        P    = aindm1 * ( Qpt(en) - 0.5*dn*(vx**2 + vy**2 + vz**2) )
        P_xz = Qpt(exz) - dn*vx*vz
        P_yz = Qpt(eyy) - dn*vy*vz
        P_zz = Qpt(ezz) - dn*vz**2

        Fpt(rh)  = Qpt(mz)
        Fpt(mx)  = Qpt(mx)*vz     + Qpt(exz)
        Fpt(my)  = Qpt(my)*vz     + Qpt(eyz)
        Fpt(mz)  = Qpt(mz)*vz + P + Qpt(ezz)
        Fpt(en)  = (Qpt(en) + P_zz)*vz + P_xz*vx + P_yz*vy

        Fpt(exx) = -c2d3cv*vz
        Fpt(eyy) = -c2d3cv*vz
        Fpt(ezz) =  c4d3cv*vz
        Fpt(exy) =  0
        Fpt(exz) =  colvis*vx
        Fpt(eyz) =  colvis*vy

        Fpt_z_1 = Fpt(:)
    end function Fpt_z_1
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! Get field var values at a cell face in body (not boundary) of spatial region
    !---------------------------------------------------------------------------
    function Qface_body(Q_ijk, bfv_r)
        real, dimension(nQ,nbasis), intent(in) :: Q_ijk
        real, dimension(nface,nQ), intent(in)  :: bfv_r
        real, dimension(nface,nQ) :: Qface_body
        integer ieq,ipnt

        do ieq = 1,nQ
            do ipnt=1,nface
                Qface_body(ipnt,ieq) = sum(bfv_r(ipnt,1:nbasis)*Q_ijk(ieq,1:nbasis))
            end do
        end do
    end function Qface_body
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Get field var values at a cell face at the edge (boundary) of spatial region
    !---------------------------------------------------------------------------
    function Qface_edge(Q_bc_ext)
        real, dimension(nface,nQ), intent(in)  :: Q_bc_ext
        real, dimension(nface,nQ) :: Qface_edge
        integer ieq,ipnt

        do ieq = 1,nQ
            do ipnt=1,nface
                Qface_edge(ipnt,ieq) = Q_bc_ext(ipnt,ieq)
            end do
        end do
    end function Qface_edge
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    subroutine calc_flux_x(Q_r, flux_x)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r
        real, dimension(nface,nx+1,ny,nz,nQ), intent(out) :: flux_x
        real, dimension(nfe,nQ) :: Qface_x,Fface_x
        real, dimension(nface,nQ) :: cfrx

        integer i,j,k,ieq,iback,i4,i4p,ipnt,ife
        real cwavex(nfe),fhllc_x(nface,5),qvin(nQ)

        do k=1,nz
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iback,Qface_x,qvin) FIRSTPRIVATE(Fface_x,cfrx,cwavex,fhllc_x)
        do j=1,ny
        do i=1,nx+1
            iback = i-1

            !--- get Qface_z --------------------------
            if (i > 1)     Qface_x(1:nface,1:nQ)     = Qface_body( Q_r(iback,j,k,:,:), bfvals_xp )
            if (i == 1)    Qface_x(1:nface,1:nQ)     = Qface_edge( Qxlo_ext(j,k,:,:) )
            if (i <  nx+1) Qface_x(nface+1:nfe,1:nQ) = Qface_body( Q_r(i,j,k,:,:),     bfvals_xm )
            if (i == nx+1) Qface_x(nface+1:nfe,1:nQ) = Qface_edge( Qxhi_ext(j,k,:,:) )

            !--- get Fface_z --------------------------
            ! call Fpts_x(Qface_x, Fface_x, nfe)
            do ife = 1,nfe
                Fface_x(ife,1:nQ) = Fpt_x( Qface_x(ife,:) )
            end do

            !--- get cfrz -----------------------------
            if (ivis == 2) then
                do i4=1,nfe
                    do ieq=1,nQ
                        qvin(ieq) = Qface_x(i4,ieq)
                    end do
                    cwavex(i4) = cfcal(qvin,1)
                end do
                do i4=1,nface
                    cfrx(i4,1:nQ) = max(cwavex(i4),cwavex(i4+nface))
                end do
            end if

            !--- get flux_z ---------------------------
            do ieq = 1,nQ
                do i4=1,nface
                    i4p = i4 + nface
                    flux_x(i4,i,j,k,ieq) = 0.5*(Fface_x(i4,ieq) + Fface_x(i4p,ieq)) &
                                         - 0.5*cfrx(i4,ieq)*(Qface_x(i4p,ieq) - Qface_x(i4,ieq))
                end do
            end do

            !--- get HLLC flux ------------------------
            if (ihllc) then ! Needs to be done for HLLC and Roe
                call flux_hllc(Qface_x,Fface_x,fhllc_x,1)
                do ieq = 1,en
                    do i4=1,nface
                        flux_x(i4,i,j,k,ieq) = fhllc_x(i4,ieq)
                    end do
                end do
            end if

        end do
        end do
        !$OMP END PARALLEL DO
        end do

    end subroutine calc_flux_x
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    subroutine calc_flux_y(Q_r, flux_y)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r
        real, dimension(nface,nx,ny+1,nz,nQ), intent(out) :: flux_y
        real, dimension(nfe,nQ) :: Qface_y,Fface_y
        real, dimension(nface,nQ) :: cfry

        integer i,j,k,ieq,jleft,i4,i4p,ipnt,ife
        real cwavey(nfe),fhllc_y(nface,5),qvin(nQ)

        do k=1,nz
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jleft,Qface_y,qvin) FIRSTPRIVATE(Fface_y,cfry,cwavey,fhllc_y)
        do j=1,ny+1
        jleft = j-1
        do i=1,nx
            !--- get Qface_z --------------------------
            if (j > 1)     Qface_y(1:nface,1:nQ)     = Qface_body( Q_r(i,jleft,k,:,:), bfvals_yp )
            if (j == 1)    Qface_y(1:nface,1:nQ)     = Qface_edge( Qylo_ext(i,k,:,:) )
            if (j <  ny+1) Qface_y(nface+1:nfe,1:nQ) = Qface_body( Q_r(i,j,k,:,:),     bfvals_ym )
            if (j == ny+1) Qface_y(nface+1:nfe,1:nQ) = Qface_edge( Qyhi_ext(i,k,:,:) )

            !--- get Fface_z --------------------------
            ! call Fpts_y(Qface_y, Fface_y, nfe)
            do ife = 1,nfe
                Fface_y(ife,1:nQ) = Fpt_y(Qface_y(ife,:))
            end do

            !--- get cfrz -----------------------------
            if (ivis == 2) then
                do i4=1,nfe
                    do ieq=1,nQ
                        qvin(ieq) = Qface_y(i4,ieq)
                    end do
                    cwavey(i4) = cfcal(qvin,2)
                end do
                do i4=1,nface
                    cfry(i4,1:nQ) = max(cwavey(i4),cwavey(i4+nface))
                end do
            end if

            !--- get flux_z ---------------------------
            do ieq = 1,nQ
                do i4=1,nface
                    i4p = i4 + nface
                    flux_y(i4,i,j,k,ieq) = 0.5*(Fface_y(i4,ieq) + Fface_y(i4p,ieq)) &
                                         - 0.5*cfry(i4,ieq)*(Qface_y(i4p,ieq) - Qface_y(i4,ieq))
                end do
            end do

            !--- get HLLC flux ------------------------
            if (ihllc) then ! Needs to be done for HLLC and Roe
                call flux_hllc(Qface_y,Fface_y,fhllc_y,2)
                do ieq = 1,en
                    do i4=1,nface
                        flux_y(i4,i,j,k,ieq) = fhllc_y(i4,ieq)
                    end do
                end do
            end if

        end do
        end do
        !$OMP END PARALLEL DO
        end do

    end subroutine calc_flux_y
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    subroutine calc_flux_z(Q_r, flux_z)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r
        real, dimension(nface,nx,ny,nz+1,nQ), intent(out) :: flux_z
        real, dimension(nfe,nQ) :: Qface_z,Fface_z
        real, dimension(nface,nQ) :: cfrz

        integer i,j,k,ieq,kdown,i4,i4p,ipnt,ife
        real cwavez(nfe),fhllc_z(nface,5),qvin(nQ)

        do k=1,nz+1
        kdown = k-1
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Qface_z,qvin) FIRSTPRIVATE(Fface_z,cfrz,cwavez,fhllc_z)
        do j=1,ny
        do i=1,nx
            !--- get Qface_z --------------------------
            if (k > 1)     Qface_z(1:nface,1:nQ)     = Qface_body( Q_r(i,j,kdown,:,:), bfvals_zp )
            if (k == 1)    Qface_z(1:nface,1:nQ)     = Qface_edge( Qzlo_ext(i,k,:,:) )
            if (k <  nz+1) Qface_z(nface+1:nfe,1:nQ) = Qface_body( Q_r(i,j,k,:,:),     bfvals_zm )
            if (k == nz+1) Qface_z(nface+1:nfe,1:nQ) = Qface_edge( Qzhi_ext(i,k,:,:) )

            !--- get Fface_z --------------------------
            ! call Fpts_z(Qface_z, Fface_z, nfe)
            do ife = 1,nfe
                Fface_z(ife,1:nQ) = Fpt_z(Qface_z(ife,:))
            end do

            !--- get cfrz -----------------------------
            if (ivis == 2) then
                do i4=1,nfe
                    do ieq=1,nQ
                        qvin(ieq) = Qface_z(i4,ieq)
                    end do
                    cwavez(i4) = cfcal(qvin,3)
                end do
                do i4=1,nface
                    cfrz(i4,1:nQ) = max(cwavez(i4),cwavez(i4+nface))
                end do
            end if

            !--- get flux_z ---------------------------
            do ieq = 1,nQ
                do i4=1,nface
                    i4p = i4 + nface
                    flux_z(i4,i,j,k,ieq) = 0.5*(Fface_z(i4,ieq) + Fface_z(i4p,ieq)) &
                                         - 0.5*cfrz(i4,ieq)*(Qface_z(i4p,ieq) - Qface_z(i4,ieq))
                end do
            end do

            !--- get HLLC flux ------------------------
            if (ihllc) then
                call flux_hllc(Qface_z,Fface_z,fhllc_z,3)
                do ieq = 1,en
                    do i4=1,nface
                        flux_z(i4,i,j,k,ieq) = fhllc_z(i4,ieq)
                    end do
                end do
            end if

        end do
        end do
        !$OMP END PARALLEL DO
        end do

    end subroutine calc_flux_z
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    subroutine flux_calc(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r

        ! NOTE: "nu" in CES' code is te*rh, which is "colvis" here; "coll" is colvis/vis
        c2d3cv = c2d3*colvis  ! NOTE: globaly defined in params.f90
        c4d3cv = c4d3*colvis  ! NOTE: globaly defined in params.f90

        call calc_flux_x(Q_r, flux_x)
        call calc_flux_y(Q_r, flux_y)
        call calc_flux_z(Q_r, flux_z)
    end subroutine flux_calc
    !---------------------------------------------------------------------------


!-------------------------------------------------------------------------------
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

        ! iparr  = mxa(ixyz)
        ! iperp1 = mya(ixyz)
        ! iperp2 = mza(ixyz)
        select case (ixyz)
        case (1)
            iparr  = mx
            iperp1 = my
            iperp2 = mz
        case (2)
            iparr  = my
            iperp1 = mz
            iperp2 = mx
        case (3)
            iparr  = mz
            iperp1 = mx
            iperp2 = my
        end select

        nhll = 5
        rhj = rh
        mxj = mx
        myj = my
        mzj = mz
        enj = en
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
            plr(k) = aindm1*(Qlr(k,enj) - 0.5*rhov(k)*qsq(k))        ! pressure
            if(ieos == 2) plr(k) = P_1*(rhov(k)**7.2 - 1.) + P_base + plr(k)
            rtrho(k) = sqrt(rhov(k))
        end do

        do k=1,nface
            k2 = k + nface
            if(ieos == 2)then
                cslr(k) = vlr(k,1) - sqrt(7.2*P_1*rhov(k)**6.2 + plr(k)*rho_i)       ! lambda_M(Q_l)
                cslr(k2) = vlr(k2,1) + sqrt(7.2*P_1*rhov(k2)**6.2 + plr(k2)*rho_i)       ! lambda_P(Q_r)
            else
                cslr(k) = vlr(k,1) - sqrt(aindex*plr(k)/rhov(k))       ! lambda_M(Q_l)
                cslr(k2) = vlr(k2,1) + sqrt(aindex*plr(k2)/rhov(k2) )       ! lambda_P(Q_r)
            end if
        end do

        if (ibatten == 1) then  ! compute wave speeds using Roe averages following Batten, 1997

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
                ctsq(k) = aindm1*(qtilde(k,4) - 0.5*qsq(k))
            end do
            if (minval(ctsq) >= 0.0) then
                ctilde = sqrt(ctsq)
                qslr(1:nface) = qtilde(1:nface,1) - ctilde(1:nface)       ! lambda_M(Q_Roe)
                qslr(nface+1:nfe) = qtilde(nface+1:nfe,1) + ctilde(nface+1:nfe)    ! lambda_P(Q_Roe)
            end if
            if (minval(ctsq) < 0.0) then
                ibatten = 0
            end if

        end if

        if (ibatten == 0) then
            do k=1,nface
                k2 = k + nface
                if(ieos == 2)then
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
        where (sm_den == 0.0) sm_den = rh_floor

        ! Calculate the wavespeed S_M of the contact discontinuity.

        do k=1,nface
            s_m(k) = sm_num(k)/sm_den(k)                                          ! Eq. (34) of Batten, 1997
            pstar(k) = rhov(k)*(vlr(k,1) - s_lr(k))*(vlr(k,1) - s_m(k)) + plr(k)  ! Eq. (36) of Batten, 1997
        end do


        ! Now, calculate Q_l* and Q_r* in order to calculate F_l* and F_r*.
        do k=1,nfe
            if (k <= nface) i4 = k
            if (k > nface) i4 = k - nface
            sm_den(1) = s_lr(k) - s_m(i4)        ! S_{L,R} - S_M

            if (sm_den(1) == 0.0) then
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
            if (s_lr(i4) > 0.0) then                               ! if S_L > 0
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = flr(i4,ivar(ieq))        ! F_HLLC = F_l
                end do
            end if
            if (s_lr(i4) <= 0.0 .and. 0.0 < s_m(i4)) then        ! if S_L <= 0 < S_M
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = fstar(i4,ieq)            ! F_HLLC = F_l*
                end do
            end if
            if (s_m(i4) <= 0.0 .and. 0.0 <= s_lr(i4+nface)) then  ! if S_M <= 0 <= S_R
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = fstar(i4+nface,ieq)      ! F_HLLC = F_r*
                end do
            end if
            if (s_lr(i4+nface) < 0.0) then                         ! if S_R < 0
                do ieq=1,nhll
                    fhllc(i4,ivar(ieq)) = flr(i4+nface,ivar(ieq))  ! F_HLLC = F_r
                end do
            end if

        end do

    end subroutine flux_hllc
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine innerintegral(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r
        real, dimension(npg,nQ)     :: Qinner, finner_x,finner_y,finner_z
        real, dimension(nbastot,nQ) :: int_r
        real sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9
        integer i,j,k,ieq,ipg,ir

        integral_r(:,:,:,:,:) = 0.

        do k = 1,nz
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Qinner,integral_r,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum8) FIRSTPRIVATE(finner_x,finner_y,finner_z,int_r)
        do j = 1,ny
        do i = 1,nx

            do ieq = 1,nQ
                do ipg = 1,npg
                    Qinner(ipg,ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                end do
            end do

            ! NOTE: this can be made more efficient since temp vars can be re-used
            do ipg = 1,npg
                finner_x(ipg,1:nQ) = Fpt_x( Qinner(ipg,:) )
                finner_y(ipg,1:nQ) = Fpt_y( Qinner(ipg,:) )
                finner_z(ipg,1:nQ) = Fpt_z( Qinner(ipg,:) )
            end do
            ! call flux_calc_pnts_r(Qinner,finner_x,1,npg)
            ! call flux_calc_pnts_r(Qinner,finner_y,2,npg)
            ! call flux_calc_pnts_r(Qinner,finner_z,3,npg)

            do ieq = 1,nQ

                ! int_r(kx,ieq) = 0.25*cbasis(kx)*dxi*sum(wgt3d(1:npg)*finner_x(1:npg,ieq))
                ! int_r(ky,ieq) = 0.25*cbasis(ky)*dyi*sum(wgt3d(1:npg)*finner_y(1:npg,ieq))
                ! int_r(kz,ieq) = 0.25*cbasis(kz)*dzi*sum(wgt3d(1:npg)*finner_z(1:npg,ieq))
                sum1 = 0.
                sum2 = 0.
                sum3 = 0.
                do ipg=1,npg
                    sum1 = sum1 + wgt3d(ipg)*finner_x(ipg,ieq)
                    sum2 = sum2 + wgt3d(ipg)*finner_y(ipg,ieq)
                    sum3 = sum3 + wgt3d(ipg)*finner_z(ipg,ieq)
                end do
                int_r(kx,ieq) = 0.25*cbasis(kx)*dxi*sum1
                int_r(ky,ieq) = 0.25*cbasis(ky)*dyi*sum2
                int_r(kz,ieq) = 0.25*cbasis(kz)*dzi*sum3


                if ( nbasis > 4 ) then

                    if ( ibitri == 1 .or. iquad > 2 ) then
                        sum1 = 0.
                        sum2 = 0.
                        sum3 = 0.
                        sum4 = 0.
                        sum5 = 0.
                        sum6 = 0.
                        do ipg=1,npg
                            sum1 = sum1 + wgt3d(ipg)*bfvals_int(ipg,kz)*finner_y(ipg,ieq)
                            sum2 = sum2 + wgt3d(ipg)*bfvals_int(ipg,ky)*finner_z(ipg,ieq)
                            sum3 = sum3 + wgt3d(ipg)*bfvals_int(ipg,kz)*finner_x(ipg,ieq)
                            sum4 = sum4 + wgt3d(ipg)*bfvals_int(ipg,kx)*finner_z(ipg,ieq)
                            sum5 = sum5 + wgt3d(ipg)*bfvals_int(ipg,ky)*finner_x(ipg,ieq)
                            sum6 = sum6 + wgt3d(ipg)*bfvals_int(ipg,kx)*finner_y(ipg,ieq)
                        end do
                        int_r(kyz,ieq) = 0.25*cbasis(kyz)*(dyi*sum1 + dzi*sum2)
                        int_r(kzx,ieq) = 0.25*cbasis(kzx)*(dxi*sum3 + dzi*sum4)
                        int_r(kxy,ieq) = 0.25*cbasis(kxy)*(dxi*sum5 + dyi*sum6)
                    end if

                    if ( iquad > 2 ) then
                        sum1 = 0.
                        sum2 = 0.
                        sum3 = 0.
                        do ipg=1,npg
                            sum1 = sum1 + wgt3d(ipg)*3.*bfvals_int(ipg,kx)*finner_x(ipg,ieq)
                            sum2 = sum2 + wgt3d(ipg)*3.*bfvals_int(ipg,ky)*finner_y(ipg,ieq)
                            sum3 = sum3 + wgt3d(ipg)*3.*bfvals_int(ipg,kz)*finner_z(ipg,ieq)
                        end do

                        int_r(kxx,ieq) = 0.25*cbasis(kxx)*(dxi*sum1)
                        int_r(kyy,ieq) = 0.25*cbasis(kyy)*(dyi*sum2)
                        int_r(kzz,ieq) = 0.25*cbasis(kzz)*(dzi*sum3)
                    end if

                    if ( ibitri == 1 .or. iquad > 3 ) then
                        sum7 = 0.
                        sum8 = 0.
                        sum9 = 0.
                        do ipg=1,npg
                            sum7 = sum7 + wgt3d(ipg)*bfvals_int(ipg,kyz)*finner_x(ipg,ieq)
                            sum8 = sum8 + wgt3d(ipg)*bfvals_int(ipg,kzx)*finner_y(ipg,ieq)
                            sum9 = sum9 + wgt3d(ipg)*bfvals_int(ipg,kxy)*finner_z(ipg,ieq)
                        end do
                        int_r(kxyz,ieq) = 0.25*cbasis(kxyz)*(dxi*sum7 + dyi*sum8 + dzi*sum9)
                    end if

                    if ( (ibitri == 1 .and. iquad > 2) .or. iquad > 3 ) then
                        int_r(kyzz,ieq) =                                                                       &
                            0.25*cbasis(kyzz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzz)*finner_y(1:npg,ieq))   &
                          + 0.25*cbasis(kyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyz)*finner_z(1:npg,ieq))
                        int_r(kzxx,ieq) =                                                                       &
                            0.25*cbasis(kzxx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxx)*finner_z(1:npg,ieq))   &
                          + 0.25*cbasis(kzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzx)*finner_x(1:npg,ieq))
                        int_r(kxyy,ieq) =                                                                       &
                            0.25*cbasis(kxyy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyy)*finner_x(1:npg,ieq))   &
                          + 0.25*cbasis(kxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxy)*finner_y(1:npg,ieq))
                        int_r(kyyz,ieq) =                                                                       &
                            0.25*cbasis(kyyz)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyy)*finner_z(1:npg,ieq))   &
                          + 0.25*cbasis(kyyz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyz)*finner_y(1:npg,ieq))
                        int_r(kzzx,ieq) =                                                                       &
                            0.25*cbasis(kzzx)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzz)*finner_x(1:npg,ieq))   &
                          + 0.25*cbasis(kzzx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzx)*finner_z(1:npg,ieq))
                        int_r(kxxy,ieq) =                                                                       &
                            0.25*cbasis(kxxy)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxx)*finner_y(1:npg,ieq))   &
                          + 0.25*cbasis(kxxy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxy)*finner_x(1:npg,ieq))
                    end if

                    if ( iquad > 3 ) then
                        int_r(kxxx,ieq) = 0.25*cbasis(kxxx)*dxi *                                       &
                            sum( wgt3d(1:npg)*(7.5*bfvals_int(1:npg,kx)**2 - 1.5)*finner_x(1:npg,ieq) )
                        int_r(kyyy,ieq) = 0.25*cbasis(kyyy)*dyi *                                       &
                            sum( wgt3d(1:npg)*(7.5*bfvals_int(1:npg,ky)**2 - 1.5)*finner_y(1:npg,ieq) )
                        int_r(kzzz,ieq) = 0.25*cbasis(kzzz)*dzi *                                       &
                            sum( wgt3d(1:npg)*(7.5*bfvals_int(1:npg,kz)**2 - 1.5)*finner_z(1:npg,ieq) )
                    end if

                    if ( ibitri == 1 .and. iquad > 2 ) then
                        int_r(kyyzz,ieq) =                                                                         &
                            0.25*cbasis(kyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzz)*finner_y(1:npg,ieq)) &
                          + 0.25*cbasis(kyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyyz)*finner_z(1:npg,ieq))
                        int_r(kzzxx,ieq) =                                                                         &
                            0.25*cbasis(kzzxx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxx)*finner_z(1:npg,ieq)) &
                          + 0.25*cbasis(kzzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzzx)*finner_x(1:npg,ieq))
                        int_r(kxxyy,ieq) =                                                                         &
                            0.25*cbasis(kxxyy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyy)*finner_x(1:npg,ieq)) &
                          + 0.25*cbasis(kxxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxxy)*finner_y(1:npg,ieq))
                        int_r(kyzxx,ieq) =                                                                         &
                            0.25*cbasis(kyzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_x(1:npg,ieq)) &
                          + 0.25*cbasis(kyzxx)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzxx)*finner_y(1:npg,ieq))    &
                          + 0.25*cbasis(kyzxx)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxxy)*finner_z(1:npg,ieq))
                        int_r(kzxyy,ieq) =                                                                         &
                            0.25*cbasis(kzxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_y(1:npg,ieq)) &
                          + 0.25*cbasis(kzxyy)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxyy)*finner_z(1:npg,ieq))    &
                          + 0.25*cbasis(kzxyy)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyyz)*finner_x(1:npg,ieq))
                        int_r(kxyzz,ieq) =                                                                         &
                            0.25*cbasis(kxyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyz)*finner_z(1:npg,ieq)) &
                          + 0.25*cbasis(kxyzz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyzz)*finner_x(1:npg,ieq))    &
                          + 0.25*cbasis(kxyzz)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzzx)*finner_y(1:npg,ieq))
                        int_r(kxyyzz,ieq) =                                                                          &
                            0.25*cbasis(kxyyzz)*dxi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kyyzz)*finner_x(1:npg,ieq))    &
                          + 0.25*cbasis(kxyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyzz)*finner_y(1:npg,ieq)) &
                          + 0.25*cbasis(kxyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxyy)*finner_z(1:npg,ieq))
                        int_r(kyzzxx,ieq) =                                                                          &
                            0.25*cbasis(kyzzxx)*dyi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kzzxx)*finner_y(1:npg,ieq))    &
                          + 0.25*cbasis(kyzzxx)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzxx)*finner_z(1:npg,ieq)) &
                          + 0.25*cbasis(kyzzxx)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyzz)*finner_x(1:npg,ieq))
                        int_r(kzxxyy,ieq) =                                                                          &
                            0.25*cbasis(kzxxyy)*dzi*sum(wgt3d(1:npg)*bfvals_int(1:npg,kxxyy)*finner_z(1:npg,ieq))    &
                          + 0.25*cbasis(kzxxyy)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxyy)*finner_x(1:npg,ieq)) &
                          + 0.25*cbasis(kzxxyy)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzxx)*finner_y(1:npg,ieq))
                        int_r(kxxyyzz,ieq) =                                                                           &
                            0.25*cbasis(kxxyyzz)*dxi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kxyyzz)*finner_x(1:npg,ieq)) &
                          + 0.25*cbasis(kxxyyzz)*dyi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kyzzxx)*finner_y(1:npg,ieq)) &
                          + 0.25*cbasis(kxxyyzz)*dzi*sum(wgt3d(1:npg)*3.*bfvals_int(1:npg,kzxxyy)*finner_z(1:npg,ieq))
                    end if

                end if

            end do

            do ieq = 1,nQ
                do ir=1,nbasis
                    integral_r(i,j,k,ieq,ir) = int_r(ir,ieq)
                end do
            end do

        end do
        end do
        !$OMP END PARALLEL DO
        end do

    end subroutine innerintegral
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine innerintegral_v0(Q_r)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_r
        integer i,j,k,ieq,ipg,ir
        real Qinner(npg,nQ),finner_x(npg,nQ), finner_y(npg,nQ), finner_z(npg,nQ), int_r(nbastot,nQ)

        integral_r(:,:,:,:,:) = 0.

        do k = 1,nz
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Qinner,integral_r) FIRSTPRIVATE(finner_x,finner_y,finner_z,int_r)
        do j = 1,ny
        do i = 1,nx

            do ieq = 1,nQ
                do ipg = 1,npg
                    Qinner(ipg,ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                end do
            end do

            ! NOTE: this can be made more efficient since temp vars can be re-used
            do ipg = 1,npg
                finner_x(ipg,1:nQ) = Fpt_x( Qinner(ipg,:) )
                finner_y(ipg,1:nQ) = Fpt_y( Qinner(ipg,:) )
                finner_z(ipg,1:nQ) = Fpt_z( Qinner(ipg,:) )
            end do
            ! call flux_calc_pnts_r(Qinner,finner_x,1,npg)
            ! call flux_calc_pnts_r(Qinner,finner_y,2,npg)
            ! call flux_calc_pnts_r(Qinner,finner_z,3,npg)

            do ieq = 1,nQ

                int_r(kx,ieq) = 0.25*cbasis(kx)*dxi*sum(wgt3d(1:npg)*finner_x(1:npg,ieq))
                int_r(ky,ieq) = 0.25*cbasis(ky)*dyi*sum(wgt3d(1:npg)*finner_y(1:npg,ieq))
                int_r(kz,ieq) = 0.25*cbasis(kz)*dzi*sum(wgt3d(1:npg)*finner_z(1:npg,ieq))

                if (nbasis > 4) then
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

                if (nbasis > 8) then
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
        !$OMP END PARALLEL DO
        end do

    end subroutine innerintegral_v0
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine glflux(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,ieq,ir

        ! if (llns) call get_region_Sflux_xyz(GRM_x, GRM_y, GRM_z, Sflux_x, Sflux_y, Sflux_z)

        !#########################################################
        ! Step 1: Calculate fluxes for boundaries of each cell
        !   --> flux_x, flux_y, flux_z  (used only in flux.f90)
        !---------------------------------------------------------
        ! call flux_calc(Q_r)
        call calc_flux_x(Q_r, flux_x)
        call calc_flux_y(Q_r, flux_y)
        call calc_flux_z(Q_r, flux_z)

        !#########################################################
        ! Step 2: Calc inner integral for each cell
        !   --> integral_r  (used only in flux.f90)
        !---------------------------------------------------------
        call innerintegral(Q_r)  ! call innerintegral_v0(Q_r)

        !#########################################################
        ! Step 3: Calc (total) "Gauss-Legendre flux" for each cell
        !   --> glflux_r  (used by advance_time_level_gl)
        !   --
        !---------------------------------------------------------
        do ieq = 1,nQ
        do k = 1,nz                        !FIRSTPRIVATE(flux_x,flux_y,flux_z)
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j = 1,ny
        do i = 1,nx
           glflux_r(i,j,k,ieq,1) = 0.25*(dxi*sum(wgt2d(1:nface)*(flux_x(1:nface,i+1,j,k,ieq) - flux_x(1:nface,i,j,k,ieq))) &
                                       + dyi*sum(wgt2d(1:nface)*(flux_y(1:nface,i,j+1,k,ieq) - flux_y(1:nface,i,j,k,ieq))) &
                                       + dzi*sum(wgt2d(1:nface)*(flux_z(1:nface,i,j,k+1,ieq) - flux_z(1:nface,i,j,k,ieq))))
            do ir=2,nbasis
                glflux_r(i,j,k,ieq,ir) = 0.25*cbasis(ir)*(  &
                    dxi*sum(wgt2d(1:nface)*(bfvals_xp(1:nface,ir)*flux_x(1:nface,i+1,j,k,ieq) - bfvals_xm(1:nface,ir)*flux_x(1:nface,i,j,k,ieq))) &
                  + dyi*sum(wgt2d(1:nface)*(bfvals_yp(1:nface,ir)*flux_y(1:nface,i,j+1,k,ieq) - bfvals_ym(1:nface,ir)*flux_y(1:nface,i,j,k,ieq))) &
                  + dzi*sum(wgt2d(1:nface)*(bfvals_zp(1:nface,ir)*flux_z(1:nface,i,j,k+1,ieq) - bfvals_zm(1:nface,ir)*flux_z(1:nface,i,j,k,ieq))) &
                  ) - integral_r(i,j,k,ieq,ir)
            end do

        end do
        end do
        end do
        !$OMP END PARALLEL DO
        end do

        ! NOTE: This was the original code before newCES
        ! do ieq = 1,nQ
        ! do k = 1,nz
        ! do j = 1,ny
        !   ! do ir=1,1
        !         do i = 1,nx
        !             glflux_r(i,j,k,ieq,1) =                                         &
        !                 0.25*( dxi*(wgt2d(1)*(flux_x(1,i+1,j,k,ieq) - flux_x(1,i,j,k,ieq)))  &
        !                      + dyi*(wgt2d(1)*(flux_y(1,i,j+1,k,ieq) - flux_y(1,i,j,k,ieq)))  &
        !                      + dzi*(wgt2d(1)*(flux_z(1,i,j,k+1,ieq) - flux_z(1,i,j,k,ieq)))  &
        !                      + dxi*(wgt2d(2)*(flux_x(2,i+1,j,k,ieq) - flux_x(2,i,j,k,ieq)))  &
        !                      + dyi*(wgt2d(2)*(flux_y(2,i,j+1,k,ieq) - flux_y(2,i,j,k,ieq)))  &
        !                      + dzi*(wgt2d(2)*(flux_z(2,i,j,k+1,ieq) - flux_z(2,i,j,k,ieq)))  &
        !                      + dxi*(wgt2d(3)*(flux_x(3,i+1,j,k,ieq) - flux_x(3,i,j,k,ieq)))  &
        !                      + dyi*(wgt2d(3)*(flux_y(3,i,j+1,k,ieq) - flux_y(3,i,j,k,ieq)))  &
        !                      + dzi*(wgt2d(3)*(flux_z(3,i,j,k+1,ieq) - flux_z(3,i,j,k,ieq)))  &
        !                      + dxi*(wgt2d(4)*(flux_x(4,i+1,j,k,ieq) - flux_x(4,i,j,k,ieq)))  &
        !                      + dyi*(wgt2d(4)*(flux_y(4,i,j+1,k,ieq) - flux_y(4,i,j,k,ieq)))  &
        !                      + dzi*(wgt2d(4)*(flux_z(4,i,j,k+1,ieq) - flux_z(4,i,j,k,ieq))) )
        !         end do
        !   ! end do ir=1
        !     do ir=2,nbasis
        !         do i = 1,nx
        !             glflux_r(i,j,k,ieq,ir) =                                        &
        !                  wgtbf_xmp(1,2,ir)*flux_x(1,i+1,j,k,ieq) + wgtbf_xmp(1,1,ir)*flux_x(1,i,j,k,ieq)  &
        !                + wgtbf_ymp(1,2,ir)*flux_y(1,i,j+1,k,ieq) + wgtbf_ymp(1,1,ir)*flux_y(1,i,j,k,ieq)  &
        !                + wgtbf_zmp(1,2,ir)*flux_z(1,i,j,k+1,ieq) + wgtbf_zmp(1,1,ir)*flux_z(1,i,j,k,ieq)  &
        !                + wgtbf_xmp(2,2,ir)*flux_x(2,i+1,j,k,ieq) + wgtbf_xmp(2,1,ir)*flux_x(2,i,j,k,ieq)  &
        !                + wgtbf_ymp(2,2,ir)*flux_y(2,i,j+1,k,ieq) + wgtbf_ymp(2,1,ir)*flux_y(2,i,j,k,ieq)  &
        !                + wgtbf_zmp(2,2,ir)*flux_z(2,i,j,k+1,ieq) + wgtbf_zmp(2,1,ir)*flux_z(2,i,j,k,ieq)  &
        !                + wgtbf_xmp(3,2,ir)*flux_x(3,i+1,j,k,ieq) + wgtbf_xmp(3,1,ir)*flux_x(3,i,j,k,ieq)  &
        !                + wgtbf_ymp(3,2,ir)*flux_y(3,i,j+1,k,ieq) + wgtbf_ymp(3,1,ir)*flux_y(3,i,j,k,ieq)  &
        !                + wgtbf_zmp(3,2,ir)*flux_z(3,i,j,k+1,ieq) + wgtbf_zmp(3,1,ir)*flux_z(3,i,j,k,ieq)  &
        !                + wgtbf_xmp(4,2,ir)*flux_x(4,i+1,j,k,ieq) + wgtbf_xmp(4,1,ir)*flux_x(4,i,j,k,ieq)  &
        !                + wgtbf_ymp(4,2,ir)*flux_y(4,i,j+1,k,ieq) + wgtbf_ymp(4,1,ir)*flux_y(4,i,j,k,ieq)  &
        !                + wgtbf_zmp(4,2,ir)*flux_z(4,i,j,k+1,ieq) + wgtbf_zmp(4,1,ir)*flux_z(4,i,j,k,ieq)  &
        !                - integral_r(i,j,k,ieq,ir)
        !         end do
        !     end do
        ! end do
        ! end do
        ! end do
    end subroutine glflux
!-------------------------------------------------------------------------------


!-----------------------------------------------------------!
!*******************calculate freezing speeds***************!
!-----------------------------------------------------------!
    real function cfcal(Qcf, cases)
        implicit none
        real, dimension(nQ), intent(in) :: Qcf
        integer, intent(in) :: cases

        real dn,dni, vx,vy,vz, P, cs

        dn = Qcf(rh)
        dni = 1./dn
        vx = Qcf(mx)*dni
        vy = Qcf(my)*dni
        vz = Qcf(mz)*dni

        !--- pressure -----------------
        select case(ivis)
            case(0)
                P = aindm1*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
            case(1)
                P = aindm1*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
            case(2)
                P = c1d3 * ( Qcf(exx) + Qcf(eyy) + Qcf(ezz) - dn*(vx**2 + vy**2 + vz**2) )
        end select

        !--- sound speed --------------
        select case(ieos)
            case(1)
                cs = sqrt(aindex*P*dni)
            case(2)
                cs = sqrt(7.2*P_1*dn**6.2)
        end select

        !--- freezing speed -----------
        select case(cases)
            case(1) !freezing speed in x direction for fluid variable
                cfcal = abs(vx) + cs

            case(2) !freezing speed in y direction for fluid variable
                cfcal = abs(vy) + cs

            case(3) !freezing speed in z direction for fluid variable
                cfcal = abs(vz) + cs
        end select

    end function cfcal
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
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
                    if (theta > 1.) then
                        theta = 1.
                    end if

                    if (theta < 0) then
                        theta = 0.
                    end if
                    do ir=2,nbasis
                        Q_r(i,j,k,rh,ir) = abs(theta)*Q_r(i,j,k,rh,ir)
                    end do

                end if

                Pave = aindm1*(Q_r(i,j,k,en,1)                                  &
                        - 0.5*( Q_r(i,j,k,mx,1)**2                              &
                              + Q_r(i,j,k,my,1)**2                              &
                              + Q_r(i,j,k,mz,1)**2 ) / Q_r(i,j,k,rh,1) )

                if (Pave < epsiP) then
                    do ir=2,nbasis
                        Q_r(i,j,k,rh:en,ir) = 0.0
                    end do
                else
                    theta = 1.
                    do ipge = 1,npge
                        do ieq = rh,en
                            Qedge(ipge,ieq) = sum(bf_faces(ipge,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
                        end do

                        dn = Qedge(ipge,rh)
                        dni = 1./dn
                        P(ipge) = aindm1*(Qedge(ipge,en)                        &
                                - 0.5*dni*( Qedge(ipge,mx)**2                   &
                                          + Qedge(ipge,my)**2                   &
                                          + Qedge(ipge,mz)**2) )

                        if (P(ipge) < epsiP) then
                            if (Pave .ne. P(ipge)) then
                                thetaj = (Pave - epsiP)/(Pave - P(ipge))
                                theta = min(theta,thetaj)
                            end if
                        end if
                    end do
                    if (theta > 1.) then
                        theta = 1.
                    end if

                    if (theta < 0.) then
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
!-------------------------------------------------------------------------------

end module flux
