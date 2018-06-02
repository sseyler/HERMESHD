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

        subroutine Fpts_ptr(Qpts, Fpts, i,j,k, npts)
            use params, only : nQ
            integer, intent(in) :: i,j,k,npts
            real, dimension(npts,nQ), intent(in)  :: Qpts
            real, dimension(npts,nQ), intent(out) :: Fpts
        end subroutine Fpts_ptr
    end interface
    !---------------------------------------------------------------------------

    !===========================================================================
    ! Initialize pointer to temporal integration subroutine
    !-----------------------------------------------------------------
    procedure (Fpt_ptr), pointer :: Fpt_x => null ()
    procedure (Fpt_ptr), pointer :: Fpt_y => null ()
    procedure (Fpt_ptr), pointer :: Fpt_z => null ()

    procedure (Fpts_ptr), pointer :: Fpts_x => null ()
    procedure (Fpts_ptr), pointer :: Fpts_y => null ()
    procedure (Fpts_ptr), pointer :: Fpts_z => null ()
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Global variables
    !-----------------------------------------------------------------
    real, dimension(nx,ny,nz,nQ,nbasis) :: glflux_r, integral_r

    ! Only used by flux_calc (flux_cal) and glflux
    real, dimension(nface,1:nx1,ny,nz,1:nQ) :: flux_x
    real, dimension(nface,nx,1:ny1,nz,1:nQ) :: flux_y
    real, dimension(nface,nx,ny,1:nz1,1:nQ) :: flux_z
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
    ! DEBUG: See the call to this subroutine in main (currently handling two situations)
    subroutine select_hydro_model(ivis, flux_x_pt, flux_y_pt, flux_z_pt, flux_x_pts, flux_y_pts, flux_z_pts)
        implicit none
        integer, intent(in) :: ivis
        procedure(Fpt_ptr),  pointer :: flux_x_pt, flux_y_pt, flux_z_pt
        procedure(Fpts_ptr), pointer :: flux_x_pts,flux_y_pts,flux_z_pts

        select case (ivis)
            case (0)
                call mpi_print(iam, 'ivis = 0: Selected Euler equations (inviscid fluid)')
                ! flux_x_pts => Fpts_x_0 ! DEBUG
                ! flux_y_pts => Fpts_y_0 ! DEBUG
                ! flux_z_pts => Fpts_z_0 ! DEBUG
            case (1)
                call mpi_print(iam, 'ivis = 1: Selected linearized 10-moment system')
                ! flux_x_pts => Fpts_x_1 ! DEBUG
                ! flux_y_pts => Fpts_y_1 ! DEBUG
                ! flux_z_pts => Fpts_z_1 ! DEBUG
                flux_x_pts => FSpts_x_1 ! DEBUG
                flux_y_pts => FSpts_y_1 ! DEBUG
                flux_z_pts => FSpts_z_1 ! DEBUG
            case (2)
                ! WARNING/TODO
                call mpi_print(iam, 'ivis = 2: Selected full (nonlinear) 10-moment system')
                flux_x_pts => Fpts_x_2 ! DEBUG
                flux_y_pts => Fpts_y_2 ! DEBUG
                flux_z_pts => Fpts_z_2 ! DEBUG
                ! call exit(-1)
            case default
                call mpi_print(iam, 'ivis not set: Defaulting to Euler equations (inviscid fluid)')
                ! flux_x_pt => Fpt_x_0
                ! flux_y_pt => Fpt_y_0
                ! flux_z_pt => Fpt_z_0
                ! flux_x_pts => Fpts_x_0 ! DEBUG
                ! flux_y_pts => Fpts_y_0 ! DEBUG
                ! flux_z_pts => Fpts_z_0 ! DEBUG
        end select
    end subroutine select_hydro_model
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Linearized 10-moment fluxes
    !---------------------------------------------------------------------------
    subroutine FSpts_x_1(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k, npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,m_x,m_y,m_z,Ener
        real :: dni, vx,vy,vz, E_xx,E_xy,E_xz, P,c2d3_Pvx
        real :: S_xx,S_xy,S_xz, EmS_xx,EmS_xy,EmS_xz

        S_xx = Sflux(npts_llns,i,j,k,1,1)
        S_xy = Sflux(npts_llns,i,j,k,1,2)
        S_xz = Sflux(npts_llns,i,j,k,1,3)

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            m_x  = Qpts(ife,mx)
            m_y  = Qpts(ife,my)
            m_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xx = Qpts(ife,exx)
            E_xy = Qpts(ife,exy)
            E_xz = Qpts(ife,exz)
            EmS_xx = E_xx - S_xx
            EmS_xy = E_xy - S_xy
            EmS_xz = E_xz - S_xz
            !-----------------
            dni = 1./dn
            vx  = m_x*dni
            vy  = m_y*dni
            vz  = m_z*dni
            v2  = vx*vx + vy*vy + vz*vz
            P   = aindm1 * ( Ener - 0.5*dn*v2 )
            if (P < P_floor) P = P_floor
            c2d3_Pvx = c2d3*P*vx
            !------------------------------------
            Fpts(ife,rh)  = m_x
            Fpts(ife,mx)  = m_x*vx + EmS_xx + P !+ c1d3*dn*v2
            Fpts(ife,my)  = m_x*vy + EmS_xy
            Fpts(ife,mz)  = m_x*vz + EmS_xz
            Fpts(ife,en)  = (Ener + P + EmS_xx)*vx + EmS_xy*vy + EmS_xz*vz
            Fpts(ife,exx) = 2*c2d3_Pvx  !  c4d3ne*vx
            Fpts(ife,eyy) = - c2d3_Pvx  ! -c2d3ne*vx
            Fpts(ife,ezz) = - c2d3_Pvx  ! -c2d3ne*vx
            Fpts(ife,exy) =       P*vy  !  nueta*vy
            Fpts(ife,exz) =       P*vz  !  nueta*vz
            Fpts(ife,eyz) =        0    !  0
        end do
    end subroutine FSpts_x_1
    !---------------------------------------------------------------------------
    subroutine FSpts_y_1(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k, npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,M_x,M_y,M_z,Ener
        real :: dni, vx,vy,vz, E_xy,E_yy,E_yz, P,c2d3_Pvy
        real :: S_xy,S_yy,S_yz, EmS_xy,EmS_yy,EmS_yz

        S_xy = Sflux(npts_llns,i,j,k,1,2)
        S_yy = Sflux(npts_llns,i,j,k,2,2)
        S_yz = Sflux(npts_llns,i,j,k,2,3)

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            m_x  = Qpts(ife,mx)
            m_y  = Qpts(ife,my)
            m_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xy = Qpts(ife,exy)
            E_yy = Qpts(ife,eyy)
            E_yz = Qpts(ife,eyz)
            EmS_xy = E_xy - S_xy
            EmS_yy = E_yy - S_yy
            EmS_yz = E_yz - S_yz
            !-----------------
            dni = 1./dn
            vx  = m_x*dni
            vy  = m_y*dni
            vz  = m_z*dni
            v2  = vx*vx + vy*vy + vz*vz
            P   = aindm1 * ( Ener - 0.5*dn*v2 )
            if (P < P_floor) P = P_floor
            c2d3_Pvy = c2d3*P*vy
            !------------------------------------
            Fpts(ife,rh)  = m_y
            Fpts(ife,mx)  = m_y*vx + EmS_xy
            Fpts(ife,my)  = m_y*vy + EmS_yy + P !+ c1d3*dn*v2
            Fpts(ife,mz)  = m_y*vz + EmS_yz
            Fpts(ife,en)  = (Ener + P + EmS_yy)*vy + EmS_yz*vz + EmS_xy*vx
            Fpts(ife,exx) = - c2d3_Pvy  ! -c2d3ne*vy
            Fpts(ife,eyy) = 2*c2d3_Pvy  !  c4d3ne*vy
            Fpts(ife,ezz) = - c2d3_Pvy  ! -c2d3ne*vy
            Fpts(ife,exy) =       P*vx  !  nueta*vx
            Fpts(ife,exz) =        0    !  0
            Fpts(ife,eyz) =       P*vz  !  nueta*vz
        end do
    end subroutine FSpts_y_1
    !---------------------------------------------------------------------------
    subroutine FSpts_z_1(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k, npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,M_x,M_y,M_z,Ener
        real :: dni, vx,vy,vz, E_xz,E_yz,E_zz, P,c2d3_Pvz
        real :: S_xz,S_yz,S_zz, EmS_xz,EmS_yz,EmS_zz

        S_xz = Sflux(npts_llns,i,j,k,1,3)
        S_yz = Sflux(npts_llns,i,j,k,2,3)
        S_zz = Sflux(npts_llns,i,j,k,3,3)

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            m_x  = Qpts(ife,mx)
            m_y  = Qpts(ife,my)
            m_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xz = Qpts(ife,exz)
            E_yz = Qpts(ife,eyz)
            E_zz = Qpts(ife,ezz)
            EmS_xz = E_xz - S_xz
            EmS_yz = E_yz - S_yz
            EmS_zz = E_zz - S_zz
            !-----------------
            dni = 1./dn
            vx  = m_x*dni
            vy  = m_y*dni
            vz  = m_z*dni
            v2  = vx*vx + vy*vy + vz*vz
            P   = aindm1 * ( Ener - 0.5*dn*v2 )
            if (P < P_floor) P = P_floor
            c2d3_Pvz = c2d3*P*vz
            !------------------------------------
            Fpts(ife,rh)  = m_z
            Fpts(ife,mx)  = m_z*vx + EmS_xz
            Fpts(ife,my)  = m_z*vy + EmS_yz
            Fpts(ife,mz)  = m_z*vz + EmS_zz + P !+ c1d3*dn*v2
            Fpts(ife,en)  = (Ener + P + Ems_zz)*vz + Ems_xz*vx + Ems_yz*vy
            Fpts(ife,exx) = - c2d3_Pvz  ! -c2d3ne*vz
            Fpts(ife,eyy) = - c2d3_Pvz  ! -c2d3ne*vz
            Fpts(ife,ezz) = 2*c2d3_Pvz  !  c4d3ne*vz
            Fpts(ife,exy) =        0    !  0
            Fpts(ife,exz) =       P*vx  !  nueta*vx
            Fpts(ife,eyz) =       P*vy  !  nueta*vy
        end do
    end subroutine FSpts_z_1
    !---------------------------------------------------------------------------

    !===========================================================================
    ! Linearized 10-moment fluxes
    !---------------------------------------------------------------------------
    subroutine Fpts_x_1(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,m_x,m_y,m_z,Ener
        real :: dni, vx,vy,vz, E_xx,E_xy,E_xz, P_xx,P_xy,P_xz, P,c2d3_Pvx

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            m_x  = Qpts(ife,mx)
            m_y  = Qpts(ife,my)
            m_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xx = Qpts(ife,exx)
            E_xy = Qpts(ife,exy)
            E_xz = Qpts(ife,exz)
            !-----------------
            dni = 1./dn
            vx  = m_x*dni
            vy  = m_y*dni
            vz  = m_z*dni
            P   = aindm1 * ( Ener - 0.5*dn*(vx*vx + vy*vy + vz*vz) )
            if (P < P_floor) P = P_floor
            c2d3_Pvx = c2d3*P*vx
            !------------------------------------
            Fpts(ife,rh)  = m_x
            Fpts(ife,mx)  = m_x*vx + E_xx + P
            Fpts(ife,my)  = m_x*vy + E_xy
            Fpts(ife,mz)  = m_x*vz + E_xz
            Fpts(ife,en)  = (Ener + P + E_xx)*vx + E_xy*vy + E_xz*vz
            Fpts(ife,exx) = 2*c2d3_Pvx  !  c4d3ne*vx
            Fpts(ife,eyy) = - c2d3_Pvx  ! -c2d3ne*vx
            Fpts(ife,ezz) = - c2d3_Pvx  ! -c2d3ne*vx
            Fpts(ife,exy) =       P*vy  !  nueta*vy
            Fpts(ife,exz) =       P*vz  !  nueta*vz
            Fpts(ife,eyz) =        0    !  0
        end do
    end subroutine Fpts_x_1
    !---------------------------------------------------------------------------
    subroutine Fpts_y_1(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,M_x,M_y,M_z,Ener
        real :: dni, vx,vy,vz, E_xy,E_yy,E_yz, P_xy,P_yy,P_yz, P,c2d3_Pvy

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            m_x  = Qpts(ife,mx)
            m_y  = Qpts(ife,my)
            m_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xy = Qpts(ife,exy)
            E_yy = Qpts(ife,eyy)
            E_yz = Qpts(ife,eyz)
            !-----------------
            dni = 1./dn
            vx  = m_x*dni
            vy  = m_y*dni
            vz  = m_z*dni
            P   = aindm1 * ( Ener - 0.5*dn*(vx*vx + vy*vy + vz*vz) )
            if (P < P_floor) P = P_floor
            c2d3_Pvy = c2d3*P*vy
            !------------------------------------
            Fpts(ife,rh)  = m_y
            Fpts(ife,mx)  = m_y*vx + E_xy
            Fpts(ife,my)  = m_y*vy + E_yy + P
            Fpts(ife,mz)  = m_y*vz + E_yz
            Fpts(ife,en)  = (Ener + P + E_yy)*vy + E_yz*vz + E_xy*vx
            Fpts(ife,exx) = - c2d3_Pvy  ! -c2d3ne*vy
            Fpts(ife,eyy) = 2*c2d3_Pvy  !  c4d3ne*vy
            Fpts(ife,ezz) = - c2d3_Pvy  ! -c2d3ne*vy
            Fpts(ife,exy) =       P*vx  !  nueta*vx
            Fpts(ife,exz) =        0    !  0
            Fpts(ife,eyz) =       P*vz  !  nueta*vz
        end do
    end subroutine Fpts_y_1
    !---------------------------------------------------------------------------
    subroutine Fpts_z_1(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,M_x,M_y,M_z,Ener
        real :: dni, vx,vy,vz, E_xz,E_yz,E_zz, P_xz,P_yz,P_zz, P,c2d3_Pvz

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            m_x  = Qpts(ife,mx)
            m_y  = Qpts(ife,my)
            m_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xz = Qpts(ife,exz)
            E_yz = Qpts(ife,eyz)
            E_zz = Qpts(ife,ezz)
            !-----------------
            dni = 1./dn
            vx  = m_x*dni
            vy  = m_y*dni
            vz  = m_z*dni
            P   = aindm1 * ( Ener - 0.5*dn*(vx*vx + vy*vy + vz*vz) )
            if (P < P_floor) P = P_floor
            c2d3_Pvz = c2d3*P*vz
            !------------------------------------
            Fpts(ife,rh)  = m_z
            Fpts(ife,mx)  = m_z*vx + E_xz
            Fpts(ife,my)  = m_z*vy + E_yz
            Fpts(ife,mz)  = m_z*vz + E_zz + P
            Fpts(ife,en)  = (Ener + P + E_zz)*vz + E_xz*vx + E_yz*vy
            Fpts(ife,exx) = - c2d3_Pvz  ! -c2d3ne*vz
            Fpts(ife,eyy) = - c2d3_Pvz  ! -c2d3ne*vz
            Fpts(ife,ezz) = 2*c2d3_Pvz  !  c4d3ne*vz
            Fpts(ife,exy) =        0    !  0
            Fpts(ife,exz) =       P*vx  !  nueta*vx
            Fpts(ife,eyz) =       P*vy  !  nueta*vy
        end do
    end subroutine Fpts_z_1
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    subroutine Spts_x_1(Qpts, Spts, i, j, k, npts)
        ! DEBUG
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Spts
        real :: dn,m_x,m_y,m_z, dni,vx,vy,vz, S_xx,S_xy,S_xz

        S_xx = Sflux_x(npts_llns,i,j,k,1,1)
        S_xy = Sflux_x(npts_llns,i,j,k,1,2)
        S_xz = Sflux_x(npts_llns,i,j,k,1,3)

        do ife = 1,npts
            dn  = Qpts(ife,rh)
            m_x = Qpts(ife,mx)
            m_y = Qpts(ife,my)
            m_z = Qpts(ife,mz)
            dni = 1./dn
            vx = m_x*dni
            vy = m_y*dni
            vz = m_z*dni

            Spts(ife,rh) = 0
            Spts(ife,mx) = S_xx
            Spts(ife,my) = S_xy
            Spts(ife,mz) = S_xz
            Spts(ife,en) = S_xx*vx + S_xy*vy + S_xz*vz
        end do
    end subroutine Spts_x_1
    !---------------------------------------------------------------------------
    subroutine Spts_y_1(Qpts, Spts, i, j, k, npts)
        ! DEBUG
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Spts
        real :: dn,m_x,m_y,m_z, dni,vx,vy,vz, S_xy,S_yy,S_yz

        S_xy = Sflux_y(npts_llns,i,j,k,1,2)
        S_yy = Sflux_y(npts_llns,i,j,k,2,2)
        S_yz = Sflux_y(npts_llns,i,j,k,2,3)

        do ife = 1,npts
            dn  = Qpts(ife,rh)
            m_x = Qpts(ife,mx)
            m_y = Qpts(ife,my)
            m_z = Qpts(ife,mz)
            dni = 1./dn
            vx = m_x*dni
            vy = m_y*dni
            vz = m_z*dni

            Spts(ife,rh) = 0
            Spts(ife,mx) = S_xy
            Spts(ife,my) = S_yy
            Spts(ife,mz) = S_yz
            Spts(ife,en) = S_yy*vy + S_yz*vz + S_xy*vx
        end do
    end subroutine Spts_y_1
    !---------------------------------------------------------------------------
    subroutine Spts_z_1(Qpts, Spts, i, j, k, npts)
        ! DEBUG
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Spts
        real :: dn,m_x,m_y,m_z, dni,vx,vy,vz ,S_xz,S_yz,S_zz

        S_xz = Sflux_z(npts_llns,i,j,k,1,3)
        S_yz = Sflux_z(npts_llns,i,j,k,2,3)
        S_zz = Sflux_z(npts_llns,i,j,k,3,3)

        do ife = 1,npts
            dn  = Qpts(ife,rh)
            m_x = Qpts(ife,mx)
            m_y = Qpts(ife,my)
            m_z = Qpts(ife,mz)
            dni = 1./dn
            vx = m_x*dni
            vy = m_y*dni
            vz = m_z*dni

            Spts(ife,rh) = 0
            Spts(ife,mx) = S_xz
            Spts(ife,my) = S_yz
            Spts(ife,mz) = S_zz
            Spts(ife,en) = S_zz*vz + S_xz*vx + S_yz*vy
        end do
    end subroutine Spts_z_1
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Nonlinear 10 moment-eqns fluxes
    !---------------------------------------------------------------------------
    subroutine Fpts_x_2(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,M_x,M_y,M_z,Ener
        real :: dni, vx,vy,vz, E_xx,E_xy,E_xz, P
        real :: c2d3_Pvx, vx2,vy2,vz2,vsq,vsqd3, Mxvx,Mxvy,Mxvz

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            M_x  = Qpts(ife,mx)
            M_y  = Qpts(ife,my)
            M_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xx = Qpts(ife,exx)
            E_xy = Qpts(ife,exy)
            E_xz = Qpts(ife,exz)
            !-----------------
            dni = 1./dn
            vx  = M_x*dni
            vy  = M_y*dni
            vz  = M_z*dni
            vx2 = vx*vx
            vy2 = vy*vy
            vz2 = vz*vz
            vsq = vx2 + vy2 + vz2
            vsqd3 = c1d3*vsq
            Mxvx = M_x*vx
            Mxvy = M_x*vy
            Mxvz = M_x*vz
            P    = aindm1 * ( Ener - 0.5*dn*vsq )
            if (P < P_floor) P = P_floor
            c2d3_Pvx = c2d3*P*vx
            !------------------------------------
            Fpts(ife,rh)  = M_x
            Fpts(ife,mx)  = Mxvx + E_xx + P
            Fpts(ife,my)  = Mxvy + E_xy
            Fpts(ife,mz)  = Mxvz + E_xz
            Fpts(ife,en)  = (Ener + P)*vx
            Fpts(ife,exx) = M_x*(vx2 - vsqd3) + 2*c2d3_Pvx
            Fpts(ife,eyy) = M_x*(vy2 - vsqd3) -   c2d3_Pvx
            Fpts(ife,ezz) = M_x*(vz2 - vsqd3) -   c2d3_Pvx
            Fpts(ife,exy) = Mxvx*vy + P*vy
            Fpts(ife,exz) = Mxvz*vx + P*vz
            Fpts(ife,eyz) = Mxvy*vz
        end do
    end subroutine Fpts_x_2
    !---------------------------------------------------------------------------
    subroutine Fpts_y_2(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,M_x,M_y,M_z,Ener
        real :: dni, vx,vy,vz, E_xy,E_yy,E_yz, P
        real :: c2d3_Pvy, vx2,vy2,vz2,vsq,vsqd3, Myvx,Myvy,Myvz

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            M_x  = Qpts(ife,mx)
            M_y  = Qpts(ife,my)
            M_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xy = Qpts(ife,exy)
            E_yy = Qpts(ife,eyy)
            E_yz = Qpts(ife,eyz)
            !-----------------
            dni = 1./dn
            vx  = M_x*dni
            vy  = M_y*dni
            vz  = M_z*dni
            vx2 = vx*vx
            vy2 = vy*vy
            vz2 = vz*vz
            vsq = vx2 + vy2 + vz2
            vsqd3 = c1d3*vsq
            Myvx = M_y*vx
            Myvy = M_y*vy
            Myvz = M_y*vz
            P    = aindm1 * ( Ener - 0.5*dn*vsq )
            if (P < P_floor) P = P_floor
            c2d3_Pvy = c2d3*P*vy
            !------------------------------------
            Fpts(ife,rh)  = M_y
            Fpts(ife,mx)  = Myvx + E_xy
            Fpts(ife,my)  = Myvy + E_yy + P
            Fpts(ife,mz)  = Myvz + E_yz
            Fpts(ife,en)  = (Ener + P)*vy
            Fpts(ife,exx) = M_y*(vx2 - vsqd3) -   c2d3_Pvy
            Fpts(ife,eyy) = M_y*(vy2 - vsqd3) + 2*c2d3_Pvy
            Fpts(ife,ezz) = M_y*(vz2 - vsqd3) -   c2d3_Pvy
            Fpts(ife,exy) = Myvy*vx + P*vx
            Fpts(ife,exz) = Myvx*vz
            Fpts(ife,eyz) = Myvz*vy + P*vz
        end do
    end subroutine Fpts_y_2
    !---------------------------------------------------------------------------
    subroutine Fpts_z_2(Qpts, Fpts, i,j,k, npts)
        ! WARNING: the stochastic terms have been excluded here!!!
        integer, intent(in) :: i,j,k,npts
        real, dimension(npts,nQ), intent(in)  :: Qpts
        real, dimension(npts,nQ), intent(out) :: Fpts
        real :: dn,M_x,M_y,M_z,Ener
        real :: dni, vx,vy,vz, E_xz,E_yz,E_zz, P
        real :: c2d3_Pvz, vx2,vy2,vz2,vsq,vsqd3, Mzvx,Mzvy,Mzvz

        do ife = 1,npts
            dn   = Qpts(ife,rh)
            M_x  = Qpts(ife,mx)
            M_y  = Qpts(ife,my)
            M_z  = Qpts(ife,mz)
            Ener = Qpts(ife,en)
            E_xz = Qpts(ife,exz)
            E_yz = Qpts(ife,eyz)
            E_zz = Qpts(ife,ezz)
            !-----------------
            dni = 1./dn
            vx  = M_x*dni
            vy  = M_y*dni
            vz  = M_z*dni
            vx2 = vx*vx
            vy2 = vy*vy
            vz2 = vz*vz
            vsq = vx2 + vy2 + vz2
            vsqd3 = c1d3*vsq
            Myvx = M_y*vx
            Myvy = M_y*vy
            Myvz = M_y*vz
            P    = aindm1 * ( Ener - 0.5*dn*vsq )
            if (P < P_floor) P = P_floor
            c2d3_Pvz = c2d3*P*vz
            !------------------------------------
            Fpts(ife,rh)  = M_z
            Fpts(ife,mx)  = Mzvx + E_xz
            Fpts(ife,my)  = Mzvy + E_yz
            Fpts(ife,mz)  = Mzvz + E_zz + P
            Fpts(ife,en)  = (Ener + P)*vz
            Fpts(ife,exx) = M_z*(vx2 - vsqd3) -   c2d3_Pvz
            Fpts(ife,eyy) = M_z*(vy2 - vsqd3) -   c2d3_Pvz
            Fpts(ife,ezz) = M_z*(vz2 - vsqd3) + 2*c2d3_Pvz
            Fpts(ife,exy) = Mzvy*vx
            Fpts(ife,exz) = Mzvx*vz + P*vx
            Fpts(ife,eyz) = Mzvz*vy + P*vy
        end do
    end subroutine Fpts_z_2
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Get field var values at a cell face in body (not boundary) of spatial region
    !---------------------------------------------------------------------------
    function Qface_body(Q_ijk, bfv_r)
        real, dimension(nQ,nbasis), intent(in)    :: Q_ijk
        real, dimension(nface,nbasis), intent(in) :: bfv_r
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
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in)  :: Q_r
        real, dimension(nface,nx1,ny,nz,nQ), intent(out) :: flux_x
        real, dimension(nfe,nQ) :: Qface_x,Fface_x,Fface_x1
        real, dimension(nfe,3,3) :: Sface_x
        real, dimension(nface,nQ) :: cfrx
        integer :: i,j,k,ieq,im1,i4,i4p,ipnt,ife
        real :: cwavex(nfe),fhllc_x(nface,5),qvin(nQ)

        ! Sface_x(:,:) = 0

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(im1,Qface_x,qvin) FIRSTPRIVATE(Fface_x,cfrx,cwavex,fhllc_x)
        do k=1,nz
        do j=1,ny
        do i=1,nx1
            im1 = i-1
            !--- get Qface_x --------------------------
            if (i > 1) then
                Qface_x(1:nface,1:nQ) = Qface_body( Q_r(im1,j,k,1:nQ,:),bfvals_xp(1:nface,:) )
            else if (i == 1) then
                Qface_x(1:nface,1:nQ) = Qface_edge( Qxlo_ext(j,k,1:nface,1:nQ) )
            end if
            if (i < nx1) then
                Qface_x(nface1:nfe,1:nQ) = Qface_body( Q_r(i,j,k,1:nQ,:),bfvals_xm(1:nface,:) )
            else if (i == nx1) then
                Qface_x(nface1:nfe,1:nQ) = Qface_edge( Qxhi_ext(j,k,1:nface,1:nQ) )
            end if
            !--- get Sface_x --------------------------
            ! if (llns) then
            !     do ipnt = 1,nfe
            !         Sface_x(ipnt,1:3,1:3) = Sflux_x(npts_llns,i,j,k,1:3,1:3)
            !     end do
            ! end if
            !--- get Fface_x --------------------------
            call FSpts_x_1(Qface_x, Fface_x, i,j,k, nfe)
            ! call Fpts_x(Qface_x, Fface_x, nfe)                    ! WARNING: UNDO
            ! if (llns) call Spts_x_1(Qface_x, Sface_x, i,j,k, nfe) ! WARNING: UNDO
            ! Fface_x = Fface_x - Sface_x                           ! WARNING: UNDO
            !--- get cfrx -----------------------------
            do i4=1,nfe
                do ieq=1,nQ
                    qvin(ieq) = Qface_x(i4,ieq)
                end do
                cwavex(i4) = cfcal(qvin,1)
            end do
            do i4=1,nface
                cfrx(i4,1:nQ) = max(cwavex(i4),cwavex(i4+nface))
            end do
            !--- get flux_x ---------------------------
            do ieq = 1,nQ
                do i4=1,nface
                    i4p = i4 + nface
                    flux_x(i4,i,j,k,ieq) = 0.5*( (Fface_x(i4,ieq) + Fface_x(i4p,ieq)) &
                                         - cfrx(i4,ieq)*(Qface_x(i4p,ieq) - Qface_x(i4,ieq)) )
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
            !------------------------------------------
        end do
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine calc_flux_x
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    subroutine calc_flux_y(Q_r, flux_y)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in)  :: Q_r
        real, dimension(nface,nx,ny1,nz,nQ), intent(out) :: flux_y
        real, dimension(nfe,nQ) :: Qface_y,Fface_y,Sface_y,Fface_y1
        real, dimension(nface,nQ) :: cfry
        integer :: i,j,k,ieq,jm1,i4,i4p,ipnt,ife
        real :: cwavey(nfe),fhllc_y(nface,5),qvin(nQ)

        Sface_y(:,:) = 0

        do k=1,nz
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jm1,Qface_y,qvin) FIRSTPRIVATE(Fface_y,cfry,cwavey,fhllc_y)
        do j=1,ny1
        jm1 = j-1
        do i=1,nx
            !--- get Qface_y --------------------------
            if (j > 1) then
                Qface_y(1:nface,1:nQ)  = Qface_body( Q_r(i,jm1,k,1:nQ,:),bfvals_yp(1:nface,:) )
            else if (j == 1) then
                Qface_y(1:nface,1:nQ)  = Qface_edge( Qylo_ext(i,k,1:nface,:) )
            end if
            if (j < ny1) then
                Qface_y(nface1:nfe,1:nQ) = Qface_body( Q_r(i,j,k,1:nQ,:),bfvals_ym(1:nface,:) )
            else if (j == ny1) then
                Qface_y(nface1:nfe,1:nQ) = Qface_edge( Qyhi_ext(i,k,1:nface,:) )
            end if
            !--- get Fface_y --------------------------
            call FSpts_y_1(Qface_y, Fface_y, i,j,k, nfe)
            ! call Fpts_y(Qface_y, Fface_y, nfe)                    ! WARNING: UNDO
            ! if (llns) call Spts_y_1(Qface_y, Sface_y, i,j,k, nfe) ! WARNING: UNDO
            ! Fface_y = Fface_y - Sface_y                           ! WARNING: UNDO
            !--- get cfry -----------------------------
            do i4=1,nfe
                do ieq=1,nQ
                    qvin(ieq) = Qface_y(i4,ieq)
                end do
                cwavey(i4) = cfcal(qvin,2)
            end do
            do i4=1,nface
                cfry(i4,1:nQ) = max(cwavey(i4),cwavey(i4+nface))
            end do
            !--- get flux_y ---------------------------
            do ieq = 1,nQ
                do i4=1,nface
                    i4p = i4 + nface
                    flux_y(i4,i,j,k,ieq) = 0.5*( (Fface_y(i4,ieq) + Fface_y(i4p,ieq)) &
                                         - cfry(i4,ieq)*(Qface_y(i4p,ieq) - Qface_y(i4,ieq)) )
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
            !------------------------------------------
        end do
        end do
        !$OMP END PARALLEL DO
        end do
    end subroutine calc_flux_y
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    subroutine calc_flux_z(Q_r, flux_z)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in)  :: Q_r
        real, dimension(nface,nx,ny,nz1,nQ), intent(out) :: flux_z
        real, dimension(nfe,nQ) :: Qface_z,Fface_z,Sface_z,Fface_z1
        real, dimension(nface,nQ) :: cfrz
        integer :: i,j,k,ieq,km1,i4,i4p,ipnt,ife
        real :: cwavez(nfe),fhllc_z(nface,5),qvin(nQ)

        Sface_z(:,:) = 0

        do k=1,nz1
        km1 = k-1
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Qface_z,qvin) FIRSTPRIVATE(Fface_z,cfrz,cwavez,fhllc_z)
        do j=1,ny
        do i=1,nx
            !--- get Qface_z --------------------------
            if (k > 1) then
                Qface_z(1:nface,1:nQ)  = Qface_body( Q_r(i,j,km1,1:nQ,:),bfvals_zp(1:nface,:) )
            else if (k == 1) then
                Qface_z(1:nface,1:nQ)  = Qface_edge( Qzlo_ext(i,j,1:nface,:) )
            end if
            if (k < nz1) then
                Qface_z(nface1:nfe,1:nQ) = Qface_body( Q_r(i,j,k,1:nQ,:),bfvals_zm(1:nface,:) )
            else if (k == nz1) then
                Qface_z(nface1:nfe,1:nQ) = Qface_edge( Qzhi_ext(i,j,1:nface,:) )
            end if
            !--- get Fface_z --------------------------
            call FSpts_z_1(Qface_z, Fface_z, i,j,k, nfe)
            ! call Fpts_z(Qface_z, Fface_z, nfe)                    ! WARNING: UNDO
            ! if (llns) call Spts_z_1(Qface_z, Sface_z, i,j,k, nfe) ! WARNING: UNDO
            ! Fface_z = Fface_z - Sface_z                           ! WARNING: UNDO
            !--- get cfrz -----------------------------
            do i4=1,nfe
                do ieq=1,nQ
                    qvin(ieq) = Qface_z(i4,ieq)
                end do
                cwavez(i4) = cfcal(qvin,3)
            end do
            do i4=1,nface
                cfrz(i4,1:nQ) = max(cwavez(i4),cwavez(i4+nface))
            end do
            !--- get flux_z ---------------------------
            do ieq = 1,nQ
                do i4=1,nface
                    i4p = i4 + nface
                    flux_z(i4,i,j,k,ieq) = 0.5*( (Fface_z(i4,ieq) + Fface_z(i4p,ieq)) &
                                         - cfrz(i4,ieq)*(Qface_z(i4p,ieq) - Qface_z(i4,ieq)) )
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
            !------------------------------------------
        end do
        end do
        !$OMP END PARALLEL DO
        end do
    end subroutine calc_flux_z
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    subroutine flux_calc(Q_io)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io

        call calc_flux_x(Q_io, flux_x)
        call calc_flux_y(Q_io, flux_y)
        call calc_flux_z(Q_io, flux_z)
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
            if (ieos == 2) plr(k) = P_1*(rhov(k)**7.2 - 1.) + P_base + plr(k)
            rtrho(k) = sqrt(rhov(k))
        end do

        do k=1,nface
            k2 = k + nface
            if (ieos == 2) then
                cslr(k) = vlr(k,1) - sqrt(7.2*P_1*rhov(k)**6.2 + plr(k)*rho_i)     ! lambda_M(Q_l)
                cslr(k2) = vlr(k2,1) + sqrt(7.2*P_1*rhov(k2)**6.2 + plr(k2)*rho_i) ! lambda_P(Q_r)
            else
                cslr(k) = vlr(k,1) - sqrt(aindex*plr(k)/rhov(k))       ! lambda_M(Q_l)
                cslr(k2) = vlr(k2,1) + sqrt(aindex*plr(k2)/rhov(k2) )       ! lambda_P(Q_r)
            end if
        end do

        if (ibatten == 1) then  ! compute wave speeds using Roe averages following Batten, 1997
          do k=1,nface
            k2 = k + nface
            rtrho_i(k) = 1.0/(rtrho(k) + rtrho(k2))
            qtilde(k,1) = ( rtrho(k)*vlr(k,1) + rtrho(k2)*vlr(k2,1) ) * rtrho_i(k)
            qtilde(k,2) = ( rtrho(k)*vlr(k,2) + rtrho(k2)*vlr(k2,2) ) * rtrho_i(k)
            qtilde(k,3) = ( rtrho(k)*vlr(k,3) + rtrho(k2)*vlr(k2,3) ) * rtrho_i(k)
            qsq(k) = qtilde(k,1)**2 + qtilde(k,2)**2 + qtilde(k,3)**2
            hlr(k)  = (Qlr(k,enj)  + plr(k)) /rhov(k)
            hlr(k2) = (Qlr(k2,enj) + plr(k2))/rhov(k2)
            qtilde(k,4) = (rtrho(k)*hlr(k) + rtrho(k2)*hlr(k2))*rtrho_i(k)
            ctsq(k) = aindm1*(qtilde(k,4) - 0.5*qsq(k))
          end do
          if (minval(ctsq) >= 0.0) then
            ctilde = sqrt(ctsq)
            qslr(1:nface) = qtilde(1:nface,1) - ctilde(1:nface)       ! lambda_M(Q_Roe)
            qslr(nface1:nfe) = qtilde(nface1:nfe,1) + ctilde(nface1:nfe)  ! lambda_P(Q_Roe)
          end if
          if (minval(ctsq) < 0.0) then
            ibatten = 0
          end if
        end if

        if (ibatten == 0) then
          do k=1,nface
            k2 = k + nface
            if (ieos == 2) then
                qslr(k) = vlr(k2,1) - sqrt(7.2*P_1*rhov(k2)**6.2 + plr(k2)*rho_i) ! lambda_M(Q_r)
                qslr(k2) = vlr(k,1) + sqrt(7.2*P_1*rhov(k)**6.2 + plr(k)*rho_i)   ! lambda_P(Q_l)
            else
                qslr(k) = vlr(k2,1) - sqrt(aindex*plr(k2)/rhov(k2))   ! lambda_M(Q_r)
                qslr(k2) = vlr(k,1) + sqrt(aindex*plr(k)/rhov(k))     ! lambda_P(Q_l)
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
            Qstar(k,1) = rhov(k)*sq_lr(k)*slrm_i(k)                             ! Eq. (35) of Batten ,1997
            Qstar(k,2) = (sq_lr(k)*Qlr(k,iparr) + pstar(i4) - plr(k))*slrm_i(k) ! Eq. (37-39) of Batten, 1997
            Qstar(k,3) = sq_lr(k)*Qlr(k,iperp1)*slrm_i(k)                       ! Eq. (37-39) of Batten, 1997
            Qstar(k,4) = sq_lr(k)*Qlr(k,iperp2)*slrm_i(k)                       ! Eq. (37-39) of Batten, 1997
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
    subroutine innerintegral(Q_in)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(in) :: Q_in
        real, dimension(npg,nQ)     :: Qinner, finner_x,finner_y,finner_z
        real, dimension(nbastot,nQ) :: int_r
        real :: sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9
        real :: wgt3d_fin_x_ieq(npg), wgt3d_fin_x_ieq_3(npg)
        real :: wgt3d_fin_y_ieq(npg), wgt3d_fin_y_ieq_3(npg)
        real :: wgt3d_fin_z_ieq(npg), wgt3d_fin_z_ieq_3(npg)
        integer :: i,j,k,ieq,ipg,ir

        integral_r(:,:,:,:,:) = 0.

        do k = 1,nz
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Qinner,integral_r,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum8) FIRSTPRIVATE(finner_x,finner_y,finner_z,int_r)
        do j = 1,ny
        do i = 1,nx

            do ieq = 1,nQ
                do ipg = 1,npg
                    Qinner(ipg,ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_in(i,j,k,ieq,1:nbasis))
                end do
            end do

            ! NOTE: this can be made more efficient since temp vars can be re-used
            call Fpts_x(Qinner, finner_x, i,j,k, npg)
            call Fpts_y(Qinner, finner_y, i,j,k, npg)
            call Fpts_z(Qinner, finner_z, i,j,k, npg)
            ! call flux_calc_pnts_r_v0(Qinner,finner_x,1,npg)
            ! call flux_calc_pnts_r_v0(Qinner,finner_y,2,npg)
            ! call flux_calc_pnts_r_v0(Qinner,finner_z,3,npg)

          do ieq = 1,nQ
            wgt3d_fin_x_ieq(1:npg) = wgt3d(1:npg) * finner_x(1:npg,ieq)
            wgt3d_fin_y_ieq(1:npg) = wgt3d(1:npg) * finner_y(1:npg,ieq)
            wgt3d_fin_z_ieq(1:npg) = wgt3d(1:npg) * finner_z(1:npg,ieq)
            wgt3d_fin_x_ieq_3(1:npg) = 3 * wgt3d_fin_x_ieq(1:npg)
            wgt3d_fin_y_ieq_3(1:npg) = 3 * wgt3d_fin_y_ieq(1:npg)
            wgt3d_fin_z_ieq_3(1:npg) = 3 * wgt3d_fin_z_ieq(1:npg)

            sum1 = 0
            sum2 = 0
            sum3 = 0
            do ipg=1,npg
                sum1 = sum1 + wgt3d_fin_x_ieq(ipg)
                sum2 = sum2 + wgt3d_fin_y_ieq(ipg)
                sum3 = sum3 + wgt3d_fin_z_ieq(ipg)
            end do
            int_r(kx,ieq) = dxid4*cbasis(kx)*sum1
            int_r(ky,ieq) = dyid4*cbasis(ky)*sum2
            int_r(kz,ieq) = dzid4*cbasis(kz)*sum3

            if ( nbasis > 4 ) then
                !=== ibitri == 1 .or. iquad > 2 ===========
                if ( ibitri == 1 .or. iquad > 2 ) then
                    sum1 = 0
                    sum2 = 0
                    sum3 = 0
                    sum4 = 0
                    sum5 = 0
                    sum6 = 0
                    do ipg=1,npg
                        sum1 = sum1 + wgt3d_fin_y_ieq(ipg) * bfvals_int(ipg,kz)
                        sum2 = sum2 + wgt3d_fin_z_ieq(ipg) * bfvals_int(ipg,ky)
                        sum3 = sum3 + wgt3d_fin_x_ieq(ipg) * bfvals_int(ipg,kz)
                        sum4 = sum4 + wgt3d_fin_z_ieq(ipg) * bfvals_int(ipg,kx)
                        sum5 = sum5 + wgt3d_fin_x_ieq(ipg) * bfvals_int(ipg,ky)
                        sum6 = sum6 + wgt3d_fin_y_ieq(ipg) * bfvals_int(ipg,kx)
                    end do
                    int_r(kyz,ieq) = cbasis(kyz)*(dyid4*sum1 + dzid4*sum2)
                    int_r(kzx,ieq) = cbasis(kzx)*(dxid4*sum3 + dzid4*sum4)
                    int_r(kxy,ieq) = cbasis(kxy)*(dxid4*sum5 + dyid4*sum6)
                end if
                !=== iquad > 2 ============================
                if ( iquad > 2 ) then
                    sum1 = 0
                    sum2 = 0
                    sum3 = 0
                    do ipg=1,npg
                        sum1 = sum1 + wgt3d_fin_x_ieq_3(ipg) * bfvals_int(ipg,kx)
                        sum2 = sum2 + wgt3d_fin_y_ieq_3(ipg) * bfvals_int(ipg,ky)
                        sum3 = sum3 + wgt3d_fin_z_ieq_3(ipg) * bfvals_int(ipg,kz)
                    end do
                    int_r(kxx,ieq) = dxid4*cbasis(kxx)*sum1
                    int_r(kyy,ieq) = dyid4*cbasis(kyy)*sum2
                    int_r(kzz,ieq) = dzid4*cbasis(kzz)*sum3
                end if
                !=== ibitri == 1 .or. iquad > 3 ===========
                if ( ibitri == 1 .or. iquad > 3 ) then
                    sum7 = 0
                    sum8 = 0
                    sum9 = 0
                    do ipg=1,npg
                        sum7 = sum7 + wgt3d_fin_x_ieq(ipg) * bfvals_int(ipg,kyz)
                        sum8 = sum8 + wgt3d_fin_y_ieq(ipg) * bfvals_int(ipg,kzx)
                        sum9 = sum9 + wgt3d_fin_z_ieq(ipg) * bfvals_int(ipg,kxy)
                    end do
                    int_r(kxyz,ieq) = cbasis(kxyz)*(dxid4*sum7 + dyid4*sum8 + dzid4*sum9)
                end if
                !=== (ibitri == 1 .and. iquad > 2) .or. iquad > 3 ========
                if ( (ibitri == 1 .and. iquad > 2) .or. iquad > 3 ) then
                  int_r(kyzz,ieq) =                                                               &
                      cbasis(kyzz)*dyid4 * sum( wgt3d_fin_y_ieq(1:npg) * bfvals_int(1:npg,kzz) )  &
                    + cbasis(kyzz)*dzid4 * sum( wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kyz) )
                  int_r(kzxx,ieq) =                                                               &
                      cbasis(kzxx)*dzid4 * sum( wgt3d_fin_z_ieq(1:npg) * bfvals_int(1:npg,kxx) )  &
                    + cbasis(kzxx)*dxid4 * sum( wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kzx) )
                  int_r(kxyy,ieq) =                                                               &
                      cbasis(kxyy)*dxid4 * sum( wgt3d_fin_x_ieq(1:npg) * bfvals_int(1:npg,kyy) )  &
                    + cbasis(kxyy)*dyid4 * sum( wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kxy) )
                  int_r(kyyz,ieq) =                                                               &
                      cbasis(kyyz)*dzid4 * sum( wgt3d_fin_z_ieq(1:npg) * bfvals_int(1:npg,kyy) )  &
                    + cbasis(kyyz)*dyid4 * sum( wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kyz) )
                  int_r(kzzx,ieq) =                                                               &
                      cbasis(kzzx)*dxid4 * sum( wgt3d_fin_x_ieq(1:npg) * bfvals_int(1:npg,kzz) )  &
                    + cbasis(kzzx)*dzid4 * sum( wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kzx) )
                  int_r(kxxy,ieq) =                                                               &
                      cbasis(kxxy)*dyid4 * sum( wgt3d_fin_y_ieq(1:npg) * bfvals_int(1:npg,kxx) )  &
                    + cbasis(kxxy)*dxid4 * sum( wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kxy) )
                end if
                !=== iquad > 3 ============================
                if ( iquad > 3 ) then
                    int_r(kxxx,ieq) = cbasis(kxxx)*dxid4 *                                        &
                        sum( wgt3d_fin_x_ieq(1:npg)*(7.5*bfvals_int(1:npg,kx)**2 - 1.5) )
                    int_r(kyyy,ieq) = cbasis(kyyy)*dyid4 *                                        &
                        sum( wgt3d_fin_y_ieq(1:npg)*(7.5*bfvals_int(1:npg,ky)**2 - 1.5) )
                    int_r(kzzz,ieq) = cbasis(kzzz)*dzid4 *                                        &
                        sum( wgt3d_fin_z_ieq(1:npg)*(7.5*bfvals_int(1:npg,kz)**2 - 1.5) )
                end if
                !=== ibitri == 1 .and. iquad > 2 ==========
                if ( ibitri == 1 .and. iquad > 2 ) then
                  int_r(kyyzz,ieq) =                                                              &
                      cbasis(kyyzz)*dyid4 * sum(wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kyzz))  &
                    + cbasis(kyyzz)*dzid4 * sum(wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kyyz))
                   int_r(kzzxx,ieq) =                                                             &
                      cbasis(kzzxx)*dzid4 * sum(wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kzxx))  &
                    + cbasis(kzzxx)*dxid4 * sum(wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kzzx))
                   int_r(kxxyy,ieq) =                                                             &
                      cbasis(kxxyy)*dxid4 * sum(wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kxyy))  &
                    + cbasis(kxxyy)*dyid4 * sum(wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kxxy))
                   int_r(kyzxx,ieq) =                                                             &
                      cbasis(kyzxx)*dxid4 * sum(wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kxyz))  &
                    + cbasis(kyzxx)*dyid4 * sum(wgt3d_fin_y_ieq(1:npg) * bfvals_int(1:npg,kzxx))  &
                    + cbasis(kyzxx)*dzid4 * sum(wgt3d_fin_z_ieq(1:npg) * bfvals_int(1:npg,kxxy))
                  int_r(kzxyy,ieq) =                                                              &
                      cbasis(kzxyy)*dyid4 * sum(wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kxyz))  &
                    + cbasis(kzxyy)*dzid4 * sum(wgt3d_fin_z_ieq(1:npg) * bfvals_int(1:npg,kxyy))  &
                    + cbasis(kzxyy)*dxid4 * sum(wgt3d_fin_x_ieq(1:npg) * bfvals_int(1:npg,kyyz))
                  int_r(kxyzz,ieq) =                                                              &
                      cbasis(kxyzz)*dzid4 * sum(wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kxyz))  &
                    + cbasis(kxyzz)*dxid4 * sum(wgt3d_fin_x_ieq(1:npg) * bfvals_int(1:npg,kyzz))  &
                    + cbasis(kxyzz)*dyid4 * sum(wgt3d_fin_y_ieq(1:npg) * bfvals_int(1:npg,kzzx))
                  int_r(kxyyzz,ieq) =                                                             &
                      cbasis(kxyyzz)*dxid4 *sum(wgt3d_fin_x_ieq(1:npg) * bfvals_int(1:npg,kyyzz)) &
                    + cbasis(kxyyzz)*dyid4 *sum(wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kxyzz)) &
                    + cbasis(kxyyzz)*dzid4 *sum(wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kzxyy))
                  int_r(kyzzxx,ieq) =                                                             &
                      cbasis(kyzzxx)*dyid4 *sum(wgt3d_fin_y_ieq(1:npg) * bfvals_int(1:npg,kzzxx)) &
                    + cbasis(kyzzxx)*dzid4 *sum(wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kyzxx)) &
                    + cbasis(kyzzxx)*dxid4 *sum(wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kxyzz))
                  int_r(kzxxyy,ieq) =                                                             &
                      cbasis(kzxxyy)*dzid4 *sum(wgt3d_fin_z_ieq(1:npg) * bfvals_int(1:npg,kxxyy)) &
                    + cbasis(kzxxyy)*dxid4 *sum(wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kzxyy)) &
                    + cbasis(kzxxyy)*dyid4 *sum(wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kyzxx))
                  int_r(kxxyyzz,ieq) =                                                             &
                      cbasis(kxxyyzz)*dxid4*sum(wgt3d_fin_x_ieq_3(1:npg)*bfvals_int(1:npg,kxyyzz)) &
                    + cbasis(kxxyyzz)*dyid4*sum(wgt3d_fin_y_ieq_3(1:npg)*bfvals_int(1:npg,kyzzxx)) &
                    + cbasis(kxxyyzz)*dzid4*sum(wgt3d_fin_z_ieq_3(1:npg)*bfvals_int(1:npg,kzxxyy))
                end if
            end if
          end do

          do ieq = 1,nQ
            do ir = 1,nbasis
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
    subroutine glflux(Q_io, dt)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        real, intent(inout) :: dt
        real :: sumx,sumy,sumz
        integer i,j,k,ieq,ir,ifa, i1,j1,k1

        !#########################################################
        ! Step 1: Calculate fluxes for boundaries of each cell
        !   --> flux_x, flux_y, flux_z  (used only in flux.f90)
        !---------------------------------------------------------
        ! NOTE: generation of Gaussian random matrices is moved to prep_advance
        !       b/c we need consistent random samples at interfaces of MPI domains
        ! call flux_calc(Q_io)
        call calc_flux_x(Q_io, flux_x)
        call calc_flux_y(Q_io, flux_y)
        call calc_flux_z(Q_io, flux_z)

        !#########################################################
        ! Step 2: Calc inner integral for each cell
        !   --> integral_r  (used only in flux.f90)
        !---------------------------------------------------------
        call innerintegral(Q_io)  ! call innerintegral_v0(Q_io)

        !#########################################################
        ! Step 3: Calc (total) "Gauss-Legendre flux" for each cell
        !   --> glflux_r  (used by advance_time_level_gl)
        !   --
        !---------------------------------------------------------
        do ieq = 1,nQ
        do k = 1,nz       !FIRSTPRIVATE(flux_x,flux_y,flux_z)  !$OMP PARALLEL DO DEFAULT(SHARED)
            k1 = k+1
        do j = 1,ny
            j1 = j+1
        do i = 1,nx
            i1 = i+1
            sumx = 0
            sumy = 0
            sumz = 0
    	    do ifa = 1,nface
        		sumx = sumx + dxid4*wgt2d(ifa)*( flux_x(ifa,i1,j,k,ieq)-flux_x(ifa,i,j,k,ieq) )
        		sumy = sumy + dyid4*wgt2d(ifa)*( flux_y(ifa,i,j1,k,ieq)-flux_y(ifa,i,j,k,ieq) )
        		sumz = sumz + dzid4*wgt2d(ifa)*( flux_z(ifa,i,j,k1,ieq)-flux_z(ifa,i,j,k,ieq) )
    	    end do
    	    glflux_r(i,j,k,ieq,1) = sumx + sumy + sumz
        end do
        do ir=2,nbasis
            do i = 1,nx
                i1 = i+1
                sumx = 0
                sumy = 0
                sumz = 0
                do ifa = 1,nface
            		sumx = sumx + wgtbf_xmp(ifa,1,ir)*flux_x(ifa,i,j,k,ieq)                     &
                                + wgtbf_xmp(ifa,2,ir)*flux_x(ifa,i1,j,k,ieq)
            		sumy = sumy + wgtbf_ymp(ifa,1,ir)*flux_y(ifa,i,j,k,ieq)                     &
                                + wgtbf_ymp(ifa,2,ir)*flux_y(ifa,i,j1,k,ieq)
            		sumz = sumz + wgtbf_zmp(ifa,1,ir)*flux_z(ifa,i,j,k,ieq)                     &
                                + wgtbf_zmp(ifa,2,ir)*flux_z(ifa,i,j,k1,ieq)
        	    end do
        	    glflux_r(i,j,k,ieq,ir) = sumx + sumy + sumz - integral_r(i,j,k,ieq,ir)
            end do
        end do
        end do
        end do
        end do   !$OMP END PARALLEL DO
    end subroutine glflux
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    ! subroutine glflux(Q_io, dt)
    !     implicit none
    !     real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
    !     real, intent(inout) :: dt
    !     real :: sumx,sumy,sumz
    !     integer i,j,k,ieq,ir,ifa, i1,j1,k1
    !
    !     if (llns) call get_region_Sflux_xyz(Sflux_x, Sflux_y, Sflux_z, dt)
    !
    !     !#########################################################
    !     ! Step 1: Calculate fluxes for boundaries of each cell
    !     !   --> flux_x, flux_y, flux_z  (used only in flux.f90)
    !     !---------------------------------------------------------
    !     ! call flux_calc(Q_io)
    !     call calc_flux_x(Q_io, flux_x)
    !     call calc_flux_y(Q_io, flux_y)
    !     call calc_flux_z(Q_io, flux_z)
    !
    !     !#########################################################
    !     ! Step 2: Calc inner integral for each cell
    !     !   --> integral_r  (used only in flux.f90)
    !     !---------------------------------------------------------
    !     call innerintegral(Q_io)  ! call innerintegral_v0(Q_io)
    !
    !     !#########################################################
    !     ! Step 3: Calc (total) "Gauss-Legendre flux" for each cell
    !     !   --> glflux_r  (used by advance_time_level_gl)
    !     !   --
    !     !---------------------------------------------------------
    !     do ieq = 1,nQ
    !     do k = 1,nz       !FIRSTPRIVATE(flux_x,flux_y,flux_z)  !$OMP PARALLEL DO DEFAULT(SHARED)
    !         k1 = k+1
    !     do j = 1,ny
    !         j1 = j+1
    !     do i = 1,nx
    !         i1 = i+1
    !         sumx = 0
    !         sumy = 0
    !         sumz = 0
    ! 	    do ifa = 1,nface
    !     		sumx = sumx + dxid4*wgt2d(ifa)*( flux_x(ifa,i1,j,k,ieq)-flux_x(ifa,i,j,k,ieq) )
    !     		sumy = sumy + dyid4*wgt2d(ifa)*( flux_y(ifa,i,j1,k,ieq)-flux_y(ifa,i,j,k,ieq) )
    !     		sumz = sumz + dzid4*wgt2d(ifa)*( flux_z(ifa,i,j,k1,ieq)-flux_z(ifa,i,j,k,ieq) )
    ! 	    end do
    ! 	    glflux_r(i,j,k,ieq,1) = sumx + sumy + sumz
    !     end do
    !     do ir=2,nbasis
    !         do i = 1,nx
    !             i1 = i+1
    !             sumx = 0
    !             sumy = 0
    !             sumz = 0
    !             do ifa = 1,nface
    !         		sumx = sumx + wgtbf_xmp(ifa,1,ir)*flux_x(ifa,i,j,k,ieq)                     &
    !                             + wgtbf_xmp(ifa,2,ir)*flux_x(ifa,i1,j,k,ieq)
    !         		sumy = sumy + wgtbf_ymp(ifa,1,ir)*flux_y(ifa,i,j,k,ieq)                     &
    !                             + wgtbf_ymp(ifa,2,ir)*flux_y(ifa,i,j1,k,ieq)
    !         		sumz = sumz + wgtbf_zmp(ifa,1,ir)*flux_z(ifa,i,j,k,ieq)                     &
    !                             + wgtbf_zmp(ifa,2,ir)*flux_z(ifa,i,j,k1,ieq)
    !     	    end do
    !     	    glflux_r(i,j,k,ieq,ir) = sumx + sumy + sumz - integral_r(i,j,k,ieq,ir)
    !         end do
    !     end do
    !     end do
    !     end do
    !     end do   !$OMP END PARALLEL DO
    ! end subroutine glflux
!-------------------------------------------------------------------------------


!-----------------------------------------------------------!
!*******************calculate freezing speeds***************!
!-----------------------------------------------------------!
    real function cfcal(Qcf, direction)
        implicit none
        real, dimension(nQ), intent(in) :: Qcf
        integer, intent(in) :: direction
        real :: dn,dni, vx,vy,vz, P, cs

        dn = Qcf(rh)
        dni = 1./dn
        vx = Qcf(mx)*dni
        vy = Qcf(my)*dni
        vz = Qcf(mz)*dni

        !=== pressure =================
        select case(ivis)
            case(0)
                P = aindm1*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
            case(1)
                P = aindm1*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
            case(2)
                P = aindm1*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))  ! NOTE: from GFV1
                ! P = c1d3 * ( Qcf(exx) + Qcf(eyy) + Qcf(ezz) - dn*(vx**2 + vy**2 + vz**2) )
        end select
        !=== sound speed ==============
        select case(ieos)
            case(1)
                cs = sqrt(aindex*P*dni)
            case(2)
                cs = sqrt(7.2*P_1*dn**6.2)
        end select
        !=== freezing speed ===========
        select case(direction)
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
    subroutine limiter(Q_io)

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_io
        integer i,j,k, ieq, ipge, minindex, ir
        real Qedge(npge,nQ),theta,Qmin(nQ), deltaQ(nQ)
        real epsi, Qrhmin, QPmin, P(npge), Pave, dn, dni, epsiP, thetaj
        real*8 a, b, c

        epsi  = rh_floor
        epsiP = rh_floor*T_floor

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            if ( Q_io(i,j,k,rh,1) < rh_floor ) then
                do ir=2,nbasis
                    Q_io(i,j,k,rh:en,ir) = 0.0
                end do
                Q_io(i,j,k,rh,1) = rh_floor
            else
                do ipge = 1,npge
                    Qedge(ipge,rh) = sum(bf_faces(ipge,1:nbasis)*Q_io(i,j,k,rh,1:nbasis))
                end do
                !--------------------------------
                Qrhmin = minval(Qedge(:,rh))
                if ( Qrhmin < epsi ) then
                    theta = (epsi-Q_io(i,j,k,rh,1))/(Qrhmin-Q_io(i,j,k,rh,1))
                    if ( theta > 1 ) then
                        theta = 1.
                    else if ( theta < 0 ) then
                        theta = 0.
                    end if
                    do ir=2,nbasis
                        Q_io(i,j,k,rh,ir) = abs(theta)*Q_io(i,j,k,rh,ir)
                    end do
                end if
                !----------------------------------------------------
                Pave = aindm1*(Q_io(i,j,k,en,1)                                  &
                        - 0.5*( Q_io(i,j,k,mx,1)**2                              &
                              + Q_io(i,j,k,my,1)**2                              &
                              + Q_io(i,j,k,mz,1)**2 ) / Q_io(i,j,k,rh,1) )
                !----------------------------------------------------
                if ( Pave < epsiP ) then
                    do ir=2,nbasis
                        Q_io(i,j,k,rh:en,ir) = 0.0
                    end do
                else
                    theta = 1.
                    do ipge = 1,npge
                      do ieq = rh,en
                        Qedge(ipge,ieq) = sum(bf_faces(ipge,1:nbasis)*Q_io(i,j,k,ieq,1:nbasis))
                      end do
                      !--------------------------------
                      dn  = Qedge(ipge,rh)
                      dni = 1./dn
                      P(ipge) = aindm1*(Qedge(ipge,en)                           &
                          - 0.5*dni*(Qedge(ipge,mx)**2 + Qedge(ipge,my)**2 + Qedge(ipge,mz)**2))
                      !--------------------------------
                      if ( P(ipge) < epsiP ) then
                        if ( Pave .ne. P(ipge) ) then
                            thetaj = (Pave - epsiP)/(Pave - P(ipge))
                            theta = min(theta,thetaj)
                        end if
                      end if
                    end do
                    if ( theta > 1 ) then
                        theta = 1.
                    else if ( theta < 0 ) then
                        theta = 0.
                    end if
                    do ir=2,nbasis
                        Q_io(i,j,k,rh:en,ir) = theta*Q_io(i,j,k,rh:en,ir)
                    end do
                end if
            end if
        end do
        end do
        end do

    end subroutine limiter
!-------------------------------------------------------------------------------


end module flux
