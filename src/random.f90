!***** RANDOM.F90 **************************************************************
! NOTE: this implementation is slow probably becaues of the dynamic
!       allocation of the arrays used for random matrices/vectors
!*******************************************************************************
module random

use MKL_VSL_TYPE
use MKL_VSL

use params
use spatial
use timestep

!===========================================================================
! MKL VSL parameters
!------------------------------------------------------------
real :: vsl_errcode
TYPE (VSL_STREAM_STATE) :: vsl_stream

integer, parameter :: vsl_brng   = VSL_BRNG_MCG31
integer, parameter :: vsl_method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
real, parameter :: vsl_mean  = 0.0
real, parameter :: vsl_sigma = 1.0
!===========================================================================

real :: dvdti, sqrt_dvdti


! NOTE: It might make sense to use global arrays for stochastic
!   stuff at some point, especially if in separate module
! real GRM_x(nface, 1:nx+1, ny,     nz,     3,3)
! real GRM_y(nface, nx,     1:ny+1, nz,     3,3)
! real GRM_z(nface, nx,     ny,     1:nz+1, 3,3)
! real Sflux_x(nface, 1:nx+1, ny,     nz,     3,3)
! real Sflux_y(nface, nx,     1:ny+1, nz,     3,3)
! real Sflux_z(nface, nx,     ny,     1:nz+1, 3,3)

contains

!-------------------------------------------------------------------------------
    subroutine random_init(seed)
        integer, intent(in) :: seed
        vsl_errcode = vslnewstream(vsl_stream, vsl_brng, seed + iam)

        dvi   = (dxi*dyi*dzi) * npg  ! NOTE: taking cell volume to be 1/(# internal quad pts)
        dvdti = dvi/dt
        sqrt_dvdti = sqrt(dvdti)
    end subroutine random_init
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine random_cleanup()
        vsl_errcode = vsldeletestream( vsl_stream )
    end subroutine random_cleanup
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine get_GRM(GRMpnts_r, npnts)

        implicit none
        integer ife,npnts,b,e,vsl_ndim
        real, dimension(npnts,3,3) :: GRMpnts_r

        real, dimension(9*npnts) :: grn
        ! real, allocatable :: grn(:)

        vsl_ndim = 9*npnts
        ! allocate(grn(vsl_ndim))

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
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine get_GRM_symmetric(GRMpnts_r, npnts)

        implicit none
        integer ife,npnts,b,e,vsl_ndim
        real, dimension(npnts,3,3) :: GRMpnts_r

        real, dimension(6*npnts) :: grn
        ! real, allocatable :: grn(:)

        vsl_ndim = 6*npnts
        ! allocate(grn(vsl_ndim))

        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)

        ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
        do ife = 1,npnts
            e = 6*ife
            b = e - 5
            GRMpnts_r(ife,1,1:3) = grn(b:b+2)
            GRMpnts_r(ife,2,2:3) = grn(b+3:b+4)
            GRMpnts_r(ife,3,3)   = grn(e)

            GRMpnts_r(ife,2,1)   = GRMpnts_r(ife,1,2)
            GRMpnts_r(ife,3,1)   = GRMpnts_r(ife,1,3)
            GRMpnts_r(ife,3,2)   = GRMpnts_r(ife,2,3)
        end do

    end subroutine get_GRM_symmetric
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine random_stresses_pnts_r(Spnts_r, npnts)

        implicit none
        integer ife,npnts
        real, dimension(npnts,3,3) :: GRMpnts_r, Spnts_r
        real Gxx,Gyy,Gzz,Gxy,Gxz,Gyz
        real eta_d,zeta_d
        real trG,trGd3,trG_zeta

        eta_d = eta_sd * sqrt_dvdti
        zeta_d = zeta_sd * sqrt_dvdti
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
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    ! subroutine random_stresses_pnts_r(Spnts_r, npnts)
    !
    !     implicit none
    !     integer ife,npnts
    !     real, dimension(npnts,3,3) :: GRMpnts_r, Spnts_r
    !     real Gxx,Gyy,Gzz,Gxy,Gxz,Gyz
    !     real eta_d,zeta_d
    !     real trG,trGd3,trG_zeta
    !
    !     eta_d = eta_sd * sqrt_dvdti
    !     zeta_d = zeta_sd * sqrt_dvdti
    !     call get_GRM_symmetric(GRMpnts_r, npnts)
    !
    !     ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
    !     do ife = 1,npnts
    !         Gxx = GRMpnts_r(ife,1,1)
    !         Gyy = GRMpnts_r(ife,2,2)
    !         Gzz = GRMpnts_r(ife,3,3)
    !
    !         Gxy = GRMpnts_r(ife,1,2)
    !         Gxz = GRMpnts_r(ife,1,3)
    !         Gyz = GRMpnts_r(ife,2,3)
    !
    !         trG = (Gxx + Gyy + Gzz)
    !         trGd3 = trG/3.0
    !         trG_zeta = zeta_sd*trG
    !
    !         Spnts_r(ife,1,1) = eta_sd*(Gxx - trGd3) + trG_zeta
    !         Spnts_r(ife,2,2) = eta_sd*(Gyy - trGd3) + trG_zeta
    !         Spnts_r(ife,3,3) = eta_sd*(Gzz - trGd3) + trG_zeta
    !
    !         Spnts_r(ife,1,2) = eta_sd*Gxy
    !         Spnts_r(ife,1,3) = eta_sd*Gxz
    !         Spnts_r(ife,2,3) = eta_sd*Gyz
    !         Spnts_r(ife,2,1) = Spnts_r(ife,1,2)
    !         Spnts_r(ife,3,1) = Spnts_r(ife,1,3)
    !         Spnts_r(ife,3,2) = Spnts_r(ife,2,3)
    !     end do
    !
    ! end subroutine random_stresses_pnts_r
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine random_heatflux_pnts_r(Hpnts_r, npnts)

        implicit none
        integer ife,npnts,b,e,vsl_ndim
        real, dimension(npnts,3) :: GRVpnts_r, Hpnts_r
        real Gx,Gy,Gz,kappa_d
        real, allocatable :: grn(:)

        vsl_ndim = npnts*3
        allocate(grn(vsl_ndim))

        kappa_d = kappa_sd * sqrt_dvdti

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
!-------------------------------------------------------------------------------


end module random
