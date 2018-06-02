!***** RANDOM.F90 **************************************************************
! NOTE: this implementation is slow probably becaues of the dynamic
!       allocation of the arrays used for random matrices/vectors
!*******************************************************************************
module random

use MKL_VSL_TYPE
use MKL_VSL

use helpers
use params
use spatial
use timestep

!===========================================================================
! MKL VSL parameters
!------------------------------------------------------------
real :: vsl_errcode
TYPE(VSL_STREAM_STATE) :: vsl_stream
! TYPE(VSL_STREAM_STATE), allocatable, dimension(:) :: vsl_stream
! TYPE NEW_TYPE
!     TYPE (VSL_STREAM_STATE) :: s1
! END TYPE NEW_TYPE
! TYPE(NEW_TYPE), allocatable, dimension(:) :: vsl_stream

integer, parameter :: vsl_brng   = VSL_BRNG_MT2203 !VSL_BRNG_MT19937 ! VSL_BRNG_MCG31
! integer, parameter :: vsl_method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
! integer, parameter :: vsl_method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
integer, parameter :: vsl_method = VSL_RNG_METHOD_GAUSSIAN_ICDF
real, parameter :: vsl_mean  = 0.0
real, parameter :: vsl_sigma = 1.0
!===========================================================================

integer, parameter :: npts_llns = 1  !npg ! number of stochastic sampling points per cell face
real :: dvdti, sqrt_dvdti


real ::   GRM(nface,1:nx1,1:ny1,1:nz1,3,3)
real :: Sflux(nface,1:nx1,1:ny1,1:nz1,3,3)

! NOTE: It might make sense to use global arrays for stochastic
!   stuff at some point, especially if in separate module
real ::   GRM_x(npts_llns, 1:nx1, ny,    nz,    3,3)
real ::   GRM_y(npts_llns, nx,    1:ny1, nz,    3,3)
real ::   GRM_z(npts_llns, nx,    ny,    1:nz1, 3,3)
real :: Sflux_x(npts_llns, 1:nx1, ny,    nz,    3,3)
real :: Sflux_y(npts_llns, nx,    1:ny1, nz,    3,3)
real :: Sflux_z(npts_llns, nx,    ny,    1:nz1, 3,3)

contains

!-------------------------------------------------------------------------------
    subroutine random_init(seed)
        ! integer, intent(in) :: npts
        integer, optional, intent(in) :: seed
        integer :: default_seed

        Sflux(:,:,:,:,:,:) = 0.0
        ! allocate(vsl_stream(numprocs))

        call mpi_print(iam, 'Selected fluctuating hydrodynamics (LLNS) model')

        if (present(seed)) then
            vsl_errcode = vslnewstream(vsl_stream, vsl_brng+iam, seed)
        else
            default_seed = 1234567 + 7654321*(1+iam)
            call mpi_print(iam, 'No seed specified. Using seed=1234567+7654321*(1+RANK)')
            vsl_errcode = vslnewstream(vsl_stream, vsl_brng+iam, default_seed)
        end if

        dvi   = (dxi*dyi*dzi) * npts_llns ! NOTE: taking cell volume as 1/(# internal quad pts)
        dvdti = dvi/dt
        sqrt_dvdti = sqrt(dvdti)
    end subroutine random_init
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine random_cleanup()
        vsl_errcode = vsldeletestream(vsl_stream)
    end subroutine random_cleanup
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine get_region_GRM(GRM)
        ! NOTE: this currently only populates a since face point!
        implicit none
        real, dimension(nface,1:nx1,1:ny1,1:nz1,3,3), intent(out) :: GRM
        real, allocatable :: grn(:)
        integer :: ipnt, vsl_ndim

        vsl_ndim = npts_llns*nx1*ny1*nz1*9
        allocate(grn(vsl_ndim))
        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)
        do ipnt = 1,nface
            ! GRM(ipnt,:,:,:,:,:) = reshape( grn, (/ npts_llns,nx1,ny1,nz1,3,3 /) )
            GRM(ipnt,:,:,:,:,:) = reshape( grn, (/ nx1,ny1,nz1,3,3 /) )
        end do
        deallocate(grn)
    end subroutine get_region_GRM
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine get_region_Sflux(Sflux, dt)
        implicit none
        real, dimension(nface,1:nx1,1:ny1,1:nz1,3,3), intent(out) :: Sflux
        real, intent(inout) :: dt
        real, dimension(nface,1:nx1,1:ny1,1:nz1,3,3) :: GRM
        integer :: i,j,k,ipnt
        real :: Gxx,Gyy,Gzz,Gxy,Gxz,Gyz
        real :: eta_d,zeta_d
        real :: trG,trGd3,trG_zeta

        dvdti = npts_llns*(dxi*dyi*dzi)/dt
        sqrt_dvdti = sqrt(dvdti)
        eta_d  = eta_sd  * sqrt_dvdti
        zeta_d = zeta_sd * sqrt_dvdti
        call get_region_GRM(GRM)

        ! NOTE: probz not necessary to generate up to nx+1, etc.
        do k=1,nz1
        do j=1,ny1
        do i=1,nx1
            do ipnt=1,npts_llns
                Gxx = GRM(ipnt,i,j,k,1,1)
                Gyy = GRM(ipnt,i,j,k,2,2)
                Gzz = GRM(ipnt,i,j,k,3,3)

                Gxy = sqrt2i*( GRM(ipnt,i,j,k,1,2) + GRM(ipnt,i,j,k,2,1) )
                Gxz = sqrt2i*( GRM(ipnt,i,j,k,1,3) + GRM(ipnt,i,j,k,3,1) )
                Gyz = sqrt2i*( GRM(ipnt,i,j,k,2,3) + GRM(ipnt,i,j,k,3,2) )

                Sflux(ipnt,i,j,k,1,2) = eta_d*Gxy
                Sflux(ipnt,i,j,k,1,3) = eta_d*Gxz
                Sflux(ipnt,i,j,k,2,3) = eta_d*Gyz
                Sflux(ipnt,i,j,k,2,1) = Sflux(ipnt,i,j,k,1,2)
                Sflux(ipnt,i,j,k,3,1) = Sflux(ipnt,i,j,k,1,3)
                Sflux(ipnt,i,j,k,3,2) = Sflux(ipnt,i,j,k,2,3)

                trG   = (Gxx + Gyy + Gzz)
                trGd3 = c1d3*trG
                trG_zeta = zeta_d*trG

                Sflux(ipnt,i,j,k,1,1) = eta_d*(Gxx - trGd3) + trG_zeta
                Sflux(ipnt,i,j,k,2,2) = eta_d*(Gyy - trGd3) + trG_zeta
                Sflux(ipnt,i,j,k,3,3) = eta_d*(Gzz - trGd3) + trG_zeta
            enddo
        enddo
        enddo
        enddo
    end subroutine get_region_Sflux
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine get_region_GRM_xyz(GRM_x, GRM_y, GRM_z)
        implicit none
        real, dimension(npts_llns, 1:nx1, ny,    nz,    3,3), intent(out) :: GRM_x
        real, dimension(npts_llns, nx,    1:ny1, nz,    3,3), intent(out) :: GRM_y
        real, dimension(npts_llns, nx,    ny,    1:nz1, 3,3), intent(out) :: GRM_z
        real, allocatable :: grn(:)
        integer vsl_ndim

        vsl_ndim = npts_llns*nx1*ny*nz*9
        allocate(grn(vsl_ndim))
        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)
        GRM_x = reshape( grn, (/ npts_llns, nx1, ny, nz, 3, 3 /) )
        deallocate(grn)

        vsl_ndim = npts_llns*nx*ny1*nz*9
        allocate(grn(vsl_ndim))
        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)
        GRM_y = reshape( grn, (/ npts_llns, nx, ny1, nz, 3, 3 /) )
        deallocate(grn)

        vsl_ndim = npts_llns*nx*ny*nz1*9
        allocate(grn(vsl_ndim))
        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)
        GRM_z = reshape( grn, (/ npts_llns, nx, ny, nz1, 3, 3 /) )
        deallocate(grn)
    end subroutine get_region_GRM_xyz
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine get_region_Sflux_xyz(Sflux_x, Sflux_y, Sflux_z, dt)
        implicit none
        real, dimension(npts_llns, 1:nx1, ny,    nz,    3,3), intent(out) :: Sflux_x
        real, dimension(npts_llns, nx,    1:ny1, nz,    3,3), intent(out) :: Sflux_y
        real, dimension(npts_llns, nx,    ny,    1:nz1, 3,3), intent(out) :: Sflux_z
        real, intent(inout) :: dt
        real, dimension(npts_llns, 1:nx1, ny,    nz,    3,3) :: GRM_x
        real, dimension(npts_llns, nx,    1:ny1, nz,    3,3) :: GRM_y
        real, dimension(npts_llns, nx,    ny,    1:nz1, 3,3) :: GRM_z
        integer :: i,j,k,ipnt
        real :: Gxx,Gyy,Gzz,Gxy,Gxz,Gyz
        real :: eta_d,zeta_d
        real :: trG,trGd3,trG_zeta

        dvdti = npts_llns*(dxi*dyi*dzi)/dt
        sqrt_dvdti = sqrt(dvdti)
        eta_d  = eta_sd  * sqrt_dvdti
        zeta_d = zeta_sd * sqrt_dvdti
        call get_region_GRM_xyz(GRM_x, GRM_y, GRM_z)

        do k=1,nz
        do j=1,ny
        do i=1,nx1
            do ipnt=1,npts_llns
                Gxx = GRM_x(ipnt,i,j,k,1,1)
                Gyy = GRM_x(ipnt,i,j,k,2,2)
                Gzz = GRM_x(ipnt,i,j,k,3,3)

                Gxy = sqrt2i*( GRM_x(ipnt,i,j,k,1,2) + GRM_x(ipnt,i,j,k,2,1) )
                Gxz = sqrt2i*( GRM_x(ipnt,i,j,k,1,3) + GRM_x(ipnt,i,j,k,3,1) )
                Gyz = sqrt2i*( GRM_x(ipnt,i,j,k,2,3) + GRM_x(ipnt,i,j,k,3,2) )

                Sflux_x(ipnt,i,j,k,1,2) = eta_d*Gxy
                Sflux_x(ipnt,i,j,k,1,3) = eta_d*Gxz
                Sflux_x(ipnt,i,j,k,2,3) = eta_d*Gyz
                Sflux_x(ipnt,i,j,k,2,1) = Sflux_x(ipnt,i,j,k,1,2)
                Sflux_x(ipnt,i,j,k,3,1) = Sflux_x(ipnt,i,j,k,1,3)
                Sflux_x(ipnt,i,j,k,3,2) = Sflux_x(ipnt,i,j,k,2,3)

                trG   = (Gxx + Gyy + Gzz)
                trGd3 = c1d3*trG
                trG_zeta = zeta_d*trG

                Sflux_x(ipnt,i,j,k,1,1) = eta_d*(Gxx - trGd3) + trG_zeta
                Sflux_x(ipnt,i,j,k,2,2) = eta_d*(Gyy - trGd3) + trG_zeta
                Sflux_x(ipnt,i,j,k,3,3) = eta_d*(Gzz - trGd3) + trG_zeta
            enddo
        enddo
        enddo
        enddo

        do k=1,nz
        do j=1,ny1
        do i=1,nx
            do ipnt=1,npts_llns
                Gxx = GRM_y(ipnt,i,j,k,1,1)
                Gyy = GRM_y(ipnt,i,j,k,2,2)
                Gzz = GRM_y(ipnt,i,j,k,3,3)

                Gxy = sqrt2i*( GRM_y(ipnt,i,j,k,1,2) + GRM_y(ipnt,i,j,k,2,1) )
                Gxz = sqrt2i*( GRM_y(ipnt,i,j,k,1,3) + GRM_y(ipnt,i,j,k,3,1) )
                Gyz = sqrt2i*( GRM_y(ipnt,i,j,k,2,3) + GRM_y(ipnt,i,j,k,3,2) )

                Sflux_y(ipnt,i,j,k,1,2) = eta_d*Gxy
                Sflux_y(ipnt,i,j,k,1,3) = eta_d*Gxz
                Sflux_y(ipnt,i,j,k,2,3) = eta_d*Gyz
                Sflux_y(ipnt,i,j,k,2,1) = Sflux_y(ipnt,i,j,k,1,2)
                Sflux_y(ipnt,i,j,k,3,1) = Sflux_y(ipnt,i,j,k,1,3)
                Sflux_y(ipnt,i,j,k,3,2) = Sflux_y(ipnt,i,j,k,2,3)

                trG   = (Gxx + Gyy + Gzz)
                trGd3 = c1d3*trG
                trG_zeta = zeta_d*trG

                Sflux_y(ipnt,i,j,k,1,1) = eta_d*(Gxx - trGd3) + trG_zeta
                Sflux_y(ipnt,i,j,k,2,2) = eta_d*(Gyy - trGd3) + trG_zeta
                Sflux_y(ipnt,i,j,k,3,3) = eta_d*(Gzz - trGd3) + trG_zeta
            enddo
        enddo
        enddo
        enddo

        do k=1,nz1
        do j=1,ny
        do i=1,nx
            do ipnt=1,npts_llns
                Gxx = GRM_z(ipnt,i,j,k,1,1)
                Gyy = GRM_z(ipnt,i,j,k,2,2)
                Gzz = GRM_z(ipnt,i,j,k,3,3)

                Gxy = sqrt2i*( GRM_z(ipnt,i,j,k,1,2) + GRM_z(ipnt,i,j,k,2,1) )
                Gxz = sqrt2i*( GRM_z(ipnt,i,j,k,1,3) + GRM_z(ipnt,i,j,k,3,1) )
                Gyz = sqrt2i*( GRM_z(ipnt,i,j,k,2,3) + GRM_z(ipnt,i,j,k,3,2) )

                Sflux_z(ipnt,i,j,k,1,2) = eta_d*Gxy
                Sflux_z(ipnt,i,j,k,1,3) = eta_d*Gxz
                Sflux_z(ipnt,i,j,k,2,3) = eta_d*Gyz
                Sflux_z(ipnt,i,j,k,2,1) = Sflux_z(ipnt,i,j,k,1,2)
                Sflux_z(ipnt,i,j,k,3,1) = Sflux_z(ipnt,i,j,k,1,3)
                Sflux_z(ipnt,i,j,k,3,2) = Sflux_z(ipnt,i,j,k,2,3)

                trG   = (Gxx + Gyy + Gzz)
                trGd3 = c1d3*trG
                trG_zeta = zeta_d*trG

                Sflux_z(ipnt,i,j,k,1,1) = eta_d*(Gxx - trGd3) + trG_zeta
                Sflux_z(ipnt,i,j,k,2,2) = eta_d*(Gyy - trGd3) + trG_zeta
                Sflux_z(ipnt,i,j,k,3,3) = eta_d*(Gzz - trGd3) + trG_zeta
            enddo
        enddo
        enddo
        enddo

    end subroutine get_region_Sflux_xyz
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    ! subroutine get_GRM(GRMpnts_r)
    !     implicit none
    !     integer ife,b,e,vsl_ndim
    !     real, dimension(npts_llns,3,3) :: GRMpnts_r
    !
    !     real, dimension(9*npts_llns) :: grn
    !     ! real, allocatable :: grn(:)
    !
    !     vsl_ndim = 9*npts_llns
    !     ! allocate(grn(vsl_ndim))
    !
    !     vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)
    !
    !     ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
    !     do ife = 1,npts_llns
    !         e = 9*ife
    !         b = e - 8
    !         GRMpnts_r(ife,1,:) = grn(b:b+2)
    !         GRMpnts_r(ife,2,:) = grn(b+3:b+5)
    !         GRMpnts_r(ife,3,:) = grn(b+6:e)
    !     end do
    !
    ! end subroutine get_GRM
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    ! subroutine random_stresses_pnts_r(Spnts_r)
    !     implicit none
    !     integer ife
    !     real, dimension(npts_llns,3,3) :: GRMpnts_r, Spnts_r
    !     real Gxx,Gyy,Gzz,Gxy,Gxz,Gyz
    !     real eta_d,zeta_d
    !     real trG,trGd3,trG_zeta
    !
    !     eta_d = eta_sd * sqrt_dvdti
    !     zeta_d = zeta_sd * sqrt_dvdti
    !     call get_GRM(GRMpnts_r)
    !
    !     ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
    !     do ife = 1,npts_llns
    !         Gxx = GRMpnts_r(ife,1,1)
    !         Gyy = GRMpnts_r(ife,2,2)
    !         Gzz = GRMpnts_r(ife,3,3)
    !
    !         Gxy = sqrt2i*( GRMpnts_r(ife,1,2) + GRMpnts_r(ife,2,1) )
    !         Gxz = sqrt2i*( GRMpnts_r(ife,1,3) + GRMpnts_r(ife,3,1) )
    !         Gyz = sqrt2i*( GRMpnts_r(ife,2,3) + GRMpnts_r(ife,3,2) )
    !
    !         Spnts_r(ife,1,2) = eta_d*Gxy
    !         Spnts_r(ife,1,3) = eta_d*Gxz
    !         Spnts_r(ife,2,3) = eta_d*Gyz
    !         Spnts_r(ife,2,1) = Spnts_r(ife,1,2)
    !         Spnts_r(ife,3,1) = Spnts_r(ife,1,3)
    !         Spnts_r(ife,3,2) = Spnts_r(ife,2,3)
    !
    !         trG = (Gxx + Gyy + Gzz)
    !         trGd3 = c1d3*trG
    !         trG_zeta = zeta_d*trG
    !
    !         Spnts_r(ife,1,1) = eta_d*(Gxx - trGd3) + trG_zeta
    !         Spnts_r(ife,2,2) = eta_d*(Gyy - trGd3) + trG_zeta
    !         Spnts_r(ife,3,3) = eta_d*(Gzz - trGd3) + trG_zeta
    !     end do
    !
    ! end subroutine random_stresses_pnts_r
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine get_GRM_symmetric(GRMpnts_r)
        implicit none
        real, dimension(npts_llns,3,3), intent(out) :: GRMpnts_r
        real, dimension(6*npts_llns) :: grn
        integer :: ife,b,e,vsl_ndim
        ! real, allocatable :: grn(:)

        vsl_ndim = 6*npts_llns
        ! allocate(grn(vsl_ndim))

        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)

        ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
        do ife = 1,npts_llns
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
    subroutine random_stresses_pnts_r(Spnts_r)
        implicit none
        real, dimension(npts_llns,3,3), intent(out) :: Spnts_r
        real, dimension(npts_llns,3,3) :: GRMpnts_r
        integer ife
        real Gxx,Gyy,Gzz,Gxy,Gxz,Gyz, eta_d,zeta_d, trG,trGd3,trG_zeta

        dvdti = npg*(dxi*dyi*dzi)/dt
        sqrt_dvdti = sqrt(dvdti)
        eta_d  = eta_sd * sqrt_dvdti
        zeta_d = zeta_sd * sqrt_dvdti
        call get_GRM_symmetric(GRMpnts_r)  ! NOTE: using nptns = 1 (FV approach)

        ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
        do ife = 1,npts_llns
            Gxx = GRMpnts_r(ife,1,1)
            Gyy = GRMpnts_r(ife,2,2)
            Gzz = GRMpnts_r(ife,3,3)

            Gxy = GRMpnts_r(ife,1,2)
            Gxz = GRMpnts_r(ife,1,3)
            Gyz = GRMpnts_r(ife,2,3)

            trG   = (Gxx + Gyy + Gzz)
            trGd3 = c1d3*trG
            trG_zeta = zeta_d*trG

            Spnts_r(ife,1,1) = eta_d*(Gxx - trGd3) + trG_zeta
            Spnts_r(ife,2,2) = eta_d*(Gyy - trGd3) + trG_zeta
            Spnts_r(ife,3,3) = eta_d*(Gzz - trGd3) + trG_zeta

            Spnts_r(ife,1,2) = eta_d*Gxy
            Spnts_r(ife,1,3) = eta_d*Gxz
            Spnts_r(ife,2,3) = eta_d*Gyz
            Spnts_r(ife,2,1) = Spnts_r(ife,1,2)
            Spnts_r(ife,3,1) = Spnts_r(ife,1,3)
            Spnts_r(ife,3,2) = Spnts_r(ife,2,3)
        end do

    end subroutine random_stresses_pnts_r
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine random_heatflux_pnts_r(Hpnts_r)
        implicit none
        real, dimension(npts_llns,3), intent(out) :: Hpnts_r
        real, dimension(npts_llns,3) :: GRVpnts_r
        integer ife,b,e,vsl_ndim
        real Gx,Gy,Gz,kappa_d
        real, allocatable :: grn(:)

        vsl_ndim = npts_llns*3
        allocate(grn(vsl_ndim))

        kappa_d = kappa_sd * sqrt_dvdti

        vsl_errcode = vsRngGaussian(vsl_method, vsl_stream, vsl_ndim, grn, vsl_mean, vsl_sigma)

        ! NOTE: There's probably a better way to do a reshape (w/o using Fortran 2003/8)
        do ife = 1,npts_llns
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
