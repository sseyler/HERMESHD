!****** BASIS_FUNCS.F90 ******************************************************************
module basis_funcs

!===========================================================================
! Parameters
!------------------------------------------------------------
use params

integer, parameter :: nvtk  = 1    ! (was 2)
integer, parameter :: nvtk2 = nvtk*nvtk
integer, parameter :: nvtk3 = nvtk*nvtk*nvtk
!===========================================================================


!===========================================================================
! Definitions
!------------------------------------------------------------
real, dimension(nvtk3,nbastot) :: bfvtk, bfvtk_dx, bfvtk_dy, bfvtk_dz

! TODO: only in init, innerintegral, source_calc
real, dimension(nbastot) :: cbasis

! TODO: only in init + innerintegral
integer ibitri  ! set by set_cbasis_3D using chosen value of iquad

! TODO: only in set_weights_3D
real wgt1d(5)    ! wgt1d: quadrature weights for 1-D integration

! TODO: only in init, set_weights_3D, glflux
real wgt2d(30)   ! wgt2d: quadrature weights for 2-D integration

! TODO: only in init, set_weights_3D, innerintegral
real wgt3d(100)  ! wgt3d: quadrature weights for 3-D integration

! TODO: only in limiter, set_face_vals_3D
real, dimension(nslim,nbastot) :: bf_faces

! TODO: only in:
!   * initialize.f90 (setup)
!   * set_face_vals_3D, set_internal_vals_3D
!   * source_calc, innerintegral
real, dimension(npg,nbastot) :: bfvals_int

! TODO: only in set_internal_vals_3D, set_face_vals_3D
real xquad(20)

real, dimension(npg,nbastot) :: bval_int_wgt  ! used in source_calc
real, dimension(nface,2,nbastot) :: wgtbf_xmp, wgtbf_ymp, wgtbf_zmp  ! used in glflux

! TODO: only in init, flux_cal, prepare_exchange, set_face_vals_3D
real, dimension(nface,nbastot) :: bfvals_zp, bfvals_zm
real, dimension(nface,nbastot) :: bfvals_yp, bfvals_ym
real, dimension(nface,nbastot) :: bfvals_xp, bfvals_xm

real :: xgrid(20) ! used in set_vtk_vals_3D

! These are used in test_basis_3D
real, dimension(nbastot,nbastot) :: cell_int,xp_int,xm_int,yp_int,ym_int,zp_int,zm_int
real, dimension(nbastot,nbastot) :: cell_int0,xp_int0,xm_int0,yp_int0,ym_int0,zp_int0,zm_int0
real, dimension(nbastot) :: cbas_xp,cbas_xm,cbas_yp,cbas_ym,cbas_zp,cbas_zm
integer, dimension(nbastot) :: ibas_x, ibas_y, ibas_z
!===========================================================================

contains

    subroutine set_bfvals_3D
        ! Defines local basis function values and weights for 1, 2, or 3-point Gaussian quadrature.
        ! Basis functions are evaluated in cell interior and on cell faces.
        implicit none

        call set_vtk_vals_3D()    ! Define basis function values on a 3D grid INTERNAL to the cell.
        call set_internal_vals_3D()     ! Define basis function values at quadrature points INTERNAL to cell.
        call set_face_vals_3D()   ! Define local basis function values at quadrature points on a cell face.
        call set_weights_3D()     ! Define weights for integral approximation using Gaussian quadrature.

        if (test_basis .and. iam == print_mpi) then
    	   print *,'Testing basis...'
    	   call test_basis_3D()
    	   print *,'Done testing basis...'
    	end if

    end subroutine set_bfvals_3D

!----------------------------------------------------------
    subroutine set_vtk_vals_3D
        ! Define basis function values on a rectangular 3D grid INTERNAL to the cell.
        ! For use in output_vtk().

        implicit none
        integer ixyz,i,igrid

        dxvtk = 1./(nvtk*dxi)
        dyvtk = 1./(nvtk*dyi)
        dzvtk = 1./(nvtk*dzi)

        if (nvtk .eq. 1) xgrid(1) = 0.

        if (nvtk .eq. 2) then
            xgrid(1) = -0.5
            xgrid(2) = 0.5
        end if
        if (nvtk .eq. 3) then
            xgrid(1) = -2./3.
            xgrid(2) = 0.
            xgrid(3) = 2./3.
        end if
        if (nvtk .eq. 4) then
            xgrid(1) = -3./4.
            xgrid(2) = -1./4.
            xgrid(3) = 1./4.
            xgrid(4) = 3./4.
        end if


        bfvtk(1,1) = 1.        ! basis function = 1
        bfvtk(1,kx) = 0.        ! basis function = x
        bfvtk(1,ky) = 0.        ! basis function = y
        bfvtk(1,kz) = 0.        ! basis function = z

        bfvtk(1:nvtk3,1) = 1.        ! basis function = 1
        do i=1,nvtk
            bfvtk((i-1)*nvtk+1:i*nvtk,ky) = xgrid(i)        ! basis function = y
            bfvtk((i-1)*nvtk+1:i*nvtk,kz) = xgrid(1:nvtk)        ! basis function = z
        end do
        do i=1,nvtk
            bfvtk((i-1)*nvtk2+1:i*nvtk2,kx) = xgrid(i)        ! basis function = x
            bfvtk((i-1)*nvtk2+1:i*nvtk2,ky) = bfvtk(1:nvtk2,ky)        ! basis function = y
            bfvtk((i-1)*nvtk2+1:i*nvtk2,kz) = bfvtk(1:nvtk2,kz)        ! basis function = z
        end do

        do i=0,2
            bfvtk(1:nvtk3,kxx+i) = 1.5*bfvtk(1:nvtk3,kx+i)**2 - 0.5    ! basis function = P_2(s)
            ! bfvtk(1:nvtk3,kxxx+i) = 2.5*bfvtk(1:nvtk3,kx+i)**3 - 1.5*bfvtk(1:nvtk3,kx+i)   ! basis function = P3(s)
        end do

        bfvtk(1:nvtk3,kyz) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kz)        ! basis function = yz
        bfvtk(1:nvtk3,kzx) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kx)        ! basis function = zx
        bfvtk(1:nvtk3,kxy) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,ky)        ! basis function = xy
        bfvtk(1:nvtk3,kxyz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kz)     ! basis function = xyz

        bfvtk(1:nvtk3,kyyz) = bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kz)     ! basis function = P2(y)z
        bfvtk(1:nvtk3,kyzz) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kzz)     ! basis function = P2(z)y
        bfvtk(1:nvtk3,kzzx) = bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kx)     ! basis function = P2(z)x
        bfvtk(1:nvtk3,kzxx) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kxx)     ! basis function = P2(x)z
        bfvtk(1:nvtk3,kxxy) = bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,ky)     ! basis function = P2(x)y
        bfvtk(1:nvtk3,kxyy) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)     ! basis function = P2(y)x
        bfvtk(1:nvtk3,kyyzz) = bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)     ! basis function = P_2(y)P_2(z)
        bfvtk(1:nvtk3,kzzxx) = bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kxx)     ! basis function = P_2(z)P_2(x)
        bfvtk(1:nvtk3,kxxyy) = bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,kyy)     ! basis function = P_2(x)P_2(y)
        bfvtk(1:nvtk3,kyzxx) = bfvtk(1:nvtk3,kyz)*bfvtk(1:nvtk3,kxx)     ! basis function = yz P_2(x)
        bfvtk(1:nvtk3,kzxyy) = bfvtk(1:nvtk3,kzx)*bfvtk(1:nvtk3,kyy)     ! basis function = zx P_2(y)
        bfvtk(1:nvtk3,kxyzz) = bfvtk(1:nvtk3,kxy)*bfvtk(1:nvtk3,kzz)     ! basis function = xy P_2(z)
        bfvtk(1:nvtk3,kxyyzz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)   ! basis function = x P_2(y)P_2(z)
        bfvtk(1:nvtk3,kyzzxx) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kzz)*bfvtk(1:nvtk3,kxx)   ! basis function = y P_2(z)P_2(x)
        bfvtk(1:nvtk3,kzxxyy) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kxx)*bfvtk(1:nvtk3,kyy)   ! basis function = z P_2(x)P_2(y)
        bfvtk(1:nvtk3,kxxyyzz) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,kyy)*bfvtk(1:nvtk3,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

        do igrid=1,nvtk3
            bfvtk_dx(igrid,1) = 0.
            bfvtk_dx(igrid,kx) = 1.
            bfvtk_dx(igrid,ky:kz) = 0.
            bfvtk_dx(igrid,kyz) = 0.
            bfvtk_dx(igrid,kzx) = bfvtk(igrid,kz)
            bfvtk_dx(igrid,kxy) = bfvtk(igrid,ky)
            bfvtk_dx(igrid,kxyz) = bfvtk(igrid,kyz)

            bfvtk_dx(igrid,kxx) = 3.*bfvtk(igrid,kx)
            bfvtk_dx(igrid,kyy:kzz) = 0.
            bfvtk_dx(igrid,kyyz) = 0.
            bfvtk_dx(igrid,kyzz) = 0.
            bfvtk_dx(igrid,kzzx) = bfvtk(igrid,kzz)
            bfvtk_dx(igrid,kzxx) = bfvtk(igrid,kz)*3.*bfvtk(igrid,kx)
            bfvtk_dx(igrid,kxxy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,ky)
            bfvtk_dx(igrid,kxyy) = bfvtk(igrid,kyy)
            bfvtk_dx(igrid,kyyzz) = 0.
            bfvtk_dx(igrid,kzzxx) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kzz)
            bfvtk_dx(igrid,kxxyy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyy)
            bfvtk_dx(igrid,kyzxx) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyz)
            bfvtk_dx(igrid,kzxyy) = bfvtk(igrid,kz)*bfvtk(igrid,kyy)
            bfvtk_dx(igrid,kxyzz) = bfvtk(igrid,ky)*bfvtk(igrid,kzz)
            bfvtk_dx(igrid,kxyyzz) = bfvtk(igrid,kyyzz)
            bfvtk_dx(igrid,kyzzxx) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyzz)
            bfvtk_dx(igrid,kzxxyy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyyz)
            bfvtk_dx(igrid,kxxyyzz) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyy)*bfvtk(igrid,kzz)

            bfvtk_dy(igrid,1) = 0.
            bfvtk_dy(igrid,ky) = 1.
            bfvtk_dy(igrid,kx) = 0.
            bfvtk_dy(igrid,kz) = 0.
            bfvtk_dy(igrid,kyz) = bfvtk(igrid,kz)
            bfvtk_dy(igrid,kzx) = 0.
            bfvtk_dy(igrid,kxy) = bfvtk(igrid,kx)
            bfvtk_dy(igrid,kxyz) = bfvtk(igrid,kzx)

            bfvtk_dy(igrid,kyy) = 3.*bfvtk(igrid,ky)
            bfvtk_dy(igrid,kxx) = 0.
            bfvtk_dy(igrid,kzz) = 0.
            bfvtk_dy(igrid,kyzz) = bfvtk(igrid,kzz)
            bfvtk_dy(igrid,kyyz) = bfvtk(igrid,kz)*3.*bfvtk(igrid,ky)
            bfvtk_dy(igrid,kzzx) = 0.
            bfvtk_dy(igrid,kzxx) = 0.
            bfvtk_dy(igrid,kxxy) = bfvtk(igrid,kxx)
            bfvtk_dy(igrid,kxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kx)
            bfvtk_dy(igrid,kyyzz) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzz)
            bfvtk_dy(igrid,kzzxx) = 0.
            bfvtk_dy(igrid,kxxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kxx)
            bfvtk_dy(igrid,kyzxx) = bfvtk(igrid,kz)*bfvtk(igrid,kxx)
            bfvtk_dy(igrid,kzxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzx)
            bfvtk_dy(igrid,kxyzz) = bfvtk(igrid,kx)*bfvtk(igrid,kzz)
            bfvtk_dy(igrid,kxyyzz) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzzx)
            bfvtk_dy(igrid,kyzzxx) = bfvtk(igrid,kzzxx)
            bfvtk_dy(igrid,kzxxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzxx)
            bfvtk_dy(igrid,kxxyyzz) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kzz)*bfvtk(igrid,kxx)


            bfvtk_dz(igrid,1) = 0.
            bfvtk_dz(igrid,kz) = 1.
            bfvtk_dz(igrid,kx) = 0.
            bfvtk_dz(igrid,ky) = 0.
            bfvtk_dz(igrid,kyz) = bfvtk(igrid,ky)
            bfvtk_dz(igrid,kzx) = bfvtk(igrid,kx)
            bfvtk_dz(igrid,kxy) = 0.
            bfvtk_dz(igrid,kxyz) = bfvtk(igrid,kxy)

            bfvtk_dz(igrid,kzz) = 3.*bfvtk(igrid,kz)
            bfvtk_dz(igrid,kxx) = 0.
            bfvtk_dz(igrid,kyy) = 0.

            bfvtk_dz(igrid,kyzz) = bfvtk(igrid,ky)*3.*bfvtk(igrid,kz)
            bfvtk_dz(igrid,kyyz) = bfvtk(igrid,kyy)
            bfvtk_dz(igrid,kzzx) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kx)
            bfvtk_dz(igrid,kzxx) = bfvtk(igrid,kxx)
            bfvtk_dz(igrid,kxxy) = 0.
            bfvtk_dz(igrid,kxyy) = 0.
            bfvtk_dz(igrid,kyyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kyy)
            bfvtk_dz(igrid,kzzxx) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxx)
            bfvtk_dz(igrid,kxxyy) = 0.
            bfvtk_dz(igrid,kyzxx) = bfvtk(igrid,ky)*bfvtk(igrid,kxx)
            bfvtk_dz(igrid,kzxyy) = bfvtk(igrid,kx)*bfvtk(igrid,kyy)
            bfvtk_dz(igrid,kxyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxy)
            bfvtk_dz(igrid,kxyyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxyy)
            bfvtk_dz(igrid,kyzzxx) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxxy)
            bfvtk_dz(igrid,kzxxyy) = bfvtk(igrid,kxxyy)
            bfvtk_dz(igrid,kxxyyzz) = 3.*bfvtk(igrid,kz)*bfvtk(igrid,kxx)*bfvtk(igrid,kyy)

            ! bfvtk_dx(igrid,kxxx) = 7.5*bfvtk(igrid,kx)**2 - 1.5
            ! bfvtk_dx(igrid,kyyy:kzzz) = 0.
            ! bfvtk_dy(igrid,kyyy) = 7.5*bfvtk(igrid,ky)**2 - 1.5
            ! bfvtk_dy(igrid,kxxx) = 0.
            ! bfvtk_dy(igrid,kzzz) = 0.
            ! bfvtk_dz(igrid,kzzz) = 7.5*bfvtk(igrid,kz)**2 - 1.5
            ! bfvtk_dz(igrid,kxxx) = 0.
            ! bfvtk_dz(igrid,kyyy) = 0.
        end do

    end subroutine set_vtk_vals_3D

!----------------------------------------------------------

!----------------------------------------------------------
    subroutine set_internal_vals_3D
        ! Define basis function values at quadrature points INTERNAL to cell.

        implicit none
        integer ixyz, i
        real c15d5,c1dsq3,xq4p,xq4m

        c15d5 = sqrt(15.)/5.
        c1dsq3 = 1./sqrt(3.)

        xq4p = sqrt(525. + 70.*sqrt(30.))/35.
        xq4m = sqrt(525. - 70.*sqrt(30.))/35.


        if (iquad .eq. 2) then
            xquad(1) = -c1dsq3
            xquad(2) = c1dsq3
        end if
        if (iquad .eq. 3) then
            xquad(1) = -c15d5
            xquad(2) = 0.
            xquad(3) = c15d5
        end if
        if (iquad .eq. 4) then
            xquad(1) = -xq4p
            xquad(2) = -xq4m
            xquad(3) = xq4m
            xquad(4) = xq4p
        end if

        if (iquad .eq. 1) then            ! 2-point Gaussian quadrature
            bfvals_int(1,1) = 1.        ! basis function = 1
            bfvals_int(1,kx) = 0.        ! basis function = x
            bfvals_int(1,ky) = 0.        ! basis function = y
            bfvals_int(1,kz) = 0.        ! basis function = z
        end if


        if (iquad .gt. 1) then
            bfvals_int(1:npg,1) = 1.        ! basis function = 1
            do i=1,nedge
                bfvals_int((i-1)*nedge+1:i*nedge,ky) = xquad(i)        ! basis function = y
                bfvals_int((i-1)*nedge+1:i*nedge,kz) = xquad(1:nedge)        ! basis function = z
            end do
            do i=1,nedge
                bfvals_int((i-1)*nface+1:i*nface,kx) = xquad(i)        ! basis function = x
                bfvals_int((i-1)*nface+1:i*nface,ky) = bfvals_int(1:nface,ky)        ! basis function = y
                bfvals_int((i-1)*nface+1:i*nface,kz) = bfvals_int(1:nface,kz)        ! basis function = z
            end do
        end if

        do i=0,2
            bfvals_int(1:npg,kxx+i) = 1.5*bfvals_int(1:npg,kx+i)**2 - 0.5    ! basis function = P_2(s)
            ! bfvals_int(1:npg,kxxx+i) = 2.5*bfvals_int(1:npg,kx+i)**3 - 1.5*bfvals_int(1:npg,kx+i)   ! basis function = P_3(s)
        end do

        bfvals_int(1:npg,kyz) = bfvals_int(1:npg,ky)*bfvals_int(1:npg,kz)        ! basis function = yz
        bfvals_int(1:npg,kzx) = bfvals_int(1:npg,kz)*bfvals_int(1:npg,kx)        ! basis function = zx
        bfvals_int(1:npg,kxy) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,ky)        ! basis function = xy
        bfvals_int(1:npg,kxyz) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,ky)*bfvals_int(1:npg,kz)     ! basis function = xyz

        bfvals_int(1:npg,kyyz) = bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kz)     ! basis function = P_2(y)z
        bfvals_int(1:npg,kyzz) = bfvals_int(1:npg,ky)*bfvals_int(1:npg,kzz)     ! basis function = P_2(z)y
        bfvals_int(1:npg,kzzx) = bfvals_int(1:npg,kzz)*bfvals_int(1:npg,kx)     ! basis function = P_2(z)x
        bfvals_int(1:npg,kzxx) = bfvals_int(1:npg,kz)*bfvals_int(1:npg,kxx)     ! basis function = P_2(x)z
        bfvals_int(1:npg,kxxy) = bfvals_int(1:npg,kxx)*bfvals_int(1:npg,ky)     ! basis function = P_2(x)y
        bfvals_int(1:npg,kxyy) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,kyy)     ! basis function = P_2(y)x
        bfvals_int(1:npg,kyyzz) = bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kzz)     ! basis function = P_2(y)P_2(z)
        bfvals_int(1:npg,kzzxx) = bfvals_int(1:npg,kzz)*bfvals_int(1:npg,kxx)     ! basis function = P_2(z)P_2(x)
        bfvals_int(1:npg,kxxyy) = bfvals_int(1:npg,kxx)*bfvals_int(1:npg,kyy)     ! basis function = P_2(x)P_2(y)
        bfvals_int(1:npg,kyzxx) = bfvals_int(1:npg,kyz)*bfvals_int(1:npg,kxx)     ! basis function = yz P_2(x)
        bfvals_int(1:npg,kzxyy) = bfvals_int(1:npg,kzx)*bfvals_int(1:npg,kyy)     ! basis function = zx P_2(y)
        bfvals_int(1:npg,kxyzz) = bfvals_int(1:npg,kxy)*bfvals_int(1:npg,kzz)     ! basis function = xy P_2(z)
        bfvals_int(1:npg,kxyyzz) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kzz)   ! basis function = x P_2(y)P_2(z)
        bfvals_int(1:npg,kyzzxx) = bfvals_int(1:npg,ky)*bfvals_int(1:npg,kzz)*bfvals_int(1:npg,kxx)   ! basis function = y P_2(z)P_2(x)
        bfvals_int(1:npg,kzxxyy) = bfvals_int(1:npg,kz)*bfvals_int(1:npg,kxx)*bfvals_int(1:npg,kyy)   ! basis function = z P_2(x)P_2(y)
        bfvals_int(1:npg,kxxyyzz) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

    end subroutine set_internal_vals_3D

!----------------------------------------------------------
    subroutine set_face_vals_3D
        ! Define local basis function values at quadrature points on a cell face.
        ! Used in flux_cal() and prepare_exchange() for interpolating from cell center onto cell face.
        ! Also used in glflux().

        implicit none
        integer ixyz,i
        real c15d5,c1dsq3,xq4p,xq4m

        c15d5 = sqrt(15.)/5.
        c1dsq3 = 1./sqrt(3.)

        xq4p = sqrt(525. + 70.*sqrt(30.))/35.
        xq4m = sqrt(525. - 70.*sqrt(30.))/35.

        if (iquad .eq. 2) then
            xquad(1) = -c1dsq3
            xquad(2) = c1dsq3
        end if
        if (iquad .eq. 3) then
            xquad(1) = -c15d5
            xquad(2) = 0.
            xquad(3) = c15d5
        end if
        if (iquad .eq. 4) then
            xquad(1) = -xq4p
            xquad(2) = -xq4m
            xquad(3) = xq4m
            xquad(4) = xq4p
        end if

        ! bfvals_rsp:  positive rs-face
        ! bfvals_rsm:  negative rs-face
        ! bfvals_rsp(,1):  value=1 on positive rs-face
        ! bfvals_rsp(,kx):  x-values on positive rs-face
        ! bfvals_rsp(,ky):  y-values on positive rs-face
        ! bfvals_rsp(,kz):  z-values on positive rs-face
        ! bfvals_rsp(,kxx):  P_2(x)-values on positive rs-face
        ! bfvals_rsp(,kyy):  P_2(y)-values on positive rs-face
        ! bfvals_rsp(,kzz):  P_2(z)-values on positive rs-face
        ! bfvals_rsp(,kyz):  yz-values on positive rs-face
        ! bfvals_rsp(,kzx):  zx-values on positive rs-face
        ! bfvals_rsp(,kxy):  xy-values on positive rs-face

        bfvals_zp(1:nface,1) = 1.
        if (iquad .eq. 1) then
            bfvals_zp(1,kx) = 0.
            bfvals_zp(1,ky) = 0.
        end if

        if (iquad .gt. 1) then
            do i=1,nedge
                bfvals_zp((i-1)*nedge+1:i*nedge,kx) = xquad(1:nedge)
                bfvals_zp((i-1)*nedge+1:i*nedge,ky) = xquad(i)
            end do
        end if

        bfvals_zp(1:nface,kz) = 1.
        bfvals_zp(1:nface,kyz) = bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kz)
        bfvals_zp(1:nface,kzx) = bfvals_zp(1:nface,kz)*bfvals_zp(1:nface,kx)
        bfvals_zp(1:nface,kxy) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,ky)
        bfvals_zp(1:nface,kxyz) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kz)     ! basis function = xyz

        bfvals_yp(1:nface,1) = 1.
        bfvals_yp(1:nface,kx) = bfvals_zp(1:nface,ky)
        bfvals_yp(1:nface,ky) = 1.
        bfvals_yp(1:nface,kz) = bfvals_zp(1:nface,kx)
        bfvals_yp(1:nface,kyz) = bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kz)
        bfvals_yp(1:nface,kzx) = bfvals_yp(1:nface,kz)*bfvals_yp(1:nface,kx)
        bfvals_yp(1:nface,kxy) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,ky)
        bfvals_yp(1:nface,kxyz) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kz)     ! basis function = xyz

        bfvals_xp(1:nface,1) = 1.
        bfvals_xp(1:nface,kx) = 1.
        bfvals_xp(1:nface,ky) = bfvals_zp(1:nface,kx)
        bfvals_xp(1:nface,kz) = bfvals_zp(1:nface,ky)
        bfvals_xp(1:nface,kyz) = bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kz)
        bfvals_xp(1:nface,kzx) = bfvals_xp(1:nface,kz)*bfvals_xp(1:nface,kx)
        bfvals_xp(1:nface,kxy) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,ky)
        bfvals_xp(1:nface,kxyz) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kz)     ! basis function = xyz

        do i=1,3
            bfvals_xp(1:nface,kxx+i-1) = 1.5*bfvals_xp(1:nface,kx+i-1)**2 - 0.5
            bfvals_yp(1:nface,kxx+i-1) = 1.5*bfvals_yp(1:nface,kx+i-1)**2 - 0.5
            bfvals_zp(1:nface,kxx+i-1) = 1.5*bfvals_zp(1:nface,kx+i-1)**2 - 0.5
            ! bfvals_xp(1:nface,kxxx+i-1) = 2.5*bfvals_xp(1:nface,kx+i-1)**3 - 1.5*bfvals_xp(1:nface,kx+i-1)
            ! bfvals_yp(1:nface,kxxx+i-1) = 2.5*bfvals_yp(1:nface,kx+i-1)**3 - 1.5*bfvals_yp(1:nface,kx+i-1)
            ! bfvals_zp(1:nface,kxxx+i-1) = 2.5*bfvals_zp(1:nface,kx+i-1)**3 - 1.5*bfvals_zp(1:nface,kx+i-1)
        end do

        bfvals_xp(1:nface,kyyz) = bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kz)     ! basis function = P_2(y)z
        bfvals_xp(1:nface,kyzz) = bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kzz)     ! basis function = P_2(z)y
        bfvals_xp(1:nface,kzzx) = bfvals_xp(1:nface,kzz)*bfvals_xp(1:nface,kx)     ! basis function = P_2(z)x
        bfvals_xp(1:nface,kzxx) = bfvals_xp(1:nface,kz)*bfvals_xp(1:nface,kxx)     ! basis function = P_2(x)z
        bfvals_xp(1:nface,kxxy) = bfvals_xp(1:nface,kxx)*bfvals_xp(1:nface,ky)     ! basis function = P_2(x)y
        bfvals_xp(1:nface,kxyy) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,kyy)     ! basis function = P_2(y)x
        bfvals_xp(1:nface,kyyzz) = bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kzz)     ! basis function = P_2(y)P_2(z)
        bfvals_xp(1:nface,kzzxx) = bfvals_xp(1:nface,kzz)*bfvals_xp(1:nface,kxx)     ! basis function = P_2(z)P_2(x)
        bfvals_xp(1:nface,kxxyy) = bfvals_xp(1:nface,kxx)*bfvals_xp(1:nface,kyy)     ! basis function = P_2(x)P_2(y)
        bfvals_xp(1:nface,kyzxx) = bfvals_xp(1:nface,kyz)*bfvals_xp(1:nface,kxx)     ! basis function = yz P_2(x)
        bfvals_xp(1:nface,kzxyy) = bfvals_xp(1:nface,kzx)*bfvals_xp(1:nface,kyy)     ! basis function = zx P_2(y)
        bfvals_xp(1:nface,kxyzz) = bfvals_xp(1:nface,kxy)*bfvals_xp(1:nface,kzz)     ! basis function = xy P_2(z)
        bfvals_xp(1:nface,kxyyzz) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kzz)   ! basis function = x P_2(y)P_2(z)
        bfvals_xp(1:nface,kyzzxx) = bfvals_xp(1:nface,ky)*bfvals_xp(1:nface,kzz)*bfvals_xp(1:nface,kxx)   ! basis function = y P_2(z)P_2(x)
        bfvals_xp(1:nface,kzxxyy) = bfvals_xp(1:nface,kz)*bfvals_xp(1:nface,kxx)*bfvals_xp(1:nface,kyy)   ! basis function = z P_2(x)P_2(y)
        bfvals_xp(1:nface,kxxyyzz) = bfvals_xp(1:nface,kx)*bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

        bfvals_xm = bfvals_xp
        bfvals_xm(1:nface,kx) = -bfvals_xp(1:nface,kx)
        bfvals_xm(1:nface,kzx) = -bfvals_xp(1:nface,kzx)
        bfvals_xm(1:nface,kxy) = -bfvals_xp(1:nface,kxy)
        bfvals_xm(1:nface,kxyz) = -bfvals_xp(1:nface,kxyz)
        bfvals_xm(1:nface,kzzx) = -bfvals_xp(1:nface,kzzx)
        bfvals_xm(1:nface,kxyy) = -bfvals_xp(1:nface,kxyy)
        ! bfvals_xm(1:nface,kxxx) = -bfvals_xp(1:nface,kxxx)
        bfvals_xm(1:nface,kzxyy) = -bfvals_xp(1:nface,kzxyy)
        bfvals_xm(1:nface,kxyzz) = -bfvals_xp(1:nface,kxyzz)
        bfvals_xm(1:nface,kxyyzz) = -bfvals_xp(1:nface,kxyyzz)

        bfvals_yp(1:nface,kyyz) = bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kz)     ! basis function = P_2(y)z
        bfvals_yp(1:nface,kyzz) = bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kzz)     ! basis function = P_2(z)y
        bfvals_yp(1:nface,kzzx) = bfvals_yp(1:nface,kzz)*bfvals_yp(1:nface,kx)     ! basis function = P_2(z)x
        bfvals_yp(1:nface,kzxx) = bfvals_yp(1:nface,kz)*bfvals_yp(1:nface,kxx)     ! basis function = P_2(x)z
        bfvals_yp(1:nface,kxxy) = bfvals_yp(1:nface,kxx)*bfvals_yp(1:nface,ky)     ! basis function = P_2(x)y
        bfvals_yp(1:nface,kxyy) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,kyy)     ! basis function = P_2(y)x
        bfvals_yp(1:nface,kyyzz) = bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kzz)     ! basis function = P_2(y)P_2(z)
        bfvals_yp(1:nface,kzzxx) = bfvals_yp(1:nface,kzz)*bfvals_yp(1:nface,kxx)     ! basis function = P_2(z)P_2(x)
        bfvals_yp(1:nface,kxxyy) = bfvals_yp(1:nface,kxx)*bfvals_yp(1:nface,kyy)     ! basis function = P_2(x)P_2(y)
        bfvals_yp(1:nface,kyzxx) = bfvals_yp(1:nface,kyz)*bfvals_yp(1:nface,kxx)     ! basis function = yz P_2(x)
        bfvals_yp(1:nface,kzxyy) = bfvals_yp(1:nface,kzx)*bfvals_yp(1:nface,kyy)     ! basis function = zx P_2(y)
        bfvals_yp(1:nface,kxyzz) = bfvals_yp(1:nface,kxy)*bfvals_yp(1:nface,kzz)     ! basis function = xy P_2(z)
        bfvals_yp(1:nface,kxyyzz) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kzz)   ! basis function = x P_2(y)P_2(z)
        bfvals_yp(1:nface,kyzzxx) = bfvals_yp(1:nface,ky)*bfvals_yp(1:nface,kzz)*bfvals_yp(1:nface,kxx)   ! basis function = y P_2(z)P_2(x)
        bfvals_yp(1:nface,kzxxyy) = bfvals_yp(1:nface,kz)*bfvals_yp(1:nface,kxx)*bfvals_yp(1:nface,kyy)   ! basis function = z P_2(x)P_2(y)
        bfvals_yp(1:nface,kxxyyzz) = bfvals_yp(1:nface,kx)*bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

        bfvals_ym = bfvals_yp
        bfvals_ym(1:nface,ky) = -bfvals_yp(1:nface,ky)
        bfvals_ym(1:nface,kyz) = -bfvals_yp(1:nface,kyz)
        bfvals_ym(1:nface,kxy) = -bfvals_yp(1:nface,kxy)
        bfvals_ym(1:nface,kxyz) = -bfvals_yp(1:nface,kxyz)
        bfvals_ym(1:nface,kyzz) = -bfvals_yp(1:nface,kyzz)
        bfvals_ym(1:nface,kxxy) = -bfvals_yp(1:nface,kxxy)
        ! bfvals_ym(1:nface,kyyy) = -bfvals_yp(1:nface,kyyy)
        bfvals_ym(1:nface,kyzxx) = -bfvals_yp(1:nface,kyzxx)
        bfvals_ym(1:nface,kxyzz) = -bfvals_yp(1:nface,kxyzz)
        bfvals_ym(1:nface,kyzzxx) = -bfvals_yp(1:nface,kyzzxx)

        bfvals_zp(1:nface,kyyz) = bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kz)     ! basis function = P_2(y)z
        bfvals_zp(1:nface,kyzz) = bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kzz)     ! basis function = P_2(z)y
        bfvals_zp(1:nface,kzzx) = bfvals_zp(1:nface,kzz)*bfvals_zp(1:nface,kx)     ! basis function = P_2(z)x
        bfvals_zp(1:nface,kzxx) = bfvals_zp(1:nface,kz)*bfvals_zp(1:nface,kxx)     ! basis function = P_2(x)z
        bfvals_zp(1:nface,kxxy) = bfvals_zp(1:nface,kxx)*bfvals_zp(1:nface,ky)     ! basis function = P_2(x)y
        bfvals_zp(1:nface,kxyy) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,kyy)     ! basis function = P_2(y)x
        bfvals_zp(1:nface,kyyzz) = bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kzz)     ! basis function = P_2(y)P_2(z)
        bfvals_zp(1:nface,kzzxx) = bfvals_zp(1:nface,kzz)*bfvals_zp(1:nface,kxx)     ! basis function = P_2(z)P_2(x)
        bfvals_zp(1:nface,kxxyy) = bfvals_zp(1:nface,kxx)*bfvals_zp(1:nface,kyy)     ! basis function = P_2(x)P_2(y)
        bfvals_zp(1:nface,kyzxx) = bfvals_zp(1:nface,kyz)*bfvals_zp(1:nface,kxx)     ! basis function = yz P_2(x)
        bfvals_zp(1:nface,kzxyy) = bfvals_zp(1:nface,kzx)*bfvals_zp(1:nface,kyy)     ! basis function = zx P_2(y)
        bfvals_zp(1:nface,kxyzz) = bfvals_zp(1:nface,kxy)*bfvals_zp(1:nface,kzz)     ! basis function = xy P_2(z)
        bfvals_zp(1:nface,kxyyzz) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kzz)   ! basis function = x P_2(y)P_2(z)
        bfvals_zp(1:nface,kyzzxx) = bfvals_zp(1:nface,ky)*bfvals_zp(1:nface,kzz)*bfvals_zp(1:nface,kxx)   ! basis function = y P_2(z)P_2(x)
        bfvals_zp(1:nface,kzxxyy) = bfvals_zp(1:nface,kz)*bfvals_zp(1:nface,kxx)*bfvals_zp(1:nface,kyy)   ! basis function = z P_2(x)P_2(y)
        bfvals_zp(1:nface,kxxyyzz) = bfvals_zp(1:nface,kx)*bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)


        bfvals_zm = bfvals_zp
        bfvals_zm(1:nface,kz) = -bfvals_zp(1:nface,kz)
        bfvals_zm(1:nface,kyz) = -bfvals_zp(1:nface,kyz)
        bfvals_zm(1:nface,kzx) = -bfvals_zp(1:nface,kzx)
        bfvals_zm(1:nface,kxyz) = -bfvals_zp(1:nface,kxyz)
        bfvals_zm(1:nface,kyyz) = -bfvals_zp(1:nface,kyyz)
        bfvals_zm(1:nface,kzxx) = -bfvals_zp(1:nface,kzxx)
        ! bfvals_zm(1:nface,kzzz) = -bfvals_zp(1:nface,kzzz)
        bfvals_zm(1:nface,kyzxx) = -bfvals_zp(1:nface,kyzxx)
        bfvals_zm(1:nface,kzxyy) = -bfvals_zp(1:nface,kzxyy)
        bfvals_zm(1:nface,kzxxyy) = -bfvals_zp(1:nface,kzxxyy)

        ! Organize local basis values on faces into 1-D vectors.
        ! Used in limiter() and max_lim().

        bf_faces(1:nslim,1) = 1.        ! basis function = 1

        do ixyz=kx,kz
            bf_faces(1:nface,ixyz) = bfvals_xm(1:nface,ixyz)        ! basis function = x,y,z
            bf_faces(nface+1:2*nface,ixyz) = bfvals_xp(1:nface,ixyz)        ! basis function = x,y,z
            bf_faces(2*nface+1:3*nface,ixyz) = bfvals_ym(1:nface,ixyz)        ! basis function = x,y,z
            bf_faces(3*nface+1:4*nface,ixyz) = bfvals_yp(1:nface,ixyz)        ! basis function = x,y,z
            bf_faces(4*nface+1:5*nface,ixyz) = bfvals_zm(1:nface,ixyz)        ! basis function = x,y,z
            bf_faces(5*nface+1:6*nface,ixyz) = bfvals_zp(1:nface,ixyz)        ! basis function = x,y,z
            bf_faces(6*nface+1:nslim,ixyz) = bfvals_int(1:npg,ixyz)        ! basis function = x,y,z
        end do

        bf_faces(1:nslim,kyz) = bf_faces(1:nslim,ky)*bf_faces(1:nslim,kz)     ! basis function = yz
        bf_faces(1:nslim,kzx) = bf_faces(1:nslim,kz)*bf_faces(1:nslim,kx)     ! basis function = zx
        bf_faces(1:nslim,kxy) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,ky)     ! basis function = xy
        bf_faces(1:nslim,kxyz) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,ky)*bf_faces(1:nslim,kz)     ! basis function = xyz

        do i=0,2
            bf_faces(1:nslim,kxx+i) = 1.5*bf_faces(1:nslim,kx+i)**2 - 0.5    ! basis function = P_2(s)
            ! bf_faces(1:nslim,kxxx+i) = 2.5*bf_faces(1:nslim,kx+i)**3 - 1.5*bf_faces(1:nslim,kx+i)   ! basis function = P_3(s)
        end do


        bf_faces(1:nslim,kyyz) = bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kz)     ! basis function = P_2(y)z
        bf_faces(1:nslim,kyzz) = bf_faces(1:nslim,ky)*bf_faces(1:nslim,kzz)     ! basis function = P_2(z)y
        bf_faces(1:nslim,kzzx) = bf_faces(1:nslim,kzz)*bf_faces(1:nslim,kx)     ! basis function = P_2(z)x
        bf_faces(1:nslim,kzxx) = bf_faces(1:nslim,kz)*bf_faces(1:nslim,kxx)     ! basis function = P_2(x)z
        bf_faces(1:nslim,kxxy) = bf_faces(1:nslim,kxx)*bf_faces(1:nslim,ky)     ! basis function = P_2(x)y
        bf_faces(1:nslim,kxyy) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,kyy)     ! basis function = P_2(y)x
        bf_faces(1:nslim,kyyzz) = bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kzz)     ! basis function = P_2(y)P_2(z)
        bf_faces(1:nslim,kzzxx) = bf_faces(1:nslim,kzz)*bf_faces(1:nslim,kxx)     ! basis function = P_2(z)P_2(x)
        bf_faces(1:nslim,kxxyy) = bf_faces(1:nslim,kxx)*bf_faces(1:nslim,kyy)     ! basis function = P_2(x)P_2(y)
        bf_faces(1:nslim,kyzxx) = bf_faces(1:nslim,kyz)*bf_faces(1:nslim,kxx)     ! basis function = yz P_2(x)
        bf_faces(1:nslim,kzxyy) = bf_faces(1:nslim,kzx)*bf_faces(1:nslim,kyy)     ! basis function = zx P_2(y)
        bf_faces(1:nslim,kxyzz) = bf_faces(1:nslim,kxy)*bf_faces(1:nslim,kzz)     ! basis function = xy P_2(z)
        bf_faces(1:nslim,kxyyzz) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kzz)   ! basis function = x P_2(y)P_2(z)
        bf_faces(1:nslim,kyzzxx) = bf_faces(1:nslim,ky)*bf_faces(1:nslim,kzz)*bf_faces(1:nslim,kxx)   ! basis function = y P_2(z)P_2(x)
        bf_faces(1:nslim,kzxxyy) = bf_faces(1:nslim,kz)*bf_faces(1:nslim,kxx)*bf_faces(1:nslim,kyy)   ! basis function = z P_2(x)P_2(y)
        bf_faces(1:nslim,kxxyyzz) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

    end subroutine set_face_vals_3D

!----------------------------------------------------------
    subroutine set_weights_3D

        integer i
        real wq4p,wq4m

        wq4p = (18. + sqrt(30.))/36.
        wq4m = (18. - sqrt(30.))/36.

        ! Define weights for integral approximation using Gaussian quadrature.

        ! Define weights for 1-D integration
        if (iquad .eq. 1) then        ! 1-point quadrature
            wgt1d(1) = 2.
        end if
        if (iquad .eq. 2) then        ! 2-point quadrature
            wgt1d(1:2) = 1.
        end if
        if (iquad .eq. 3) then        ! 3-point quadrature
            wgt1d(1) = 5./9.
            wgt1d(2) = 8./9.
            wgt1d(3) = 5./9.
        end if
        if (iquad .eq. 4) then        ! 4-point quadrature
            wgt1d(1) = wq4m
            wgt1d(2) = wq4p
            wgt1d(3) = wq4p
            wgt1d(4) = wq4m
        end if

        ! Define weights for 2-D integration
        if (iquad .eq. 1) then        ! 1-point quadrature
            wgt2d(1) = 4.
        end if
        if (iquad .eq. 2) then        ! 2-point quadrature
            wgt2d(1:4) = 1.
        end if

        if (iquad .ge. 3) then
            do i= 1,nedge
                wgt2d((i-1)*nedge+1:i*nedge) = wgt1d(1:nedge)*wgt1d(i)
            end do
        end if

        ! Define weights for 3-D integration
        if (iquad .eq. 1) then        ! 1-point quadrature
            wgt3d(1) = 8.
        end if
        if (iquad .eq. 2) then        ! 2-point quadrature
            wgt3d(1:8) = 1.
        end if

        if (iquad .ge. 3) then
            do i= 1,nedge
                wgt3d((i-1)*nface+1:i*nface) = wgt2d(1:nface)*wgt1d(i)
            end do
        end if

    end subroutine set_weights_3D


    !===========================================================================
    ! This routine tests whether the basis functions are orthonormal over the
    ! cell volume and over each cell face.
    !---------------------------------------------------------------------------
    subroutine test_basis_3D
        implicit none
        integer ir,jr,ix1,iy1,iz1,ipass
        real sq3,sq5,sq7,epsm

        sq3 = sqrt(3.)
        sq5 = sqrt(5.)
        sq7 = sqrt(7.)

        ! Indices of basis elements projected onto the x-faces (yz-plane):
        ibas_x(1) = 1
        ibas_x(kx) = 1
        ibas_x(ky) = ky
        ibas_x(kz) = kz
        ibas_x(kxy) = ky
        ibas_x(kyz) = kyz
        ibas_x(kzx) = kz
        ibas_x(kxx) = 1
        ibas_x(kzz) = kzz
        ibas_x(kyy) = kyy
        ibas_x(kxxy) = ky
        ibas_x(kyyz) = kyyz
        ibas_x(kzzx) = kzz
        ibas_x(kxyy) = kyy
        ibas_x(kyzz) = kyzz
        ibas_x(kzxx) = kz
        ibas_x(kxyz) = kyz
        ibas_x(kxxyy) = kyy
        ibas_x(kyyzz) = kyyzz
        ibas_x(kzzxx) = kzz
        ibas_x(kxyzz) = kyzz
        ibas_x(kyzxx) = kyz
        ibas_x(kzxyy) = kyyz
        ibas_x(kxyyzz) = kyyzz
        ibas_x(kyzzxx) = kyzz
        ibas_x(kzxxyy) = kyyz
        ibas_x(kxxyyzz) = kyyzz

        ! Indices of basis elements projected onto the y-faces (zx-plane):
        ibas_y(1) = 1
        ibas_y(kx) = kx
        ibas_y(ky) = 1
        ibas_y(kz) = kz
        ibas_y(kxy) = kx
        ibas_y(kyz) = kz
        ibas_y(kzx) = kzx
        ibas_y(kxx) = kxx
        ibas_y(kzz) = kzz
        ibas_y(kyy) = 1
        ibas_y(kxxy) = kxx
        ibas_y(kyyz) = kz
        ibas_y(kzzx) = kzzx
        ibas_y(kxyy) = kx
        ibas_y(kyzz) = kzz
        ibas_y(kzxx) = kzxx
        ibas_y(kxyz) = kzx
        ibas_y(kxxyy) = kxx
        ibas_y(kyyzz) = kzz
        ibas_y(kzzxx) = kzzxx
        ibas_y(kxyzz) = kzzx
        ibas_y(kyzxx) = kzxx
        ibas_y(kzxyy) = kzx
        ibas_y(kxyyzz) = kzzx
        ibas_y(kyzzxx) = kzzxx
        ibas_y(kzxxyy) = kzxx
        ibas_y(kxxyyzz) = kzzxx

        ! Indices of basis elements projected onto the z-faces (xy-plane):
        ibas_z(1) = 1
        ibas_z(kx) = kx
        ibas_z(ky) = ky
        ibas_z(kz) = 1
        ibas_z(kxy) = kxy
        ibas_z(kyz) = ky
        ibas_z(kzx) = kx
        ibas_z(kxx) = kxx
        ibas_z(kzz) = 1
        ibas_z(kyy) = kyy
        ibas_z(kxxy) = kxxy
        ibas_z(kyyz) = kyy
        ibas_z(kzzx) = kx
        ibas_z(kxyy) = kxyy
        ibas_z(kyzz) = ky
        ibas_z(kzxx) = kxx
        ibas_z(kxyz) = kxy
        ibas_z(kxxyy) = kxxyy
        ibas_z(kyyzz) = kyy
        ibas_z(kzzxx) = kxx
        ibas_z(kxyzz) = kxy
        ibas_z(kyzxx) = kxxy
        ibas_z(kzxyy) = kxyy
        ibas_z(kxyyzz) = kxyy
        ibas_z(kyzzxx) = kxxy
        ibas_z(kzxxyy) = kxxyy
        ibas_z(kxxyyzz) = kxxyy

        ! Normalization of basis elements over the  positive x-face of a cell:
        cbas_xp(1) = 1.
        cbas_xp(kx) = 1.
        cbas_xp(ky) = sq3
        cbas_xp(kz) = sq3
        cbas_xp(kxy) = sq3
        cbas_xp(kyz) = 3.
        cbas_xp(kzx) = sq3
        cbas_xp(kxx) = 1.
        cbas_xp(kzz) = sq5
        cbas_xp(kyy) = sq5
        cbas_xp(kxxy) = sq3
        cbas_xp(kyyz) = sq5*sq3
        cbas_xp(kzzx) = sq5
        cbas_xp(kxyy) = sq5
        cbas_xp(kyzz) = sq3*sq5
        cbas_xp(kzxx) = sq3
        cbas_xp(kxyz) = sq3*sq3
        cbas_xp(kxxyy) = sq5
        cbas_xp(kyyzz) = sq5*sq5
        cbas_xp(kzzxx) = sq5
        cbas_xp(kxyzz) = sq3*sq5
        cbas_xp(kyzxx) = sq3*sq3
        cbas_xp(kzxyy) = sq3*sq5
        cbas_xp(kxyyzz) = sq5*sq5
        cbas_xp(kyzzxx) = sq3*sq5
        cbas_xp(kzxxyy) = sq3*sq5
        cbas_xp(kxxyyzz) = sq5*sq5

        ! Normalization of basis elements over the  positive y-face of a cell:
        cbas_yp(1) = 1.
        cbas_yp(kx) = sq3
        cbas_yp(ky) = 1.
        cbas_yp(kz) = sq3
        cbas_yp(kxy) = sq3
        cbas_yp(kyz) = sq3
        cbas_yp(kzx) = 3.
        cbas_yp(kxx) = sq5
        cbas_yp(kzz) = sq5
        cbas_yp(kyy) = 1.
        cbas_yp(kxxy) = sq5
        cbas_yp(kyyz) = sq3
        cbas_yp(kzzx) = sq5*sq3
        cbas_yp(kxyy) = sq3
        cbas_yp(kyzz) = sq5
        cbas_yp(kzxx) = sq3*sq5
        cbas_yp(kxyz) = sq3*sq3
        cbas_yp(kxxyy) = sq5
        cbas_yp(kyyzz) = sq5
        cbas_yp(kzzxx) = sq5*sq5
        cbas_yp(kxyzz) = sq3*sq5
        cbas_yp(kyzxx) = sq3*sq5
        cbas_yp(kzxyy) = sq3*sq3
        cbas_yp(kxyyzz) = sq3*sq5
        cbas_yp(kyzzxx) = sq5*sq5
        cbas_yp(kzxxyy) = sq3*sq5
        cbas_yp(kxxyyzz) = sq5*sq5

        ! Normalization of basis elements over the  positive z-face of a cell:
        cbas_zp(1) = 1.
        cbas_zp(kx) = sq3
        cbas_zp(ky) = sq3
        cbas_zp(kz) = 1.
        cbas_zp(kxy) = 3.
        cbas_zp(kyz) = sq3
        cbas_zp(kzx) = sq3
        cbas_zp(kxx) = sq5
        cbas_zp(kzz) = 1.
        cbas_zp(kyy) = sq5
        cbas_zp(kxxy) = sq5*sq3
        cbas_zp(kyyz) = sq5
        cbas_zp(kzzx) = sq3
        cbas_zp(kxyy) = sq3*sq5
        cbas_zp(kyzz) = sq3
        cbas_zp(kzxx) = sq5
        cbas_zp(kxyz) = sq3*sq3
        cbas_zp(kxxyy) = sq5*sq5
        cbas_zp(kyyzz) = sq5
        cbas_zp(kzzxx) = sq5
        cbas_zp(kxyzz) = sq3*sq3
        cbas_zp(kyzxx) = sq3*sq5
        cbas_zp(kzxyy) = sq3*sq5
        cbas_zp(kxyyzz) = sq3*sq5
        cbas_zp(kyzzxx) = sq3*sq5
        cbas_zp(kzxxyy) = sq5*sq5
        cbas_zp(kxxyyzz) = sq5*sq5

        ! Normalization of basis elements over the  negative x-face of a cell:
        cbas_xm = cbas_xp
        cbas_xm(kx) = -1.
        cbas_xm(kxy) = -sq3
        cbas_xm(kzx) = -sq3
        cbas_xm(kxyy) = -sq5
        cbas_xm(kzzx) = -sq5
        cbas_xm(kxyz) = -sq3*sq3
        cbas_xm(kxyzz) = -sq3*sq5
        cbas_xm(kzxyy) = -sq3*sq5
        cbas_xm(kxyyzz) = -sq5*sq5

        ! Normalization of basis elements over the  negative y-face of a cell:
        cbas_ym = cbas_yp
        cbas_ym(ky) = -1.
        cbas_ym(kxy) = -sq3
        cbas_ym(kyz) = -sq3
        cbas_ym(kyzz) = -sq5
        cbas_ym(kxxy) = -sq5
        cbas_ym(kxyz) = -sq3*sq3
        cbas_ym(kxyzz) = -sq3*sq5
        cbas_ym(kyzxx) = -sq3*sq5
        cbas_ym(kyzzxx) = -sq5*sq5

        ! Normalization of basis elements over the  negative z-face of a cell:
        cbas_zm = cbas_zp
        cbas_zm(kz) = -1.
        cbas_zm(kyz) = -sq3
        cbas_zm(kzx) = -sq3
        cbas_zm(kzxx) = -sq5
        cbas_zm(kyyz) = -sq5
        cbas_zm(kxyz) = -sq3*sq3
        cbas_zm(kyzxx) = -sq3*sq5
        cbas_zm(kzxyy) = -sq3*sq5
        cbas_zm(kzxxyy) = -sq5*sq5


        ! Below, we specify what the scalar product of any two basis elements should
        ! evaluate to.

        ! cell_int0: scalar product over cell volume
        ! xp_int0: scalar product over positive x-face
        ! xm_int0: scalar product over negative x-face
        ! yp_int0: scalar product over positive y-face
        ! ym_int0: scalar product over negative y-face
        ! zp_int0: scalar product over positive z-face
        ! zm_int0: scalar product over negative z-face

        xp_int0 = 0.
        xm_int0 = 0.
        yp_int0 = 0.
        ym_int0 = 0.
        zp_int0 = 0.
        zm_int0 = 0.

        do jr=1,nbasis
            do ir=1,nbasis
                if (ir .ne. jr) cell_int0(ir,jr) = 0.
                if (ir .eq. jr) cell_int0(ir,jr) = 1.
                ix1 = 0
                iy1 = 0
                iz1 = 0

                if (ibas_x(ir) .eq. 1 .and. ibas_x(jr) .eq. 1) ix1 = 1
                if (ibas_x(ir) .eq. ky .and. ibas_x(jr) .eq. ky) ix1 = 1
                if (ibas_x(ir) .eq. kyy .and. ibas_x(jr) .eq. kyy) ix1 = 1
                if (ibas_x(ir) .eq. kz .and. ibas_x(jr) .eq. kz) ix1 = 1
                if (ibas_x(ir) .eq. kzz .and. ibas_x(jr) .eq. kzz) ix1 = 1
                if (ibas_x(ir) .eq. kyz .and. ibas_x(jr) .eq. kyz) ix1 = 1
                if (ibas_x(ir) .eq. kyyz .and. ibas_x(jr) .eq. kyyz) ix1 = 1
                if (ibas_x(ir) .eq. kyzz .and. ibas_x(jr) .eq. kyzz) ix1 = 1
                if (ibas_x(ir) .eq. kyyzz .and. ibas_x(jr) .eq. kyyzz) ix1 = 1

                if (ibas_y(ir) .eq. 1 .and. ibas_y(jr) .eq. 1) iy1 = 1
                if (ibas_y(ir) .eq. kx .and. ibas_y(jr) .eq. kx) iy1 = 1
                if (ibas_y(ir) .eq. kxx .and. ibas_y(jr) .eq. kxx) iy1 = 1
                if (ibas_y(ir) .eq. kz .and. ibas_y(jr) .eq. kz) iy1 = 1
                if (ibas_y(ir) .eq. kzz .and. ibas_y(jr) .eq. kzz) iy1 = 1
                if (ibas_y(ir) .eq. kzx .and. ibas_y(jr) .eq. kzx) iy1 = 1
                if (ibas_y(ir) .eq. kzxx .and. ibas_y(jr) .eq. kzxx) iy1 = 1
                if (ibas_y(ir) .eq. kzzx .and. ibas_y(jr) .eq. kzzx) iy1 = 1
                if (ibas_y(ir) .eq. kzzxx .and. ibas_y(jr) .eq. kzzxx) iy1 = 1

                if (ibas_z(ir) .eq. 1 .and. ibas_z(jr) .eq. 1) iz1 = 1
                if (ibas_z(ir) .eq. kx .and. ibas_z(jr) .eq. kx) iz1 = 1
                if (ibas_z(ir) .eq. kxx .and. ibas_z(jr) .eq. kxx) iz1 = 1
                if (ibas_z(ir) .eq. ky .and. ibas_z(jr) .eq. ky) iz1 = 1
                if (ibas_z(ir) .eq. kyy .and. ibas_z(jr) .eq. kyy) iz1 = 1
                if (ibas_z(ir) .eq. kxy .and. ibas_z(jr) .eq. kxy) iz1 = 1
                if (ibas_z(ir) .eq. kxxy .and. ibas_z(jr) .eq. kxxy) iz1 = 1
                if (ibas_z(ir) .eq. kxyy .and. ibas_z(jr) .eq. kxyy) iz1 = 1
                if (ibas_z(ir) .eq. kxxyy .and. ibas_z(jr) .eq. kxxyy) iz1 = 1

                if (ix1 .eq. 1) then
                    xp_int0(ir,jr) = 1.  !0.5*3.
                    xm_int0(ir,jr) = 1.  !0.5*3.
                end if
                if (iy1 .eq. 1) then
                    yp_int0(ir,jr) = 1.  !0.5*3.
                    ym_int0(ir,jr) = 1.  !0.5*3.
                end if
                if (iz1 .eq. 1) then
                    zp_int0(ir,jr) = 1.  !0.5*3.
                    zm_int0(ir,jr) = 1.  !0.5*3.
                end if
           end do
        end do


        ! Below, we evaluate the scalar product of all pairs of basis elements and
        ! test whether it's equal to the expected value computed above.

        ! cell_int: scalar product over cell volume
        ! xp_int: scalar product over positive x-face
        ! xm_int: scalar product over negative x-face
        ! yp_int: scalar product over positive y-face
        ! ym_int: scalar product over negative y-face
        ! zp_int: scalar product over positive z-face
        ! zm_int: scalar product over negative z-face

        ipass = 1

        do jr=1,nbasis
            do ir=1,nbasis
                cell_int(ir,jr) = 0.125*cbasis(ir)*sum(wgt3d(1:npg)*bfvals_int(1:npg,ir)*bfvals_int(1:npg,jr))
                xp_int(ir,jr) = 0.25*cbas_xp(ir)*cbas_xp(jr)*sum(wgt2d(1:nface)*bfvals_xp(1:nface,ir)*bfvals_xp(1:nface,jr))
                xm_int(ir,jr) = 0.25*cbas_xm(ir)*cbas_xm(jr)*sum(wgt2d(1:nface)*bfvals_xm(1:nface,ir)*bfvals_xm(1:nface,jr))
                yp_int(ir,jr) = 0.25*cbas_yp(ir)*cbas_yp(jr)*sum(wgt2d(1:nface)*bfvals_yp(1:nface,ir)*bfvals_yp(1:nface,jr))
                ym_int(ir,jr) = 0.25*cbas_ym(ir)*cbas_ym(jr)*sum(wgt2d(1:nface)*bfvals_ym(1:nface,ir)*bfvals_ym(1:nface,jr))
                zp_int(ir,jr) = 0.25*cbas_zp(ir)*cbas_zp(jr)*sum(wgt2d(1:nface)*bfvals_zp(1:nface,ir)*bfvals_zp(1:nface,jr))
                zm_int(ir,jr) = 0.25*cbas_zm(ir)*cbas_zm(jr)*sum(wgt2d(1:nface)*bfvals_zm(1:nface,ir)*bfvals_zm(1:nface,jr))

                epsm = 1.e-5

         !     if (iam .eq. print_mpi) then
            if (abs(cell_int(ir,jr) - cell_int0(ir,jr)) .gt. epsm) then
                print *,'incorrect cell integral: bases',ir,jr
                print *,'computed value: ',cell_int(ir,jr),'  expected value: ',cell_int0(ir,jr)
                ipass = 0
            end if
                if (abs(xp_int(ir,jr) - xp_int0(ir,jr)) .gt. epsm) then
                print *,'incorrect +x-face integral: bases',ir,jr
                print *,'computed value: ',xp_int(ir,jr),'  expected value: ',xp_int0(ir,jr)
                ipass = 0
            end if
                if (abs(xm_int(ir,jr) - xm_int0(ir,jr)) .gt. epsm) then
                print *,'incorrect -x-face integral: bases',ir,jr
                print *,'computed value: ',xm_int(ir,jr),'  expected value: ',xm_int0(ir,jr)
                ipass = 0
            end if
                if (abs(yp_int(ir,jr) - yp_int0(ir,jr)) .gt. epsm) then
                print *,'incorrect +y-face integral: bases',ir,jr
                print *,'computed value: ',yp_int(ir,jr),'  expected value: ',yp_int0(ir,jr)
                ipass = 0
            end if
                if (abs(ym_int(ir,jr) - ym_int0(ir,jr)) .gt. epsm) then
                print *,'incorrect -y-face integral: bases',ir,jr
                print *,'computed value: ',ym_int(ir,jr),'  expected value: ',ym_int0(ir,jr)
                ipass = 0
            end if
                if (abs(zp_int(ir,jr) - zp_int0(ir,jr)) .gt. epsm) then
                print *,'incorrect +z-face integral: bases',ir,jr
                print *,'computed value: ',zp_int(ir,jr),'  expected value: ',zp_int0(ir,jr)
                ipass = 0
            end if
                if (abs(zm_int(ir,jr) - zm_int0(ir,jr)) .gt. epsm) then
                print *,'incorrect -z-face integral: bases',ir,jr
                print *,'computed value: ',zm_int(ir,jr),'  expected value: ',zm_int0(ir,jr)
                ipass = 0
            end if
          !    end if
           end do
        end do

        ! if (ipass .eq. 0) call exit(-1)
    end subroutine test_basis_3D
    !---------------------------------------------------------------------------

end module basis_funcs
