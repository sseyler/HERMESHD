module basis_funcs

use parameters

contains


    subroutine set_bfvals_3D
        ! Defines local basis function values and weights for 1, 2, or 3-point Gaussian quadrature.
        ! Basis functions are evaluated in cell interior and on cell faces.

        implicit none

        call set_vtk_vals_3D()    ! Define basis function values on a 3D grid INTERNAL to the cell.
        call set_internal_vals_3D()     ! Define basis function values at quadrature points INTERNAL to cell.
        call set_face_vals_3D()   ! Define local basis function values at quadrature points on a cell face.
        call set_weights_3D()     ! Define weights for integral approximation using Gaussian quadrature.

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

end module basis_funcs
