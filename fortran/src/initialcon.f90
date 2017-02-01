module initialcon

use parameters
use helpers
use custom_boundary
use boundaries

contains

    subroutine set_ic(Q_r, id)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        ! character(*), intent(in) :: name
        integer i,j,k,id
        real wtev

        wtev = T_floor

        Q_r(:,:,:,:,:) = 0.0
        Q_r(:,:,:,rh,1) = rh_floor
        Q_r(:,:,:,en,1) = T_floor*rh_floor/(aindex - 1.)
        ! QMask(:,:,:) = .false.
        ! MMask(:,:,:) = .false.

        select case(id)
            case(0)
                call mpi_print(iam, 'Setting up 2D hydrodynamic jet (old version)...')
                call fill_fluid2(Q_r)
            case(1)
                call mpi_print(iam, 'Setting up 2D hydrodynamic jet...')
                call hydro_jet(Q_r)
            case(2)
                call mpi_print(iam, 'Setting up 2D isentropic vortex...')
                call isentropic_vortex(Q_r)
            case(3)
                call mpi_print(iam, 'Setting up 1D Sod Shock Tube v1...')
                call sod_shock_tube_1d(Q_r, 0)
            case(4)
                call mpi_print(iam, 'Setting up 1D Sod Shock Tube v2...')
                call sod_shock_tube_1d(Q_r, 1)
            case(5)
                call mpi_print(iam, 'Setting up 2D cylinder-in-pipe (laminar)...')
                call pipe_cylinder_2d(Q_r, 0)
            case(6)
                call mpi_print(iam, 'Setting up 2D cylinder-in-pipe (periodic)...')
                call pipe_cylinder_2d(Q_r, 1)
        end select

    end subroutine set_ic

!-------------------------------------------------------

    subroutine hydro_jet(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,ieq
        real x,y,wtev,rnum,jet_strength,rand_num,rh_fluid
        real qquad(npg,nQ),xcc,ycc,zcc  ! bfint(npg,nbasis),qquadv(npg)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        wtev = T_floor
        jet_strength = 1.0
        rh_fluid = 1.0

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

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
        real rh_amb,vx_amb,vy_amb,vz_amb,pr_amb,T_amb,beta
        real xctr,yctr,zctr,xp,yp,r2,delta_vx,delta_vy,delta_T
        real dn,vx,vy,vz,temp

        beta = 5.0             ! vortex strength
        rh_amb = 1.0           ! ambient density
        vx_amb = 1.0           ! ambient x-velocity
        vy_amb = 0.0           ! ambient y-velocity
        vz_amb = 0.0           ! ambient z-velocity
        T_amb  = 1.0           ! ambient temperature
        pr_amb = T_amb*rh_amb  ! ambient pressure

        xctr = 0               ! vortex center in x-direction
        yctr = 0               ! vortex center in y-direction
        zctr = 0               ! vortex center in z-direction

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx

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

    subroutine sod_shock_tube_1d(Q_r, version)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,version
        real rh_hi,rh_lo,pr_hi,pr_lo,xctr,yctr,xp,yp,dn,vx,vy,vz,pr

        rh_hi = 1.0e-4
        pr_hi = 1.0*P_base  ! atmospheric pressure
        rh_lo = 0.125e-4
        pr_lo = 0.1*P_base
        vx = 0.0
        vy = 0.0
        vz = 0.0

        xctr = 0
        yctr = 0  ! only needed if shock normal not aligned w/ x-direction

        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            xp = xc(i) - xctr
            if ( version .eq. 1 ) then
                yp = yc(j) - yctr
            else
                yp = 0
            endif

            if ( xp+yp .le. 0 ) then
                dn = rh_hi
                pr = pr_hi
            endif
            if ( xp+yp .gt. 0 ) then
                dn = rh_lo
                pr = pr_lo
            endif

            Q_r(i,j,k,rh,1) = dn
            Q_r(i,j,k,mx,1) = dn*vx
            Q_r(i,j,k,my,1) = dn*vy
            Q_r(i,j,k,mz,1) = dn*vz
            Q_r(i,j,k,en,1) = pr/aindm1 + 0.5*dn*(vx**2 + vy**2 + vz**2)
        end do
        end do
        end do

    end subroutine sod_shock_tube_1d


    !===========================================================================
    ! 2D pipe flow around a cylinder
    !   * Re ~ 20 for laminar case (version 0)
    !   * Re ~ 100 for periodic case (version 1)
    !   * Need outflow BCs on right wall: nu d u/d eta - p*eta = 0
    !   * No-slip BCs everywhere else
    ! NOTE: This subroutine sets BCs for the problem using custom_boundary
    !------------------------------------------------------------
    subroutine pipe_cylinder_2d(Q_r, version)
        use custom_boundary

        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        logical, dimension(nx,ny,nz):: Qmask

        integer i,j,k,i4,version
        real dn,pr,vx,vy,vz,ux_amb,cyl_x0,cyl_y0,cyl_rad

        !-------------------------------------------------------
        ! Definitions
        dn = 1.0
        pr = 1.0*P_base  ! can be adjusted to get stable results
        vx = 0.0
        vy = 0.0
        vz = 0.0

        cyl_x0 = 0.2
        cyl_y0 = 0.2
        cyl_rad = 0.05

        select case(version)
            case(0)
                ux_amb = 1.5
            case(1)
                ux_amb = 0.3
        end select

        !-------------------------------------------------------
        ! Set initial conditions
        do k = 1,nz
        do j = 1,ny
        do i = 1,nx
            Q_r(i,j,k,rh,1) = dn
            Q_r(i,j,k,mx,1) = dn*vx
            Q_r(i,j,k,my,1) = dn*vy
            Q_r(i,j,k,mz,1) = dn*vz
            Q_r(i,j,k,en,1) = pr/aindm1 + 0.5*dn*(vx**2 + vy**2 + vz**2)
        end do
        end do
        end do

        !-------------------------------------------------------
        ! Generate cylinder mask --> zero out grid quantities in mask
        !   NOTE: might be an offset b/c mask is constructed at faces, not cells
        Qmask(:,:,:) = cyl_in_2d_pipe_mask(cyl_x0, cyl_y0, cyl_rad)

        do i4=1,nface
            where (Qmask)
                Q_r(:,:,:,i4,1) = 0.0
            end where
        end do

        ! NOTE: Qxlow_ext_custom, Qcyl_ext, and QMask should already be initialized!
        ! call add_custom_boundaries(icname)
        call set_cyl_in_2d_pipe_boundaries(Qxlow_ext_custom, Qcyl_ext, Qmask, ux_amb)

    end subroutine pipe_cylinder_2d
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    subroutine fill_fluid2(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,igrid,ieq
        real x,y,wtev,rnum,rand_num,rh_fluid
        real qquad(npg,nQ),xcc,ycc,zcc  ! bfint(npg,nbasis),qquadv(npg)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        rh_fluid = 1.0
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


end module initialcon
