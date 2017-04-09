!***** INITIALIZE.F90 ********************************************************************
module initialize

use input
use params
use helpers
use basis_funcs!, only: wgt1d, wgt2d, wgt3d, ibitri, cbasis, set_bfvals_3D

use initialcon
use random  ! TODO: commented to get working w/o MKL

implicit none

integer, parameter :: iseed = 123456789  ! 1317345*mpi_P + 5438432*mpi_Q + 38472613*mpi_R

real :: cflm
integer :: nout

contains

    !===========================================================================
    ! initializer : perform general setup & variable initialization
    !------------------------------------------------------------
    subroutine initializer(t, dt, nout, comm)
        implicit none
        real, intent(out) :: t, dt
        integer, intent(out) :: nout

        integer :: comm

        !-------------------------
        ! Create output directory
        call system('mkdir -p '//outdir)  ! NOTE: transfer to scratch in BW batch script!!!

        !-----------------------------------------
        ! Initialize grid sizes and local lengths
        ibitri = set_ibitri(iquad, nbasis)
        cbasis = set_cbasis_3D(ibitri)
        cflm = set_cflm(iquad, ibitri)

        !-----------------------------------------
        ! Initialize various parameters
        call setup_MPI(comm, mpi_nx, mpi_ny, mpi_nz, mpi_P, mpi_Q, mpi_R)
        call init_spatial_params(dx,dy,dz, dxi,dyi,dzi, loc_lxd,loc_lyd,loc_lzd)
        ! call set_mxa_mya_mza(mxa, mya, mza)  ! not needed anymore probably
        call init_temporal_params(t, dt, nout)
        !-----------------------------------------

        !-----------------------------------------
        ! Evaluate local cell values of basis functions on cell interior and faces
        call set_bfvals_3D  ! This is done for 1, 2, or 3 point Gaussian quadrature

        call init_bf_weights(bval_int_wgt, wgtbf_xmp, wgtbf_ymp, wgtbf_zmp)

        ! call init_random_seed(iam, iseed)

        call random_init(iseed)  ! initialize MKL random number generator

        call print_startup_info

    end subroutine initializer
    !---------------------------------------------------------------------------


    !===========================================================================
    ! integer function setup_output(iquad, nbasis)
    !     implicit none
    !     integer, intent(in) :: iquad, nbasis
    !
    !     nstdout
    !     nstldout
    !     nstvout
    !     nsttout
    !     nsteiout
    !     nstenout
    !     nstesout
    !     nstpout
    !     nststout
    !     nstvrout
    !
    !     setup_output
    !     return
    ! end function setup_output
    !---------------------------------------------------------------------------


    !===========================================================================
    integer function set_ibitri(iquad, nbasis)
        implicit none
        integer, intent(in) :: iquad, nbasis

        if (iquad == 2 .and. nbasis == 4 ) set_ibitri = 0
        if (iquad == 3 .and. nbasis == 10) set_ibitri = 0
        if (iquad == 4 .and. nbasis == 20) set_ibitri = 0
        if (iquad == 2 .and. nbasis == 8 ) set_ibitri = 1
        if (iquad == 3 .and. nbasis == 27) set_ibitri = 1
        if (iquad == 4 .and. nbasis == 64) set_ibitri = 1
        return
    end function set_ibitri
    !---------------------------------------------------------------------------


    !===========================================================================
    real function set_cflm(iquad, ibitri)
        implicit none
        integer, intent(in) :: iquad, ibitri

        ! if (nbasis <= 8)  cflm = 0.14  ! nbasis = 4   0.28 is unstable for hydro
        ! if (nbasis == 27) cflm = 0.1   ! nbasis = 10  0.15 is unstable for hydro
        ! if (nbasis == 64) cflm = 0.08  ! nbasis = 20  0.1  is unstable for hydro
        select case (ibitri)
            case (0)
        	    if (iquad == 2) set_cflm = 0.14
        	    if (iquad == 3) set_cflm = 0.08
        	    if (iquad == 4) set_cflm = 0.05
        	case (1)  ! coefficients for basis functions {P2(x)P2(y)P2(z)}
        	    if (iquad == 2) set_cflm = 0.14
        	    if (iquad == 3) set_cflm = 0.08
        	    if (iquad == 4) set_cflm = 0.05
        end select
        return
    end function set_cflm
    !---------------------------------------------------------------------------

    !===========================================================================
    subroutine set_basis_function_flags(ibitri)
        implicit none
        integer, intent(in) :: ibitri
        integer klim
        ! Using bi/tri elements:
            ! iquad=2, nbasis=8:  { 1,x,y,z, yz,zx,xy, xyz }
            ! iquad=3, nbasis=27: { 1,x,y,z, yz,zx,xy, xyz, P2(x),P2(y),P2(z),
            !     yP2(z), zP2(x), xP2(y), P2(y)z, P2(z)x, P2(x)y,
            !     P2(y)P2(z), P2(z)P2(x), P2(x)P2(y), yzP2(x),zxP2(y),xyP2(z),
            !     xP2(y)P2(z),yP2(z)P2(x),zP2(x)P2(y), P2(x)P2(y)P2(z) }
        ! Not using bi/tri elements:
            ! iquad=2, nbasis=4:  { 1,x,y,z }
            ! iquad=3, nbasis=10: { 1,x,y,z, yz,zx,xy, P2(x),P2(y),P2(z) }
            ! iquad=4, nbasis=20: { 1,x,y,z, yz,zx,xy, P2(x),P2(y),P2(z), xyz,
            !     yP2(z), zP2(x), xP2(y), P2(y)z, P2(z)x, P2(x)y,
            !     P3(x), P3(y), P3(z) }
        kx  = 2
        ky  = 3
        kz  = 4
        kyz = 5
        kzx = 6
        kxy = 7

        if ( ibitri == 0 ) klim = kxy
        if ( ibitri == 1 ) then
            kxyz = 8
            klim = kxyz
        end if
        kxx = klim + 1
        kyy = klim + 2
        kzz = klim + 3
        if ( ibitri == 0 ) then
            kxyz = klim + 4
            klim = kxyz
        end if

        if ( ibitri == 1 ) klim = kzz
            kyzz = klim + 1
            kzxx = klim + 2
            kxyy = klim + 3
            kyyz = klim + 4
            kzzx = klim + 5
            kxxy = klim + 6
        if ( ibitri == 1 ) klim = kxxy
        if ( ibitri == 0 ) then
            kxxx = klim + 7
            kyyy = klim + 8
            kzzz = klim + 9
            klim = kzzz
        end if

        ! If ibitri = 1.......
        kyyzz   = klim + 1
        kzzxx   = klim + 2
        kxxyy   = klim + 3
        kyzxx   = klim + 4
        kzxyy   = klim + 5
        kxyzz   = klim + 6
        kxyyzz  = klim + 7
        kyzzxx  = klim + 8
        kzxxyy  = klim + 9
        kxxyyzz = klim + 10
        if ( ibitri == 1 ) then
            kxxx = klim + 11
            kyyy = klim + 12
            kzzz = klim + 13
            klim = kzzz
        end if
    end subroutine set_basis_function_flags
    !---------------------------------------------------------------------------


    !===========================================================================
    function set_cbasis_3D(ibitri)
        implicit none
        integer, intent(in) :: ibitri
        real, dimension(nbastot) :: set_cbasis_3D
        real, dimension(nbastot) :: cbasis

        call set_basis_function_flags(ibitri)
        cbasis(1)             = 1.   ! basis func coeff   {1}
        cbasis(kx:kz)         = 3.   ! basis funcs coeffs {x,y,z}
        cbasis(kyz:kxy)       = 9.   ! basis funcs coeffs {yz,zx,xy}
        cbasis(kxyz)          = 27.  ! basis func coeff   {xyz}
        cbasis(kxx:kzz)       = 5.   ! basis funcs coeffs {P2(x),   P2(y),  P2(z)}
        cbasis(kyzz:kxyy)     = 15.  ! basis funcs coeffs {yP2(z), zP2(x), xP2(y)}
        cbasis(kyyz:kxxy)     = 15.  ! basis funcs coeffs {P2(y)z, P2(z)y, P2(z)x}
        cbasis(kyyzz:kxxyy)   = 25.  ! basis funcs coeffs {P2(y)P2(z), P2(z)P2(x), P2(x)P2(y)}
        cbasis(kyzxx:kxyzz)   = 45.  ! basis funcs coeffs {yzP_2(x),   zxP_2(y),   xyP_2(z)}
        cbasis(kxyyzz:kzxxyy) = 75.  ! basis funcs coeffs {xP2(y)P2(z),yP2(z)P2(x),zP2(x)P2(y)}
        cbasis(kxxyyzz)       = 125. ! basis funcs coeffs {P2(x)P2(y)P2(z)}
        cbasis(kxxx:kzzz)     = 7.   ! basis funcs coeffs {P3(x), P3(y), P3(z)}
        set_cbasis_3D(:) = cbasis(:)
        return
    end function set_cbasis_3D
    !---------------------------------------------------------------------------


    !===========================================================================
    subroutine setup_MPI(comm, mpi_nx, mpi_ny, mpi_nz, mpi_P, mpi_Q, mpi_R)
        implicit none
        integer, intent(inout) :: mpi_nx, mpi_ny
        integer, intent(out)   :: mpi_nz
        integer, intent(out)   :: mpi_P, mpi_Q, mpi_R

        integer :: comm

        integer, dimension(3) :: dims, coords, periods ! only used in init for MPI things
        integer reorder

        ! NOTE: USES the following global parameters
        !   * MPI_COMM_WORLD, numprocs, ierr, cartcomm, nbrs
        !   * EAST, WEST, NORTH, SOUTH, UP, DOWN
        !   * mpi_nx, mpi_ny (set by user)
        !   * clt, tf, ntout, xhibc, yhibc, zhibc (set by user)
        ! NOTE: SETS the following global parameters
        !   * mpi_nz
        !   * mpi_P, mpi_Q, mpi_R

        ! call MPI_Init ( ierr )
        ! call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
        call MPI_COMM_SIZE(comm, numprocs, ierr)

        mpi_nz = numprocs/(mpi_nx*mpi_ny)

        dims(1) = mpi_nx
        dims(2) = mpi_ny
        dims(3) = mpi_nz

        periods(:) = 0
        if (xhibc == 'periodic') then
            periods(1) = 1
        end if
        if (yhibc == 'periodic') then
            periods(2) = 1
        end if
        if (zhibc == 'periodic') then
            periods(3) = 1
        end if
        reorder = 1

        call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, reorder,cartcomm, ierr)
        call MPI_COMM_RANK (cartcomm, iam, ierr )
        call MPI_CART_COORDS(cartcomm, iam, 3, coords, ierr)

        call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(WEST), nbrs(EAST), ierr)
        call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(SOUTH), nbrs(NORTH), ierr)
        call MPI_CART_SHIFT(cartcomm, 2, 1, nbrs(DOWN), nbrs(UP), ierr)

        mpi_P = coords(1) + 1
        mpi_Q = coords(2) + 1
        mpi_R = coords(3) + 1
    end subroutine setup_MPI
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Currently only needed for HLLC in flux.f90
    !------------------------------------------------------------
    subroutine set_mxa_mya_mza(mxa, mya, mza)
        implicit none
        integer, dimension(3), intent(out) :: mxa, mya, mza
        mxa = (/ mx, my, mz /)
        mya = (/ my, mz, mx /)
        mza = (/ mz, mx, my /)
    end subroutine set_mxa_mya_mza
    !---------------------------------------------------------------------------


    !===========================================================================
    !   loc_l[]d: Set the starting x,y,z coords for domain of this MPI process
    !   NOTE: the center of the computational grid is the origin (0,0,0)
    !------------------------------------------------------------
    subroutine init_spatial_params(dx,dy,dz, dxi,dyi,dzi, loc_lxd,loc_lyd,loc_lzd)
        implicit none
        real, intent(out) :: dx,dy,dz, dxi,dyi,dzi
        real, intent(out) :: loc_lxd, loc_lyd, loc_lzd

        ! NOTE: USES the following global parameters
        !   * lx, ly, lz (set by user)
        !   * mpi_nx, mpi_ny, mpi_nz, mpi_P, mpi_Q, mpi_R (depends on setup_MPI)
        ! NOTE: SETS the following global parameters
        !   * dx,dy,dz, dxi,dyi,dzi
        !   * loc_lxd, loc_lyd, loc_lzd

        lxd = -(lx/2.0)
        lxu =  (lx/2.0)
        lyd = -(ly/2.0)
        lyu =  (ly/2.0)
        lzd = -(lz/2.0)
        lzu =  (lz/2.0)

        dxi = (nx*mpi_nx)/(lxu-lxd)
        dyi = (ny*mpi_ny)/(lyu-lyd)
        dzi = (nz*mpi_nz)/(lzu-lzd)
        dx = 1./dxi
        dy = 1./dyi
        dz = 1./dzi

        loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
        loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
        loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz
    end subroutine init_spatial_params
    !---------------------------------------------------------------------------


    !===========================================================================
    subroutine init_temporal_params(t, dt, nout)
        implicit none
        real, intent(out) :: t, dt
        integer, intent(out) :: nout
        ! NOTE: USES the following global parameters
        !   * dx (set by init_spatial_params)
        !   * cflm (set by set_cflm)
        !   * clt, ntout (set by user)
        ! NOTE: SETS the following global parameters
        !   * t, dt, nout
        t = 0.0
        dt = cflm*dx/clt
        nout = 0
    end subroutine init_temporal_params
    !---------------------------------------------------------------------------


    !===========================================================================
    subroutine init_bf_weights(bval_int_wgt, wgtbf_xmp, wgtbf_ymp, wgtbf_zmp)
        implicit none
        real, dimension(npg,nbastot), intent(out) :: bval_int_wgt
        real, dimension(nface,2,nbastot), intent(out) :: wgtbf_xmp, wgtbf_ymp, wgtbf_zmp

        real, dimension(nface,nbastot) :: wgtbfvals_xp, wgtbfvals_xm
        real, dimension(nface,nbastot) :: wgtbfvals_yp, wgtbfvals_ym
        real, dimension(nface,nbastot) :: wgtbfvals_zp, wgtbfvals_zm
        integer ir,ipg

        ! NOTE: USES the following global parameters
        !   * dxi, dyi, dzi (set by init_spatial_params)
        !   * wgt2d, wgt3d (set by set_weights_3D)
        !   * nbasis, npg, nface (set implicitly by user)
        !   * cbasis (set by set_cbasis_3D)

        do ir=1,nbasis
            do ipg=1,npg
                bval_int_wgt(ipg,ir) = wgt3d(ipg)*bfvals_int(ipg,ir)
            end do
        end do

        do ir=1,nbasis
            wgtbfvals_xp(1:nface,ir) = wgt2d(1:nface)*bfvals_xp(1:nface,ir)
            wgtbfvals_yp(1:nface,ir) = wgt2d(1:nface)*bfvals_yp(1:nface,ir)
            wgtbfvals_zp(1:nface,ir) = wgt2d(1:nface)*bfvals_zp(1:nface,ir)
            wgtbfvals_xm(1:nface,ir) = wgt2d(1:nface)*bfvals_xm(1:nface,ir)
            wgtbfvals_ym(1:nface,ir) = wgt2d(1:nface)*bfvals_ym(1:nface,ir)
            wgtbfvals_zm(1:nface,ir) = wgt2d(1:nface)*bfvals_zm(1:nface,ir)
        end do

        do ir=1,nbasis
            wgtbf_xmp(1:nface,1,ir) = -0.25*cbasis(ir)*dxi*wgtbfvals_xm(1:nface,ir)
            wgtbf_ymp(1:nface,1,ir) = -0.25*cbasis(ir)*dyi*wgtbfvals_ym(1:nface,ir)
            wgtbf_zmp(1:nface,1,ir) = -0.25*cbasis(ir)*dzi*wgtbfvals_zm(1:nface,ir)
            wgtbf_xmp(1:nface,2,ir) =  0.25*cbasis(ir)*dxi*wgtbfvals_xp(1:nface,ir)
            wgtbf_ymp(1:nface,2,ir) =  0.25*cbasis(ir)*dyi*wgtbfvals_yp(1:nface,ir)
            wgtbf_zmp(1:nface,2,ir) =  0.25*cbasis(ir)*dzi*wgtbfvals_zp(1:nface,ir)
        end do
    end subroutine init_bf_weights
    !---------------------------------------------------------------------------


    !===========================================================================
    ! set_ic_from_file : initialize simulation from checkpoint file
    !------------------------------------------------------------
    subroutine set_ic_from_file(Q_r, t, dt, dtout, nout)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
        real, intent(inout) :: t, dt, dtout
        integer, intent(inout) :: nout

        real t_p,dt_p,dtout_p
        integer nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p
        ! This applies only if the initial data are being read from an input file.
        ! - If resuming a run, keep the previous clock (i.e., t at nout) running.
        ! - If not resuming a run, treat input as initial conditions at t=0, nout=0.
        call readQ(fpre,iam,iread,Q_r,t_p,dt_p,nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p)

        if (resuming) then
            t = t_p
            dt = dt_p
            nout = nout_p
        end if
        ! Note, nout=1 corresponds to t=dt, but nout=2 corresponds to t~dtout, etc.
        if (nout > 1) then
            dtout_p = t_p/(nout_p-1)
        else  ! Automatically pass consistency check
            dtout_p = dtout
        end if
        if (iam == print_mpi) then
            print *, 'resuming = ', resuming
            print *, 't = ', t
            print *, 'dt = ', dt
            print *, 'nout = ', nout
            print *, 'dtout_p = ', dtout_p, ' dtout = ', dtout
            print *, 'mpi_nx_p = ', mpi_nx_p, ' mpi_nx = ', mpi_nx
            print *, 'mpi_ny_p = ', mpi_ny_p, ' mpi_ny = ', mpi_ny
            print *, 'mpi_nz_p = ', mpi_nz_p, ' mpi_nz = ', mpi_nz
        end if
            ! Quit if dtout is incompatible with input t/(nout-1)
        if (abs(dtout_p-dtout)/dt_p > 1.01) then
            call mpi_print(iam, 'Bad restart, non-matching dtout')
            call exit(-1)
        end if
        if ((mpi_nx_p /= mpi_nx) .or. (mpi_ny_p /= mpi_ny) .or. (mpi_nz_p /= mpi_nz)) then
            call mpi_print(iam, 'Bad restart, non-matching mpi_nx, mpi_ny, or mpi_nz')
            call exit(-1)
        end if
    end subroutine set_ic_from_file
    !---------------------------------------------------------------------------


    !===========================================================================
    ! print_startup_info : print initial information about simulation
    !------------------------------------------------------------
    subroutine print_startup_info()

        if (iam == print_mpi) then
            print *, ''
            print *, '---------------------------------------------------------'
            print *, 'Starting simulation...'
            print *, '---------------------------------------------------------'
            write(*,'(A13,I10,I7,I7)') ' total dim = ', mpi_nx*nx, mpi_ny*ny, mpi_nz*nz
            write(*,'(A13,I10,I7,I7)') ' mpi dim   = ', mpi_nx,    mpi_ny,    mpi_nz
            write(*,'(A13,ES10.3)')    ' te0 is    = ', te0
            write(*,'(A13,ES10.3)')    ' dx is     = ', ly/(ny*mpi_ny)*L0
            write(*,'(A13,I10)')       ' iquad is  = ', iquad
            write(*,'(A13,I10)')       ' nbasis is = ', nbasis
            print *, '----------------------------------------------'
            write(*,'(A16,A8,A13,A8)') ' X BC:  lower = ', xlobc, '  |  upper = ', xhibc
            write(*,'(A16,A8,A13,A8)') ' Y BC:  lower = ', ylobc, '  |  upper = ', yhibc
            write(*,'(A16,A8,A13,A8)') ' Z BC:  lower = ', zlobc, '  |  upper = ', zhibc
            print *, '---------------------------------------------------------'
            print *, ''
        end if

    end subroutine print_startup_info
    !---------------------------------------------------------------------------


    !===========================================================================
    ! writeQ : Write a checkpoint file
    !------------------------------------------------------------
    subroutine writeQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,              &
                      mpi_nxnow,mpi_nynow,mpi_nznow)

        implicit none
        real :: Qin(nx,ny,nz,nQ,nbasis),tnow,dtnow
        integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow,nump,numd,qq,k,j,i,ir
        character (4) :: fprefix,pname,dname
        character (5) :: pname1,dname1
        character (30) :: fname2

        nump = iam + 10000

        write(pname1,'(i5)')nump
        pname=pname1(2:5)
        pname = trim(pname)
        pname = adjustr(pname)

        numd = iddump + 10000

        write(dname1,'(i5)')numd
        dname=dname1(2:5)
        dname = trim(dname)
        dname = adjustr(dname)

        fname2 = 'data/'//fprefix//'_p'//pname//'_d'//dname//'.dat'
        ! print *,'fname2 ',fname2

        open(unit=3,file=fname2)

        ! open(unit = 10, file = 'data/perseus_t'//dname//'_p'//pname//'.bin',form = 'unformatted',access = 'stream')

        do ir=1,nbasis
            do qq=1,nQ
                do k=1,nz
                    do j=1,ny
                        write(3,*) (Qin(i,j,k,qq,ir),i=1,nx)
                    enddo
                enddo
            enddo
        enddo

        write(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow
        close(3)

    end subroutine writeQ
    !---------------------------------------------------------------------------


    !===========================================================================
    ! Read a checkpoint file (set iread to nonzero integer)
    !------------------------------------------------------------
    subroutine readQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,               &
                     mpi_nxnow,mpi_nynow,mpi_nznow)
        implicit none
        real :: Qin(nx,ny,nz,nQ,nbasis),tnow,dtnow
        integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow,nump,numd,qq,k,j,i,ir
        character (4) :: fprefix,pname,dname
        character (5) :: pname1,dname1
        character (30) :: fname2

        nump = irank + 10000

        write(pname1,'(i5)')nump
        pname=pname1(2:5)
        pname = trim(pname)
        pname = adjustr(pname)

        numd = iddump + 10000

        write(dname1,'(i5)')numd
        dname=dname1(2:5)
        dname = trim(dname)
        dname = adjustr(dname)

        fname2 = 'data/'//fpre//'_p'//pname//'_d'//dname//'.dat'

        open(unit=3,file=fname2,action='read')

        do ir=1,nbasis
        do qq=1,nQ
            do k=1,nz
                do j=1,ny
                    read(3,*) (Qin(i,j,k,qq,ir),i=1,nx)
                enddo
            enddo
        enddo
        enddo

        read(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow
        close(3)

    end subroutine readQ
    !---------------------------------------------------------------------------

end module initialize
