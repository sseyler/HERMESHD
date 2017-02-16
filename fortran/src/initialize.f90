module initialize

use parameters
use helpers

use flux, only : mxa,mya,mza, set_mxa_mya_mza

use initialcon
use basis_funcs

implicit none

real :: cflm
integer :: nout

contains

    !===========================================================================
    ! perform_setup : perform general setup & variable initialization
    !------------------------------------------------------------
    subroutine perform_setup
        implicit none
        integer :: ir,ipg

        !-------------------------
        ! Create output directory
        call system('mkdir -p '//outdir)

        !-----------------------------------------
        ! Initialize grid sizes and local lengths
        ibitri = set_ibitri(iquad, nbasis)
        cbasis = set_cbasis_3D(ibitri)

        cflm = set_cflm(iquad, ibitri)

        call setup_MPI

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
        dVi = dxi*dyi*dzi

        ! Set the starting x,y,z coords for the domain of this MPI process
        !   Note: the center of the computational grid is the origin (0,0,0)
        loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
        loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
        loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz

        call set_mxa_mya_mza(mxa, mya, mza)

        t = 0.
        dt = cflm*dx/clt
        dtoriginal = dt
        nout = 0
        dtout = tf/ntout

        ! Evaluate local cell values of basis functions on cell interior and faces.
        ! This is done for 1, 2, or 3 point Gaussian quadrature.
        call set_bfvals_3D

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


        call init_random_seed(iam, 123456789)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 38472613*mpi_R

        ! Initialize MKL random number generator
        vsl_errcode = vslnewstream(vsl_stream, vsl_brng, iseed)

        call print_startup_info

    end subroutine perform_setup
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
    integer function set_cflm(iquad, ibitri)
        implicit none
        integer, intent(in) :: iquad, ibitri

        ! if (nbasis <= 8)  cflm = 0.14  ! nbasis = 4   0.28 is unstable for hydro
        ! if (nbasis == 27) cflm = 0.1   ! nbasis = 10  0.15 is unstable for hydro
        ! if (nbasis == 64) cflm = 0.08  ! nbasis = 20  0.1  is unstable for hydro
        select case (ibitri)
            case (0)
        	    if (iquad == 2) set_cflm = 0.2
        	    if (iquad == 3) set_cflm = 0.12
        	    if (iquad == 4) set_cflm = 0.08
        	case (1)  ! coefficients for basis functions {P2(x)P2(y)P2(z)}
        	    if (iquad == 2) set_cflm = 0.14
        	    if (iquad == 3) set_cflm = 0.10
        	    if (iquad == 4) set_cflm = 0.07
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
        real, intent(out) :: cbasis(nbastot)

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
        return
    end function set_cbasis_3D
    !---------------------------------------------------------------------------


    !===========================================================================
    subroutine setup_MPI
        integer, dimension(3) :: dims, coords, periods ! only used in init for MPI things

        call MPI_Init ( ierr )
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

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
        mpi_P = coords(1) + 1
        mpi_Q = coords(2) + 1
        mpi_R = coords(3) + 1
        call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(WEST), nbrs(EAST), ierr)
        call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(SOUTH), nbrs(NORTH), ierr)
        call MPI_CART_SHIFT(cartcomm, 2, 1, nbrs(DOWN), nbrs(UP), ierr)
    end subroutine setup_MPI
    !---------------------------------------------------------------------------


    !===========================================================================
    ! initialize_from_file : initialize simulation from checkpoint file
    !------------------------------------------------------------
    subroutine initialize_from_file(Q_r)
        implicit none
        real, dimension(nx,ny,nz,nQ,nbasis), intent(inout) :: Q_r
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
        if (nout .gt. 1) then
            dtout_p = t_p/(nout_p-1)
        else  ! Automatically pass consistency check
            dtout_p = dtout
        end if
        if (iam .eq. print_mpi) then
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
            if (iam .eq. print_mpi) then
                print *, 'Bad restart, non-matching dtout'
            end if
            call exit(-1)
        end if
        if ((mpi_nx_p .ne. mpi_nx) .or. (mpi_ny_p .ne. mpi_ny) .or. (mpi_nz_p .ne. mpi_nz)) then
            if (iam .eq. print_mpi) then
                print *, 'Bad restart, non-matching mpi_nx, mpi_ny, or mpi_nz'
            end if
            call exit(-1)
        end if
    end subroutine initialize_from_file
    !---------------------------------------------------------------------------


    !===========================================================================
    ! print_startup_info : print initial information about simulation
    !------------------------------------------------------------
    subroutine print_startup_info()

        if (iam .eq. print_mpi) then
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
