module initialize

use parameters
use helpers

use ic_mod
use basis_funcs_mod
use io_mod

implicit none

integer :: nout


!===============================================================================
contains

    subroutine perform_setup()

        implicit none
        integer ir,ipg
        real half_length

        if (nbasis .le. 8) cflm = 0.14
        if (nbasis .eq. 27) cflm = 0.1
        if (nbasis .eq. 64) cflm = 0.08

        ! Initialize grid sizes and local lengths

        cbasis(1) = 1.              ! coeff for basis func {1}
        cbasis(kx:kz) = 3.          ! coeff for basis funcs {x,y,z}
        cbasis(kyz:kxy) = 9.        ! coeff for basis funcs {yz,zx,xy}
        cbasis(kxyz) = 27.          ! coeff for basis func {xyz}

        cbasis(kxx:kzz) = 5.        ! coeff for basis funcs {P2(x),P2(y),P2(z)}
        cbasis(kyzz:kxyy) = 15.     ! coeff for basis funcs {yP2(z),zP2(x),xP2(y)}
        cbasis(kyyz:kxxy) = 15.     ! coeff for basis funcs {P2(y)z,P2(z)y,P2(z)x}
        cbasis(kyyzz:kxxyy) = 25.   ! coeff for basis funcs {P2(y)P2(z),P2(z)P2(x),P2(x)P2(y)}
        cbasis(kyzxx:kxyzz) = 45.   ! coeff for basis funcs {yzP_2(x),zxP_2(y),xyP_2(z)}
        cbasis(kxyyzz:kzxxyy) = 75. ! coeff for basis funcs {xP2(y)P2(z),yP2(z)P2(x),zP2(x)P2(y)}
        cbasis(kxxyyzz) = 125.      ! coeff for basis funcs {P2(x)P2(y)P2(z)}

        call MPI_Init ( ierr )
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

        mpi_nz=numprocs/(mpi_nx*mpi_ny)

        dims(1)=mpi_nx
        dims(2)=mpi_ny
        dims(3)=mpi_nz

        periods(:)=0
        if (xhbc .eq. 2) then
            periods(1)=1
        end if
        if (yhbc .eq. 2) then
            periods(2)=1
        end if
        if (zhbc .eq. 2) then
            periods(3)=1
        end if
        reorder = 1

        call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, reorder,cartcomm, ierr)
        call MPI_COMM_RANK (cartcomm, iam, ierr )
        call MPI_CART_COORDS(cartcomm, iam, 3, coords, ierr)
        mpi_P=coords(1)+1
        mpi_Q=coords(2)+1
        mpi_R=coords(3)+1
        call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(WEST), nbrs(EAST), ierr)
        call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(SOUTH), nbrs(NORTH), ierr)
        call MPI_CART_SHIFT(cartcomm, 2, 1, nbrs(DOWN), nbrs(UP), ierr)

        half_length = lx/2.0
        lxd = -half_length
        lxu = half_length
        half_length = ly/2.0
        lyd = -half_length
        lyu = half_length
        half_length = lz/2.0
        lzd = -half_length
        lzu = half_length

        dxi = (nx*mpi_nx)/(lxu-lxd)
        dyi = (ny*mpi_ny)/(lyu-lyd)
        dzi = (nz*mpi_nz)/(lzu-lzd)

        dx = 1./dxi
        dy = 1./dyi
        dz = 1./dzi

        ! Set the starting x,y,z coords for the domain of this MPI process
        !   Note: the center of the computational grid is the origin (0,0,0)
        loc_lxd = lxd + (mpi_P-1)*(lxu-lxd)/mpi_nx
        loc_lyd = lyd + (mpi_Q-1)*(lyu-lyd)/mpi_ny
        loc_lzd = lzd + (mpi_R-1)*(lzu-lzd)/mpi_nz

        rh_fluid = 1.

        ! indices used in poynting() for computing Poynting Maxwell flux

        mxa(1) = mx
        mxa(2) = my
        mxa(3) = mz
        mya(1) = my
        mya(2) = mz
        mya(3) = mx
        mza(1) = mz
        mza(2) = mx
        mza(3) = my

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
            wgtbf_xmp(1:nface,2,ir) = 0.25*cbasis(ir)*dxi*wgtbfvals_xp(1:nface,ir)
            wgtbf_ymp(1:nface,2,ir) = 0.25*cbasis(ir)*dyi*wgtbfvals_yp(1:nface,ir)
            wgtbf_zmp(1:nface,2,ir) = 0.25*cbasis(ir)*dzi*wgtbfvals_zp(1:nface,ir)
        end do


        call init_random_seed(123456789)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 38472613*mpi_R

        ! Initialize MKL random number generator
        vsl_errcode = vslnewstream(vsl_stream, vsl_brng, iseed)


        call print_startup_info()

        t_start = get_clock_time()

        if (iread .eq. 0) then
            call set_ic(Q_r0, icid)
        else
            call initialize_from_file()
        endif

    end subroutine perform_setup


    subroutine initialize_from_file()

        implicit none

        real t_p,dt_p,dtout_p
        integer nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p
        ! This applies only if the initial data are being read from an input file.
        ! - If resuming a run, keep the previous clock (i.e., t at nout) running.
        ! - If not resuming a run, treat input as initial conditions at t=0, nout=0.
        call readQ(fpre,iam,iread,Q_r0,t_p,dt_p,nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p)

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


    subroutine print_startup_info()
        if (iam .eq. print_mpi) then
            print *, ''
            print *, '----------------------------------------------'
            print *, 'Starting simulation...'
            print *, '----------------------------------------------'
            print *, 'total dim= ',mpi_nx*nx,mpi_ny*ny,mpi_nz*nz
            print *, 'mpi dim= ',mpi_nx,mpi_ny,mpi_nz
            print *, 'te0 is: ', te0
            print *, 'dx is: ', ly/(ny*mpi_ny)*L0
            print *, 'iquad is: ',  iquad
            print *, 'nbasis is: ', nbasis
            print *, '----------------------------------------------'
            print *, ''
        end if
    end subroutine print_startup_info

end module initialize
