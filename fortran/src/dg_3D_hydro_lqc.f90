!      DG-Hydro in 3 dimensions: Extension of the finite volume PERSEUS Algorithm (Physics of the Extended-mhd Relaxation System using an Efficient Upwind Scheme, 
!                                by Matt Martin and Charles Seyler) to the Discontinuous Galerkin Method specialized to hydrodynamics.
!                
!                DG Algorithm by Xuan Zhao, Nat Hamlin and Charles Seyler. 
!
!                Solves the 3D compressible Euler equations
!                 
!                Variables:    There are 5 dependent variables: density (rh), velocity (vx,vy,vz), energy (en)
!
!                Units:        A value of unity for each variable or parameter corresponds to the following dimensional units:
!
!                              Length              L0
!                              Time                t0
!                              Number density      n0 
!                              Velocity            v0
!                              Temperature         te0
!
!-------------------------------------------------------------------------------------------------------------------------------------------  

!********************  SUMMARY OF TEST RESULTS USING VARIUOS NUMBER OF BASIS FUNCTIONS AND QUADRATURE POINTS  ******************************
!
!   The nbasis = 8 and iquad = 2 case gives extremely close results to nbasis = 10 and iquad = 3.  The (8,2) case takes 21 seconds per iteration and (10,3) takes 35 seconds.      
!   Close examination of the output indicates that the (10,3) case has slightly better resolution/sharpness/definition as seen in the transitions from high to low 
!   or low to high vorticity.  The (27,3) case has sightly better resolution/sharpness/definition than the (10,3) case in about the same proportion as the (10,3) case does compared
!   to the (8,2) case.  However, the time for an iteration is about 190 sec, which is about 4x higher than the (10,3) case.  It would appear for this test there is no significant
!   advantage of using 27 elements over the 10 elements. There of ourse could be applications in which the 27 elements does provide significant improvement of accuracy.    
!   Finally in comparing the cases (3,27) to (4,20) there does not appear to be any real differences, only that the (4,20) case was about 15% faster.

    program main  
    use lib_vtk_io 
    implicit none	

    include 'mpif.h'  
	integer, parameter :: rh=1, mx=2, my=3, mz=4, en=5, exx=6, eyy=7, ezz=8, exy=9, exz=10, eyz=11, nQ=11

	integer, parameter :: nx = 40, ny = 40, nz = 1, ngu = 0, nbasis = 27, nbastot = 30

	        ! nbasis = 4: {1, x, y, z}   (use iquad = 2)
		! nbasis = 8: {1, x, y, z, yz, zx, xy, xyz}    (use iquad = 2)
	        ! nbasis = 10: {1, x, y, z, yz, zx, xy, P_2(x), P_2(y), P_2(z)}   (use iquad = 3)
	        ! nbasis = 20: nbasis10 + {xyz, xP2(y), yP2(x), xP2(z), zP2(x), yP2(z), zP2(y), P3(x), P3(y), P3(z)}  (use iquad = 4)
		! nbasis = 27: {1, x, y, z, yz, zx, xy, xyz, P2(x), P2(y), P2(z), yP2(z), zP2(x), xP2(y), P2(y)z, P2(z)x, P2(x)y,
		!               P2(y)P2(z), P2(z)P2(x), P2(x)P2(y), yzP2(x), zxP2(y), xyP2(z), xP2(y)P2(z), yP2(z)P2(x), zP2(x)P2(y), P2(x)P2(y)P2(z)}  (use iquad = 3)

        ! The jump in accuracy between the linear basis (nbasis = 4) and quadratic basis (nbasis = 10)
        ! is much greater than the jump in accuracy between quadratic and cubic (nbasis = 20) bases.

	integer, parameter :: iquad=3, nedge = iquad
	! iquad: number of Gaussian quadrature points per direction. iquad should not be less than ipoly, where ipoly is the maximum Legendre polynomial order used,
        ! otherwise instability results. iquad should not be larger than ipoly + 1, which is an exact Gaussian quadrature for the Legendre polynomial used. Thus there are only
        ! two cases of iquad for a given nbasis.  Both cases give similar results although the iquad = ipoly + 1 case is formally more accurate. 
        ! Note: In 3D iquad = 2 corresponds to 8 internal points and 4 face points, iquad = 3 corresponds to 27 internal points and 9 face points,
        ! iquad = 4 corresponds to 64 internal points and 16 face points.

    	integer, parameter :: nface = iquad*iquad, npg = nface*iquad, nfe = 2*nface, npge = 6*nface, nslim=npg+6*nface
	! nface: number of quadrature points per cell face.
	! npg: number of internal points per cell.

	integer, parameter :: xlbc = 2, xhbc = 2, ylbc = 2, yhbc = 2, zlbc = 2, zhbc = 2

	integer, parameter :: ntout = 100, iorder = 2

	integer, parameter :: ihllc = 1, iroe = 0, ivis = -1
	! Choose Riemann solver for computing fluxes.  Set chosen solver = 1.
	! If all of them = 0, then LLF is used for fluxes.
        ! LLF is very diffusive for the hydro problem.  Roe and HLLC are much less diffusive than LLF and give very similar results with similar cpu overhead.

	! to restart from a checkpoint, set iread to 1 or 2 (when using the odd/even scheme)
	integer, parameter :: iread = 0, iwrite = 0
	character (4), parameter :: fpre = 'Qout'
	logical, parameter :: resuming = .false.

	real, parameter :: lx = 300., ly = 300., lz = 1.25 

	real, parameter :: tf = 2000.
                 
        real, parameter :: pi = 4.0*atan(1.0)

        real, parameter :: aindex = 5./3., mu = 22.

	real, parameter :: aindm1 = aindex - 1.0, cp = aindex/(aindex - 1.0), clt = 2.

    	! dimensional units (expressed in MKS)
        real, parameter :: L0=1.0e4, t0=1.0e2, n0=6.e18

    	! derived units 
        real, parameter :: v0=L0/t0, p0=mu*1.67e-27*n0*v0**2, te0=p0/n0/1.6e-19, vis=1.e-3

	real, parameter :: T_floor = 0.01/te0, rh_floor = 1.e-5, P_floor = T_floor*rh_floor, rh_mult = 1.01

 
	real, dimension(nx,ny,nz,nQ,nbasis) ::  Q_r0, Q_r1, Q_r2, Q_r3	
	real, dimension(nx,ny,nz,nQ,nbasis) ::  glflux_r, source_r
	real, dimension(nx,ny,nz,nQ,nbasis) :: integral_r 
	
 	real eta(nx,ny,nz,npg),den0(nx,ny,nz),Ez0,Zdy(nx,ny,nz,npg),coll
	real flux_x(nface,1:nx+1,ny,nz,1:nQ), flux_y(nface,nx,1:ny+1,nz,1:nQ), flux_z(nface,nx,ny,1:nz+1,1:nQ)
        real cfrx(nface,nQ),cfry(nface,nQ),cfrz(nface,nQ)
	logical MMask(nx,ny,nz),BMask(nx,ny,nz)   
	real xcell(npg), ycell(npg), zcell(npg), xface(npge), yface(npge), zface(npge)
	integer ticks, count_rate, count_max
	real t1, t2, t3, t4, elapsed_time, t_start, t_stop, dtoriginal
        real t,dt,tout,dtout,vf,dxi,dyi,dzi,loc_lxd,loc_lyd,loc_lzd,check_Iz,sl, dz, dy, dx
        real lxd,lxu,lyd,lyu,lzd,lzu,pin_rad,pin_height,foil_thickness,rh_foil,rh_fluid,pin_rad_in,pin_rad_out,rim_rad
        real disk_rad,disk_height,foil_rad,buf_rad,buf_z,dish_height,foil_height
        real gpz_rad,rh_gpz,kappa
	real Qxhigh_ext(ny,nz,nface,nQ), Qxlow_int(ny,nz,nface,nQ), Qxlow_ext(ny,nz,nface,nQ), Qxhigh_int(ny,nz,nface,nQ) 
	real Qyhigh_ext(nx,nz,nface,nQ), Qylow_int(nx,nz,nface,nQ), Qylow_ext(nx,nz,nface,nQ), Qyhigh_int(nx,nz,nface,nQ)	
	real Qzhigh_ext(nx,ny,nface,nQ), Qzlow_int(nx,ny,nface,nQ), Qzlow_ext(nx,ny,nface,nQ), Qzhigh_int(nx,ny,nface,nQ)
	integer mxa(3),mya(3),mza(3),kroe(nface),niter

	integer ir,ipg,nout
	real half_length

!   Variables for resuming from saved data.
	real t_p,dt_p,dtout_p
	integer nout_p,mpi_nx_p,mpi_ny_p,mpi_nz_p,ioe


	! Parameters relating to quadratures and basis functions.

	real wgt1d(5), wgt2d(30), wgt3d(100), cbasis(nbastot)
	! wgt1d: quadrature weights for 1-D integration
	! wgt2d: quadrature weights for 2-D integration
	! wgt3d: quadrature weights for 3-D integration	

	
	real, dimension(nface,nbastot) :: bfvals_zp, bfvals_zm, bfvals_yp, bfvals_ym, bfvals_xp, bfvals_xm
	real bf_faces(nslim,nbastot), bfvals_int(npg,nbastot),xquad(20)
        real wgtbf_xmp(nface,2,nbastot),wgtbf_ymp(nface,2,nbastot),wgtbf_zmp(nface,2,nbastot)
        real bval_int_wgt(npg,nbastot),wgtbfvals_xp(nface,nbastot),wgtbfvals_yp(nface,nbastot),wgtbfvals_zp(nface,nbastot)
        real wgtbfvals_xm(nface,nbastot),wgtbfvals_ym(nface,nbastot),wgtbfvals_zm(nface,nbastot)

        integer kx,ky,kz,kyz,kzx,kxy,kxyz,kxx,kyy,kzz,kyzz,kzxx,kxyy
        integer kyyz,kzzx,kxxy,kyyzz,kzzxx,kxxyy,kyzxx,kzxyy,kxyzz,kxyyzz,kyzzxx,kzxxyy,kxxyyzz
        integer kxxx,kyyy,kzzz

	real, dimension(nbastot,nbastot) :: cell_int,xp_int,xm_int,yp_int,ym_int,zp_int,zm_int
	real, dimension(nbastot,nbastot) :: cell_int0,xp_int0,xm_int0,yp_int0,ym_int0,zp_int0,zm_int0
	real, dimension(nbastot) :: cbas_xp,cbas_xm,cbas_yp,cbas_ym,cbas_zp,cbas_zm	
	integer ibas_x(nbastot), ibas_y(nbastot), ibas_z(nbastot), ibitri

	! Parameters for VTK output.

	integer, parameter :: nvtk=1
	integer, parameter :: nvtk2=nvtk*nvtk, nvtk3=nvtk*nvtk*nvtk
	integer(I4P), parameter :: nnx=nx*nvtk, nny=ny*nvtk, nnz=nz*nvtk
	real, dimension(nvtk3,nbastot) :: bfvtk, bfvtk_dx, bfvtk_dy, bfvtk_dz
	real xgrid(20),dxvtk,dyvtk,dzvtk
	
			
    ! MPI definitions

	integer :: mpi_nx=6, mpi_ny=6, print_mpi=0  
       
    	integer iam,ierr,mpi_nz,numprocs,reorder,cartcomm,mpi_P,mpi_Q,mpi_R
    	integer dims(3),coords(3),periods(3),nbrs(6),reqs(4),stats(MPI_STATUS_SIZE,4)
    	integer,parameter:: NORTH=1,SOUTH=2,EAST=3,WEST=4,UP=5,DOWN=6,MPI_TT=MPI_REAL4

        real cflm
   
    ! Initialize grid sizes and local lengths

	call set_cbasis_3D

!	if(nbasis .eq. 4)  cflm = 0.25     !  nbasis = 4   0.28 is unstable for hydro
!	if(nbasis .eq. 10) cflm = 0.12     !  nbasis = 10  0.15 is unstable for hydro
!	if(nbasis .eq. 20) cflm = 0.08     !  nbasis = 20  0.1 is unstable for hydro

	if (ibitri .eq. 0) then
	   if (iquad .eq. 2) cflm = 0.2
	   if (iquad .eq. 3) cflm = 0.12 
	   if (iquad .eq. 4) cflm = 0.08
	end if
	if (ibitri .eq. 1) then
	   if (iquad .eq. 2) cflm = 0.14
	   if (iquad .eq. 3) cflm = 0.10 
	   if (iquad .eq. 4) cflm = 0.07
	end if	! coefficients for basis functions {P2(x)P2(y)P2(z)}
    
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
	niter = 0
	dtout = tf/ntout

	! Evaluate local cell values of basis functions on cell interior and faces.
	! This is done for 1, 2, 3, or 4 point Gaussian quadrature.
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
    
    if (iam .eq. print_mpi) then
    	print *,'total dim= ',mpi_nx*nx,mpi_ny*ny,mpi_nz*nz
    	print *,'mpi dim= ',mpi_nx,mpi_ny,mpi_nz
	print *, 'te0 is: ', te0
	print *, 'dx is: ', ly/(ny*mpi_ny)*L0
        print *, 'iquad is: ',  iquad
        print *, 'nbasis is: ', nbasis
    end if
    	
	call system_clock(ticks, count_rate, count_max)
	t_start = ticks*1./count_rate
	
	if (iread .eq. 0) then

		call initial_condition

	else

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

	end if

	
	call system_clock( ticks, count_rate, count_max )
    	t1 = 1.*ticks / count_rate
	call output_vtk(Q_r0,nout,iam)	
    	
	do while (t<tf)

!	   if (mod(niter,200) .eq. 0 .and. iam .eq. print_mpi) print *,'niter,t,dt = ',niter,t,dt,dtout*nout
	   niter = niter + 1
 	call get_min_dt(dt)
		
	if (iorder .eq. 2) then
		
		call prep_advance(Q_r0)
		call advance_time_level_gl(Q_r0,Q_r1)

		call prep_advance(Q_r1)
		call advance_time_level_gl(Q_r1,Q_r2)

		Q_r0 = 0.5*(Q_r0 + Q_r2)
	end if

	if (iorder .eq. 3) then

		call prep_advance(Q_r0)
		call advance_time_level_gl(Q_r0,Q_r1)

		call prep_advance(Q_r1)
		call advance_time_level_gl(Q_r1,Q_r2)

		Q_r3 = 0.75*Q_r0 + 0.25*Q_r2

		call prep_advance(Q_r3)		
		call advance_time_level_gl(Q_r3,Q_r2)	

                Q_r0 = (Q_r0 + 2.*Q_r2)/3.
			
	end if

		t = t+dt	

	if (t .gt. dtout*nout) then
		nout = nout+1
		if (iam .eq. print_mpi) then
			call system_clock(ticks, count_rate, count_max)
			t2 = 1.*ticks/count_rate
			print *, 'Iteration time', (t2-t1), 'seconds'
			t1 = t2
			print *, 'nout = ', nout
		        print *, "        t= ",t*100.,"         dt= ",dt
		end if	
		
		call MPI_BARRIER(cartcomm,ierr)
		call output_vtk(Q_r0,nout,iam)	

		! write checkpoint files; assign an odd/even id to ensure last two sets are kept
		if (iwrite .eq. 1) then
			ioe = 2 - mod(nout,2) 
			call writeQ(fpre,iam,ioe,Q_r0,t,dt,nout,mpi_nx,mpi_ny,mpi_nz)
		end if
		call MPI_BARRIER(cartcomm,ierr)

		if (iam .eq. print_mpi) then
			call system_clock(ticks, count_rate, count_max)
			t2 = ticks/count_rate
			print *, 'Output time', (t2-t1), 'seconds'
			t1 = t2
		end if	

	end if

	end do 

    	call MPI_Finalize (ierr)

    contains

!------------------------------------------        
	subroutine init_random_seed(iseed)
		implicit none
		integer :: i, n, clock,iseed
		integer, dimension(:), allocatable :: seed
		
		call random_seed(size = n)
		allocate(seed(n))
		  
		call system_clock(count=clock)
		  
		if (iseed .eq. 0) then 
			seed =  clock*(iam+1) + 37 * (/ (i - 1, i = 1, n) /)
		else
			seed =  iseed*(iam+1)
		endif
		call random_seed(put = seed)
		  
	!	print *,seed(1)
		deallocate(seed)
	end subroutine

!------------------------------------------------------------------------------- 

       subroutine prep_advance(Q_ri)

	    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri

		call prepare_exchange(Q_ri)
		call set_bc
		call flux_cal(Q_ri)
		call innerintegral(Q_ri)
		call glflux
!		call source_calc(Q_ri,t)

       end subroutine prep_advance

!------------------------------------------------------------------------------- 
 
       subroutine get_min_dt(dt)
       real dt,dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0,vex,vey,vez,vem,vem0,dni,dn,vx,vy,vz,Pr,sqdni,vacc,vacc0,cs
       integer :: i,j,k,main_proc=0,mpi_size=1,loc_reqs(numprocs-1),loc_stats(MPI_STATUS_SIZE,numprocs-1)

       vmag = 0.
		
       do k=1,nz
       do j=1,ny
       do i=1,nx
              
             dn = Q_r0(i,j,k,rh,1)
             dni = 1./dn
             vx = Q_r0(i,j,k,mx,1)*dni*dxi
             vy = Q_r0(i,j,k,my,1)*dni*dyi
             vz = Q_r0(i,j,k,mz,1)*dni*dzi
             cs = sqrt(aindex*(Q_r0(i,j,k,en,1)*dni - 0.5*(vx**2 + vy**2 + vz**2)))

             vmag0 = max(abs(vx),abs(vy),abs(vz))
             if(vmag0 > vmag .and. dn > rh_mult*rh_floor)vmag = vmag0

       end do
       end do  
       end do

       vmax = (vmag + cs)*dxi
       dt_min = cflm/vmax  !max(5.,vmax)
       
		call MPI_BARRIER(cartcomm,ierr)
		if (iam.eq.main_proc) then
		do i=1,numprocs-1
        	call MPI_IRecv(dt_val(i),mpi_size,MPI_TT,i,0,cartcomm,loc_reqs(i),ierr)        			
		enddo
		call MPI_WaitAll(numprocs-1,loc_reqs,loc_stats,ierr)
		do i=1,numprocs-1
        	dt_min=min(dt_min,dt_val(i))
		enddo
		do i=1,numprocs-1
        	call MPI_ISend(dt_min,mpi_size,MPI_TT,i,0,cartcomm,loc_reqs(i),ierr)        			
		enddo
		call MPI_WaitAll(numprocs-1,loc_reqs,loc_stats,ierr)
		else
        	call MPI_ISend(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
		call MPI_Wait(reqs(1),stats(:,1),ierr)
        	call MPI_IRecv(dt_min,mpi_size,MPI_TT,main_proc,0,cartcomm,reqs(1),ierr)
		call MPI_Wait(reqs(1),stats(:,1),ierr)
		endif
		
		dt = dt_min
		
	end subroutine get_min_dt

!----------------------------------------------------------------------------------------------

    subroutine initial_condition
    implicit none
    integer i,j,k
    real coeff
    real den, Pre, wtev

        wtev = T_floor
	
	Q_r0(:,:,:,:,:) = 0.0
	Q_r0(:,:,:,rh,1) = rh_floor
	Q_r0(:,:,:,en,1) = T_floor*rh_floor/(aindex - 1.)	
	MMask(:,:,:) = .false.

        call fill_fluid

    end subroutine initial_condition

!-------------------------------------------------------

    subroutine fill_fluid3
        
!-------------------------------------------------------

        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,igrid,ieq,inds(nbasis)
        real wirer,x,y,x0,y0,rhline,wtev,rx,ry,rnum,w0,P,vx,vy,vz,rho
	real qquad(npg,nQ),xcc,ycc,zcc,rand_number
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        wtev = T_floor
	w0 = 0.3
        coll = rh_fluid*wtev/vis

!       test problem is an unstable flow jet in x with velocity perturbations in y

        do i = 1,nx
           do j = 1,ny
              do k = 1,nz

              call random_number(rand_number)
              rnum = (rand_number - 0.5)

	         Q_r0(i,j,k,rh,1) = rh_floor
	         Q_r0(i,j,k,en,1) = wtev*Q_r0(i,j,k,rh,1)/(aindex - 1.)

		 xcc = xc(i)	
		 ycc = yc(j) 	
		 zcc = zc(k) 
               
                 Q_r0(i,j,k,rh,1) = rh_fluid
		 Q_r0(i,j,k,my,1) = 0.0001*rnum/cosh(20*yc(j)/lyu)**2 
		 Q_r0(i,j,k,mz,1) = 0.  
		 Q_r0(i,j,k,mx,1) = 1.0*Q_r0(i,j,k,rh,1)/cosh(20*ycc/lyu)/1.
	         Q_r0(i,j,k,en,1) = (1. + 0.001*rnum)*wtev*Q_r0(i,j,k,rh,1)/(aindex - 1.)  &
                                  + 0.5*(Q_r0(i,j,k,mx,1)**2 + Q_r0(i,j,k,my,1)**2 + Q_r0(i,j,k,mz,1)**2)/Q_r0(i,j,k,rh,1)

                 rho = rh_fluid
                 vx = Q_r0(i,j,k,mx,1)/rh_fluid
                 vy = Q_r0(i,j,k,my,1)/rh_fluid
                 vz = Q_r0(i,j,k,mz,1)/rh_fluid
                 P = wtev*Q_r0(i,j,k,rh,1)

                 Q_r0(i,j,k,exx,1) = P + rho*vx**2 
                 Q_r0(i,j,k,eyy,1) = P + rho*vy**2
                 Q_r0(i,j,k,ezz,1) = P + rho*vz**2
                 Q_r0(i,j,k,exy,1) = rho*vx*vy 
                 Q_r0(i,j,k,exz,1) = rho*vx*vz
                 Q_r0(i,j,k,eyz,1) = rho*vy*vz

              end do     
           end do
        end do
           
    end subroutine fill_fluid3

!-------------------------------------------------------

    subroutine fill_fluid
        
!-------------------------------------------------------

        integer i,j,k,ir,izw(4),ixw(4),iyw(4),iw,iseed,igrid,ieq,inds(nbasis)
        real wirer,x,y,x0,y0,rhline,wtev,rx,ry,rnum,w0,rand_number
	real qquad(npg,nQ),xcc,ycc,zcc,rho,vx,vy,vz,P  ! bfint(npg,nbasis),qquadv(npg)
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R


        wtev = 2*T_floor
	w0 = 0.3
        coll = rh_fluid*wtev/vis


        do i = 1,nx
           do j = 1,ny
              do k = 1,nz

              call random_number(rand_number)
              rnum = (rand_number - 0.5)

	      qquad(:,:) = 0.

	      do igrid=1,npg

	         qquad(igrid,rh) = rh_floor
	         qquad(igrid,en) = wtev*qquad(igrid,rh)/(aindex - 1.)

		 xcc = xc(i) + bfvals_int(igrid,kx)*0.5/dxi	
		 ycc = yc(j) + bfvals_int(igrid,ky)*0.5/dyi	
		 zcc = zc(k) + bfvals_int(igrid,kz)*0.5/dzi
               
                 qquad(igrid,rh) = rh_fluid
		 qquad(igrid,my) = 0.0001*rnum  
		 qquad(igrid,mz) = 0.  
		 qquad(igrid,mx) = 1.0*qquad(igrid,rh)/cosh(20*ycc/lyu)
	         qquad(igrid,en) = (1. + 0.001*rnum)*wtev*qquad(igrid,rh)/(aindex - 1.) + 0.5*(qquad(igrid,mx)**2 + qquad(igrid,my)**2 + qquad(igrid,mz)**2)/qquad(igrid,rh)

                 rho = rh_fluid
                 vx = qquad(igrid,mx)/rh_fluid
                 vy = qquad(igrid,my)/rh_fluid
                 vz = qquad(igrid,mz)/rh_fluid
                 P = wtev*qquad(igrid,rh)

	         if (ivis .ge. 0) then
                    qquad(igrid,exx) = P + rho*vx**2 
                    qquad(igrid,eyy) = P + rho*vy**2
                    qquad(igrid,ezz) = P + rho*vz**2
                    qquad(igrid,exy) = rho*vx*vy 
                    qquad(igrid,exz) = rho*vx*vz
                    qquad(igrid,eyz) = rho*vy*vz
		 end if

	      end do

	      do ieq=rh,en
		 do ir=1,nbasis  
		    Q_r0(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(wgt3d(1:npg)*bfvals_int(1:npg,ir)*qquad(1:npg,ieq))
		 end do
	      end do

              end do     
           end do
        end do
           
    end subroutine fill_fluid

!-------------------------------------------------------

    subroutine fill_fluid2
        
!-------------------------------------------------------

        integer i,j,k,izw(4),ixw(4),iyw(4),iw,iseed
        real wirer,x,y,x0,y0,rhline,wtev,rx,ry,rnum,w0,rand_number
        iseed = 1317345*mpi_P + 5438432*mpi_Q + 3338451*mpi_R

        wtev = 2*T_floor
	w0 = 0.3

        do i = 1,nx
           do j = 1,ny
              do k = 1,nz
              call random_number(rand_number)

                rnum = (rand_number - 0.5)
                
                Q_r0(i,j,k,rh,1) = rh_fluid!*exp(-zc(k)/Hscale)
		Q_r0(i,j,k,my,1) = 0.005*rnum/cosh(20*yc(j)/lyu)**2
                rnum = (rand_number - 0.5)
!		Q_r0(i,j,k,mz,1) = 0.005*rnum 	
!		Q_0(i,j,k,mx) = 1.0*Q_0(i,j,k,rh)/cosh(10*(yc(j)-0.25*lyu)/lyu) - 1.0*Q_0(i,j,k,rh)/cosh(10*(yc(j)-0.75*lyu)/lyu) &
!                              - 1.0*Q_0(i,j,k,rh)/cosh(10*(yc(j)+0.25*lyu)/lyu) + 1.0*Q_0(i,j,k,rh)/cosh(10*(yc(j)+0.75*lyu)/lyu)
		Q_r0(i,j,k,mx,1) = 1.0*Q_r0(i,j,k,rh,1)/cosh(20*yc(j)/lyu)
!               Q_0(i,j,k,my) = rh_fluid*sin(4*pi*xc(i)/lx)*cos(4*pi*yc(j)/ly)
!               Q_0(i,j,k,mx) = -rh_fluid*sin(4*pi*yc(j)/ly)  !*cos(4*pi*xc(i)/lx)  
!	        Q_0(i,j,k,en) = grav*Hscale*Q_0(i,j,k,rh)/(aindex - 1.) + 0.5*(Q_0(i,j,k,mx)**2 + Q_0(i,j,k,my)**2)/Q_0(i,j,k,rh)

!               rnum = (rand_number - 0.5)
!		Q_0(i,j,k,mx) = 3*rnum
!                rnum = (rand_number - 0.5)
!		Q_0(i,j,k,my) = 3*rnum
!                rnum = (rand_number - 0.5)
!		Q_0(i,j,k,mz) = 3*rnum

	        Q_r0(i,j,k,en,1) = wtev*Q_r0(i,j,k,rh,1)/(aindex - 1.) + 0.5*(Q_r0(i,j,k,mx,1)**2 + Q_r0(i,j,k,my,1)**2 + Q_r0(i,j,k,mz,1)**2)/Q_r0(i,j,k,rh,1)

                end do     
            end do
        end do
           
    end subroutine fill_fluid2

!----------------------------------------------------------------------------------------------

    subroutine set_bc
    implicit none
    real P1, P2, den, vz(nface)
    integer i,j,k,ne,ieq,l,i4


	if (mpi_P .eq. 1 .and. xlbc .ne. 2) then
		! Set B.C.'s at bottom x-boundary.

	do k = 1,nz	
	do j = 1,ny

		do i4=1,nface
		   Qxlow_ext(j,k,i4,:) = Qxlow_int(j,k,i4,:)
		end do


		if (maxval(Qxlow_ext(j,k,1:nface,mx)) .gt. 0.) then
		   Qxlow_ext(j,k,1:nface,mx) = 0.
		end if
			
	end do
	end do

	end if

!---------------------------------------------------------

	if (mpi_P .eq. mpi_nx .and. xlbc .ne. 2) then
		! Set B.C.'s at top x-boundary.

	do k = 1,nz	
	do j = 1,ny

		do i4=1,nface
		   Qxhigh_ext(j,k,i4,:) = Qxhigh_int(j,k,i4,:)
		end do
		
		if (minval(Qxhigh_ext(j,k,1:nface,mx)) .lt. 0.) then
		   Qxhigh_ext(j,k,1:nface,mx) = 0.
		end if		
	
	end do 
	end do

        end if

!----------------------------------------------------

	if (mpi_Q .eq. 1 .and. ylbc .ne. 2) then
		! Set B.C.'s at bottom y-boundary.
	
	do k = 1,nz	
	do i = 1,nx

		do i4=1,nface
		   Qylow_ext(i,k,i4,:) = Qylow_int(i,k,i4,:)
		end do
		
		if (maxval(Qylow_ext(i,k,1:nface,my)) .gt. 0.) then
		   Qylow_ext(i,k,1:nface,my) = 0.
		end if		
			
	end do
	end do

	end if

!------------------------------------------------------

	if (mpi_Q .eq. mpi_ny .and. ylbc .ne. 2) then
		! Set B.C.'s at top y-boundary.

	do k = 1,nz	
	do i = 1,nx

		do i4=1,nface
		   Qyhigh_ext(i,k,i4,:) = Qyhigh_int(i,k,i4,:)
		end do

		if (minval(Qyhigh_ext(i,k,1:nface,my)) .lt. 0.) then
		   Qyhigh_ext(i,k,1:nface,my) = 0.
		end if

	end do 
	end do

	end if

!--------------------------------------------------------

        if (mpi_R .eq. 1 .and. zlbc .ne. 2) then
		! Set B.C.'s at bottom z-boundary.
	
	do j = 1,ny	
	do i = 1,nx

            if (r(i,j) .ge. disk_rad .and. r(i,j) .le. lxu - 2.*dx) then

	      do i4=1,nface
		Qzlow_ext(i,j,i4,:) = Qzlow_int(i,j,i4,:)
	      end do
	    end if

		if (maxval(Qzlow_ext(i,j,1:nface,mz)) .gt. 0.) then
		   Qzlow_ext(i,j,1:nface,mz) = 0.
		end if
	
	end do
	end do

       end if

!-----------------------------------------------------------

	if (mpi_R .eq. mpi_nz .and. zlbc .ne. 2) then
		! Set B.C.'s at top z-boundary.

	do j = 1,ny	
	do i = 1,nx
	
                Qzhigh_ext(i,j,1:nface,:) = Qzhigh_int(i,j,1:nface,:)

		if (minval(Qzhigh_ext(i,j,1:nface,mz)) .lt. 0.) then
		   Qzhigh_ext(i,j,1:nface,mz) = 0
		end if	
			
	end do 
	end do

	end if

    end subroutine set_bc

!----------------------------------------------------------------------------------------------

    subroutine source_calc(Q_ri,t)
    implicit none
    integer i,j,k,ieq,ipg,ir
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri 	
    real t,source(npg,nQ),Qin(nQ),dn,dni,Zin,vx,vy,vz,alpha,temp,dne,eta_a,Teev,Tiev,etaJ2,nuei,TemR,Tev,vmax
    real Tcoef,fac,en_floor,gyro,Pres,rh_buf,oth

	source_r(:,:,:,:,:) = 0.0
	source(:,:) = 0.0
        en_floor = P_floor/(aindex - 1.)
        oth = 1./3.
	
	do k = 1,nz
	do j = 1,ny
	do i = 1,nx

	do ipg = 1,npg
	
	do ieq = 1,nQ
	   Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,k,ieq,1:nbasis))  
	end do
		
        dn = Qin(rh)
	dni = 1./Qin(rh)
        vx = Qin(mx)*dni
        vy = Qin(my)*dni
        vz = Qin(mz)*dni
        Pres = (aindex - 1.)*(Qin(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
        Tev = Pres*dni

        source(ipg,mx) = 0 
        source(ipg,my) = 0 
        source(ipg,mz) = 0  
	source(ipg,en) = 0

        if(ivis .eq. 1)then
           source(ipg,exx) = 0 
           source(ipg,eyy) = 0 
           source(ipg,ezz) = 0     
           source(ipg,exy) = 0   
           source(ipg,exz) = 0    
           source(ipg,eyz) = 0 
        end if
	if (ivis .eq. 0) then
           source(ipg,exx) = coll*dn*(2*vx**2 - vy**2 - vz**2)*oth
           source(ipg,eyy) = coll*dn*(2*vy**2 - vz**2 - vx**2)*oth
           source(ipg,ezz) = coll*dn*(2*vz**2 - vx**2 - vy**2)*oth 
           source(ipg,exy) = coll*dn*vx*vy  
           source(ipg,exz) = coll*dn*vx*vz  
           source(ipg,eyz) = coll*dn*vy*vz 
        end if

	end do

	do ieq = 1,nQ
	  do ir=1,nbasis
	     source_r(i,j,k,ieq,ir) = 0.125*cbasis(ir)*sum(wgt3d(1:npg)*bfvals_int(1:npg,ir)*source(1:npg,ieq))
	  end do
	end do

        source_r(i,j,k,rh:en,kxx:kyy) = -1.e-1*Q_ri(i,j,k,rh:en,kxx:kyy)
	
	end do
	end do
	end do
	
    end subroutine source_calc

!----------------------------------------------------------------------------------------------

    subroutine advance_time_level_gl(Q_ri,Q_rp)
    implicit none
    integer i,j,k,ieq,ir
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_ri, Q_rp
    real Q_xtemp, Q_ytemp, Q_ztemp, dti,oth
    real c1d3, c1d5

	oth = 1./3.
        dti = 1./dt

	do k = 1,nz
	do j = 1,ny
	do i = 1,nx

	do ieq = 1,nQ
	  do ir=1,nbasis
	    Q_rp(i,j,k,ieq,ir) = Q_ri(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir) + dt*source_r(i,j,k,ieq,ir)
	  end do
	end do

        if(ivis .eq. 0)then

	do ir=1,nbasis

        Q_rp(i,j,k,exx,ir) = (Q_ri(i,j,k,exx,ir) - dt*glflux_r(i,j,k,exx,ir) + dt*source_r(i,j,k,exx,ir)                                 &
                           + coll*dt*(Q_ri(i,j,k,exx,ir) + Q_ri(i,j,k,eyy,ir) + Q_ri(i,j,k,ezz,ir))*oth                                  &
                           - coll*dt**2*(glflux_r(i,j,k,exx,ir) + glflux_r(i,j,k,eyy,ir) + glflux_r(i,j,k,ezz,ir))*oth                   &
                           + coll*dt**2*(source_r(i,j,k,exx,ir) + source_r(i,j,k,eyy,ir) + source_r(i,j,k,ezz,ir))*oth)/(1. + dt*coll)

        Q_rp(i,j,k,eyy,ir) = (Q_ri(i,j,k,eyy,ir) - dt*glflux_r(i,j,k,eyy,ir) + dt*source_r(i,j,k,eyy,ir)                                 &
                           + coll*dt*(Q_ri(i,j,k,exx,ir) + Q_ri(i,j,k,eyy,ir) + Q_ri(i,j,k,ezz,ir))*oth                                  &
                           - coll*dt**2*(glflux_r(i,j,k,exx,ir) + glflux_r(i,j,k,eyy,ir) + glflux_r(i,j,k,ezz,ir))*oth                   &
                           + coll*dt**2*(source_r(i,j,k,exx,ir) + source_r(i,j,k,eyy,ir) + source_r(i,j,k,ezz,ir))*oth)/(1. + dt*coll)

        Q_rp(i,j,k,ezz,ir) = (Q_ri(i,j,k,ezz,ir) - dt*glflux_r(i,j,k,ezz,ir) + dt*source_r(i,j,k,ezz,ir)                                 &
                           + coll*dt*(Q_ri(i,j,k,exx,ir) + Q_ri(i,j,k,eyy,ir) + Q_ri(i,j,k,ezz,ir))*oth                                  &
                           - coll*dt**2*(glflux_r(i,j,k,exx,ir) + glflux_r(i,j,k,eyy,ir) + glflux_r(i,j,k,ezz,ir))*oth                   &
                           + coll*dt**2*(source_r(i,j,k,exx,ir) + source_r(i,j,k,eyy,ir) + source_r(i,j,k,ezz,ir))*oth)/(1. + dt*coll)

        end do

	do ieq = exy,nQ
	  do ir=1,nbasis
	    Q_rp(i,j,k,ieq,ir) = (Q_ri(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir) + dt*source_r(i,j,k,ieq,ir))/(1. + dt*coll)
	  end do
	end do

        end if

        if(ivis .eq.1)then

	do ieq = exx,nQ
	  do ir=1,nbasis
	    Q_rp(i,j,k,ieq,ir) = (Q_ri(i,j,k,ieq,ir) - dt*glflux_r(i,j,k,ieq,ir) + dt*source_r(i,j,k,ieq,ir))/(1. + dt*coll)
	  end do
	end do

        end if

        do ieq = 1,nQ
        if ( Q_rp(i,j,k,ieq,1) .ne. Q_rp(i,j,k,ieq,1)) then
	        print *,'NaN. Bailing out...','  xc  =',xc(i),'  yc  =',yc(j),'  zc  =',zc(k),'  ieq  =',ieq
        call exit(-1)
        endif
        end do

	end do
	end do
	end do
	
    end subroutine advance_time_level_gl

!----------------------------------------------------------------------------------------------

    subroutine flux_calc_pnts_r5(Qpnts_r,fpnts_r,ixyz,npnts)

	! Calculate the flux "fpnts_r" in direction "ixyz" (x, y, or z) at a set of 
	! points corresponding to conserved quantities "Qpnts_r".

	! ixyz=1: x-direction
	! ixyz=2: y-direction
	! ixyz=3: z-direction	

    implicit none
    integer ife, ixyz,npnts
    real, dimension(npnts,nQ):: Qpnts_r, fpnts_r
    real dn,dni,vr,P,asqr,fac,Pre,dnei,Psol,dx2,Tem,smsq

	do ife = 1,npnts

	dn = Qpnts_r(ife,rh)
	dni = 1./dn

	smsq = Qpnts_r(ife,mx)**2 + Qpnts_r(ife,my)**2 + Qpnts_r(ife,mz)**2
	vr = Qpnts_r(ife,mxa(ixyz))*dni

	P = (aindex - 1.)*(Qpnts_r(ife,en) - 0.5*dni*smsq) 
	if (P < P_floor) P = P_floor

	fpnts_r(ife,rh) = Qpnts_r(ife,mxa(ixyz))
	fpnts_r(ife,mx:mz) = Qpnts_r(ife,mx:mz)*vr
	fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vr

	fpnts_r(ife,mxa(ixyz)) = fpnts_r(ife,mxa(ixyz)) + P 	 
		
	end do

    end subroutine

!----------------------------------------------------------------------------------------------

    subroutine flux_calc_pnts_r(Qpnts_r,fpnts_r,ixyz,npnts)

	! Calculate the flux "fpnts_r" in direction "ixyz" (x, y, or z) at a set of 
	! points corresponding to conserved quantities "Qpnts_r".

	! ixyz=1: x-direction
	! ixyz=2: y-direction
	! ixyz=3: z-direction	

	! For ivis = 0 the fluxes for the 10-moment equations are computed based on Grad's 13-moment approximation. 
	! The energy tensor components are E_ij = S_ij + rho u_i u_j + delta_ij P, where S_ij is the deviatoric stress tensor
	! and P is the scalar pressure.

	! For ivis = 1 the fluxes are for the linearized stress tensor in which case the equations solved are:
	! partial_t S_ij + P_0 (partial_x_i u_j + partial_x_j u_i - (2/3)delta_ij partial_x_k u_k) = -coll S_ij, where
	! P_0 is the avarage pressure and coll = P_0/mu, where mu is the shear viscosity.

    implicit none
    integer ife, ixyz,npnts
    real, dimension(npnts,nQ):: Qpnts_r, fpnts_r
    real dn,dni,vx,vy,vz,P,asqr,fac,Pre,dnei,Psol,dx2,Tem,nu,c2d3,c4d3
    real E_xx,E_yy,E_zz,E_xy,E_xz,E_yz,P_10,P_5
    real ampx,ampy,ampz,ampd

        nu = rh_fluid*T_floor

	c2d3 = 2./3.
	c4d3 = 4./3.

	do ife = 1,npnts

	dn = Qpnts_r(ife,rh) 
	dni = 1./dn
	vx = Qpnts_r(ife,mx)*dni
	vy = Qpnts_r(ife,my)*dni
	vz = Qpnts_r(ife,mz)*dni
        E_xx = Qpnts_r(ife,exx)
        E_yy = Qpnts_r(ife,eyy)
        E_zz = Qpnts_r(ife,ezz)
        E_xy = Qpnts_r(ife,exy)
        E_xz = Qpnts_r(ife,exz)
        E_yz = Qpnts_r(ife,eyz)
	
	P = (aindex - 1.)*(Qpnts_r(ife,en) - 0.5*dn*(vx**2 + vy**2 + vz**2)) 
        P_10 = (E_xx + E_yy + E_zz - dn*(vx**2 + vy**2 + vz**2))/3.
	if (P < P_floor) P = P_floor
	if (P_10 < P_floor) P_10 = P_floor
        P_5 = P

	if (ivis .lt. 0) then
           if(ixyz .eq. 1)then
	      fpnts_r(ife,rh) = Qpnts_r(ife,mx) 
	      fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vx + P
	      fpnts_r(ife,my) = Qpnts_r(ife,my)*vx    
	      fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vx 
	      fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vx
	   end if
           if(ixyz .eq. 2)then
	      fpnts_r(ife,rh) = Qpnts_r(ife,my) 
	      fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vy
	      fpnts_r(ife,my) = Qpnts_r(ife,my)*vy + P    
	      fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vy 
	      fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vy
	   end if
           if(ixyz .eq. 3)then
	      fpnts_r(ife,rh) = Qpnts_r(ife,mz) 
	      fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vz
	      fpnts_r(ife,my) = Qpnts_r(ife,my)*vz    
	      fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vz + P 
	      fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vz
	   end if
	end if


	if (ivis .ge. 0) then
        if(ixyz .eq. 1)then

	fpnts_r(ife,rh) = Qpnts_r(ife,mx) 

        if(ivis .eq. 1)then

	fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vx + E_xx + P
	fpnts_r(ife,my) = Qpnts_r(ife,my)*vx + E_xy    
	fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vx + E_xz 

        else

	fpnts_r(ife,mx) = E_xx - P_10 + P_5
	fpnts_r(ife,my) = E_xy   
	fpnts_r(ife,mz) = E_xz   

        end if

	fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vx !+ (vx*E_xx + vy*E_xy + vz*E_xz) 

        if(ivis .eq. 1)then

	fpnts_r(ife,exx) =  c4d3*nu*vx 
	fpnts_r(ife,eyy) = -c2d3*nu*vx 
	fpnts_r(ife,ezz) = -c2d3*nu*vx 

	fpnts_r(ife,exy) = nu*vy
	fpnts_r(ife,exz) = nu*vz
	fpnts_r(ife,eyz) = 0

        else

	fpnts_r(ife,exx) = vx*(3*E_xx - 2*vx*Qpnts_r(ife,mx))                                       ! term 1
	fpnts_r(ife,eyy) = 2*vy*(E_xy - vy*Qpnts_r(ife,mx)) + vx*E_yy                               ! term 4
	fpnts_r(ife,ezz) = 2*vz*(E_xz - vz*Qpnts_r(ife,mx)) + vx*E_zz                               ! term 7

	fpnts_r(ife,exy) = 2*vx*(E_xy - vy*Qpnts_r(ife,mx)) + vy*E_xx                               ! term 10
	fpnts_r(ife,exz) = 2*vx*(E_xz - vz*Qpnts_r(ife,mx)) + vz*E_xx                               ! term 13
	fpnts_r(ife,eyz) = vx*E_yz + vy*E_xz + vz*E_xy - 2*vy*vz*Qpnts_r(ife,mx)                    ! term 16

        end if

        end if

        if(ixyz .eq. 2)then

	fpnts_r(ife,rh) = Qpnts_r(ife,mxa(ixyz)) 

        if(ivis .eq. 1)then

	fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vy + E_xy    
	fpnts_r(ife,my) = Qpnts_r(ife,my)*vy + E_yy + P
	fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vy + E_yz 

        else

	fpnts_r(ife,mx) = E_xy     
	fpnts_r(ife,my) = E_yy - P_10 + P_5
	fpnts_r(ife,mz) = E_yz 

        end if

	fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vy !+ (vy*E_yy + vx*E_xy + vz*E_yz) 

        if(ivis .eq. 1)then

	fpnts_r(ife,exx) = -c2d3*nu*vy
	fpnts_r(ife,eyy) =  c4d3*nu*vy
	fpnts_r(ife,ezz) = -c2d3*nu*vy

	fpnts_r(ife,exy) = nu*vx
	fpnts_r(ife,exz) = 0
	fpnts_r(ife,eyz) = nu*vz

        else

	fpnts_r(ife,exx) = 2*vx*(E_xy - vx*Qpnts_r(ife,my)) + vy*E_xx                            ! term 2
	fpnts_r(ife,eyy) = vy*(3*E_yy - 2*vy*Qpnts_r(ife,my))                                    ! term 5
	fpnts_r(ife,ezz) = 2*vz*(E_yz - vz*Qpnts_r(ife,my)) + vy*E_zz                            ! term 8

	fpnts_r(ife,exy) = 2*vy*(E_xy - vx*Qpnts_r(ife,my)) + vx*E_yy                            ! term 11
	fpnts_r(ife,exz) = vx*E_yz + vy*E_xz + vz*E_xy - 2*vx*vz*Qpnts_r(ife,my)                 ! term 14
	fpnts_r(ife,eyz) = 2*vy*(E_yz - vz*Qpnts_r(ife,my)) + vz*E_yy                            ! term 17

        end if

        end if

        if(ixyz .eq. 3)then

	fpnts_r(ife,rh) = Qpnts_r(ife,mz) 

        if(ivis .eq. 1)then

	fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vz + E_xz   
	fpnts_r(ife,my) = Qpnts_r(ife,my)*vz + E_yz    
	fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vz + E_zz + P 

        else

	fpnts_r(ife,mx) = E_xz    
	fpnts_r(ife,my) = E_yz    
	fpnts_r(ife,mz) = E_zz - P_10 + P_5

        end if

	fpnts_r(ife,en) = (Qpnts_r(ife,en) + P)*vz !+ (vz*E_zz + vx*E_xz + vy* E_yz) 

        if(ivis .eq. 1)then

	fpnts_r(ife,exx) = -c2d3*nu*vz
	fpnts_r(ife,eyy) = -c2d3*nu*vz
	fpnts_r(ife,ezz) =  c4d3*nu*vz

	fpnts_r(ife,exy) = 0.
	fpnts_r(ife,exz) = nu*vx
	fpnts_r(ife,eyz) = nu*vy

        else

	fpnts_r(ife,exx) = 2*vx*(E_xz - vx*Qpnts_r(ife,mz)) + vz*E_xx                               ! term 3
	fpnts_r(ife,eyy) = 2*vy*(E_yz - vy*Qpnts_r(ife,mz)) + vz*E_yy                               ! term 6
	fpnts_r(ife,ezz) = vz*(3*E_zz - 2*vz*Qpnts_r(ife,mz))                                       ! term 9

	fpnts_r(ife,exy) = vx*E_yz + vy*E_xz + vz*E_xy - 2*vx*vy*Qpnts_r(ife,mz)                    ! term 12
	fpnts_r(ife,exz) = 2*vz*(E_xz - vx*Qpnts_r(ife,mz)) + vx*E_zz                               ! term 15
	fpnts_r(ife,eyz) = 2*vz*(E_yz - vy*Qpnts_r(ife,mz)) + vy*E_zz                               ! term 18

        end if

        end if
	end if
		
	end do

    end subroutine

!----------------------------------------------------------------------------------------------

    subroutine flux_cal(Q_r)
    implicit none
    integer i,j,k,ieq,iback,jleft,kdown,i4,i4p,ipnt
    real Qface_x(nfe,nQ), Qface_y(nfe,nQ),Qface_z(nfe,nQ),fface_x(nfe,nQ), fface_y(nfe,nQ), fface_z(nfe,nQ)
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
    real cwavex(nfe),cwavey(nfe),cwavez(nfe),Pre,dni,B2, cfbound, cfx(nface,nQ),cfy(nface,nQ),cfz(nface,nQ)
    real fhllc_x(nface,5),fhllc_y(nface,5),fhllc_z(nface,5),fs(nface,nQ),qvin(nQ)

!--------------------------------------------------------------

        kroe(:) = 1

	do k=1,nz
	do j=1,ny
	do i=1,nx+1

	iback = i-1
	
	if (i .gt. 1) then
	do ieq = 1,nQ
        do ipnt=1,nface	
		Qface_x(ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_r(iback,j,k,ieq,1:nbasis))
        end do
	end do
	end if

	if (i .eq. 1) then
	do ieq = 1,nQ	
		Qface_x(1:nface,ieq) = Qxlow_ext(j,k,1:nface,ieq)
	end do
	end if

	if (i .lt. nx+1) then
	do ieq = 1,nQ	
	  do ipnt=1,nface
		Qface_x(ipnt+nface,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
	  end do
	end do
	end if

	if (i .eq. nx+1) then
	do ieq = 1,nQ	
		Qface_x(nface+1:nfe,ieq) = Qxhigh_ext(j,k,1:nface,ieq)
	end do
	end if

	call flux_calc_pnts_r(Qface_x,fface_x,1,nfe)

	do i4=1,nfe
	do ieq=1,nQ
	      qvin(ieq) = Qface_x(i4,ieq)
	end do
        cwavex(i4) = cfcal(qvin,1)
        end do

        do i4=1,nface
           cfrx(i4,rh:en) = max(cwavex(i4),cwavex(i4+nface))
        end do

	do ieq = 1,nQ
	 do i4=1,nface
	   i4p = i4 + nface
	   flux_x(i4,i,j,k,ieq) = 0.5*(fface_x(i4,ieq) + fface_x(i4p,ieq)) - 0.5*cfrx(i4,ieq)*(Qface_x(i4p,ieq) - Qface_x(i4,ieq))
	end do
	end do

	kroe(1:nface) = 1

	if (ihllc .eq. 1) call flux_hllc(Qface_x,fface_x,fhllc_x,1)
	if (iroe .eq. 1) call flux_roe(Qface_x,fface_x,fhllc_x,1)

        if (ihllc .eq. 1 .or. iroe .eq. 1) then

	do ieq = 1,en
	   do i4=1,nface
	   if (kroe(i4) .gt. 0) then
	     flux_x(i4,i,j,k,ieq) = fhllc_x(i4,ieq)
	   end if
	   end do 
	end do

        end if

	end do
	end do
	end do

!----------------------------------------------------

	do k=1,nz
	do j=1,ny+1
		jleft = j-1
	do i=1,nx

	if (j .gt. 1) then
	do ieq = 1,nQ
		do ipnt=1,nface
	          Qface_y(ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_r(i,jleft,k,ieq,1:nbasis))
		end do
	end do
	end if
	if (j .eq. 1) then
	do ieq = 1,nQ
	   Qface_y(1:nface,ieq) = Qylow_ext(i,k,1:nface,ieq)
	end do
	end if

	if (j .lt. ny+1) then
	do ieq = 1,nQ
	    do ipnt=1,nface
		Qface_y(ipnt+nface,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
	    end do
	end do
	end if
	if (j .eq. ny+1) then
	do ieq = 1,nQ
	   Qface_y(nface+1:nfe,ieq) = Qyhigh_ext(i,k,1:nface,ieq)
	end do
	end if

	   call flux_calc_pnts_r(Qface_y,fface_y,2,nfe)

	   do i4=1,nfe
	      do ieq=1,nQ
	         qvin(ieq) = Qface_y(i4,ieq)
	      end do
              cwavey(i4) = cfcal(qvin,2)
       end do

           do i4=1,nface
              cfry(i4,rh:en) = max(cwavey(i4),cwavey(i4+nface))            
           end do

	do ieq = 1,nQ
	 do i4=1,nface
	   i4p = i4 + nface
	   flux_y(i4,i,j,k,ieq) = 0.5*(fface_y(i4,ieq) + fface_y(i4p,ieq)) - 0.5*cfry(i4,ieq)*(Qface_y(i4p,ieq) - Qface_y(i4,ieq))
	end do
	end do

	kroe(1:nface) = 1

	if (ihllc .eq. 1) call flux_hllc(Qface_y,fface_y,fhllc_y,2)
	if (iroe .eq. 1) call flux_roe(Qface_y,fface_y,fhllc_y,2)

        if (ihllc .eq. 1 .or. iroe .eq. 1) then

	do ieq = 1,en
	   do i4=1,nface
	   if (kroe(i4) .gt. 0) then
	     flux_y(i4,i,j,k,ieq) = fhllc_y(i4,ieq)
	   end if
	   end do 
	end do

        end if

	end do
	end do
	end do

!-------------------------------------------------------

	do k=1,nz+1
	kdown = k-1
	do j=1,ny
	do i=1,nx


	if (k .gt. 1) then
	do ieq = 1,nQ
		do ipnt=1,nface
			Qface_z(ipnt,ieq) = sum(bfvals_zp(ipnt,1:nbasis)*Q_r(i,j,kdown,ieq,1:nbasis))
		end do
	end do
	end if
	if (k .eq. 1) then
	do ieq = 1,nQ
	   Qface_z(1:nface,ieq) = Qzlow_ext(i,j,1:nface,ieq)
	end do
	end if

	if (k .lt. nz+1) then
	do ieq = 1,nQ
		do ipnt=1,nface
			Qface_z(ipnt+nface,ieq) = sum(bfvals_zm(ipnt,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
		end do
	end do
	end if
	if (k .eq. nz+1) then
	do ieq = 1,nQ
	   Qface_z(nface+1:nfe,ieq) = Qzhigh_ext(i,j,1:nface,ieq)
	end do
	end if

	   call flux_calc_pnts_r(Qface_z,fface_z,3,nfe)

	do i4=1,nfe
	   do ieq=1,nQ
	      qvin(ieq) = Qface_z(i4,ieq)
	   end do
           cwavez(i4) = cfcal(qvin,3)
        end do

        do i4=1,nface           
           cfrz(i4,rh:en) = max(cwavez(i4),cwavez(i4+nface))
        end do

	do ieq = 1,nQ
	 do i4=1,nface
	   i4p = i4 + nface
	   flux_z(i4,i,j,k,ieq) = 0.5*(fface_z(i4,ieq) + fface_z(i4p,ieq)) - 0.5*cfrz(i4,ieq)*(Qface_z(i4p,ieq) - Qface_z(i4,ieq))
	end do
	end do

	kroe(1:nface) = 1

	if (ihllc .eq. 1) call flux_hllc(Qface_z,fface_z,fhllc_z,3)
	if (iroe .eq. 1) call flux_roe(Qface_z,fface_z,fhllc_z,3)

        if (ihllc .eq. 1 .or. iroe .eq. 1) then

	do ieq = 1,en
	   do i4=1,nface
	   if (kroe(i4) .gt. 0) then
	     flux_z(i4,i,j,k,ieq) = fhllc_z(i4,ieq)
	   end if
	   end do 
	end do

        end if
		
	end do
	end do
	end do

    end subroutine flux_cal

!----------------------------------------------------------------------------------------------

    subroutine flux_hllc(Qlr,flr,fhllc,ixyz)

!	Compute ion or electron fluxes (density, momentum, energy) using 
!	nonrelativistic HD HLLC approximate Riemann solver developed by Batten, 1997,
!	*****    "On the Choice of Wavespeeds for the HLLC Riemann Solver"

!	This takes into account the left and right propagating shocks along with 
!	the contact (tangential) discontinuity.
	
	implicit none
	real Qlr(nfe,nQ),flr(nfe,nQ),fhllc(nface,5)
	real sm_num(nface),sm_den(nface),qtilde(nface,5),rtrho(nfe),rtrho_i(nface),qsq(nfe)
	real s_lr(nfe),ctilde(nface),hlr(nfe),cslr(nfe),ctsq(nface),Zi,mfact
	real aq(nfe),bq(nfe),Qstar(nfe,6),fstar(nfe,6),pstar(nface),s_m(nface)
	real rhov(nfe),vlr(nfe,3),plr(nfe),slsm_i(nface),rho_i,qslr(nfe),sq_lr(nfe),slrm_i(nfe),B2(nfe),cf(nfe)
	integer ixyz,i4,i4p1,nr,jie,k,k2,ieq,iparr,iperp1,iperp2,ibatten
	integer rhj,mxj,myj,mzj,enj,psj,ivar(5),ipassive,nhll,ib1,ib2

        nhll = 5
	rhj = rh
	mxj = mx
	myj = my
	mzj = mz
	enj = en

	iparr  = mxa(ixyz)
	iperp1 = mya(ixyz)
	iperp2 = mza(ixyz)
	
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
	      vlr(k,1) = Qlr(k,iparr)*rho_i		! velocity parallel to direction of flux computation
	      vlr(k,2) = Qlr(k,iperp1)*rho_i		! velocity in perpendicular direction 1
	      vlr(k,3) = Qlr(k,iperp2)*rho_i		! velocity in perpendicular direction 2
	      qsq(k) = vlr(k,1)**2 + vlr(k,2)**2 + vlr(k,3)**2
	      plr(k) = (aindex - 1.0)*(Qlr(k,enj) - 0.5*rhov(k)*qsq(k))		! pressure
	      rtrho(k) = sqrt(rhov(k))
	   end do

	   do k=1,nface
	      k2 = k + nface
	      cslr(k) = vlr(k,1) - sqrt(aindex*plr(k)/rhov(k))   	! lambda_M(Q_l)
	      cslr(k2) = vlr(k2,1) + sqrt(aindex*plr(k2)/rhov(k2) )   	! lambda_P(Q_r)
	   end do

	   if (ibatten .eq. 1) then		! compute wave speeds using Roe averages following Batten, 1997

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
	      ctsq(k) = (aindex - 1.0)*(qtilde(k,4) - 0.5*qsq(k))
	    end do
	      if (minval(ctsq) .ge. 0.0) then
		ctilde = sqrt(ctsq)
		qslr(1:nface) = qtilde(1:nface,1) - ctilde(1:nface)   	! lambda_M(Q_Roe)
	        qslr(nface+1:nfe) = qtilde(nface+1:nfe,1) + ctilde(nface+1:nfe)    ! lambda_P(Q_Roe)
	      end if
	      if (minval(ctsq) .lt. 0.0) then
		ibatten = 0
	      end if
	   end if

	   if (ibatten .eq. 0) then
	    do k=1,nface
	      k2 = k + nface
	      qslr(k) = vlr(k2,1) - sqrt(aindex*plr(k2)/rhov(k2))   	! lambda_M(Q_r)
	      qslr(k2) = vlr(k,1) + sqrt(aindex*plr(k)/rhov(k))   	! lambda_P(Q_l)
	    end do
	   end if

!	  Calculate the slow and fast wavespeeds S_L, S_R of the left, right propagating shocks.

	    do k=1,nface
	      k2 = k + nface
	      s_lr(k) = min(cslr(k),qslr(k))         ! S_L = min(lambda_M(Q_l),lambda_M(Q_r or Q_Roe))
	      s_lr(k2) = max(cslr(k2),qslr(k2))	     ! S_R = max(lambda_P(Q_r),lambda_P(Q_l or Q_Roe))
	      sm_num(k) = rhov(k2)*vlr(k2,1)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*vlr(k,1)*(s_lr(k) - vlr(k,1))
	      sm_num(k) = sm_num(k) + plr(k) - plr(k2)
	      sm_den(k) = rhov(k2)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*(s_lr(k) - vlr(k,1))
	    end do
	   where (sm_den .eq. 0.0) sm_den = rh_floor

!	   Calculate the wavespeed S_M of the contact discontinuity.
	
	   do k=1,nface
	      s_m(k) = sm_num(k)/sm_den(k)	     	! Eq. (34) of Batten, 1997
	      pstar(k) = rhov(k)*(vlr(k,1) - s_lr(k))*(vlr(k,1) - s_m(k)) + plr(k) 		
							! Eq. (36) of Batten, 1997
	   end do


!	   Now, calculate Q_l* and Q_r* in order to calculate F_l* and F_r*.

	   do k=1,nfe
	      if (k .le. nface) i4 = k
	      if (k .gt. nface) i4 = k - nface
	      sm_den(1) = s_lr(k) - s_m(i4)		! S_{L,R} - S_M

	      if (sm_den(1) .eq. 0.0) then
	         sm_den(1) = rh_floor
	      end if

	      slrm_i(k) = 1.0/sm_den(1)
	      sq_lr(k) = s_lr(k) - vlr(k,1)			! S_{L,R} - q_{l,r}
	      Qstar(k,1) = rhov(k)*sq_lr(k)*slrm_i(k)		! Eq. (35) of Batten ,1997			
	      Qstar(k,2) = (sq_lr(k)*Qlr(k,iparr) + pstar(i4) - plr(k))*slrm_i(k) 	! Eq. (37-39) of Batten, 1997
	      Qstar(k,3) = sq_lr(k)*Qlr(k,iperp1)*slrm_i(k)		! Eq. (37-39) of Batten, 1997
	      Qstar(k,4) = sq_lr(k)*Qlr(k,iperp2)*slrm_i(k)		! Eq. (37-39) of Batten, 1997
	      Qstar(k,5) = (sq_lr(k)*Qlr(k,enj) - plr(k)*vlr(k,1) + pstar(i4)*s_m(i4))*slrm_i(k) 	
								! Eq. (40) of Batten, 1997

	      do ieq=1,nhll
	        fstar(k,ieq) = flr(k,ivar(ieq)) + s_lr(k)*(Qstar(k,ieq) - Qlr(k,ivar(ieq)))		
								! Eq. (29) of Batten, 1997
	      end do

	   end do

	! Finally, calculate the HLLC fluxes from F_l, F_l*, F_r*, and F_r...

	do i4 = 1,nface						! Use Eq. (26) of Batten ,1997
	    if (s_lr(i4) .gt. 0.0) then				! if S_L > 0
	      do ieq=1,nhll
		fhllc(i4,ivar(ieq)) = flr(i4,ivar(ieq))		! F_HLLC = F_l
	      end do
	    end if
	    if (s_lr(i4) .le. 0.0 .and. 0.0 .lt. s_m(i4)) then		! if S_L <= 0 < S_M
	      do ieq=1,nhll
	         fhllc(i4,ivar(ieq)) = fstar(i4,ieq)		! F_HLLC = F_l*
	      end do
	    end if
	    if (s_m(i4) .le. 0.0 .and. 0.0 .le. s_lr(i4+nface)) then    	! if S_M <= 0 <= S_R
	      do ieq=1,nhll
		fhllc(i4,ivar(ieq)) = fstar(i4+nface,ieq)		! F_HLLC = F_r*
	      end do
	    end if
	    if (s_lr(i4+nface) .lt. 0.0) then			! if S_R < 0
	      do ieq=1,nhll
		fhllc(i4,ivar(ieq)) = flr(i4+nface,ivar(ieq))	! F_HLLC = F_r
	      end do
	    end if

	end do

	end subroutine flux_hllc

	!----------------------------------------------------------------

    subroutine flux_roe(Qlr,flr,froef,ixyz)

    implicit none
    integer ixyz,i9,j9,ilr,iparr,iperp
    real Qlr(nfe,nQ),flr(nfe,nQ),froef(nface,5)
    real dflux(nface,5),dwr(nface,5),dsk(nface),Qri
    real evec(5,5),swr(2,nface,5),skr(2,nface),asq,lam9(5),a9(5),sk9
    real a9e,ea2i,vsq,Pri,vels(3),dnii,en_floor
    real vc(3),hc,dele(nface),delrho(nface),delmom1(nface),delmom2(nface),delmom3(nface),sc9,vcsq,vcdel,csq,gmi,bsq


	! approximate Roe solver for non-relativistic two-fluid, from Eulderink and Mellema, 1994 

	! Starting from HD equations for a given fluid (either electrons or ions), this routine 
	! computes "dflux_i = (cQ)_i+1 - (cQ)_i" terms by characteristic wave decomposition using 
	! a Roe matrix.

	! INPUT is Qr(,1)=conserved density, Qr(,2)=energy, Qr(,3:5)=momentum
	!	   Qpr(,1)=density, Qpr(,2)=pressure, Qpr(,3:5)=velocity
	!	   n = number of spatial points in given direction (x or y) 

	! 	   ixyz = 1: flux is being computed in x-direction
	! 	   ixyz = 2: flux is being computed in y-direction


	! OUTPUT is  "dflux_i = (cQ)_i+1 - (cQ)_i"
	!	dflux(,1)=density term, dflux(,2)=energy term, dflux(,3:5)=momentum terms



	! evec(,1) through evec(,5) are the eigenvectors of the Roe matrix.

	en_floor = P_floor*aindm1
        gmi = 1./(aindex - 1.)
	
	evec(:,:) = 0.
	evec(1,1) = 1.
	evec(1,2) = 1.
	evec(1,3) = 1.
	evec(4,4) = 1.
	evec(5,5) = 1.


	do i9=1,nface     
	  kroe(i9) = 1
	  do j9=1,2
	    if (j9 .eq. 1) ilr = i9
	    if (j9 .eq. 2) ilr = i9 + nface		

		! Make sure density, pressure, and energy are at least at floor values.

	     dnii = 1./Qlr(ilr,rh)
	     vels(1) = Qlr(ilr,mxa(ixyz))*dnii
	     vels(2) = Qlr(ilr,mya(ixyz))*dnii
	     vels(3) = Qlr(ilr,mza(ixyz))*dnii
	     vsq = vels(1)**2 + vels(2)**2 + vels(3)**2
	     Pri = aindm1*(Qlr(ilr,en) - 0.5*Qlr(ilr,rh)*vsq)

             skr(j9,i9) = sqrt(Qlr(ilr,rh))	 
	     swr(j9,i9,1) = skr(j9,i9)*vels(1)  ! sqrt(rho) * v_x
	     swr(j9,i9,2) = skr(j9,i9)*vels(2)  ! sqrt(rho) * v_y
	     swr(j9,i9,3) = skr(j9,i9)*vels(3)  ! sqrt(rho) * v_z
	     swr(j9,i9,4) = 0.5*skr(j9,i9)*(vsq + cp*Pri/Qlr(ilr,rh))
	   end do
	end do

	do i9=1,nface
	     Qri = 1./Qlr(i9,rh)	! Increments in conserved quantities are normalized w.r.t. density.

             delrho(i9) = Qlr(i9+nface,rh)*Qri - 1.	! delrho = increment in conserved density
             dele(i9) = (Qlr(i9+nface,en) - Qlr(i9,en))*Qri    ! *one_mime(jie)
			! dele = increment in conserved energy

             delmom1(i9) = (Qlr(i9+nface,mxa(ixyz)) - Qlr(i9,mxa(ixyz)))*Qri	! del1 = increment in x-momentum
             delmom2(i9) = (Qlr(i9+nface,mya(ixyz)) - Qlr(i9,mya(ixyz)))*Qri	! del2 = increment in y-momentum
             delmom3(i9) = (Qlr(i9+nface,mza(ixyz)) - Qlr(i9,mza(ixyz)))*Qri	! del3 = increment in z-momentum

	!    dwr(i,1:3) = 0.5*[sqrt(rho_i)v_{1:3,i} + sqrt(rho_{i+1})v_{1:3,i+1}]
	!    dwr(i,4) = 0.5*[sqrt(rho_i) enthalpy(i)/rho(i) + sqrt(rho_{i+1}) enthalpy(i+1)/rho(i+1)]

	   do j9=1,4
	     dwr(i9,j9) = 0.5*(swr(2,i9,j9) + swr(1,i9,j9))
	   enddo

	   dsk(i9) = 0.5*(skr(2,i9) + skr(1,i9))  

	enddo

	! The Roe average of a quantity is the arithmetic average 
	! between neighboring cells weighted by the square root of density.
	! For example, for "v_x" at position "i" the Roe average "v_{cx}" is

	!    v_{cx} = [sqrt(rho_i)v_{xi} + sqrt(rho_{i+1}) v_{x,i+1}]/[sqrt(rho_i) + sqrt(rho_{i+1})]

	do i9=1,nface	
	  vc(1) = dwr(i9,1)/dsk(i9)	! component 1 of Roe-averaged velocity (x if jie=1, y if jie=2)
	  vc(2) = dwr(i9,2)/dsk(i9)	! component 2 of Roe-averaged velocity (y if jie=1, x if jie=2)
	  vc(3) = dwr(i9,3)/dsk(i9)	! component 3 of Roe-averaged velocity (z-component)
	  hc = dwr(i9,4)/dsk(i9)	! Roe-averaged enthalpy/density

	  vcsq = vc(1)*vc(1) + vc(2)*vc(2) + vc(3)*vc(3)
	  asq = aindm1*(hc - 0.5*vcsq)	! squared sound speed
	  if (asq .le. 0.0) then
		kroe(i9) = 0    ! asq = (aindex - 1.0d0)*hc
	  end if
	  sc9 = sqrt(asq)	! sound speed


	! Define the characteristic speeds (eigenvalues of the Roe matrix).
             lam9(1) = abs(vc(1) - sc9)
	     lam9(2) = abs(vc(1) + sc9)
	     lam9(3) = abs(vc(1))
	     lam9(4) = lam9(3)
	     lam9(5) = lam9(3)

	! Define the eigenvectors evec(,1)...evec(,5) of the Roe matrix.
	  evec(2,1) = hc - sc9*vc(1)
	  evec(3,1) = vc(1) - sc9
	  evec(4,1) = vc(2)
	  evec(5,1) = vc(3)

	  evec(2,2) = hc + sc9*vc(1)
	  evec(3,2) = vc(1) + sc9
	  evec(4,2) = vc(2)
	  evec(5,2) = vc(3)

	  evec(2,3) = 0.5*vcsq
	  evec(3,3) = vc(1)
	  evec(4,3) = vc(2)
	  evec(5,3) = vc(3)
	  evec(2,4) = vc(2)
	  evec(2,5) = vc(3)

	! Define a few intermediate variables needed for computing the expansion coefficients.

	  ea2i = aindm1/asq
	  vcdel = vc(1)*delmom1(i9) + vc(2)*delmom2(i9) + vc(3)*delmom3(i9) - dele(i9)
	  a9e = 0.5*ea2i*(0.5*vcsq*delrho(i9) - vcdel)
	  sk9 = 0.5*(delmom1(i9) - vc(1)*delrho(i9))/sc9

	 ! Define the expansion coefficients a9_1...a9_5 such that   
	 ! 	Delta Q = {delrho,dele,del1,del2,del3} = sum [a9_j evec(,j)]

	  a9(1) = a9e - sk9
	  a9(2) = a9e + sk9
	  a9(3) = ea2i*((hc - vcsq)*delrho(i9) + vcdel)
	  a9(4) = delmom2(i9) - vc(2)*delrho(i9)
	  a9(5) = delmom3(i9) - vc(3)*delrho(i9)

	 ! The flux increments "dflux" are now given by   Delta F = sum [a9_j lam_j evec(,j)]

	  dflux(i9,1:5) = 0.
          do j9=1,5
	     dflux(i9,1) = dflux(i9,1) + a9(j9)*lam9(j9)*evec(1,j9)
	     dflux(i9,2) = dflux(i9,2) + a9(j9)*lam9(j9)*evec(2,j9)
	     dflux(i9,3) = dflux(i9,3) + a9(j9)*lam9(j9)*evec(3,j9)
	     dflux(i9,4) = dflux(i9,4) + a9(j9)*lam9(j9)*evec(4,j9)
	     dflux(i9,5) = dflux(i9,5) + a9(j9)*lam9(j9)*evec(5,j9)
	  enddo
	  dflux(i9,1) = dflux(i9,1)*Qlr(i9,rh)	 ! flux increment in density
	  dflux(i9,2) = dflux(i9,2)*Qlr(i9,rh)	 ! flux increment in energy    
	  dflux(i9,3) = dflux(i9,3)*Qlr(i9,rh)	 ! flux increment in parallel momentum
	  dflux(i9,4) = dflux(i9,4)*Qlr(i9,rh)	 ! flux increment in perpendicular momentum
	  dflux(i9,5) = dflux(i9,5)*Qlr(i9,rh)	 ! flux increment in z-momentum

	  if (kroe(i9) .gt. 0) then
            froef(i9,rh) = 0.5*(flr(i9,rh) + flr(i9+nface,rh) - dflux(i9,1))
            froef(i9,en) = 0.5*(flr(i9,en) + flr(i9+nface,en) - dflux(i9,2))
            froef(i9,mxa(ixyz)) = 0.5*(flr(i9,mxa(ixyz)) + flr(i9+nface,mxa(ixyz)) - dflux(i9,3))
            froef(i9,mya(ixyz)) = 0.5*(flr(i9,mya(ixyz)) + flr(i9+nface,mya(ixyz)) - dflux(i9,4))
            froef(i9,mza(ixyz)) = 0.5*(flr(i9,mza(ixyz)) + flr(i9+nface,mza(ixyz)) - dflux(i9,5))
	  end if

	end do 
  
    end subroutine flux_roe 

!----------------------------------------------------------------------------------------------

    subroutine innerintegral(Q_r)
    implicit none
    real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r
    integer i,j,k,ieq,ipg,ir
    real Qinner(npg,nQ),finner_x(npg,nQ), finner_y(npg,nQ), finner_z(npg,nQ), int_r(nbastot,nQ)
    real sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9

	integral_r(:,:,:,:,:) = 0.
	
	do k = 1,nz
	do j = 1,ny
	do i = 1,nx
	
	   do ieq = 1,nQ
	      do ipg = 1,npg
	         Qinner(ipg,ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_r(i,j,k,ieq,1:nbasis))
	      end do
	   end do

	   call flux_calc_pnts_r(Qinner,finner_x,1,npg)
	   call flux_calc_pnts_r(Qinner,finner_y,2,npg)
	   call flux_calc_pnts_r(Qinner,finner_z,3,npg)
		
	   do ieq = 1,nQ

!	      int_r(kx,ieq) = 0.25*cbasis(kx)*dxi*sum(wgt3d(1:npg)*finner_x(1:npg,ieq)) 
!	      int_r(ky,ieq) = 0.25*cbasis(ky)*dyi*sum(wgt3d(1:npg)*finner_y(1:npg,ieq)) 
!	      int_r(kz,ieq) = 0.25*cbasis(kz)*dzi*sum(wgt3d(1:npg)*finner_z(1:npg,ieq)) 

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


	   if (nbasis .gt. 4) then

	      if (ibitri .eq. 1 .or. iquad .gt. 2) then
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

	         int_r(kyz,ieq)  = 0.25*cbasis(kyz)*(dyi*sum1 + dzi*sum2)
	         int_r(kzx,ieq)  = 0.25*cbasis(kzx)*(dxi*sum3 + dzi*sum4)
	         int_r(kxy,ieq)  = 0.25*cbasis(kxy)*(dxi*sum5 + dyi*sum6)
	      end if
	      if (iquad .gt. 2) then
                 sum1 = 0.
                 sum2 = 0.
                 sum3 = 0.

                 do ipg=1,npg
                    sum1 = sum1 + wgt3d(ipg)*3.*bfvals_int(ipg,kx)*finner_x(ipg,ieq)
                    sum2 = sum2 + wgt3d(ipg)*3.*bfvals_int(ipg,ky)*finner_y(ipg,ieq)
                    sum3 = sum3 + wgt3d(ipg)*3.*bfvals_int(ipg,kz)*finner_z(ipg,ieq)
	         end do

	         int_r(kxx,ieq)  = 0.25*cbasis(kxx)*(dxi*sum1)
	         int_r(kyy,ieq)  = 0.25*cbasis(kyy)*(dyi*sum2)
	         int_r(kzz,ieq)  = 0.25*cbasis(kzz)*(dzi*sum3)
	      end if

	      if (ibitri .eq. 1 .or. iquad .gt. 3) then
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

	      if ((ibitri .eq. 1 .and. iquad .gt. 2) .or. iquad .gt. 3) then
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
	      end if

	      if (iquad .gt. 3) then
	         int_r(kxxx,ieq) = 0.25*cbasis(kxxx)*dxi*sum(wgt3d(1:npg)*(7.5*bfvals_int(1:npg,kx)**2 - 1.5)*finner_x(1:npg,ieq))
	         int_r(kyyy,ieq) = 0.25*cbasis(kyyy)*dyi*sum(wgt3d(1:npg)*(7.5*bfvals_int(1:npg,ky)**2 - 1.5)*finner_y(1:npg,ieq))
	         int_r(kzzz,ieq) = 0.25*cbasis(kzzz)*dzi*sum(wgt3d(1:npg)*(7.5*bfvals_int(1:npg,kz)**2 - 1.5)*finner_z(1:npg,ieq))
	      end if

	      if (ibitri .eq. 1 .and. iquad .gt. 2) then
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

	      end if   
           end do

	   do ieq = 1,nQ
	      do ir=1,nbasis
	         integral_r(i,j,k,ieq,ir) = int_r(ir,ieq) 
	      end do 
           end do
		
	end do
	end do
	end do
	
    end subroutine innerintegral

!----------------------------------------------------------------------------------------------

    subroutine glflux2
    implicit none
    real sumx,sumy,sumz
    integer i,j,k,ieq,ir,iqfa

	do ieq = 1,nQ
	do k = 1,nz
	do j = 1,ny

	do i = 1,nx 

            sumx = 0.
            sumy = 0.
            sumz = 0.

	    do iqfa = 1,nface
		sumx = sumx + 0.25*dxi*wgt2d(iqfa)*(flux_x(iqfa,i+1,j,k,ieq) - flux_x(iqfa,i,j,k,ieq))
		sumy = sumy + 0.25*dyi*wgt2d(iqfa)*(flux_y(iqfa,i,j+1,k,ieq) - flux_y(iqfa,i,j,k,ieq))
		sumz = sumz + 0.25*dzi*wgt2d(iqfa)*(flux_z(iqfa,i,j,k+1,ieq) - flux_z(iqfa,i,j,k,ieq))
	    end do

	    glflux_r(i,j,k,ieq,1) = sumx + sumy + sumz 

	end do

	do ir=2,nbasis
	do i = 1,nx

            sumx = 0.
            sumy = 0.
            sumz = 0.

	    do iqfa = 1,nface
		sumx = sumx + wgtbf_xmp(iqfa,1,ir)*flux_x(iqfa,i,j,k,ieq) + wgtbf_xmp(iqfa,2,ir)*flux_x(iqfa,i+1,j,k,ieq)
		sumy = sumy + wgtbf_ymp(iqfa,1,ir)*flux_y(iqfa,i,j,k,ieq) + wgtbf_ymp(iqfa,2,ir)*flux_y(iqfa,i,j+1,k,ieq)
		sumz = sumz + wgtbf_zmp(iqfa,1,ir)*flux_z(iqfa,i,j,k,ieq) + wgtbf_zmp(iqfa,2,ir)*flux_z(iqfa,i,j,k+1,ieq)
	    end do

	    glflux_r(i,j,k,ieq,ir) = sumx + sumy + sumz - integral_r(i,j,k,ieq,ir)

	end do  
	end do

	end do
	end do
	end do

    end subroutine glflux2

!----------------------------------------------------------------------------------------------

    subroutine glflux
    implicit none
    integer i,j,k,ieq,ir
	
	do ieq = 1,nQ
	do k = 1,nz
	do j = 1,ny
	do i = 1,nx


	   glflux_r(i,j,k,ieq,1) =  0.25*(dxi*sum(wgt2d(1:nface)*(flux_x(1:nface,i+1,j,k,ieq) - flux_x(1:nface,i,j,k,ieq))) &
				        + dyi*sum(wgt2d(1:nface)*(flux_y(1:nface,i,j+1,k,ieq) - flux_y(1:nface,i,j,k,ieq))) &
				        + dzi*sum(wgt2d(1:nface)*(flux_z(1:nface,i,j,k+1,ieq) - flux_z(1:nface,i,j,k,ieq))))
					  
	do ir=2,nbasis

	    glflux_r(i,j,k,ieq,ir) =  0.25*cbasis(ir)*  &
			(dxi*sum(wgt2d(1:nface)*(bfvals_xp(1:nface,ir)*flux_x(1:nface,i+1,j,k,ieq) - bfvals_xm(1:nface,ir)*flux_x(1:nface,i,j,k,ieq))) &
		       + dyi*sum(wgt2d(1:nface)*(bfvals_yp(1:nface,ir)*flux_y(1:nface,i,j+1,k,ieq) - bfvals_ym(1:nface,ir)*flux_y(1:nface,i,j,k,ieq))) & 
		       + dzi*sum(wgt2d(1:nface)*(bfvals_zp(1:nface,ir)*flux_z(1:nface,i,j,k+1,ieq) - bfvals_zm(1:nface,ir)*flux_z(1:nface,i,j,k,ieq)))) &
				      - integral_r(i,j,k,ieq,ir) 
	end do  

	end do
	end do
	end do
	end do

    end subroutine glflux


!-----------------------------------------------------------!
!*******************calculate freezing speeds***************!
!-----------------------------------------------------------!        
   real function cfcal(Qcf,cases)
        implicit none
        integer cases
        real Qcf(nQ)
     	real Pi, Pe, P, B2, ne
    	real dn,dni,vx,vy,vz,hx,hy,hz,dnei,va2,vf1,lil02,va,fac
      
        dn = Qcf(rh)
        dni = 1./dn
        vx = Qcf(mx)*dni
        vy = Qcf(my)*dni 
        vz = Qcf(mz)*dni
        if(ivis .ne. 0) P = (aindex - 1.)*(Qcf(en) - 0.5*dn*(vx**2 + vy**2 + vz**2))
        if(ivis .eq. 0) P = (Qcf(exx) + Qcf(eyy) + Qcf(ezz) - dn*(vx**2 + vy**2 + vz**2))/3.
	 
	select case (cases)

   	case (1) !freezing speed in x direction for fluid variable 
              cfcal = abs(vx) + sqrt(aindex*P*dni)

	case (2) !freezing speed in y direction for fluid variable
              cfcal = abs(vy) + sqrt(aindex*P*dni)

	case (3) !freezing speed in z direction for fluid variable
              cfcal = abs(vz) + sqrt(aindex*P*dni)

        end select
        
    end function cfcal 

!----------------------------------------------------------------------------------------------

    subroutine prepare_exchange(Q_r) 
	integer ieq, i, j, k, ipnt
        real, dimension(nx,ny,nz,nQ,nbasis) :: Q_r

	do ieq = 1,nQ
	do j = 1,ny
	do i = 1,nx
      do ipnt=1,nface	
		 Qzlow_int(i,j,ipnt,ieq) = sum(bfvals_zm(ipnt,1:nbasis)*Q_r(i,j,1,ieq,1:nbasis))
		 Qzhigh_int(i,j,ipnt,ieq) = sum(bfvals_zp(ipnt,1:nbasis)*Q_r(i,j,nz,ieq,1:nbasis))
	  end do
	end do
	end do
	end do
	
	do ieq = 1,nQ
	do k = 1,nz
	do i = 1,nx	
      do ipnt=1,nface	
		Qylow_int(i,k,ipnt,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_r(i,1,k,ieq,1:nbasis))
		Qyhigh_int(i,k,ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_r(i,ny,k,ieq,1:nbasis))
	  end do
	end do
	end do
	end do	
	
	do ieq = 1,nQ
	do k = 1,nz
	do j = 1,ny	
      do ipnt=1,nface	
		Qxlow_int(j,k,ipnt,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_r(1,j,k,ieq,1:nbasis))
		Qxhigh_int(j,k,ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_r(nx,j,k,ieq,1:nbasis))
	  end do		
	end do
	end do
	end do	
	call exchange_flux	

    end subroutine

!----------------------------------------------------------------------------------------------


    subroutine exchange_flux
        integer mpi_size 
    
	call MPI_BARRIER(cartcomm,ierr)

 	mpi_size=ny*nz*nface*nQ 

        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Qxhigh_int,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(1),ierr)
	endif

        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Qxlow_ext,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(2),ierr)
	endif
	
        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(1),stats(:,1),ierr)
           call MPI_IRecv(Qxhigh_ext,mpi_size,MPI_TT,nbrs(EAST),0,cartcomm,reqs(3),ierr)
	endif

        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(2),stats(:,2),ierr)
           call MPI_ISend(Qxlow_int,mpi_size,MPI_TT,nbrs(WEST),0,cartcomm,reqs(4),ierr)        			
	endif
		
        if (nbrs(EAST) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(3),stats(:,3),ierr)
	endif

        if (nbrs(WEST) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(4),stats(:,4),ierr)
	endif

	if (mpi_P .eq. 1 .and. xlbc .eq. 0) then
		Qxlow_ext = Qxlow_int
	end if	
	
	if (mpi_P .eq. mpi_nx .and. xhbc .eq. 0) then
		Qxhigh_ext = Qxhigh_int
	end if

        mpi_size=nface*nx*nz*nQ
       
        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
           call MPI_ISend(Qyhigh_int,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(1),ierr)
	endif

        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
           call MPI_IRecv(Qylow_ext,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(2),ierr)
	endif

        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(1),stats(:,1),ierr)
           call MPI_IRecv(Qyhigh_ext,mpi_size,MPI_TT,nbrs(NORTH),0,cartcomm,reqs(3),ierr)
	endif

        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(2),stats(:,2),ierr)
           call MPI_ISend(Qylow_int,mpi_size,MPI_TT,nbrs(SOUTH),0,cartcomm,reqs(4),ierr)        			
	endif
		
        if (nbrs(NORTH) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(3),stats(:,3),ierr)
	endif

        if (nbrs(SOUTH) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(4),stats(:,4),ierr)
	endif
	
		
	if (mpi_Q .eq. 1 .and. ylbc .eq. 0) then
		Qylow_ext = Qylow_int
	end if	
	
	if (mpi_Q .eq. mpi_ny .and. yhbc .eq. 0) then
		Qyhigh_ext = Qyhigh_int
	end if
	
        
        mpi_size=nface*nx*ny*nQ
       
        if (nbrs(UP) .ne. MPI_PROC_NULL) then
           call MPI_ISend(Qzhigh_int,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(1),ierr)
	endif
        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
           call MPI_IRecv(Qzlow_ext,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(2),ierr)
	endif
        if (nbrs(UP) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(1),stats(:,1),ierr)
           call MPI_IRecv(Qzhigh_ext,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(3),ierr)
	endif
        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(2),stats(:,2),ierr)
           call MPI_ISend(Qzlow_int,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(4),ierr)        			
	endif
        if (nbrs(UP) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(3),stats(:,3),ierr)
	endif

        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(4),stats(:,4),ierr)
	endif

		
	if (mpi_R .eq. 1 .and. zlbc .eq. 0) then
		Qzlow_ext = Qzlow_int
	end if	
	
	if (mpi_R .eq. mpi_nz .and. zhbc .eq. 0) then
		Qzhigh_ext = Qzhigh_int
	end if
		
    end subroutine
	

!-----------------------------------------------------------

        real function xc(i)
        integer i
        xc = loc_lxd + (i - 0.5)*dx
        end function xc

!-----------------------------------------------------------
        
        real function yc(j)
        integer j
        yc = loc_lyd + (j - 0.5)*dy
        end function yc

!-----------------------------------------------------------
        
        real function zc(k)
        integer k
        zc = loc_lzd + (k - 0.5)*dz
        end function zc

!-----------------------------------------------------------

        real function rz(i,j)
        integer i,j
        rz = sqrt(yc(j)**2)
        end function rz

!-----------------------------------------------------------

        real function r(i,j)
        integer i,j,k
        r = sqrt(xc(i)**2 + yc(j)**2)
        end function r

!-----------------------------------------------------------

        real function theta(i,j)
        integer i,j
        theta = atan2(yc(j),xc(i))
        end function theta

!-----------------------------------------------------------

        real function xvtk(i)
        integer i
        xvtk = loc_lxd + (i - 0.5)*dxvtk
        end function xvtk

!-----------------------------------------------------------
        
        real function yvtk(j)
        integer j
        yvtk = loc_lyd + (j - 0.5)*dyvtk
        end function yvtk

!-----------------------------------------------------------
        
        real function zvtk(k)
        integer k
        zvtk = loc_lzd + (k - 0.5)*dzvtk
        end function zvtk

!-----------------------------------------------------------

        real function rvtk(i,j)
        integer i,j,k
        rvtk = sqrt(xvtk(i)**2 + yvtk(j)**2)
        end function rvtk

!-----------------------------------------------------------

        real function thetavtk(i,j)
        integer i,j
        thetavtk = atan2(yvtk(j),xvtk(i))
        end function thetavtk


!-----------------------------------------------------------
       
        subroutine output_vtk0(Qin,nout,iam)
        implicit none
        real Qin(nx,ny,nz,nQ,nbasis)
	integer nout
	integer(I4P) , parameter :: nb=1*ngu,sprd=0*1
!	integer(I4P), parameter:: nnx=nx-2*nb,nny=ny-2*nb,nnz=nz-2*nb 
	real(R4P), dimension(nnx+1):: x_xml_rect
	real(R4P), dimension(nny+1):: y_xml_rect
	real(R4P), dimension(nnz+1):: z_xml_rect
	real(R4P), dimension(nnx*nny*nnz):: var_xml_val_x
	real(R4P), dimension(nnx*nny*nnz):: var_xml_val_y
	real(R4P), dimension(nnx*nny*nnz):: var_xml_val_z
	real P, vx, vy, vz,  dni
	integer(I4P):: E_IO,i,j,k,l,num,iam
        character (50) :: out_name
        character (4) :: tname
        character (5) :: tname1
        character (4) :: pname
        character (5) :: pname1
        
       
	num=nout+10000

        write(tname1,'(i5)')num
	tname=tname1(2:5) 
        tname = trim(tname)
        tname = adjustr(tname) 

	num=iam+10000
        
        write(pname1,'(i5)')num
	pname=pname1(2:5) 
        pname = trim(pname)
        pname = adjustr(pname) 
	out_name='/data/data8/perseus_p'//pname//'_t'//tname//'.vtr'
!	out_name='data/perseus_p'//pname//'_t'//tname//'.vtr'
!        print *, out_name
        out_name = trim(out_name)
        out_name = adjustr(out_name) 
			
		E_IO = VTK_INI_XML(output_format = 'BINARY',              &
		                   filename      = out_name, &
		                   mesh_topology = 'RectilinearGrid',     &
		                   nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1)

	
		do i=1+nb,nx-nb+1
			x_xml_rect(i)=(xc(i) - 0.5/dxi)*1e-3
		enddo
		do i=1+nb,ny-nb+1
			y_xml_rect(i)=(yc(i) - 0.5/dyi)*1e-3
		enddo
		do i=1+nb,nz-nb+1
			z_xml_rect(i)=(zc(i) - 0.5/dzi)*1e-3
		enddo


		E_IO = VTK_GEO_XML(nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1, &
			                 X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)
		
		E_IO = VTK_DAT_XML(var_location     = 'cell', &
		                   var_block_action = 'OPEN')

		do i=1+nb,nx-nb
			do j=1+nb,ny-nb
				do k=1+nb,nz-nb
					l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
					var_xml_val_x(l)=log(Qin(i,j,k,rh,1)*n0)/log(10.)

				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Log Density',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nx-nb
			do j=1+nb,ny-nb
				do k=1+nb,nz-nb
					l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
					var_xml_val_x(l)=Qin(i,j,k,rh,1)*n0
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Density',                    &
		                   var     = var_xml_val_x)
		                   
		do i=1+nb,nx-nb
			do j=1+nb,ny-nb
				do k=1+nb,nz-nb
					l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
					dni = 1./Qin(i,j,k,rh,1)					
					vx = Qin(i,j,k,mx,1)*dni
					vy = Qin(i,j,k,my,1)*dni
					vz = Qin(i,j,k,mz,1)*dni
					P = (aindex - 1)*(Qin(i,j,k,en,1) - 0.5*Qin(i,j,k,rh,1)*(vx**2 + vy**2 + vz**2))*dni
					var_xml_val_x(l)=P*te0
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Temperature',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nx-nb
			do j=1+nb,ny-nb
				do k=1+nb,nz-nb
					l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
					dni = 1./Qin(i,j,k,rh,1)					
					vx = Qin(i,j,k,mx,1)*dni
					vy = Qin(i,j,k,my,1)*dni
					vz = Qin(i,j,k,mz,1)*dni
					P = (aindex - 1)*(Qin(i,j,k,en,1) - 0.5*Qin(i,j,k,rh,1)*(vx**2 + vy**2 + vz**2))
					var_xml_val_x(l)=P
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'pressure',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nx-nb
			do j=1+nb,ny-nb
				do k=1+nb,nz-nb
					l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
					var_xml_val_x(l)=v0*Qin(i,j,k,mx,1)/Qin(i,j,k,rh,1)
					var_xml_val_y(l)=v0*Qin(i,j,k,my,1)/Qin(i,j,k,rh,1)
					var_xml_val_z(l)=v0*Qin(i,j,k,mz,1)/Qin(i,j,k,rh,1)
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Velocity',                    &
		                   varX     = var_xml_val_x, varY     = var_xml_val_y, varZ     = var_xml_val_Z )

		do i=1+nb,nx-nb
			do j=1+nb,ny-nb
				do k=1+nb,nz-nb
					l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
				        var_xml_val_x(l)= -2*(Q_r2(i,j,k,my,2)-Q_r2(i,j,k,mx,3))/Q_r2(i,j,k,rh,1)*dxi/t0
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Vorticity',                    &
		                   var     = var_xml_val_x)
		                   
		                   
		E_IO = VTK_DAT_XML(var_location     = 'cell', &
		                   var_block_action = 'Close')
		E_IO = VTK_GEO_XML()
		E_IO = VTK_END_XML()

        end subroutine output_vtk0

!-----------------------------------------------------------
       
        subroutine output_vtk(Qin,nout,iam)
        implicit none
        real Qin(nx,ny,nz,nQ,nbasis)
	integer nout
	integer(I4P) , parameter :: nb=1*ngu,sprd=0*1

	real(R4P), dimension(nnx+1) :: x_xml_rect
	real(R4P), dimension(nny+1) :: y_xml_rect
	real(R4P), dimension(nnz+1) :: z_xml_rect
	real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_x
	real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_y
	real(R4P), dimension(nnx*nny*nnz) :: var_xml_val_z
	real(R4P), dimension(nnx,nny,nnz,nQ) :: qvtk
	real(R4P), dimension(nnx,nny,nnz) :: qvtk_dxvy,qvtk_dyvx
	real P, vx, vy, vz,  dni, dxrh,dyrh,dxmy,dymx
	integer(I4P):: E_IO,i,j,k,l,num,iam,igrid,ir,jr,kr,ib,jb,kb,ieq
        character (50) :: out_name
        character (4) :: tname
        character (5) :: tname1
        character (4) :: pname
        character (5) :: pname1

       
	num=nout+10000

        write(tname1,'(i5)')num
	tname=tname1(2:5) 
        tname = trim(tname)
        tname = adjustr(tname) 

	num=iam+10000
        
        write(pname1,'(i5)')num
	pname=pname1(2:5) 
        pname = trim(pname)
        pname = adjustr(pname) 
	out_name='/data/data7/perseus_p'//pname//'_t'//tname//'.vtr'
!	out_name='data3/perseus_p'//pname//'_t'//tname//'.vtr'
!        print *, out_name
        out_name = trim(out_name)
        out_name = adjustr(out_name) 
			
		E_IO = VTK_INI_XML(output_format = 'BINARY',              &
		                   filename      = out_name, &
		                   mesh_topology = 'RectilinearGrid',     &
		                   nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1)



	
		do i=1+nb,nnx-nb+1
			x_xml_rect(i-nb)=(xvtk(i) - 0.5*dxvtk)*1e-3
		enddo
		do j=1+nb,nny-nb+1
			y_xml_rect(j-nb)=(yvtk(j) - 0.5*dyvtk)*1e-3
		enddo
		do k=1+nb,nnz-nb+1
			z_xml_rect(k-nb)=(zvtk(k) - 0.5*dzvtk)*1e-3
		enddo


		E_IO = VTK_GEO_XML(nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=1,nz2=nnz+1, &
			                 X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)
		
		E_IO = VTK_DAT_XML(var_location     = 'cell', &
		                   var_block_action = 'OPEN')


		do i=1+nb,nnx-nb
					ir = int(i/nvtk) + 1
				    ib = i - nvtk*int(i/nvtk)
					if (ib .eq. 0) then
						ib = nvtk
						ir = ir - 1
					end if
			do j=1+nb,nny-nb
					jr = int(j/nvtk) + 1
				    jb = j - nvtk*int(j/nvtk)
					if (jb .eq. 0) then
						jb = nvtk
						jr = jr - 1
					end if
				do k=1+nb,nnz-nb
					kr = int(k/nvtk) + 1
				    kb = k - nvtk*int(k/nvtk)
					if (kb .eq. 0) then
						kb = nvtk
						kr = kr - 1
					end if
					igrid = nvtk2*(ib-1) + nvtk*(jb-1) + kb
					do ieq=1,nQ
	      			   qvtk(i,j,k,ieq) = sum(bfvtk(igrid,1:nbasis)*Qin(ir,jr,kr,ieq,1:nbasis))
					end do

					   dxrh = sum(bfvtk_dx(igrid,1:nbasis)*Qin(ir,jr,kr,rh,1:nbasis))
					   dyrh = sum(bfvtk_dy(igrid,1:nbasis)*Qin(ir,jr,kr,rh,1:nbasis))
					   dxmy = sum(bfvtk_dx(igrid,1:nbasis)*Qin(ir,jr,kr,my,1:nbasis))
					   dymx = sum(bfvtk_dy(igrid,1:nbasis)*Qin(ir,jr,kr,mx,1:nbasis))
	      			   	qvtk_dxvy(i,j,k) = (qvtk(i,j,k,rh)*dxmy - qvtk(i,j,k,my)*dxrh)/qvtk(i,j,k,rh)**2
	      			   	qvtk_dyvx(i,j,k) = (qvtk(i,j,k,rh)*dymx - qvtk(i,j,k,mx)*dyrh)/qvtk(i,j,k,rh)**2	
				end do
			end do
		end do

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
				do k=1+nb,nnz-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
					var_xml_val_x(l)=log(qvtk(i,j,k,rh)*n0)/log(10.)
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Log Density',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
				do k=1+nb,nnz-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
					var_xml_val_x(l)=qvtk(i,j,k,rh)*n0
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Density',                    &
		                   var     = var_xml_val_x)	

	                   
		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
				do k=1+nb,nnz-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
					dni = 1./qvtk(i,j,k,rh)					
					vx = qvtk(i,j,k,mx)*dni
					vy = qvtk(i,j,k,my)*dni
					vz = qvtk(i,j,k,mz)*dni
					P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))*dni
					var_xml_val_x(l)=P*te0
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Temperature',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
				do k=1+nb,nnz-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
					dni = 1./qvtk(i,j,k,rh)					
					vx = qvtk(i,j,k,mx)*dni
					vy = qvtk(i,j,k,my)*dni
					vz = qvtk(i,j,k,mz)*dni
					P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))
					var_xml_val_x(l)=P
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'pressure',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
				do k=1+nb,nnz-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
					var_xml_val_x(l)=v0*qvtk(i,j,k,mx)/qvtk(i,j,k,rh)
					var_xml_val_y(l)=v0*qvtk(i,j,k,my)/qvtk(i,j,k,rh)
					var_xml_val_z(l)=v0*qvtk(i,j,k,mz)/qvtk(i,j,k,rh)
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Velocity',                    &
		                   varX     = var_xml_val_x, varY     = var_xml_val_y, varZ     = var_xml_val_Z )

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
				do k=1+nb,nnz-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
				        var_xml_val_x(l)= -2.*(qvtk_dxvy(i,j,k)-qvtk_dyvx(i,j,k))*dxi/t0
				enddo
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
		                   varname = 'Vorticity',                    &
		                   var     = var_xml_val_x)
		                   
		                   
		E_IO = VTK_DAT_XML(var_location     = 'cell', &
		                   var_block_action = 'Close')
		E_IO = VTK_GEO_XML()
		E_IO = VTK_END_XML()

        end subroutine output_vtk

!--------------------------------------------------------------------------------

	subroutine writeQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow)

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
!	print *,'fname2 ',fname2

	open(unit=3,file=fname2)

!	open(unit = 10, file = 'data/perseus_t'//dname//'_p'//pname//'.bin',form = 'unformatted',access = 'stream')

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

!------------------------------------------------------------------------------

	subroutine readQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow,mpi_nznow)

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

!--------------------------------------------------------------------------
	subroutine set_cbasis_3D

	! Using bi/tri elements:
		! iquad = 2, nbasis = 8: {1,x,y,z, yz, zx, xy, xyz}
		! iquad = 3, nbasis = 27: {1,x,y,z, yz, zx, xy, xyz, P2(x),P2(y),P2(z), yP2(z), zP2(x), xP2(y), P2(y)z, P2(z)x, P2(x)y,
		!    P2(y)P2(z), P2(z)P2(x), P2(x)P2(y), yzP2(x), zxP2(y), xyP2(z), xP2(y)P2(z), yP2(z)P2(x), zP2(x)P2(y), P2(x)P2(y)P2(z)}

	! Not using bi/tri elements:
		! iquad = 2, nbasis = 4: {1,x,y,z}
		! iquad = 3, nbasis = 10: {1,x,y,z,yz,zx,xy,P2(x),P2(y),P2(z)}
		! iquad = 4, nbasis = 20: {1,x,y,z,yz,zx,xy,P2(x),P2(y),P2(z),
!						xyz,yP2(z), zP2(x), xP2(y), P2(y)z, P2(z)x, P2(x)y,P3(x),P3(y),P3(z)}

	integer klim

	ibitri = 0
	if (iquad .eq. 2 .and. nbasis .eq. 4) ibitri = 0
	if (iquad .eq. 3 .and. nbasis .eq. 10) ibitri = 0
	if (iquad .eq. 4 .and. nbasis .eq. 20) ibitri = 0
	if (iquad .eq. 2 .and. nbasis .eq. 8) ibitri = 1
	if (iquad .eq. 3 .and. nbasis .eq. 27) ibitri = 1
	if (iquad .eq. 4 .and. nbasis .eq. 64) ibitri = 1

        kx = 2
	ky = 3
	kz = 4
	kyz = 5
	kzx = 6
	kxy = 7
	if (ibitri .eq. 0) klim = kxy
	if (ibitri .eq. 1) then
	   kxyz = 8
	   klim = kxyz
	end if

	kxx = klim + 1
	kyy = klim + 2
	kzz = klim + 3

	if (ibitri .eq. 0) then
	   kxyz = klim + 4
	   klim = kxyz
	end if
	if (ibitri .eq. 1) klim = kzz
	kyzz = klim + 1
	kzxx = klim + 2
	kxyy = klim + 3
    	kyyz = klim + 4
	kzzx = klim + 5
	kxxy = klim + 6
	if (ibitri .eq. 1) klim = kxxy
	if (ibitri .eq. 0) then
	   kxxx = klim + 7
	   kyyy = klim + 8
	   kzzz = klim + 9
	   klim = kzzz
	end if

!	If ibitri = 1.......
	kyyzz = klim + 1
	kzzxx = klim + 2
	kxxyy = klim + 3
	kyzxx = klim + 4
	kzxyy = klim + 5
	kxyzz = klim + 6
    	kxyyzz = klim + 7
	kyzzxx = klim + 8
	kzxxyy = klim + 9
	kxxyyzz = klim + 10
	if (ibitri .eq. 1) then
	   kxxx = klim + 11
	   kyyy = klim + 12
	   kzzz = klim + 13
	   klim = kzzz
	end if


	cbasis(1) = 1.			! coefficient for basis function {1}
	cbasis(kx:kz) = 3.		! coefficients for basis functions {x,y,z}
	cbasis(kyz:kxy) = 9.		! coefficients for basis functions {yz,zx,xy}
	cbasis(kxyz) = 27. 		! coefficient for basis function {xyz}

	cbasis(kxx:kzz) = 5.		! coefficients for basis functions {P2(x),P2(y),P2(z)}
	cbasis(kyzz:kxyy) = 15.		! coefficients for basis functions {yP2(z),zP2(x),xP2(y)}
	cbasis(kyyz:kxxy) = 15.		! coefficients for basis functions {P2(y)z,P2(z)y,P2(z)x}
	cbasis(kyyzz:kxxyy) = 25.	! coefficients for basis functions {P2(y)P2(z),P2(z)P2(x),P2(x)P2(y)}
	cbasis(kyzxx:kxyzz) = 45.	! coefficients for basis functions {yzP_2(x),zxP_2(y),xyP_2(z)}
	cbasis(kxyyzz:kzxxyy) = 75.	! coefficients for basis functions {xP2(y)P2(z),yP2(z)P2(x),zP2(x)P2(y)}
	cbasis(kxxyyzz) = 125.		! coefficients for basis functions {P2(x)P2(y)P2(z)}
	cbasis(kxxx:kzzz) = 7.		! coefficients for basis functions {P3(x),P3(y),P3(z)}
    
	end subroutine

!--------------------------------------------------------------------------------
	
	subroutine set_bfvals_3D
		! Defines local basis function values and weights for 1, 2, or 3-point Gaussian quadrature.
		! Basis functions are evaluated in cell interior and on cell faces.

	implicit none
	
        call set_vtk_vals_3D()    ! Define basis function values on a 3D grid INTERNAL to the cell.
        call set_internal_vals_3D()     ! Define basis function values at quadrature points INTERNAL to cell. 
        call set_face_vals_3D()   ! Define local basis function values at quadrature points on a cell face.   	
        call set_weights_3D()     ! Define weights for integral approximation using Gaussian quadrature.

	if (iam .eq. print_mpi) then
	   print *,'testing basis...'
	   call test_basis_3D()
	   print *,'done testing basis...'
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


	    bfvtk(1,1) = 1.		! basis function = 1
	    bfvtk(1,kx) = 0.		! basis function = x		
	    bfvtk(1,ky) = 0.		! basis function = y	
	    bfvtk(1,kz) = 0.		! basis function = z	
	
	
		bfvtk(1:nvtk3,1) = 1.		! basis function = 1
		do i=1,nvtk
	           bfvtk((i-1)*nvtk+1:i*nvtk,ky) = xgrid(i)		! basis function = y	   		
	           bfvtk((i-1)*nvtk+1:i*nvtk,kz) = xgrid(1:nvtk)		! basis function = z
		end do
		do i=1,nvtk
		   bfvtk((i-1)*nvtk2+1:i*nvtk2,kx) = xgrid(i)		! basis function = x
		   bfvtk((i-1)*nvtk2+1:i*nvtk2,ky) = bfvtk(1:nvtk2,ky)		! basis function = y
		   bfvtk((i-1)*nvtk2+1:i*nvtk2,kz) = bfvtk(1:nvtk2,kz)		! basis function = z
		end do		

	do i=0,2		
		bfvtk(1:nvtk3,kxx+i) = 1.5*bfvtk(1:nvtk3,kx+i)**2 - 0.5    ! basis function = P_2(s)
	!	bfvtk(1:nvtk3,kxxx+i) = 2.5*bfvtk(1:nvtk3,kx+i)**3 - 1.5*bfvtk(1:nvtk3,kx+i)   ! basis function = P3(s)
	end do
	
	bfvtk(1:nvtk3,kyz) = bfvtk(1:nvtk3,ky)*bfvtk(1:nvtk3,kz)		! basis function = yz
	bfvtk(1:nvtk3,kzx) = bfvtk(1:nvtk3,kz)*bfvtk(1:nvtk3,kx)		! basis function = zx
	bfvtk(1:nvtk3,kxy) = bfvtk(1:nvtk3,kx)*bfvtk(1:nvtk3,ky)		! basis function = xy
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


		!   bfvtk_dx(igrid,kxxx) = 7.5*bfvtk(igrid,kx)**2 - 1.5
		!   bfvtk_dx(igrid,kyyy:kzzz) = 0.
		!   bfvtk_dy(igrid,kyyy) = 7.5*bfvtk(igrid,ky)**2 - 1.5
		!   bfvtk_dy(igrid,kxxx) = 0.
		!   bfvtk_dy(igrid,kzzz) = 0.
		!   bfvtk_dz(igrid,kzzz) = 7.5*bfvtk(igrid,kz)**2 - 1.5
		!   bfvtk_dz(igrid,kxxx) = 0.
		!   bfvtk_dz(igrid,kyyy) = 0.			   
	   end do


	end subroutine set_vtk_vals_3D

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

	if (iquad .eq. 1) then			! 2-point Gaussian quadrature
	    bfvals_int(1,1) = 1.		! basis function = 1
	    bfvals_int(1,kx) = 0.		! basis function = x		
	    bfvals_int(1,ky) = 0.		! basis function = y	
	    bfvals_int(1,kz) = 0.		! basis function = z	
	end if	
	

	if (iquad .gt. 1) then
		bfvals_int(1:npg,1) = 1.		! basis function = 1
		do i=1,nedge
	       bfvals_int((i-1)*nedge+1:i*nedge,ky) = xquad(i)		! basis function = y	   		
	       bfvals_int((i-1)*nedge+1:i*nedge,kz) = xquad(1:nedge)		! basis function = z
		end do
		do i=1,nedge
		   bfvals_int((i-1)*nface+1:i*nface,kx) = xquad(i)		! basis function = x
		   bfvals_int((i-1)*nface+1:i*nface,ky) = bfvals_int(1:nface,ky)		! basis function = y
		   bfvals_int((i-1)*nface+1:i*nface,kz) = bfvals_int(1:nface,kz)		! basis function = z
		end do		
	end if

	do i=0,2		
		bfvals_int(1:npg,kxx+i) = 1.5*bfvals_int(1:npg,kx+i)**2 - 0.5    ! basis function = P_2(s)
	!	bfvals_int(1:npg,kxxx+i) = 2.5*bfvals_int(1:npg,kx+i)**3 - 1.5*bfvals_int(1:npg,kx+i)   ! basis function = P_3(s)
	end do
	
	bfvals_int(1:npg,kyz) = bfvals_int(1:npg,ky)*bfvals_int(1:npg,kz)		! basis function = yz
	bfvals_int(1:npg,kzx) = bfvals_int(1:npg,kz)*bfvals_int(1:npg,kx)		! basis function = zx
	bfvals_int(1:npg,kxy) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,ky)		! basis function = xy
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
	bfvals_int(1:npg,kxxyyzz) = bfvals_int(1:npg,kxx)*bfvals_int(1:npg,kyy)*bfvals_int(1:npg,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

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
	!	bfvals_xp(1:nface,kxxx+i-1) = 2.5*bfvals_xp(1:nface,kx+i-1)**3 - 1.5*bfvals_xp(1:nface,kx+i-1)
	!	bfvals_yp(1:nface,kxxx+i-1) = 2.5*bfvals_yp(1:nface,kx+i-1)**3 - 1.5*bfvals_yp(1:nface,kx+i-1)
	!	bfvals_zp(1:nface,kxxx+i-1) = 2.5*bfvals_zp(1:nface,kx+i-1)**3 - 1.5*bfvals_zp(1:nface,kx+i-1)
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
	bfvals_xp(1:nface,kxxyyzz) = bfvals_xp(1:nface,kxx)*bfvals_xp(1:nface,kyy)*bfvals_xp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)
	
	bfvals_xm = bfvals_xp
	bfvals_xm(1:nface,kx) = -bfvals_xp(1:nface,kx)
	bfvals_xm(1:nface,kzx) = -bfvals_xp(1:nface,kzx)
	bfvals_xm(1:nface,kxy) = -bfvals_xp(1:nface,kxy)
	bfvals_xm(1:nface,kxyz) = -bfvals_xp(1:nface,kxyz)
	bfvals_xm(1:nface,kzzx) = -bfvals_xp(1:nface,kzzx)
	bfvals_xm(1:nface,kxyy) = -bfvals_xp(1:nface,kxyy)
	!bfvals_xm(1:nface,kxxx) = -bfvals_xp(1:nface,kxxx)	
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
	bfvals_yp(1:nface,kxxyyzz) = bfvals_yp(1:nface,kxx)*bfvals_yp(1:nface,kyy)*bfvals_yp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)
	
	bfvals_ym = bfvals_yp
	bfvals_ym(1:nface,ky) = -bfvals_yp(1:nface,ky)
	bfvals_ym(1:nface,kyz) = -bfvals_yp(1:nface,kyz)
	bfvals_ym(1:nface,kxy) = -bfvals_yp(1:nface,kxy)
	bfvals_ym(1:nface,kxyz) = -bfvals_yp(1:nface,kxyz)
	bfvals_ym(1:nface,kyzz) = -bfvals_yp(1:nface,kyzz)
	bfvals_ym(1:nface,kxxy) = -bfvals_yp(1:nface,kxxy)
	!bfvals_ym(1:nface,kyyy) = -bfvals_yp(1:nface,kyyy)
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
	bfvals_zp(1:nface,kxxyyzz) = bfvals_zp(1:nface,kxx)*bfvals_zp(1:nface,kyy)*bfvals_zp(1:nface,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

	
	bfvals_zm = bfvals_zp
	bfvals_zm(1:nface,kz) = -bfvals_zp(1:nface,kz)
	bfvals_zm(1:nface,kyz) = -bfvals_zp(1:nface,kyz)
	bfvals_zm(1:nface,kzx) = -bfvals_zp(1:nface,kzx)
	bfvals_zm(1:nface,kxyz) = -bfvals_zp(1:nface,kxyz)
	bfvals_zm(1:nface,kyyz) = -bfvals_zp(1:nface,kyyz)	
	bfvals_zm(1:nface,kzxx) = -bfvals_zp(1:nface,kzxx)
	!bfvals_zm(1:nface,kzzz) = -bfvals_zp(1:nface,kzzz)
	bfvals_zm(1:nface,kyzxx) = -bfvals_zp(1:nface,kyzxx)
	bfvals_zm(1:nface,kzxyy) = -bfvals_zp(1:nface,kzxyy)
	bfvals_zm(1:nface,kzxxyy) = -bfvals_zp(1:nface,kzxxyy)
	
	! Organize local basis values on faces into 1-D vectors.
	! Used in limiter() and max_lim().

	bf_faces(1:nslim,1) = 1.		! basis function = 1
	
	do ixyz=kx,kz
		bf_faces(1:nface,ixyz) = bfvals_xm(1:nface,ixyz)		! basis function = x,y,z
		bf_faces(nface+1:2*nface,ixyz) = bfvals_xp(1:nface,ixyz)		! basis function = x,y,z
		bf_faces(2*nface+1:3*nface,ixyz) = bfvals_ym(1:nface,ixyz)		! basis function = x,y,z
		bf_faces(3*nface+1:4*nface,ixyz) = bfvals_yp(1:nface,ixyz)		! basis function = x,y,z
		bf_faces(4*nface+1:5*nface,ixyz) = bfvals_zm(1:nface,ixyz)		! basis function = x,y,z
		bf_faces(5*nface+1:6*nface,ixyz) = bfvals_zp(1:nface,ixyz)		! basis function = x,y,z
		bf_faces(6*nface+1:nslim,ixyz) = bfvals_int(1:npg,ixyz)		! basis function = x,y,z
	end do

	bf_faces(1:nslim,kyz) = bf_faces(1:nslim,ky)*bf_faces(1:nslim,kz)     ! basis function = yz
	bf_faces(1:nslim,kzx) = bf_faces(1:nslim,kz)*bf_faces(1:nslim,kx)     ! basis function = zx
	bf_faces(1:nslim,kxy) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,ky)     ! basis function = xy
	bf_faces(1:nslim,kxyz) = bf_faces(1:nslim,kx)*bf_faces(1:nslim,ky)*bf_faces(1:nslim,kz)     ! basis function = xyz
	
	do i=0,2		
		bf_faces(1:nslim,kxx+i) = 1.5*bf_faces(1:nslim,kx+i)**2 - 0.5    ! basis function = P_2(s)
	!	bf_faces(1:nslim,kxxx+i) = 2.5*bf_faces(1:nslim,kx+i)**3 - 1.5*bf_faces(1:nslim,kx+i)   ! basis function = P_3(s)
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
	bf_faces(1:nslim,kxxyyzz) = bf_faces(1:nslim,kxx)*bf_faces(1:nslim,kyy)*bf_faces(1:nslim,kzz)  ! basis function = P_2(x)P_2(y)P_2(z)

	end subroutine set_face_vals_3D

!----------------------------------------------------------
	subroutine set_weights_3D

	integer i
	real wq4p,wq4m

	wq4p = (18. + sqrt(30.))/36.
	wq4m = (18. - sqrt(30.))/36.

	! Define weights for integral approximation using Gaussian quadrature.

	!  Define weights for 1-D integration

	 if (iquad .eq. 1) then		! 1-point quadrature
	    wgt1d(1) = 2.
	 end if
	 if (iquad .eq. 2) then		! 2-point quadrature
	    wgt1d(1:2) = 1.
	 end if
	 if (iquad .eq. 3) then		! 3-point quadrature
	    wgt1d(1) = 5./9.
	    wgt1d(2) = 8./9.
	    wgt1d(3) = 5./9.
	 end if
	 if (iquad .eq. 4) then		! 4-point quadrature
	    wgt1d(1) = wq4m
	    wgt1d(2) = wq4p
	    wgt1d(3) = wq4p
	    wgt1d(4) = wq4m
	 end if
	
	!  Define weights for 2-D integration
	
	 if (iquad .eq. 1) then		! 1-point quadrature
	    wgt2d(1) = 4.
	 end if
	 if (iquad .eq. 2) then		! 2-point quadrature
	    wgt2d(1:4) = 1.
	 end if

	if (iquad .ge. 3) then
	   do i= 1,nedge
		  wgt2d((i-1)*nedge+1:i*nedge) = wgt1d(1:nedge)*wgt1d(i)
	   end do
	end if
	
	!  Define weights for 3-D integration

	 if (iquad .eq. 1) then		! 1-point quadrature
	    wgt3d(1) = 8.
	 end if
	 if (iquad .eq. 2) then		! 2-point quadrature
	    wgt3d(1:8) = 1.
	 end if

	if (iquad .ge. 3) then
	   do i= 1,nedge
		  wgt3d((i-1)*nface+1:i*nface) = wgt2d(1:nface)*wgt1d(i)
	   end do
	end if


	end subroutine set_weights_3D

!-----------------------------------------------------------------------------

	subroutine test_basis_3D

!	This routine tests whether the basis functions are orthonormal over the cell volume and over each cell face. 

	implicit none
	integer ir,jr,ix1,iy1,iz1,ipass
	real sq3,sq5,sq7,epsm

	sq3 = sqrt(3.)
	sq5 = sqrt(5.)
	sq7 = sqrt(7.)


!	Indices of basis elements projected onto the x-faces (yz-plane):

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

!	Indices of basis elements projected onto the y-faces (zx-plane):

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

!	Indices of basis elements projected onto the z-faces (xy-plane):

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

!	Normalization of basis elements over the  positive x-face of a cell:

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

!	Normalization of basis elements over the  positive y-face of a cell:

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

!	Normalization of basis elements over the  positive z-face of a cell:

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

!	Normalization of basis elements over the  negative x-face of a cell:

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

!	Normalization of basis elements over the  negative y-face of a cell:

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

!	Normalization of basis elements over the  negative z-face of a cell:

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

!	Below, we specify what the scalar product of any two basis elements should evaluate to.

!	cell_int0: scalar product over cell volume
!	xp_int0: scalar product over positive x-face
!	xm_int0: scalar product over negative x-face
!	yp_int0: scalar product over positive y-face
!	ym_int0: scalar product over negative y-face
!	zp_int0: scalar product over positive z-face
!	zm_int0: scalar product over negative z-face

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
		

!	Below, we evaluate the scalar product of all pairs of basis elements and test whether it's equal to the 
!	expected value computed above.

!	cell_int: scalar product over cell volume
!	xp_int: scalar product over positive x-face
!	xm_int: scalar product over negative x-face
!	yp_int: scalar product over positive y-face
!	ym_int: scalar product over negative y-face
!	zp_int: scalar product over positive z-face
!	zm_int: scalar product over negative z-face

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

!	if (ipass .eq. 0) call exit(-1)

	end subroutine test_basis_3D

!--------------------------------------------------------------------------------

end program


