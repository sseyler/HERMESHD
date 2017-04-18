!      DG-PERSEUS on 3 dimension: Extension of the finite volume PERSEUS Algorithm (Physics of the Extended-mhd Relaxation System using an Efficient Upwind Scheme, by Matt Martin and Charles Seyler) to the Discontinuous Galerkin Method.
!                
!                DG Algorithm by Xuan Zhao, Charles Seyler, and Yang Yang. 
!
!                Solves the two-fluid equations in Generalized Ohm's Law (GOL) form in 1D cartesian geometry
!                 
!                Variables:    There are 14 dependent variables: density (rh), velocity (vx,vy,vz), energy (en), magnetic field (bx,by,bz), electric field (ex,ey,ez), current density (jx,jy,jz)
!
!                Units:        A value of unity for each variable or parameter corresponds to the following dimensional units:
!
!                              Length              1 mm
!                              Time                100 ns
!                              Number density      6 x 10^22 cm^-3  
!                              Velocity            10^4 ms^-1
!                              Magnetic field      580 T
!                              Electric field      5.8 x 10^6 Vm^-1
!                              Current density     4.6 x 10^11 Am^-2
!                              Temperature         14 eV
!                              Resistivity         1.3 x 10^-5 Ohm-m
!    
!                Features: 
!                              Written in Fortran 90 in highly structured form
!                              GOL is treated as an evolution equation 
!                              Electric field is advanced using displacement current in Ampere-Maxwell law
!                              E and J are treated implicitly as source terms (i.e. no spatial derivatives an only a 3x3 solve)
!                              Second order discontinuous Galerkin method using positivity-preserving limiting for spatial discretization
!                              Div B = 0 preserved by Local-Structure-preserving DG method
!                              Freezing speed - relaxation method used instead of a Riemann solver (more diffusive but faster and simpler)
!                              No guard cells are needed for Boundary Condition  
!
!                Capabilities:
!                              Can handle about 8 orders of magnitude in density variation in single precision
!                              Can study Hall phenomena on long time scales
!                              No difficulties with "vacuum currents" - handles density floor very well
!                              Scales well with increasing number of processors
!                              Parallel electric field physics is much better than resistive MHD
!
!                Limitations:               
!                              Numerical speed of light cannot be too low otherwise displacement current dominates J (er2 <= 10 is ok)
!-------------------------------------------------------------------------------------------------------------------------------------------                                        

    program main  
    use lib_vtk_io 
    implicit none
    include 'mpif.h'  
!    implicit none
	integer, parameter :: rh=1, mx=2, my=3, mz=4, en=5, exx=6, eyy=7, ezz=8, exy=9, exz=10, eyz=11, nQ=11

	integer, parameter :: iquad=4
	    ! iquad: number of Gaussian quadrature points per direction.

	integer, parameter :: nx = 40, ny = 40, ngu = 0, nbasis = 16, nbastot = 16
        integer nbasv(nQ)	

    	integer, parameter :: nedge = iquad, npg = nedge*iquad, nfe = 2*nedge, nslim = npg+4*nedge
		! nedge: number of quadrature points per cell edge.
		! npg: number of internal points per cell.
		! nslim: number of quadrature points on cell edges and interior; for use in limiters.

	integer, parameter :: xlbc = 2, xhbc = 2, ylbc = 2, yhbc = 2

	integer, parameter :: ntout = 100, iorder = 2

	integer, parameter :: ihllc = 1, iroe = 0, ivis = 0
	! Choose Riemann solver for computing fluxes.  Set chosen solver = 1.
	! If all of them = 0, then LLF is used for fluxes.
        ! LLF is very diffusive for the hydro problem.  Roe and HLLC are much less diffusive than LLF and give very similar results with similar cpu overhead.

	! to restart from a checkpoint, set iread to 1 or 2 (when using the odd/even scheme)
	integer, parameter :: iread = 0, iwrite = 0
	character (4), parameter :: fpre = 'Qout'
	logical, parameter :: resuming = .false.

	real, parameter :: lx = 300., ly = 300.

	integer :: mpi_nx=6, print_mpi=0 

	real, parameter :: tf = 2000.
                 
        real, parameter :: pi = 4.0*atan(1.0)

        real, parameter :: aindex = 5./3., mu = 22.

	real, parameter :: aindm1 = aindex - 1.0, cp = aindex/(aindex - 1.0), clt = 2.

    	! dimensional units (expressed in MKS)
        real, parameter :: L0=1.0e4, t0=1.0e2, n0=6.e18

    	! derived units 
        real, parameter :: v0=L0/t0, p0=mu*1.67e-27*n0*v0**2, te0=p0/n0/1.6e-19, vis=1.e-2

	real, parameter :: T_floor = 0.01/te0, rh_floor = 1.e-5, P_floor = T_floor*rh_floor, rh_mult = 1.01

	real, dimension(nx,ny,nQ,nbasis) ::  Q_r0, Q_r1, Q_r2, Q_r3	
	real, dimension(nx,ny,nQ,nbasis) ::  glflux_r, source_r
	real, dimension(nx,ny,nQ,nbasis) :: integral_r 

	real cflm
	real flux_x(nedge,1:nx+1,ny,1:nQ), flux_y(nedge,nx,1:ny+1,1:nQ)
        real cfrx(nedge,nQ),cfry(nedge,nQ)
	logical MMask(nx,ny),BMask(nx,ny)   
	integer ticks, count_rate, count_max
	real t1, t2, t3, t4, elapsed_time, t_start, t_stop, dtoriginal
        real t,dt,tout,dtout,tlast,vf,dxi,dyi,loc_lxd,loc_lyd,check_Iz,sl, dy, dx
        real lxd,lxu,lyd,lyu,pin_rad,pin_height,foil_thickness,rh_foil,rh_solid,pin_rad_in,pin_rad_out,rim_rad,wire_rad
        real disk_rad,disk_height,foil_rad,buf_rad,dish_height,foil_height,T_init
        real gpz_rad,rh_gpz,kappa,rh_fluid,coll
	real Qxhigh_ext(ny,nedge,nQ), Qxlow_int(ny,nedge,nQ), Qxlow_ext(ny,nedge,nQ), Qxhigh_int(ny,nedge,nQ) 
	real Qyhigh_ext(nx,nedge,nQ), Qylow_int(nx,nedge,nQ), Qylow_ext(nx,nedge,nQ), Qyhigh_int(nx,nedge,nQ)
	
	integer mxa(3),mya(3),mza(3),kroe(nedge),niter,nlim_calls

    ! Parameters for liquid-vapor EOS

	integer, parameter :: nte0=50
	integer nte
			
    ! MPI definitions 
       
    	integer iam,ierr,mpi_ny,numprocs,reorder,cartcomm,mpi_P,mpi_Q
    	integer dims(2),coords(2),periods(2),nbrs(4),reqs(4),stats(MPI_STATUS_SIZE,4)
    	integer,parameter:: UP=1,DOWN=2,LEFT=3,RIGHT=4,MPI_TT=MPI_REAL4

	! Parameters relating to quadratures and basis functions.

	real wgt1d(5), wgt2d(30), wgt3d(100), cbasis(nbastot)
		! wgt1d: quadrature weights for 1-D integration
		! wgt2d: quadrature weights for 2-D integration
		! wgt3d: quadrature weights for 3-D integration	

	
	real, dimension(nedge,nbastot) :: bfvals_yp, bfvals_ym, bfvals_xp, bfvals_xm
	real bf_edges(nslim,nbastot), bfvals_int(npg,nbastot)
	real bfint_dx(npg,nbastot), bfint_dy(npg,nbastot)

    	integer ibitri,kx,ky,kxx,kyy,kxy,kxxx,kyyy,kxxy,kxyy,kxxyy

	! Parameters for VTK output.

	integer, parameter :: nvtk=1
	integer, parameter :: nvtk2=nvtk*nvtk
	integer(I4P), parameter :: nnx=nx*nvtk, nny=ny*nvtk, nnz=1
	real, dimension(nvtk2,nbastot) :: bfvtk, bfvtk_dx, bfvtk_dy
	real xgrid(20),dxvtk,dyvtk
	integer nout, ioe,mpi_nx_p,mpi_ny_p,nout_p,e_betainv
        real t_p,dt_p,dtout_p


	ibitri = 0
	if (iquad .eq. 2 .and. nbasis .eq. 3) ibitri = 0
	if (iquad .eq. 3 .and. nbasis .eq. 6) ibitri = 0
	if (iquad .eq. 4 .and. nbasis .eq. 10) ibitri = 0

	if (iquad .eq. 2 .and. nbasis .eq. 4) ibitri = 1
	if (iquad .eq. 3 .and. nbasis .eq. 9) ibitri = 1
	if (iquad .eq. 4 .and. nbasis .eq. 16) ibitri = 1

	kx = 2
	ky = 3
	kxy = 4
	kxx = 5
	kyy = 6
	kxxy = 7
	kxyy = 8

	if (ibitri .eq. 1) then
	   kxxyy = 9
	   kxxx = 10
	   kyyy = 11
	end if	
	if (ibitri .eq. 0) then
	   kxxx = 9
	   kyyy = 10
	   kxxyy = 11
	end if

	cbasis(1) = 1.			! coefficient for basis function {1}
	cbasis(kx:ky) = 3.		! coefficients for basis functions {x,y}
	cbasis(kxy) = 9.		! coefficients for basis functions {yz,zx,xy}
	cbasis(kxx:kyy) = 5.		! coefficients for basis functions {P_2(x),P_2(y)}
	cbasis(kxxy:kxyy) = 15.		! coefficients for basis functions {P2(x)y,P2(y)x}
	cbasis(kxxyy) = 25.		! coefficient for basis functions P2(x)P2(y)
	cbasis(kxxx:kyyy) = 7.		! coefficients for basis functions {P3(x),P3(y)}

	if (ibitri .eq. 0) then
	   if (iquad .eq. 2) cflm = 0.14
	   if (iquad .eq. 3) cflm = 0.09 
	   if (iquad .eq. 4) cflm = 0.05
	end if
	if (ibitri .eq. 1) then
	   if (iquad .eq. 2) cflm = 0.14
	   if (iquad .eq. 3) cflm = 0.07 
	   if (iquad .eq. 4) cflm = 0.05
	end if	! coefficients for basis functions {P2(x)P2(y)P2(z)}

   
    ! Initialize grid sizes and local lengths
    
	call MPI_Init ( ierr )
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	
	mpi_ny=numprocs/mpi_nx

	dims(1)=mpi_nx
	dims(2)=mpi_ny
	
	periods(:)=0
	if (xhbc .eq. 2) then
     	   periods(1)=1
	end if
	if (yhbc .eq. 2) then
     	   periods(2)=1
	end if
	
	reorder = 1
	
	call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, periods, reorder,cartcomm, ierr)
    	call MPI_COMM_RANK (cartcomm, iam, ierr )
    	call MPI_CART_COORDS(cartcomm, iam, 2, coords, ierr)
    	mpi_P=coords(1)+1
    	mpi_Q=coords(2)+1
    	call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(LEFT), nbrs(RIGHT), ierr)
    	call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(DOWN), nbrs(UP), ierr)	

    	lyd = -lx/2.
        lxd = -ly/2.

    	lxu = lxd + lx
    	lyu = lyd + ly 
    	
    	dxi = (nx*mpi_nx)/lx
    	dyi = (ny*mpi_ny)/ly

	dx = 1./dxi
	dy = 1./dyi

    	loc_lxd = lxd + (mpi_P-1)*lx/mpi_nx
    	loc_lyd = lyd + (mpi_Q-1)*ly/mpi_ny

        rh_fluid = 1.

	nbasv(:) = nbasis

	mxa(1) = mx
	mxa(2) = my
	mxa(3) = mz
	mya(1) = my
	mya(2) = mx
	mya(3) = mz
	mza(1) = mz
	mza(2) = mz
	mza(3) = mz

	nlim_calls = 0
	t = 0.
	dt = cflm*dx/clt
	dtoriginal = dt
	nout = 0
	niter = 0
	dtout = tf/ntout

	! Evaluate local cell values of basis functions on cell interior and edges.
	! Also defines cell values for sub-sampled output_vtk grid.
	! This is done for 1, 2, 3 or 4 point Gaussian quadrature.
	call set_bfvals
    	
	call system_clock(ticks, count_rate, count_max)
	t_start = ticks*1./count_rate
	
	if (iread .eq. 0) then
		call initial_condition

	else

		! This applies only if the initial data are being read from an input file.
		! - If resuming a run, keep the previous clock (i.e., t at nout) running.
		! - If not resuming a run, treat input as initial conditions at t=0, nout=0.

		call readQ(fpre,iam,iread,Q_r0,t_p,dt_p,nout_p,mpi_nx_p,mpi_ny_p)

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
		end if
 		! Quit if dtout is incompatible with input t/(nout-1)
		if (abs(dtout_p-dtout)/dt_p > 1.01) then
			if (iam .eq. print_mpi) then
				print *, 'Bad restart, non-matching dtout'
			end if
			call exit(-1)
		end if
		if ((mpi_nx_p .ne. mpi_nx) .or. (mpi_ny_p .ne. mpi_ny)) then
			if (iam .eq. print_mpi) then
				print *, 'Bad restart, non-matching mpi_nx or mpi_ny'
			end if
			call exit(-1)
		end if

	end if

        if (iam .eq. print_mpi) then
    	print *,'total dim= ',mpi_nx*nx,mpi_ny*ny
    	print *,'mpi dim= ',mpi_nx,mpi_ny
      	print *, 'v0 is: ', v0
	print *, 'te0 is: ', te0
	print *, 'dx is: ', ly/(ny*mpi_ny)*L0
	print *, 'coll is: ', coll
      	print *, 'aindex is: ', aindex
        end if

	tlast = t

	call system_clock( ticks, count_rate, count_max )
    	t1 = 1.*ticks / count_rate
	call output_vtk(Q_r0,nout,iam)	
    	
	do while (t<tf)

!	   if (mod(niter,2000) .eq. 0 .and. iam .eq. print_mpi) print *,'niter,t,dt = ',niter,t,dt,dtout*nout
	   niter = niter + 1
 	call get_min_dt(dt)
		
	if (iorder .eq. 2) then
		
		call prep_advance(Q_r0)
		call advance_time_level_gl(Q_r0,Q_r1)
                if(ivis .eq. 0)call implicit_source(Q_r1,Q_r1)

		call prep_advance(Q_r1)
		call advance_time_level_gl(Q_r1,Q_r2)
                if(ivis .eq. 0)call implicit_source(Q_r2,Q_r2)

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

	if (t-tlast .gt. dtout) then

		tlast = t
		nout = nout + 1

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
	!		call writeQ(fpre,iam,ioe,Q_r0,t,dt,nout,mpi_nx,mpi_ny)
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

!------------------------------------------------------------------------------- 

       subroutine prep_advance(Q_ri)

	    real, dimension(nx,ny,nQ,nbasis) :: Q_ri

		call prepare_exchange(Q_ri)
		call flux_cal(Q_ri)
		call innerintegral(Q_ri)
		call glflux

       end subroutine prep_advance

!------------------------------------------------------------------------------- 
 
       subroutine get_min_dt(dt)
       real dt,dt_min,dt_val(numprocs-1),tt,cfl,vmax,vmag,valf,vmag0,valf0,vex,vey,vez,vem,vem0,dni,dn,vx,vy,vz,Pr,sqdni,vacc,vacc0,cs
       integer :: i,j,main_proc=0,mpi_size=1,loc_reqs(numprocs-1),loc_stats(MPI_STATUS_SIZE,numprocs-1)

       vmag0 = 0.
       vmag = 0.
		
       do j=1,ny
       do i=1,nx
              
             dn = Q_r0(i,j,rh,1)
             dni = 1./dn
             sqdni = sqrt(dni)
             vx = Q_r0(i,j,mx,1)*dni
             vy = Q_r0(i,j,my,1)*dni
             vz = Q_r0(i,j,mz,1)*dni
             Pr = 0.6667*(Q_r0(i,j,en,1)*dni - 0.5*(vx**2 + vy**2 + vz**2))
             if(Pr < 0.) Pr = 0.
             cs = sqrt(aindex*Pr)
             vmag0 = max(abs(vx),abs(vy),abs(vz)) + cs
             if(vmag0 > vmag .and. dn > rh_mult*rh_floor) vmag = vmag0

       end do
       end do  

       vmax = max(clt,vmag)
       dt_min = cflm/(dxi*vmax)
       
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
	
	Q_r0(:,:,:,:) = 0.0
	Q_r0(:,:,rh,1) = rh_floor
	Q_r0(:,:,en,1) = T_floor*rh_floor/(aindex - 1.)	

        call fill_fluid

    end subroutine initial_condition

!-------------------------------------------------------

    subroutine fill_fluid
        
!-------------------------------------------------------

        integer i,j,ir,izw(4),ixw(4),iyw(4),iw,iseed,igrid,ieq,inds(nbasis)
        real wirer,x,y,x0,y0,rhline,wtev,rx,ry,rnum,w0,P,vx,vy,vz,rho
	real qquad(npg,nQ),xcc,ycc,zcc,rand_number
        iseed = 1317345*mpi_P + 5438432*mpi_Q 

        wtev = 2*T_floor
	w0 = 0.3
        coll = rh_fluid*wtev/vis

!       test problem is an unstable flow jet in x with velocity perturbations in y

        do i = 1,nx
           do j = 1,ny

              rnum = (ran(iseed) - 0.5)

		 xcc = xc(i)	
		 ycc = yc(j) 	
               
                 Q_r0(i,j,rh,1) = rh_fluid
		 Q_r0(i,j,my,1) = 0.0001*rnum/cosh(20*yc(j)/lyu)**2 
		 Q_r0(i,j,mz,1) = 0.  
		 Q_r0(i,j,mx,1) = 1.0*Q_r0(i,j,rh,1)/cosh(20*ycc/lyu)/1.
	         Q_r0(i,j,en,1) = (1. + 0.001*rnum)*wtev*Q_r0(i,j,rh,1)/(aindex - 1.)  &
                                + 0.5*(Q_r0(i,j,mx,1)**2 + Q_r0(i,j,my,1)**2 + Q_r0(i,j,mz,1)**2)/Q_r0(i,j,rh,1)

                 rho = rh_fluid
                 vx = Q_r0(i,j,mx,1)/rh_fluid
                 vy = Q_r0(i,j,my,1)/rh_fluid
                 vz = Q_r0(i,j,mz,1)/rh_fluid
                 P = wtev*Q_r0(i,j,rh,1)

                 Q_r0(i,j,exx,1) = P + rho*vx**2 
                 Q_r0(i,j,eyy,1) = P + rho*vy**2
                 Q_r0(i,j,ezz,1) = P + rho*vz**2
                 Q_r0(i,j,exy,1) = rho*vx*vy 
                 Q_r0(i,j,exz,1) = rho*vx*vz
                 Q_r0(i,j,eyz,1) = rho*vy*vz

           end do     
        end do
           
    end subroutine fill_fluid


!----------------------------------------------------------------------------------------------

    subroutine source_calc(Q_ri,t)
    implicit none
    integer i,j,ieq,ipg,ir
    real, dimension(nx,ny,nQ,nbasis) :: Q_ri 	
    real t,source(npg,nQ),Qin(nQ),dn,dni,Zin,vx,vy,vz,alpha,temp,dne,eta_a,Teev,Tiev,etaJ2,nuei,TemR,Tev,vmax,P
    real Tcoef,fac,en_floor,gyro,Pres,rh_buf,oth

	source_r(:,:,:,:) = 0.0
	source(:,:) = 0.0
        en_floor = P_floor/(aindex - 1.)
        oth = 1./3.
	
	do j = 1,ny
	do i = 1,nx

	do ipg = 1,npg
	
	do ieq = 1,nQ
	   Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,ieq,1:nbasis))  
	end do

        dn = Qin(rh)
        dni = 1./dn
        vx = Qin(mx)*dni
        vy = Qin(my)*dni
        vz = Qin(mz)*dni

!       Sources for the fluid variables. Nonzero if randomly forced.

        source(ipg,rh) = 0 
        source(ipg,mx) = 0 
        source(ipg,my) = 0 
        source(ipg,mz) = 0  
	source(ipg,en) = 0 

	! For ivis = 1 the sources for the deviatoric stress are zero for the implicit solve used in advance_time_level.
	! For ivis = 0 the sources are for the 6 equations for the energy tensor E_ij with the implicit solve implemented in advance_time_level.

        if(ivis .eq. 1)then
           source(ipg,exx) = 0 
           source(ipg,eyy) = 0 
           source(ipg,ezz) = 0     
           source(ipg,exy) = 0   
           source(ipg,exz) = 0    
           source(ipg,eyz) = 0 
        end if

        if(ivis .eq. 0)then
           source(ipg,exx) = coll*dn*(2*vx**2 - vy**2 - vz**2)*oth
           source(ipg,eyy) = coll*dn*(2*vy**2 - vz**2 - vx**2)*oth
           source(ipg,ezz) = coll*dn*(2*vz**2 - vx**2 - vy**2)*oth 
           source(ipg,exy) = coll*dn*vx*vy  
           source(ipg,exz) = coll*dn*vx*vz  
           source(ipg,eyz) = coll*dn*vy*vz 
        end if

	end do

	do ieq = exx,nQ
	  do ir=1,nbasis
	      source_r(i,j,ieq,ir) = 0.25*cbasis(ir)*sum(wgt2d(1:npg)*bfvals_int(1:npg,ir)*source(1:npg,ieq))
	  end do
	end do
	
	end do
	end do
	
    end subroutine source_calc

!----------------------------------------------------------------------------------------------

    subroutine implicit_source(Q_ri,Q_rip)
    implicit none
    integer i,j,k,ieq,ipg,ir
    real, dimension(nx,ny,nQ,nbasis) :: Q_ri,Q_rip 
    real Qout(npg,nQ)	
    real t,source(npg,nQ),Qin(nQ),dn,dni,Zin,vx,vy,vz,alpha,temp,dne,eta_a,Teev,Tiev,etaJ2,nuei,TemR,Tev,vmax
    real Tcoef,fac,en_floor,gyro,P,rh_buf,oth

        oth = 1./3.
        fac = 1./(1. + coll*dt)

	do j = 1,ny
	do i = 1,nx

	do ipg = 1,npg
	
	do ieq = 1,nQ
	   Qin(ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_ri(i,j,ieq,1:nbasis)) 
	end do
		
        dn = Qin(rh)
	dni = 1./Qin(rh)
        vx = Qin(mx)*dni
        vy = Qin(my)*dni
        vz = Qin(mz)*dni
	P = (aindex - 1.)*(Qin(en) - 0.5*dn*(vx**2 + vy**2 + vz**2)) 

        Qout(ipg,exx) = (Qin(exx) + coll*dt*(P + dn*vx**2))*fac
        Qout(ipg,eyy) = (Qin(eyy) + coll*dt*(P + dn*vy**2))*fac
        Qout(ipg,ezz) = (Qin(ezz) + coll*dt*(P + dn*vz**2))*fac
        Qout(ipg,exy) = (Qin(exy) + coll*dt*dn*vx*vy)*fac
        Qout(ipg,exz) = (Qin(exz) + coll*dt*dn*vx*vz)*fac
        Qout(ipg,eyz) = (Qin(eyz) + coll*dt*dn*vy*vz)*fac

	end do

	do ieq=exx,nQ	
	   do ir=1,nbasis
	      Q_rip(i,j,ieq,ir) = 0.25*cbasis(ir)*sum(wgt2d(1:npg)*bfvals_int(1:npg,ir)*Qout(1:npg,ieq))
	   end do
	end do
	
	end do
	end do
	
    end subroutine implicit_source

!----------------------------------------------------------------------------------------------

    subroutine advance_time_level_gl(Q_ri,Q_rp)
    implicit none
    integer i,j,k,ieq,ir
    real, dimension(nx,ny,nQ,nbasis) :: Q_ri, Q_rp
    real Q_xtemp, Q_ytemp, Q_ztemp, dti, oth, fac
    real c1d3, c1d5

	oth = 1./3.
        dti = 1./dt
        fac = 1.
        if(ivis .eq. 1)fac = 1./(1. + coll*dt)

	do j = 1,ny
	do i = 1,nx

	do ieq = 1,en
	  do ir=1,nbasis
	    Q_rp(i,j,ieq,ir) = Q_ri(i,j,ieq,ir) - dt*glflux_r(i,j,ieq,ir) + dt*source_r(i,j,ieq,ir)
	  end do
	end do

	do ieq = exx,nQ
	  do ir=1,nbasis
	    Q_rp(i,j,ieq,ir) = (Q_ri(i,j,ieq,ir) - dt*glflux_r(i,j,ieq,ir) + dt*source_r(i,j,ieq,ir))*fac
	  end do
	end do

        do ieq = 1,nQ
        if ( Q_rp(i,j,ieq,1) .ne. Q_rp(i,j,ieq,1)) then
	        print *,'NaN. Bailing out...','  xc  =',xc(i),'  yc  =',yc(j),'  ieq  =',ieq
        call exit(-1)
        endif
        end do

	end do
	end do


    end subroutine advance_time_level_gl

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
    real E_xx,E_yy,E_zz,E_xy,E_xz,E_yz,P_10,P_5,Pxx,Pyy,Pzz,Pxy,Pxz,Pyz
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
        Pxx = E_xx - dn*vx**2
        Pyy = E_yy - dn*vy**2
        Pzz = E_zz - dn*vz**2
        Pxy = E_xy - dn*vx*vy
        Pxz = E_xz - dn*vx*vz
        Pyz = E_yz - dn*vy*vz
	
	P = (aindex - 1.)*(Qpnts_r(ife,en) - 0.5*dn*(vx**2 + vy**2 + vz**2)) 
        P_10 = (Pxx + Pyy + Pzz)/3.
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

	end if

	if(ivis .ge. 0)then

        if(ixyz .eq. 1)then

	fpnts_r(ife,rh) = Qpnts_r(ife,mx)
	fpnts_r(ife,en) = (Qpnts_r(ife,en) + Pxx)*vx + vy*Pxy + vz*Pxz
 
        if(ivis .eq. 1)then
  	   fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vx + E_xx + P
	   fpnts_r(ife,my) = Qpnts_r(ife,my)*vx + E_xy    
	   fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vx + E_xz
	   fpnts_r(ife,exx) =  c4d3*nu*vx 
	   fpnts_r(ife,eyy) = -c2d3*nu*vx 
	   fpnts_r(ife,ezz) = -c2d3*nu*vx 
	   fpnts_r(ife,exy) = nu*vy
	   fpnts_r(ife,exz) = nu*vz
	   fpnts_r(ife,eyz) = 0
        else
	   fpnts_r(ife,mx) = E_xx - P_10 + P_5
	   fpnts_r(ife,my) = E_xy   
	   fpnts_r(ife,mz) = E_xz 
	   fpnts_r(ife,exx) = vx*(dn*vx**2 + 3*Pxx)                                        ! term 1
	   fpnts_r(ife,eyy) = vx*(dn*vy**2 + Pyy) + 2*vy*Pxy                               ! term 4
	   fpnts_r(ife,ezz) = vx*(dn*vz**2 + Pzz) + 2*vz*Pxz                               ! term 7
	   fpnts_r(ife,exy) = vy*(dn*vx**2 + Pxx) + 2*vx*Pxy                               ! term 10
	   fpnts_r(ife,exz) = vz*(dn*vx**2 + Pxx) + 2*vx*Pxz                               ! term 13
	   fpnts_r(ife,eyz) = dn*vx*vy*vz + vx*Pyz + vy*Pxz * vz*Pxy                       ! term 16  
        end if

        end if

        if(ixyz .eq. 2)then

	fpnts_r(ife,rh) = Qpnts_r(ife,mxa(ixyz)) 
	fpnts_r(ife,en) = (Qpnts_r(ife,en) + Pyy)*vy + vx*Pxy + vz*Pyz

        if(ivis .eq. 1)then
	   fpnts_r(ife,mx) = Qpnts_r(ife,mx)*vy + E_xy    
	   fpnts_r(ife,my) = Qpnts_r(ife,my)*vy + E_yy + P 
	   fpnts_r(ife,mz) = Qpnts_r(ife,mz)*vy + E_yz 
	   fpnts_r(ife,exx) = -c2d3*nu*vy
	   fpnts_r(ife,eyy) =  c4d3*nu*vy
	   fpnts_r(ife,ezz) = -c2d3*nu*vy
	   fpnts_r(ife,exy) = nu*vx
	   fpnts_r(ife,exz) = 0
	   fpnts_r(ife,eyz) = nu*vz
        else
	   fpnts_r(ife,mx) = E_xy    
	   fpnts_r(ife,my) = E_yy - P_10 + P_5
	   fpnts_r(ife,mz) = E_yz 
	   fpnts_r(ife,exx) = vy*(dn*vx**2 + Pxx) + 2*vx*Pxy                            ! term 2
	   fpnts_r(ife,eyy) = vy*(dn*vy**2 + 3*Pyy)                                     ! term 5
	   fpnts_r(ife,ezz) = vy*(dn*vz**2 + Pzz) + 2*vz*Pyz                            ! term 8
	   fpnts_r(ife,exy) = vx*(dn*vy**2 + Pyy) + 2*vy*Pxy                            ! term 11
	   fpnts_r(ife,exz) = dn*vx*vy*vz + vx*Pyz + vy*Pxz + vz*Pxy                    ! term 14
	   fpnts_r(ife,eyz) = vz*(dn*vy**2 + Pyy) + 2*vy*Pyz                            ! term 17
        end if

        end if

        end if
		
	end do

    end subroutine

!----------------------------------------------------------------------------------------------

    subroutine flux_cal(Q_r)
    implicit none
    integer i,j,ieq,iback,jleft,i4,i4p,ipnt
    real Qedge_x(nfe,nQ), Qedge_y(nfe,nQ),fedge_x(nfe,nQ), fedge_y(nfe,nQ)
    real, dimension(nx,ny,nQ,nbasis) :: Q_r
    real cwavex(nfe),cwavey(nfe),cwavez(nfe),Pre,dni,B2, cfbound, cfx(nedge,nQ),cfy(nedge,nQ)
    real fhllc_x(nedge,5),fhllc_y(nedge,5),fs(nedge,nQ),qvin(nQ)

!--------------------------------------------------------------

	do j=1,ny
	do i=1,nx+1

	iback = i-1
	
	if (i .gt. 1) then
	   do ieq = 1,nQ
           do ipnt=1,nedge	
	      Qedge_x(ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_r(iback,j,ieq,1:nbasis))
           end do
	   end do
	end if

	if (i .eq. 1) then
	   do ieq = 1,nQ	
	      Qedge_x(1:nedge,ieq) = Qxlow_ext(j,1:nedge,ieq)
	   end do
	end if

	if (i .lt. nx+1) then
	   do ieq = 1,nQ	
	   do ipnt=1,nedge
	      Qedge_x(ipnt+nedge,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_r(i,j,ieq,1:nbasis))
	   end do
	   end do
	end if

	if (i .eq. nx+1) then
	   do ieq = 1,nQ	
	      Qedge_x(nedge+1:nfe,ieq) = Qxhigh_ext(j,1:nedge,ieq)
	   end do
	end if

	   call flux_calc_pnts_r(Qedge_x,fedge_x,1,nfe)

	   do i4=1,nfe
	      do ieq=1,nQ
	         qvin(ieq) = Qedge_x(i4,ieq)
	      end do
              cwavex(i4) = cfcal(qvin,1)
           end do

           do i4=1,nedge
              cfrx(i4,rh:en) = max(cwavex(i4),cwavex(i4+nedge))
           end do

	   do ieq = 1,nQ
	      do i4=1,nedge
	         i4p = i4 + nedge
	         flux_x(i4,i,j,ieq) = 0.5*(fedge_x(i4,ieq) + fedge_x(i4p,ieq)) - 0.5*cfrx(i4,ieq)*(Qedge_x(i4p,ieq) - Qedge_x(i4,ieq))
	      end do
	   end do

	   if (ihllc .eq. 1) then

              call flux_hllc(Qedge_x,fedge_x,fhllc_x,1)

	      do ieq = 1,en
	         do i4=1,nedge
	            flux_x(i4,i,j,ieq) = fhllc_x(i4,ieq)
	         end do 
	      end do

           end if

	end do
	end do

!-----------------------------------------------------------------------------------------------

	do j=1,ny+1
	   jleft = j-1
	do i=1,nx

	if (j .gt. 1) then
	   do ieq = 1,nQ
	   do ipnt=1,nedge
	      Qedge_y(ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_r(i,jleft,ieq,1:nbasis))
	   end do
	   end do
	end if

	if (j .eq. 1) then
	   do ieq = 1,nQ
	      Qedge_y(1:nedge,ieq) = Qylow_ext(i,1:nedge,ieq)
	   end do
	end if

	if (j .lt. ny+1) then
	   do ieq = 1,nQ
	   do ipnt=1,nedge
	      Qedge_y(ipnt+nedge,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_r(i,j,ieq,1:nbasis))
	   end do
	   end do
	end if

	if (j .eq. ny+1) then
	   do ieq = 1,nQ
	      Qedge_y(nedge+1:nfe,ieq) = Qyhigh_ext(i,1:nedge,ieq)
	   end do
	end if

	   call flux_calc_pnts_r(Qedge_y,fedge_y,2,nfe)

	   do i4=1,nfe
	      do ieq=1,nQ
	         qvin(ieq) = Qedge_y(i4,ieq)
	      end do
              cwavey(i4) = cfcal(qvin,2)
           end do

           do i4=1,nedge
              cfry(i4,rh:en) = max(cwavey(i4),cwavey(i4+nedge))            
           end do

	   do ieq = 1,nQ
	   do i4=1,nedge
	      i4p = i4 + nedge
	      flux_y(i4,i,j,ieq) = 0.5*(fedge_y(i4,ieq) + fedge_y(i4p,ieq)) - 0.5*cfry(i4,ieq)*(Qedge_y(i4p,ieq) - Qedge_y(i4,ieq))
	   end do
	   end do

	   if (ihllc .eq. 1) then

              call flux_hllc(Qedge_y,fedge_y,fhllc_y,2)

	      do ieq = 1,en
	         do i4=1,nedge
	            flux_y(i4,i,j,ieq) = fhllc_y(i4,ieq)
	         end do 
	      end do

    	   end if

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
	real Qlr(nfe,nQ),flr(nfe,nQ),fhllc(nedge,5)
	real sm_num(nedge),sm_den(nedge),qtilde(nedge,5),rtrho(nfe),rtrho_i(nedge),qsq(nfe)
	real s_lr(nfe),ctilde(nedge),hlr(nfe),cslr(nfe),ctsq(nedge),Zi,mfact,qfloor(nedge)
	real aq(nfe),bq(nfe),Qstar(nfe,6),fstar(nfe,6),pstar(nedge),s_m(nedge)
	real rhov(nfe),vlr(nfe,3),plr(nfe),slsm_i(nedge),rho_i,qslr(nfe),sq_lr(nfe),slrm_i(nfe),B2(nfe),cs(nfe),cs2(nfe)
        real dn,u2,pe,Temp(nfe),targ,theta_t,dnnu1,dnnu2,P
	integer ixyz,i4,i4p1,nr,jie,k,k2,ieq,iparr,iperp1,iperp2,ibatten
	integer rhj,mxj,myj,mzj,enj,psj,ivar(5),ipassive,nhll,ib1,ib2

        nhll = 5
	rhj = rh
	mxj = mx
	myj = my
	mzj = mz
	enj = en

	iparr = mxa(ixyz)
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
              dn = rhov(k)
              u2 = qsq(k)
              pe = Qlr(k,enj)
              plr(k) = 0.6667*(pe - 0.5*dn*u2)
              if(plr(k) < 0.) plr(k) = 0.
              cs(k) = sqrt(1.6667*plr(k)/dn)
	   end do

	   do k=1,nedge
	      k2 = k + nedge
	         cslr(k) = vlr(k,1)   - cs(k)	       
	         cslr(k2) = vlr(k2,1) + cs(k2)	
	         qslr(k) = vlr(k2,1) - cs(k2)
	         qslr(k2) = vlr(k,1) + cs(k)
	   end do

!	  Calculate the slow and fast wavespeeds S_L, S_R of the left, right propagating shocks.

	    do k=1,nedge
	      k2 = k + nedge
	      s_lr(k) = min(cslr(k),qslr(k))         ! S_L = min(lambda_M(Q_l),lambda_M(Q_r or Q_Roe))
	      s_lr(k2) = max(cslr(k2),qslr(k2))	     ! S_R = max(lambda_P(Q_r),lambda_P(Q_l or Q_Roe))
	      sm_num(k) = rhov(k2)*vlr(k2,1)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*vlr(k,1)*(s_lr(k) - vlr(k,1))
	      sm_num(k) = sm_num(k) + plr(k) - plr(k2)
	      sm_den(k) = rhov(k2)*(s_lr(k2) - vlr(k2,1)) - rhov(k)*(s_lr(k) - vlr(k,1))
	    end do
	   where (sm_den .eq. 0.0) sm_den = rh_floor

!	   Calculate the wavespeed S_M of the contact discontinuity.
	
	   do k=1,nedge
	      s_m(k) = sm_num(k)/sm_den(k)	     	! Eq. (34) of Batten, 1997
	      pstar(k) = rhov(k)*(vlr(k,1) - s_lr(k))*(vlr(k,1) - s_m(k)) + plr(k) 		
							! Eq. (36) of Batten, 1997
	   end do


!	   Now, calculate Q_l* and Q_r* in order to calculate F_l* and F_r*.

	   do k=1,nfe
	      if (k .le. nedge) i4 = k
	      if (k .gt. nedge) i4 = k - nedge
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

	do i4 = 1,nedge						! Use Eq. (26) of Batten ,1997
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
	    if (s_m(i4) .le. 0.0 .and. 0.0 .le. s_lr(i4+nedge)) then    	! if S_M <= 0 <= S_R
	      do ieq=1,nhll
		fhllc(i4,ivar(ieq)) = fstar(i4+nedge,ieq)		! F_HLLC = F_r*
	      end do
	    end if
	    if (s_lr(i4+nedge) .lt. 0.0) then			! if S_R < 0
	      do ieq=1,nhll
		fhllc(i4,ivar(ieq)) = flr(i4+nedge,ivar(ieq))	! F_HLLC = F_r
	      end do
	    end if

	end do

	end subroutine flux_hllc


!----------------------------------------------------------------------------------------------

    subroutine innerintegral(Q_r)
    implicit none
    real, dimension(nx,ny,nQ,nbasis) :: Q_r
    integer i,j,ieq,ipg,ir,ebz(2),jr
    real Qinner(npg,nQ),finner_x(npg,nQ), finner_y(npg,nQ),rinv,rwgt(npg),int_r(nbastot)

	integral_r(:,:,:,:) = 0.
	rwgt(1:npg) = 1.0
	rinv = 1.0
	
	do j = 1,ny
	do i = 1,nx
	
	   do ieq = 1,nQ
	   do ipg = 1,npg
	      Qinner(ipg,ieq) = sum(bfvals_int(ipg,1:nbasis)*Q_r(i,j,ieq,1:nbasis))
	   end do
	   end do

	   call flux_calc_pnts_r(Qinner,finner_x,1,npg)
	   call flux_calc_pnts_r(Qinner,finner_y,2,npg)

		
	   do ieq = 1,nQ

	      int_r(kx) = 0.5*cbasis(kx)*dxi*sum(rwgt(1:npg)*wgt2d(1:npg)*finner_x(1:npg,ieq))
	      int_r(ky) = 0.5*cbasis(ky)*dyi*sum(rwgt(1:npg)*wgt2d(1:npg)*finner_y(1:npg,ieq))

		 if (nbasis .gt. 3) then
	          int_r(kxy) = 0.5*cbasis(kxy)*sum(rwgt(1:npg)*wgt2d(1:npg)*(bfvals_int(1:npg,ky)*dxi*finner_x(1:npg,ieq) &
                             + bfvals_int(1:npg,kx)*dyi*finner_y(1:npg,ieq)))
		 end if

		 if (iquad .ge. 3) then
		  int_r(kxx) = 0.5*cbasis(kxx)*dxi*sum(rwgt(1:npg)*wgt2d(1:npg)*3.*bfvals_int(1:npg,kx)*finner_x(1:npg,ieq))
	          int_r(kyy) = 0.5*cbasis(kyy)*dyi*sum(rwgt(1:npg)*wgt2d(1:npg)*3.*bfvals_int(1:npg,ky)*finner_y(1:npg,ieq))
	         end if	

		 if (nbasis .gt. 6) then

		   int_r(kxxy) = 0.5*cbasis(kxxy)*sum(rwgt(1:npg)*wgt2d(1:npg)*(bfvals_int(1:npg,kxx)*dyi*finner_y(1:npg,ieq) &
                               + 3.*bfvals_int(1:npg,kxy)*dxi*finner_x(1:npg,ieq)))	
		   int_r(kxyy) = 0.5*cbasis(kxyy)*sum(rwgt(1:npg)*wgt2d(1:npg)*(bfvals_int(1:npg,kyy)*dxi*finner_x(1:npg,ieq) &
                               + 3.*bfvals_int(1:npg,kxy)*dyi*finner_y(1:npg,ieq)))

		   if (ibitri .eq. 1) then
		      int_r(kxxyy) = 0.5*cbasis(kxxyy)*sum(rwgt(1:npg)*wgt2d(1:npg)*(3.*bfvals_int(1:npg,kxyy)*dxi*finner_x(1:npg,ieq) &
                                   + 3.*bfvals_int(1:npg,kxxy)*dyi*finner_y(1:npg,ieq)))
		   end if

		   if (iquad .ge. 4) then
		     int_r(kxxx) = 0.5*cbasis(kxxx)*dxi*sum(rwgt(1:npg)*wgt2d(1:npg)*(7.5*bfvals_int(1:npg,kx)**2 - 1.5)*finner_x(1:npg,ieq))	
		     int_r(kyyy) = 0.5*cbasis(kyyy)*dyi*sum(rwgt(1:npg)*wgt2d(1:npg)*(7.5*bfvals_int(1:npg,ky)**2 - 1.5)*finner_y(1:npg,ieq))
		   end if

		 end if

		 do ir=1,nbasis
			integral_r(i,j,ieq,ir) = int_r(ir)
		 end do
           end do
		
		
	end do
	end do
	
    end subroutine innerintegral

!----------------------------------------------------------------------------------------------

    subroutine glflux
    implicit none
    real rbp,rbm,rinv,redge(nedge)
    integer i,j,ieq,ir
	
	rbp = 1.0
	rbm = 1.0
	rinv = 1.0
	redge(1:nedge) = 1.0

	
	do ieq = 1,nQ
	do j = 1,ny
	do i = 1,nx
	   
	   glflux_r(i,j,ieq,1) =  0.5*(rinv*dxi*sum(wgt1d(1:nedge)*(rbp*flux_x(1:nedge,i+1,j,ieq) - rbm*flux_x(1:nedge,i,j,ieq))) &
				     + rinv*dyi*sum(redge(1:nedge)*wgt1d(1:nedge)*(flux_y(1:nedge,i,j+1,ieq) - flux_y(1:nedge,i,j,ieq))))
					  
	do ir=2,nbasis

	    glflux_r(i,j,ieq,ir) =  0.5*cbasis(ir)*(rinv*dxi*sum(wgt1d(1:nedge)*(rbp*bfvals_xp(1:nedge,ir)*flux_x(1:nedge,i+1,j,ieq) &
                                 - rbm*bfvals_xm(1:nedge,ir)*flux_x(1:nedge,i,j,ieq)))                                               &
		                 + rinv*dyi*sum(redge(1:nedge)*wgt1d(1:nedge)*(bfvals_yp(1:nedge,ir)*flux_y(1:nedge,i,j+1,ieq)       &
                                 - bfvals_ym(1:nedge,ir)*flux_y(1:nedge,i,j,ieq)))) - integral_r(i,j,ieq,ir) 
	end do

 
	end do
	end do
	end do


    end subroutine glflux

!-----------------------------------------------------------!
!*******************calculate freezing speeds***************!
!-----------------------------------------------------------!        
        real function cfcal(Qin,cases)
        implicit none
        integer cases
        real Qin(nQ)
    	real dn,dni,vx,vy,vz,u2,pe,B2,P,Temp
            
        dn = Qin(rh)
        dni = 1./dn
        vx = Qin(mx)*dni 
        vy = Qin(my)*dni
        vz = Qin(mz)*dni
        u2 = vx**2 + vy**2 + vz**2
        pe = Qin(en)
        P = 0.6667*(pe - 0.5*dn*u2)
        if(P < 0.) P = 0.

	select case (cases)

   	case (1) !freezing speed in x direction for fluid variable 
	      cfcal = abs(vx) + sqrt(1.6667*P*dni)
	case (2) !freezing speed in y direction for fluid variable
      	      cfcal = abs(vy) + sqrt(1.6667*P*dni)

        end select
        
        end function cfcal 

!----------------------------------------------------------------------------------------------

    subroutine limit_flow(Q_r)
    implicit none
    integer i,j,ieq,ipg,ir
    real, dimension(nx,ny,nQ,nbasis) :: Q_r
    real dn,en_floor,jmax,tev,vx,vy,vz,dni

        en_floor = T_floor/(aindex-1.)

	do j = 1,ny
	do i = 1,nx

	   dn = Q_r(i,j,rh,1)
           dni = 1./dn
	   vx = Q_r(i,j,mx,1)*dni
	   vy = Q_r(i,j,my,1)*dni
	   vz = Q_r(i,j,mz,1)*dni

       	   if (dn < rh_mult*rh_floor .and. .true.) then
            	Q_r(i,j,en,1) = en_floor*rh_floor
		Q_r(i,j,mx:mz,1) = 0.0
		do ir=2,nbasis
	           Q_r(i,j,mx:en,ir) = 0.0
		end do
	   end if

	end do
	end do

    end subroutine

!----------------------------------------------------------------------------------

    subroutine limiter(Q_r)
    implicit none
    integer i, j, ieq, ipge, minindex, ir
    real, dimension(nx,ny,nQ,nbasis) :: Q_r
    real Qedge(nslim,nQ),theta,Qmin(nQ), deltaQ(nQ)
    real epsi, Qrhmin, QPmin, P(nslim), Pave, dn, dni, epsiP, thetaj
    real*8 a, b, c

    epsi = rh_floor 
    epsiP = rh_floor*T_floor

    do j = 1,ny
    do i = 1,nx
        dn = Q_r(i,j,rh,1)
	if (Q_r(i,j,rh,1) < rh_floor) then

		do ir=2,nbasis
		   Q_r(i,j,rh:en,ir) = 0.0	
	    	end do	   

        	Q_r(i,j,rh,1) = rh_floor

	else	
		do ipge = 1,nslim
		   Qedge(ipge,rh) = sum(bf_edges(ipge,1:nbasv(rh))*Q_r(i,j,rh,1:nbasv(rh)))
		end do
		
		Qrhmin = minval(Qedge(:,rh))
		if (Qrhmin < epsi) then
			theta = (epsi-Q_r(i,j,rh,1))/(Qrhmin-Q_r(i,j,rh,1))

			if (theta .gt. 1.) theta = 1.
			if (theta .lt. 0) theta = 0.

			do ir=2,nbasis
				Q_r(i,j,rh,ir) = abs(theta)*Q_r(i,j,rh,ir)
			end do

		end if	

                if(dn .lt. 0.9)then
			
		Pave = (aindex-1.)*(Q_r(i,j,en,1) - 0.5*(Q_r(i,j,mx,1)**2 + Q_r(i,j,my,1)**2 + Q_r(i,j,mz,1)**2)/Q_r(i,j,rh,1))
	
		if (Pave < epsiP) then

			do ir=2,nbasis
			   Q_r(i,j,rh:en,ir) = 0.0	
			end do
	
		else
			theta = 1.
			do ipge = 1,nslim 
			   do ieq = 1,5
			     Qedge(ipge,ieq) = sum(bf_edges(ipge,1:nbasis)*Q_r(i,j,ieq,1:nbasis))
			   end do
			
			   dn = Qedge(ipge,rh)
        		   dni = 1./dn
			   P(ipge) = (aindex - 1.)*(Qedge(ipge,en) - 0.5*(Qedge(ipge,mx)**2+Qedge(ipge,my)**2+Qedge(ipge,mz)**2)*dni)

			   if (P(ipge) < epsiP) then
                  	      if (Pave .ne. P(ipge)) then
			         thetaj = (Pave - epsiP)/(Pave - P(ipge))
				 theta = min(theta,thetaj)
			      end if
		           end if
			end do

			if (theta .gt. 1.) theta = 1.
			if (theta .lt. 0.) theta = 0.

			do ir=2,nbasis
			   Q_r(i,j,rh:en,ir) = theta*Q_r(i,j,rh:en,ir) 
			end do

		end if
                end if
	end if
    end do
    end do


    end subroutine limiter

!----------------------------------------------------------------------------------

    subroutine limiter_vtk(Q_r)
    implicit none
    integer i, j, ieq, ipge, minindex, ir
    real, dimension(nx,ny,nQ,nbasis) :: Q_r
    real Qedge(nslim,nQ),theta,Qmin(nQ), deltaQ(nQ)
    real epsi, Qrhmin, QPmin, P(nslim), Pave, dn, dni, epsiP, thetaj
    real*8 a, b, c

    epsi = rh_floor 
    epsiP = rh_floor*T_floor

    do j = 1,ny
    do i = 1,nx
	if (Q_r(i,j,rh,1) < rh_floor) then

		do ir=2,nbasis
		   Q_r(i,j,rh:en,ir) = 0.0	
	    	end do	   

        	Q_r(i,j,rh,1) = rh_floor

	else	
		do ipge = 1,nvtk2
		   Qedge(ipge,rh) = sum(bfvtk(ipge,1:nbasv(rh))*Q_r(i,j,rh,1:nbasv(rh)))
		end do
		
		Qrhmin = minval(Qedge(1:nvtk2,rh))
		if (Qrhmin < epsi) then
			theta = (epsi-Q_r(i,j,rh,1))/(Qrhmin-Q_r(i,j,rh,1))

			if (theta .gt. 1.) theta = 1.
			if (theta .lt. 0) theta = 0.

			do ir=2,nbasis
				Q_r(i,j,rh,ir) = abs(theta)*Q_r(i,j,rh,ir)
			end do

		end if	
			
		Pave = (aindex-1.)*(Q_r(i,j,en,1) - 0.5*(Q_r(i,j,mx,1)**2 + Q_r(i,j,my,1)**2 + Q_r(i,j,mz,1)**2)/Q_r(i,j,rh,1))
	
		if (Pave < epsiP) then

			do ir=2,nbasis
			   Q_r(i,j,rh:en,ir) = 0.0	
			end do
	
		else
			theta = 1.
			do ipge = 1,nvtk2 
			   do ieq = 1,5
			     Qedge(ipge,ieq) = sum(bfvtk(ipge,1:nbasis)*Q_r(i,j,ieq,1:nbasis))
			   end do
			
			   dn = Qedge(ipge,rh)
        		   dni = 1./dn
			   P(ipge) = (aindex - 1.)*(Qedge(ipge,en) - 0.5*(Qedge(ipge,mx)**2+Qedge(ipge,my)**2+Qedge(ipge,mz)**2)*dni)

			   if (P(ipge) < epsiP) then
                  	      if (Pave .ne. P(ipge)) then
			         thetaj = (Pave - epsiP)/(Pave - P(ipge))
				 theta = min(theta,thetaj)
			      end if
		           end if
			end do

			if (theta .gt. 1.) theta = 1.
			if (theta .lt. 0.) theta = 0.

			do ir=2,nbasis
			   Q_r(i,j,rh:en,ir) = theta*Q_r(i,j,rh:en,ir) 
			end do

		end if
	end if
    end do
    end do


    end subroutine limiter_vtk

!----------------------------------------------------------------------------------------------

    subroutine prepare_exchange(Q_r) 
	integer ieq, i, j, ipnt
        real, dimension(nx,ny,nQ,nbasis) :: Q_r

	
	do ieq = 1,nQ
	do i = 1,nx	
      do ipnt=1,nedge	
		Qylow_int(i,ipnt,ieq) = sum(bfvals_ym(ipnt,1:nbasis)*Q_r(i,1,ieq,1:nbasis))
		Qyhigh_int(i,ipnt,ieq) = sum(bfvals_yp(ipnt,1:nbasis)*Q_r(i,ny,ieq,1:nbasis))
	  end do
	end do
	end do
	
	
	do ieq = 1,nQ
	do j = 1,ny	
      do ipnt=1,nedge	
		Qxlow_int(j,ipnt,ieq) = sum(bfvals_xm(ipnt,1:nbasis)*Q_r(1,j,ieq,1:nbasis))
		Qxhigh_int(j,ipnt,ieq) = sum(bfvals_xp(ipnt,1:nbasis)*Q_r(nx,j,ieq,1:nbasis))
	  end do		
	end do
	end do	
	call exchange_flux	

    end subroutine

!----------------------------------------------------------------------------------------------

    subroutine exchange_flux
        integer mpi_size 
    
	call MPI_BARRIER(cartcomm,ierr)

        mpi_size=nedge*nx*nQ
       
        if (nbrs(UP) .ne. MPI_PROC_NULL) then
           call MPI_ISend(Qyhigh_int,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(1),ierr)
	endif

        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
           call MPI_IRecv(Qylow_ext,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(2),ierr)
	endif

        if (nbrs(UP) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(1),stats(:,1),ierr)
           call MPI_IRecv(Qyhigh_ext,mpi_size,MPI_TT,nbrs(UP),0,cartcomm,reqs(3),ierr)
	endif

        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(2),stats(:,2),ierr)
           call MPI_ISend(Qylow_int,mpi_size,MPI_TT,nbrs(DOWN),0,cartcomm,reqs(4),ierr)        			
	endif
		
        if (nbrs(UP) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(3),stats(:,3),ierr)
	endif

        if (nbrs(DOWN) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(4),stats(:,4),ierr)
	endif


 	mpi_size=ny*nedge*nQ 

        if (nbrs(RIGHT) .ne. MPI_PROC_NULL) then
            call MPI_ISend(Qxhigh_int,mpi_size,MPI_TT,nbrs(RIGHT),0,cartcomm,reqs(1),ierr)
	endif

        if (nbrs(LEFT) .ne. MPI_PROC_NULL) then
            call MPI_IRecv(Qxlow_ext,mpi_size,MPI_TT,nbrs(LEFT),0,cartcomm,reqs(2),ierr)
	endif
	
        if (nbrs(RIGHT) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(1),stats(:,1),ierr)
           call MPI_IRecv(Qxhigh_ext,mpi_size,MPI_TT,nbrs(RIGHT),0,cartcomm,reqs(3),ierr)
	endif

        if (nbrs(LEFT) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(2),stats(:,2),ierr)
           call MPI_ISend(Qxlow_int,mpi_size,MPI_TT,nbrs(LEFT),0,cartcomm,reqs(4),ierr)        			
	endif
		
        if (nbrs(RIGHT) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(3),stats(:,3),ierr)
	endif

        if (nbrs(LEFT) .ne. MPI_PROC_NULL) then
           call MPI_Wait(reqs(4),stats(:,4),ierr)
	endif

		
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

        real function rz(i,j)
        integer i,j
        rz = sqrt(yc(j)**2)
        end function rz

!-----------------------------------------------------------

        real function r(i,j)
        integer i,j
        r = sqrt(xc(i)**2 + yc(j)**2)
        end function r

!-----------------------------------------------------------

        real function theta(i,j)
        integer i,j
        theta = atan2(yc(j),xc(i))
        end function theta


!-----------------------------------------------------------

        real function x_edge(i,ipnt)
        integer i,ipnt,ix
        x_edge = xc(i) + bfvals_yp(ipnt,kx)*0.5*dx
        end function x_edge

!-----------------------------------------------------------

        real function y_edge(j,ipnt)
        integer j,ipnt
        y_edge = yc(j) + bfvals_xp(ipnt,ky)*0.5*dy
        end function y_edge

!-----------------------------------------------------------

        real function x_int(i,ipnt)
        integer i,ipnt
        x_int = xc(i) + bfvals_int(ipnt,kx)*0.5*dx
        end function x_int

!-----------------------------------------------------------

        real function y_int(j,ipnt)
        integer j,ipnt
        y_int = yc(j) + bfvals_int(ipnt,ky)*0.5*dy
        end function y_int

!-----------------------------------------------------------

        real function r_xy(xq,yq)
        real xq,yq
        r_xy= sqrt(xq**2 + yq**2)
        end function r_xy

!-----------------------------------------------------------

        real function theta_xy(xq,yq)
        real xq,yq
        theta_xy = atan2(yq,xq)
        end function theta_xy


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
       
        subroutine output_vtk(Qin,nout,iam)
        implicit none
        real Qin(nx,ny,nQ,nbasis)
	integer nout
	integer(I4P) , parameter :: nb=1*ngu,sprd=0*1
	real, dimension(:),allocatable :: x_xml_rect
	real, dimension(:),allocatable:: y_xml_rect
	real, dimension(1):: z_xml_rect
	real, dimension(:),allocatable :: var_xml_val_x
	real, dimension(:),allocatable :: var_xml_val_y
	real, dimension(:),allocatable :: var_xml_val_z
	real qvtk(nnx,nny,nQ),qvtk_dxmx(nnx,nny),qvtk_dymy(nnx,nny),dxmy,dymx,dxrh,dyrh,qvtk_dxvy(nnx,nny),qvtk_dyvx(nnx,nny)
	real P, vx, vy, vz,  dni, dn, u2, pe, dnt,Temp,targ,theta_t, Zi
	integer :: mil=1,miu=1,mju=1,mjl=1!,mpi_P=1,mpi_Q=1
	integer(I4P):: E_IO,i,j,l,num,iam,ira(nnx),jra(nny),iba(nnx),jba(nny),igrid,ieq
        character (50) :: out_name
        character (4) :: tname
        character (5) :: tname1
        character (4) :: pname
        character (5) :: pname1
        
       
       allocate(x_xml_rect(nnx+1))
       allocate(y_xml_rect(nny+1))
       allocate(var_xml_val_x(nnx*nny))
       allocate(var_xml_val_y(nnx*nny))
       allocate(var_xml_val_z(nnx*nny))
       
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
!---------------------------------------------------------------------
	out_name='/data/data6/perseus_p'//pname//'_t'//tname//'.vtr'
!	out_name='data/perseus_p'//pname//'_t'//tname//'.vtr'
!---------------------------------------------------------------------
        out_name = trim(out_name)
        out_name = adjustr(out_name) 

		call limiter_vtk(Qin)
		call limit_flow(Qin)
		
		
		E_IO = VTK_INI_XML(output_format = 'BINARY',              &
		                   filename      = out_name, &
		                   mesh_topology = 'RectilinearGrid',     &
		                   nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=0,nz2=0)

		do i=1+nb,nnx-nb+1
			x_xml_rect(i-nb)=(xvtk(i) - 0.5*dxvtk)*L0
		enddo
		do j=1+nb,nny-nb+1
			y_xml_rect(j-nb)=(yvtk(j) - 0.5*dyvtk)*L0
		enddo
		z_xml_rect(:)=0.

		E_IO = VTK_GEO_XML(nx1=1,nx2=nnx+1,ny1=1,ny2=nny+1,nz1=0,nz2=0, &
			                 X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)
		
		E_IO = VTK_DAT_XML(var_location     = 'cell', &
		                   var_block_action = 'OPEN')

			do i=1,nnx
				ira(i) = int(i/nvtk) + 1
				iba(i) = i - nvtk*int(i/nvtk)
				if (iba(i) .eq. 0) then
					ira(i) = ira(i) - 1
					iba(i) = nvtk
				end if
			end do
			do j=1,nny
				jra(j) = int(j/nvtk) + 1
				jba(j) = j - nvtk*int(j/nvtk)
				if (jba(j) .eq. 0) then
					jra(j) = jra(j) - 1
					jba(j) = nvtk
				end if
			end do

			do i=1,nnx
			   do j=1,nny
				igrid = nvtk*(iba(i) - 1) + jba(j)
				do ieq=1,nQ
				   qvtk(i,j,ieq) = sum(bfvtk(igrid,1:nbasis)*Qin(ira(i),jra(j),ieq,1:nbasis))
				   dxrh = sum(bfvtk_dx(igrid,1:nbasis)*Qin(ira(i),jra(j),ieq,1:nbasis))
				   dyrh = sum(bfvtk_dy(igrid,1:nbasis)*Qin(ira(i),jra(j),ieq,1:nbasis))
				   dxmy = sum(bfvtk_dx(igrid,1:nbasis)*Qin(ira(i),jra(j),my,1:nbasis))
				   dymx = sum(bfvtk_dy(igrid,1:nbasis)*Qin(ira(i),jra(j),mx,1:nbasis))
				end do
	      			qvtk_dxvy(i,j) = (qvtk(i,j,rh)*dxmy - qvtk(i,j,my)*dxrh)/qvtk(i,j,rh)**2
	      			qvtk_dyvx(i,j) = (qvtk(i,j,rh)*dymx - qvtk(i,j,mx)*dyrh)/qvtk(i,j,rh)**2
			   end do
			end do

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)
					var_xml_val_x(l)=log(abs(qvtk(i,j,rh))*n0)/log(10.)
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny, &
		                   varname = 'Log Density',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)
					var_xml_val_x(l)=qvtk(i,j,rh)*n0
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny, &
		                   varname = 'Density',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)
                                        dn = qvtk(i,j,rh)
					dni = 1./dn					
					vx = qvtk(i,j,mx)*dni
					vy = qvtk(i,j,my)*dni
					vz = qvtk(i,j,mz)*dni
                                        u2 = vx**2 + vy**2 + vz**2
                                        pe = qvtk(i,j,en)
					var_xml_val_x(l)=0.6667*(pe - 0.5*dn*u2)*p0
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny, &
		                   varname = 'Pressure',                    &
		                   var     = var_xml_val_x)

		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)
					var_xml_val_x(l)=v0*qvtk(i,j,mx)/qvtk(i,j,rh)
					var_xml_val_y(l)=v0*qvtk(i,j,my)/qvtk(i,j,rh)
					var_xml_val_z(l)=v0*qvtk(i,j,mz)/qvtk(i,j,rh)
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny, &
		                   varname = 'Velocity',                    &
		                   varX     = var_xml_val_x, varY     = var_xml_val_y, varZ     = var_xml_val_Z )


		do i=1+nb,nnx-nb
			do j=1+nb,nny-nb
					l=(i-nb)+(j-nb-1)*(nnx-2*nb)
					var_xml_val_x(l)=-2.*(qvtk_dxvy(i,j)-qvtk_dyvx(i,j))*dxi/t0
			enddo
		enddo
		E_IO = VTK_VAR_XML(NC_NN   = nnx*nny, &
		                   varname = 'Vorticity',                    &
		                   var     = var_xml_val_x)

		var_xml_val_x(:)=-8   
                
		                   
		E_IO = VTK_DAT_XML(var_location     = 'cell', &
		                   var_block_action = 'Close')
		E_IO = VTK_GEO_XML()
		E_IO = VTK_END_XML()

        end subroutine output_vtk

!--------------------------------------------------------------------------------

	subroutine writeQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow)

	implicit none
	real :: Qin(nx,ny,nQ,nbasis),tnow,dtnow
	integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,nump,numd,qq,j,i,ir
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

	do qq=1,nQ
    do ir=1,nbasv(qq)
			do j=1,ny
				write(3,*) (Qin(i,j,qq,ir),i=1,nx)
			enddo
	enddo
	enddo

	write(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow

	close(3)

	end subroutine writeQ

!------------------------------------------------------------------------------

	subroutine readQ(fprefix,irank,iddump,Qin,tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow)

	implicit none
	real :: Qin(nx,ny,nQ,nbasis),tnow,dtnow
	integer :: irank,iddump,noutnow,mpi_nxnow,mpi_nynow,nump,numd,qq,j,i,ir
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

	do qq=1,nQ
    do ir=1,nbasv(qq)
			do j=1,ny
				read(3,*) (Qin(i,j,qq,ir),i=1,nx)
			enddo
	enddo
	enddo

	read(3,*) tnow,dtnow,noutnow,mpi_nxnow,mpi_nynow

	close(3)

	end subroutine readQ

!--------------------------------------------------------------------------------
	
	subroutine set_bfvals
		! Defines local basis function values and weights for 1, 2, or 3-point Gaussian quadrature.
		! Basis functions are evaluated in cell interior and on cell edges.

	implicit none
	
    	call set_vtk_vals_2D()
   	call set_internal_vals_2D()     ! Define basis function values at quadrature points INTERNAL to cell. 
    	call set_edge_vals()   ! Define local basis function values at quadrature points on a cell edge.   	
	call set_weights()     ! Define weights for integral approximation using Gaussian quadrature.

	end subroutine set_bfvals

!----------------------------------------------------------
	subroutine set_vtk_vals_2D

	! Define basis function values on a rectangular 3D grid INTERNAL to the cell.  
	! For use in output_vtk().

	implicit none
	integer ixyz,i,igrid

	dxvtk = 1./(nvtk*dxi)
	dyvtk = 1./(nvtk*dyi)

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
		
		bfvtk(1:nvtk2,1) = 1.		! basis function = 1
		do i=1,nvtk
	           bfvtk((i-1)*nvtk+1:i*nvtk,kx) = xgrid(i)		! basis function = y	   		
	           bfvtk((i-1)*nvtk+1:i*nvtk,ky) = xgrid(1:nvtk)	! basis function = z
		end do	
		
		bfvtk(1:nvtk2,kxx) = 1.5*bfvtk(1:nvtk2,kx)**2 - 0.5    ! basis function = P_2(x)
		bfvtk(1:nvtk2,kxxx) = 2.5*bfvtk(1:nvtk2,kx)**3 - 1.5*bfvtk(1:nvtk2,kx)   ! basis function = P3(x)
		bfvtk(1:nvtk2,kyy) = 1.5*bfvtk(1:nvtk2,ky)**2 - 0.5    ! basis function = P_2(y)
		bfvtk(1:nvtk2,kyyy) = 2.5*bfvtk(1:nvtk2,ky)**3 - 1.5*bfvtk(1:nvtk2,ky)   ! basis function = P3(y)
	
	bfvtk(1:nvtk2,kxy) = bfvtk(1:nvtk2,kx)*bfvtk(1:nvtk2,ky)		! basis function = xy
	bfvtk(1:nvtk2,kxxy) = bfvtk(1:nvtk2,kxx)*bfvtk(1:nvtk2,ky)     ! basis function = P2(x)y
	bfvtk(1:nvtk2,kxyy) = bfvtk(1:nvtk2,kx)*bfvtk(1:nvtk2,kyy)     ! basis function = P2(y)x
	bfvtk(1:nvtk2,kxxyy) = bfvtk(1:nvtk2,kxx)*bfvtk(1:nvtk2,kyy)     ! basis function = P2(x)P_2(y)

	   do igrid=1,nvtk2
		   bfvtk_dx(igrid,1) = 0.
		   bfvtk_dx(igrid,kx) = 1.	
		   bfvtk_dx(igrid,ky) = 0.
		   bfvtk_dx(igrid,kxx) = 3.*bfvtk(igrid,kx)
		   bfvtk_dx(igrid,kyy) = 0.
		   bfvtk_dx(igrid,kxy) = bfvtk(igrid,ky)
		   bfvtk_dx(igrid,kxxx) = 7.5*bfvtk(igrid,kx)**2 - 1.5
		   bfvtk_dx(igrid,kyyy) = 0.

		   bfvtk_dx(igrid,kxxy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,ky)
		   bfvtk_dx(igrid,kxyy) = bfvtk(igrid,kyy)
		   bfvtk_dx(igrid,kxxyy) = 3.*bfvtk(igrid,kx)*bfvtk(igrid,kyy)

		   bfvtk_dy(igrid,1) = 0.
		   bfvtk_dy(igrid,ky) = 1.	
		   bfvtk_dy(igrid,kx) = 0.
		   bfvtk_dy(igrid,kyy) = 3.*bfvtk(igrid,ky)
		   bfvtk_dy(igrid,kxx) = 0.
		   bfvtk_dy(igrid,kxy) = bfvtk(igrid,kx)
		   bfvtk_dy(igrid,kyyy) = 7.5*bfvtk(igrid,ky)**2 - 1.5
		   bfvtk_dy(igrid,kxxx) = 0.
		   bfvtk_dy(igrid,kxxy) = bfvtk(igrid,kxx)
		   bfvtk_dy(igrid,kxyy) = 3.*bfvtk(igrid,ky)*bfvtk(igrid,kx)
		   bfvtk_dy(igrid,kxxyy) = 3.*bfvtk(igrid,kxx)*bfvtk(igrid,ky)	   
	   end do


	end subroutine set_vtk_vals_2D

!----------------------------------------------------------
	subroutine set_internal_vals_2D	

	! Define basis function values at quadrature points INTERNAL to cell.

	implicit none
	integer ixyz, i
	real c15d5,c1dsq3,xq4p,xq4m,xquad(20)

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
	end if	
	

	if (iquad .gt. 1) then
	     bfvals_int(1:npg,1) = 1.		! basis function = 1
	     do i=1,nedge
	       bfvals_int((i-1)*nedge+1:i*nedge,ky) = xquad(1:nedge)	! basis function = y
	       bfvals_int((i-1)*nedge+1:i*nedge,kx) = xquad(i)		! basis function = x	   		
	     end do	
	end if

		
	bfvals_int(1:npg,kxx) = 1.5*bfvals_int(1:npg,kx)**2 - 0.5    ! basis function = P_2(x)
	bfvals_int(1:npg,kxxx) = 2.5*bfvals_int(1:npg,kx)**3 - 1.5*bfvals_int(1:npg,kx)   ! basis function = P_3(x)
	bfvals_int(1:npg,kyy) = 1.5*bfvals_int(1:npg,ky)**2 - 0.5    ! basis function = P_2(y)
	bfvals_int(1:npg,kyyy) = 2.5*bfvals_int(1:npg,ky)**3 - 1.5*bfvals_int(1:npg,ky)   ! basis function = P_3(y)

	

	bfvals_int(1:npg,kxy) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,ky)		! basis function = xy
	bfvals_int(1:npg,kxxy) = bfvals_int(1:npg,kxx)*bfvals_int(1:npg,ky)     ! basis function = P_2(x)y
	bfvals_int(1:npg,kxyy) = bfvals_int(1:npg,kx)*bfvals_int(1:npg,kyy)     ! basis function = P_2(y)x
	bfvals_int(1:npg,kxxyy) = bfvals_int(1:npg,kxx)*bfvals_int(1:npg,kyy)     ! basis function = P_2(x)P_2(y)

	end subroutine set_internal_vals_2D	

!----------------------------------------------------------
	subroutine set_edge_vals

	! Define local basis function values at quadrature points on a cell edge.
	! Used in flux_cal() and prepare_exchange() for interpolating from cell center onto cell edge.
	! Also used in glflux().

	implicit none
	integer ixyz,i
	real c15d5,c1dsq3,xquad(20),xq4p,xq4m

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


	! bfvals_rp:  positive r-edge
	! bfvals_rm:  negative r-edge
	! bfvals_rp(,1):  value=1 on positive r-edge
	! bfvals_rp(,kx):  x-values on positive r-edge 
	! bfvals_rp(,ky):  y-values on positive r-edge
	! bfvals_rp(,kxx):  P_2(x)-values on positive r-edge	
	! bfvals_rp(,kyy):  P_2(y)-values on positive r-edge
	! bfvals_rp(,kxy):  xy-values on positive r-edge

	bfvals_xp(1:nedge,1) = 1.
	bfvals_yp(1:nedge,1) = 1.
	if (iquad .eq. 1) then
		bfvals_yp(1,kx) = 0.
		bfvals_yp(1,ky) = 0.
	end if	

	if (iquad .gt. 1) then
	do i=1,nedge
	    bfvals_yp(1:nedge,kx) = xquad(1:nedge)  
	    bfvals_yp(1:nedge,ky) = 1.
	    bfvals_xp(1:nedge,ky) = xquad(1:nedge) 
	    bfvals_xp(1:nedge,kx) = 1. 
	end do
	end if	


	bfvals_xp(1:nedge,kxx) = 1.5*bfvals_xp(1:nedge,kx)**2 - 0.5
	bfvals_xp(1:nedge,kyy) = 1.5*bfvals_xp(1:nedge,ky)**2 - 0.5
	bfvals_xp(1:nedge,kxy) = bfvals_xp(1:nedge,kx)*bfvals_xp(1:nedge,ky)     ! basis function = xy

	bfvals_yp(1:nedge,kxx) = 1.5*bfvals_yp(1:nedge,kx)**2 - 0.5
	bfvals_yp(1:nedge,kyy) = 1.5*bfvals_yp(1:nedge,ky)**2 - 0.5
	bfvals_yp(1:nedge,kxy) = bfvals_yp(1:nedge,kx)*bfvals_yp(1:nedge,ky)     ! basis function = xy

	bfvals_xp(1:nedge,kxxx) = 2.5*bfvals_xp(1:nedge,kx)**3 - 1.5*bfvals_xp(1:nedge,kx)
	bfvals_xp(1:nedge,kyyy) = 2.5*bfvals_xp(1:nedge,ky)**3 - 1.5*bfvals_xp(1:nedge,ky)
	bfvals_xp(1:nedge,kxxy) = bfvals_xp(1:nedge,kxx)*bfvals_xp(1:nedge,ky)     ! basis function = P_2(x)y
	bfvals_xp(1:nedge,kxyy) = bfvals_xp(1:nedge,kx)*bfvals_xp(1:nedge,kyy)     ! basis function = P_2(y)x
	bfvals_xp(1:nedge,kxxyy) = bfvals_xp(1:nedge,kxx)*bfvals_xp(1:nedge,kyy)     ! basis function = P_2(x)P_2(y)


	bfvals_yp(1:nedge,kxxx) = 2.5*bfvals_yp(1:nedge,kx)**3 - 1.5*bfvals_yp(1:nedge,kx)
	bfvals_yp(1:nedge,kyyy) = 2.5*bfvals_yp(1:nedge,ky)**3 - 1.5*bfvals_yp(1:nedge,ky)
	bfvals_yp(1:nedge,kxxy) = bfvals_yp(1:nedge,kxx)*bfvals_yp(1:nedge,ky)     ! basis function = P_2(x)y
	bfvals_yp(1:nedge,kxyy) = bfvals_yp(1:nedge,kx)*bfvals_yp(1:nedge,kyy)     ! basis function = P_2(y)x
	bfvals_yp(1:nedge,kxxyy) = bfvals_yp(1:nedge,kxx)*bfvals_yp(1:nedge,kyy)     ! basis function = P_2(x)P_2(y)


	bfvals_xm = bfvals_xp
	bfvals_xm(1:nedge,kx) = -bfvals_xp(1:nedge,kx)
	bfvals_xm(1:nedge,kxy) = -bfvals_xp(1:nedge,kxy)
	bfvals_xm(1:nedge,kxyy) = -bfvals_xp(1:nedge,kxyy)
	bfvals_xm(1:nedge,kxxx) = -bfvals_xp(1:nedge,kxxx)	

	
	bfvals_ym = bfvals_yp
	bfvals_ym(1:nedge,ky) = -bfvals_yp(1:nedge,ky)
	bfvals_ym(1:nedge,kxy) = -bfvals_yp(1:nedge,kxy)
	bfvals_ym(1:nedge,kxxy) = -bfvals_yp(1:nedge,kxxy)
	bfvals_ym(1:nedge,kyyy) = -bfvals_yp(1:nedge,kyyy)

	
	! Organize local basis values on faces into 1-D vectors.
	! Used in limiter() and max_lim().

	bf_edges(1:nslim,1) = 1.		! basis function = 1
	
	do ixyz=kx,ky
		bf_edges(1:nedge,ixyz) = bfvals_xm(1:nedge,ixyz)		! basis function = x,y
		bf_edges(nedge+1:2*nedge,ixyz) = bfvals_xp(1:nedge,ixyz)		! basis function = x,y
		bf_edges(2*nedge+1:3*nedge,ixyz) = bfvals_ym(1:nedge,ixyz)		! basis function = x,y
		bf_edges(3*nedge+1:4*nedge,ixyz) = bfvals_yp(1:nedge,ixyz)		! basis function = x,y
		bf_edges(4*nedge+1:nslim,ixyz) = bfvals_int(1:npg,ixyz)		! basis function = x,y
	end do

	bf_edges(1:nslim,kxy) = bf_edges(1:nslim,kx)*bf_edges(1:nslim,ky)     ! basis function = xy

	
	do i=0,1		
		bf_edges(1:nslim,kxx+i) = 1.5*bf_edges(1:nslim,kx+i)**2 - 0.5    ! basis function = P_2(s)
		bf_edges(1:nslim,kxxx+i) = 2.5*bf_edges(1:nslim,kx+i)**3 - 1.5*bf_edges(1:nslim,kx+i)   
										! basis function = P_3(s)
	end do

	bf_edges(1:nslim,kxxy) = bf_edges(1:nslim,kxx)*bf_edges(1:nslim,ky)     ! basis function = P_2(x)y
	bf_edges(1:nslim,kxyy) = bf_edges(1:nslim,kx)*bf_edges(1:nslim,kyy)     ! basis function = P_2(y)x
	bf_edges(1:nslim,kxxyy) = bf_edges(1:nslim,kxx)*bf_edges(1:nslim,kyy)     ! basis function = P_2(x)P_2(y)

	end subroutine set_edge_vals

!----------------------------------------------------------
	subroutine set_weights

	implicit none
	integer i,nedge2
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

	 nedge2 = nedge*nedge

	 if (iquad .eq. 1) then		! 1-point quadrature
	    wgt3d(1) = 8.
	 end if
	 if (iquad .eq. 2) then		! 2-point quadrature
	    wgt3d(1:8) = 1.
	 end if

	if (iquad .ge. 3) then
	   do i= 1,nedge
		  wgt3d((i-1)*nedge2+1:i*nedge2) = wgt2d(1:nedge2)*wgt1d(i)
	   end do
	end if


	end subroutine set_weights


!--------------------------------------------------------------------------------

end program


