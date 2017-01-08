module io_mod

use lib_vtk_io

contains


subroutine output_vtk0(Qin,nout,iam)

implicit none
real Qin(nx,ny,nz,nQ,nbasis)
integer nout
integer(I4P) , parameter :: nb=1*ngu,sprd=0*1
!    integer(I4P), parameter:: nnx=nx-2*nb,nny=ny-2*nb,nnz=nz-2*nb
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
out_name='data/data3/perseus_p'//pname//'_t'//tname//'.vtr'
! out_name='data/data4/perseus_p'//pname//'_t'//tname//'.vtr'
print *, out_name
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

            ! from "viscosity" version
            P = (aindex - 1)*(Qin(i,j,k,en,1) - 0.5*Qin(i,j,k,rh,1)*(vx**2 + vy**2 + vz**2))
            ! if(ieos .eq. 1)P = (aindex - 1)*(Qin(i,j,k,en,1) - 0.5*Qin(i,j,k,rh,1)*(vx**2 + vy**2 + vz**2))
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
real P, vx, vy, vz, dni, dxrh,dyrh,dxmy,dymx
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

! "outdir" is a global variable specifying the output directory
out_name=''//outdir//'/perseus_p'//pname//'_t'//tname//'.vtr'

if (iam .eq. print_mpi) print *, out_name
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
            P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))
            var_xml_val_x(l)=P*te0*dni
        enddo
    enddo
enddo
E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                   varname = 'Temperature',                    &
                   var     = var_xml_val_x)

do i=1+nb,nnx-nb
    do j=1+nb,nny-nb
        do k=1+nb,nnz-nb
            l = (i-nb)+(j-nb-1)*(nnx-2*nb)+(k-nb-1)*(nnx-2*nb)*(nny-2*nb)
            dni = 1./qvtk(i,j,k,rh)
            vx = qvtk(i,j,k,mx)*dni
            vy = qvtk(i,j,k,my)*dni
            vz = qvtk(i,j,k,mz)*dni
            ! from "viscosity" version
            P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))
            var_xml_val_x(l)=P
            ! if(ieos .eq. 1) P = (aindex - 1.)*(qvtk(i,j,k,en) - 0.5*qvtk(i,j,k,rh)*(vx**2 + vy**2 + vz**2))
            ! if(ieos .eq. 2)then
            !     P = P_1*(qvtk(i,j,k,rh)**7.2 - 1.) + P_base
            !     if(P < P_floor) P= P_floor
            ! end if
            ! var_xml_val_x(l) = P*P0
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

! from "viscosity" version
do i=1+nb,nnx-nb
    do j=1+nb,nny-nb
        do k=1+nb,nnz-nb
            l=(i-nb)+(j-nb-1)*(nx-2*nb)+(k-nb-1)*(nx-2*nb)*(ny-2*nb)
                var_xml_val_x(l)= -2*(Q_r2(i,j,k,my,2)-Q_r2(i,j,k,mx,3))/Q_r2(i,j,k,rh,1)*dxi/t0
        enddo
    enddo
enddo
E_IO = VTK_VAR_XML(NC_NN   = nnx*nny*nnz, &
                  varname = 'Vorticity2',                    &
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


end module io_mod
