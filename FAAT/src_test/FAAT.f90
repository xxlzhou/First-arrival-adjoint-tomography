!****************************************************************************
! First-arrival adjoint tomography 
! Author: Xiaole Zhou
! Xiaole Zhou, Yu Zhang, Shijie Hao et al. 2025. Near-surface site characterization 
! at Bishan-AMK Park, Singapore, using first-arrival adjoint tomography
! via the joint inversion of absolute and differential traveltimes. GJI
!****************************************************************************
module globalp

 real*8::                dx,dy,dsou,drecx,W_t,W_pr,W_ps
 integer::             nx,ny,ns,nr,nop,nw1,nw2
 real*8,allocatable::    x(:,:),y(:,:),pa(:,:),pb(:,:),pc(:,:)
 integer,allocatable:: sorci(:),sorcj(:),recnum(:),reci(:,:),recj(:,:)
 integer,allocatable:: sorc1i(:),sorc1j(:),recnum1(:),rec1i(:,:),rec1j(:,:)

end module
 
program FTDET2d
 use globalp
 include 'optim_type.h'
 type(optim_type)::   optim
 character(len=4)::   FLAG
 character(len=160):: gridgeo,gridveli,sorcfile,recfile,obsfile,obs1file,acqui,acqui1,gridvel1
 integer::            i,j,k,igrad,niter,optimalgo,nitermax
 real::               perc,alpha,cost 
 real,allocatable::   vel1D(:),grad(:),pgrad(:)
 real*8,allocatable::   veli(:,:),vel(:,:),vel1(:,:)
 real*8,allocatable::   Trobs(:,:),douobs(:,:,:)
 integer,parameter::  lbytes=4,ibound=0
 real*8,parameter::     pi=3.141592654,thresconv=1e-8,ngradlbfgs=5
 logical,parameter::  debug=.false.

 open(unit=11,file='input.in',status='old') 
 read(11,*)nx, ny, dx, dy
 read(11,*)W_t,W_pr
 read(11,*)gridveli
 read(11,*)gridvel1
 read(11,*)acqui
 read(11,*)obsfile
 read(11,*)obs1file
 read(11,*)nitermax
 read(11,*)perc
 read(11,*)optimalgo
 read(11,*)nw1
 read(11,*)nw2
 close(11)
 
 nop=nx*ny
 
 allocate(x(nx,ny),y(nx,ny))
 allocate(pa(nx,ny),pb(nx,ny),pc(nx,ny))
 allocate(veli(nx,ny),vel(nx,ny),vel1(nx,ny))
 allocate(vel1D(nop))
 allocate(grad(nop))
 allocate(pgrad(nop))
 
 open(unit=11,file=gridveli,status='old')
 do j=1, ny
 do i=1, nx
    read(11,*)x(i,j),y(i,j),veli(i,j)
 enddo
 enddo
 close(11)

 open(unit=11,file=gridvel1,status='old')
 do j=1, ny
 do i=1, nx
    read(11,*)x(i,j),y(i,j),vel1(i,j)
 enddo
 enddo
 close(11)

  open(unit=11,file=acqui,status='old')
 read(11,*)ns
 read(11,*)nr
 
 allocate(sorci(1:ns),sorcj(1:ns),recnum(1:ns))
 allocate(reci(1:ns,1:nr),recj(1:ns,1:nr))
 
 do i=1, ns
   read(11,*)
   read(11,*)sorci(i),sorcj(i)
   read(11,*)recnum(i)
   do k=1, recnum(i)
      read(11,*)reci(i,k),recj(i,k)
   enddo
 enddo
 close(11)
 
 allocate(Trobs(ns,nr))
 Trobs=0
 
 open(unit=14,file=obsfile,status='old')
 do i=1, ns
    read(14,*)
    do j=1, recnum(i)
       read(14,*)Trobs(i,j)
    enddo
 enddo
 close(14) 

 allocate(douobs(ns,nr,nr))
 douobs=0
 
 open(unit=16,file=obs1file,status='old')
 do i=1, ns
    read(16,*)
    do j=1, recnum(i)
    do k=1,recnum(i)
        read(16,*)douobs(i,j,k)
    enddo
    enddo
 enddo
 close(16)

   call subgradient(Trobs,douobs,veli,vel1D,grad,cost)
   
 if (debug) then
    open(15,file='vel0.bin',access='direct',recl=nop*lbytes)
    write(15,rec=1) vel1D
    close(15)
 endif
 
 call scalegradient(grad,cost,vel1D,perc,alpha)
 
 optim%niter_max = nitermax
 optim%conv = thresconv
 optim%print_flag = 1
 optim%bound = ibound
 optim%l = ngradlbfgs
 optim%debug = .true.
 
 FLAG='INIT'       
 niter = 0
 igrad=0
 
 do while ( (FLAG .ne. 'CONV') .and. (FLAG .ne. 'FAIL'))  
          
          call LBFGS(nop,vel1D,cost,grad,optim,FLAG)
         
       if (FLAG.EQ.'GRAD') then
          igrad = igrad + 1
 
          call subm2v(vel1D,vel)
          
            call subgradient(Trobs,douobs,vel,vel1D,grad,cost)
            
          grad(:)=grad(:)*alpha
          cost=cost*alpha

       else if (FLAG.eq.'NSTE' .OR. FLAG.eq.'CONV') then
          niter =  niter + 1 
          
          call subm2v(vel1D,vel)

       endif

 enddo
 
  open(unit=11,file='vel_anomaly.dat')
 do j=1, ny
 do i=1, nx
     write(11,*)x(i,j),y(i,j),vel(i,j)-vel1(i,j)
 enddo
 enddo
 close(11)
 
  open(unit=11,file='vel.dat')
 do j=1, ny
 do i=1, nx
     write(11,*)x(i,j),y(i,j),vel(i,j)
 enddo
 enddo
 close(11)
 
 deallocate(x,y)
 deallocate(pa,pb,pc)
 deallocate(veli,vel,vel1)
 deallocate(vel1D,grad,pgrad)
 deallocate(sorci,sorcj,reci,recj)
 deallocate(Trobs,douobs)
 
end program FTDET2d
    
