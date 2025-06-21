
module globaladj
 
 real*8,parameter:: big=1.0d8
 integer:: imin(4),imax(4),istep(4),jmin(4),jmax(4),jstep(4)   
 real*8,allocatable:: T(:,:),u(:,:)
 integer,allocatable:: isw(:,:)
 
end module
 
subroutine subadjoint_rec(k,Tt,residuals_time_rec,residuals_dou_rec,adj)
 
 use globalp
 use globaladj

 integer:: i,j,k,l,ite
 real*8:: Tt(nx,ny),adj(nx,ny),residuals_time_rec(nr),residuals_dou_rec(nr,nr)
 real*8:: cvg_err,eps,maxu
 real*8,allocatable:: u0(:,:)

 allocate(T(-2:nx+2,-2:ny+3),u(-2:nx+3,-2:ny+3),u0(-2:nx+3,-2:ny+3),isw(nx,ny))
 
 T=big 
 do i=1, nx
 do j=1, ny
    T(i,j)=Tt(i,j)
 enddo
 enddo
 
 call init_sweep2_rec
 
 call init_receiver_rec(k,residuals_time_rec,residuals_dou_rec)

 eps=1.0d-12
 ite=0   
 cvg_err=100000

 do while((cvg_err.gt.eps).and.(ite.le.100))
    u0=u
    call Tadjoint_sweep_rec(k)
    ite=ite+1
    cvg_err=sum(abs(u(1:nx,1:ny)-u0(1:nx,1:ny)))/dble((nx-1)*(ny-1))   
 enddo
    
do j=1, ny
do i=1, nx
   adj(i,j)=u(i,j)
enddo
enddo
 
 deallocate(T,u,u0,isw)

end subroutine subadjoint_rec                               

subroutine init_sweep2_rec
 use globalp
 use globaladj
 
 imin(1)=1
 imax(1)=nx
 istep(1)=1
 jmin(1)=1
 jmax(1)=ny
 jstep(1)=1

 imin(2)=1
 imax(2)=nx
 istep(2)=1
 jmin(2)=ny
 jmax(2)=1
 jstep(2)=-1

 imin(3)=nx
 imax(3)=1
 istep(3)=-1
 jmin(3)=ny
 jmax(3)=1
 jstep(3)=-1

 imin(4)=nx
 imax(4)=1
 istep(4)=-1
 jmin(4)=1
 jmax(4)=ny
 jstep(4)=1
 
end subroutine        

subroutine init_receiver_rec(ishot,residuals_time_rec,residuals_dou_rec)
 use globalp
 use globaladj
 integer:: i,j,k,m,ishot                                        
 real*8:: Tq,Tr,par1,par2
 real*8:: nnq,nnr
 real*8::residuals_time_rec(nr),residuals_dou_rec(nr,nr)
 
 u=0
 isw=0

 do k=1, recnum(ishot)
    i=reci(ishot,k)
    j=recj(ishot,k)

    isw(i,j)=1

    if(j.eq.ny)then
    Tr=(T(i,ny-2)-4*T(i,ny-1)+3*T(i,ny))/(2*dy)
   elseif(j.eq.1)then
    Tr=(-T(i,3)+4*T(i,2)-3*T(i,1))/(2*dy)
   else
     Tr=(T(i,j+1)-T(i,j-1))/(2*dy)
   endif

   if(i.eq.nx)then
     Tq=(T(nx-2,j)-4*T(nx-1,j)+3*T(nx,j))/(2*dx)
   elseif(i.eq.1)then
     Tq=(-T(3,j)+4*T(2,j)-3*T(1,j))/(2*dx)
   else
     Tq=(T(i+1,j)-T(i-1,j))/(2*dx)
   endif
    
    nnr=1
    nnq=0
    
    par1=0
    do m=1, recnum(ishot)
       if(m.ne.k)then
          par1=par1-residuals_dou_rec(k,m)/(recnum(ishot)-1)
       endif
    enddo
    
    par1=par1+residuals_time_rec(k)
    par2=nnq*Tq+nnr*Tr
    
    if(abs(par1).gt.(1.0d-10).and.abs(par2).gt.(1.0d-6))then
        u(i,j)=par1/par2
    endif

 enddo

end subroutine

subroutine Tadjoint_sweep_rec(l)
 use globalp
 use globaladj
 integer:: i,j,k,l
 real*8:: alphaw,alphae,betas,betan
 real*8:: ac,aw,ae,bc,bs,bn,cc,cw,ce,cs,cn
 real*8:: par1,par2,par3,par4,par5
 real*8,external:: pos,neg
           
 do k=1, 4                     
    do j=jmin(k), jmax(k), jstep(k)
    do i=imin(k), imax(k), istep(k) 

       if(isw(i,j).eq.0)then
    
         if(T(i,j).ge.big.or.T(i-1,j).ge.big.or.T(i+1,j).ge.big.or.T(i,j+1).ge.big.or.T(i,j-1).ge.big)then
            u(i,j)=0
         elseif(i.eq.sorci(l).and.j.eq.sorcj(l))then
            u(i,j)=0
         else          

           alphaw=-(T(i,j)-T(i-1,j))/dx
           alphae=-(T(i+1,j)-T(i,j))/dx
           betas=-(T(i,j)-T(i,j-1))/dy
           betan=-(T(i,j+1)-T(i,j))/dy

           par1=pos(alphaw)*u(i-1,j)
           par2=neg(alphae)*u(i+1,j)
           par3=pos(betas)*u(i,j-1)
           par4=neg(betan)*u(i,j+1)
           par5=(pos(alphae)-neg(alphaw))/dx+(pos(betan)-neg(betas))/dy   

           if(abs(par5).gt.1.0d-9)then
             u(i,j)=((par1-par2)/dx+(par3-par4)/dy)/par5
             if(abs(u(i,j)).lt.1.0d-9)u(i,j)=0
           else
             u(i,j)=0
           endif

         endif
       endif        
    enddo
    enddo
 enddo

end subroutine

!********************************************************************
real*8 function pos(par)
 real*8:: par

 pos=(par+abs(par))/2.0

end function

!********************************************************************
real*8 function neg(par)
 real*8:: par

 neg=(par-abs(par))/2.0

end function

