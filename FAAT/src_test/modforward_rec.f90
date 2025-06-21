 
 module globalfwd_rec 
 real*8,parameter:: big=1.0d10
 integer:: imin(4),imax(4),istep(4),jmin(4),jmax(4),jstep(4)   
 real*8,allocatable:: tau(:,:),tau0(:,:),t0(:,:),t0dx(:,:),t0dy(:,:)
 integer,allocatable:: isw(:,:)
 real*8::ssa,ssb,ssc
 real*8,external:: t_anad2d_rec
 
end module
 
subroutine subforward_rec(k,vel,u)

 use globalp
 use globalfwd_rec
 integer:: i,j,k,ite
 real*8:: vel(nx,ny),u(nx,ny)
 real*8:: cvg_err,eps,tc  

 allocate(tau(-2:nx+3,-2:ny+3),tau0(-2:nx+3,-2:ny+3),isw(nx,ny))
 allocate(t0(-2:nx+3,-2:ny+3),t0dx(-2:nx+3,-2:ny+3),t0dy(-2:nx+3,-2:ny+3))
 
 call init_sweep_rec
 
 call init_source_rec(k,vel)                                                
 
 eps=1.0d-12
 ite=0   
 cvg_err=100000
 
 do while((cvg_err.gt.eps).and.(ite.le.20))
    tau0=tau    
    call factored_sweep_rec(vel)                   
    ite=ite+1                                                                                         
    cvg_err=sum(abs(tau(1:nx,1:ny)-tau0(1:nx,1:ny)))/dble((nx-1)*(ny-1))                               
 enddo                                                                                     
 
 u=big
 do j=1, ny
 do i=1, nx
    u(i,j)=t0(i,j)*tau(i,j)
 enddo
 enddo

 deallocate(tau,tau0,t0,t0dx,t0dy,isw)

end subroutine                                                  

subroutine init_sweep_rec
 use globalp
 use globalfwd_rec
 
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
 
subroutine init_source_rec(k,vel)
 use globalp
 use globalfwd_rec

 integer:: i,j,k,isx,isy
 real*8:: vel(nx,ny)

 isx=sorci(k)
 isy=sorcj(k)

 ssa=vel(isx,isy)**2                       
 ssb=vel(isx,isy)**2
 ssc=0.
 
 do i=-2, nx+3
 do j=-2, ny+3
    t0(i,j)=t_anad2d_rec(t0dx(i,j),t0dy(i,j),i,j,dx,dy,isx,isy,ssa,ssb,ssc)    
 enddo
 enddo
 
 isw=0      
 tau=big                                    
 isw(isx,isy)=1                                    
 tau(isx,isy)=1                                                        

end subroutine
                                                           
 subroutine factored_sweep_rec(vel)
 use globalp
 use globalfwd_rec  
 real*8:: vel(nx,ny)
 real*8:: tauxmin,tauymin,taubar,slown,sa,sb,sc,h              
 integer:: i,j,k,m,sx(4),sy(4)                   
 real*8:: txc,tyc,t0c,aa,bb,cc                                  
 real*8:: taux(4),tauy(4),t0x(4),t0y(4),ubar(4)    
 real*8:: tauc1,tauc2,taub1,taub2,cau11,cau12,caub1,caub2,cau21,cau22  
 
 h=dx

 do k=1, 4                                                      
  
    do j=jmin(k),jmax(k),jstep(k)                        
    do i=imin(k),imax(k),istep(k)                
       if(isw(i,j).eq.0)then                                 

   sa = 1                     
   sb = 1
   slown = 1.0/vel(i,j)**2

   t0c = t0(i,j)                     
   txc = t0dx(i,j)
   tyc = t0dy(i,j)

   taux(1) = tau(i-1,j)
   tauy(1) = tau(i,j-1)
   taux(2) = tau(i-1,j)
   tauy(2) = tau(i,j+1)
   taux(3) = tau(i+1,j)
   tauy(3) = tau(i,j+1)
   taux(4) = tau(i+1,j)
   tauy(4) = tau(i,j-1)                                 

   t0x(1) = t0(i-1,j)  
   t0y(1) = t0(i,j-1)  
   t0x(2) = t0(i-1,j)  
   t0y(2) = t0(i,j+1)  
   t0x(3) = t0(i+1,j)  
   t0y(3) = t0(i,j+1) 
   t0x(4) = t0(i+1,j) 
   t0y(4) = t0(i,j-1)         

   sx = (/ 1,  1, -1, -1 /)
   sy = (/ 1, -1, -1,  1 /)

   do m = 1, 4
      aa = sa*txc**2 + sb*tyc**2 + (sa + sb)*(t0c/h)**2 + 2*t0c/h*(txc*sa*sx(m) + tyc*sb*sy(m))

      bb = 2*t0c/h*(sx(m)*(-sa*txc)*taux(m) + sy(m)*(-sb*tyc)*tauy(m)) - 2*(sa*taux(m) + sb*tauy(m))*(t0c/h)**2

      cc = t0c**2 / (h*h) * (sa*taux(m)**2 + sb*tauy(m)**2) - slown**2

      disc = bb**2 - 4*aa*cc
      if(disc .ge. 0)then
         tauc1 = (-bb + sqrt(disc)) / (2*aa)
         tauc2 = (-bb - sqrt(disc)) / (2*aa)

         cau11 = sa * (sx(m)*(tauc1 - taux(m))/h * t0c + tauc1*txc)
         cau12 = sb * (sy(m)*(tauc1 - tauy(m))/h * t0c + tauc1*tyc)
         cau21 = sa * (sx(m)*(tauc2 - taux(m))/h * t0c + tauc2*txc)
         cau22 = sb * (sy(m)*(tauc2 - tauy(m))/h * t0c + tauc2*tyc)

         if(sx(m)*cau11.ge.0.and.sy(m)*cau12.ge.0.and.sx(m)*cau21.ge.0.and.sy(m)*cau22.ge.0)then
            ubar(m) = min(tauc1*t0c, tauc2*t0c)
         elseif(sx(m)*cau11.ge.0.and.sy(m)*cau12.ge.0)then
            ubar(m) = tauc1*t0c
         elseif(sx(m)*cau21.ge.0.and.sy(m)*cau22.ge.0)then
            ubar(m) = tauc2*t0c
         else                                                                        
            taub1 = (t0c*taux(m) + h*sqrt(sb / (sa*sb))) / (t0c + txc*h*sx(m))                             
            taub2 = (t0c*tauy(m) + h*sqrt(sa / (sa*sb))) / (t0c + tyc*h*sy(m))
            caub1 = taub1*t0c - taux(m)*t0x(m)
            caub2 = taub2*t0c - tauy(m)*t0y(m)
            if(caub1.ge.0.and.caub2.ge.0)then
               ubar(m) = min(taub1*t0c, taub2*t0c)
            elseif(caub1.ge.0.and.caub2.lt.0)then
               ubar(m) = taub1*t0c
            elseif(caub1.lt.0.and.caub2.ge.0)then
               ubar(m) = taub2*t0c                         
            endif
         endif
      else                                         
         taub1 = (t0c*taux(m) + h*sqrt(sb / (sa*sb))) / (t0c + txc*h*sx(m))                             
         taub2 = (t0c*tauy(m) + h*sqrt(sa / (sa*sb))) / (t0c + tyc*h*sy(m))
         caub1 = taub1*t0c - taux(m)*t0x(m)
         caub2 = taub2*t0c - tauy(m)*t0y(m)
         if(caub1.ge.0.and.caub2.ge.0)then
            ubar(m) = min(taub1*t0c, taub2*t0c)
         elseif(caub1.ge.0.and.caub2.lt.0)then
            ubar(m) = taub1*t0c
         elseif(caub1.lt.0.and.caub2.ge.0)then
            ubar(m) = taub2*t0c                         
         endif
      endif 
   enddo

   umin = min(ubar(1), ubar(2), ubar(3), ubar(4)) 
   taubar = umin / t0c
   tau(i,j) = min(taubar, tau(i,j)) 
endif
    enddo
    enddo                                                          
 enddo

end subroutine 

 real*8 function t_anad2d_rec(txc,tyc,i,j,dx,dy,xsa,ysa,aa,bb,cc)
 integer::i,j,xsa,ysa
 real*8:: dx,dy,vzero
 real*8:: d0,aa,bb,cc
 real*8:: txc,tyc
 
 d0=bb*((i-xsa)*dx)**2+2*cc*(i-xsa)*dx*(j-ysa)*dy+aa*((j-ysa)*dy)**2

 t_anad2d_rec=sqrt(d0/(aa*bb-cc**2))

 if(d0.gt.0)then
   txc=(bb*(i-xsa)*dx+cc*(j-ysa)*dy)/sqrt((aa*bb-cc**2)*d0)
   tyc=(cc*(i-xsa)*dx+aa*(j-ysa)*dy)/sqrt((aa*bb-cc**2)*d0)
 else    
   txc=0.d0
   tyc=0.d0
 endif

end function t_anad2d_rec
