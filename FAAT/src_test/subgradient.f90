
subroutine subgradient(Trobs,douobs,vel,vel1D,grad,cost)
 
 use globalp
 integer:: i,j,k,l,adjn
 real::vel1D(nop),grad(nop),pgrad(nop)
 real::cost
 real*8:: Trobs(ns,nr),douobs(ns,nr,nr),vel(nx,ny),damp,hess(nop)
 real*8,allocatable:: T_rec(:,:),adj_rec(:,:,:),padj_rec(:,:,:)
 real*8,allocatable:: Trcal(:,:),doucal(:,:,:),residuals_time_rec(:,:),residuals_dou_rec(:,:,:),res_deri_slope_rec(:,:)
 real*8,allocatable:: gradt_rec(:,:,:),hesst_rec(:,:,:),costt_rec(:)
 real*8,allocatable::gradtt_rec(:,:),hesstt_rec(:,:)
 real*8,allocatable:: gradtt(:,:),hesstt(:,:),sgradtt(:,:),pgradtt(:,:),spgradtt(:,:)
 
 
 allocate(T_rec(nx,ny),adj_rec(ns,nx,ny),padj_rec(ns,nx,ny))
 allocate(Trcal(ns,nr),doucal(ns,nr,nr),residuals_time_rec(ns,nr),residuals_dou_rec(ns,nr,nr),res_deri_slope_rec(ns,nr))
 allocate(gradt_rec(ns,nx,ny),hesst_rec(ns,nx,ny),costt_rec(ns))
 allocate(gradtt_rec(nx,ny),hesstt_rec(nx,ny))
 allocate(gradtt(nx,ny),hesstt(nx,ny),sgradtt(nx,ny),pgradtt(nx,ny),spgradtt(nx,ny))
 
 call subv2m(vel,vel1D)
 
 gradt_rec=0
 
 do k=1, ns
    
    call subforward_rec(k,vel,T_rec)
    call subdata_rec(k,T_rec,Trcal(k,:),doucal(k,:,:))
    call subresiduals_rec(k,Trcal(k,:),Trobs(k,:),doucal(k,:,:),douobs(k,:,:),&
        &residuals_time_rec(k,:),residuals_dou_rec(k,:,:))  
    call submisfit_rec(k,Trcal(k,:),Trobs(k,:),doucal(k,:,:),douobs(k,:,:),costt_rec(k)) 
    call subadjoint_rec(k,T_rec,residuals_time_rec(k,:),residuals_dou_rec(k,:,:),adj_rec(k,:,:))

    do i=1, nx
    do j=1, ny
       gradt_rec(k,i,j)=-adj_rec(k,i,j)/vel(i,j)**3
    enddo
    enddo
    
 enddo
 
  cost=sum(costt_rec(1:ns))

 do i=1, nx
 do j=1, ny
    gradtt(i,j)=sum(gradt_rec(1:ns,i,j))
 enddo
 enddo

   call gauss_smooth(gradtt,sgradtt)
 
 k=0
 do i=1, nx
 do j=1, ny
    k=k+1
    grad(k)=sgradtt(i,j)
 enddo
 enddo
 
 deallocate(adj_rec,padj_rec,T_rec)
 deallocate(Trcal,doucal,residuals_time_rec,residuals_dou_rec,res_deri_slope_rec)
 deallocate(gradt_rec,hesst_rec,costt_rec)
 deallocate(gradtt_rec,hesstt_rec)
 deallocate(gradtt,hesstt,sgradtt,pgradtt,spgradtt)

end subroutine subgradient

!**********************************************************************
! smooth the gradient
!**********************************************************************
subroutine gauss_smooth(g,g1)
 use globalp
 real*8:: g(1:nx,1:ny),g1(1:nx,1:ny),w(-nw1:nw1,-nw2:nw2)
 real*8:: sigma,sigma2,mywet,mysum,sigma0_1,sigma0_2,sigma_1,sigma_2,sigma1
 integer:: i,j,k,i1,j1,k1,i2,j2,k2,iter

 sigma0_1=dble(nw1)/3.
 sigma0_2=dble(nw2)/3.
 sigma_1=sigma0_1
 sigma_2=sigma0_2
 sigma1=sigma_1*sigma_1
 sigma2=sigma_2*sigma_2

 do j=-nw1, nw1
 do k=-nw2, nw2
    w(j,k)=exp(-0.5d0*(dble(j*j)/sigma1+dble(k*k)/sigma2))
 enddo
 enddo

 w=w/sum(w)
 g1=g

 do j=1, ny
 do i=1, nx
    mywet=0.0
    mysum=0.0
    do i1=-nw1, nw1
    do j1=-nw2, nw2
       k1=i+i1
       k2=j+j1
       if(k1.lt.1)  k1=1-k1
       if(k1.gt.nx) k1=2*nx-k1+1
       if(k2.lt.1)  k2=1-k2
       if(k2.gt.ny) k2=2*ny-k2+1
       mywet=mywet+w(i1,j1)
       mysum=mysum+w(i1,j1)*g(k1,k2)
    enddo
    enddo
    g1(i,j)=mysum/mywet
 enddo
 enddo

end subroutine

