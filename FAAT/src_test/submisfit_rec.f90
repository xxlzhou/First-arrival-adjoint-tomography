subroutine submisfit_rec(ishot,trcal,trobs,doucal,douobs,cost)
 use globalp
 integer:: i,j,k,ishot
 real*8:: cost,cost1,cost2
 real*8::trcal(nr),trobs(nr),doucal(nr,nr),douobs(nr,nr)
 
 cost1=0.
 do i=1,recnum(ishot)
    cost1=cost1+0.5*W_t*(trcal(i)-trobs(i))**2
 enddo

 cost2=0
 do j=1, recnum(ishot)-1
 do k=j+1, recnum(ishot)
     cost2=cost2+0.5*W_pr*(doucal(j,k)-douobs(j,k))**2
 enddo
 enddo
 
 cost=cost1+cost2
 
end subroutine submisfit_rec
