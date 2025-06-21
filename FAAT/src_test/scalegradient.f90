  subroutine scalegradient(grad,cost,m,perc,alpha)
 use globalp
 integer:: i,imax
 real:: perc,cost,gradmax,alpha,pgradnorm,gradnorm
 real:: m(nop),grad(nop)

 gradmax=-1.
 imax=-1
 
 do i=1,nop
    if (abs(grad(i)).gt.gradmax) then
       gradmax=abs(grad(i))
       imax=i
    endif
 enddo
 
 alpha=m(imax)*perc/gradmax
 grad(:)=grad(:)*alpha
 cost=cost*alpha
 
 end subroutine scalegradient

