subroutine subresiduals_rec(ishot,trcal,trobs,doucal,douobs,residuals_time_rec,&
            &residuals_dou_rec)
 use globalp
 integer::i,j,k,ishot
 real*8::trcal(nr),trobs(nr),doucal(nr,nr),douobs(nr,nr)
 real*8::residuals_time_rec(nr),residuals_dou_rec(nr,nr)

 residuals_time_rec=0
 do i=1, recnum(ishot)
    residuals_time_rec(i)=W_t*(trcal(i)-trobs(i))
 enddo

 residuals_dou_rec=0
 do j=1, recnum(ishot)
 do k=1, recnum(ishot)
    residuals_dou_rec(j,k)=W_pr*(doucal(j,k)-douobs(j,k))
 enddo
 enddo

 end subroutine subresiduals_rec
