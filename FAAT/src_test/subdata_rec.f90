
subroutine subdata_rec(ishot,T,Tr,dou)
 use globalp
 integer:: i,j,k,ishot
 real*8:: T(nx,ny),Tr(nr),dou(nr,nr)
 
 Tr=0
 do k=1, recnum(ishot)
    Tr(k)=T(reci(ishot,k),recj(ishot,k))
 enddo

dou=0
do i=1, recnum(ishot)
do j=1, recnum(ishot)
   dou(i,j)=T(reci(ishot,j),recj(ishot,j))-T(reci(ishot,i),recj(ishot,i))
enddo
enddo 
   
end subroutine subdata_rec
