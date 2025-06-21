 subroutine subm2v(m,v)
 use globalp
 integer:: i,j,k
 real:: m(nop)
 real*8::v(nx,ny)

 k=0
 do i=1, nx
 do j=1, ny
    k=k+1
    v(i,j)=m(k)
 enddo
 enddo
 
 end subroutine subm2v
 
 
 subroutine subv2m(v,m)
 use globalp
 integer:: i,j,k
 real:: m(nop)
 real*8::v(nx,ny)
 
 k=0
 do i=1, nx
 do j=1, ny
    k=k+1
    m(k)=v(i,j)
 enddo
 enddo
 
 end subroutine subv2m
