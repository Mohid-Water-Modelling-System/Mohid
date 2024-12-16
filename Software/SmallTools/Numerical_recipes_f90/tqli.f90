	SUBROUTINE tqli(d,e,z)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : pythag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
	REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
	INTEGER(I4B) :: i,iter,l,m,n,ndum
	REAL(SP) :: b,c,dd,f,g,p,r,s
	REAL(SP), DIMENSION(size(e)) :: ff
	n=assert_eq(size(d),size(e),'tqli: n')
	if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
	e(:)=eoshift(e(:),1)
	do l=1,n
		iter=0
		iterate: do
			do m=l,n-1
				dd=abs(d(m))+abs(d(m+1))
				if (abs(e(m))+dd == dd) exit
			end do
			if (m == l) exit iterate
			if (iter == 30) call nrerror('too many iterations in tqli')
			iter=iter+1
			g=(d(l+1)-d(l))/(2.0_sp*e(l))
			r=pythag(g,1.0_sp)
			g=d(m)-d(l)+e(l)/(g+sign(r,g))
			s=1.0
			c=1.0
			p=0.0
			do i=m-1,l,-1
				f=s*e(i)
				b=c*e(i)
				r=pythag(f,g)
				e(i+1)=r
				if (r == 0.0) then
					d(i+1)=d(i+1)-p
					e(m)=0.0
					cycle iterate
				end if
				s=f/r
				c=g/r
				g=d(i+1)-p
				r=(d(i)-g)*s+2.0_sp*c*b
				p=s*r
				d(i+1)=g+p
				g=c*r-b
				if (present(z)) then
					ff(1:n)=z(1:n,i+1)
					z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
					z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
				end if
			end do
			d(l)=d(l)-p
			e(l)=g
			e(m)=0.0
		end do iterate
	end do
	END SUBROUTINE tqli
