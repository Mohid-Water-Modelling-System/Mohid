	SUBROUTINE hqr(a,wr,wi)
	USE nrtype; USE nrutil, ONLY : assert_eq,diagadd,nrerror,upper_triangle
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: wr,wi
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B) :: i,its,k,l,m,n,nn,mnnk
	REAL(SP) :: anorm,p,q,r,s,t,u,v,w,x,y,z
	REAL(SP), DIMENSION(size(a,1)) :: pp
	n=assert_eq(size(a,1),size(a,2),size(wr),size(wi),'hqr')
	anorm=sum(abs(a),mask=upper_triangle(n,n,extra=2))
	nn=n
	t=0.0
	do
		if (nn < 1) exit
		its=0
		iterate: do
			do l=nn,2,-1
				s=abs(a(l-1,l-1))+abs(a(l,l))
				if (s == 0.0) s=anorm
				if (abs(a(l,l-1))+s == s) exit
			end do
			x=a(nn,nn)
			if (l == nn) then
				wr(nn)=x+t
				wi(nn)=0.0
				nn=nn-1
				exit iterate
			end if
			y=a(nn-1,nn-1)
			w=a(nn,nn-1)*a(nn-1,nn)
			if (l == nn-1) then
				p=0.5_sp*(y-x)
				q=p**2+w
				z=sqrt(abs(q))
				x=x+t
				if (q >= 0.0) then
					z=p+sign(z,p)
					wr(nn)=x+z
					wr(nn-1)=wr(nn)
					if (z /= 0.0) wr(nn)=x-w/z
					wi(nn)=0.0
					wi(nn-1)=0.0
				else
					wr(nn)=x+p
					wr(nn-1)=wr(nn)
					wi(nn)=z
					wi(nn-1)=-z
				end if
				nn=nn-2
				exit iterate
			end if
			if (its == 30) call nrerror('too many iterations in hqr')
			if (its == 10 .or. its == 20) then
				t=t+x
				call diagadd(a(1:nn,1:nn),-x)
				s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
				x=0.75_sp*s
				y=x
				w=-0.4375_sp*s**2
			end if
			its=its+1
			do m=nn-2,l,-1
				z=a(m,m)
				r=x-z
				s=y-z
				p=(r*s-w)/a(m+1,m)+a(m,m+1)
				q=a(m+1,m+1)-z-r-s
				r=a(m+2,m+1)
				s=abs(p)+abs(q)+abs(r)
				p=p/s
				q=q/s
				r=r/s
				if (m == l) exit
				u=abs(a(m,m-1))*(abs(q)+abs(r))
				v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
				if (u+v == v) exit
			end do
			do i=m+2,nn
				a(i,i-2)=0.0
				if (i /= m+2) a(i,i-3)=0.0
			end do
			do k=m,nn-1
				if (k /= m) then
					p=a(k,k-1)
					q=a(k+1,k-1)
					r=0.0
					if (k /= nn-1) r=a(k+2,k-1)
					x=abs(p)+abs(q)+abs(r)
					if (x /= 0.0) then
						p=p/x
						q=q/x
						r=r/x
					end if
				end if
				s=sign(sqrt(p**2+q**2+r**2),p)
				if (s /= 0.0) then
					if (k == m) then
						if (l /= m) a(k,k-1)=-a(k,k-1)
					else
						a(k,k-1)=-s*x
					end if
					p=p+s
					x=p/s
					y=q/s
					z=r/s
					q=q/p
					r=r/p
					pp(k:nn)=a(k,k:nn)+q*a(k+1,k:nn)
					if (k /= nn-1) then
						pp(k:nn)=pp(k:nn)+r*a(k+2,k:nn)
						a(k+2,k:nn)=a(k+2,k:nn)-pp(k:nn)*z
					end if
					a(k+1,k:nn)=a(k+1,k:nn)-pp(k:nn)*y
					a(k,k:nn)=a(k,k:nn)-pp(k:nn)*x
					mnnk=min(nn,k+3)
					pp(l:mnnk)=x*a(l:mnnk,k)+y*a(l:mnnk,k+1)
					if (k /= nn-1) then
						pp(l:mnnk)=pp(l:mnnk)+z*a(l:mnnk,k+2)
						a(l:mnnk,k+2)=a(l:mnnk,k+2)-pp(l:mnnk)*r
					end if
					a(l:mnnk,k+1)=a(l:mnnk,k+1)-pp(l:mnnk)*q
					a(l:mnnk,k)=a(l:mnnk,k)-pp(l:mnnk)
				end if
			end do
		end do iterate
	end do
	END SUBROUTINE hqr
