	SUBROUTINE svdcmp_sp(a,w,v)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerprod
	USE nr, ONLY : pythag
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: w
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
	INTEGER(I4B) :: i,its,j,k,l,m,n,nm
	REAL(SP) :: anorm,c,f,g,h,s,scale,x,y,z
	REAL(SP), DIMENSION(size(a,1)) :: tempm
	REAL(SP), DIMENSION(size(a,2)) :: rv1,tempn
	m=size(a,1)
	n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_sp')
	g=0.0
	scale=0.0
	do i=1,n
		l=i+1
		rv1(i)=scale*g
		g=0.0
		scale=0.0
		if (i <= m) then
			scale=sum(abs(a(i:m,i)))
			if (scale /= 0.0) then
				a(i:m,i)=a(i:m,i)/scale
				s=dot_product(a(i:m,i),a(i:m,i))
				f=a(i,i)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,i)=f-g
				tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
				a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
				a(i:m,i)=scale*a(i:m,i)
			end if
		end if
		w(i)=scale*g
		g=0.0
		scale=0.0
		if ((i <= m) .and. (i /= n)) then
			scale=sum(abs(a(i,l:n)))
			if (scale /= 0.0) then
				a(i,l:n)=a(i,l:n)/scale
				s=dot_product(a(i,l:n),a(i,l:n))
				f=a(i,l)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,l)=f-g
				rv1(l:n)=a(i,l:n)/h
				tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
				a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
				a(i,l:n)=scale*a(i,l:n)
			end if
		end if
	end do
	anorm=maxval(abs(w)+abs(rv1))
	do i=n,1,-1
		if (i < n) then
			if (g /= 0.0) then
				v(l:n,i)=(a(i,l:n)/a(i,l))/g
				tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
				v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
			end if
			v(i,l:n)=0.0
			v(l:n,i)=0.0
		end if
		v(i,i)=1.0
		g=rv1(i)
		l=i
	end do
	do i=min(m,n),1,-1
		l=i+1
		g=w(i)
		a(i,l:n)=0.0
		if (g /= 0.0) then
			g=1.0_sp/g
			tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
			a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
			a(i:m,i)=a(i:m,i)*g
		else
			a(i:m,i)=0.0
		end if
		a(i,i)=a(i,i)+1.0_sp
	end do
	do k=n,1,-1
		do its=1,30
			do l=k,1,-1
				nm=l-1
				if ((abs(rv1(l))+anorm) == anorm) exit
				if ((abs(w(nm))+anorm) == anorm) then
					c=0.0
					s=1.0
					do i=l,k
						f=s*rv1(i)
						rv1(i)=c*rv1(i)
						if ((abs(f)+anorm) == anorm) exit
						g=w(i)
						h=pythag(f,g)
						w(i)=h
						h=1.0_sp/h
						c= (g*h)
						s=-(f*h)
						tempm(1:m)=a(1:m,nm)
						a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
						a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
					end do
					exit
				end if
			end do
			z=w(k)
			if (l == k) then
				if (z < 0.0) then
					w(k)=-z
					v(1:n,k)=-v(1:n,k)
				end if
				exit
			end if
			if (its == 30) call nrerror('svdcmp_sp: no convergence in svdcmp')
			x=w(l)
			nm=k-1
			y=w(nm)
			g=rv1(nm)
			h=rv1(k)
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
			g=pythag(f,1.0_sp)
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
			c=1.0
			s=1.0
			do j=l,nm
				i=j+1
				g=rv1(i)
				y=w(i)
				h=s*g
				g=c*g
				z=pythag(f,h)
				rv1(j)=z
				c=f/z
				s=h/z
				f= (x*c)+(g*s)
				g=-(x*s)+(g*c)
				h=y*s
				y=y*c
				tempn(1:n)=v(1:n,j)
				v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
				v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
				z=pythag(f,h)
				w(j)=z
				if (z /= 0.0) then
					z=1.0_sp/z
					c=f*z
					s=h*z
				end if
				f= (c*g)+(s*y)
				x=-(s*g)+(c*y)
				tempm(1:m)=a(1:m,j)
				a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
				a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
			end do
			rv1(l)=0.0
			rv1(k)=f
			w(k)=x
		end do
	end do
	END SUBROUTINE svdcmp_sp

	SUBROUTINE svdcmp_dp(a,w,v)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerprod
	USE nr, ONLY : pythag
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(DP), DIMENSION(:), INTENT(OUT) :: w
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
	INTEGER(I4B) :: i,its,j,k,l,m,n,nm
	REAL(DP) :: anorm,c,f,g,h,s,scale,x,y,z
	REAL(DP), DIMENSION(size(a,1)) :: tempm
	REAL(DP), DIMENSION(size(a,2)) :: rv1,tempn
	m=size(a,1)
	n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_dp')
	g=0.0
	scale=0.0
	do i=1,n
		l=i+1
		rv1(i)=scale*g
		g=0.0
		scale=0.0
		if (i <= m) then
			scale=sum(abs(a(i:m,i)))
			if (scale /= 0.0) then
				a(i:m,i)=a(i:m,i)/scale
				s=dot_product(a(i:m,i),a(i:m,i))
				f=a(i,i)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,i)=f-g
				tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
				a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
				a(i:m,i)=scale*a(i:m,i)
			end if
		end if
		w(i)=scale*g
		g=0.0
		scale=0.0
		if ((i <= m) .and. (i /= n)) then
			scale=sum(abs(a(i,l:n)))
			if (scale /= 0.0) then
				a(i,l:n)=a(i,l:n)/scale
				s=dot_product(a(i,l:n),a(i,l:n))
				f=a(i,l)
				g=-sign(sqrt(s),f)
				h=f*g-s
				a(i,l)=f-g
				rv1(l:n)=a(i,l:n)/h
				tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
				a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
				a(i,l:n)=scale*a(i,l:n)
			end if
		end if
	end do
	anorm=maxval(abs(w)+abs(rv1))
	do i=n,1,-1
		if (i < n) then
			if (g /= 0.0) then
				v(l:n,i)=(a(i,l:n)/a(i,l))/g
				tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
				v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
			end if
			v(i,l:n)=0.0
			v(l:n,i)=0.0
		end if
		v(i,i)=1.0
		g=rv1(i)
		l=i
	end do
	do i=min(m,n),1,-1
		l=i+1
		g=w(i)
		a(i,l:n)=0.0
		if (g /= 0.0) then
			g=1.0_dp/g
			tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
			a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
			a(i:m,i)=a(i:m,i)*g
		else
			a(i:m,i)=0.0
		end if
		a(i,i)=a(i,i)+1.0_dp
	end do
	do k=n,1,-1
		do its=1,30
			do l=k,1,-1
				nm=l-1
				if ((abs(rv1(l))+anorm) == anorm) exit
				if ((abs(w(nm))+anorm) == anorm) then
					c=0.0
					s=1.0
					do i=l,k
						f=s*rv1(i)
						rv1(i)=c*rv1(i)
						if ((abs(f)+anorm) == anorm) exit
						g=w(i)
						h=pythag(f,g)
						w(i)=h
						h=1.0_dp/h
						c= (g*h)
						s=-(f*h)
						tempm(1:m)=a(1:m,nm)
						a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
						a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
					end do
					exit
				end if
			end do
			z=w(k)
			if (l == k) then
				if (z < 0.0) then
					w(k)=-z
					v(1:n,k)=-v(1:n,k)
				end if
				exit
			end if
			if (its == 30) call nrerror('svdcmp_dp: no convergence in svdcmp')
			x=w(l)
			nm=k-1
			y=w(nm)
			g=rv1(nm)
			h=rv1(k)
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
			g=pythag(f,1.0_dp)
			f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
			c=1.0
			s=1.0
			do j=l,nm
				i=j+1
				g=rv1(i)
				y=w(i)
				h=s*g
				g=c*g
				z=pythag(f,h)
				rv1(j)=z
				c=f/z
				s=h/z
				f= (x*c)+(g*s)
				g=-(x*s)+(g*c)
				h=y*s
				y=y*c
				tempn(1:n)=v(1:n,j)
				v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
				v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
				z=pythag(f,h)
				w(j)=z
				if (z /= 0.0) then
					z=1.0_dp/z
					c=f*z
					s=h*z
				end if
				f= (c*g)+(s*y)
				x=-(s*g)+(c*y)
				tempm(1:m)=a(1:m,j)
				a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
				a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
			end do
			rv1(l)=0.0
			rv1(k)=f
			w(k)=x
		end do
	end do
	END SUBROUTINE svdcmp_dp
