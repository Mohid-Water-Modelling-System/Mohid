	SUBROUTINE tred2(a,d,e,novectors)
	USE nrtype; USE nrutil, ONLY : assert_eq,outerprod
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
	LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
	INTEGER(I4B) :: i,j,l,n
	REAL(SP) :: f,g,h,hh,scale
	REAL(SP), DIMENSION(size(a,1)) :: gg
	LOGICAL(LGT), SAVE :: yesvec=.true.
	n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
	if (present(novectors)) yesvec=.not. novectors
	do i=n,2,-1
		l=i-1
		h=0.0
		if (l > 1) then
			scale=sum(abs(a(i,1:l)))
			if (scale == 0.0) then
				e(i)=a(i,l)
			else
				a(i,1:l)=a(i,1:l)/scale
				h=sum(a(i,1:l)**2)
				f=a(i,l)
				g=-sign(sqrt(h),f)
				e(i)=scale*g
				h=h-f*g
				a(i,l)=f-g
				if (yesvec) a(1:l,i)=a(i,1:l)/h
				do j=1,l
					e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
					+dot_product(a(j+1:l,j),a(i,j+1:l)))/h
				end do
				f=dot_product(e(1:l),a(i,1:l))
				hh=f/(h+h)
				e(1:l)=e(1:l)-hh*a(i,1:l)
				do j=1,l
					a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
				end do
			end if
		else
			e(i)=a(i,l)
		end if
		d(i)=h
	end do
	if (yesvec) d(1)=0.0
	e(1)=0.0
	do i=1,n
		if (yesvec) then
			l=i-1
			if (d(i) /= 0.0) then
				gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
				a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
			end if
			d(i)=a(i,i)
			a(i,i)=1.0
			a(i,1:l)=0.0
			a(1:l,i)=0.0
		else
			d(i)=a(i,i)
		end if
	end do
	END SUBROUTINE tred2
