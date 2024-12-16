	SUBROUTINE balanc(a)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), PARAMETER :: RADX=radix(a),SQRADX=RADX**2
	INTEGER(I4B) :: i,last,ndum
	REAL(SP) :: c,f,g,r,s
	ndum=assert_eq(size(a,1),size(a,2),'balanc')
	do
		last=1
		do i=1,size(a,1)
			c=sum(abs(a(:,i)))-a(i,i)
			r=sum(abs(a(i,:)))-a(i,i)
			if (c /= 0.0 .and. r /= 0.0) then
				g=r/RADX
				f=1.0
				s=c+r
				do
					if (c >= g) exit
					f=f*RADX
					c=c*SQRADX
				end do
				g=r*RADX
				do
					if (c <= g) exit
					f=f/RADX
					c=c/SQRADX
				end do
				if ((c+r)/f < 0.95_sp*s) then
					last=0
					g=1.0_sp/f
					a(i,:)=a(i,:)*g
					a(:,i)=a(:,i)*f
				end if
			end if
		end do
		if (last /= 0) exit
	end do
	END SUBROUTINE balanc
