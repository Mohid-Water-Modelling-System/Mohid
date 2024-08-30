	FUNCTION recur2(a,b,c)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a,b,c
	REAL(SP), DIMENSION(size(a)) :: recur2
	INTEGER(I4B) :: n
	REAL(SP), DIMENSION(size(a)-1) :: a1,a2,u1,u2
	REAL(SP), DIMENSION(size(a)-2) :: b11,b12,b21,b22
	n=assert_eq(size(a),size(b)+2,size(c)+2,'recur2')
	a1(1)=a(1)
	a2(1)=a(2)
	a1(2:n-1)=0.0
	a2(2:n-1)=a(3:n)
	b11(1:n-2)=0.0
	b12(1:n-2)=1.0
	b21(1:n-2)=c(1:n-2)
	b22(1:n-2)=b(1:n-2)
	call recur1_v(a1,a2,b11,b12,b21,b22,u1,u2)
	recur2(1:n-1)=u1(1:n-1)
	recur2(n)=u2(n-1)
	CONTAINS
!BL
	RECURSIVE SUBROUTINE recur1_v(a1,a2,b11,b12,b21,b22,u1,u2)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: a1,a2,b11,b12,b21,b22
	REAL(SP), DIMENSION(:), INTENT(OUT) :: u1,u2
	INTEGER(I4B), PARAMETER :: NPAR_RECUR2=8
	INTEGER(I4B) :: n,j,nn,nn1
	REAL(SP), DIMENSION(size(a1)/2) :: aa1,aa2
	REAL(SP), DIMENSION(size(a1)/2-1) :: bb11,bb12,bb21,bb22
	n=assert_eq((/size(a1),size(a2),size(b11)+1,size(b12)+1,size(b21)+1,&
		size(b22)+1,size(u1),size(u2)/),'recur1_v')
	u1(1)=a1(1)
	u2(1)=a2(1)
	if (n < NPAR_RECUR2) then
		do j=2,n
			u1(j)=a1(j)+b11(j-1)*u1(j-1)+b12(j-1)*u2(j-1)
			u2(j)=a2(j)+b21(j-1)*u1(j-1)+b22(j-1)*u2(j-1)
		end do
	else
		nn=n/2
		nn1=nn-1
		aa1(1:nn)=a1(2:n:2)+b11(1:n-1:2)*a1(1:n-1:2)+&
			b12(1:n-1:2)*a2(1:n-1:2)
		aa2(1:nn)=a2(2:n:2)+b21(1:n-1:2)*a1(1:n-1:2)+&
				b22(1:n-1:2)*a2(1:n-1:2)
		bb11(1:nn1)=b11(3:n-1:2)*b11(2:n-2:2)+&
				b12(3:n-1:2)*b21(2:n-2:2)
		bb12(1:nn1)=b11(3:n-1:2)*b12(2:n-2:2)+&
				b12(3:n-1:2)*b22(2:n-2:2)
		bb21(1:nn1)=b21(3:n-1:2)*b11(2:n-2:2)+&
				b22(3:n-1:2)*b21(2:n-2:2)
		bb22(1:nn1)=b21(3:n-1:2)*b12(2:n-2:2)+&
				b22(3:n-1:2)*b22(2:n-2:2)
		call recur1_v(aa1,aa2,bb11,bb12,bb21,bb22,u1(2:n:2),u2(2:n:2))
		u1(3:n:2)=a1(3:n:2)+b11(2:n-1:2)*u1(2:n-1:2)+&
				b12(2:n-1:2)*u2(2:n-1:2)
		u2(3:n:2)=a2(3:n:2)+b21(2:n-1:2)*u1(2:n-1:2)+&
				b22(2:n-1:2)*u2(2:n-1:2)
	end if
	END SUBROUTINE recur1_v
	END FUNCTION recur2
