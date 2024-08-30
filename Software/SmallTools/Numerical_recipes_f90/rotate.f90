	SUBROUTINE rotate(r,qt,i,a,b)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
	INTEGER(I4B), INTENT(IN) :: i
	REAL(SP), INTENT(IN) :: a,b
	REAL(SP), DIMENSION(size(r,1)) :: temp
	INTEGER(I4B) :: n
	REAL(SP) :: c,fact,s
	n=assert_eq(size(r,1),size(r,2),size(qt,1),size(qt,2),'rotate')
	if (a == 0.0) then
		c=0.0
		s=sign(1.0_sp,b)
	else if (abs(a) > abs(b)) then
		fact=b/a
		c=sign(1.0_sp/sqrt(1.0_sp+fact**2),a)
		s=fact*c
	else
		fact=a/b
		s=sign(1.0_sp/sqrt(1.0_sp+fact**2),b)
		c=fact*s
	end if
	temp(i:n)=r(i,i:n)
	r(i,i:n)=c*temp(i:n)-s*r(i+1,i:n)
	r(i+1,i:n)=s*temp(i:n)+c*r(i+1,i:n)
	temp=qt(i,:)
	qt(i,:)=c*temp-s*qt(i+1,:)
	qt(i+1,:)=s*temp+c*qt(i+1,:)
	END SUBROUTINE rotate
