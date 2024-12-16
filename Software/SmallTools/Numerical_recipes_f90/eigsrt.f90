	SUBROUTINE eigsrt(d,v)
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: d
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: v
	INTEGER(I4B) :: i,j,n
	n=assert_eq(size(d),size(v,1),size(v,2),'eigsrt')
	do i=1,n-1
		j=imaxloc(d(i:n))+i-1
		if (j /= i) then
			call swap(d(i),d(j))
			call swap(v(:,i),v(:,j))
		end if
	end do
	END SUBROUTINE eigsrt
