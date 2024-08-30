	SUBROUTINE sort_pick(arr)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,j,n
	REAL(SP) :: a
	n=size(arr)
	do j=2,n
		a=arr(j)
		do i=j-1,1,-1
			if (arr(i) <= a) exit
			arr(i+1)=arr(i)
		end do
		arr(i+1)=a
	end do
	END SUBROUTINE sort_pick
