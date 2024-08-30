	SUBROUTINE sort_shell(arr)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,j,inc,n
	REAL(SP) :: v
	n=size(arr)
	inc=1
	do
		inc=3*inc+1
		if (inc > n) exit
	end do
	do
		inc=inc/3
		do i=inc+1,n
			v=arr(i)
			j=i
			do
				if (arr(j-inc) <= v) exit
				arr(j)=arr(j-inc)
				j=j-inc
				if (j <= inc) exit
			end do
			arr(j)=v
		end do
		if (inc <= 1) exit
	end do
	END SUBROUTINE sort_shell
