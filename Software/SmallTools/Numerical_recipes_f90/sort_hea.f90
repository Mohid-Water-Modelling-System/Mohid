	SUBROUTINE sort_heap(arr)
	USE nrtype
	USE nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,n
	n=size(arr)
	do i=n/2,1,-1
		call sift_down(i,n)
	end do
	do i=n,2,-1
		call swap(arr(1),arr(i))
		call sift_down(1,i-1)
	end do
	CONTAINS
!BL
	SUBROUTINE sift_down(l,r)
	INTEGER(I4B), INTENT(IN) :: l,r
	INTEGER(I4B) :: j,jold
	REAL(SP) :: a
	a=arr(l)
	jold=l
	j=l+l
	do
		if (j > r) exit
		if (j < r) then
			if (arr(j) < arr(j+1)) j=j+1
		end if
		if (a >= arr(j)) exit
		arr(jold)=arr(j)
		jold=j
		j=j+j
	end do
	arr(jold)=a
	END SUBROUTINE sift_down
	END SUBROUTINE sort_heap
