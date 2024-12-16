	SUBROUTINE select_heap(arr,heap)
	USE nrtype; USE nrutil, ONLY : nrerror,swap
	USE nr, ONLY : sort
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	REAL(SP), DIMENSION(:), INTENT(OUT) :: heap
	INTEGER(I4B) :: i,j,k,m,n
	m=size(heap)
	n=size(arr)
	if (m > n/2 .or. m < 1) call nrerror('probable misuse of select_heap')
	heap=arr(1:m)
	call sort(heap)
	do i=m+1,n
		if (arr(i) > heap(1)) then
			heap(1)=arr(i)
			j=1
			do
				k=2*j
				if (k > m) exit
				if (k /= m) then
					if (heap(k) > heap(k+1)) k=k+1
				end if
				if (heap(j) <= heap(k)) exit
				call swap(heap(k),heap(j))
				j=k
			end do
		end if
	end do
	END SUBROUTINE select_heap
