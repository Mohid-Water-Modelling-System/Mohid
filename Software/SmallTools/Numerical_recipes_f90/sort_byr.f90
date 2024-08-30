	SUBROUTINE sort_byreshape(arr)
	USE nrtype; USE nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP), DIMENSION(:,:), ALLOCATABLE :: tab
	REAL(SP), PARAMETER :: big=huge(arr)
	INTEGER(I4B) :: inc,n,m
	n=size(arr)
	inc=1
	do
		inc=2*inc+1
		if (inc > n) exit
	end do
	do
		inc=inc/2
		m=(n+inc-1)/inc
		allocate(tab(inc,m))
		tab=reshape(arr, (/inc,m/) , (/big/) )
		do
			call swap(tab(:,1:m-1:2),tab(:,2:m:2), &
				tab(:,1:m-1:2)>tab(:,2:m:2))
			call swap(tab(:,2:m-1:2),tab(:,3:m:2), &
				tab(:,2:m-1:2)>tab(:,3:m:2))
			if (all(tab(:,1:m-1) <= tab(:,2:m))) exit
		end do
		arr=reshape(tab,shape(arr))
		deallocate(tab)
		if (inc <= 1) exit
	end do
	END SUBROUTINE sort_byreshape
