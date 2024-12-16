	RECURSIVE SUBROUTINE sort_bypack(arr)
	USE nrtype; USE nrutil, ONLY : array_copy,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,nl,nerr
	INTEGER(I4B), SAVE :: level=0
	LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: mask
	REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: temp
	n=size(arr)
	if (n <= 1) RETURN
	k=(1+n)/2
	call swap(arr(1),arr(k),arr(1)>arr(k))
	call swap(arr(k),arr(n),arr(k)>arr(n))
	call swap(arr(1),arr(k),arr(1)>arr(k))
	if (n <= 3) RETURN
	level=level+1
	if (level == 1) allocate(mask(n),temp(n))
	a=arr(k)
	mask(1:n) = (arr <= a)
	mask(k) = .false.
	call array_copy(pack(arr,mask(1:n)),temp,nl,nerr)
	mask(k) = .true.
	temp(nl+2:n)=pack(arr,.not. mask(1:n))
	temp(nl+1)=a
	arr=temp(1:n)
	call sort_bypack(arr(1:nl))
	call sort_bypack(arr(nl+2:n))
	if (level == 1) deallocate(mask,temp)
	level=level-1
	END SUBROUTINE sort_bypack
