	FUNCTION select_bypack(k,arr)
	USE nrtype; USE nrutil, ONLY : array_copy,assert,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP) :: select_bypack
	LOGICAL, DIMENSION(size(arr)) :: mask
	REAL(SP), DIMENSION(size(arr)) :: temp
	INTEGER(I4B) :: i,r,j,l,n,nl,nerr
	REAL(SP) :: a
	n=size(arr)
	call assert(k >= 1, k <= n, 'select_bypack args')
	l=1
	r=n
	do
		if (r-l <= 1) exit
		i=(l+r)/2
		call swap(arr(l),arr(i),arr(l)>arr(i))
		call swap(arr(i),arr(r),arr(i)>arr(r))
		call swap(arr(l),arr(i),arr(l)>arr(i))
		a=arr(i)
		mask(l:r) = (arr(l:r) <= a)
		mask(i) = .false.
		call array_copy(pack(arr(l:r),mask(l:r)),temp(l:),nl,nerr)
		j=l+nl
		mask(i) = .true.
		temp(j+1:r)=pack(arr(l:r),.not. mask(l:r))
		temp(j)=a
		arr(l:r)=temp(l:r)
		if (k > j) then
			l=j+1
		else if (k < j) then
			r=j-1
		else
			l=j
			r=j
		end if
	end do
	if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r))
	select_bypack=arr(k)
	END FUNCTION select_bypack
