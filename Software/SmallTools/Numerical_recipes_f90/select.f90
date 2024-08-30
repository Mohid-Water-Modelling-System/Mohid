	FUNCTION select(k,arr)
	USE nrtype; USE nrutil, ONLY : assert,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: k
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	REAL(SP) :: select
	INTEGER(I4B) :: i,r,j,l,n
	REAL(SP) :: a
	n=size(arr)
	call assert(k >= 1, k <= n, 'select args')
	l=1
	r=n
	do
		if (r-l <= 1) then
			if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r))
			select=arr(k)
			RETURN
		else
			i=(l+r)/2
			call swap(arr(i),arr(l+1))
			call swap(arr(l),arr(r),arr(l)>arr(r))
			call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
			call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
			i=l+1
			j=r
			a=arr(l+1)
			do
				do
					i=i+1
					if (arr(i) >= a) exit
				end do
				do
					j=j-1
					if (arr(j) <= a) exit
				end do
				if (j < i) exit
				call swap(arr(i),arr(j))
			end do
			arr(l+1)=arr(j)
			arr(j)=a
			if (j >= k) r=j-1
			if (j <= k) l=i
		end if
	end do
	END FUNCTION select
