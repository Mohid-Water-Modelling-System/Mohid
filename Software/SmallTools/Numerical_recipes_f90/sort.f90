	SUBROUTINE sort(arr)
	USE nrtype; USE nrutil, ONLY : swap,nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,i,j,jstack,l,r
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	n=size(arr)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				a=arr(j)
				do i=j-1,l,-1
					if (arr(i) <= a) exit
					arr(i+1)=arr(i)
				end do
				arr(i+1)=a
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(arr(k),arr(l+1))
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
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('sort: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	END SUBROUTINE sort
