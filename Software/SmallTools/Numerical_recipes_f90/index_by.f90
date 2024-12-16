	RECURSIVE SUBROUTINE index_bypack(arr,index,partial)
	USE nrtype; USE nrutil, ONLY : array_copy,arth,assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: index
	INTEGER, OPTIONAL, INTENT(IN) :: partial
	REAL(SP) :: a
	INTEGER(I4B) :: n,k,nl,indext,nerr
	INTEGER(I4B), SAVE :: level=0
	LOGICAL, DIMENSION(:), ALLOCATABLE, SAVE :: mask
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE, SAVE :: temp
	if (present(partial)) then
		n=size(index)
	else
		n=assert_eq(size(index),size(arr),'indexx_bypack')
		index=arth(1,1,n)
	end if
	if (n <= 1) RETURN
	k=(1+n)/2
	call icomp_xchg(index(1),index(k))
	call icomp_xchg(index(k),index(n))
	call icomp_xchg(index(1),index(k))
	if (n <= 3) RETURN
	level=level+1
	if (level == 1) allocate(mask(n),temp(n))
	indext=index(k)
	a=arr(indext)
	mask(1:n) = (arr(index) <= a)
	mask(k) = .false.
	call array_copy(pack(index,mask(1:n)),temp,nl,nerr)
	mask(k) = .true.
	temp(nl+2:n)=pack(index,.not. mask(1:n))
	temp(nl+1)=indext
	index=temp(1:n)
	call index_bypack(arr,index(1:nl),partial=1)
	call index_bypack(arr,index(nl+2:n),partial=1)
	if (level == 1) deallocate(mask,temp)
	level=level-1
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (arr(j) < arr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE index_bypack
