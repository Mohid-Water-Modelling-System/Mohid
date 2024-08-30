	FUNCTION eclass(lista,listb,n)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: lista,listb
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), DIMENSION(n) :: eclass
	INTEGER :: j,k,l,m
	m=assert_eq(size(lista),size(listb),'eclass')
	eclass(1:n)=arth(1,1,n)
	do l=1,m
		j=lista(l)
		do
			if (eclass(j) == j) exit
			j=eclass(j)
		end do
		k=listb(l)
		do
			if (eclass(k) == k) exit
			k=eclass(k)
		end do
		if (j /= k) eclass(j)=k
	end do
	do j=1,n
		do
			if (eclass(j) == eclass(eclass(j))) exit
			eclass(j)=eclass(eclass(j))
		end do
	end do
	END FUNCTION eclass
