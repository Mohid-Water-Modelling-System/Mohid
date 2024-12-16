	FUNCTION eclazz(equiv,n)
	USE nrtype; USE nrutil, ONLY : arth
	IMPLICIT NONE
	INTERFACE
		FUNCTION equiv(i,j)
		USE nrtype
		IMPLICIT NONE
		LOGICAL(LGT) :: equiv
		INTEGER(I4B), INTENT(IN) :: i,j
		END FUNCTION equiv
	END INTERFACE
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B), DIMENSION(n) :: eclazz
	INTEGER :: i,j
	eclazz(1:n)=arth(1,1,n)
	do i=2,n
		do j=1,i-1
			eclazz(j)=eclazz(eclazz(j))
			if (equiv(i,j)) eclazz(eclazz(eclazz(j)))=i
		end do
	end do
	do i=1,n
		eclazz(i)=eclazz(eclazz(i))
	end do
	END FUNCTION eclazz
