	FUNCTION decchk(string,ch)
	USE nrtype; USE nrutil, ONLY : ifirstloc
	IMPLICIT NONE
	CHARACTER(1), DIMENSION(:), INTENT(IN) :: string
	CHARACTER(1), INTENT(OUT) :: ch
	LOGICAL(LGT) :: decchk
	INTEGER(I4B) :: i,j,k,m
	INTEGER(I4B) :: ip(0:9,0:7) = reshape((/ &
		0,1,2,3,4,5,6,7,8,9,1,5,7,6,2,8,3,0,9,4,&
		5,8,0,3,7,9,6,1,4,2,8,9,1,6,0,4,3,5,2,7,9,4,5,3,1,2,6,8,7,0,&
		4,2,8,6,5,7,3,9,0,1,2,7,9,3,8,0,6,4,1,5,7,0,4,6,9,1,3,2,5,8 /),&
		(/ 10,8 /) )
	INTEGER(I4B) :: ij(0:9,0:9) = reshape((/ &
		0,1,2,3,4,5,6,7,8,9,1,2,3,4,0,9,5,6,7,8,2,3,4,0,1,8,9,5,6,&
		7,3,4,0,1,2,7,8,9,5,6,4,0,1,2,3,6,7,8,9,5,5,6,7,8,9,0,1,2,3,&
		4,6,7,8,9,5,4,0,1,2,3,7,8,9,5,6,3,4,0,1,2,8,9,5,6,7,2,3,4,0,&
		1,9,5,6,7,8,1,2,3,4,0 /),(/ 10,10 /))
	k=0
	m=0
	do j=1,size(string)
		i=ichar(string(j))
		if (i >= 48 .and. i <= 57) then
			k=ij(k,ip(mod(i+2,10),mod(m,8)))
			m=m+1
		end if
	end do
	decchk=logical(k == 0,kind=lgt)
	i=mod(m,8)
	i=ifirstloc(ij(k,ip(0:9,i)) == 0)-1
	ch=char(i+48)
	END FUNCTION decchk
