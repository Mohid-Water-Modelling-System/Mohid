	FUNCTION chebpc(c)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: c
	REAL(SP), DIMENSION(size(c)) :: chebpc
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(c)) :: dd,sv
	n=size(c)
	chebpc=0.0
	dd=0.0
	chebpc(1)=c(n)
	do j=n-1,2,-1
		sv(2:n-j+1)=chebpc(2:n-j+1)
		chebpc(2:n-j+1)=2.0_sp*chebpc(1:n-j)-dd(2:n-j+1)
		dd(2:n-j+1)=sv(2:n-j+1)
		sv(1)=chebpc(1)
		chebpc(1)=-dd(1)+c(j)
		dd(1)=sv(1)
	end do
	chebpc(2:n)=chebpc(1:n-1)-dd(2:n)
	chebpc(1)=-dd(1)+0.5_sp*c(1)
	END FUNCTION chebpc
