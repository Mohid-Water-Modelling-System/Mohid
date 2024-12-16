	SUBROUTINE kstwo(data1,data2,d,prob)
	USE nrtype; USE nrutil, ONLY : cumsum
	USE nr, ONLY : probks,sort2
	IMPLICIT NONE
	REAL(SP), INTENT(OUT) :: d,prob
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	INTEGER(I4B) :: n1,n2
	REAL(SP) :: en1,en2,en
	REAL(SP), DIMENSION(size(data1)+size(data2)) :: dat,org
	n1=size(data1)
	n2=size(data2)
	en1=n1
	en2=n2
	dat(1:n1)=data1
	dat(n1+1:)=data2
	org(1:n1)=0.0
	org(n1+1:)=1.0
	call sort2(dat,org)
	d=maxval(abs(cumsum(org)/en2-cumsum(1.0_sp-org)/en1))
	en=sqrt(en1*en2/(en1+en2))
	prob=probks((en+0.12_sp+0.11_sp/en)*d)
	END SUBROUTINE kstwo
