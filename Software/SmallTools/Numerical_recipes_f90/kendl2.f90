	SUBROUTINE kendl2(tab,tau,z,prob)
	USE nrtype; USE nrutil, ONLY : cumsum
	USE nr, ONLY : erfcc
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: tab
	REAL(SP), INTENT(OUT) :: tau,z,prob
	REAL(SP), DIMENSION(size(tab,1),size(tab,2)) :: cum,cumt
	INTEGER(I4B) :: i,j,ii,jj
	REAL(SP) :: sc,sd,en1,en2,points,var
	ii=size(tab,1)
	jj=size(tab,2)
	do i=1,ii
		cumt(i,jj:1:-1)=cumsum(tab(i,jj:1:-1))
	end do
	en2=sum(tab(1:ii,1:jj-1)*cumt(1:ii,2:jj))
	do j=1,jj
		cum(ii:1:-1,j)=cumsum(cumt(ii:1:-1,j))
	end do
	points=cum(1,1)
	sc=sum(tab(1:ii-1,1:jj-1)*cum(2:ii,2:jj))
	do j=1,jj
		cum(1:ii,j)=cumsum(cumt(1:ii,j))
	end do
	sd=sum(tab(2:ii,1:jj-1)*cum(1:ii-1,2:jj))
	do j=1,jj
		cumt(ii:1:-1,j)=cumsum(tab(ii:1:-1,j))
	end do
	en1=sum(tab(1:ii-1,1:jj)*cumt(2:ii,1:jj))
	tau=(sc-sd)/sqrt((en1+sc+sd)*(en2+sc+sd))
	var=(4.0_sp*points+10.0_sp)/(9.0_sp*points*(points-1.0_sp))
	z=tau/sqrt(var)
	prob=erfcc(abs(z)/SQRT2)
	END SUBROUTINE kendl2
