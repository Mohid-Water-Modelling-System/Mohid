	FUNCTION wwghts(n,h,kermom)
	USE nrtype; USE nrutil, ONLY : geop
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: h
	REAL(SP), DIMENSION(n) :: wwghts
	INTERFACE
		FUNCTION kermom(y,m)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: y
		INTEGER(I4B), INTENT(IN) :: m
		REAL(DP), DIMENSION(m) :: kermom
		END FUNCTION kermom
	END INTERFACE
	INTEGER(I4B) :: j
	REAL(DP) :: hh,hi,c,a,b
	REAL(DP), DIMENSION(4) :: wold,wnew,w
	hh=h
	hi=1.0_dp/hh
	wwghts(1:n)=0.0
	wold(1:4)=kermom(0.0_dp,4)
	if (n >= 4) then
		b=0.0
		do j=1,n-3
			c=j-1
			a=b
			b=a+hh
			if (j == n-3) b=(n-1)*hh
			wnew(1:4)=kermom(b,4)
			w(1:4)=(wnew(1:4)-wold(1:4))*geop(1.0_dp,hi,4)
			wwghts(j:j+3)=wwghts(j:j+3)+(/&
				((c+1.0_dp)*(c+2.0_dp)*(c+3.0_dp)*w(1)&
				-(11.0_dp+c*(12.0_dp+c*3.0_dp))*w(2)&
					+3.0_dp*(c+2.0_dp)*w(3)-w(4))/6.0_dp,&
				(-c*(c+2.0_dp)*(c+3.0_dp)*w(1)&
				+(6.0_dp+c*(10.0_dp+c*3.0_dp))*w(2)&
					-(3.0_dp*c+5.0_dp)*w(3)+w(4))*0.50_dp,&
				(c*(c+1.0_dp)*(c+3.0_dp)*w(1)&
				-(3.0_dp+c*(8.0_dp+c*3.0_dp))*w(2)&
					+(3.0_dp*c+4.0_dp)*w(3)-w(4))*0.50_dp,&
				(-c*(c+1.0_dp)*(c+2.0_dp)*w(1)&
				+(2.0_dp+c*(6.0_dp+c*3.0_dp))*w(2)&
				-3.0_dp*(c+1.0_dp)*w(3)+w(4))/6.0_dp /)
			wold(1:4)=wnew(1:4)
		end do
	else if (n == 3) then
		wnew(1:3)=kermom(hh+hh,3)
		w(1:3)= (/ wnew(1)-wold(1), hi*(wnew(2)-wold(2)),&
			hi**2*(wnew(3)-wold(3)) /)
		wwghts(1:3)= (/ w(1)-1.50_dp*w(2)+0.50_dp*w(3),&
			2.0_dp*w(2)-w(3), 0.50_dp*(w(3)-w(2)) /)
	else if (n == 2) then
		wnew(1:2)=kermom(hh,2)
		wwghts(2)=hi*(wnew(2)-wold(2))
		wwghts(1)=wnew(1)-wold(1)-wwghts(2)
	end if
	END FUNCTION wwghts
