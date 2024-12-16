	SUBROUTINE daub4(a,isign)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: isign
	REAL(SP), DIMENSION(size(a)) :: wksp
	REAL(SP), PARAMETER :: C0=0.4829629131445341_sp,&
		C1=0.8365163037378079_sp,C2=0.2241438680420134_sp,&
		C3=-0.1294095225512604_sp
	INTEGER(I4B) :: n,nh,nhp,nhm
	n=size(a)
	if (n < 4) RETURN
	nh=n/2
	nhp=nh+1
	nhm=nh-1
	if (isign >= 0) then
		wksp(1:nhm) = C0*a(1:n-3:2)+C1*a(2:n-2:2) &
			+C2*a(3:n-1:2)+C3*a(4:n:2)
		wksp(nh)=C0*a(n-1)+C1*a(n)+C2*a(1)+C3*a(2)
		wksp(nhp:n-1) = C3*a(1:n-3:2)-C2*a(2:n-2:2) &
			+C1*a(3:n-1:2)-C0*a(4:n:2)
		wksp(n)=C3*a(n-1)-C2*a(n)+C1*a(1)-C0*a(2)
	else
		wksp(1)=C2*a(nh)+C1*a(n)+C0*a(1)+C3*a(nhp)
		wksp(2)=C3*a(nh)-C0*a(n)+C1*a(1)-C2*a(nhp)
		wksp(3:n-1:2) = C2*a(1:nhm)+C1*a(nhp:n-1) &
			+C0*a(2:nh)+C3*a(nh+2:n)
		wksp(4:n:2) = C3*a(1:nhm)-C0*a(nhp:n-1) &
			+C1*a(2:nh)-C2*a(nh+2:n)
	end if
	a(1:n)=wksp(1:n)
	END SUBROUTINE daub4
