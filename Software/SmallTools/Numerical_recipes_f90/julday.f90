	FUNCTION julday(mm,id,iyyy)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: mm,id,iyyy
	INTEGER(I4B) :: julday
	INTEGER(I4B), PARAMETER :: IGREG=15+31*(10+12*1582)
	INTEGER(I4B) :: ja,jm,jy
	jy=iyyy
	if (jy == 0) call nrerror('julday: there is no year zero')
	if (jy < 0) jy=jy+1
	if (mm > 2) then
		jm=mm+1
	else
		jy=jy-1
		jm=mm+13
	end if
	julday=int(365.25_sp*jy)+int(30.6001_sp*jm)+id+1720995
	if (id+31*(mm+12*iyyy) >= IGREG) then
		ja=int(0.01_sp*jy)
		julday=julday+2-ja+int(0.25_sp*ja)
	end if
	END FUNCTION julday
