	PROGRAM badluk
	USE nrtype
	USE nr, ONLY : flmoon,julday
	IMPLICIT NONE
	INTEGER(I4B) :: ic,icon,idwk,ifrac,im,iyyy,jd,jday,n
	INTEGER(I4B) :: iybeg=1900,iyend=2000
	REAL(SP) :: frac
	REAL(SP), PARAMETER :: TIMZON=-5.0_sp/24.0_sp
	write (*,'(1x,a,i5,a,i5)') 'Full moons on Friday the 13th from',&
		iybeg,' to',iyend
	do iyyy=iybeg,iyend
		do im=1,12
			jday=julday(im,13,iyyy)
			idwk=mod(jday+1,7)
			if (idwk == 5) then
				n=12.37_sp*(iyyy-1900+(im-0.5_sp)/12.0_sp)
				icon=0
				do
					call flmoon(n,2,jd,frac)
					ifrac=nint(24.0_sp*(frac+TIMZON))
					if (ifrac < 0) then
						jd=jd-1
						ifrac=ifrac+24
					end if
					if (ifrac > 12) then
						jd=jd+1
						ifrac=ifrac-12
					else
						ifrac=ifrac+12
					end if
					if (jd == jday) then
						write (*,'(/1x,i2,a,i2,a,i4)') im,'/',13,'/',iyyy
						write (*,'(1x,a,i2,a)') 'Full moon ',ifrac,&
							' hrs after midnight (EST).'
						exit
					else
						ic=isign(1,jday-jd)
						if (ic == -icon) exit
						icon=ic
						n=n+ic
					end if
				end do
			end if
		end do
	end do
	END PROGRAM badluk
