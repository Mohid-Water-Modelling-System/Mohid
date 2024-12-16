	SUBROUTINE arcode(ich,codep,lcd,isign,acode)
	USE nrtype; USE nrutil, ONLY : nrerror,reallocate
	USE arcode_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(INOUT) :: ich,lcd
	INTEGER(I4B), INTENT(IN) :: isign
	CHARACTER(1), DIMENSION(:), POINTER :: codep
	TYPE(arithcode) :: acode
	INTEGER(I4B) :: ihi,j,ja,jh,jl,m
	if (isign == 0) then
		acode%jdif=acode%nrad-1
		acode%ilob(:)=0
		acode%iupb(:)=acode%nrad-1
		do j=NWK,1,-1
			acode%nc=j
			if (acode%jdif > acode%minint) RETURN
			acode%jdif=(acode%jdif+1)*acode%nrad-1
		end do
		call nrerror('NWK too small in arcode')
	else
		if (isign > 0) then
			if (ich > acode%nch .or. ich < 0) call nrerror('bad ich in arcode')
		else
			ja=ichar(codep(lcd))-acode%ilob(acode%nc)
			do j=acode%nc+1,NWK
				ja=ja*acode%nrad+(ichar(codep(j+lcd-acode%nc))-acode%ilob(j))
			end do
			ich=0
			ihi=acode%nch+1
			do
				if (ihi-ich <= 1) exit
				m=(ich+ihi)/2
				if (ja >= jtry(acode%jdif,acode%ncumfq(m+1),acode%ncum)) then
					ich=m
				else
					ihi=m
				end if
			end do
			if (ich == acode%nch) RETURN
		end if
		jh=jtry(acode%jdif,acode%ncumfq(ich+2),acode%ncum)
		jl=jtry(acode%jdif,acode%ncumfq(ich+1),acode%ncum)
		acode%jdif=jh-jl
		call arcsum(acode%ilob,acode%iupb,jh,NWK,acode%nrad,acode%nc)
		call arcsum(acode%ilob,acode%ilob,jl,NWK,acode%nrad,acode%nc)
		do j=acode%nc,NWK
			if (ich /= acode%nch .and. acode%iupb(j) /= acode%ilob(j)) exit
			if (acode%nc > size(codep)) codep=>reallocate(codep,2*size(codep))
			if (isign > 0) codep(lcd)=char(acode%ilob(j))
			lcd=lcd+1
		end do
		if (j > NWK) RETURN
		acode%nc=j
		j=0
		do
			if (acode%jdif >= acode%minint) exit
			j=j+1
			acode%jdif=acode%jdif*acode%nrad
		end do
		if (acode%nc-j < 1) call nrerror('NWK too small in arcode')
		if (j /= 0) then
			acode%iupb((acode%nc-j):(NWK-j))=acode%iupb(acode%nc:NWK)
			acode%ilob((acode%nc-j):(NWK-j))=acode%ilob(acode%nc:NWK)
		end if
		acode%nc=acode%nc-j
		acode%iupb((NWK-j+1):NWK)=0
		acode%ilob((NWK-j+1):NWK)=0
	end if
	CONTAINS
!BL
	FUNCTION jtry(m,n,k)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: m,n,k
	INTEGER(I4B) :: jtry
	jtry=int((real(m,dp)*real(n,dp))/real(k,dp))
	END FUNCTION jtry
!BL
	SUBROUTINE arcsum(iin,iout,ja,nwk,nrad,nc)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iin
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: iout
	INTEGER(I4B), INTENT(IN) :: nwk,nrad,nc
	INTEGER(I4B), INTENT(INOUT) :: ja
	INTEGER(I4B) :: j,jtmp,karry
	karry=0
	do j=nwk,nc+1,-1
		jtmp=ja
		ja=ja/nrad
		iout(j)=iin(j)+(jtmp-ja*nrad)+karry
		if (iout(j) >= nrad) then
			iout(j)=iout(j)-nrad
			karry=1
		else
			karry=0
		end if
	end do
	iout(nc)=iin(nc)+ja+karry
	END SUBROUTINE arcsum
	END SUBROUTINE arcode
