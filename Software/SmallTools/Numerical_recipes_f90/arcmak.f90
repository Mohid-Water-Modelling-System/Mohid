MODULE arcode_info
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NWK=20
	TYPE arithcode
		INTEGER(I4B), DIMENSION(:), POINTER :: ilob,iupb,ncumfq
		INTEGER(I4B) :: jdif,nc,minint,nch,ncum,nrad
	END TYPE arithcode
CONTAINS
	SUBROUTINE arcode_allocate(acode,mc)
	USE nrtype
	IMPLICIT NONE
	TYPE(arithcode) :: acode
	INTEGER(I4B) :: mc
	allocate(acode%ilob(NWK),acode%iupb(NWK),acode%ncumfq(mc+2))
	END SUBROUTINE arcode_allocate
!BL
	SUBROUTINE arcode_deallocate(acode)
	USE nrtype
	IMPLICIT NONE
	TYPE(arithcode) :: acode
	deallocate(acode%ncumfq,acode%iupb,acode%ilob)
	nullify(acode%ilob)
	nullify(acode%iupb)
	nullify(acode%ncumfq)
	END SUBROUTINE arcode_deallocate
END MODULE arcode_info

	SUBROUTINE arcmak(nfreq,nradd,acode)
	USE nrtype; USE nrutil, ONLY : cumsum,nrerror
	USE arcode_info
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: nradd
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nfreq
	TYPE(arithcode) :: acode
	INTEGER(I4B), PARAMETER :: MAXINT=huge(nradd)
	if (nradd > 256) call nrerror('output radix may not exceed 256 in arcmak')
	acode%minint=MAXINT/nradd
	acode%nch=size(nfreq)
	acode%nrad=nradd
	call arcode_allocate(acode,acode%nch)
	acode%ncumfq(1)=0
	acode%ncumfq(2:acode%nch+1)=cumsum(max(nfreq(1:acode%nch),1))
	acode%ncumfq(acode%nch+2)=acode%ncumfq(acode%nch+1)+1
	acode%ncum=acode%ncumfq(acode%nch+2)
	END SUBROUTINE arcmak
