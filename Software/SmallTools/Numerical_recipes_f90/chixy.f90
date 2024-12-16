MODULE chixyfit
	USE nrtype; USE nrutil, ONLY : nrerror
	REAL(SP), DIMENSION(:), POINTER :: xxp,yyp,sxp,syp,wwp
	REAL(SP) :: aa,offs
CONTAINS
!BL
	FUNCTION chixy(bang)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: bang
	REAL(SP) :: chixy
	REAL(SP), PARAMETER :: BIG=1.0e30_sp
	REAL(SP) :: avex,avey,sumw,b
	if (.not. associated(wwp)) call nrerror("chixy: bad pointers")
	b=tan(bang)
	wwp(:)=(b*sxp(:))**2+syp(:)**2
	where (wwp(:) < 1.0/BIG)
		wwp(:)=BIG
	elsewhere
		wwp(:)=1.0_sp/wwp(:)
	end where
	sumw=sum(wwp)
	avex=dot_product(wwp,xxp)/sumw
	avey=dot_product(wwp,yyp)/sumw
	aa=avey-b*avex
	chixy=sum(wwp(:)*(yyp(:)-aa-b*xxp(:))**2)-offs
	END FUNCTION chixy
END MODULE chixyfit
