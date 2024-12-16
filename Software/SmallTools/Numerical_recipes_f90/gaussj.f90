	SUBROUTINE gaussj(a,b)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
	INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
	LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
	REAL(SP) :: pivinv
	REAL(SP), DIMENSION(size(a,1)) :: dumc
	INTEGER(I4B), TARGET :: irc(2)
	INTEGER(I4B) :: i,l,n
	INTEGER(I4B), POINTER :: irow,icol
	n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
	irow => irc(1)
	icol => irc(2)
	ipiv=0
	do i=1,n
		lpiv = (ipiv == 0)
		irc=maxloc(abs(a),outerand(lpiv,lpiv))
		ipiv(icol)=ipiv(icol)+1
		if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
		if (irow /= icol) then
			call swap(a(irow,:),a(icol,:))
			call swap(b(irow,:),b(icol,:))
		end if
		indxr(i)=irow
		indxc(i)=icol
		if (a(icol,icol) == 0.0) &
			call nrerror('gaussj: singular matrix (2)')
		pivinv=1.0_sp/a(icol,icol)
		a(icol,icol)=1.0
		a(icol,:)=a(icol,:)*pivinv
		b(icol,:)=b(icol,:)*pivinv
		dumc=a(:,icol)
		a(:,icol)=0.0
		a(icol,icol)=pivinv
		a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
		b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
		a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
		b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
	end do
	do l=n,1,-1
		call swap(a(:,indxr(l)),a(:,indxc(l)))
	end do
	END SUBROUTINE gaussj
