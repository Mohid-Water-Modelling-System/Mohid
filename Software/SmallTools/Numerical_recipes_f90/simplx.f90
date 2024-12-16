	SUBROUTINE simplx(a,m1,m2,m3,icase,izrov,iposv)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,ifirstloc,imaxloc,&
		nrerror,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2,m3
	INTEGER(I4B), INTENT(OUT) :: icase
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: izrov,iposv
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	INTEGER(I4B) :: ip,k,kh,kp,nl1,m,n
	INTEGER(I4B), DIMENSION(size(a,2)) :: l1
	INTEGER(I4B), DIMENSION(m2) :: l3
	REAL(SP) :: bmax
	LOGICAL(LGT) :: init
	m=assert_eq(size(a,1)-2,size(iposv),'simplx: m')
	n=assert_eq(size(a,2)-1,size(izrov),'simplx: n')
	if (m /= m1+m2+m3) call nrerror('simplx: bad input constraint counts')
	if (any(a(2:m+1,1) < 0.0)) call nrerror('bad input tableau in simplx')
	nl1=n
	l1(1:n)=arth(1,1,n)
	izrov(:)=l1(1:n)
	iposv(:)=n+arth(1,1,m)
	init=.true.
	phase1: do
		if (init) then
			if (m2+m3 == 0) exit phase1
			init=.false.
			l3(1:m2)=1
			a(m+2,1:n+1)=-sum(a(m1+2:m+1,1:n+1),dim=1)
		end if
		if (nl1 > 0) then
			kp=l1(imaxloc(a(m+2,l1(1:nl1)+1)))
			bmax=a(m+2,kp+1)
		else
			bmax=0.0
		end if
		phase1a: do
			if (bmax <= EPS .and. a(m+2,1) < -EPS) then
				icase=-1
				RETURN
			else if (bmax <= EPS .and. a(m+2,1) <= EPS) then
				do ip=m1+m2+1,m
					if (iposv(ip) == ip+n) then
						if (nl1 > 0) then
							kp=l1(imaxloc(abs(a(ip+1,l1(1:nl1)+1))))
							bmax=a(ip+1,kp+1)
						else
							bmax=0.0
						end if
						if (bmax > EPS) exit phase1a
					end if
				end do
				where (spread(l3(1:m2),2,n+1) == 1) &
					a(m1+2:m1+m2+1,1:n+1)=-a(m1+2:m1+m2+1,1:n+1)
				exit phase1
			end if
			call simp1
			if (ip == 0) then
				icase=-1
				RETURN
			end if
			exit phase1a
		end do phase1a
		call simp2(m+1,n)
		if (iposv(ip) >= n+m1+m2+1) then
			k=ifirstloc(l1(1:nl1) == kp)
			nl1=nl1-1
			l1(k:nl1)=l1(k+1:nl1+1)
		else
			kh=iposv(ip)-m1-n
			if (kh >= 1) then
				if (l3(kh) /= 0) then
					l3(kh)=0
					a(m+2,kp+1)=a(m+2,kp+1)+1.0_sp
					a(1:m+2,kp+1)=-a(1:m+2,kp+1)
				end if
			end if
		end if
		call swap(izrov(kp),iposv(ip))
	end do phase1
	phase2: do
		if (nl1 > 0) then
			kp=l1(imaxloc(a(1,l1(1:nl1)+1)))
			bmax=a(1,kp+1)
		else
			bmax=0.0
		end if
		if (bmax <= EPS) then
			icase=0
			RETURN
		end if
		call simp1
		if (ip == 0) then
			icase=1
			RETURN
		end if
		call simp2(m,n)
		call swap(izrov(kp),iposv(ip))
	end do phase2
	CONTAINS
!BL
	SUBROUTINE simp1
	IMPLICIT NONE
	INTEGER(I4B) :: i,k
	REAL(SP) :: q,q0,q1,qp
	ip=0
	i=ifirstloc(a(2:m+1,kp+1) < -EPS)
	if (i > m) RETURN
	q1=-a(i+1,1)/a(i+1,kp+1)
	ip=i
	do i=ip+1,m
		if (a(i+1,kp+1) < -EPS) then
			q=-a(i+1,1)/a(i+1,kp+1)
			if (q < q1) then
				ip=i
				q1=q
			else if (q == q1) then
				do k=1,n
					qp=-a(ip+1,k+1)/a(ip+1,kp+1)
					q0=-a(i+1,k+1)/a(i+1,kp+1)
					if (q0 /= qp) exit
				end do
				if (q0 < qp) ip=i
			end if
		end if
	end do
	END SUBROUTINE simp1
!BL
	SUBROUTINE simp2(i1,k1)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: i1,k1
	INTEGER(I4B) :: ip1,kp1
	REAL(SP) :: piv
	INTEGER(I4B), DIMENSION(k1) :: icol
	INTEGER(I4B), DIMENSION(i1) :: irow
	INTEGER(I4B), DIMENSION(max(i1,k1)+1) :: itmp
	ip1=ip+1
	kp1=kp+1
	piv=1.0_sp/a(ip1,kp1)
	itmp(1:k1+1)=arth(1,1,k1+1)
	icol=pack(itmp(1:k1+1),itmp(1:k1+1) /= kp1)
	itmp(1:i1+1)=arth(1,1,i1+1)
	irow=pack(itmp(1:i1+1),itmp(1:i1+1) /= ip1)
	a(irow,kp1)=a(irow,kp1)*piv
	a(irow,icol)=a(irow,icol)-outerprod(a(irow,kp1),a(ip1,icol))
	a(ip1,icol)=-a(ip1,icol)*piv
	a(ip1,kp1)=piv
	END SUBROUTINE simp2
	END SUBROUTINE simplx
