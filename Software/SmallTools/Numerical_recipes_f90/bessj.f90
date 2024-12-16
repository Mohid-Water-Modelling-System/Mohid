	FUNCTION bessj_s(n,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessj0,bessj1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: bessj_s
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
	INTEGER(I4B) :: j,jsum,m
	REAL(SP) :: ax,bj,bjm,bjp,summ,tox
	call assert(n >= 2, 'bessj_s args')
	ax=abs(x)
	if (ax*ax <= 8.0_sp*tiny(x)) then
		bessj_s=0.0
	else if (ax > real(n,sp)) then
		tox=2.0_sp/ax
		bjm=bessj0(ax)
		bj=bessj1(ax)
		do j=1,n-1
			bjp=j*tox*bj-bjm
			bjm=bj
			bj=bjp
		end do
		bessj_s=bj
	else
		tox=2.0_sp/ax
		m=2*((n+int(sqrt(real(IACC*n,sp))))/2)
		bessj_s=0.0
		jsum=0
		summ=0.0
		bjp=0.0
		bj=1.0
		do j=m,1,-1
			bjm=j*tox*bj-bjp
			bjp=bj
			bj=bjm
			if (exponent(bj) > IEXP) then
				bj=scale(bj,-IEXP)
				bjp=scale(bjp,-IEXP)
				bessj_s=scale(bessj_s,-IEXP)
				summ=scale(summ,-IEXP)
			end if
			if (jsum /= 0) summ=summ+bj
			jsum=1-jsum
			if (j == n) bessj_s=bjp
		end do
		summ=2.0_sp*summ-bj
		bessj_s=bessj_s/summ
	end if
	if (x < 0.0 .and. mod(n,2) == 1) bessj_s=-bessj_s
	END FUNCTION bessj_s

	FUNCTION bessj_v(n,xx)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : bessj0,bessj1
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	REAL(SP), DIMENSION(size(xx)) :: bessj_v
	INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(xx)/2
	REAL(SP), DIMENSION(size(xx)) :: ax
	LOGICAL(LGT), DIMENSION(size(xx)) :: mask,mask0
	REAL(SP), DIMENSION(:), ALLOCATABLE :: x,bj,bjm,bjp,summ,tox,bessjle
	LOGICAL(LGT), DIMENSION(:), ALLOCATABLE :: renorm
	INTEGER(I4B) :: j,jsum,m,npak
	call assert(n >= 2, 'bessj_v args')
	ax=abs(xx)
	mask = (ax <= real(n,sp))
	mask0 = (ax*ax <= 8.0_sp*tiny(xx))
	bessj_v=bessjle_v(n,ax,logical(mask .and. .not.mask0, kind=lgt))
	bessj_v=merge(bessjgt_v(n,ax,.not. mask),bessj_v,.not. mask)
	where (mask0) bessj_v=0.0
	where (xx < 0.0 .and. mod(n,2) == 1) bessj_v=-bessj_v
	CONTAINS
!BL
!BL
	FUNCTION bessjgt_v(n,xx,mask)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	LOGICAL(LGT), DIMENSION(size(xx)), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(xx)) :: bessjgt_v
	npak=count(mask)
	if (npak == 0) RETURN
	allocate(x(npak),bj(npak),bjm(npak),bjp(npak),tox(npak))
	x=pack(xx,mask)
	tox=2.0_sp/x
	bjm=bessj0(x)
	bj=bessj1(x)
	do j=1,n-1
		bjp=j*tox*bj-bjm
		bjm=bj
		bj=bjp
	end do
	bessjgt_v=unpack(bj,mask,0.0_sp)
	deallocate(x,bj,bjm,bjp,tox)
	END FUNCTION bessjgt_v
!BL
	FUNCTION bessjle_v(n,xx,mask)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP), DIMENSION(:), INTENT(IN) :: xx
	LOGICAL(LGT), DIMENSION(size(xx)), INTENT(IN) :: mask
	REAL(SP), DIMENSION(size(xx)) :: bessjle_v
	npak=count(mask)
	if (npak == 0) RETURN
	allocate(x(npak),bj(npak),bjm(npak),bjp(npak),summ(npak), &
		bessjle(npak),tox(npak),renorm(npak))
	x=pack(xx,mask)
	tox=2.0_sp/x
	m=2*((n+int(sqrt(real(IACC*n,sp))))/2)
	bessjle=0.0
	jsum=0
	summ=0.0
	bjp=0.0
	bj=1.0
	do j=m,1,-1
		bjm=j*tox*bj-bjp
		bjp=bj
		bj=bjm
		renorm = (exponent(bj)>IEXP)
		bj=merge(scale(bj,-IEXP),bj,renorm)
		bjp=merge(scale(bjp,-IEXP),bjp,renorm)
		bessjle=merge(scale(bessjle,-IEXP),bessjle,renorm)
		summ=merge(scale(summ,-IEXP),summ,renorm)
		if (jsum /= 0) summ=summ+bj
		jsum=1-jsum
		if (j == n) bessjle=bjp
	end do
	summ=2.0_sp*summ-bj
	bessjle=bessjle/summ
	bessjle_v=unpack(bessjle,mask,0.0_sp)
	deallocate(x,bj,bjm,bjp,summ,bessjle,tox,renorm)
	END FUNCTION bessjle_v
	END FUNCTION bessj_v
