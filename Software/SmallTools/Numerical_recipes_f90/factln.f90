	FUNCTION factln_s(n)
	USE nrtype; USE nrutil, ONLY : arth,assert
	USE nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: factln_s
	INTEGER(I4B), PARAMETER :: TMAX=100
	REAL(SP), DIMENSION(TMAX), SAVE :: a
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		a(1:TMAX)=gammln(arth(1.0_sp,1.0_sp,TMAX))
		init=.false.
	end if
	call assert(n >= 0, 'factln_s arg')
	if (n < TMAX) then
		factln_s=a(n+1)
	else
		factln_s=gammln(n+1.0_sp)
	end if
	END FUNCTION factln_s


	FUNCTION factln_v(n)
	USE nrtype; USE nrutil, ONLY : arth,assert
	USE nr, ONLY : gammln
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	REAL(SP), DIMENSION(size(n)) :: factln_v
	LOGICAL(LGT), DIMENSION(size(n)) :: mask
	INTEGER(I4B), PARAMETER :: TMAX=100
	REAL(SP), DIMENSION(TMAX), SAVE :: a
	LOGICAL(LGT), SAVE :: init=.true.
	if (init) then
		a(1:TMAX)=gammln(arth(1.0_sp,1.0_sp,TMAX))
		init=.false.
	end if
	call assert(all(n >= 0), 'factln_v arg')
	mask = (n >= TMAX)
	factln_v=unpack(gammln(pack(n,mask)+1.0_sp),mask,0.0_sp)
	where (.not. mask) factln_v=a(n+1)
	END FUNCTION factln_v
