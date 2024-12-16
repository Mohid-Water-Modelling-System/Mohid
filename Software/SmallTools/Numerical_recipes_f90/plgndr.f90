	FUNCTION plgndr_s(l,m,x)
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: l,m
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: plgndr_s
	INTEGER(I4B) :: ll
	REAL(SP) :: pll,pmm,pmmp1,somx2
	call assert(m >= 0, m <= l, abs(x) <= 1.0, 'plgndr_s args')
	pmm=1.0
	if (m > 0) then
		somx2=sqrt((1.0_sp-x)*(1.0_sp+x))
		pmm=product(arth(1.0_sp,2.0_sp,m))*somx2**m
		if (mod(m,2) == 1) pmm=-pmm
	end if
	if (l == m) then
		plgndr_s=pmm
	else
		pmmp1=x*(2*m+1)*pmm
		if (l == m+1) then
			plgndr_s=pmmp1
		else
			do ll=m+2,l
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
				pmm=pmmp1
				pmmp1=pll
			end do
			plgndr_s=pll
		end if
	end if
	END FUNCTION plgndr_s


	FUNCTION plgndr_v(l,m,x)
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: l,m
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: plgndr_v
	INTEGER(I4B) :: ll
	REAL(SP), DIMENSION(size(x)) :: pll,pmm,pmmp1,somx2
	call assert(m >= 0, m <= l, all(abs(x) <= 1.0), 'plgndr_v args')
	pmm=1.0
	if (m > 0) then
		somx2=sqrt((1.0_sp-x)*(1.0_sp+x))
		pmm=product(arth(1.0_sp,2.0_sp,m))*somx2**m
		if (mod(m,2) == 1) pmm=-pmm
	end if
	if (l == m) then
		plgndr_v=pmm
	else
		pmmp1=x*(2*m+1)*pmm
		if (l == m+1) then
			plgndr_v=pmmp1
		else
			do ll=m+2,l
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
				pmm=pmmp1
				pmmp1=pll
			end do
			plgndr_v=pll
		end if
	end if
	END FUNCTION plgndr_v
