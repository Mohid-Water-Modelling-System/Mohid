	SUBROUTINE atimes(x,r,itrnsp)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : sprsax,sprstx
	USE xlinbcg_data
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(:), INTENT(OUT) :: r
	INTEGER(I4B), INTENT(IN) :: itrnsp
	INTEGER(I4B) :: n
	n=assert_eq(size(x),size(r),'atimes')
	if (itrnsp == 0) then
		call sprsax(sa,x,r)
	else
		call sprstx(sa,x,r)
	end if
	END SUBROUTINE atimes
