	SUBROUTINE svbksb_sp(u,w,v,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: u,v
	REAL(SP), DIMENSION(:), INTENT(IN) :: w,b
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B) :: mdum,ndum
	REAL(SP), DIMENSION(size(x)) :: tmp
	mdum=assert_eq(size(u,1),size(b),'svbksb_sp: mdum')
	ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),&
		'svbksb_sp: ndum')
	where (w /= 0.0)
		tmp=matmul(b,u)/w
	elsewhere
		tmp=0.0
	end where
	x=matmul(v,tmp)
	END SUBROUTINE svbksb_sp

	SUBROUTINE svbksb_dp(u,w,v,b,x)
	USE nrtype; USE nrutil, ONLY : assert_eq
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: u,v
	REAL(DP), DIMENSION(:), INTENT(IN) :: w,b
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B) :: mdum,ndum
	REAL(DP), DIMENSION(size(x)) :: tmp
	mdum=assert_eq(size(u,1),size(b),'svbksb_dp: mdum')
	ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),&
		'svbksb_dp: ndum')
	where (w /= 0.0)
		tmp=matmul(b,u)/w
	elsewhere
		tmp=0.0
	end where
	x=matmul(v,tmp)
	END SUBROUTINE svbksb_dp
