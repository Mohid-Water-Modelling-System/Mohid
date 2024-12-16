	SUBROUTINE lfit(x,y,sig,a,maska,covar,chisq,funcs)
	USE nrtype; USE nrutil, ONLY : assert_eq,diagmult,nrerror
	USE nr, ONLY :covsrt,gaussj
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
	REAL(SP), INTENT(OUT) :: chisq
	INTERFACE
		SUBROUTINE funcs(x,arr)
		USE nrtype
		IMPLICIT NONE
		REAL(SP),INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
		END SUBROUTINE funcs
	END INTERFACE
	INTEGER(I4B) :: i,j,k,l,ma,mfit,n
	REAL(SP) :: sig2i,wt,ym
	REAL(SP), DIMENSION(size(maska)) :: afunc
	REAL(SP), DIMENSION(size(maska),1) :: beta
	n=assert_eq(size(x),size(y),size(sig),'lfit: n')
	ma=assert_eq(size(maska),size(a),size(covar,1),size(covar,2),'lfit: ma')
	mfit=count(maska)
	if (mfit == 0) call nrerror('lfit: no parameters to be fitted')
	covar(1:mfit,1:mfit)=0.0
	beta(1:mfit,1)=0.0
	do i=1,n
		call funcs(x(i),afunc)
		ym=y(i)
		if (mfit < ma) ym=ym-sum(a(1:ma)*afunc(1:ma), mask=.not. maska)
		sig2i=1.0_sp/sig(i)**2
		j=0
		do l=1,ma
			if (maska(l)) then
				j=j+1
				wt=afunc(l)*sig2i
				k=count(maska(1:l))
				covar(j,1:k)=covar(j,1:k)+wt*pack(afunc(1:l),maska(1:l))
				beta(j,1)=beta(j,1)+ym*wt
			end if
		end do
	end do
	call diagmult(covar(1:mfit,1:mfit),0.5_sp)
	covar(1:mfit,1:mfit)= &
		covar(1:mfit,1:mfit)+transpose(covar(1:mfit,1:mfit))
	call gaussj(covar(1:mfit,1:mfit),beta(1:mfit,1:1))
	a(1:ma)=unpack(beta(1:ma,1),maska,a(1:ma))
	chisq=0.0
	do i=1,n
		call funcs(x(i),afunc)
		chisq=chisq+((y(i)-dot_product(a(1:ma),afunc(1:ma)))/sig(i))**2
	end do
	call covsrt(covar,maska)
	END SUBROUTINE lfit
