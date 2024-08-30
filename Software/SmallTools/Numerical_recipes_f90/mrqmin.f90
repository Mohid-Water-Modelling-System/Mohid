	SUBROUTINE mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
	USE nrtype; USE nrutil, ONLY : assert_eq,diagmult
	USE nr, ONLY : covsrt,gaussj
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
	REAL(SP), INTENT(OUT) :: chisq
	REAL(SP), INTENT(INOUT) :: alamda
	LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
	INTERFACE
		SUBROUTINE funcs(x,a,yfit,dyda)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x,a
		REAL(SP), DIMENSION(:), INTENT(OUT) :: yfit
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: dyda
		END SUBROUTINE funcs
	END INTERFACE
	INTEGER(I4B) :: ma,ndata
	INTEGER(I4B), SAVE :: mfit
	call mrqmin_private
	CONTAINS
!BL
	SUBROUTINE mrqmin_private
	REAL(SP), SAVE :: ochisq
	REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: atry,beta
	REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: da
	ndata=assert_eq(size(x),size(y),size(sig),'mrqmin: ndata')
	ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2),&
		size(alpha,1),size(alpha,2)/),'mrqmin: ma')
	mfit=count(maska)
	if (alamda < 0.0) then
		allocate(atry(ma),beta(ma),da(ma,1))
		alamda=0.001_sp
		call mrqcof(a,alpha,beta)
		ochisq=chisq
		atry=a
	end if
	covar(1:mfit,1:mfit)=alpha(1:mfit,1:mfit)
	call diagmult(covar(1:mfit,1:mfit),1.0_sp+alamda)
	da(1:mfit,1)=beta(1:mfit)
	call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1))
	if (alamda == 0.0) then
		call covsrt(covar,maska)
		call covsrt(alpha,maska)
		deallocate(atry,beta,da)
		RETURN
	end if
	atry=a+unpack(da(1:mfit,1),maska,0.0_sp)
	call mrqcof(atry,covar,da(1:mfit,1))
	if (chisq < ochisq) then
		alamda=0.1_sp*alamda
		ochisq=chisq
		alpha(1:mfit,1:mfit)=covar(1:mfit,1:mfit)
		beta(1:mfit)=da(1:mfit,1)
		a=atry
	else
		alamda=10.0_sp*alamda
		chisq=ochisq
	end if
	END SUBROUTINE mrqmin_private
!BL
	SUBROUTINE mrqcof(a,alpha,beta)
	REAL(SP), DIMENSION(:), INTENT(IN) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
	INTEGER(I4B) :: j,k,l,m
	REAL(SP), DIMENSION(size(x),size(a)) :: dyda
	REAL(SP), DIMENSION(size(x)) :: dy,sig2i,wt,ymod
	call funcs(x,a,ymod,dyda)
	sig2i=1.0_sp/(sig**2)
	dy=y-ymod
	j=0
	do l=1,ma
		if (maska(l)) then
			j=j+1
			wt=dyda(:,l)*sig2i
			k=0
			do m=1,l
				if (maska(m)) then
					k=k+1
					alpha(j,k)=dot_product(wt,dyda(:,m))
					alpha(k,j)=alpha(j,k)
				end if
			end do
			beta(j)=dot_product(dy,wt)
		end if
	end do
	chisq=dot_product(dy**2,sig2i)
	END SUBROUTINE mrqcof
	END SUBROUTINE mrqmin
