	SUBROUTINE wtn(a,nn,isign,wtstep)
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nn
	INTEGER(I4B), INTENT(IN) :: isign
	INTERFACE
		SUBROUTINE wtstep(a,isign)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
		INTEGER(I4B), INTENT(IN) :: isign
		END SUBROUTINE wtstep
	END INTERFACE
	INTEGER(I4B) :: i1,i2,i3,idim,n,ndim,nnew,nprev,nt,ntot
	REAL(SP), DIMENSION(:), ALLOCATABLE :: wksp
	call assert(iand(nn,nn-1)==0, 'each dimension must be a power of 2 in wtn')
	allocate(wksp(maxval(nn)))
	ndim=size(nn)
	ntot=product(nn(:))
	nprev=1
	do idim=1,ndim
		n=nn(idim)
		nnew=n*nprev
		if (n > 4) then
			do i2=0,ntot-1,nnew
				do i1=1,nprev
					i3=i1+i2
					wksp(1:n)=a(arth(i3,nprev,n))
					i3=i3+n*nprev
					if (isign >= 0) then
						nt=n
						do
							if (nt < 4) exit
							call wtstep(wksp(1:nt),isign)
							nt=nt/2
						end do
					else
						nt=4
						do
							if (nt > n) exit
							call wtstep(wksp(1:nt),isign)
							nt=nt*2
						end do
					end if
					i3=i1+i2
					a(arth(i3,nprev,n))=wksp(1:n)
					i3=i3+n*nprev
				end do
			end do
		end if
		nprev=nnew
	end do
	deallocate(wksp)
	END SUBROUTINE wtn
