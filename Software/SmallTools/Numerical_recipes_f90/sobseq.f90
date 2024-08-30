	SUBROUTINE sobseq(x,init)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
	INTEGER(I4B), PARAMETER :: MAXBIT=30,MAXDIM=6
	REAL(SP), SAVE :: fac
	INTEGER(I4B) :: i,im,ipp,j,k,l
	INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE:: iu
	INTEGER(I4B), SAVE :: in
	INTEGER(I4B), DIMENSION(MAXDIM), SAVE :: ip,ix,mdeg
	INTEGER(I4B), DIMENSION(MAXDIM*MAXBIT), SAVE :: iv
	DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
	DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
	if (present(init)) then
		ix=0
		in=0
		if (iv(1) /= 1) RETURN
		fac=1.0_sp/2.0_sp**MAXBIT
		allocate(iu(MAXDIM,MAXBIT))
		iu=reshape(iv,shape(iu))
		do k=1,MAXDIM
			do j=1,mdeg(k)
				iu(k,j)=iu(k,j)*2**(MAXBIT-j)
			end do
			do j=mdeg(k)+1,MAXBIT
				ipp=ip(k)
				i=iu(k,j-mdeg(k))
				i=ieor(i,i/2**mdeg(k))
				do l=mdeg(k)-1,1,-1
					if (btest(ipp,0)) i=ieor(i,iu(k,j-l))
					ipp=ipp/2
				end do
				iu(k,j)=i
			end do
		end do
		iv=reshape(iu,shape(iv))
		deallocate(iu)
	else
		im=in
		do j=1,MAXBIT
			if (.not. btest(im,0)) exit
			im=im/2
		end do
		if (j > MAXBIT) call nrerror('MAXBIT too small in sobseq')
		im=(j-1)*MAXDIM
		j=min(size(x),MAXDIM)
		ix(1:j)=ieor(ix(1:j),iv(1+im:j+im))
		x(1:j)=ix(1:j)*fac
		in=in+1
	end if
	END SUBROUTINE sobseq
