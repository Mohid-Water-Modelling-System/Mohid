	RECURSIVE SUBROUTINE miser(func,regn,ndim,npts,dith,ave,var)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP) :: func
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	REAL(SP), DIMENSION(:), INTENT(IN) :: regn
	INTEGER(I4B), INTENT(IN) :: ndim,npts
	REAL(SP), INTENT(IN) :: dith
	REAL(SP), INTENT(OUT) :: ave,var
	REAL(SP), PARAMETER :: PFAC=0.1_sp,TINY=1.0e-30_sp,BIG=1.0e30_sp
	INTEGER(I4B), PARAMETER :: MNPT=15,MNBS=60
	REAL(SP), DIMENSION(:), ALLOCATABLE :: regn_temp
	INTEGER(I4B) :: j,jb,n,ndum,npre,nptl,nptr
	INTEGER(I4B), SAVE :: iran=0
	REAL(SP) :: avel,varl,fracl,fval,rgl,rgm,rgr,&
		s,sigl,siglb,sigr,sigrb,sm,sm2,sumb,sumr
	REAL(SP), DIMENSION(:), ALLOCATABLE :: fmaxl,fmaxr,fminl,fminr,pt,rmid
	ndum=assert_eq(size(regn),2*ndim,'miser')
	allocate(pt(ndim))
	if (npts < MNBS) then
		sm=0.0
		sm2=0.0
		do n=1,npts
			call ranpt(pt,regn)
			fval=func(pt)
			sm=sm+fval
			sm2=sm2+fval**2
		end do
		ave=sm/npts
		var=max(TINY,(sm2-sm**2/npts)/npts**2)
	else
		npre=max(int(npts*PFAC),MNPT)
		allocate(rmid(ndim),fmaxl(ndim),fmaxr(ndim),fminl(ndim),fminr(ndim))
		fminl(:)=BIG
		fminr(:)=BIG
		fmaxl(:)=-BIG
		fmaxr(:)=-BIG
		do j=1,ndim
			iran=mod(iran*2661+36979,175000)
			s=sign(dith,real(iran-87500,sp))
			rmid(j)=(0.5_sp+s)*regn(j)+(0.5_sp-s)*regn(ndim+j)
		end do
		do n=1,npre
			call ranpt(pt,regn)
			fval=func(pt)
			where (pt <= rmid)
				fminl=min(fminl,fval)
				fmaxl=max(fmaxl,fval)
			elsewhere
				fminr=min(fminr,fval)
				fmaxr=max(fmaxr,fval)
			end where
		end do
		sumb=BIG
		jb=0
		siglb=1.0
		sigrb=1.0
		do j=1,ndim
			if (fmaxl(j) > fminl(j) .and. fmaxr(j) > fminr(j)) then
				sigl=max(TINY,(fmaxl(j)-fminl(j))**(2.0_sp/3.0_sp))
				sigr=max(TINY,(fmaxr(j)-fminr(j))**(2.0_sp/3.0_sp))
				sumr=sigl+sigr
				if (sumr <= sumb) then
					sumb=sumr
					jb=j
					siglb=sigl
					sigrb=sigr
				end if
			end if
		end do
		deallocate(fminr,fminl,fmaxr,fmaxl)
		if (jb == 0) jb=1+(ndim*iran)/175000
		rgl=regn(jb)
		rgm=rmid(jb)
		rgr=regn(ndim+jb)
		fracl=abs((rgm-rgl)/(rgr-rgl))
		nptl=(MNPT+(npts-npre-2*MNPT)*fracl*siglb/ &
			(fracl*siglb+(1.0_sp-fracl)*sigrb))
		nptr=npts-npre-nptl
		allocate(regn_temp(2*ndim))
		regn_temp(:)=regn(:)
		regn_temp(ndim+jb)=rmid(jb)
		call miser(func,regn_temp,ndim,nptl,dith,avel,varl)
		regn_temp(jb)=rmid(jb)
		regn_temp(ndim+jb)=regn(ndim+jb)
		call miser(func,regn_temp,ndim,nptr,dith,ave,var)
		deallocate(regn_temp)
		ave=fracl*avel+(1-fracl)*ave
		var=fracl*fracl*varl+(1-fracl)*(1-fracl)*var
		deallocate(rmid)
	end if
	deallocate(pt)
	CONTAINS
!BL
	SUBROUTINE ranpt(pt,region)
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(OUT) :: pt
	REAL(SP), DIMENSION(:), INTENT(IN) :: region
	INTEGER(I4B) :: n
	call ran1(pt)
	n=size(pt)
	pt(1:n)=region(1:n)+(region(n+1:2*n)-region(1:n))*pt(1:n)
	END SUBROUTINE ranpt
	END SUBROUTINE miser
