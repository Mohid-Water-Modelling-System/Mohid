	SUBROUTINE anneal(x,y,iorder)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,swap
	USE nr, ONLY : ran1
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	INTEGER(I4B), DIMENSION(6) :: n
	INTEGER(I4B) :: i1,i2,j,k,nlimit,ncity,nn,nover,nsucc
	REAL(SP) :: de,harvest,path,t,tfactr
	LOGICAL(LGT) :: ans
	ncity=assert_eq(size(x),size(y),size(iorder),'anneal')
	nover=100*ncity
	nlimit=10*ncity
	tfactr=0.9_sp
	t=0.5_sp
	path=sum(alen_v(x(iorder(1:ncity-1)),x(iorder(2:ncity)),&
		y(iorder(1:ncity-1)),y(iorder(2:ncity))))
	i1=iorder(ncity)
	i2=iorder(1)
	path=path+alen(x(i1),x(i2),y(i1),y(i2))
	do j=1,100
		nsucc=0
		do k=1,nover
			do
				call ran1(harvest)
				n(1)=1+int(ncity*harvest)
				call ran1(harvest)
				n(2)=1+int((ncity-1)*harvest)
				if (n(2) >= n(1)) n(2)=n(2)+1
				nn=1+mod((n(1)-n(2)+ncity-1),ncity)
				if (nn >= 3) exit
			end do
			call ran1(harvest)
			if (harvest < 0.5_sp) then
				call ran1(harvest)
				n(3)=n(2)+int(abs(nn-2)*harvest)+1
				n(3)=1+mod(n(3)-1,ncity)
				call trncst(x,y,iorder,n,de)
				call metrop(de,t,ans)
				if (ans) then
					nsucc=nsucc+1
					path=path+de
					call trnspt(iorder,n)
				end if
			else
				call revcst(x,y,iorder,n,de)
				call metrop(de,t,ans)
				if (ans) then
					nsucc=nsucc+1
					path=path+de
					call revers(iorder,n)
				end if
			end if
			if (nsucc >= nlimit) exit
		end do
		write(*,*)
		write(*,*) 'T =',t,' Path Length =',path
		write(*,*) 'Successful Moves: ',nsucc
		t=t*tfactr
		if (nsucc == 0) RETURN
	end do
	CONTAINS
!BL
	FUNCTION alen(x1,x2,y1,y2)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,y1,y2
	REAL(SP) :: alen
	alen=sqrt((x2-x1)**2+(y2-y1)**2)
	END FUNCTION alen
!BL
	FUNCTION alen_v(x1,x2,y1,y2)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1,x2,y1,y2
	REAL(SP), DIMENSION(size(x1)) :: alen_v
	alen_v=sqrt((x2-x1)**2+(y2-y1)**2)
	END FUNCTION alen_v
!BL
	SUBROUTINE metrop(de,t,ans)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: de,t
	LOGICAL(LGT), INTENT(OUT) :: ans
	call ran1(harvest)
	ans=(de < 0.0) .or. (harvest < exp(-de/t))
	END SUBROUTINE metrop
!BL
	SUBROUTINE revcst(x,y,iorder,n,de)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: n
	REAL(SP), INTENT(OUT) :: de
	INTEGER(I4B) :: ncity
	REAL(SP), DIMENSION(4) :: xx,yy
	ncity=size(x)
	n(3)=1+mod((n(1)+ncity-2),ncity)
	n(4)=1+mod(n(2),ncity)
	xx(1:4)=x(iorder(n(1:4)))
	yy(1:4)=y(iorder(n(1:4)))
	de=-alen(xx(1),xx(3),yy(1),yy(3))&
		-alen(xx(2),xx(4),yy(2),yy(4))&
		+alen(xx(1),xx(4),yy(1),yy(4))&
		+alen(xx(2),xx(3),yy(2),yy(3))
	END SUBROUTINE revcst
!BL
	SUBROUTINE revers(iorder,n)
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	INTEGER(I4B) :: j,k,l,nn,ncity
	ncity=size(iorder)
	nn=(1+mod(n(2)-n(1)+ncity,ncity))/2
	do j=1,nn
		k=1+mod((n(1)+j-2),ncity)
		l=1+mod((n(2)-j+ncity),ncity)
		call swap(iorder(k),iorder(l))
	end do
	END SUBROUTINE revers
!BL
	SUBROUTINE trncst(x,y,iorder,n,de)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: n
	REAL(SP), INTENT(OUT) :: de
	INTEGER(I4B) :: ncity
	REAL(SP), DIMENSION(6) :: xx,yy
	ncity=size(x)
	n(4)=1+mod(n(3),ncity)
	n(5)=1+mod((n(1)+ncity-2),ncity)
	n(6)=1+mod(n(2),ncity)
	xx(1:6)=x(iorder(n(1:6)))
	yy(1:6)=y(iorder(n(1:6)))
	de=-alen(xx(2),xx(6),yy(2),yy(6))&
		-alen(xx(1),xx(5),yy(1),yy(5))&
		-alen(xx(3),xx(4),yy(3),yy(4))&
		+alen(xx(1),xx(3),yy(1),yy(3))&
		+alen(xx(2),xx(4),yy(2),yy(4))&
		+alen(xx(5),xx(6),yy(5),yy(6))
	END SUBROUTINE trncst
!BL
	SUBROUTINE trnspt(iorder,n)
	IMPLICIT NONE
	INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: iorder
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: n
	INTEGER(I4B) :: m1,m2,m3,nn,ncity
	INTEGER(I4B), DIMENSION(size(iorder)) :: jorder
	ncity=size(iorder)
	m1=1+mod((n(2)-n(1)+ncity),ncity)
	m2=1+mod((n(5)-n(4)+ncity),ncity)
	m3=1+mod((n(3)-n(6)+ncity),ncity)
	jorder(1:m1)=iorder(1+mod((arth(1,1,m1)+n(1)-2),ncity))
	nn=m1
	jorder(nn+1:nn+m2)=iorder(1+mod((arth(1,1,m2)+n(4)-2),ncity))
	nn=nn+m2
	jorder(nn+1:nn+m3)=iorder(1+mod((arth(1,1,m3)+n(6)-2),ncity))
	iorder(1:ncity)=jorder(1:ncity)
	END SUBROUTINE trnspt
	END SUBROUTINE anneal
