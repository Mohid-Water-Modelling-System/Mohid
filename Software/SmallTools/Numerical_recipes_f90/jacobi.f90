	SUBROUTINE jacobi(a,d,v,nrot)
	USE nrtype; USE nrutil, ONLY : assert_eq,get_diag,nrerror,unit_matrix,&
		upper_triangle
	IMPLICIT NONE
	INTEGER(I4B), INTENT(OUT) :: nrot
	REAL(SP), DIMENSION(:), INTENT(OUT) :: d
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: v
	INTEGER(I4B) :: i,ip,iq,n
	REAL(SP) :: c,g,h,s,sm,t,tau,theta,tresh
	REAL(SP), DIMENSION(size(d)) :: b,z
	n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),'jacobi')
	call unit_matrix(v(:,:))
	b(:)=get_diag(a(:,:))
	d(:)=b(:)
	z(:)=0.0
	nrot=0
	do i=1,50
		sm=sum(abs(a),mask=upper_triangle(n,n))
		if (sm == 0.0) RETURN
		tresh=merge(0.2_sp*sm/n**2,0.0_sp, i < 4 )
		do ip=1,n-1
			do iq=ip+1,n
				g=100.0_sp*abs(a(ip,iq))
				if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
					.and. (abs(d(iq))+g == abs(d(iq)))) then
					a(ip,iq)=0.0
				else if (abs(a(ip,iq)) > tresh) then
					h=d(iq)-d(ip)
					if (abs(h)+g == abs(h)) then
						t=a(ip,iq)/h
					else
						theta=0.5_sp*h/a(ip,iq)
						t=1.0_sp/(abs(theta)+sqrt(1.0_sp+theta**2))
						if (theta < 0.0) t=-t
					end if
					c=1.0_sp/sqrt(1+t**2)
					s=t*c
					tau=s/(1.0_sp+c)
					h=t*a(ip,iq)
					z(ip)=z(ip)-h
					z(iq)=z(iq)+h
					d(ip)=d(ip)-h
					d(iq)=d(iq)+h
					a(ip,iq)=0.0
					call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
					call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
					call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
					call jrotate(v(:,ip),v(:,iq))
					nrot=nrot+1
				end if
			end do
		end do
		b(:)=b(:)+z(:)
		d(:)=b(:)
		z(:)=0.0
	end do
	call nrerror('too many iterations in jacobi')
	CONTAINS
!BL
	SUBROUTINE jrotate(a1,a2)
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: a1,a2
	REAL(SP), DIMENSION(size(a1)) :: wk1
	wk1(:)=a1(:)
	a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
	a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
	END SUBROUTINE jrotate
	END SUBROUTINE jacobi
