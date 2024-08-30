	MODULE xlinbcg_data
	USE nrtype
	TYPE(sprs2_dp) :: sa
	END MODULE xlinbcg_data

	PROGRAM xlinbcg
!	driver for routine linbcg
	USE nrtype
	USE nr
	USE xlinbcg_data
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=20,ITOL=1,ITMAX=75
	REAL(DP), PARAMETER :: TOL=1.e-9_dp
	INTEGER(I4B) :: i,iter
	REAL(DP) :: err
	REAL(DP), DIMENSION(NP) :: b,bcmp,x
	REAL(DP), DIMENSION(NP,NP) :: afull
	afull=0.0
	do i=1,NP
		if (i /= 1) afull(i,i-1)=-2.0
		afull(i,i)=1.0
		if (i /= NP) afull(i,i+1)=2.0
	end do
	call sprsin(afull,1.e-6_dp,sa)
	x=0.0
	b=1.0
	b(1)=3.0
	b(NP)=-1.0
	call linbcg(b,x,ITOL,TOL,ITMAX,iter,err)
	write(*,'(/1x,a,e15.6)') 'Estimated error:',err
	write(*,'(/1x,a,i6)') 'Iterations needed:',iter
	write(*,'(/1x,a)') 'Solution vector:'
	write(*,'(1x,5f12.6)') x
	call sprsax(sa,x,bcmp)
!	this is a double precision version of sprsax
	write(*,'(/1x,a)') 'press RETURN to continue...'
	read(*,*)
	write(*,'(1x,a/t8,a,t22,a)') 'Test of solution vector:','a*x','b'
	do i=1,NP
		write(*,'(1x,2f12.6)') bcmp(i),b(i)
	end do
	END PROGRAM xlinbcg