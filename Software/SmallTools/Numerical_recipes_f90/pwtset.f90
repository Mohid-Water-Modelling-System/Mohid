MODULE pwtcom
	USE nrtype
	INTEGER(I4B), SAVE :: ncof=0,ioff,joff
	REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: cc,cr
END MODULE pwtcom

	SUBROUTINE pwtset(n)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE pwtcom
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: sig
	REAL(SP), PARAMETER :: &
		c4(4)=(/&
		0.4829629131445341_sp, 0.8365163037378079_sp, &
		0.2241438680420134_sp,-0.1294095225512604_sp /), &
		c12(12)=(/&
		0.111540743350_sp, 0.494623890398_sp, 0.751133908021_sp, &
		0.315250351709_sp,-0.226264693965_sp,-0.129766867567_sp, &
		0.097501605587_sp, 0.027522865530_sp,-0.031582039318_sp, &
		0.000553842201_sp, 0.004777257511_sp,-0.001077301085_sp /), &
		c20(20)=(/&
		0.026670057901_sp, 0.188176800078_sp, 0.527201188932_sp, &
		0.688459039454_sp, 0.281172343661_sp,-0.249846424327_sp, &
		-0.195946274377_sp, 0.127369340336_sp, 0.093057364604_sp, &
		-0.071394147166_sp,-0.029457536822_sp, 0.033212674059_sp, &
		0.003606553567_sp,-0.010733175483_sp, 0.001395351747_sp, &
		0.001992405295_sp,-0.000685856695_sp,-0.000116466855_sp, &
		0.000093588670_sp,-0.000013264203_sp /)
	if (allocated(cc)) deallocate(cc)
	if (allocated(cr)) deallocate(cr)
	allocate(cc(n),cr(n))
	ncof=n
	ioff=-n/2
	joff=-n/2
	sig=-1.0
	select case(n)
		case(4)
			cc=c4
		case(12)
			cc=c12
		case(20)
			cc=c20
		case default
			call nrerror('unimplemented value n in pwtset')
	end select
	cr(n:1:-1) = cc
	cr(n:1:-2) = -cr(n:1:-2)
	END SUBROUTINE pwtset
