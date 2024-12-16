	SUBROUTINE hypdrv(s,ry,rdyds)
	USE nrtype
	USE hypgeo_info
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: s
	REAL(SP), DIMENSION(:), INTENT(IN) :: ry
	REAL(SP), DIMENSION(:), INTENT(OUT) :: rdyds
	COMPLEX(SPC), DIMENSION(2) :: y,dyds
	COMPLEX(SPC) :: z
	y=cmplx(ry(1:4:2),ry(2:4:2),kind=spc)
	z=hypgeo_z0+s*hypgeo_dz
	dyds(1)=y(2)*hypgeo_dz
	dyds(2)=((hypgeo_aa*hypgeo_bb)*y(1)-(hypgeo_cc-&
		((hypgeo_aa+hypgeo_bb)+1.0_sp)*z)*y(2))*hypgeo_dz/(z*(1.0_sp-z))
	rdyds(1:4:2)=real(dyds)
	rdyds(2:4:2)=aimag(dyds)
	END SUBROUTINE hypdrv
