MODULE hypgeo_info
	USE nrtype
	COMPLEX(SPC) :: hypgeo_aa,hypgeo_bb,hypgeo_cc,hypgeo_dz,hypgeo_z0
END MODULE hypgeo_info


	FUNCTION hypgeo(a,b,c,z)
	USE nrtype
	USE hypgeo_info
	USE nr, ONLY : bsstep,hypdrv,hypser,odeint
	IMPLICIT NONE
	COMPLEX(SPC), INTENT(IN) :: a,b,c,z
	COMPLEX(SPC) :: hypgeo
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	COMPLEX(SPC), DIMENSION(2) :: y
	REAL(SP), DIMENSION(4) :: ry
	if (real(z)**2+aimag(z)**2 <= 0.25) then
		call hypser(a,b,c,z,hypgeo,y(2))
		RETURN
	else if (real(z) < 0.0) then
		hypgeo_z0=cmplx(-0.5_sp,0.0_sp,kind=spc)
	else if (real(z) <= 1.0) then
		hypgeo_z0=cmplx(0.5_sp,0.0_sp,kind=spc)
	else
		hypgeo_z0=cmplx(0.0_sp,sign(0.5_sp,aimag(z)),kind=spc)
	end if
	hypgeo_aa=a
	hypgeo_bb=b
	hypgeo_cc=c
	hypgeo_dz=z-hypgeo_z0
	call hypser(hypgeo_aa,hypgeo_bb,hypgeo_cc,hypgeo_z0,y(1),y(2))
	ry(1:4:2)=real(y)
	ry(2:4:2)=aimag(y)
	call odeint(ry,0.0_sp,1.0_sp,EPS,0.1_sp,0.0001_sp,hypdrv,bsstep)
	y=cmplx(ry(1:4:2),ry(2:4:2),kind=spc)
	hypgeo=y(1)
	END FUNCTION hypgeo
