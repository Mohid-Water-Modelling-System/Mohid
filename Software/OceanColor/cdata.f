      subroutine cdata

c  This subroutines initializes several constants for general use in 
c  navigation routines.  The constants are passed in the common /gconst/.
c
c  Calling Arguments:  None
c
c  Variables in common /gconst/
c
c  	Name		Type 	Description
c
c  	pi		R*8	Pi
c  	radeg		R*8	Radians to degrees conversion factor
c  	re		R*8	Earth equatorial radius (km)
c  	rem		R*8	Earth mean radius (km)
c  	f		R*8	Earth flattening factor
c  	omf2		R*8 	(1-f)**2
c	omegae		R*8 	Earth rotation rate (radians/second)
c
c
c	Program written by:	Frederick S. Patt
c				General Sciences Corporation
c				May 10, 1993
c
c	Modification History:

	real*8 pi,radeg,re,rem,f,omf2,omegae
	common /gconst/pi,radeg,re,rem,f,omf2,omegae

	pi = dacos(-1.0D0)
	re = 6378.137d0 
	rem = 6371.d0
	radeg = 180.d0/pi
	f = 1.d0/298.257d0
	omf2 = (1.d0-f)**2
        omegae = 7.29211585494d-5

	return
	end
