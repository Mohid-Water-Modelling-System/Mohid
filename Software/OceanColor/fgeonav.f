	SUBROUTINE FGEONAV(pos,rm,coef,sun,nsta,ninc,npix,latitude,
     &  longitude)
c
      REAL*4, dimension(3)   :: pos   
	REAL*4, dimension(6)   :: coef  
	REAL*4, dimension(3,3) :: rm    
	REAL*4, dimension(3)   :: sun
      
      INTEGER                :: nsta  [REFERENCE]
      INTEGER                :: ninc  [REFERENCE] 
      INTEGER                :: npix  [REFERENCE] 
      
      REAL*4, dimension(npix) :: latitude
      REAL*4, dimension(npix) :: longitude 
      
      REAL*4 b
      INTEGER i 
c
c      
c	
      real geovec(3),no(3),up(3),ea(3),rmtq(3)
      real*8 pi,radeg,re,rem,f,omf2,omegae
      real*8 sinc,elev,sinl,cosl,h
	logical first/.true./
	common /navicom/sina(1285),cosa(1285),sinl,cosl
      common /gconst/pi,radeg,re,rem,f,omf2,omegae
	data sinc/0.0015911d0/
      data ea/0.0,0.0,0.0/
c
c
c  If first call, store array of sines and cosines of scan angles for 
c  future use
   	if (first) then
	  first = .false.
          call cdata
c
c  Compute elevation (out-of-plane) angle
	  elev = sinc*1.2
	  sinl = sin(elev)
	  cosl = cos(elev)
	  do i=1,1285
	    sina(i) = sin((i-643)*sinc)*cosl
	    cosa(i) = cos((i-643)*sinc)*cosl
	  end do
	end if
c
c  Compute correction factor for out-of-plane angle
	h = (rm(2,1)*pos(1)+rm(2,2)*pos(2)+rm(2,3)*pos(3)/omf2)*2.d0
c
c  Compute sensor-to-surface vectors for all scan angles
	do i=1,npix
	  in = ninc*(i-1) + nsta
	  a = coef(1)*cosa(in)*cosa(in)+coef(2)*cosa(in)*sina(in)
     *		+coef(3)*sina(in)*sina(in)
	  b = coef(4)*cosa(in)+coef(5)*sina(in)
	  c = coef(6)
	  r = b*b-4.d0*c*a  !begin solve quadratic equation
c
c  Check for scan past edge of Earth
	  if (r.lt.0.) then
	    latitude(i) = 999.
	    longitude(i) = 999.
	  else
c  Solve for magnitude of sensor-to-pixel vector and compute components
	    q = (-b-sqrt(r))/(2.d0*a)

c  Add out-of-plane correction 
	    q = q*(1.d0 + sinl*h/sqrt(r))
	    
	    Qx = q*cosa(in)
	    Qy = q*sinl
	    Qz = q*sina(in)
c
c  Transform vector from sensor to geocentric frame
	    do j=1,3
	      rmtq(j) = Qx*rm(1,j) + Qy*rm(2,j) + Qz*rm(3,j)
	      geovec(j) = rmtq(j) + pos(j)
	    end do
c
c    Compute geodetic latitude and longitude
          tmp = sqrt(geovec(1)*geovec(1)+geovec(2)*geovec(2))*omf2
	    latitude(i) = radeg*atan2(geovec(3),tmp)
	    longitude(i) = radeg*atan2(geovec(2),geovec(1))
c         
         endif
      enddo
c	
      END