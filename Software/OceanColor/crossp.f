        subroutine crossp(v1,v2,v3)
c
c  crossp(v1,v2,v3)
c
c  Purpose: Computes vector cross product v3 = v1 x v2
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  v1(3)        R*4      I      Input Vector
c  v2(3)        R*4      I      Input Vector
c  v3(3)        R*4      O      Output Vector
c
c  By: Fred Patt, GSC
c
c  Notes:  
c
c  Modification History:
c
      implicit none
c
	real v1(3),v2(3),v3(3)
c
c
	v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
	v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
	v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
c
	return
	end
