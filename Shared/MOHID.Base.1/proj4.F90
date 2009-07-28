! Copyright 2004, Magnus Hagdorn
! 
! This file is part of proj4.
! 
! proj4 is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
! 
! proj4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with proj4; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

module proj4
  include 'proj4.inc'
  integer, parameter :: PRJ90_NOERR = PRJF_NOERR

  type prj90_projection
     private
     integer :: prj
  end type prj90_projection

  interface prj90_fwd
     module procedure prj90_fwd_pt, prj90_fwd_array
  end interface

  interface prj90_inv
     module procedure prj90_inv_pt, prj90_inv_array
  end interface

contains
  function prj90_strerrno(prj_errno)
    implicit none
    character(len=80) :: prj90_strerrno
    integer, intent(in) :: prj_errno

    prj90_strerrno = prjf_strerrno(prj_errno)
  end function prj90_strerrno

  function prj90_init(prj,args)
    implicit none
    integer :: prj90_init
    type(prj90_projection) :: prj
    character(len=*),dimension(:) :: args

    prj90_init = prjf_init(prj%prj,size(args),len(args),args)
  end function prj90_init

  function prj90_free(prj)
    implicit none
    integer :: prj90_free
    type(prj90_projection) :: prj

    prj90_free =  prjf_free(prj%prj)
  end function prj90_free

  function prj90_fwd_pt(prj,lam,phi,x,y)
    implicit none
    integer :: prj90_fwd_pt
    type(prj90_projection) :: prj
    real(kind=kind(1.0d0)), intent(in) :: lam, phi
    real(kind=kind(1.0d0)), intent(out) :: x,y

    prj90_fwd_pt = prjf_fwd(prj%prj,lam,phi,x,y)
  end function prj90_fwd_pt
    
  function prj90_fwd_array(prj,lam,phi,x,y)
    implicit none
    integer :: prj90_fwd_array
    type(prj90_projection) :: prj
    real(kind=kind(1.0d0)), dimension(:), intent(in) :: lam, phi
    real(kind=kind(1.0d0)), dimension(:), intent(out) :: x,y

    integer, dimension(size(lam)) :: res
    integer i

    do i=1,size(lam)
       res(i) = prj90_fwd_pt(prj,lam(i),phi(i),x(i),y(i))
    end do

    if (any(res.ne.PRJ90_NOERR)) then
       do i=1,size(lam)
          if (res(i).ne.PRJ90_NOERR) then
             prj90_fwd_array = res(i)
          end if
       end do
    else
       prj90_fwd_array = PRJ90_NOERR
    end if
  end function prj90_fwd_array

  function prj90_inv_pt(prj,x,y,lam,phi)
    implicit none
    integer :: prj90_inv_pt
    type(prj90_projection) :: prj
    real(kind=kind(1.0d0)), intent(in) :: x,y
    real(kind=kind(1.0d0)), intent(out) :: lam, phi

    prj90_inv_pt = prjf_inv(prj%prj,x,y,lam,phi)
  end function prj90_inv_pt  

  
  function prj90_inv_array(prj,x,y,lam,phi)
    implicit none
    integer :: prj90_inv_array
    type(prj90_projection) :: prj
    real(kind=kind(1.0d0)), dimension(:), intent(in) :: x,y
    real(kind=kind(1.0d0)), dimension(:), intent(out) :: lam, phi

    integer, dimension(size(x)) :: res
    integer i

    do i=1,size(x)    
       res(i) = prj90_inv_pt(prj,x(i),y(i),lam(i),phi(i))
    end do

    
    if (any(res.ne.PRJ90_NOERR)) then
       do i=1,size(x)
          if (res(i).ne.PRJ90_NOERR) then
             prj90_inv_array = res(i)
          end if
       end do
    else
       prj90_inv_array = PRJ90_NOERR
    end if
  end function prj90_inv_array

end module proj4
