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

!     global variables
      integer prjf_noerr
      parameter (prjf_noerr=0)

!     routines
      character*80 prjf_strerrno
      external prjf_strerrno

      integer prjf_init
      external prjf_init

      integer prjf_free
      external prjf_free

      integer prjf_fwd
      external prjf_fwd

      integer prjf_inv
      external prjf_inv
