# ModuleMPImanagement

!------------------------------------------------------------------------------ ! IST/MARETEC, Water Modelling Group, Mohid modelling system !------------------------------------------------------------------------------

!

! TITLE  : moduleMPImanagement

! PROJECT  : Example of the actor model using Fortran and MPI

! URL  : [http://www.mohid.com](http://www.mohid.com)

! AFFILIATION  : IST/MARETEC, Marine Modelling Group

! DATE  : Nov 2014

! REVISION  : Rciardo Miranda - v1.0

!

!------------------------------------------------------------------------------

!

!This program is free software; you can redistribute it and/or

!modify it under the terms of the GNU General Public License

!version 2, as published by the Free Software Foundation.

!

!This program is distributed in the hope that it will be useful,

!but WITHOUT ANY WARRANTY; without even the implied warranty of

!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the

!GNU General Public License for more details.

!

!You should have received a copy of the GNU General Public License

!along with this program; if not, write to the Free Software

!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

!

!------------------------------------------------------------------------------

module moduleMPImanagement

```
   implicit none
```

```
   private
```

```
   integer, parameter :: msgPlayBall1Tag = 11111
   integer, parameter :: msgPlayBall2Tag = 11122
   integer, parameter :: msgPlayBall3Tag = 33333
   integer, parameter :: msgPlayBall4Tag = 44444
   integer, parameter :: msgMyMPITag     = 22222
   integer, parameter :: msgEndGameTag   = 34232
   integer, parameter :: msgMsgIDTag     = 34555
```

```
   !Subroutines---------------------------------------------------------------
```

```
   public  :: getMsgPlayBall1Tag
   public  :: getMsgPlayBall2Tag
   public  :: getMsgPlayBall3Tag
   public  :: getMsgPlayBall4Tag
   public  :: getMsgMyMPITag
   public  :: getMsgEndGameTag
   public  :: getMsgIDTag
```

```
   contains
```

```
   integer function getMsgPlayBall1Tag()
       getMsgPlayBall1Tag = msgPlayBall1Tag
   end function getMsgPlayBall1Tag
```

```
   integer function getMsgPlayBall2Tag()
       getMsgPlayBall2Tag = msgPlayBall2Tag
   end function getMsgPlayBall2Tag
```

```
   integer function getMsgPlayBall3Tag()
       getMsgPlayBall3Tag = msgPlayBall3Tag
   end function getMsgPlayBall3Tag
```

```
   integer function getMsgPlayBall4Tag()
       getMsgPlayBall4Tag = msgPlayBall4Tag
   end function getMsgPlayBall4Tag
```

```
   integer function getMsgMyMPITag()
       getMsgMyMPITag = msgMyMPITag
   end function getMsgMyMPITag
```

```
   integer function getMsgEndGameTag()
       getMsgEndGameTag = msgEndGameTag
   end function getMsgEndGameTag
```

```
   integer function getMsgIDTag()
       getMsgIDTag = msgMsgIDTag
   end function getMsgIDTag
```

end module moduleMPImanagement
