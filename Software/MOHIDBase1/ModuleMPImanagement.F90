!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : ModuleMPImanagement
! PROJECT       : Example of the actor model using Fortran and MPI
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Nov 2014
! REVISION      : Rciardo Miranda - v1.0
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

module ModuleMPImanagement

    implicit none

    private

    integer, parameter :: msgSendTaskTag      = 11111
    integer, parameter :: msgPoisonPillTag    = 11122
    integer, parameter :: msgIdleWorkerTag    = 33333
    integer, parameter :: msgTaskCompletedTag = 44444
    integer, parameter :: msgEndGameTag       = 34232
    integer, parameter :: msgMsgRankTag       = 34555

   !Functions-----------------------------------------------------------------

    public  :: getMsgSendTaskTag
    public  :: getMsgPoisonPillTag
    public  :: getMsgIdleWorkerTag
    public  :: getMsgTaskCompletedTag
    public  :: getMsgRankTag

    contains

    integer function getMsgSendTaskTag()
        getMsgSendTaskTag = msgSendTaskTag
    end function getMsgSendTaskTag

    integer function getMsgPoisonPillTag()
        getMsgPoisonPillTag = msgPoisonPillTag
    end function getMsgPoisonPillTag

    integer function getMsgIdleWorkerTag()
        getMsgIdleWorkerTag = msgIdleWorkerTag
    end function getMsgIdleWorkerTag

    integer function getMsgTaskCompletedTag()
        getMsgTaskCompletedTag = msgTaskCompletedTag
    end function getMsgTaskCompletedTag

    integer function getMsgRankTag()
        getMsgRankTag = msgMsgRankTag
    end function getMsgRankTag

end module ModuleMPImanagement
