!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : TwoWayDischarges
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jan 2019
! REVISION      : Joao Sobrinho - v1.0
! DESCRIPTION   : Module which contains routines to build twoway discharges
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

Module ModuleTwoWayDischarges

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData,        only : GetData

    implicit none

    private

    !types---------------------------------------------------------------------


    !Subroutines-----------------------------------------------------------------

    public  :: BuildDischargesMatrix
    
    !begin-----------------------------------------------------------------------
    contains

    subroutine BuildDischargesMatrix(Connections, WaterPoints3D, FatherWaterPoints3D, TwoWayMomentumFlow)
    
    end subroutine BuildDischargesMatrix
    



    end module ModuleTwoWayDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon.
!----------------------------------------------------------------------------------------------------------
