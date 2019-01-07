!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Discharges
! PROJECT       : Mohid Base 2
! MODULE        : UpscallingDischarges
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jan 2019
! REVISION      : Joao Sobrinho - v1.0
! DESCRIPTION   : Module which contains routines to build discharges set by a child domain
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

Module ModuleUpscallingDischarges

    use ModuleGlobalData

    implicit none

    private

    !types---------------------------------------------------------------------


    !Subroutines-----------------------------------------------------------------

    public  :: BuildDischargesMatrix
    public  :: SearchForDischarges
    
    !begin-----------------------------------------------------------------------
    contains
    
    subroutine SearchForDischarges(Connections, SonWaterPoints3D, FatherWaterPoints3D, SonLandPoints2D, &
                                   FatherLandPoints2D, SizeSon, SizeFather, Present)
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer   :: Connections
        integer, dimension(:, :), pointer   :: SonWaterPoints3D, FatherWaterPoints3D, SonLandPoints2D, &
                                               FatherLandPoints2D
        type (T_Size2D)                     :: SizeSon, SizeFather
        logical                             :: Present
        integer                             :: i, iFather, jFather, iSon, jSon, k, ConnectionsUB
        !Local-----------------------------------------------------------------
        
        !----------------------------------------------------------------------
        ConnectionsUB = SizeSon%IUB * SizeSon%JUB
        do i = 1, ConnectionsUB
            iFather = Connections(i, 1)
            jFather = Connections(i, 2)
            
            k = i
            !find next father cell in Connections matrix
            do
                k = k + 1
                if ((Connections(k, 1) * Connections(k, 2)) .NE. aux)then
                    exit
                endif
            enddo
            
            if (Connections(k, 1) == iFather)then
                
                if (FatherWaterPoints3D(iFather, jFather) .and. FatherLandPoints2D(iFather, jFather + 1))then
                    aux = iFather * (jFather + 1)
                    aux2 = aux
                    do while (aux2 == aux)
                        k = k + 1
                        if (k == ConnectionsUB - 1)then
                            exit
                        endif
                        
                        if (
                        
                        iSon = Connections(k, 3)
                        jSon = Connections(k, 4)
                        
                        if (SonWaterPoints3D(iSon, jSon) == 1 .and. SonWaterPoints3D(iSon, jSon - 1) == 1) then
                            Flag1 = .true.
                        endif
                        

                        aux2 = Connections(k, 1) * Connections(k, 2)                
                    enddo

                endif

            endif
            
            Waterpoint3D(ison, json) * FatherWaterPoints3D(iFather, jFather)
            Waterpoint3D(ison, json) * Waterpoint3D(ison, json)
            
        enddo
        
    
    end subroutine SearchForDischarges


    subroutine BuildDischargesMatrix(Connections, WaterPoints3D, FatherWaterPoints3D, MomentumDischargesMatrix)
        !Arguments-------------------------------------------------------------
        integer                            :: FatherTwoWayID, TwoWayID
        integer, dimension(:), pointer     :: Connections
        integer, dimension(:,:), pointer   :: WaterPoints3D, FatherWaterPoints3D
        real, dimension(:,:,:), pointer    :: 
        !Local-----------------------------------------------------------------
        integer                            :: ready_, ILB, IUB, JLB, JUB, KLB, KUB, STAT_CALL
        !----------------------------------------------------------------------
        
        
    
    end subroutine BuildDischargesMatrix
    



    end module ModuleUpscallingDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon.
!----------------------------------------------------------------------------------------------------------
