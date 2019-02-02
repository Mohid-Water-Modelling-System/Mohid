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

    !public  :: BuildDischargesMatrix
    public  :: SearchDischargeFace
    
    
    private :: SearchFace
    
    !begin-----------------------------------------------------------------------
    contains
    
    subroutine SearchDischargeFace(ConnectionMatrix, SonWaterPoints2D, FatherWaterPoints2D, SonSize2D, &
                                   ICell, JCell, CellsToAllocate)
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :, :, :), pointer   :: ConnectionMatrix
        integer, dimension(:, :), pointer         :: SonWaterPoints2D, FatherWaterPoints2D
        type (T_Size2D)                           :: SonSize2D
        integer                                   :: ICell, JCell
        integer, intent(OUT), optional            :: CellsToAllocate
        !Local-----------------------------------------------------------------
        integer                                   :: di, dj
        !-------------------------------------------------------------------------
        
        if (present(CellsToAllocate)) then
            
            Size = size(ConnectionMatrix)            
            
            !If father cell to the north is land, check for son cells(compare with adjacent southern cell)
            if (FatherWaterPoints2D(Icell + 1, JCell) == 0) then
                IFather = Icell + 1
                JFather = JCell
                di      = -1 !only going to search southwards
                dj      = 0
                call SearchFace (ConnectionMatrix, Size, IFather, JFather, di, dj, SonWaterPoints2D, &
                                 CellsToAllocate, link = Ilink)
            endif
            if (FatherWaterPoints2D(Icell - 1, JCell) == 0) then
                IFather = Icell - 1
                JFather = JCell
                di      = 1 !only going to search northwards
                dj      = 0
                call SearchFace (ConnectionMatrix, Size, IFather, JFather, di, dj, SonWaterPoints2D, &
                                 CellsToAllocate, link = Ilink)
            endif 
            if (FatherWaterPoints2D(Icell, JCell + 1) == 0) then
                IFather = Icell - 1
                JFather = JCell
                di      = 0 
                dj      = -1 !only going to search westward
                call SearchFace (ConnectionMatrix, Size, IFather, JFather, di, dj, SonWaterPoints2D, &
                                 CellsToAllocate, link = Jlink)
            endif
            if (FatherWaterPoints2D(Icell, JCell - 1) == 0) then
                IFather = Icell - 1
                JFather = JCell
                di      = 0 !only going to search eastward
                dj      = 1
                call SearchFace (ConnectionMatrix, Size, IFather, JFather, di, dj, SonWaterPoints2D, &
                                 CellsToAllocate, link = Jlink)
            endif
            

        else
            
            
            
        endif
        
        
        
        
        
        
        
        
        !----------------------------------------------------------------------
        !ConnectionsUB = SizeSon%IUB * SizeSon%JUB
        !do i = 1, ConnectionsUB
        !    iFather = Connections(i, 1)
        !    jFather = Connections(i, 2)
        !    
        !    k = i
        !    !find next father cell in Connections matrix
        !    do
        !        k = k + 1
        !        if ((Connections(k, 1) * Connections(k, 2)) .NE. aux)then
        !            exit
        !        endif
        !    enddo
        !    
        !    if (Connections(k, 1) == iFather)then
        !        
        !        if (FatherWaterPoints3D(iFather, jFather) .and. FatherLandPoints2D(iFather, jFather + 1))then
        !            aux = iFather * (jFather + 1)
        !            aux2 = aux
        !            do while (aux2 == aux)
        !                k = k + 1
        !                if (k == ConnectionsUB - 1)then
        !                    exit
        !                endif
        !                
        !                if (
        !                
        !                iSon = Connections(k, 3)
        !                jSon = Connections(k, 4)
        !                
        !                if (SonWaterPoints3D(iSon, jSon) == 1 .and. SonWaterPoints3D(iSon, jSon - 1) == 1) then
        !                    Flag1 = .true.
        !                endif
        !                
        !
        !                aux2 = Connections(k, 1) * Connections(k, 2)                
        !            enddo
        !
        !        endif
        !
        !    endif
        !    
        !    Waterpoint3D(ison, json) * FatherWaterPoints3D(iFather, jFather)
        !    Waterpoint3D(ison, json) * Waterpoint3D(ison, json)
        !    
        !enddo
        !
    
    end subroutine SearchDischargeFace


    subroutine SearchFace(Connection, Size, IFather, JFather, di, dj, SonWaterPoints, n, link)
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :, :, :), pointer   :: Connection
        integer, dimension(:, :), pointer         :: SonWaterPoints, link
        integer                                   :: IFather, JFather, n, di, dj, Size
        
        !Local-----------------------------------------------------------------
        integer                         :: Aux, Aux2, StartIndex, i, ISon, JSon, ISonAdjacent, JSonAdjacent, IJFather
        !----------------------------------------------------------------------
        
        !Find index of matrix where connections to cell (IFather, JCell) begin
        do i = 1, Size
            if (Connection(i, 1) == IFather)then
                if (Connection(i, 2) == JFather)then
                    StartIndex = i
                    Aux = IFather * JFather
                    exit
                endif
            endif
        enddo
        
        if (di /=0) IJFather = IFather
        if (dj /=0) IJFather = JFather
        
        !Check if northern face needs to be considered for the discharge velocity
        do while (Aux2 == Aux)
            i = StartIndex
            ISon         = Connection(i, 3)
            JSon         = Connection(i, 4)
            ISonAdjacent = Connection(i, 3) + di
            JSonAdjacent = Connection(i, 3) + dj
                    
            if (SonWaterPoints(ISon, JSon) == 1)then
                if (link(ISonAdjacent, JSonAdjacent) == IJFather)then
                    if (SonWaterPoints2D(ISonAdjacent, JSonAdjacent) == 1)then
                        n = n + 1 ! Found a discharge face
                        exit
                    endif

                endif
            endif
            i = 1 + 1
            Aux2 = Connection(i, 1) * Connection(i, 2)
                         
        enddo     
    
    end subroutine SearchFace
    



    end module ModuleUpscallingDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon.
!----------------------------------------------------------------------------------------------------------
