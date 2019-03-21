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

    public  :: SearchDischargeFace
    
    
    private :: SearchFace
    
    !begin-----------------------------------------------------------------------
    contains
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches discharge faces of father cell, and provides it to the son domain
    !>@param[in] ConnectionMatrix, SonWaterPoints, FatherWaterPoints, ICell, JCell, CellsU, CellsV
    subroutine SearchDischargeFace(ConnectionMatrix, SonWaterPoints, FatherWaterPoints, ICell, JCell, n_U, n_V, IZ, JZ)
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer, intent(IN)         :: ConnectionMatrix
        integer, dimension(:, :), pointer, intent(IN)         :: SonWaterPoints, FatherWaterPoints, IZ, JZ
        integer, intent(IN)                                   :: ICell, JCell
        integer, intent(INOUT)                                :: n_U, n_V !Number of son cells in U and V direction
        !Local-----------------------------------------------------------------
        integer                                               :: di, dj, MaxSize, IFather, JFather
        !-------------------------------------------------------------------------
         
        MaxSize = size(ConnectionMatrix, 1)            
            
        !If father cell to the north is land, check for son cells(compare with adjacent southern cell)
        if (FatherWaterPoints(Icell + 1, JCell) == 0) then
            IFather = Icell + 1
            JFather = JCell
            di      = -1 !only going to search southwards
            dj      = 0
            call SearchFace (ConnectionMatrix, MaxSize, IFather, JFather, di, dj, SonWaterPoints, n_V, link = IZ)
        endif
        if (FatherWaterPoints(Icell - 1, JCell) == 0) then
            IFather = Icell - 1
            JFather = JCell
            di      = 1 !only going to search northwards
            dj      = 0
            call SearchFace (ConnectionMatrix, MaxSize, IFather, JFather, di, dj, SonWaterPoints, n_V, link = IZ)
        endif 
        if (FatherWaterPoints(Icell, JCell + 1) == 0) then
            IFather = Icell - 1
            JFather = JCell
            di      = 0 
            dj      = -1 !only going to search westward
            call SearchFace (ConnectionMatrix, MaxSize, IFather, JFather, di, dj, SonWaterPoints, n_U, link = JZ)
        endif
        if (FatherWaterPoints(Icell, JCell - 1) == 0) then
            IFather = Icell - 1
            JFather = JCell
            di      = 0 !only going to search eastward
            dj      = 1
            call SearchFace (ConnectionMatrix, MaxSize, IFather, JFather, di, dj, SonWaterPoints, n_U, link = JZ)
        endif        
    
    end subroutine SearchDischargeFace


    subroutine SearchFace(Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, n, link)
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer         :: Connection
        integer, dimension(:, :), pointer         :: SonWaterPoints, link
        integer                                   :: IFather, JFather, di, dj, MaxSize
        integer, intent(INOUT)                    :: n
        
        !Local-----------------------------------------------------------------
        integer                         :: Aux, Aux2, StartIndex, i, ISon, JSon, ISonAdjacent, JSonAdjacent, IJFather
        !----------------------------------------------------------------------
        
        !Find index of matrix where connections to cell (IFather, JCell) begin
        !columns in connection(:) : 1 - IFather; 2 - JFather; 3 - ISon; 4 - JSon
        do i = 1, MaxSize
            if (Connection(i, 1) == IFather)then
                if (Connection(i, 2) == JFather)then
                    StartIndex = i
                    Aux = IFather * JFather
                    exit
                endif
            endif
        enddo
        
        if (di /=0) IJFather = IFather ! means we are searching the north/South direction
        if (dj /=0) IJFather = JFather ! means we are searching the west/east direction
        
        !Check if current face needs to be considered for the discharge velocity
        do while (Aux2 == Aux)
            i = StartIndex
            ISon         = Connection(i, 3)
            JSon         = Connection(i, 4)
            ISonAdjacent = Connection(i, 3) + di
            JSonAdjacent = Connection(i, 4) + dj
                    
            if (SonWaterPoints(ISon, JSon) == 1)then
                !Check if adjacent cell of son domain is inside Father dicharge cell
                if (link(ISonAdjacent, JSonAdjacent) == (IJFather - 1))then
                    if (SonWaterPoints(ISonAdjacent, JSonAdjacent) == 1)then
                        n = n + 1 ! Found a discharge face
                    endif
                endif
            endif
            i = 1 + 1
            Aux2 = Connection(i, 1) * Connection(i, 2)            
        enddo     
    
    end subroutine SearchFace
    !
    !
    end module ModuleUpscallingDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon.
!----------------------------------------------------------------------------------------------------------
