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
    public  :: UpsDischargesLinks
    public  :: ComputeUpscalingVelocity
    
    private :: SearchFace
    
    !begin-----------------------------------------------------------------------
    contains
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches discharge faces of father cell, and provides it to the son domain
    !>@param[in] Connection, SonWaterPoints, FatherWaterPoints, ICell, JCell, n_U, n_V, IZ, JZ
    subroutine SearchDischargeFace(Connection, SonWaterPoints, FatherWaterPoints, ICell, JCell, IZ, JZ, n_U, n_V)
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer, intent(IN)         :: Connection !Connection beteen father-Son Z cells
        integer, dimension(:, :), pointer, intent(IN)         :: SonWaterPoints, FatherWaterPoints, IZ, JZ
        integer, intent(IN)                                   :: ICell, JCell
        integer                                               :: n_U, n_V !Number of son cells in U and V direction
        !Local-----------------------------------------------------------------
        integer                                               :: di, dj, MaxSize, IFather, JFather
        !-------------------------------------------------------------------------
         
        MaxSize = size(Connection, 1)            
            
        !If father cell to the north is land, check for son cells(compare with adjacent southern cell)
        if (FatherWaterPoints(Icell + 1, JCell) == 0) then
            IFather = Icell + 1
            JFather = JCell
            di      = -1 !only going to search southwards
            dj      = 0
            call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, IZ, n = n_V)
        endif
        if (FatherWaterPoints(Icell - 1, JCell) == 0) then
            IFather = Icell - 1
            JFather = JCell
            di      = 1 !only going to search northwards
            dj      = 0
            call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, IZ, n = n_V)
        endif 
        if (FatherWaterPoints(Icell, JCell + 1) == 0) then
            IFather = Icell - 1
            JFather = JCell
            di      = 0 
            dj      = -1 !only going to search westward
            call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, JZ, n = n_U)
        endif
        if (FatherWaterPoints(Icell, JCell - 1) == 0) then
            IFather = Icell - 1
            JFather = JCell
            di      = 0 !only going to search eastward
            dj      = 1
            call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, JZ, n = n_U)
        endif        
    
    end subroutine SearchDischargeFace

    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches discharge faces of father cell, and its son cells responsible for an upscaling discharge.
    !> Will also update the connection matrix of the upscaling discharges (one matrix for all discharges)
    !>@param[in] Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, link, tracer, n, Cells
    subroutine SearchFace(Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, link, tracer, n, Cells)
        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer, intent(IN)    :: Connection
        integer, dimension(:, :), pointer, intent(IN)    :: SonWaterPoints, link
        integer, intent(IN)                              :: IFather, JFather, di, dj, MaxSize
        integer, optional                                :: n, tracer
        Integer, dimension(:, :), optional               :: Cells
        
        !Local-----------------------------------------------------------------
        integer                                          :: Aux, Aux2, StartIndex, i, ISon, JSon, ISonAdjacent, &
                                                            JSonAdjacent, IJFather
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
        if (present(Cells)) then
            do while (Aux2 == Aux)
                i = StartIndex
                ISon         = Connection(i, 3)
                JSon         = Connection(i, 4)
                ISonAdjacent = Connection(i, 3) + di
                JSonAdjacent = Connection(i, 4) + dj
                    
                if (SonWaterPoints(ISon, JSon) == 1)then
                    !Check if adjacent cell of son domain is inside Father dicharge cell
                    !this will need to change if a radious of search is used, as a son cell will &
                    !belong to more than 1 father cell
                    if (link(ISonAdjacent, JSonAdjacent) == (IJFather - 1))then
                        if (SonWaterPoints(ISonAdjacent, JSonAdjacent) == 1)then
                            tracer = tracer + 1
                            Cells(tracer, 1) = Connection(i, 1)
                            Cells(tracer, 2) = Connection(i, 2)
                            Cells(tracer, 3) = Connection(i, 3)
                            Cells(tracer, 4) = Connection(i, 4)
                        endif
                    endif
                endif
                i = 1 + 1
                Aux2 = Connection(i, 1) * Connection(i, 2)            
            enddo
        elseif (present(n)) then
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
            
        endif
        
    end subroutine SearchFace
    
    !-------------------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches discharge faces of father cell, and saves the upscaling discharge cell links between father-son cells
    !>@param[in] n, Cells, Connection, SonWaterPoints, FatherWaterPoints, ICell, JCell, IZ, JZ, I   
    subroutine UpsDischargesLinks (Connection, SonWaterPoints, FatherWaterPoints, ICell, JCell, IZ, JZ, VelID, tracer, Cells)
        !Arguments------------------------------------------------------------------------------------------------
        integer, dimension(:, :), pointer, intent(IN)         :: Connection !Connection beteen father-Son Z cells
        integer, dimension(:, :), pointer, intent(IN)         :: SonWaterPoints, FatherWaterPoints, IZ, JZ
        integer, intent(IN)                                   :: ICell, JCell, VelID, tracer
        integer, dimension(:, :)                              :: Cells
        !Local----------------------------------------------------------------------------------------------------
        integer                                               :: di, dj, MaxSize, IFather, JFather
        !---------------------------------------------------------------------------------------------------------
    
        MaxSize = size(Connection, 1)
        !If father cell to the north is land, check for son cells(compare with adjacent southern cell)
        
        if (VelID == VelocityU_) then
            if (FatherWaterPoints(Icell, JCell + 1) == 0) then
                IFather = Icell - 1
                JFather = JCell
                di      = 0 
                dj      = -1 !only going to search westward
                call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, JZ, Cells = Cells, &
                                 tracer = tracer)
            endif
            if (FatherWaterPoints(Icell, JCell - 1) == 0) then
                IFather = Icell - 1
                JFather = JCell
                di      = 0 !only going to search eastward
                dj      = 1
                call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, JZ, Cells = Cells, &
                                 tracer = tracer)
            endif
        else
            if (FatherWaterPoints(Icell + 1, JCell) == 0) then
                IFather = Icell + 1
                JFather = JCell
                di      = -1 !only going to search southwards
                dj      = 0
                call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, IZ, Cells = Cells, &
                                 tracer = tracer)
            endif
            if (FatherWaterPoints(Icell - 1, JCell) == 0) then
                IFather = Icell - 1
                JFather = JCell
                di      = 1 !only going to search northwards
                dj      = 0
                call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, IZ, Cells = Cells, &
                                 tracer = tracer)
            endif
        endif
        
    end subroutine UpsDischargesLinks
    !-------------------------------------------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches discharge faces of father cell, and saves the upscaling discharge cell links between father-son cells
    !>@param[in] DischargeVel, SonVel, Cells, Area, VelocityID     
    subroutine ComputeUpscalingVelocity(DischargeVel, SonVel, DLink, AreaSon, AreaFather, SonComputeFaces, KFloor, &
		                                KUBSon, KUBFather, VelocityID)
        !Arguments----------------------------------------------------------------------------------
        real, dimension(:, :, :),          intent(INOUT) :: DischargeVel
        real, dimension(:, :, :), pointer, intent(IN)    :: AreaSon, AreaFather, SonVel
        integer, dimension(:, :), pointer, intent(IN)    :: SonComputeFaces, KFloor
        integer, intent(IN)                              :: VelocityID
        integer, dimension(:, :), intent(IN)             :: DLink
        !Local--------------------------------------------------------------------------------------
        integer                                          :: i, j, Maxlines, line, KUBSon, KUBFather, k, k2, i, i2
        !-------------------------------------------------------------------------------------------
        Maxlines = size(CellsToUse, 1)

        if (VelocityID = VelocityU_) then
			
            do line = 1, Maxlines
                i = DLink(line, 1)
                j = DLink(line, 2)
				i2 = DLink(line, 3)
				j2 = DLink(line, 4)
				KBottom = KFloor(i, j)
				do k = Kbottom, KUBFather
					k2 = k - (KUBFather - KUBSon) 
					! I am assuming the vertical discretization will be the same even if the son has less layers.
					AcumulatedVel(i, j, k) = AcumulatedVel(i, j, k) + &
					                        SonVel(i2, j2, k2) * AreaSon(i2, j2, k2) / AreaFather(i, j, k)
				enddo
                  
			enddo
			
			
			
			
			
        else
        
        endif
    
    
    end subroutine ComputeUpscalingVelocity
    !
    !
    end module ModuleUpscallingDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon.
!----------------------------------------------------------------------------------------------------------