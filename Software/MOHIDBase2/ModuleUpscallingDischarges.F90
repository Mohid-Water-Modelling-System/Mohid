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
    public  :: ComputeDischargeVolume
    public  :: UpdateDischargeConnections
    
    private :: SearchFace
    private :: Update_n_Z
    
    !begin-----------------------------------------------------------------------
    contains
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Searches discharge faces of father cell, and provides it to the son domain
    !>@param[in] Connection, SonWaterPoints, FatherWaterPoints, ICell, JCell, n_U, n_V, IZ, JZ
    subroutine SearchDischargeFace(Connection, SonWaterPoints, FatherWaterPoints, ICell, JCell, IZ, JZ, n_U, n_V, n_Z,&
                                   Kfloor, KUB)
        !Arguments------------------------------------------------------------- 
        integer, dimension(:, :), pointer, intent(IN)         :: Connection !Connection beteen father-Son Z cells
        integer, dimension(:, :), pointer, intent(IN)         :: SonWaterPoints, FatherWaterPoints, IZ, JZ, Kfloor
        integer, intent(IN)                                   :: ICell, JCell
        integer                                               :: n_U, n_V, n_Z, KUB !Number of son cells in U/V directions
        !Local-----------------------------------------------------------------
        integer                                               :: di, dj, MaxSize, IFather, JFather
        !-------------------------------------------------------------------------
         
        MaxSize = size(Connection, 1)
        
        call Update_n_Z(n_Z, Kfloor, KUB, ICell, JCell)
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
            IFather = Icell
            JFather = JCell + 1
            di      = 0 
            dj      = -1 !only going to search westward
            call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, JZ, n = n_U)
        endif
        if (FatherWaterPoints(Icell, JCell - 1) == 0) then
            IFather = Icell
            JFather = JCell - 1
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
        Aux2 = Aux
        i = StartIndex
        !Check if current face needs to be considered for the discharge velocity
        if (present(Cells)) then
            do while (Aux2 == Aux)
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
                i = i + 1
                Aux2 = Connection(i, 1) * Connection(i, 2)            
            enddo
        elseif (present(n)) then
            do while (Aux2 == Aux)
                ISon         = Connection(i, 3)
                JSon         = Connection(i, 4)
                ISonAdjacent = Connection(i, 3) + di
                JSonAdjacent = Connection(i, 4) + dj
                    
                if (SonWaterPoints(ISon, JSon) == 1)then
                    !Check if adjacent cell of son domain is inside Father dicharge cell
                    if (link(ISonAdjacent, JSonAdjacent) == (IJFather - 1))then
                        if (SonWaterPoints(ISonAdjacent, JSonAdjacent) == 1)then
                            n = n + 1 ! Found a discharge face   AQUI
                        endif
                    endif
                endif
                i = i + 1
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
                IFather = Icell
                JFather = JCell + 1
                di      = 0 
                dj      = -1 !only going to search westward
                call SearchFace (Connection, MaxSize, IFather, JFather, di, dj, SonWaterPoints, JZ, Cells = Cells, &
                                 tracer = tracer)
            endif
            if (FatherWaterPoints(Icell, JCell - 1) == 0) then
                IFather = Icell
                JFather = JCell - 1
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
    !> Computes area averaged velocity of son cells nested in a father cell
    !>@param[in] DischargeVel, SonVel, DLink, SonArea, FatherArea, SonComputeFaces, AcumulatedVel, AuxArea, &
    !> KFloor, KUBSon, KUBFather 
    subroutine ComputeUpscalingVelocity(DischargeVel, SonVel, DLink, SonArea, FatherArea, SonComputeFaces, &
                                        AcumulatedVel, AuxArea, KFloor, KUBSon, KUBFather)
        !Arguments----------------------------------------------------------------------------------
        real, dimension(:, :, :),          intent(INOUT)    :: DischargeVel
        real, dimension(:, :, :),          intent(INOUT)    :: AcumulatedVel, AuxArea 
        real, dimension(:, :, :),    pointer, intent(IN)    :: SonArea, FatherArea, SonVel
        integer, dimension(:, :),    pointer, intent(IN)    :: KFloor
        integer, dimension(:, :, :), pointer, intent(IN)    :: SonComputeFaces
        integer, dimension(:, :),             intent(IN)    :: DLink
        !Local--------------------------------------------------------------------------------------
        integer                                             :: i, j, k, Maxlines, line, KUBSon, KUBFather, &
                                                               i2, j2, k2, Kbottom
        !-------------------------------------------------------------------------------------------
        Maxlines = size(DLink, 1)
        !Not worth paralelizing... too litle work for each thread.
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
                                         SonVel(i2, j2, k2) * SonArea(i2, j2, k2) * SonComputeFaces(i2, j2, k2)
                    
                AuxArea(i, j, k) = AuxArea(i, j, k) + SonArea(i2, j2, k2)
            enddo
                  
        enddo
            
        do line = 1, Maxlines
            i = DLink(line, 1)
            j = DLink(line, 2)
            KBottom = KFloor(i, j)
            do k = Kbottom, KUBFather
                DischargeVel(i, j, k) = AcumulatedVel(i, j, k) / FatherArea(i, j, k)
            enddo
                  
        enddo            
    
    end subroutine ComputeUpscalingVelocity
    !--------------------------------------------------------------------------------------
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Computes volume to be added or removed due to upscaling discharge
    !>@param[in] FatherU_old, FatherU, FatherV_old, FatherV, Volume, KUB, KFloorU, KFloorV, &
    !>                                 AreaU, AreaV, CellsZ                                        
                                      
    subroutine ComputeDischargeVolume(FatherU_old, FatherU, FatherV_old, FatherV, AreaU, AreaV, &
        UpscaleFlow, DischargeConnection)
        !Arguments--------------------------------------------------------------------------
        real,    dimension(:, :, :), pointer, intent(IN)     :: FatherU, FatherV, AreaU, AreaV
        real,    dimension(:, :, :), allocatable, intent(IN) :: FatherU_old, FatherV_old
        real(8),    dimension(:)            , intent(INOUT)  :: UpscaleFlow
        integer, dimension(:, :)            , intent(IN)     :: DischargeConnection
        !Local-------------------------------------------------------------------------------
        integer                                              :: line, i, j, k, MaxSize
        real                                                 :: F_East, F_West, F_South, F_North
        !------------------------------------------------------------------------------------
        MaxSize = size(UpscaleFlow)
        do line = 1, MaxSize
            i = DischargeConnection(line, 1)
            j = DischargeConnection(line, 2)
            k = DischargeConnection(line, 3)
            
            F_West = (FatherU_old(i, j  , k) - FatherU(i, j  , k)) * AreaU(i, j  , k)
            F_East = -(FatherU_old(i, j+1, k) - FatherU(i, j+1, k)) * AreaU(i, j+1, k)
            
            F_South = (FatherV_old(i  , j, k) - FatherV(i  ,j , k)) * AreaV(i  , j, k)
            F_North = -(FatherV_old(i+1, j, k) - FatherV(i+1,j , k)) * AreaV(i+1, j, k)
                
            UpscaleFlow(line) = F_South + F_North + F_East + F_West
                
            UpscaleFlow(line) = F_South + F_North + F_East + F_West
        enddo
    end subroutine ComputeDischargeVolume
        
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Updates n_Z matrix
    !>@param[in] n_Z, Kfloor, KUB, I, J         
    subroutine Update_n_Z(n_Z, KFloorZ, KUB, I, J)
        !Arguments--------------------------------------------------------------------------
        integer, dimension(:, :), pointer, intent(IN)        :: KFloorZ
        integer                          , intent(IN)        :: KUB, I, J
        integer                          , intent(INOUT)     :: n_Z
        !Local-------------------------------------------------------------------------------
        integer                                              :: k, KBottom
        !------------------------------------------------------------------------------------
        KBottom = KFloorZ(I, J)
        do k = KBottom, KUB
            n_Z = n_Z + 1
        enddo
    end subroutine Update_n_Z
    !----------------------------------------------------------------------------------------
    
    !>@author Joao Sobrinho Maretec
    !>@Brief
    !> Updates discharge matrix table with the I, J and K  Parent cells
    !>@param[in] n_Z, Kfloor, KUB, I, J       
    subroutine UpdateDischargeConnections(CurrentZ, Matrix, KFloor, KUB, I, J)
        !Arguments--------------------------------------------------------------------------
        integer, dimension(:, :), allocatable , intent(INOUT)     :: Matrix
        integer,                           intent(IN)             :: I, J, KUB
        integer, dimension(:, :), pointer, intent(IN)             :: KFloor
        integer,                           intent(INOUT)          :: CurrentZ
        !Local-------------------------------------------------------------------------------
        integer                                                   :: k, KBottom
        !------------------------------------------------------------------------------------
        KBottom = KFloor(I, J)
        do k = KBottom, KUB
            CurrentZ = CurrentZ + 1
            Matrix(CurrentZ, 1) = I
            Matrix(CurrentZ, 2) = J
            Matrix(CurrentZ, 3) = k
        enddo
    end subroutine UpdateDischargeConnections
    !
    !
    end module ModuleUpscallingDischarges

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon.
!----------------------------------------------------------------------------------------------------------