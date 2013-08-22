!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Gauge
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to calculate the waterlevel along the open boundaries
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

Module ModuleOpenBoundary

    use ModuleGlobalData
    use ModuleTime  
    use ModuleTriangulation          
    use ModuleHorizontalGrid
    use ModuleHorizontalMap,    only : GetBoundaryFaces, UnGetHorizontalMap,GetBoundaries
    use ModuleGauge
    use ModuleFunctions,        only : RodaXY, CHUNK_J
    use ModuleStopWatch,        only : StartWatch, StopWatch
    
    implicit none
    
    private

    !Subroutines---------------------------------------------------------------
    
    !Constructor
    public  :: ConstructOpenBoundary
    private ::      AllocateInstance

    !Selector
    public  :: GetImposedElevation
    public  :: GetOldImposedElevation
    public  :: GetAverageImposedElevation
    public  :: GetBoundaryReferenceLevel
    public  :: GetOpenBoundParameter
    public  :: GetImposedVelocity


    public  :: UnGetOpenBoundary


    !Modifier
    public  :: Modify_OpenBoundary


    !Destructor
    public  :: KillOpenBoundary
    private ::      DeallocateInstance


    !Management
    private ::  Ready
    private ::      LocateObjOpenBoundary

    !Interfaces----------------------------------------------------------------

    private :: UnGetOpenBoundary2D
    interface  UnGetOpenBoundary
        module procedure UnGetOpenBoundary2D
    end interface  UnGetOpenBoundary

    !Parameter
    integer, parameter :: DirectionX_ = 1, DirectionY_ = 2

    !Type----------------------------------------------------------------------
    type T_Station
        real, dimension(:), pointer     :: Elevation      => null() !inicialization: Carina
        real, dimension(:), pointer     :: OpenPoints     => null() !inicialization: Carina
        real, dimension(:), pointer     :: ReferenceLevel => null() !inicialization: Carina
        !Velocities measured by a current meter
        real, dimension(:), pointer     :: VelocityU      => null() !inicialization: Carina
        real, dimension(:), pointer     :: VelocityV      => null() !inicialization: Carina
        !Station position 
        real, dimension(:), pointer     :: MetricX        => null() !inicialization: Carina
        real, dimension(:), pointer     :: MetricY        => null() !inicialization: Carina
    end type T_Station


    type      T_OpenBoundary
        integer                             :: InstanceID   = null_int !inicialization: Carina

        type (T_Time)                       :: StartTime
        type (T_Size2D)                     :: Size, WorkSize

        type (T_Station)                    :: Station

        !Elevation along the boundary
        real, pointer, dimension(:, :)      :: ImposedElevation        => null() !inicialization: Carina

        !Old Elevation along the boundary
        real, pointer, dimension(:, :)      :: OldImposedElevation     => null() !inicialization: Carina
        
        !Average Elevation imposed in the open boundary
        real                                :: AverageImposedElevation = null_real !inicialization: Carina

        !Boundary reference level
        real, pointer, dimension(:, :)      :: BoundaryReferenceLevel  => null() !inicialization: Carina

        !Initial Reference level common to the interior and to the boundary
        real                                :: InitialReferenceLevel  = null_real !inicialization: Carina

        !Velocity along the open boundary
        real, pointer, dimension(:, :)      :: VelocityU => null() !inicialization: Carina
        real, pointer, dimension(:, :)      :: VelocityV => null() !inicialization: Carina
        
        !Compute Tide
        logical                             :: Compute_tide = .false.  !inicialization: Carina

        !Impose the inverted barometer effect 
        logical                             :: InvertBarometer      = .false.  !inicialization: Carina
        logical                             :: InvertBaromSomeBound = .false.  !inicialization: Carina
        real, pointer, dimension(:,:)       :: InvertBarometerCells => null() !inicialization: Carina     

        !This coefficient is zero if the water level is the real one
        !Any slow start is consider in this case
        real                                :: SlowStartCoef = null_real !inicialization: Carina     
        
        logical                             :: TriangGaugesON = .true.

        !Instance of other modules
        integer                             :: ObjTime           = 0
        integer                             :: ObjHorizontalGrid = 0
        integer                             :: ObjHorizontalMap  = 0
        integer                             :: ObjGauge          = 0
        integer                             :: ObjTriangulation  = 0

        type (T_OpenBoundary), pointer      :: Next

    end type T_OpenBoundary

    !Global Module Variables
    type (T_OpenBoundary), pointer           :: FirstOpenBoundary => null() !inicialization: Carina   
    type (T_OpenBoundary), pointer           :: Me    => null() !inicialization: Carina   

    !--------------------------------------------------------------------------
    
    contains    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructOpenBoundary(OpenBoundaryID,                          &
                                     HorizontalGridID,                        &
                                     HorizontalMapID,                         &
                                     TimeID,                                  & 
                                     Compute_Tide,                            & 
                                     InitialReferenceLevel,                   &
                                     SlowStartCoef,                           &
                                     InvertBarometer,                         &
                                     InvertBaromSomeBound,                    &
                                     InvertBarometerCells,                    &
                                     STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: OpenBoundaryID
        integer                                     :: HorizontalGridID
        integer                                     :: HorizontalMapID
        integer                                     :: TimeID
        logical, intent(IN)                         :: Compute_Tide
        real,    intent(IN)                         :: InitialReferenceLevel
        real,    intent(IN)                         :: SlowStartCoef
        logical, intent(IN)                         :: InvertBarometer
        logical, intent(IN)                         :: InvertBaromSomeBound
        real,    pointer, dimension(:,:)            :: InvertBarometerCells                                      
        integer, optional, intent(OUT)              :: STAT     

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        integer                                     :: NGauges
        logical                                     :: TriangGaugesON

        type(T_Time) :: CurrentTime  

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mOpenBoundary_)) then
            nullify (FirstOpenBoundary)
            call RegisterModule (mOpenBoundary_) 
        endif

        call Ready(OpenBoundaryID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Associates other instances
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )


            !Nullify Water level variables
            nullify(Me%Station%Elevation)
            nullify(Me%Station%ReferenceLevel)
            nullify(Me%Station%MetricX)
            nullify(Me%Station%MetricY)
            nullify(Me%Station%OpenPoints)
            nullify(Me%Station%VelocityU)
            nullify(Me%Station%VelocityV)

            !Imposed elevation 
            nullify(Me%ImposedElevation)

            nullify(Me%OldImposedElevation)

            nullify(Me%VelocityU)
            nullify(Me%VelocityV)
            
            nullify(Me%BoundaryReferenceLevel)


            !Time Properties - Actualizes CurrentTime
            call GetComputeCurrentTime(Me%ObjTime, CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR17'


            !Gets the size from the Bathymetry
            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                             &
                                       Size     = Me%Size,                               &
                                       WorkSize = Me%WorkSize,                           &
                                       STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR07'

          
            Me%Compute_Tide          = Compute_Tide 
            Me%InitialReferenceLevel = InitialReferenceLevel
            Me%SlowStartCoef         = SlowStartCoef
            Me%StartTime             = CurrentTime
            
            Me%InvertBarometer       = InvertBarometer
            Me%InvertBaromSomeBound  = InvertBaromSomeBound
            Me%InvertBarometerCells => InvertBarometerCells


cd1:        if (Me%Compute_Tide) then   

                !Initialize the Gauges
                call ConstructGauges(GaugeID          = Me%ObjGauge,                     &
                                     HorizontalGridID = Me%ObjHorizontalGrid,            &
                                     TimeID           = Me%ObjTime,                      &
                                     STAT = STAT_CALL)

                if (STAT_CALL .NE. 0) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR04'

                !Get the number of gauges in use
                call GetNGauges(Me%ObjGauge, NGauges, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR05'

                !Allocates Station
                allocate(Me%Station%Elevation(NGauges), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR06'

                allocate(Me%Station%MetricX(NGauges), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR06a'

                allocate(Me%Station%MetricY(NGauges), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR06b'

                allocate(Me%Station%ReferenceLevel(NGauges), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR06c'

                allocate(Me%Station%OpenPoints(NGauges), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR06d'


                allocate(Me%Station%VelocityU(NGauges), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR06e'

                allocate(Me%Station%VelocityV(NGauges), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR06f'
                
                
                call GetTriangGaugesON(Me%ObjGauge,                                     &
                                       TriangGaugesON,                                  &
                                       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR13a'                

                !Starts Triangulation
                if (NGauges > 2 .and. TriangGaugesON) then
                    Me%TriangGaugesON = .true.
                else
                     Me%TriangGaugesON = .false.
                endif
                
                if (Me%TriangGaugesON) then
                    call GetMetricGaugeLocation(Me%ObjGauge,                            &
                                                Me%Station%MetricX,                     &
                                                Me%Station%MetricY,                     &
                                                STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR13'

                    call ConstructTriangulation(Me%ObjTriangulation,                    &
                                                NGauges,                                &
                                                Me%Station%MetricX,                     &
                                                Me%Station%MetricY,                     &
                                                WriteTriangles = ON,                    &
                                                STAT= STAT_CALL)
                    if  (STAT_CALL /= SUCCESS_      ) then
                        stop 'ConstructOpenBoundary - ModuleOpenBoundary - ERR14'
                    endif
                    
                    call WriteTriangles
                    
                endif

            endif cd1

            !Allocates ImposedElevation
            allocate(Me%ImposedElevation    (Me%Size%ILB: Me%Size%IUB,  Me%Size%JLB:Me%Size%JUB))
            allocate(Me%OldImposedElevation (Me%Size%ILB: Me%Size%IUB,  Me%Size%JLB:Me%Size%JUB))
            allocate(Me%VelocityU           (Me%Size%ILB: Me%Size%IUB,  Me%Size%JLB:Me%Size%JUB))
            allocate(Me%VelocityV           (Me%Size%ILB: Me%Size%IUB,  Me%Size%JLB:Me%Size%JUB))

            !By default the exterior level is imposed equal to zero
            Me%ImposedElevation = 0.

            !By default the exterior level is imposed equal to zero
            Me%OldImposedElevation = 0.

            !By default the exterior velocity is imposed equal to zero
            Me%VelocityU = 0.
            Me%VelocityV = 0.


            call ComputeReferenceLevel

            !Returns ID
            OpenBoundaryID = Me%InstanceID
            STAT_          = SUCCESS_

        else 
            
            stop 'ModuleOpenBoundary - ConstructOpenBoundary - ERR99' 

        end if 

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructOpenBoundary

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_OpenBoundary), pointer             :: NewObjOpenBoundary
        type (T_OpenBoundary), pointer             :: PreviousObjOpenBoundary


        !Allocates new instance
        allocate (NewObjOpenBoundary)
        nullify  (NewObjOpenBoundary%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstOpenBoundary)) then
            FirstOpenBoundary            => NewObjOpenBoundary
            Me                           => NewObjOpenBoundary
        else
            PreviousObjOpenBoundary      => FirstOpenBoundary
            Me                           => FirstOpenBoundary%Next
            do while (associated(Me))
                PreviousObjOpenBoundary  => Me
                Me                       => Me%Next
            enddo
            Me                           => NewObjOpenBoundary
            PreviousObjOpenBoundary%Next => NewObjOpenBoundary
        endif

        Me%InstanceID = RegisterNewInstance (mOpenBoundary_)

    end subroutine AllocateInstance
    
    !--------------------------------------------------------------------------
    
    subroutine WriteTriangles

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_CALL, nNodes
        integer                         :: UnitNumber, iT, nTriangles
        real,    dimension(:), pointer  :: XT, YT, ZT
        integer, dimension(:), pointer  :: V1, V2, V3
        character(len=PathLength)       :: TrianglesFileName = "TideTriangulationNetwork.xy"

        !Begin-----------------------------------------------------------------

        !Get the number of triangles
        call GetNumberOfTriangles   (Me%ObjTriangulation, nTriangles)

        !Allocates space for the Triangle vertices and gets them
        allocate(V1(nTriangles))
        allocate(V2(nTriangles))
        allocate(V3(nTriangles))

        call GetTriangleList (Me%ObjTriangulation, v1, v2, v3, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - ModuleOpenBoundary - ERR10'


        !Gets nodes effictive used and the reordered nodes 
        call GetNumberOfNodes (Me%ObjTriangulation, nNodes, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - ModuleOpenBoundary - ERR20'

        allocate(XT(nNodes))
        allocate(YT(nNodes))
        allocate(ZT(nNodes))

        call GetNodesList   (Me%ObjTriangulation, XT, YT, ZT, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - ModuleOpenBoundary - ERR30'


        call UnitsManager (UnitNumber, OPEN_FILE, STAT = STAT_CALL)
        open (unit=UnitNumber, status = 'unknown', file = TrianglesFileName)
        do iT = 1, nTriangles
            write(UnitNumber,*)'<beginpolygon>'
            write(UnitNumber,*)XT(V1(iT)), YT(V1(iT))
            write(UnitNumber,*)XT(V2(iT)), YT(V2(iT))
            write(UnitNumber,*)XT(V3(iT)), YT(V3(iT))
            write(UnitNumber,*)XT(V1(iT)), YT(V1(iT))
            write(UnitNumber,*)'<endpolygon>'
        enddo
        call UnitsManager (UnitNumber, CLOSE_FILE, STAT = STAT_CALL)

        deallocate(XT, YT, ZT, v1, v2, v3)


    end subroutine WriteTriangles

    subroutine ComputeReferenceLevel

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------       
        integer                             :: STAT_CALL
        integer, dimension(:,:), pointer    :: BoundaryPoints2D
        integer, dimension(:,:), pointer    :: BoundaryFacesU2D, BoundaryFacesV2D

        !Local-----------------------------------------------------------------
        integer                             :: i, j, NGauges
        real, dimension(:,:), pointer       :: CoordX, CoordY
        real                                :: PX, PY
        logical                             :: FoundBound

        !----------------------------------------------------------------------                         
        
        if (.not. associated(Me%BoundaryReferenceLevel)) then

            !Allocates ImposedElevation
            allocate(Me%BoundaryReferenceLevel(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        endif


        !Gets BoundaryPoints from the HorizontalMap
        call GetBoundaries(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR10'

         call GetBoundaryFaces(Me%ObjHorizontalMap, BoundaryFacesU = BoundaryFacesU2D, &
                                                    BoundaryFacesV = BoundaryFacesV2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR20'



cd2:    if (Me%Compute_Tide) then

            !Get the number of gauges in use
            call GetNGauges(Me%ObjGauge, NGauges, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR10'


            !Get the current elevation at the gauges
            call GetReferenceLevel(Me%ObjGauge,                         &
                                   Me%Station%ReferenceLevel,           &
                                   STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR20'


            !If there are less then three gauges, only the first is considered
cd3:        if (NGauges < 3) then

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

!                    if ( BoundaryPoints2D(i, j) == Boundary ) then
                    if ( BoundaryFacesU2D(i  , j  ) == Boundary .or. &
                         BoundaryFacesU2D(i  , j+1) == Boundary .or. &
                         BoundaryFacesV2D(i  , j  ) == Boundary .or. &
                         BoundaryFacesV2D(i+1, j  ) == Boundary) then
                        Me%BoundaryReferenceLevel(i, j) = Me%Station%ReferenceLevel(1)
                    else
                        Me%BoundaryReferenceLevel(i, j) = 0.
                    endif

                enddo
                enddo
            
                if (NGauges == 2) write(*,*)'Warning - Second gauge ignored'

            else cd3

                !Gets Center cell coordinates
                call GetZCoordinates(Me%ObjHorizontalGrid,  CoordX = CoordX, CoordY = CoordY, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR30'


                !Interpolates Elevation at the boundary points
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    FoundBound = .false.

                    if (Me%TriangGaugesON) then
                        if (BoundaryFacesU2D(i  , j  ) == Boundary .or.                     &
                            BoundaryFacesU2D(i  , j+1) == Boundary .or.                     &
                            BoundaryFacesV2D(i  , j  ) == Boundary .or.                     &
                            BoundaryFacesV2D(i+1, j  ) == Boundary) FoundBound = .true.

                    else
                        if (BoundaryPoints2D(i, j) == Boundary) FoundBound = .true.
                    endif
                         
cd4:                if (FoundBound) then

                        !Points where to interpolate
                        PX = CoordX(i,j)
                        PY = CoordY(i,j)
                            
                        call GetIJReferenceLevel(Me%ObjGauge, Me%ObjHorizontalGrid,      &
                                                 i, j,                                   &
                                                 Me%BoundaryReferenceLevel(i, j),        &
                                                 STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= NOT_FOUND_ERR_)     &
                            stop "ComputeReferenceLevel - ModuleOpenBoundary - ERR40"

cd24:                   if (STAT_CALL == NOT_FOUND_ERR_) then 
                          
                            if (Me%TriangGaugesON) then                               
                                call SetHeightValues (Me%ObjTriangulation,  &
                                                      Me%Station%ReferenceLevel, &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR60'

                                !Interpolation
                                Me%BoundaryReferenceLevel(i, j) = InterPolation(             &
                                                                  Me%ObjTriangulation,       &
                                                                  PX, PY,                    &
                                                                  FillOutsidePoints=.true.,  &
                                                                  STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR70'
                            else
                                write(*,*) 'Triangulation is OFF'
                                write(*,*) 'It is missing a gauge in cell (i,j)',i,j,' coordinates (x,y)', PX, PY
                                stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR80'
                            endif
                            
                        endif cd24

                    else cd4

                        Me%BoundaryReferenceLevel(i, j) = 0.

                    endif cd4
                enddo
                enddo

                !Ungets CoordX and CoordY
                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, stat = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR80'

                call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, stat = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR90'

            endif cd3

        else cd2

            Me%BoundaryReferenceLevel(:,:) = Me%InitialReferenceLevel

        endif cd2

        !Unget boundary points 2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR13'

        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesU2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR100'

        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesV2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeReferenceLevel - ModuleOpenBoundary - ERR110'


        !----------------------------------------------------------------------

    end subroutine ComputeReferenceLevel



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetImposedElevation(OpenBoundaryID, ImposedElevation, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: OpenBoundaryID  
        real, dimension(:,:), pointer      :: ImposedElevation
        integer, optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_
         
        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mOPENBOUNDARY_, Me%InstanceID)
            ImposedElevation => Me%ImposedElevation

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetImposedElevation

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine GetImposedVelocity(OpenBoundaryID, VelocityUV, Direction, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: OpenBoundaryID  
        real, dimension(:,:), pointer      :: VelocityUV
        integer                            :: Direction
        integer, optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_
         
        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mOPENBOUNDARY_, Me%InstanceID)

            if (Direction == DirectionX_) then

                VelocityUV => Me%VelocityU

                STAT_ = SUCCESS_

            else if (Direction == DirectionY_) then

                VelocityUV => Me%VelocityV

                STAT_ = SUCCESS_

            endif 

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

    end subroutine GetImposedVelocity

    !--------------------------------------------------------------------------

    subroutine GetOpenBoundParameter(OpenBoundaryID, DirectionX, DirectionY, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: OpenBoundaryID  
        integer, optional                  :: DirectionX, DirectionY
        integer, optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                            :: STAT_, ready_
         
        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(DirectionX)) DirectionX = DirectionX_
            if (present(DirectionY)) DirectionY = DirectionY_

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetOpenBoundParameter

    !--------------------------------------------------------------------------

    subroutine GetOldImposedElevation(OpenBoundaryID, OldImposedElevation, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: OpenBoundaryID  
        real, dimension(:,:), pointer      :: OldImposedElevation
        integer, optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                            :: STAT_, ready_
         
        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mOPENBOUNDARY_, Me%InstanceID)
            OldImposedElevation => Me%OldImposedElevation

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetOldImposedElevation

    !--------------------------------------------------------------------------

    subroutine GetBoundaryReferenceLevel(OpenBoundaryID, BoundaryReferenceLevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: OpenBoundaryID  
        real, dimension(:,:), pointer      :: BoundaryReferenceLevel
        integer, optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_
         
        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mOPENBOUNDARY_, Me%InstanceID)
            BoundaryReferenceLevel => Me%BoundaryReferenceLevel

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryReferenceLevel

    !--------------------------------------------------------------------------

    subroutine GetAverageImposedElevation(OpenBoundaryID, AverageImposedElevation, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: OpenBoundaryID  
        real                               :: AverageImposedElevation
        integer, optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_, ready_
         
        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            AverageImposedElevation = Me%AverageImposedElevation

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) &
            STAT = STAT_

    end subroutine GetAverageImposedElevation

    !--------------------------------------------------------------------------

    subroutine UngetOpenBoundary2D(OpenBoundaryID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: OpenBoundaryID  
        real, pointer, dimension(:,:)      :: Array
        integer, optional, intent(OUT)     :: STAT
         
        !Local-----------------------------------------------------------------
        integer                            :: ready_, STAT_           

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mOPENBOUNDARY_, Me%InstanceID, "UngetOpenBoundary2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetOpenBoundary2D


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Modify_OpenBoundary(OpenBoundaryID, time_, AtmosphericPressure,          &
                                   AtmosphericCoef, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: OpenBoundaryID  
        type(T_Time), intent(IN)            :: time_
        real,  dimension(:,:), pointer      :: AtmosphericPressure
        real   , intent(IN)                 :: AtmosphericCoef
        integer, optional, intent(OUT)      :: STAT  

        !External--------------------------------------------------------------
        integer                             :: ready_             
        integer                             :: STAT_CALL
        integer, dimension(:,:), pointer    :: BoundaryPoints2D
        integer, dimension(:,:), pointer    :: BoundaryFacesU2D, BoundaryFacesV2D

        !Local-----------------------------------------------------------------
        integer                             :: STAT_ 
        integer                             :: i, j, ii, jj, NGauges, NOpen
        real, dimension(:,:), pointer       :: CoordX, CoordY
        real                                :: PX, PY, Counter, SumX
        real                                :: WaterLevel, DT_Run, Coef,                 &
                                               RefLevel
        real, dimension(:), pointer         :: AuxElevation, AuxRefLevel,                &
                                               AuxMetricX  , AuxMetricY,                 &
                                               AuxVelocityU, AuxVelocityV 
        logical                             :: FoundBound 
        !$ integer                             :: CHUNK  

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !$ CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                    
            Me%OldImposedElevation = Me%ImposedElevation
            
            Me%ImposedElevation(:,:) = 0.

            !Gets BoundaryPoints from the HorizontalMap
            call GetBoundaries(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR01'

             call GetBoundaryFaces(Me%ObjHorizontalMap, BoundaryFacesU = BoundaryFacesU2D, &
                                                        BoundaryFacesV = BoundaryFacesV2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR01'

            if (MonitorPerformance)                                     &
                call StartWatch ("ModuleOpenBoundary", "Modify_OpenBoundary")

cd2:        if (Me%Compute_Tide) then

                DT_Run = Time_ - Me%StartTime 

                if (DT_Run < Me%SlowStartCoef) then

                    Coef = (DT_Run / Me%SlowStartCoef) ** 2. 

                else

                    Coef = 1.

                endif


                !Get the number of gauges in use
                call GetNGauges(Me%ObjGauge, NGauges, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR02'


                !Get the current elevation at the gauges
                call GaugeLevel(Me%ObjGauge,                                &
                                Me%Station%Elevation,                       &
                                Me%Station%OpenPoints,                      &
                                Time_,                                                   &    
                                ReferenceLevel = Me%Station%ReferenceLevel, &
                                VelocityU      = Me%Station%VelocityU,      &
                                VelocityV      = Me%Station%VelocityV,      &
                                STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR03'

                !If there are less then three gauges, only the first is considered
cd3:            if (NGauges < 3) then
                    
                    !$OMP PARALLEL PRIVATE(i,j)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

!                        if ( BoundaryPoints2D(i, j) == Boundary ) then
                        if ( BoundaryFacesU2D(i  , j  ) == Boundary .or. &
                             BoundaryFacesU2D(i  , j+1) == Boundary .or. &
                             BoundaryFacesV2D(i  , j  ) == Boundary .or. &
                             BoundaryFacesV2D(i+1, j  ) == Boundary) then


                            Me%ImposedElevation(i, j) =                     &
                                Me%Station%Elevation(1)      *  Coef +      &
                                Me%Station%ReferenceLevel(1) 

                            Me%VelocityU(i, j) =                            &
                                Me%Station%VelocityU(1)      *  Coef 

                            Me%VelocityV(i, j) =                            &
                                Me%Station%VelocityV(1)      *  Coef 


                        else

                            Me%ImposedElevation(i, j) = 0.
                            Me%VelocityU       (i, j) = 0.
                            Me%VelocityV       (i, j) = 0.

                        endif
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL
            
                    if (NGauges == 2) stop 'Warning - Second gauge ignored'

                else cd3
                
                    !Gets Center cell coordinates
                    call GetZCoordinates(Me%ObjHorizontalGrid,  CoordX = CoordX, CoordY = CoordY, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR06'

                    NOpen = Sum(Me%Station%OpenPoints(1:NGauges))
                    if (NOpen < 3) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR07'

                    allocate (AuxElevation (1 : NOpen)) 
                    allocate (AuxRefLevel  (1 : NOpen)) 
                    allocate (AuxMetricX   (1 : NOpen)) 
                    allocate (AuxMetricY   (1 : NOpen)) 
                    allocate (AuxVelocityU (1 : NOpen)) 
                    allocate (AuxVelocityV (1 : NOpen)) 

                    jj = 1

                    do ii = 1, NGauges
                        
cd4:                    if (Me%Station%OpenPoints(ii) == 1) then

                            AuxElevation (jj) = Me%Station%Elevation     (ii)

                            AuxRefLevel  (jj) = Me%Station%ReferenceLevel(ii)

                            AuxMetricX   (jj) = Me%Station%MetricX       (ii)

                            AuxMetricY   (jj) = Me%Station%MetricY       (ii)
                        
                            jj = jj + 1

                        endif cd4

                    enddo

cd5:                if (NOpen < NGauges .and. Me%TriangGaugesON) then

                        call KillTriangulation(Me%ObjTriangulation, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR08'

                        call ConstructTriangulation(Me%ObjTriangulation,                 &
                                                    NOpen,                               &
                                                    AuxMetricX, AuxMetricY,              &
                                                    WriteTriangles = ON,                 &
                                                    STAT = STAT_CALL)
                        if  (STAT_CALL /= SUCCESS_) then
                            stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR09'
                        endif
                        
                        call WriteTriangles

                    endif cd5
                    
                    if (Me%TriangGaugesON) then
                    
                        call SetHeightValues (Me%ObjTriangulation,  &
                                              AuxElevation,                      &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR11'
                                    
                    endif
                  
                    !Interpolates Elevation at the boundary points

                    !$OMP PARALLEL PRIVATE( i,j,FoundBound,PX,PY,WaterLevel,STAT_CALL)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        FoundBound = .false.

                        if (Me%TriangGaugesON) then
                            if (BoundaryFacesU2D(i  , j  ) == Boundary .or.             &
                                BoundaryFacesU2D(i  , j+1) == Boundary .or.             &
                                BoundaryFacesV2D(i  , j  ) == Boundary .or.             &
                                BoundaryFacesV2D(i+1, j  ) == Boundary) FoundBound = .true.

                        else
                            if (BoundaryPoints2D(i, j) == Boundary) FoundBound = .true.
                        endif
                        
cd6:                    if (FoundBound) then

                            !Points where to interpolate
                            PX = CoordX(i,j)
                            PY = CoordY(i,j)

                            STAT_CALL = UNKNOWN_
                            
                            call GetIJWaterLevel_ThreadSafe(Me%ObjGauge, Me%ObjHorizontalGrid, i, j, WaterLevel, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= NOT_FOUND_ERR_) then
                                write(*,*) ' STAT_CALL is ', STAT_CALL
                                stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR10'
                            endif

cd23:                       if (STAT_CALL == NOT_FOUND_ERR_) then

                                if (Me%TriangGaugesON) then

                                    STAT_CALL = UNKNOWN_

                                    WaterLevel = InterPolation_ThreadSafe(Me%ObjTriangulation,  &
                                                               PX, PY, FillOutsidePoints=.true., &
                                                               STAT = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_) then
                                        write(*,*) ' STAT_CALL is ', STAT_CALL
                                        stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR12'
                                    endif

                                else
                                    write(*,*) 'Triangulation is OFF'
                                    write(*,*) 'It is missing a gauge in cell (i,j)',i,j,' coordinates (x,y)', PX, PY
                                    stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR80'
                                endif                                

                            endif cd23

                            !Interpolation
                            Me%ImposedElevation(i, j) = WaterLevel * Coef

                        endif cd6

                    enddo
                    enddo
                    !$OMP END DO NOWAIT
                    !$OMP END PARALLEL

                    if (Me%TriangGaugesON) then

                        call SetHeightValues (  Me%ObjTriangulation,        &
                                                AuxRefLevel,                &
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR14'

                    endif

                    !$OMP PARALLEL PRIVATE(i,j,FoundBound,PX,PY,RefLevel,STAT_CALL)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        FoundBound = .false.

                        if (Me%TriangGaugesON) then
                            if (BoundaryFacesU2D(i  , j  ) == Boundary .or.             &
                                BoundaryFacesU2D(i  , j+1) == Boundary .or.             &
                                BoundaryFacesV2D(i  , j  ) == Boundary .or.             &
                                BoundaryFacesV2D(i+1, j  ) == Boundary) FoundBound = .true.

                        else
                            if (BoundaryPoints2D(i, j) == Boundary) FoundBound = .true.
                        endif
                        
                        if (FoundBound) then

                            !Points where to interpolate
                            PX = CoordX(i,j)
                            PY = CoordY(i,j)

                            STAT_CALL = UNKNOWN_
                            
                            call GetIJReferenceLevel_ThreadSafe(Me%ObjGauge, Me%ObjHorizontalGrid, i, j, RefLevel, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= NOT_FOUND_ERR_) then                                
                                 write(*,*) 'STAT_CALL is ', STAT_CALL
                                stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR13'
                            endif

                          if (STAT_CALL == NOT_FOUND_ERR_) then                               

                                if (Me%TriangGaugesON) then

                                    STAT_CALL = UNKNOWN_

                                    RefLevel = InterPolation_ThreadSafe(Me%ObjTriangulation, &
                                                             PX, PY, FillOutsidePoints=.true., STAT = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_) then
                                        write(*,*) 'STAT_CALL is ', STAT_CALL
                                        stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR15'
                                    endif
                                    
                                else
                                    write(*,*) 'Triangulation is OFF'
                                    write(*,*) 'It is missing a gauge in cell (i,j)',i,j,' coordinates (x,y)', PX, PY
                                    stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR80'
                                endif                                

                            endif

                            !Interpolation
                            Me%ImposedElevation(i, j) = Me%ImposedElevation(i, j) + RefLevel

                        else

                            Me%ImposedElevation(i, j) = 0.

                        endif
                        
                    enddo
                    enddo
                    !$OMP END DO NOWAIT                    
                    !$OMP END PARALLEL

                    deallocate (AuxElevation) 
                    deallocate (AuxRefLevel ) 
                    deallocate (AuxMetricX  ) 
                    deallocate (AuxMetricY  ) 
                    deallocate (AuxVelocityU) 
                    deallocate (AuxVelocityV) 

                    !Ungets CoordX and CoordY
                    call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, stat = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR16'

                    call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, stat = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR17'

                endif cd3
            
            endif cd2

            Counter = 0.
            SumX    = 0.

            !$OMP PARALLEL PRIVATE(i,j,FoundBound)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK) REDUCTION(+:Counter)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                FoundBound = .false.
                
                if (Me%TriangGaugesON) then
                    if (BoundaryFacesU2D(i  , j  ) == Boundary .or.                     &
                        BoundaryFacesU2D(i  , j+1) == Boundary .or.                     &
                        BoundaryFacesV2D(i  , j  ) == Boundary .or.                     &
                        BoundaryFacesV2D(i+1, j  ) == Boundary) FoundBound = .true.

                else
                    if (BoundaryPoints2D(i, j) == Boundary) FoundBound = .true.
                endif

                if (FoundBound) then

                    !$OMP CRITICAL (MOP_RED01)
                    SumX    = SumX + Me%ImposedElevation(i, j)
                    !$OMP END CRITICAL (MOP_RED01)
                    Counter = Counter + 1.

                endif
                    
                
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL   

            if (MonitorPerformance)                                     &
                call StopWatch ("ModuleOpenBoundary", "Modify_OpenBoundary")


            if (Counter == 0) then
            
                Me%AverageImposedElevation = 0.

            else

                Me%AverageImposedElevation = SumX/Counter

            endif
            
            if (Me%InvertBarometer .and. associated(AtmosphericPressure))               &
                    call ImposeInvertBarometer(AtmosphericPressure, AtmosphericCoef)

            !Unget boundary points 2D
            call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR18'

            call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesU2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR18'

            call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesV2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Modify_OpenBoundary - ModuleOpenBoundary - ERR19'



            STAT_ = SUCCESS_
        else               

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine Modify_OpenBoundary

    !--------------------------------------------------------------------------

    subroutine ImposeInvertBarometer(AtmosphericPressure, AtmosphericCoef)

        !Arguments-------------------------------------------------------------
        real,      dimension(:,:), pointer      :: AtmosphericPressure
        real                                    :: AtmosphericCoef

        !Local-----------------------------------------------------------------
        integer,   dimension(:,:), pointer      :: BoundaryFacesU2D, BoundaryFacesV2D
        
        integer                                 :: i, j, STAT_CALL
        
        !$ integer                              :: CHUNK

        !----------------------------------------------------------------------                         

        call GetBoundaryFaces(Me%ObjHorizontalMap, BoundaryFacesU = BoundaryFacesU2D, &
                                                   BoundaryFacesV = BoundaryFacesV2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ImposeInvertBarometer - ModuleOpenBoundary - ERR10'

        !$ CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j)

        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if ( BoundaryFacesU2D(i  , j  ) == Boundary .or. &
                 BoundaryFacesU2D(i  , j+1) == Boundary .or. &
                 BoundaryFacesV2D(i  , j  ) == Boundary .or. &
                 BoundaryFacesV2D(i+1, j  ) == Boundary) then
                                 
                if (Me%InvertBaromSomeBound) then 
                    if(Me%InvertBarometerCells(i, j) > 0. ) cycle
                endif
                    
                !Inverted barometer effect
                Me%ImposedElevation(i, j) = Me%ImposedElevation(i, j) + AtmosphericCoef * &
                                            (AtmPressSeaLevelReference - AtmosphericPressure(i,j)) /( 1.e3 * Gravity) 
                             
            endif

        enddo
        enddo
        !$OMP END DO

        !$OMP END PARALLEL

        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesU2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ImposeInvertBarometer - ModuleOpenBoundary - ERR20'

        call UnGetHorizontalMap(Me%ObjHorizontalMap, BoundaryFacesV2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ImposeInvertBarometer - ModuleOpenBoundary - ERR30'
        !----------------------------------------------------------------------

    end subroutine ImposeInvertBarometer


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillOpenBoundary(OpenBoundaryID, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: OpenBoundaryID
        integer, optional, intent(OUT)  :: STAT    

        !Local-----------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: STAT_, ready_, nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(OpenBoundaryID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mOPENBOUNDARY_,  Me%InstanceID)

            if (nUsers == 0) then

                !Deassociates External Instances
                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime          )
                if (nUsers == 0) stop 'KillOpenBoundary - OpenBoundary - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap )
                if (nUsers == 0) stop 'KillOpenBoundary - OpenBoundary - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillOpenBoundary - OpenBoundary - ERR03'


cd4:            if (Me%Compute_Tide) then 

                    !Deallocates ObjGauge
cd3 :               if (Me%ObjGauge /= 0) then
                        call KillGauge(Me%ObjGauge, STAT = STAT_CALL) 
                        if (STAT_CALL /= SUCCESS_) stop 'KillOpenBoundary - OpenBoundary - ERR04'
                    end if cd3


                    !Deallocates Triangulation
                    if (Me%ObjTriangulation /= 0) then
                        call KillTriangulation(Me%ObjTriangulation, STAT = STAT_CALL) 
                        if (STAT_CALL /= SUCCESS_) stop 'KillOpenBoundary - OpenBoundary - ERR05'
                    endif


                    deallocate(Me%Station%Elevation)
                    deallocate(Me%Station%ReferenceLevel)
                    deallocate(Me%Station%MetricX)
                    deallocate(Me%Station%MetricY)
                    deallocate(Me%Station%VelocityU)
                    deallocate(Me%Station%VelocityV)

                    nullify(Me%Station%Elevation)
                    nullify(Me%Station%ReferenceLevel)
                    nullify(Me%Station%MetricX)
                    nullify(Me%Station%MetricY)
                    nullify(Me%Station%VelocityU)
                    nullify(Me%Station%VelocityV)


                endif cd4


                deallocate(Me%ImposedElevation)
                nullify(Me%ImposedElevation)

                deallocate(Me%OldImposedElevation)
                nullify(Me%OldImposedElevation)

                deallocate(Me%BoundaryReferenceLevel)
                nullify(Me%BoundaryReferenceLevel)

                call DeallocateInstance

                OpenBoundaryID = 0
                STAT_          = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


        !----------------------------------------------------------------------

    end subroutine KillOpenBoundary

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_OpenBoundary), pointer          :: AuxObjOpenBoundary
        type (T_OpenBoundary), pointer          :: PreviousObjOpenBoundary

        !Updates pointers
        if (Me%InstanceID == FirstOpenBoundary%InstanceID) then
            FirstOpenBoundary => FirstOpenBoundary%Next
        else
            PreviousObjOpenBoundary => FirstOpenBoundary
            AuxObjOpenBoundary      => FirstOpenBoundary%Next
            do while (AuxObjOpenBoundary%InstanceID /= Me%InstanceID)
                PreviousObjOpenBoundary => AuxObjOpenBoundary
                AuxObjOpenBoundary      => AuxObjOpenBoundary%Next
            enddo

            !Now update linked list
            PreviousObjOpenBoundary%Next => AuxObjOpenBoundary%Next

        endif
            
        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

    end subroutine DeallocateInstance


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (ObjOpenBoundary_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjOpenBoundary_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjOpenBoundary_ID > 0) then
            call LocateObjOpenBoundary(ObjOpenBoundary_ID)
            ready_ = VerifyReadLock (mOpenBoundary_,  Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjOpenBoundary (ObjOpenBoundaryID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjOpenBoundaryID

        !Local-----------------------------------------------------------------

        Me => FirstOpenBoundary
        do while (associated (Me))
            if (Me%InstanceID == ObjOpenBoundaryID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleOpenBoundary - LocateObjOpenBoundary - ERR01'

    end subroutine LocateObjOpenBoundary

    !--------------------------------------------------------------------------

end Module ModuleOpenBoundary

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
