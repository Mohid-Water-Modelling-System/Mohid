!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : RunOff
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jan 2004
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which calculates the Surface RunOff
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

Module ModuleRunOff

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHDF5
    use ModuleFunctions         ,only : TimeToString, SetMatrixValue, ChangeSuffix,      &
                                        CHUNK_J, LinearInterpolation
    use ModuleHorizontalGrid    ,only : GetHorizontalGridSize, GetHorizontalGrid,        &
                                        UnGetHorizontalGrid, WriteHorizontalGrid,        &
                                        GetGridCellArea, GetXYCellZ,                     &
                                        GetCellZInterceptByLine,                         &
                                        GetCellZInterceptByPolygon
    use ModuleHorizontalMap     ,only : GetBoundaries, UngetHorizontalMap
    use ModuleGridData          ,only : GetGridData, UngetGridData, WriteGridData
    use ModuleBasinGeometry     ,only : GetBasinPoints, GetRiverPoints, GetCellSlope,    &
                                        GetDrainageDirection, TargetPoint,               &
                                        UnGetBasin
    use ModuleStopWatch         ,only : StartWatch, StopWatch
    use ModuleFillMatrix        ,only : ConstructFillMatrix, ModifyFillMatrix,           &
                                        KillFillMatrix
    use ModuleDrainageNetwork   ,only : GetChannelsWaterLevel, GetChannelsSurfaceWidth,  &
                                        GetChannelsBankSlope, GetChannelsNodeLength,     &
                                        GetChannelsBottomLevel, UnGetDrainageNetwork,    &
                                        GetChannelsID,GetChannelsVolume,                 &
                                        GetChannelsMaxVolume, GetChannelsActiveState,    &
                                        GetChannelsTopArea
    use ModuleDischarges        ,only : Construct_Discharges, GetDischargesNumber,       &
                                        GetDischargesGridLocalization,                   &
                                        GetDischargeWaterFlow, GetDischargesIDName,      &
                                        TryIgnoreDischarge, GetDischargeSpatialEmission, &
                                        CorrectsCellsDischarges, Kill_Discharges
                                        
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  ::  ConstructRunOff
    private ::      AllocateInstance
    private ::      ReadDataFile
    private ::      AllocateVariables
    private ::      ConstructOverLandCoefficient
    private ::      ConstructStormWaterDrainage
    private ::      ConstructHDF5Output

    !Selector
    public  ::  GetOverLandFlow
    public  ::  GetManning
    public  ::  GetManningDelta
    public  ::  GetFlowToChannels
    !public  ::  GetFlowAtBoundary
    public  ::  GetBoundaryImposed
    public  ::  GetRouteDFour
    public  ::  GetRouteDFourCells
    public  ::  GetRouteDFourNeighbours
    public  ::  GetRouteDFourFlux
    public  ::  GetBoundaryFlux
    public  ::  GetBoundaryCells
    public  ::  GetFlowDischarge
    public  ::  GetRunoffWaterLevel 
    public  ::  GetRunoffWaterColumn        !Final WaterColumn 
    public  ::  GetRunoffWaterColumnOld     !Initial WaterColumn
    public  ::  GetRunoffWaterColumnAT      !WaterColumn After Transport (For RP) 
    public  ::  GetRunoffCenterVelocity
    public  ::  GetRunoffTotalStoredVolume
    public  ::  GetRunOffStoredVolumes
    public  ::  GetRunOffBoundaryFlowVolume
    public  ::  GetMassError
    public  ::  GetNextRunOffDT
    public  ::  SetBasinColumnToRunoff
    public  ::  UnGetRunOff
    

    !Modifier
    public  ::  ModifyRunOff
    private ::      LocalWaterColumn
    private ::      IntegrateFlow
    private ::      ComputeNextDT
    private ::      RunOffOutput
    private ::  AdjustSlope

    !Destructor
    public  ::  KillRunOff                                                     
    
    !Management
    private ::  ReadLockExternalVar
    private ::  ReadUnLockExternalVar
    private ::  Ready
    private ::      LocateObjRunOff 

    !Interfaces----------------------------------------------------------------
    private :: UnGetRunOff2D_R4
    private :: UnGetRunOff2D_R8
    interface  UnGetRunOff
        module procedure UnGetRunOff2D_R4
        module procedure UnGetRunOff2D_R8
    end interface  UnGetRunOff
    
    !Parameters----------------------------------------------------------------
    integer, parameter                              :: KinematicWave_   = 1
    integer, parameter                              :: DiffusionWave_   = 2
    integer, parameter                              :: DynamicWave_     = 3      

    integer, parameter                              :: UnitMax          = 80
    
    !water column computation in faces
    integer, parameter                              :: WCMaxBottom_     = 1
    integer, parameter                              :: WCAverageBottom_ = 2
    
    !Boundary flux
    integer, parameter                              :: ComputeFlow_       = 1
    integer, parameter                              :: InstantaneousFlow_ = 2

    !Route D4 flux
    integer, parameter                              :: Celerity_          = 1
    integer, parameter                              :: Manning_           = 2    
    
    !Types---------------------------------------------------------------------
    type T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime              => null()
         integer                                    :: NextOutPut           = 1
         logical                                    :: Yes                  = .false.
         type (T_Time), dimension(:), pointer       :: RestartOutTime       => null()
        logical                                     :: WriteRestartFile     = .false.
        logical                                     :: RestartOverwrite     = .false.
        integer                                     :: NextRestartOutput    = 1         
    end type T_OutPut


    type T_Files
        character(PathLength)                       :: DataFile             = null_str
        character(PathLength)                       :: InitialFile          = null_str
        character(PathLength)                       :: FinalFile            = null_str
        character(PathLength)                       :: TransientHDF         = null_str
    end type T_Files    

    type T_ExtVar
        integer, dimension(:,:), pointer            :: BasinPoints              => null()
        real(8), dimension(:,:), pointer            :: WaterColumn              => null()
        real   , dimension(:,:), pointer            :: GridCellArea             => null()
        real   , dimension(:,:), pointer            :: DUX, DVY                 => null()
        real   , dimension(:,:), pointer            :: DXX, DYY                 => null()
        real   , dimension(:,:), pointer            :: DZX, DZY                 => null()
        real   , dimension(:,:), pointer            :: XX2D_Z, YY2D_Z           => null()
        real   , dimension(:,:), pointer            :: Topography               => null()
        integer, dimension(:,:), pointer            :: BoundaryPoints2D         => null()
        integer, dimension(:,:), pointer            :: RiverPoints              => null()
        real   , dimension(:,:), pointer            :: CellSlope                => null()
        type (T_Time)                               :: Now
        real                                        :: DT                       = null_real
    end type T_ExtVar

    type T_Converge
        integer                                     :: MinIterations                = 1               
        integer                                     :: MaxIterations                = 1024
        logical                                     :: Stabilize                    = .false.
        real                                        :: StabilizeFactor              = 0.01        
        real                                        :: DTFactorUp                   = 1.25
        real                                        :: DTFactorDown                 = 1.25
        real                                        :: StabilizeHardCutLimit        = 128
        real                                        :: DTSplitFactor                = 2.0               
        real                                        :: CurrentDT                    = null_real  
        real                                        :: NextDT                       = null_real
        integer                                     :: LastGoodNiteration           = 1
        integer                                     :: NextNiteration               = 1               
        logical                                     :: LimitDTCourant               = .false.        
        real                                        :: MaxCourant                   = 1.0  
        integer                                     :: MinToRestart                 = 0  
        real                                        :: MinimumValueToStabilize      = 0.001
        logical                                     :: CheckDecreaseOnly            = .false.        
    end type T_Converge
  
    type  T_RunOff
        integer                                     :: InstanceID               = 0
        character(len=StringLength)                 :: ModelName                = null_str
        integer                                     :: ObjBasinGeometry         = 0
        integer                                     :: ObjTime                  = 0
        integer                                     :: ObjHorizontalGrid        = 0
        integer                                     :: ObjHorizontalMap         = 0
        integer                                     :: ObjGridData              = 0
        integer                                     :: ObjHDF5                  = 0
        integer                                     :: ObjDrainageNetwork       = 0
        integer                                     :: ObjDischarges            = 0
        integer                                     :: ObjEnterData             = 0
        type (T_OutPut   )                          :: OutPut
        type (T_ExtVar)                             :: ExtVar
        type (T_Files)                              :: Files
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime
        real(8), dimension(:,:), pointer            :: myWaterLevel             => null()
        real(8), dimension(:,:), pointer            :: myWaterColumn            => null()
        real,    dimension(:,:), pointer            :: InitialWaterColumn       => null()
        real(8), dimension(:,:), pointer            :: myWaterVolume            => null() 
        real(8), dimension(:,:), pointer            :: myWaterColumnOld         => null() !OldColumn from Basin
        real(8), dimension(:,:), pointer            :: myWaterColumnAfterTransport => null() !for property transport
        real(8), dimension(:,:), pointer            :: myWaterVolumePred        => null() !to avoid negative collumns
        real(8), dimension(:,:), pointer            :: myWaterVolumeOld         => null()
        real,    dimension(:,:), pointer            :: lFlowToChannels          => null() !Instantaneous Flow To Channels
        real,    dimension(:,:), pointer            :: iFlowToChannels          => null() !Integrated    Flow
        real,    dimension(:,:), pointer            :: iFlowRouteDFour          => null() !Integrated Route D4 flux
        real,    dimension(:,:), pointer            :: lFlowBoundary            => null() !Instantaneous Flow to impose BC
        real,    dimension(:,:), pointer            :: iFlowBoundary            => null() !Integrated    Flow to impose BC
        real,    dimension(:,:), pointer            :: lFlowDischarge           => null() !Instantaneous Flow of discharges
        real,    dimension(:,:), pointer            :: iFlowDischarge           => null() !Integrated    Flow of discharges
        real(8), dimension(:,:), pointer            :: lFlowX, lFlowY           => null() !Instantaneous OverLandFlow (LocalDT   )
        real(8), dimension(:,:), pointer            :: iFlowX, iFlowY           => null() !Integrated    OverLandFlow (AfterSumDT)
        real(8), dimension(:,:), pointer            :: FlowXOld, FlowYOld       => null() !Flow From previous iteration
        real(8), dimension(:,:), pointer            :: InitialFlowX, InitialFlowY => null() !Initial Flow of convergence
        real,    dimension(:,:), pointer            :: AreaU, AreaV             => null()
        integer, dimension(:,:), pointer            :: ComputeFaceU             => null()
        integer, dimension(:,:), pointer            :: ComputeFaceV             => null()
        real,    dimension(:,:), pointer            :: OverLandCoefficient      => null() !Manning or Chezy
        real,    dimension(:,:), pointer            :: OverLandCoefficientDelta => null() !For erosion/deposition
        real,    dimension(:,:), pointer            :: OverLandCoefficientX     => null() !Manning or Chezy
        real,    dimension(:,:), pointer            :: OverLandCoefficientY     => null() !Manning or Chezy
        real,    dimension(:,:), pointer            :: StormWaterDrainageCoef   => null() !Sewer System Percentagem (area)
        real,    dimension(:,:), pointer            :: StormWaterVolume         => null() !Volume of storm water stored in each cell
        real,    dimension(:,:), pointer            :: StormWaterFlowX          => null() !Auxilizary Var for explicit routing
        real,    dimension(:,:), pointer            :: StormWaterFlowY          => null() !Auxilizary Var for explicit routing
        real,    dimension(:,:), pointer            :: StormWaterCenterFlowX    => null() !Output
        real,    dimension(:,:), pointer            :: StormWaterCenterFlowY    => null() !Output
        real,    dimension(:,:), pointer            :: StormWaterCenterModulus  => null() !Output 
        real,    dimension(:,:), pointer            :: BuildingsHeight          => null() !Height of building in cell
        real,    dimension(:,:), pointer            :: StormWaterInteraction    => null() !Points where interaction with SWMM occurs
        real,    dimension(:,:), pointer            :: StreetGutterLength       => null() !Length of Stret Gutter in a given cell
        real,    dimension(:,:), pointer            :: MassError                => null() !Contains mass error
        real,    dimension(:,:), pointer            :: CenterFlowX              => null()
        real,    dimension(:,:), pointer            :: CenterFlowY              => null()
        real,    dimension(:,:), pointer            :: CenterVelocityX          => null()
        real,    dimension(:,:), pointer            :: CenterVelocityY          => null()
        real,    dimension(:,:), pointer            :: FlowModulus              => null()
        real,    dimension(:,:), pointer            :: VelocityModulus          => null()
        integer, dimension(:,:), pointer            :: LowestNeighborI          => null() !Lowest Neighbor in the surroundings
        integer, dimension(:,:), pointer            :: LowestNeighborJ          => null() !Lowest Neighbor in the surroundings       
        integer, dimension(:,:), pointer            :: DFourSinkPoint           => null() !Point which can't drain with in X/Y only
        integer, dimension(:,:), pointer            :: StabilityPoints          => null() !Points where models check stability
        type(T_PropertyID)                          :: OverLandCoefficientID
        logical                                     :: StormWaterModel          = .false.          !If connected to SWMM
        real,    dimension(:,:), pointer            :: StormWaterModelFlow      => null() !Flow from SWMM
        real,    dimension(:,:), pointer            :: StreetGutterFlow         => null() !Flow through "street gutters"
        real,    dimension(:,:), pointer            :: SewerInflow              => null() !Integrated inflow at sewer points 
                                                                                          !(potential)
        real,    dimension(:,:), pointer            :: StormInteractionFlow     => null() !Interaction Flow 
                                                                                          !(at gutters + sewer points) (real)
        integer, dimension(:,:), pointer            :: StreetGutterTargetI      => null() !Sewer interaction point...
        integer, dimension(:,:), pointer            :: StreetGutterTargetJ      => null() !...where street gutter drains to
        real                                        :: MinSlope              = null_real
        logical                                     :: AdjustSlope           = .false.
        logical                                     :: Stabilize             = .false.
        logical                                     :: Discharges            = .false.
        logical                                     :: RouteDFourPoints      = .false.
        logical                                     :: RouteDFourPointsOnDN  = .false.
        integer                                     :: RouteDFourMethod      = null_int
        logical                                     :: StormWaterDrainage    = .false.
        real                                        :: StormWaterInfiltrationVelocity  = 1.4e-5  !~50mm/h
        real                                        :: StormWaterFlowVelocity          = 0.2     !velocity in pipes
        logical                                     :: Buildings             = .false.
!        real                                        :: StabilizeFactor       = null_real
!        real                                        :: StabilizeHardCutLimit = 0.1
        integer                                     :: HydrodynamicApproximation = DiffusionWave_
        logical                                     :: CalculateAdvection    = .true.
        logical                                     :: CalculateCellMargins  = .true.
        logical                                     :: ImposeMaxVelocity     = .false.
        real                                        :: ImposedMaxVelocity    = 0.1
!        integer                                     :: LastGoodNiter         = 1
!        integer                                     :: NextNiter             = 1
!        real                                        :: InternalTimeStepSplit = 1.5
!        integer                                     :: MinIterations         = 1
!        integer                                     :: MinToRestart     = 0
        real                                        :: MinimumWaterColumn    = null_real
        real                                        :: MinimumWaterColumnAdvection = null_real
!        real                                        :: MinimumWaterColumnStabilize = null_real
!        real                                        :: NextDT               = null_real
!        real                                        :: DTFactor             = null_real
!        real                                        :: DTFactorUp           = null_real
!        real                                        :: DTFactorDown         = null_real
!        real                                        :: CurrentDT            = null_real
!        logical                                     :: LimitDTCourant       = .false.
!        logical                                     :: LimitDTVariation     = .true.
!        real                                        :: MaxCourant           = 1.0        
        logical                                     :: ImposeBoundaryValue  = .false.
        real                                        :: BoundaryValue        = null_real
        real                                        :: MaxDtmForBoundary    = null_real
        integer                                     :: BoundaryMethod       = null_int
        integer, dimension(:,:), pointer            :: BoundaryCells        => null()
!        integer                                     :: MaxIterations        = 5
        logical                                     :: SimpleChannelInteraction = .false.
        logical                                     :: LimitToCriticalFlow  = .true.
        integer                                     :: FaceWaterColumn      = WCMaxBottom_
!        real                                        :: MaxVariation         = null_real
        integer                                     :: OverlandChannelInteractionMethod = null_int
!        logical                                     :: CheckDecreaseOnly    = .false.

        type(T_Converge)                            :: CV !Convergence options
        
        real(8)                                     :: BoundaryFlowVolume        = 0.0 !m3 => positive if flow is towards boundary.          
        real(8)                                     :: VolumeStoredInSurface     = 0.0
        real(8)                                     :: VolumeStoredInStormSystem = 0.0

        logical                                     :: WriteMaxFlowModulus  = .false.
        character(Pathlength)                       :: MaxFlowModulusFile   = null_str
        real, dimension(:,:), pointer               :: MaxFlowModulus       => null()

        logical                                     :: WriteMaxWaterColumn  = .false.        
        character(Pathlength)                       :: MaxWaterColumnFile   = null_str
        real, dimension(:,:), pointer               :: MaxWaterColumn       => null()
        
        logical                                     :: Continuous          = .false.
        logical                                     :: StopOnWrongDate     = .true.
        
        real(8)                                     :: TotalStoredVolume  = 0.
        integer                                     :: BasinCellsCount    = 0

        !Grid size
        type (T_Size2D)                             :: Size
        type (T_Size2D)                             :: WorkSize

        type(T_RunOff), pointer                     :: Next                 => null()
    end type  T_RunOff


    !Global Module Variables
    type (T_RunOff), pointer                        :: FirstObjRunOff       => null()
    type (T_RunOff), pointer                        :: Me                   => null()

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructRunOff(ModelName,                                       &            
                               RunOffID,                                        &
                               ComputeTimeID,                                   &
                               HorizontalGridID,                                &
                               HorizontalMapID,                                 &
                               GridDataID,                                      &
                               BasinGeometryID,                                 &
                               DrainageNetworkID,                               &
                               DischargesID,                                    &
                               STAT)

        !Arguments---------------------------------------------------------------
        character(len=*)                                :: ModelName
        integer                                         :: RunOffID
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: GridDataID
        integer                                         :: BasinGeometryID
        integer                                         :: DrainageNetworkID
        integer, optional, intent(OUT)                  :: STAT     
        integer, intent (OUT)                           :: DischargesID

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------
        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mRunOff_)) then
            nullify (FirstObjRunOff)
            call RegisterModule (mRunOff_) 
        endif

        call Ready(RunOffID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%ModelName = ModelName
            
            !Associates External Instances
            Me%ObjTime            = AssociateInstance (mTIME_           , ComputeTimeID     )
            Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_ , HorizontalGridID  )
            Me%ObjHorizontalMap   = AssociateInstance (mHORIZONTALMAP_  , HorizontalMapID   )
            Me%ObjGridData        = AssociateInstance (mGRIDDATA_       , GridDataID        )
            Me%ObjBasinGeometry   = AssociateInstance (mBASINGEOMETRY_  , BasinGeometryID   )

            if (DrainageNetworkID /= 0) then
                Me%ObjDrainageNetwork   = AssociateInstance (mDRAINAGENETWORK_, DrainageNetworkID)
            endif
            
            !Time Stuff
            call GetComputeTimeLimits   (Me%ObjTime, BeginTime = Me%BeginTime,           &
                                         EndTime = Me%EndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR010'

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR011'
        
            Me%CV%NextNiteration = 1
            Me%CV%CurrentDT = Me%ExtVar%DT

            call ReadLockExternalVar (StaticOnly = .false.)


            !Gets the size of the grid
            call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                        Size     = Me%Size,                              &
                                        WorkSize = Me%WorkSize,                          &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR020'
            
            call AllocateVariables
            
            call ReadDataFile

            call InitializeVariables

            call ConstructOverLandCoefficient

            !Checks if River Network is consistent with the one previously constructed
            if (DrainageNetworkID /= 0) then
                call CheckRiverNetWorkConsistency
            endif

            !Constructs Discharges
            if (Me%Discharges) then
                call ConstructDischarges
            endif
           
            !Constructs StormWaterDrainage
            if (Me%StormWaterDrainage .or. Me%StormWaterModel) then
                call ConstructStormWaterDrainage
            endif
            
            !Constructs Boundary Cells
            if (Me%ImposeBoundaryValue) then
                call CheckBoundaryCells
            endif
            
            !Reads conditions from previous run
            if (Me%Continuous) call ReadInitialFile
            
            if (Me%OutPut%Yes) then
                call ConstructHDF5Output
            endif


            call CalculateTotalStoredVolume

            !Output Results
            if (Me%OutPut%Yes) then
                call ComputeCenterValues                   
                call RunOffOutput
            endif


            call ReadUnLockExternalVar (StaticOnly = .false.)

            !Returns ID
            RunOffID          = Me%InstanceID
            DischargesID      = Me%ObjDischarges

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleRunOff - ConstructRunOff - ERR030' 

        end if cd0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructRunOff
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_RunOff), pointer                         :: NewObjRunOff
        type (T_RunOff), pointer                         :: PreviousObjRunOff


        !Allocates new instance
        allocate (NewObjRunOff)
        nullify  (NewObjRunOff%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjRunOff)) then
            FirstObjRunOff          => NewObjRunOff
            Me                      => NewObjRunOff
        else
            PreviousObjRunOff       => FirstObjRunOff
            Me                      => FirstObjRunOff%Next
            do while (associated(Me))
                PreviousObjRunOff   => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjRunOff
            PreviousObjRunOff%Next  => NewObjRunOff
        endif

        Me%InstanceID = RegisterNewInstance (mRUNOFF_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadDataFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        type(T_PropertyID)                          :: InitialWaterColumnID
        type(T_PropertyID)                          :: OverLandCoefficientDeltaID
        type(T_PropertyID)                          :: StormWaterDrainageID
        type(T_PropertyID)                          :: BuildingsHeightID
        type(T_PropertyID)                          :: StormWaterInteractionID
        type(T_PropertyID)                          :: StreetGutterLengthID
        integer                                     :: iflag, ClientNumber
        logical                                     :: BlockFound
        integer                                     :: i, j
        logical                                     :: DynamicAdjustManning
        real                                        :: dummy

        !Reads the name of the data file from nomfich
        call ReadFileName ('RUNOFF_DATA', Me%Files%DataFile, "RunOff Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR010'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('RUNOFF_HDF', Me%Files%TransientHDF, "RunOff HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR020'

        call ReadFileName('RUNOFF_FIN', Me%Files%FinalFile,                              &
                           Message = "RunOff Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR030'

        !Constructs the DataFile
        call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR040'

        !Initial Water Column
        call GetData(dummy,                                                              &
                     Me%ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'INITIAL_WATER_COLUMN',                              &
                     default      = 0.0,                                                 & 
                     ClientModule = 'ModuleRunOff',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR050'
        
        if (iflag /= 0) then
            write(*,*)'The keyword INITIAL_WATER_COLUMN is obselete.'
            write(*,*)'Please use the block <BeginInitialWaterColumn> / <EndInitialWaterColumn>'
            stop 'ReadDataFile - ModuleRunOff - ERR060'
        endif        

        !Gets Block 
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                          &
                                    '<BeginInitialWaterColumn>',                      &
                                    '<EndInitialWaterColumn>', BlockFound,            &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR070'

        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = InitialWaterColumnID,      &
                                        EnterDataID      = Me%ObjEnterData,                 &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%InitialWaterColumn,        &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR080'

            call KillFillMatrix(InitialWaterColumnID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR090'

        else
            write(*,*)'Missing Block <BeginInitialWaterColumn> / <EndInitialWaterColumn>' 
            stop      'ReadDataFile - ModuleRunOff - ERR100'
        endif

         !Gets Minimum Slope 
        call GetData(Me%MinSlope,                                               &
                     Me%ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'MIN_SLOPE',                                &
                     default      = 0.0,                                        &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0110'

        if (Me%MinSlope < 0.0 .or. Me%MinSlope >= 1.) then
            write (*,*) 'Invalid Minimum Slope [MIN_SLOPE]'
            stop 'ReadDataFile - ModuleRunOff - ERR0120'
        end if

        !Adjusts Slope according to
        !http://www.hkh-friend.net.np/rhdc/training/lectures/HEGGEN/Tc_3.pdf
        call GetData(Me%AdjustSlope,                                            &
                     Me%ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'ADJUST_SLOPE',                             &
                     default      = .true.,                                     &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0130'


        !Gets Routing method
        call GetData(Me%HydrodynamicApproximation,                              &
                     Me%ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'HYDRODYNAMIC_APROX',                       &
                     default      = DiffusionWave_,                             &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0140'

        if (Me%HydrodynamicApproximation /= KinematicWave_ .and.                &
            Me%HydrodynamicApproximation /= DiffusionWave_ .and.                &
            Me%HydrodynamicApproximation /= DynamicWave_) then
            write (*,*) 'Invalid Hydrodynamic Approximation [HYDRODYNAMIC_APROX]'
            stop 'ReadDataFile - ModuleRunOff - ERR0150'
        end if
        
        if (Me%HydrodynamicApproximation == DynamicWave_) then

            !Gets if advection is to be calculated
            call GetData(Me%CalculateAdvection,                                 &
                         Me%ObjEnterData, iflag,                                   &
                         SearchType   = FromFile,                               &
                         keyword      = 'ADVECTION',                            &
                         default      = .true.,                                 &
                         ClientModule = 'ModuleRunOff',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0160'

            
            if (Me%CalculateAdvection) then    
                
                !Minimum Water Column for advection computation
                call GetData(Me%MinimumWaterColumnAdvection,                                     &
                             Me%ObjEnterData, iflag,                                                &
                             SearchType   = FromFile,                                            &
                             keyword      = 'MIN_WATER_COLUMN_ADVECTION',                        &
                             default      = 0.0,                                                 &
                             ClientModule = 'ModuleRunOff',                                      &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0170'
            endif
        endif

        !Method for computing water column in the face (1 - Using max level and max bottom; 
        !2- using max level and average of bottom)
        call GetData(Me%FaceWaterColumn,                                    &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'WATER_COLUMN_FACE',                    &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = WCMaxBottom_,                           &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR171'
        
        if (Me%FaceWaterColumn /= WCMaxBottom_ .and. Me%FaceWaterColumn /= WCAverageBottom_) then
            write(*,*) 'Unknown option for WATER_COLUMN_FACE'
            stop 'ReadDataFile - ModuleRunOff - ERR172'
        endif        
        
        if (Me%FaceWaterColumn == WCMaxBottom_) then
            !Gets if compute "margins" aside of adjacent cells that produce friction
            call GetData(Me%CalculateCellMargins,                               &
                         Me%ObjEnterData, iflag,                                   &
                         SearchType   = FromFile,                               &
                         keyword      = 'HYDRAULIC_RADIUS_MARGINS',             &
                         default      = .true.,                                 &
                         ClientModule = 'ModuleRunOff',                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR180'        
        endif
        
        !Gets if solution is limited by an maximum velocity
        call GetData(Me%ImposeMaxVelocity,                                      &
                     Me%ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'IMPOSE_MAX_VELOCITY',                      &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR190'

        if (Me%ImposeMaxVelocity) then
        
            !Gets if solution is limited by an maximum velocity
            call GetData(Me%ImposedMaxVelocity,                                     &
                         Me%ObjEnterData, iflag,                                       &
                         SearchType   = FromFile,                                   &
                         keyword      = 'MAX_VELOCITY',                             &
                         default      = 0.1,                                        &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR200'
        
        endif


        !Gets if Manning Coeficient is increased with water depth
        call GetData(DynamicAdjustManning,                                      &
                     Me%ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'DYNAMIC_ADJUST_MANNING',                   &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR210'

        if (iflag > 0 .and. .not. DynamicAdjustManning) then
            write(*,*)'The option DynamicAdjustManning (DYNAMIC_ADJUST_MANNING) has been removed.'
            write(*,*)'Please review your runoff data file!'
        endif

        if (DynamicAdjustManning) then
            write(*,*)'The option DynamicAdjustManning (DYNAMIC_ADJUST_MANNING) has been removed.'
            write(*,*)'Please review your runoff data file!'
            stop
        endif

        !Minimum Water Column for overland flow
        call GetData(Me%MinimumWaterColumn,                                              &
                     Me%ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MIN_WATER_COLUMN',                                  &
!                     default      = 0.001,                                              &
                     ClientModule = 'ModuleRunOff',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR220'
        if (iflag == 0) then
            write(*,*)'MIN_WATER_COLUMN must be defined in module Runoff'
            stop 'ReadDataFile - ModuleRunOff - ERR230'
        endif

        !Continuous Computation
        call GetData(Me%Continuous,                                                      &
                     Me%ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONTINUOUS',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleRunoff',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR240'

        if (Me%Continuous) then
            call ReadFileName('RUNOFF_INI', Me%Files%InitialFile,                         &
                               Message = "Runoff Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0250'
        endif

        call GetData(Me%StopOnWrongDate,                                    &
                     Me%ObjEnterData, iflag,                                   &
                     SearchType   = FromFile,                               &
                     keyword      = 'STOP_ON_WRONG_DATE',                   &
                     default      = .true.,                                 &
                     ClientModule = 'ModuleBasin',                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR260'
!
!        !Factor for DT Prediction
!        call GetData(Me%DTFactor,                                           &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'DT_FACTOR',                            &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = 1.05,                                   &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR270'        
!
!        if (Me%DTFactor <= 1.0) then
!            write (*,*)'Invalid DT Factor [DT_FACTOR]'
!            write (*,*)'Value must be greater then 1.0'
!            stop 'ReadDataFile - ModuleRunOff - ERR275'              
!        endif
!        
!        call GetData(Me%CV%DTFactorUp,                                         &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'DT_FACTOR_UP',                         &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = Me%DTFactor,                            &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR280'
!        if (iflag /= 1) then
!            write(*,*) 'Assumed a value of ', Me%DTFactor, ' for RunOff DT_FACTOR_UP' 
!        endif
!        if (Me%DTFactorUp <= 1.0) then
!            write (*,*)'Invalid DT Factor Up [DT_FACTOR_UP]'
!            write (*,*)'Value must be greater then 1.0'
!            stop 'ReadDataFile - ModuleRunOff - ERR281'              
!        endif
!                
!        call GetData(Me%CV%DTFactorDown,                                       &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'DT_FACTOR_DOWN',                       &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = Me%DTFactor,                            &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR285'
!        if (iflag /= 1) then
!            write(*,*) 'Assumed a value of ', Me%DTFactor, ' for RunOff DT_FACTOR_DOWN'
!        endif
!        if (Me%DTFactorDown <= 1.0) then
!            write (*,*)'Invalid DT Factor Down [DT_FACTOR_DOWN]'
!            write (*,*)'Value must be greater then 1.0'
!            stop 'ReadDataFile - ModuleRunOff - ERR286'              
!        endif      
!
!        call GetData(Me%CV%MaxIterations,                                      &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'MAX_ITERATIONS',                       &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = 5,                                      &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR310'

!        !Stabilize Solution
!        call GetData(Me%CV%Stabilize,                                          &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'STABILIZE',                            &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = .true.,                                 &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR290'        
        
!        if (Me%Stabilize) then
!            call GetData(Me%CV%StabilizeFactor,                                &
!                         Me%ObjEnterData, iflag,                               &  
!                         keyword      = 'STABILIZE_FACTOR',                 &
!                         ClientModule = 'ModuleRunOff',                     &
!                         SearchType   = FromFile,                           &
!                         Default      = 0.1,                                &
!                         STAT         = STAT_CALL)                                  
!            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR300' 
!            
!            !Minimum Water Column for checking stabilize
!            call GetData(Me%MinimumWaterColumnStabilize,                    &
!                         Me%ObjEnterData, iflag,                               &
!                         SearchType   = FromFile,                           &
!                         keyword      = 'MIN_WATER_COLUMN_STABILIZE',       &
!                         default      = Me%MinimumWaterColumn,              &
!                         ClientModule = 'ModuleRunOff',                     &
!                         STAT         = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR301'
!            
!            !Minimum Water Column for checking stabilize
!            call GetData(dummy,                                             &
!                         Me%ObjEnterData, iflag,                               &
!                         SearchType   = FromFile,                           &
!                         keyword      = 'STABILIZE_HARD_CUT_FACTOR',        &
!                         default      = 0.1,                                &
!                         ClientModule = 'ModuleRunOff',                     &
!                         STAT         = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR302' 
!            if (dummy <= 0.0 .or. dummy > 1.0) stop 'ReadDataFile - ModuleRunOff - ERR303'                         
!            Me%CV%StabilizeHardCutLimit = dummy * Me%MaxIterations
!        endif
!
!        !Internal Time Step Split
!        call GetData(Me%InternalTimeStepSplit,                              &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'DT_SPLIT_FACTOR',                      &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = 1.25,                                   &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR311'        
!        if (Me%InternalTimeStepSplit <= 1.0) then
!            write (*,*)'Invalid DT Factor [DT_SPLIT_FACTOR]'
!            write (*,*)'Value must be greater then 1.0'
!            stop 'ReadDataFile - ModuleRunOff - ERR312'              
!        endif        
!        
!        call GetData(Me%MinIterations,                                      &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'MIN_ITERATIONS',                       &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = 1,                                      &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR313'     
!        if (Me%MinIterations < 1) then
!            write (*,*) 'MIN_ITERATIONS must be greater or equal to 1'
!            stop 'ReadDataFile - ModuleRunOff - ERR313a'
!        endif
!        
!       call GetData(dummy,                                                  &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'PERCENT_TO_RESTART',                   &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = 0.,                                     &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR313b'
!        if (dummy <= 0.) then
!            Me%MinToRestart = 0
!        else
!            call CountBasinPoints(dummy)
!        endif          
!        
!        call GetData(Me%CheckDecreaseOnly,                                  &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'CHECK_DEC_ONLY',                       &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = .false.,                                &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR314'                 
!        
!        !Gets flag of DT is limited by the courant number
!        call GetData(Me%CV%LimitDTCourant,                                  &
!                     Me%ObjEnterData, iflag,                                &  
!                     keyword      = 'LIMIT_DT_COURANT',                     &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = .false.,                                &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR320'        
!
!        if (Me%LimitDTCourant) then
!
!            !Gets Maximum allowed Courant Number
!            call GetData(Me%MaxCourant,                                         &
!                         Me%ObjEnterData, iflag,                                   &  
!                         keyword      = 'MAX_COURANT',                          &
!                         ClientModule = 'ModuleRunOff',                         &
!                         SearchType   = FromFile,                               &
!                         Default      = 1.0,                                    &
!                         STAT         = STAT_CALL)                                  
!            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR330'        
!
!        endif
!
!        !Gets flag of DT is limited by the volume variation
!        call GetData(Me%LimitDTVariation,                                   &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'LIMIT_DT_VARIATION',                   &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = .true.,                                 &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR340'        
!        
        
        !Impose Boundary Value
        call GetData(Me%ImposeBoundaryValue,                                    &
                     Me%ObjEnterData, iflag,                                       &  
                     keyword      = 'IMPOSE_BOUNDARY_VALUE',                    &
                     ClientModule = 'ModuleRunOff',                             &
                     SearchType   = FromFile,                                   &
                     Default      = .false.,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR350'        
        
        if (Me%ImposeBoundaryValue) then
            call GetData(Me%BoundaryValue,                                      &
                         Me%ObjEnterData, iflag,                                   &  
                         keyword      = 'BOUNDARY_VALUE',                       &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR360'        

            if (iflag == 0) then
                write(*,*)'BOUNDARY_VALUE must be defined in module Runoff'
                stop 'ReadDataFile - ModuleRunOff - ERR0361'
            endif


            call GetData(Me%MaxDtmForBoundary,                                  &
                         Me%ObjEnterData, iflag,                                   &  
                         keyword      = 'MAX_DTM_FOR_BOUNDARY',                 &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR362'        

            if (iflag == 0) then
                write(*,*)'MAX_DTM_FOR_BOUNDARY must be defined in module Runoff'
                stop 'ReadDataFile - ModuleRunOff - ERR0363'
            endif

            call GetData(Me%BoundaryMethod,                                     &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'BOUNDARY_METHOD',                      &
                         Default      = ComputeFlow_,                           &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR360'        

            if (Me%BoundaryMethod /= ComputeFlow_ .and. Me%BoundaryMethod /= InstantaneousFlow_) then
                write(*,*)'BOUNDARY_METHOD must be or 1 - Compute Flow or 2 - Instantaneous FlowOut'
                stop 'ReadDataFile - ModuleRunOff - ERR0363.5'
            endif
            
        endif
        
          ! this keywords were removed because two methods create instabilities and
        ! negative volumes (instanataneous mixing method 1 and 3) and because 
        ! method 4 is just the same as 2 just the maxflow is different (based on 
        ! instant mixing) that is not justified (instantaneous = instabilities)
        ! This is left just for debugging      
!        call GetData(Me%OverlandChannelInteractionMethod,                   &
!                     Me%ObjEnterData, iflag,                                   &  
!                     keyword      = 'OVERLAND_CHANNEL_INTERACTION_METHOD',  &
!                     ClientModule = 'ModuleRunOff',                         &
!                     SearchType   = FromFile,                               &
!                     Default      = 1,                                      &
!                     STAT         = STAT_CALL)                                  
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR364'          
!        if (iflag == 0) then               
!            call GetData(Me%OverlandChannelInteractionMethod,                   &
!                         Me%ObjEnterData, iflag,                                   &  
!                         keyword      = 'CHANNEL_LINK_METHOD',                  &
!                         ClientModule = 'ModuleRunOff',                         &
!                         SearchType   = FromFile,                               &
!                         Default      = 1,                                      &
!                         STAT         = STAT_CALL)                                  
!            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR365'         
!        else
!            write (*,*) 'OVERLAND_CHANNEL_INTERACTION_METHOD keyword is deprecated.'
!            write (*,*) 'Use CHANNEL_LINK_METHOD instead.'
!            stop 'ReadDataFile - ModuleRunOff - ERR366'
!        endif        
        
        !Discharges
        call GetData(Me%Discharges,                                         &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'DISCHARGES',                           &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR370'        

        !Discharges
        call GetData(Me%SimpleChannelInteraction,                           &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'SIMPLE_CHANNEL_FLOW',                  &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR380'        



        !Routes D4 Points
        call GetData(Me%RouteDFourPoints,                                   &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'ROUTE_D4',                             &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR390'
        
        if (Me%RouteDFourPoints) then
            !Routes D4 Points
            call GetData(Me%RouteDFourPointsOnDN,                               &
                         Me%ObjEnterData, iflag,                                   &  
                         keyword      = 'ROUTE_D4_ON_DN',                       &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = .false.,                                &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR391'

            call GetData(Me%RouteDFourMethod,                                   &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'ROUTE_D4_METHOD',                      &
 !                        Default      = Celerity_,                              &
                         Default      = Manning_,                               &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR392'        

            if (Me%RouteDFourMethod /= Celerity_ .and. Me%RouteDFourMethod /= Manning_) then
                write(*,*)'ROUTE_D4_METHOD must be or 1 - Celerity based or 2 - Manning Equation'
                stop 'ReadDataFile - ModuleRunOff - ERR0393'
            endif
            
        endif

        !Limits Flow to critical
        call GetData(Me%LimitToCriticalFlow,                                &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'LIMIT_TO_CRITICAL_FLOW',               &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .true.,                                 &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR394'

        

        !Storm Water Drainage
        call GetData(Me%StormWaterDrainage,                                 &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'STORM_WATER',                          &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR400'

        if (Me%StormWaterDrainage) then
        
            !Storm Water Infiltration Velocity
            call GetData(Me%StormWaterInfiltrationVelocity,                     &
                         Me%ObjEnterData, iflag,                                   &  
                         keyword      = 'STORM_WATER_INF_VEL',                  &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         default      = 1.4e-5,                                 & !~50mm/h
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR410'

            !Storm Water Transfer Coeficient
            call GetData(Me%StormWaterFlowVelocity,                             &
                         Me%ObjEnterData, iflag,                                   &  
                         keyword      = 'STORM_WATER_FLOW_VEL',                 &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         default      = 0.2,                                    & !0.2m/s
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR420'

        endif

        !If Buildings are to be simulated (flow ocuation in urban areas)
        call GetData(Me%Buildings,                                          &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'BUILDINGS',                            &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR430'
        
        !If Connected to a StormWater model
        call GetData(Me%StormWaterModel,                                    &
                     Me%ObjEnterData, iflag,                                   &  
                     keyword      = 'STORM_WATER_MODEL_LINK',               &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR440'
        
                   

        !Gets Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                                 &
                           CurrentTime = Me%ExtVar%Now,                                  &
                           EndTime     = Me%EndTime,                                     &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%Yes,                                  &
                           STAT        = STAT_CALL)
        Me%OutPut%NextOutPut = 1    


        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                                &
                           CurrentTime  = Me%ExtVar%Now,                                &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoff - ERR450'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                      &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleBasin',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR460'



        call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0470'

        !Gets Block for OverLand Coef
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                 &
                                    '<BeginOverLandCoefficient>',               &
                                    '<EndOverLandCoefficient>', BlockFound,     &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0480'
        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = Me%OverLandCoefficientID,     &
                                        EnterDataID      = Me%ObjEnterData,                 &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%OverLandCoefficient,       &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0490'

            call KillFillMatrix(Me%OverLandCoefficientID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0500'


            !Check that manning values entered are not zero or negative
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    if (.not. Me%OverLandCoefficient(i,j) .gt. 0.0) then
                        write(*,*) 'Found Manning Overland coefficient zero or negative in input'
                        write(*,*) 'in cell', i, j
                        stop 'ReadDataFile - Module Runoff - ERR0510'
                    endif
                
                
                endif
                
            enddo
            enddo

        else
            write(*,*)'Missing Block <BeginOverLandCoefficient> / <EndOverLandCoefficient>' 
            stop      'ReadDataFile - ModuleRunOff - ERR0520'
        endif
        

        
        !Gets Block for OverLand Coef Difference 
        !To compute overland resistance in bottom for shear computation (erosion/deposition).
        !This process was created to remove from manning the resistance given by 
        !aerial vegetation parts that affect flow but do not affect bottom shear. Without that,
        !a manning increase (e.g. forestation) in one cell increases water depth (and reduces velocity)
        !but may increase shear stress (because water height increase is transformed in bottom resistance 
        !using manning - chezy see module runoff properties)
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                    '<BeginOverLandCoefficientDelta>',           &
                                    '<EndOverLandCoefficientDelta>', BlockFound, &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0530'
        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = OverLandCoefficientDeltaID,   &
                                        EnterDataID      = Me%ObjEnterData,                 &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%OverLandCoefficientDelta,  &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0540'

            call KillFillMatrix(OverLandCoefficientDeltaID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0550'

            !Check that final manning values are not zero or negative
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    if (.not. (Me%OverLandCoefficient(i,j) - Me%OverLandCoefficientDelta(i,j)) .gt. 0.0) then
                        write(*,*) 'Manning Overland coefficient delta found zero or negative in input'
                        write(*,*) 'in cell', i, j
                        stop 'ReadDataFile - Module Runoff - ERR0560'
                    endif
                
                
                endif
                
            enddo
            enddo

        else
            !Do not remove aerial vegetation effect from manning 
            Me%OverLandCoefficientDelta(:,:) = 0.0
        endif
        
        
        !Looks for StormWater DrainageCoef
        if (Me%StormWaterDrainage) then

            allocate(Me%StormWaterDrainageCoef (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%StormWaterDrainageCoef  = null_real
            

            call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR570'

            !Gets Flag with Sewer Points
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                          &
                                        '<BeginStormWaterDrainage>',                         &
                                        '<EndStormWaterDrainage>', BlockFound,               &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR580'

            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = StormWaterDrainageID,      &
                                            EnterDataID      = Me%ObjEnterData,                 &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%StormWaterDrainageCoef,    &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR590'

                call KillFillMatrix(StormWaterDrainageID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR600'

            else
                write(*,*)'Missing Block <BeginStormWaterDrainage> / <EndStormWaterDrainage>' 
                stop      'ReadDataFile - ModuleRunOff - ERR0610'
            endif
            
        endif

        allocate(Me%BuildingsHeight(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%BuildingsHeight = 0.0

        if (Me%Buildings) then
        
            call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR620'

            !Gets Flag with Sewer Points
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                          &
                                        '<BeginBuildingsHeight>',                          &
                                        '<EndBuildingsHeight>', BlockFound,                &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR630'

            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = BuildingsHeightID,          &
                                            EnterDataID      = Me%ObjEnterData,                 &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%BuildingsHeight,         &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR640'

                call KillFillMatrix(BuildingsHeightID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR650'

            else
                write(*,*)'Missing Block <BeginBuildingsHeight> / <EndBuildingsHeight>' 
                stop      'ReadDataFile - ModuleRunOff - ERR0670'
            endif
        
        endif
        
        !Looks for StormWater Interaction Points
        if (Me%StormWaterModel) then
        
            allocate(Me%StormWaterInteraction(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StreetGutterLength   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR680'

            !Gets Flag with Sewer Points
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                          &
                                        '<BeginStormWaterInteraction>',                      &
                                        '<EndStormWaterInteraction>', BlockFound,            &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR690'

            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = StormWaterInteractionID,      &
                                            EnterDataID      = Me%ObjEnterData,                 &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%StormWaterInteraction,     &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR691'

                call KillFillMatrix(StormWaterInteractionID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR692'

            else
                write(*,*)'Missing Block <BeginStormWaterInteraction> / <EndStormWaterInteraction>' 
                stop      'ReadDataFile - ModuleRunOff - ERR693'
            endif
            
            !Gets Street Gutter Length in each grid cell
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                          &
                                        '<BeginStreetGutterLength>',                         &
                                        '<EndStreetGutterLength>', BlockFound,               &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR694'

            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = StreetGutterLengthID,         &
                                            EnterDataID      = Me%ObjEnterData,                 &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%StreetGutterLength,        &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR700'

                call KillFillMatrix(StreetGutterLengthID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR710'

            else
                write(*,*)'Missing Block <BeginStreetGutterLength> / <EndStreetGutterLength>' 
                stop      'ReadDataFile - ModuleRunOff - ERR711'
            endif
            
        endif
                
        

        !Write Max Flow Modulus File 
        call GetData(Me%WriteMaxFlowModulus,                                    &
                     Me%ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'WRITE_MAX_FLOW_FILE',                      &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR720'

        if(Me%WriteMaxFlowModulus) then
            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%MaxFlowModulusFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR730'
            Me%MaxFlowModulusFile = trim(adjustl(Me%MaxFlowModulusFile))//"MaxRunOff.dat"
        end if

        !Write Max Channels Level  
        call GetData(Me%WriteMaxWaterColumn,                                    &
                     Me%ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'WRITE_MAX_WATER_COLUMN',                   &
                     default      = .true.,                                     &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR740'

        if(Me%WriteMaxWaterColumn) then
            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%MaxWaterColumnFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0750'
            Me%MaxWaterColumnFile = trim(adjustl(Me%MaxWaterColumnFile))//"MaxWaterColumn.dat"
        end if

        call ReadConvergenceParameters

        !Closes Data File
        call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR760'


    end subroutine ReadDataFile

    !--------------------------------------------------------------------------
    
    subroutine ReadConvergenceParameters
    
        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL,               &
                                                       iflag,                   &
                                                       MIN_WATER_COLUMN_STABILIZE_flag    
                                                            
        real                                        :: dummy_real
        
        !----------------------------------------------------------------------    
        
        !----------------------------------------------------------------------
        !Find deprecated keywords in data file
        !----------------------------------------------------------------------
        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, MIN_WATER_COLUMN_STABILIZE_flag,          &
                     SearchType     = FromFile,                                 &
                     keyword        ='MIN_WATER_COLUMN_STABILIZE',              &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR010")

        if (MIN_WATER_COLUMN_STABILIZE_flag > 0) then
            
            write (*,*) '======================================================================='
            write (*,*) 'The following deprecated keywords were found in RunOff data file:'
            write (*,*) ''
            
            if (MIN_WATER_COLUMN_STABILIZE_flag > 0) &
                write(*,*) 'MIN_WATER_COLUMN_STABILIZE: Use STABILIZE_MIN_WATER_COLUMN instead.'
                
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR070")                              
        endif

        !----------------------------------------------------------------------
        !Read convergence options
        !----------------------------------------------------------------------  
        call GetData(Me%CV%Stabilize,                                           &
                     Me%ObjEnterData, iflag,                                    &  
                     keyword      = 'STABILIZE',                                &
                     ClientModule = 'ModuleRunOff',                             &
                     SearchType   = FromFile,                                   &
                     Default      = .false.,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) & 
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR080")        
        if (iflag <= 0) then 
            write(*,*) 'WARNING: Missing STABILIZE keyword in RunOff input data file.'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR081")        
            
        endif
        if (Me%CV%Stabilize) then                
            !Maximun change of water content (in %) allowed in one time step.
            call GetData(Me%CV%StabilizeFactor,                                     &
                         Me%ObjEnterData, iflag,                                    &  
                         keyword      = 'STABILIZE_FACTOR',                         &
                         ClientModule = 'ModuleRunOff',                             &
                         SearchType   = FromFile,                                   &
                         Default      = 0.1,                                        &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) & 
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR082")

            if (Me%CV%StabilizeFactor < 0.0 .or. Me%CV%StabilizeFactor > 1.0) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR083")
                
            call GetData(Me%CV%MinimumValueToStabilize,                     &
                         Me%ObjEnterData, iflag,                            &
                         SearchType   = FromFile,                           &
                         keyword      = 'STABILIZE_MIN_WATER_COLUMN',       &
                         default      = Me%MinimumWaterColumn,              &
                         ClientModule = 'ModuleRunOff',                     &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR084")
            if (Me%CV%MinimumValueToStabilize < Me%MinimumWaterColumn) then
                write (*,*)'Invalid Minimun Water Column to Stabilize value [STABILIZE_MIN_WATER_COLUMN]'
                write (*,*)'Value must be greater than MIN_WATER_COLUMN'            
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR085")
            endif      
            
            call GetData(dummy_real,                                            &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'STABILIZE_RESTART_FACTOR',             &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = 0.,                                     &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR086")
            if (dummy_real <= 0.) then
                Me%CV%MinToRestart = 0
            else
                call CountDomainPoints(dummy_real)
            endif  
            
            call GetData(Me%CV%CheckDecreaseOnly,                               &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'CHECK_DEC_ONLY',                       &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = .false.,                                &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR087")
        endif        

       !Number of iterations threshold for starting to ask for a lower DT 
        call GetData(Me%CV%MinIterations,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MIN_ITERATIONS',                          &
                     Default        = 1,                                        &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR090")
        if (Me%CV%MinIterations < 1) then
            write (*,*)'Invalid Minimun Iterations value [MIN_ITERATIONS]'
            write (*,*)'Value must be greater than 0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR091")
        endif                                 

        !Number of iterations threshold that causes the model to stop
        call GetData(Me%CV%MaxIterations,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MAX_ITERATIONS',                          &
                     Default        = 1024,                                     &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR100")
        if (Me%CV%MaxIterations < Me%CV%MinIterations) then
            write (*,*)'Invalid Maximun Iterations value [MAX_ITERATIONS]'
            write (*,*)'Value must be greater than the value of MIN_ITERATIONS'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR101")              
        endif
                            
        !% of the maximun iterations that causes the DT to be cut to the value of one internal time step
        call GetData(dummy_real,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'DT_CUT_FACTOR',                    &
                     default      = 0.1,                                &
                     ClientModule = 'ModuleRunOff',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR110") 
        if (dummy_real <= 0.0 .or. dummy_real > 1.0) then
            write (*,*)'Invalid DT Cut Factor [DT_CUT_FACTOR]'
            write (*,*)'Value must be >= 0.0 and < 1.0'        
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR111") 
        endif
        Me%CV%StabilizeHardCutLimit = dummy_real * Me%CV%MaxIterations
        
       !Internal Time Step Split
        call GetData(Me%CV%DTSplitFactor,                                   &
                     Me%ObjEnterData, iflag,                                &
                     keyword      = 'DT_SPLIT_FACTOR',                      &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = 2.0,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadConvergenceParameters - ModuleRunOff - ERR120'        
        if (Me%CV%DTSplitFactor <= 1.0) then
            write (*,*)'Invalid DT Split Factor [DT_SPLIT_FACTOR]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR121")              
        endif            

        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR',                               &
                     Default        = 1.25,                                     &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR130")             
        if (dummy_real <= 1.0) then
            write (*,*)'Invalid DT Factor [DT_FACTOR]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR131")              
        endif            
        
        call GetData(Me%CV%DTFactorUp,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR_UP',                            &
                     Default        = dummy_real,                               &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR140")  
        if (Me%CV%DTFactorUp <= 1.0) then
            write (*,*)'Invalid DT Factor Up [DT_FACTOR_UP]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR141")              
        endif                  
                
        call GetData(Me%CV%DTFactorDown,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR_DOWN',                          &
                     Default        = dummy_real,                               &
                     ClientModule   ='ModuleRunOff',                            &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR150")  
        if (Me%CV%DTFactorDown <= 1.0) then
            write (*,*)'Invalid DT Factor Down [DT_FACTOR_DOWN]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR151")
        endif                                           
        
        call GetData(Me%CV%LimitDTCourant,                                  &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'LIMIT_DT_COURANT',                     &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &                     
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR180") 
        if (iflag <= 0) then
            write(*,*) 'WARNING: Missing LIMIT_DT_COURANT keyword in RunOff input data file.'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR181")
        endif                    
        if (Me%CV%LimitDTCourant) then
            !Gets Maximum allowed Courant Number
            call GetData(Me%CV%MaxCourant,                                      &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'MAX_COURANT',                          &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = 1.0,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModuleRunOff - ERR181")        
        endif
        
        !----------------------------------------------------------------------
    
    end subroutine ReadConvergenceParameters
    
    !--------------------------------------------------------------------------    

    subroutine CheckBoundaryCells
        
        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                      :: CHUNK, i, j, di, dj
        real                                         :: Sum
        !Begin-----------------------------------------------------------------

   
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,di,dj,Sum)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j)  == BasinPoint) then

                !Check if near a boundary point (no diagonal)
do3:            do dj = -1, 1
do4:            do di = -1, 1
                    Sum = dj + di
                    if ((Me%ExtVar%BasinPoints(i+di, j+dj) == 0) .and. (Sum .eq. -1 .or. Sum .eq. 1)) then
                        Me%BoundaryCells(i,j) = BasinPoint
                        exit do3 
                    endif
                enddo do4
                enddo do3
                
            endif
        enddo do2
        enddo do1
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
    end subroutine CheckBoundaryCells
    
    !--------------------------------------------------------------------------   
    
    subroutine CountDomainPoints (percent)
    
        !Arguments-------------------------------------------------------------
        real                                        :: percent
        
        !Local----------------------------------------------------------------- 
        integer                                     :: i, j
        integer                                     :: count
        
        !Begin-----------------------------------------------------------------       
                
        count = 0
        
        !Initializes Water Column
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                count = count + 1
            endif

        enddo
        enddo
        
        Me%CV%MinToRestart = max(int(count * percent), 0)
    
    end subroutine CountDomainPoints
    
    !-------------------------------------------------------------------------

    subroutine InitializeVariables

        !Arguments-------------------------------------------------------------
        
        !Local----------------------------------------------------------------- 
        integer                                     :: i, j
        integer                                     :: di, dj
        real                                        :: lowestValue
        
        !Begin-----------------------------------------------------------------       
        
        !Initializes Water Column
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                Me%myWaterLevel(i, j)       = Me%InitialWaterColumn(i,j) + Me%ExtVar%Topography(i, j)
                Me%MyWaterColumn(i, j)      = Me%InitialWaterColumn(i,j)
                Me%MyWaterColumnOld(i, j)   = Me%MyWaterColumn(i,j)
                Me%StabilityPoints(i, j)    = 1
            endif

        enddo
        enddo
        
               
        if (Me%Buildings) then

            !Checks Building Height
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    if (Me%BuildingsHeight(i, j) .lt. 0.0) then
                        write(*,*)'Buildings Height must be greater then 0', i, j
                        stop 'InitializeVariables - ModuleRunOff - ERR01'                        
                    endif
                endif
            enddo
            enddo
           
        endif
        
        
        !Finds lowest neighbor for from D8
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                !Finds lowest neighbour
                lowestValue = Me%ExtVar%Topography(i, j)
                do dj = -1, 1
                do di = -1, 1
                    
                    if (dj /= 0 .and. di /= 0 .and. Me%ExtVar%BasinPoints(i+di, j+dj) == BasinPoint) then
                    
                        !Checks lowest neighbor
                        if (Me%ExtVar%Topography(i + di, j + dj) < lowestValue) then
                            
                            lowestValue = Me%ExtVar%Topography(i + di, j + dj)
                            Me%LowestNeighborI(i, j) = i + di
                            Me%LowestNeighborJ(i, j) = j + dj

                        endif
                        
                    endif
                    
                enddo
                enddo        

            endif
        
        enddo
        enddo

        !If drainage network module is associated and simple interaction, then don't apply stability
        !to river points
        if (Me%ObjDrainageNetwork /= 0 .and. Me%SimpleChannelInteraction) then
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%ExtVar%RiverPoints (i, j) == BasinPoint) then
                    Me%StabilityPoints(i, j)    =  0
                endif
                
            enddo
            enddo
        endif
        
        if (Me%RouteDFourPoints) then
            !Checks if a given point is a DFourSink Point -> No point in the four direction is lower then the current point
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                    if (((Me%ExtVar%BasinPoints(i+1, j) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i+1, j) >= Me%ExtVar%Topography(i, j)) .or.  &  
                          Me%ExtVar%BasinPoints(i+1, j) /= BasinPoint)                .and. &
                        ((Me%ExtVar%BasinPoints(i-1, j) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i-1, j) >= Me%ExtVar%Topography(i, j)) .or.  &
                          Me%ExtVar%BasinPoints(i-1, j) /= BasinPoint)                .and. &
                        ((Me%ExtVar%BasinPoints(i, j+1) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i, j+1) >= Me%ExtVar%Topography(i, j)) .or.  &
                          Me%ExtVar%BasinPoints(i, j+1) /= BasinPoint)                .and. &
                        ((Me%ExtVar%BasinPoints(i, j-1) == BasinPoint                 .and. &
                          Me%ExtVar%Topography (i, j-1) >= Me%ExtVar%Topography(i, j)) .or.  &
                          Me%ExtVar%BasinPoints(i, j-1) /= BasinPoint)) then
                        
                        if (Me%LowestNeighborI(i, j) /= i .or. Me%LowestNeighborJ(i, j) /= j) then

                            Me%DFourSinkPoint(i, j) = BasinPoint
                            
                            !D 4 Sink Points are not points where stability criteria is verified
                            Me%StabilityPoints(i, j)= 0
                        
                        endif
                        
                    endif

                endif
            
            enddo
            enddo
        
            !If drainage network modules is associated, then don't apply D4 on drainage network point
            if (Me%ObjDrainageNetwork /= 0) then
            
                if (.not. Me%RouteDFourPointsOnDN) then
                    do j = Me%Size%JLB, Me%Size%JUB
                    do i = Me%Size%ILB, Me%Size%IUB
                    
                        !Source Point is a DNet Point
                        if (Me%ExtVar%RiverPoints (i, j) == BasinPoint) then
                            Me%DFourSinkPoint(i, j) = 0
                        endif
                    
                    enddo
                    enddo
                endif       
            endif
        
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%DFourSinkPoint(i, j) == BasinPoint) then
                
                    if (Me%LowestNeighborI(i, j) /= null_int) then
                
                        !Neighbors of D 4 Sink Points are not points where stability criteria is verified
                        Me%StabilityPoints(Me%LowestNeighborI(i, j), Me%LowestNeighborJ(i, j)) = 0
                
                    endif

                endif

            enddo
            enddo        
        
        endif
                
    end subroutine InitializeVariables

    !--------------------------------------------------------------------------

    subroutine CheckRiverNetWorkConsistency

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real   , dimension(:, :), pointer           :: ChannelsNodeLength 


        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR01'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) then
                
                    if (ChannelsNodeLength(i, j) < 0.0) then
                        write(*,*)'Inconsistent River Network', i, j
                        stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR02'
                    endif
                
                else
                
                    if (ChannelsNodeLength(i, j) > 0.0) then
                        write(*,*)'Inconsistent River Network', i, j
                        stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR03'
                    endif
                
                endif

            endif

        enddo
        enddo
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRiverNetWorkConsistency - ModuleRunOff - ERR04'
    
    
    end subroutine CheckRiverNetWorkConsistency

    !--------------------------------------------------------------------------
    
    subroutine ConstructDischarges

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        character(len=StringLength)                 :: DischargeName
        real                                        :: CoordinateX, CoordinateY
        logical                                     :: CoordinatesON, IgnoreOK
        integer                                     :: Id, Jd, dn, DischargesNumber
        integer                                     :: STAT_CALL

        call Construct_Discharges(Me%ObjDischarges, Me%ObjTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR01' 
                
        call GetDischargesNumber(Me%ObjDischarges, DischargesNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR02' 

        do dn = 1, DischargesNumber

            call GetDischargesGridLocalization(Me%ObjDischarges, dn,            &
                                               CoordinateX   = CoordinateX,     &
                                               CoordinateY   = CoordinateY,     & 
                                               CoordinatesON = CoordinatesON,   &
                                               STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR03' 
                    
            call GetDischargesIDName (Me%ObjDischarges, dn, DischargeName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR03' 

            if (CoordinatesON) then
                
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordinateX, CoordinateY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR04' 

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreDischarge(Me%ObjDischarges, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR05' 

                    if (IgnoreOK) then
                        write(*,*) 'Discharge outside the domain - ',trim(DischargeName),' - ',trim(Me%ModelName)
                        cycle
                    else
                        stop 'ModuleRunOff - ConstructDischarges - ERR06' 
                    endif

                endif

                call CorrectsCellsDischarges(Me%ObjDischarges, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR07' 
                    
            endif

            !
            !FROM HERE TO THE END NOT TESTED CODE
            !

!
!            call GetDischargeSpatialEmission(Me%ObjDischarges, dn, LineX, PolygonX, &
!                                             SpatialEmission, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR08' 
!                    
!            if (SpatialEmission == DischPoint_) then
! 
!                call GetDischargesGridLocalization(Me%ObjDischarges, dn,            &
!                                                   Igrid         = Id,              &
!                                                   JGrid         = Jd,              &
!                                                   STAT          = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR09' 
!
!                if (Me%ExtVar%BasinPoints(Id,Jd) /= WaterPoint) then
!                    call TryIgnoreDischarge(Me%ObjDischarges, dn, IgnoreOK, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR10' 
!
!                    write(*,*) 'Discharge outside the domain I=',Id,' J=',Jd,'Model name=',trim(Me%ModelName)
!
!                    if (IgnoreOK) then
!                        write(*,*) 'Discharge in a land cell - ',trim(DischargeName),' - ',trim(Me%ModelName)
!                        cycle
!                    else
!                        stop 'ModuleRunOff - ConstructDischarges - ERR11' 
!                    endif
!                endif
!
!                nCells    = 1
!                allocate(VectorI(nCells), VectorJ(nCells))
!                VectorJ(nCells) = Jd
!                VectorI(nCells) = Id
!
!            else
!
!                if (SpatialEmission == DischLine_) then
!                    call GetCellZInterceptByLine(Me%ObjHorizontalGrid, LineX,       &
!                                                 Me%ExtVar%BasinPoints, VectorI, VectorJ,   &
!                                                 nCells, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR12' 
!
!                    if (nCells < 1) then
!                        write(*,*) 'Discharge line intercept 0 cells'       
!                        stop 'ModuleRunOff - ConstructDischarges - ERR13' 
!                    endif
!
!                endif 
!
!
!                if (SpatialEmission == DischPolygon_) then
!                    call GetCellZInterceptByPolygon(Me%ObjHorizontalGrid, PolygonX, &
!                                                 Me%ExtVar%BasinPoints, VectorI, VectorJ,   &
!                                                 nCells, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR14' 
!
!                    if (nCells < 1) then
!                        write(*,*) 'Discharge contains 0 center cells'       
!                        write(*,*) 'Or the polygon is to small and is best to a discharge in a point or'
!                        write(*,*) 'the polygon not define properly'
!                        stop 'ModuleRunOff - ConstructDischarges - ERR15' 
!                    endif
!
!                endif
!
!
!            endif
                        
!            if (SpatialEmission /= DischPoint_) then
!
!
!                call SetLocationCellsZ (Me%ObjDischarges, dn, nCells, VectorI, VectorJ, VectorK, STAT= STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructDischarges - ERR16' 
!
!            else  i4
!                if (DischVertical == DischBottom_ .or. DischVertical == DischSurf_) then
!                    call SetLayer (Me%ObjDischarges, dn, VectorK(nCells), STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Sub_Modules - ModuleHydrodynamic - ERR220' 
!                endif
!                deallocate(VectorI, VectorJ, VectorK)
!            endif i4

        enddo

   
    end subroutine
    
    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------
        !Local----------------------------------------------------------------- 
        !Begin-----------------------------------------------------------------       

        allocate(Me%iFlowToChannels  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowToChannels  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%iFlowToChannels      = 0.0   !Sets values initially to zero, so 
        Me%lFlowToChannels      = 0.0   !model can run without DNet
        
        allocate(Me%lFlowBoundary    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%iFlowBoundary    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%lFlowBoundary        = 0.0   !Sets values initially to zero, so 
        Me%iFlowBoundary        = 0.0   !model can run without BC

        allocate(Me%lFlowDischarge    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%iFlowDischarge    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%lFlowDischarge        = 0.0   !Sets values initially to zero, so 
        Me%iFlowDischarge        = 0.0   !model can run without Dis

        allocate(Me%iFlowRouteDFour  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%iFlowRouteDFour       = 0.0   !Sets values initially to zero, so  
        
        allocate(Me%BoundaryCells     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%BoundaryCells = 0
        
        allocate(Me%myWaterLevel         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterColumn        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterColumnAfterTransport (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolumePred   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%myWaterVolumePred = null_real
        allocate(Me%InitialWaterColumn   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolume        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterColumnOld     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolumeOld     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%MassError            (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%myWaterLevel            = null_real
        Me%myWaterColumn           = null_real
        Me%InitialWaterColumn      = null_real
        Me%myWaterVolume           = 0.0        !For OpenMI
        Me%myWaterColumnOld        = null_real
        Me%myWaterVolumeOld        = null_real
        Me%MassError               = 0
        Me%myWaterColumnAfterTransport = null_real

        allocate(Me%iFlowX               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%iFlowY               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowX               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowY               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%FlowXOld             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%FlowYOld             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%InitialFlowX         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%InitialFlowY         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AreaU                (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AreaV                (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ComputeFaceU         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ComputeFaceV         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        
        
        
        Me%iFlowX               = 0.0
        Me%iFlowY               = 0.0
        Me%lFlowX               = 0.0
        Me%lFlowY               = 0.0
        Me%FlowXOld             = 0.0
        Me%FlowYOld             = 0.0
        Me%InitialFlowX         = 0.0
        Me%InitialFlowY         = 0.0
        Me%AreaU                = 0.0
        Me%AreaV                = 0.0
        Me%ComputeFaceU         = 0
        Me%ComputeFaceV         = 0

        allocate(Me%OverLandCoefficient  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%OverLandCoefficientDelta (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%OverLandCoefficientX (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%OverLandCoefficientY (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%OverLandCoefficient       = null_real
        Me%OverLandCoefficientDelta  = null_real
        Me%OverLandCoefficientX      = null_real
        Me%OverLandCoefficientY      = null_real
        
       
        allocate (Me%CenterFlowX    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%CenterFlowY    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%FlowModulus    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%CenterVelocityX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%CenterVelocityY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%VelocityModulus(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        allocate (Me%MaxFlowModulus (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%MaxWaterColumn (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))      

        Me%MaxFlowModulus = null_real
        Me%MaxWaterColumn = null_real

        allocate (Me%LowestNeighborI (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%LowestNeighborJ (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%DFourSinkPoint  (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate (Me%StabilityPoints (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                    
        Me%LowestNeighborI = null_int
        Me%LowestNeighborJ = null_int
        Me%DFourSinkPoint  = 0
        Me%StabilityPoints = 0


    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ConstructOverLandCoefficient

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB    
        integer                                             :: i, j

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !TODO: OpenMP - Missing implementation
        do j = JLB, JUB + 1
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j) + Me%ExtVar%BasinPoints(i, j-1) == 2) then !Two Basin Points
            
                Me%OverlandCoefficientX(i, j) = (Me%ExtVar%DUX(i, j  ) * Me%OverlandCoefficient(i, j-1  )  + &
                                                 Me%ExtVar%DUX(i, j-1) * Me%OverlandCoefficient(i, j)) / &
                                                 (Me%ExtVar%DUX(i, j) + Me%ExtVar%DUX(i, j-1))
            endif

        enddo
        enddo

        do j = JLB, JUB
        do i = ILB, IUB + 1

            if (Me%ExtVar%BasinPoints(i, j) + Me%ExtVar%BasinPoints(i-1, j) == 2) then !Two Basin Points
            
                Me%OverlandCoefficientY(i, j) =     (Me%ExtVar%DVY(i, j  ) * Me%OverlandCoefficient(i-1, j  )  + &
                                                     Me%ExtVar%DVY(i-1, j) * Me%OverlandCoefficient(i, j)) / &
                                                     (Me%ExtVar%DVY(i, j) + Me%ExtVar%DVY(i-1, j))
            endif

        enddo
        enddo


    end subroutine ConstructOverLandCoefficient

    !--------------------------------------------------------------------------

    subroutine ConstructStormWaterDrainage
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB    
        integer                                             :: i, j
        logical                                             :: nearestfound
        integer                                             :: dij, lowestI, lowestJ
        integer                                             :: iAux, jAux
        real                                                :: lowestValue
        logical                                             :: IgnoreTopography

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
                
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Checks consistency
        if (Me%StormWaterDrainage .and. Me%StormWaterModel) then
            write(*,*)'It is not possible to activate a simplifed Storm Water model and SWMM at the same time'
            stop 'ConstructStormWaterDrainage - ModuleRunOff - ERR01'
        endif
        
        !Simplified Storm Water Drainage
        if (Me%StormWaterDrainage) then
            allocate(Me%StormWaterVolume       (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StormWaterFlowX        (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StormWaterFlowY        (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StormWaterCenterFlowX  (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StormWaterCenterFlowY  (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StormWaterCenterModulus(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            Me%StormWaterVolume        = null_real
            Me%StormWaterFlowX         = 0.0
            Me%StormWaterFlowY         = 0.0
            Me%StormWaterCenterFlowX   = 0.0
            Me%StormWaterCenterFlowY   = 0.0
            Me%StormWaterCenterModulus = 0.0

            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then 
                    Me%StormWaterVolume(i, j)     = 0.0
                endif

            enddo
            enddo                
        endif            

        !Model link like SMWM
        if (Me%StormWaterModel) then
        
            allocate(Me%StormWaterModelFlow    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StreetGutterFlow       (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SewerInflow            (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StreetGutterTargetI    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StreetGutterTargetJ    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%StormInteractionFlow   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%StormWaterModelFlow    = 0.0
            Me%StreetGutterFlow       = 0.0
            Me%SewerInflow            = 0.0
            Me%StormInteractionFlow   = 0.0
            
            Me%StreetGutterTargetI    = null_int
            Me%StreetGutterTargetJ    = null_int
            
            !Algorithm to find the nearest sewer interaction point near the street gutter. 
            !Point must be lower equal current point
            !Algorithm is not quiet correct, since it does not search in circles, but in rectangles
            !Algorithm is not very eficent, since it should look only to the border points of the rectangls. But we are in the constructor
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint .and. Me%StreetGutterLength(i, j) > AllmostZero) then 
                    
                    nearestfound = .false.
                    dij = 0
                    IgnoreTopography = .false.
                    do while (.not. nearestfound)

                        lowestValue = Me%ExtVar%Topography(i, j)
                        lowestI     = null_int
                        lowestJ     = null_int

                        !Left
                        jAux = Max(j-dij, JLB)
                        do iAux = Max(i-dij, ILB), Min(i+dij, IUB)

                            if (Me%ExtVar%BasinPoints(iAux, jAux) == OpenPoint) then
                                if ((IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= Me%ExtVar%Topography(i, j)) .and. &
                                     Me%StormWaterInteraction(iAux, jAux) > AllmostZero) then
                                    nearestfound = .true.
                                    if (IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= lowestValue) then
                                        lowestValue = Me%ExtVar%Topography(iAux, jAux)
                                        lowestI     = iAux
                                        lowestJ     = jAux
                                    endif
                                endif
                            endif

                        enddo

                        !Right
                        jAux = Min(j+dij, JUB)
                        do iAux = Max(i-dij, ILB), Min(i+dij, IUB)

                            if (Me%ExtVar%BasinPoints(iAux, jAux) == OpenPoint) then
                                if ((IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= Me%ExtVar%Topography(i, j)) .and. &
                                    Me%StormWaterInteraction(iAux, jAux) > AllmostZero) then
                                    nearestfound = .true.
                                    if (IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= lowestValue) then
                                        lowestValue = Me%ExtVar%Topography(iAux, jAux)
                                        lowestI     = iAux
                                        lowestJ     = jAux
                                    endif
                                endif
                            endif

                        enddo
 
                        !Bottom
                        iAux = Max(i-dij, ILB)
                        do jAux = Max(j-dij, JLB), Min(j+dij, JUB)

                            if (Me%ExtVar%BasinPoints(iAux, jAux) == OpenPoint) then
                                if ((IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= Me%ExtVar%Topography(i, j)) .and. &
                                     Me%StormWaterInteraction(iAux, jAux) > AllmostZero) then
                                    nearestfound = .true.
                                    if (IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= lowestValue) then
                                        lowestValue = Me%ExtVar%Topography(iAux, jAux)
                                        lowestI     = iAux
                                        lowestJ     = jAux
                                    endif
                                endif
                            endif
                        enddo
 
                        !Top
                        iAux = Min(i+dij, IUB)
                        do jAux = Max(j-dij, JLB), Min(j+dij, JUB)

                            if (Me%ExtVar%BasinPoints(iAux, jAux) == OpenPoint) then
                                if ((IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= Me%ExtVar%Topography(i, j)) .and. &
                                    Me%StormWaterInteraction(iAux, jAux) > AllmostZero) then
                                    nearestfound = .true.
                                    if (IgnoreTopography .or. Me%ExtVar%Topography(iAux, jAux) <= lowestValue) then
                                        lowestValue = Me%ExtVar%Topography(iAux, jAux)
                                        lowestI     = iAux
                                        lowestJ     = jAux
                                    endif
                                endif
                            endif
                        enddo
 
!
!                        
!                        !Efficiency is somethink else.... we should only travel arround the outer cells
!                        do jAux = Max(j-dij, JLB), Min(j+dij, JUB)
!                        do iAux = Max(i-dij, ILB), Min(i+dij, IUB)
!
!                            if (Me%ExtVar%BasinPoints(iAux, jAux) == OpenPoint) then
!                                if (Me%ExtVar%Topography(iAux, jAux) <= Me%ExtVar%Topography(i, j) .and. &
!                                    Me%StormWaterInteraction(iAux, jAux) > AllmostZero) then
!                                    nearestfound = .true.
!                                    if (Me%ExtVar%Topography(iAux, jAux) <= lowestValue) then
!                                        lowestValue = Me%ExtVar%Topography(iAux, jAux)
!                                        lowestI     = iAux
!                                        lowestJ     = jAux
!                                    endif
!                                endif
!                            endif
!
!                        enddo
!                        enddo

                        dij = dij + 1
                        if (dij > IUB .and. dij > JUB) then
                            if (.not. IgnoreTopography) then
                                IgnoreTopography = .true.
                                dij = 0
                                write(*,*)'Topography Ignored for', i, j 
                            else
                                write(*,*)'Internal Error locating Street Gutter Target', i, j 
                            endif
                        endif


                    enddo
                    
                    !Sets Link
                    Me%StreetGutterTargetI(i, j)  = lowestI
                    Me%StreetGutterTargetJ(i, j)  = lowestJ
                    
                endif

            enddo
            enddo                
            
            
        endif            

    end subroutine

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Output

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR04'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR05'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",                   &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR06'

        !Writes the River Points
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "RiverPoints", "-",                   &
                              Array2D = Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR07'


        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR08'


    end subroutine ConstructHDF5Output

    !--------------------------------------------------------------------------

    subroutine ReadInitialFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: InitialFile
        type (T_Time)                               :: BeginTime, EndTimeFile, EndTime
        real                                        :: DT_error
        integer                                     :: STAT_CALL, i, j

        !----------------------------------------------------------------------

        call UnitsManager(InitialFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleRunoff - ERR01'

        open(Unit = InitialFile, File = Me%Files%InitialFile, Form = 'UNFORMATTED',     &
             status = 'OLD', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleRunoff - ERR02'

        !Reads Date
        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleRunoff - ERR03'
        
        DT_error = EndTimeFile - BeginTime

        !Avoid rounding erros - Frank 08-2001
        if (abs(DT_error) >= 0.01) then
            
            write(*,*) 'The end time of the previous run is different from the start time of this run'
            write(*,*) 'Date in the file'
            write(*,*) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
            write(*,*) 'DT_error', DT_error
            if (Me%StopOnWrongDate) stop 'ReadInitialFile - ModuleRunoff - ERR04'   

        endif

        read(InitialFile)Me%myWaterColumn
        
        if (Me%StormWaterDrainage) then
            read(InitialFile)Me%StormWaterVolume
        endif
        
   10   continue

        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleRunoff - ERR05'  
        
        !Updates Volume & Level from Column        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
            
                Me%myWaterLevel(i, j)   = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                Me%myWaterVolume(i, j)  = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                
            endif

        enddo
        enddo      
      
    end subroutine ReadInitialFile
 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
    subroutine GetOverLandFlow (ObjRunOffID, FlowX, FlowY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: FlowX, FlowY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowX => Me%iFlowX

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowY => Me%iFlowY

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetOverLandFlow

    !--------------------------------------------------------------------------

    subroutine GetFlowToChannels (ObjRunOffID, FlowToChannels, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: FlowToChannels
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowToChannels => Me%iFlowToChannels

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetFlowToChannels

    !--------------------------------------------------------------------------

    subroutine GetBoundaryImposed (ObjRunOffID, BoundaryOpen, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        logical, intent(OUT)                            :: BoundaryOpen
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryOpen = Me%ImposeBoundaryValue

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetBoundaryImposed
    
    !--------------------------------------------------------------------------

    subroutine GetRouteDFour (ObjRunOffID, RouteD4, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        logical, intent(OUT)                            :: RouteD4
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            RouteD4 = Me%RouteDFourPoints

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetRouteDFour

    !--------------------------------------------------------------------------   

    subroutine GetRouteDFourCells (ObjRunOffID, RouteD4Cells, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        integer, pointer, dimension (:,:)               :: RouteD4Cells
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            RouteD4Cells => Me%DFourSinkPoint

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetRouteDFourCells

    !--------------------------------------------------------------------------   

    subroutine GetRouteDFourNeighbours (ObjRunOffID, RouteD4LowerI, RouteD4LowerJ, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        integer, pointer, dimension (:,:)               :: RouteD4LowerI
        integer, pointer, dimension (:,:)               :: RouteD4LowerJ
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            RouteD4LowerI => Me%LowestNeighborI
            RouteD4LowerJ => Me%LowestNeighborJ

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetRouteDFourNeighbours

    !-------------------------------------------------------------------------- 

    subroutine GetRouteDFourFlux (ObjRunOffID, DFourFlow, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: DFourFlow
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            DFourFlow => Me%iFlowRouteDFour

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRouteDFourFlux
    
    !--------------------------------------------------------------------------

    subroutine GetBoundaryFlux (ObjRunOffID, FlowAtBoundary, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: FlowAtBoundary
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowAtBoundary => Me%iFlowBoundary

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryFlux
    
    !--------------------------------------------------------------------------

    subroutine GetBoundaryCells (ObjRunOffID, BoundaryCells, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        integer,   pointer, dimension(:,:)              :: BoundaryCells
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryCells => Me%BoundaryCells

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryCells

    !--------------------------------------------------------------------------

    subroutine GetFlowDischarge (ObjRunOffID, FlowDischarge, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: FlowDischarge
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            FlowDischarge => Me%iFlowDischarge

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetFlowDischarge    
    
    !--------------------------------------------------------------------------

        subroutine GetRunoffWaterLevel (ObjRunOffID, Waterlevel, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterLevel
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            WaterLevel => Me%MyWaterLevel

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterLevel

    !--------------------------------------------------------------------------
    
    subroutine GetRunoffWaterColumn (ObjRunOffID, WaterColumn, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumn
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            WaterColumn => Me%MyWaterColumn

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterColumn

    !--------------------------------------------------------------------------

    subroutine GetRunoffWaterColumnOld (ObjRunOffID, WaterColumnOld, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumnOld
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            WaterColumnOld => Me%MyWaterColumnOld

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterColumnOld

    !--------------------------------------------------------------------------
    
    subroutine GetRunoffWaterColumnAT (ObjRunOffID, WaterColumn, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumn
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            !because runoff water may go all to the river in one time step and Runoff Conc would be zero
            !but not flow, transport has to be separated from drainage network interaction
            !and explicit/implicit transport is only evaluated with water column after transport            
            WaterColumn => Me%MyWaterColumnAfterTransport

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffWaterColumnAT    

    !--------------------------------------------------------------------------

    subroutine GetRunoffCenterVelocity (ObjRunOffID, VelX, VelY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: VelX, VelY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            VelX => Me%CenterVelocityX

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            VelY => Me%CenterVelocityY

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffCenterVelocity

    !--------------------------------------------------------------------------
    
    subroutine GetManning (ObjRunOffID, Manning, ManningX, ManningY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer, optional        :: Manning, ManningX, ManningY
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(Manning)) then
                call Read_Lock(mRUNOFF_, Me%InstanceID)
                Manning => Me%OverlandCoefficient
            endif

            if (present(ManningX)) then
                call Read_Lock(mRUNOFF_, Me%InstanceID)
                ManningX => Me%OverlandCoefficientX
            endif

            if (present(ManningY)) then
                call Read_Lock(mRUNOFF_, Me%InstanceID)
                ManningY => Me%OverlandCoefficientY
            endif
            
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetManning

    !--------------------------------------------------------------------------    


    subroutine GetManningDelta (ObjRunOffID, ManningDelta, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: ManningDelta
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mRUNOFF_, Me%InstanceID)
            ManningDelta => Me%OverlandCoefficientDelta

           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetManningDelta

    !--------------------------------------------------------------------------   

    subroutine GetMassError (ObjRunOffID, MassError, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, dimension(:, :), pointer                  :: MassError
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjRunOffID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mRUNOFF_, Me%InstanceID)
            MassError => Me%MassError

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetMassError
    
    !--------------------------------------------------------------------------
    
    subroutine GetRunoffTotalStoredVolume (ObjRunoffID, TotalStoredVolume, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffID
        real(8)                                         :: TotalStoredVolume
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjRunoffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TotalStoredVolume = Me%TotalStoredVolume

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetRunoffTotalStoredVolume

    !-------------------------------------------------------------------------
    
    subroutine GetRunOffStoredVolumes (ID, Surface, StormSystem, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: ID
        real(8), intent(OUT), optional                  :: Surface
        real(8), intent(OUT), optional                  :: StormSystem
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        !Begin-----------------------------------------------------------------
        call Ready(ID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Surface))     Surface     = Me%VolumeStoredInSurface
            if (present(StormSystem)) StormSystem = Me%VolumeStoredInStormSystem

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_            
        !----------------------------------------------------------------------
    
    end subroutine GetRunOffStoredVolumes
    
    !-------------------------------------------------------------------------
    
    subroutine GetRunOffBoundaryFlowVolume (ID, BoundaryFlowVolume, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: ID
        real(8), intent(OUT)                            :: BoundaryFlowVolume
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        !Begin-----------------------------------------------------------------
        call Ready(ID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryFlowVolume = Me%BoundaryFlowVolume

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_            
        !----------------------------------------------------------------------    
    
    end subroutine GetRunOffBoundaryFlowVolume
    
    !-------------------------------------------------------------------------
        
    subroutine GetNextRunOffDT (ObjRunOffID, DT, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real, intent(OUT)                               :: DT
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ObjRunOffID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DT        = Me%CV%NextDT

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNextRunOffDT

    !--------------------------------------------------------------------------

    subroutine SetBasinColumnToRunoff(ObjRunOffID, WaterColumnOld, WaterColumn, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: WaterColumnOld
        real(8), dimension(:, :), pointer               :: WaterColumn
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL, ready_
        integer                                         :: i, j
        integer                                         :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(ObjRunOffID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then
        
            !Actualizes water column, water level and water volume
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR01'

            !Gets a pointer to Topography
            call GetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR10'
        
            call GetGridCellArea  (Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea,             &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR020'
        
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                    Me%myWaterColumnOld(i, j) = WaterColumnOld(i, j)
                    Me%myWaterColumn(i, j)    = WaterColumn(i, j)

                    Me%myWaterLevel (i, j) = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                    !Here the water column is the uniformly distributed one. Inside 
                    Me%myWaterVolume(i, j) = WaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)

!                    Me%myWaterVolume(i, j) = WaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
!                    Me%myWaterColumn(i, j) = Me%myWaterVolume(i, j) / Me%FreeGridCellArea(i, j)
!                    Me%myWaterLevel (i, j) = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                   
                endif
            enddo
            enddo            

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR030'

            !Ungets the Topography
            call UngetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR40'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinColumnToRunoff - ModuleRunOff - ERR050'

            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
                    
    end subroutine SetBasinColumnToRunoff

    !--------------------------------------------------------------------------

    subroutine UnGetRunOff2D_R4(ObjRunOffID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunOffID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRUNOFF_, Me%InstanceID, "UnGetRunOff2D_R4")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunOff2D_R4

    !--------------------------------------------------------------------------

    subroutine UnGetRunOff2D_R8(ObjRunOffID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunOffID
        real(8), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunOffID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRUNOFF_, Me%InstanceID, "UnGetRunOff2D_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunOff2D_R8
        
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine ModifyRunOff(RunOffID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL
        real                                        :: SumDT
        logical                                     :: Restart
        integer                                     :: Niter, iter
        integer                                     :: n_restart

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunOffID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ModifyRunOff")

            !Time Stuff
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunOff - ModuleRunOff - ERR010'

            !Stores initial values = from basin
            call SetMatrixValue(Me%myWaterColumnOld, Me%Size, Me%myWaterColumn)
            call SetMatrixValue(Me%InitialFlowX,     Me%Size, Me%iFlowX)
            call SetMatrixValue(Me%InitialFlowY,     Me%Size, Me%iFlowY)            

            !Adds Flow from SEWER OverFlow to Water Column OLD
            if (Me%StormWaterModel) then
                call ReadLockExternalVar   (StaticOnly = .true.)
                call AddFlowFromStormWaterModel
                call ReadUnLockExternalVar (StaticOnly = .true.)
            endif
            
            Restart = .true.
            n_restart = 0
            
            if (Me%CV%NextNiteration > 1 .and. Me%ExtVar%DT < (Me%CV%CurrentDT * Me%CV%NextNiteration)) then
                Me%CV%NextNiteration = max(aint(Me%ExtVar%DT / Me%CV%CurrentDT), 1.0)
            endif            
            
            do while (Restart)
            
                !Calculates local Watercolumn
                call ReadLockExternalVar   (StaticOnly = .true.)
                call LocalWaterColumn      (Me%myWaterColumnOld)
                call ReadUnLockExternalVar (StaticOnly = .true.)

                SumDT        = 0.0
                Restart      = .false.                                 
                iter         = 1
                Niter        = Me%CV%NextNiteration    !DB
                Me%CV%CurrentDT = Me%ExtVar%DT / Niter

                if (Niter > 1) then                
                    call WriteDTLog_ML ('ModuleRunOff', Niter, Me%CV%CurrentDT)
                endif
                
                call SetMatrixValue(Me%iFlowX, Me%Size, dble(0.0))
                call SetMatrixValue(Me%iFlowY, Me%Size, dble(0.0))
                call SetMatrixValue(Me%lFlowX, Me%Size, Me%InitialFlowX)
                call SetMatrixValue(Me%lFlowY, Me%Size, Me%InitialFlowY)
                call SetMatrixValue(Me%iFlowToChannels, Me%Size, 0.0)
                call SetMatrixValue(Me%iFlowBoundary, Me%Size, 0.0)
                call SetMatrixValue(Me%iFlowRouteDFour, Me%Size, 0.0)
                
doIter:         do while (iter <= Niter)

                    !Gets ExternalVars
                    call ReadLockExternalVar (StaticOnly = .false.)

                    !Stores WaterVolume for convergence test
                    call SetMatrixValue(Me%myWaterVolumeOld, Me%Size, Me%myWaterVolume)

                    call SetMatrixValue(Me%FlowXOld,         Me%Size, Me%lFlowX)
                    call SetMatrixValue(Me%FlowYOld,         Me%Size, Me%lFlowY)

                    !Inputs Water from discharges
                    if (Me%Discharges) then
                        call ModifyWaterDischarges  (Me%CV%CurrentDT)                
                    endif

                    !Updates Geometry
                    call ModifyGeometryAndMapping
                    
                    !save most recent water volume to predict if negative occur. in that case flux will be
                    !limited to water volume and next fluxes will be zero
                    if (.not. Me%LimitToCriticalFlow) then
                        call SetMatrixValue (Me%myWaterVolumePred, Me%Size, Me%myWaterVolume)
                    endif        
                    
                    select case (Me%HydrodynamicApproximation)
                        case (KinematicWave_)
                            call KinematicWave  ()            !Slope based on topography
                        case (DiffusionWave_)
                            call KinematicWave  ()            !Slope based on surface
                        case (DynamicWave_)
                            call DynamicWaveXX    (Me%CV%CurrentDT)   !Consider Advection, Friction and Pressure
                            call DynamicWaveYY    (Me%CV%CurrentDT)
                    end select

                    
                    !Updates waterlevels, based on fluxes
                    call UpdateWaterLevels(Me%CV%CurrentDT)
                    
                    !Interaction with channels
                    if (Me%ObjDrainageNetwork /= 0 .and. .not. Me%SimpleChannelInteraction) then
                        call FlowIntoChannels       (Me%CV%CurrentDT)
                    endif

                    !Boundary Condition
!                    if (Me%ImposeBoundaryValue) then
!                        call ImposeBoundaryValue    (Me%CV%CurrentDT)
!                    endif

                    
                    call CheckStability(Restart) 
                    
                    call ReadUnLockExternalVar (StaticOnly = .false.)
                    
                    if (Restart) then
                        exit doIter
                    endif

                    call IntegrateFlow     (Me%CV%CurrentDT, SumDT)  
                        
                    SumDT = SumDT + Me%CV%CurrentDT
                    iter  = iter  + 1                                        
                    
                enddo doIter
                                            
            enddo
            
!            !DB
!            if (Niter <= Me%LastGoodNiter) then
!                Me%CV%NextNiteration = max (min(int(Niter / Me%InternalTimeStepSplit), NIter - 1), 1)
!            else
!                Me%CV%NextNiteration = Niter
!            endif     
            
            !save water column before removes from next processes
            !important for property transport if river cells get out of water, conc has to be computed
            !after transport and not zero because there was no water left
            call SetMatrixValue(Me%myWaterColumnAfterTransport, Me%Size, Me%myWaterColumn)
            
            !Gets ExternalVars
            call ReadLockExternalVar (StaticOnly = .false.)

            !Flow through street gutter
            if (Me%StormWaterModel) then
                call StreetGutterFlow
            endif
                    
            !StormWater Drainage
            if (Me%StormWaterDrainage) then
                call StormWaterDrainage
            endif


            !Calculates flow from channels to land
            if (Me%ObjDrainageNetwork /= 0 .and. .not. Me%SimpleChannelInteraction) then
                call FlowFromChannels 
            endif

            !Calculates flow from channels to land and the other way round. New approach
            if (Me%ObjDrainageNetwork /=0 .and. Me%SimpleChannelInteraction) then
                call OverLandChannelInteraction_2
                
!                !call OverLandChannelInteraction_New
!                select case (Me%OverlandChannelInteractionMethod)
!                case (1)
!                    call OverLandChannelInteraction_1
!                case (2)
!                    call OverLandChannelInteraction_2
!                case (3)
!                    call OverLandChannelInteraction_3
!                case (4)
!                    call OverLandChannelInteraction_4
!                case default
!                    stop 'ModifyRunOff - ModuleRunOff - ERR020'
!                endselect
            endif
            
            !Routes Ponded levels which occour due to X/Y direction (Runoff does not route in D8)
            !the defaul method was celerity (it was corrected) but it ccould create high flow changes. Manning method is stabler
            !because of resistance. However in both methods the area used is not consistent (regular faces flow
            !already used all the cell vertical areas and the route D4 will overlapp areas - review this in the future
            if (Me%RouteDFourPoints) then
                if (Me%RouteDFourMethod == Manning_) then
                    call RouteDFourPoints
                elseif (Me%RouteDFourMethod == Celerity_) then
                    call RouteDFourPoints_v3
                endif
            endif

            !Boundary Condition
            !Only compute if case of waterlevel higher than boundary (overflow)
            !the default method was instantaneous flow (instantaneous go to boundary level)
            !but it was changed to compute flow (based on celerity) to be more consistent
            !with a free drop to boundary level (that can be much lower than topography)
            if (Me%ImposeBoundaryValue) then
                if (Me%BoundaryMethod == ComputeFlow_) then
                    call ImposeBoundaryValue
                elseif (Me%BoundaryMethod == InstantaneousFlow_) then
                    call ImposeBoundaryValue_v2
                endif
            endif

            !Calculates center flow and velocities (for output and next DT)
            call ComputeCenterValues

            call ComputeNextDT (Niter)                                    
            
            !Output Results
            if (Me%OutPut%Yes) then                   
                call RunOffOutput
            endif

            if (Me%WriteMaxWaterColumn) &
                call OutputMaxWaterColumn

            call CalculateTotalStoredVolume

            !Restart Output
            if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%Now == Me%EndTime)) then
                if(Me%ExtVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    call WriteFinalFile
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                endif
            endif

            !Ungets external variables
            call ReadUnLockExternalVar (StaticOnly = .false.)

            STAT_ = SUCCESS_
            if (MonitorPerformance) call StopWatch ("ModuleRunOff", "ModifyRunOff")

        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyRunOff
    
    !---------------------------------------------------------------------------

    subroutine ModifyWaterDischarges (LocalDT)

        !Arguments--------------------------------------------------------------
        real                                    :: LocalDT

        !Local------------------------------------------------------------------
        integer                                 :: iDis, nDischarges
        integer                                 :: i, j, k
        real                                    :: SurfaceElevation    
        real                                    :: Flow, MaxFlow   
        integer                                 :: STAT_CALL

        !Sets to 0
        call SetMatrixValue(Me%lFlowDischarge, Me%Size, 0.0)

        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR01'

        do iDis = 1, nDischarges

            call GetDischargesGridLocalization(Me%ObjDischarges,                        &
                                               DischargeIDNumber = iDis,                &
                                               Igrid = i,                               &
                                               JGrid = j,                               &
                                               KGrid = k,                               &
                                               STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR02'
            
            if (k == 0) then
                
                !real(8) to real as expected in GetDischargeWaterFlow
                SurfaceElevation = Me%myWaterLevel(i, j)
                call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                        Me%ExtVar%Now, iDis,                            &
                                        SurfaceElevation,                               &
                                        Flow, STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOff - ModifyWaterDischarges - ERR04'
                
                !each additional flow can remove all water column left
                if (Flow .lt. 0.0) then
                    !m3/s = m3 /s
                    MaxFlow = - Me%myWaterVolume(i, j) / LocalDT
                  
                    if (abs(Flow) .gt. abs(MaxFlow)) then
                        Flow = MaxFlow 
                    endif
                endif

                Me%lFlowDischarge(i, j)     = Me%lFlowDischarge(i, j) + Flow

                !Updates Water Volume
                Me%myWaterVolume(i, j)      = Me%myWaterVolume(i, j) + Flow * LocalDT

                !Updates Water Column
                Me%myWaterColumn  (i, j)    = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                !Updates Water Level
                Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j) + Me%ExtVar%Topography(i, j)

                !if (Me%CheckMass) Me%TotalInputVolume = Me%TotalInputVolume + Me%DischargesFlow(iDis) * LocalDT
 
            endif
           
        enddo

    end subroutine ModifyWaterDischarges  
    
    !--------------------------------------------------------------------------

    subroutine ModifyGeometryAndMapping

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: WCL, WCR, WCA, Bottom
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, WCL, WCR, WCA, Bottom)


        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j-1) == BasinPoint .and. Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                if (Me%FaceWaterColumn == WCMaxBottom_) then
                    !Maximum Bottom Level
                    Bottom = max(Me%ExtVar%Topography(i, j-1), Me%ExtVar%Topography(i, j))
                elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                    !Average Bottom Level
                    Bottom = (Me%ExtVar%Topography(i,j) + Me%ExtVar%Topography(i,j-1)) / 2.0                
                endif
                
                !Water Column Left (above MaxBottom)
                WCL       = max(Me%myWaterLevel(i, j-1) + Me%BuildingsHeight(i, j-1) - Bottom, dble(0.0))
            
                !Water Column Right (above MaxBottom)
                WCR       = max(Me%myWaterLevel(i, j  ) + Me%BuildingsHeight(i, j)   - Bottom, dble(0.0))

                !In the case of kinematic wave, always consider the "upstream" area, otherwise the average above "max bottom"
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                    if (Me%ExtVar%Topography(i, j-1) > Me%ExtVar%Topography(i, j)) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                else
                    !Average Water Column
                    !WCA = (WCL + WCR) / 2.0
                    if (Me%myWaterLevel(i, j-1) + Me%BuildingsHeight(i, j-1) >           &
                        Me%myWaterLevel(i, j) + Me%BuildingsHeight(i, j)) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                endif
                
                !Area  = Water Column * Side lenght of cell
                Me%AreaU(i, j) = WCA * Me%ExtVar%DYY(i, j)
                
                if (WCA > Me%MinimumWaterColumn) then
                    Me%ComputeFaceU(i, j) = 1
                else
                    Me%ComputeFaceU(i, j) = 0
                endif

            endif
            
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%BasinPoints(i-1, j) == BasinPoint .and. Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                if (Me%FaceWaterColumn == WCMaxBottom_) then
                    !Maximum Bottom Level
                    Bottom = max(Me%ExtVar%Topography(i-1, j), Me%ExtVar%Topography(i, j))
                elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                    !Average Bottom Level
                    Bottom = (Me%ExtVar%Topography(i,j) +  Me%ExtVar%Topography(i-1,j)) / 2.0                
                endif
                
                !Water Column Left
                WCL       = max(Me%myWaterLevel(i-1, j) + Me%BuildingsHeight(i-1, j) - Bottom, dble(0.0))
            
                !Water Column Right
                WCR       = max(Me%myWaterLevel(i, j  ) + Me%BuildingsHeight(i, j)   - Bottom, dble(0.0))
               
                !In the case of kinematic wave, always consider the "upstream" area, otherwise the average above "max bottom"
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                    if (Me%ExtVar%Topography(i-1, j) > Me%ExtVar%Topography(i, j)) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                else
                    !Average Water Column
                    !WCA = (WCL + WCR) / 2.0
                    if (Me%myWaterLevel(i-1, j) + Me%BuildingsHeight(i-1, j) >           &
                        Me%myWaterLevel(i, j) + Me%BuildingsHeight(i, j)) then
                        WCA = WCL
                    else
                        WCA = WCR
                    endif
                endif

                !Area  = Water Column * Side lenght of cell
                Me%AreaV(i, j) = WCA * Me%ExtVar%DXX(i, j)
                
                if (WCA > Me%MinimumWaterColumn) then
                    Me%ComputeFaceV(i, j) = 1
                else
                    Me%ComputeFaceV(i, j) = 0
                endif

            endif
            
        enddo
        enddo    
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
    
    end subroutine ModifyGeometryAndMapping

    !--------------------------------------------------------------------------

    subroutine KinematicWave ()
    
        !Arguments-------------------------------------------------------------
        !real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Slope
        real                                        :: level_left, level_right
        real                                        :: level_bottom, level_top
        real                                        :: HydraulicRadius, WetPerimeter
        real                                        :: Margin1, Margin2
        real                                        :: WaterDepth, MaxBottom
        integer                                     :: CHUNK, di, dj
        real(8)                                     :: MaxFlow
        !character(len=StringLength)                 :: Direction

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, Slope, level_left, level_right, level_bottom, level_top, &
        !$OMP HydraulicRadius, MaxFlow, Margin1, Margin2, WaterDepth, MaxBottom, WetPerimeter, di, dj)
        
        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaceU(i, j) == Compute) then
            
                !Adds to the final level the height of the buidings, if any
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                
                    if (Me%Buildings) then
                        level_left  = Me%ExtVar%Topography(i, j-1) + Me%BuildingsHeight(i, j-1)
                        level_right = Me%ExtVar%Topography(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_left  = Me%ExtVar%Topography(i, j-1)
                        level_right = Me%ExtVar%Topography(i, j)
                    endif
                    
                else if (Me%HydrodynamicApproximation == DiffusionWave_) then
                
                    if (Me%Buildings) then
                        level_left  = Me%myWaterLevel(i, j-1) + Me%BuildingsHeight(i, j-1)
                        level_right = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_left  = Me%myWaterLevel(i, j-1)
                        level_right = Me%myWaterLevel(i, j)
                    endif

                else
                
                    write(*,*)'Internal error'
                
                endif
                    
                !Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_left - level_right) / Me%ExtVar%DZX(i, j-1))
                else
                    Slope           = (level_left - level_right) / Me%ExtVar%DZX(i, j-1)
                endif
                
                !Hydraulic Radius
!                Direction = "X"
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_left,level_right)
                !Wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DYY(i, j)
                
                !only compute in water column as MaxBottom (topography stairs descritization)
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Water Depth consistent with AreaU computed (only water above max bottom)
                    WaterDepth = Me%AreaU(i,j) / Me%ExtVar%DYY(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))
                    
                    !to check wich cell to use to use since areaU depends on higher water level and max bottom
                    if (level_left .gt. level_right) then
                        dj = -1
                    else
                        dj = 0
                    endif
                   
                    !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i+1, j + dj) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i-1, j + dj) - MaxBottom

                    !if positive than there is a “margin” on the side and friction occurs at wet length
                    !If not basin points than result will be negative.
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                HydraulicRadius = Me%AreaU(i, j) / WetPerimeter
                             
                
                !
                !MANNING'S EQUATION -  KINEMATIC WAVE
                !
                !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                if (Slope >= 0.0) then
                    Me%lFlowX(i, j) = Me%AreaU(i, j) * HydraulicRadius**(2./3.) * sqrt(Slope)          &
                                      / Me%OverlandCoefficientX(i,j)
                else
                    Me%lFlowX(i, j) = - Me%AreaU(i, j) * HydraulicRadius**(2./3.) * sqrt(-1.0 * Slope) &
                                      / Me%OverlandCoefficientX(i,j)
                endif
                
                
                !Limits Velocity to celerity if a free drop exists
                if (Me%HydrodynamicApproximation == DiffusionWave_ .and. Me%LimitToCriticalFlow) then
                    if ((level_left .lt. Me%ExtVar%Topography(i,j)) .or. (level_right .lt. Me%ExtVar%Topography(i,j-1))) then
                        
                        !already defined in shorter
                        !WaterDepth = max (level_left, level_right) - max(Me%ExtVar%Topography(i, j-1), Me%ExtVar%Topography(i, j))
                        WaterDepth      = Me%AreaU(i, j)/Me%ExtVar%DYY(i,j)
                        MaxFlow         = Me%AreaU(i, j) * sqrt(Gravity * WaterDepth)
                        Me%lFlowX(i, j) = Min (MaxFlow, Me%lFlowX(i, j))       
                                    
                    endif
    
                endif
                
            else
                
                Me%lFlowX(i, j) = 0.0
            
            endif
                
        enddo
        enddo        
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ComputeFaceV(i, j) == Compute) then
            
                !Adds to the final level the height of the buidings, if any
                if (Me%HydrodynamicApproximation == KinematicWave_) then
                
                    !Adds to the final level the height of the buidings, if any
                    if (Me%Buildings) then
                        level_bottom = Me%ExtVar%Topography(i-1, j) + Me%BuildingsHeight(i-1, j)
                        level_top    = Me%ExtVar%Topography(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_bottom = Me%ExtVar%Topography(i-1, j)
                        level_top    = Me%ExtVar%Topography(i, j)
                    endif
                    
                else if (Me%HydrodynamicApproximation == DiffusionWave_) then

                    !Adds to the final level the height of the buidings, if any
                    if (Me%Buildings) then
                        level_bottom = Me%myWaterLevel(i-1, j) + Me%BuildingsHeight(i-1, j)
                        level_top    = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_bottom = Me%myWaterLevel(i-1, j)
                        level_top    = Me%myWaterLevel(i, j)
                    endif
                
                else
                
                    write(*,*)'Internal error'
                
                endif

                
                !Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_bottom - level_top) / Me%ExtVar%DZY(i-1, j))
                else
                    Slope           = (level_bottom - level_top) / Me%ExtVar%DZY(i-1, j)
                endif
                
                !Hydraulic Radius
!                Direction = "Y"
!               !This function produced an overhead in openmp and the simulation took
!               !double the time so it was abandoned
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_bottom,level_top)                
                !Wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DXX(i, j)
                
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Water Depth consistent with AreaV computed (only water above max bottom)
                    WaterDepth = Me%AreaV(i,j) / Me%ExtVar%DXX(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))

                    !to check wich cell to use since areaV depends on higher water level
                    if (level_bottom .gt. level_top) then
                        di = -1
                    else
                        di = 0
                    endif

                    !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i + di,j+1) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i + di,j-1) - MaxBottom

                    !if positive than there is a “margin” on the side and friction occurs at wet length
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                !m = m2 / m
                HydraulicRadius = Me%AreaV(i, j) / WetPerimeter
                
                                
                !
                !MANNING'S EQUATION -  KINEMATIC WAVE
                !
                !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                if (Slope >= 0.0) then
                    Me%lFlowY(i, j) = Me%AreaV(i, j) * HydraulicRadius**(2./3.) * sqrt(Slope)            &
                                      / Me%OverlandCoefficientY(i,j)
                else
                    Me%lFlowY(i, j) = - Me%AreaV(i, j) * HydraulicRadius**(2./3.) * sqrt(-1.0 * Slope)   &
                                      / Me%OverlandCoefficientY(i,j)
                endif
                
                !Limits Velocity to reasonable values
                if (Me%HydrodynamicApproximation == DiffusionWave_ .and. Me%LimitToCriticalFlow) then

                    if ((level_bottom .lt. Me%ExtVar%Topography(i,j)) .or. (level_top .lt. Me%ExtVar%Topography(i-1,j))) then
                        
                        !already defined in shorter
                        !WaterDepth = max (level_bottom, level_top) - max(Me%ExtVar%Topography(i-1,j), Me%ExtVar%Topography(i, j))
                        WaterDepth      = Me%AreaV(i, j)/Me%ExtVar%DXX(i,j)
                        MaxFlow         = Me%AreaV(i, j) * sqrt(Gravity * WaterDepth)
                        Me%lFlowY(i, j) = Min (MaxFlow, Me%lFlowY(i, j))
                    
                    endif

                
                
                endif
                
            else
            
                Me%lFlowY(i, j) = 0.0
            
            endif

        enddo
        enddo              
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    
    end subroutine KinematicWave
    
    !--------------------------------------------------------------------------
    
    subroutine DynamicWaveXX (LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Slope
        real                                        :: level_left, level_right
        real                                        :: HydraulicRadius
        real                                        :: Friction
        real                                        :: Pressure
        !real                                        :: upAdv, downAdv, 
        real                                        :: XLeftAdv, XRightAdv, YBottomAdv, YTopAdv
        real                                        :: Advection, Qf, WetPerimeter
        real(8)                                     :: CriticalFlow
        real                                        :: Margin1, Margin2
        integer                                     :: CHUNK, dj
        real                                        :: MaxBottom, WaterDepth
        !character(len=StringLength)                 :: Direction


        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "DynamicWaveXX")


        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, Slope, level_left, level_right, &
        !$OMP HydraulicRadius, Friction, Pressure, XLeftAdv, XRightAdv, YBottomAdv, YTopAdv, Advection, Qf, &
        !$OMP CriticalFlow, Margin1, Margin2, MaxBottom, WaterDepth, dj, WetPerimeter)

        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%ComputeFaceU(i, j) == Compute) then
            
                !Adds to the final level the height of the buidings, if any
                if (Me%Buildings) then
                    level_left  = Me%myWaterLevel(i, j-1) + Me%BuildingsHeight(i, j-1)
                    level_right = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                else
                    level_left  = Me%myWaterLevel(i, j-1)
                    level_right = Me%myWaterLevel(i, j)
                endif
                    
                !!Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_left - level_right) / Me%ExtVar%DZX(i, j-1))
                else
                    Slope           = (level_left - level_right) / Me%ExtVar%DZX(i, j-1)
                endif
                 
                !!Hydraulic Radius
!                Direction = "X"
!                !This function produced an overhead in openmp and the simulation took 
!                !double the time so it was abandoned
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_left,level_right)
                !wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DYY(i, j)
                
                !only compute margins if water column method is MaxBottom (topography discretization by "stairs")
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Then, is checked if "margins" occur on the cell of the highest water level
                    !water depth consistent with AreaU computed (only water above max bottom)
                    WaterDepth = Me%AreaU(i,j) / Me%ExtVar%DYY(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))
                    
                    !to check which cell to use since areaU depends on higher water level
                    if (level_left .gt. level_right) then
                        dj = -1
                    else
                        dj = 0
                    endif
                    
                    !bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i+1, j + dj) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i-1, j + dj) - MaxBottom

                    !if positive, than there is a “margin” on the side and friction occurs at wet length
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                HydraulicRadius = Me%AreaU(i, j) / WetPerimeter
       
                !
                !Sant Venant
                !

                !Pressure
                !m3/s             = s  * m/s2    * m2   * m/m
                Pressure          = LocalDT * Gravity * Me%AreaU(i, j) * Slope


                !FRICTION - semi-implicit -----------------------------------------------
                !   -    =  (s * m.s-2  * m3.s-1 * (s.m(-1/3))^2) / (m2 * m(4/3)) = m(10/3) / m(10/3)
                Friction = LocalDT * Gravity * abs(Me%FlowXOld(i, j)) * Me%OverlandCoefficientX(i,j)** 2. &
                         / ( Me%AreaU(i, j) * HydraulicRadius ** (4./3.) ) 
        
                
                !Advection (may be limited to water column height)
                if ((Me%CalculateAdvection) .and. (Me%myWaterColumn(i,j) .gt. Me%MinimumWaterColumnAdvection)   &
                     .and. (Me%myWaterColumn(i,j-1) .gt. Me%MinimumWaterColumnAdvection)) then
                    
                    !Face XU(i,j+1). Z U Faces have to be open
                    if ((Me%ComputeFaceU(i, j) +  Me%ComputeFaceU(i, j+1) == 2)) then 

                        !OLD Version
                        !Theold formulation had a problem when flows in adjacent reaches
                        !had opposite directions. Flow was the average and velocity would be
                        !in opposite  direction of average flow.
                       
                        !New Version 
                        !The new formulation, in case of opposite directions, in adjacent reaches does not compute
                        !advection. In case of same direction, is hard-upwind meaning that it will use flow and 
                        !velocity from the upwind reach. This option may be more stable than soft-upwind 
                        !(average flow and velocity from upwind reach) or central differences (average flow 
                        !and average velocity).
                        !if flows in same direction, advection is computed                        
                        if ((Me%FlowXOld(i, j) * Me%FlowXOld(i, j+1)).ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j) + Me%FlowXOld(i, j+1)) / 2.0

                            if (Qf > 0.0) then
                                XRightAdv = Me%FlowXOld(i, j)   * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            else
                                XRightAdv = Me%FlowXOld(i, j+1) * Me%FlowXOld(i, j+1) / Me%AreaU(i, j+1)
                            endif
                        else
                            XRightAdv = 0.0
                        endif
                        
                    else
                        XRightAdv = 0.0
                    endif       
                    
                    !Face XU(i,j). Z U Faces have to be open
                    if ((Me%ComputeFaceU(i, j-1) + Me%ComputeFaceU(i, j) == 2)) then  
                        
                        !New Version
                        if ((Me%FlowXOld(i, j-1) * Me%FlowXOld(i, j)) .ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j-1) + Me%FlowXOld(i, j)) / 2.0

                            if (Qf > 0.0) then
                                XLeftAdv = Me%FlowXOld(i, j-1) * Me%FlowXOld(i, j-1) / Me%AreaU(i, j-1)
                            else
                                XLeftAdv = Me%FlowXOld(i, j) * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            endif
                        else
                            XLeftAdv = 0.0
                        endif
                        
                        
                    else
                        XLeftAdv = 0.0
                    endif       
                    
                    !Faces of U(i,j) that were not being accounted (in 2D need to be accounted)
                    !Face YU(i+1,j)
!                    if (Me%ComputeFaceV(i+1, j-1) +  Me%ComputeFaceV(i+1, j)     &
!                        + Me%ComputeFaceU(i, j) +  Me%ComputeFaceU(i+1, j) == 4) then
                    if (Me%ComputeFaceV(i+1, j-1) +  Me%ComputeFaceV(i+1, j) .ge. 1) then

                       !if flows in same direction, advection is computed                        
                        if ((Me%FlowYOld(i+1, j-1) * Me%FlowYOld(i+1, j)).ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i+1, j-1) + Me%FlowYOld(i+1, j)) / 2.0
                            
                            if ((Qf > 0.0)) then
                                YTopAdv = Qf   * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            elseif ((Qf < 0.0) .and. (Me%ComputeFaceU(i+1,j) == Compute)) then
                                YTopAdv = Qf   * Me%FlowXOld(i+1, j) / Me%AreaU(i+1, j)
                            else
                                YTopAdv = 0.0
                            endif
                        else
                            YTopAdv = 0.0
                        endif
                        
                    else
                        YTopAdv = 0.0
                    endif       
                    
                    !Faces of U(i,j) that were not being accounted (in 2D need to be accounted)
                    !Face YU(i,j)
!                    if (Me%ComputeFaceV(i, j-1) +  Me%ComputeFaceV(i, j)      &
!                        + Me%ComputeFaceU(i, j) +  Me%ComputeFaceU(i-1, j) == 4) then
                    if (Me%ComputeFaceV(i, j-1) +  Me%ComputeFaceV(i, j) .ge. 1) then

                       !if flows in same direction, advection is computed                        
                        if ((Me%FlowYOld(i, j-1) * Me%FlowYOld(i, j)).ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i, j-1) + Me%FlowYOld(i, j)) / 2.0

                            if ((Qf > 0.0) .and. (Me%ComputeFaceU(i-1,j) == Compute)) then
                                YBottomAdv =  Qf   * Me%FlowXOld(i-1, j) / Me%AreaU(i-1, j)
                            elseif ((Qf < 0.0)) then
                                YBottomAdv = Qf   * Me%FlowXOld(i, j) / Me%AreaU(i, j)
                            else
                                YBottomAdv = 0.0
                            endif
                        else
                            YBottomAdv = 0.0
                        endif
                        
                    else
                        YBottomAdv = 0.0
                    endif       

                    !Advection = (upAdv - downAdv) * LocalDT / Me%ExtVar%DUX(i, j)
                    Advection = (XLeftAdv - XRightAdv) * LocalDT / Me%ExtVar%DZX(i, j-1)     &
                                + (YBottomAdv - YTopAdv) * LocalDT / Me%ExtVar%DYY(i, j)
                                
                else
                
                    Advection = 0.0
                    
                endif
                
                Me%lFlowX(i, j) = (Me%FlowXOld(i, j) + Pressure + Advection) / (1.0 + Friction)

                if (Me%LimitToCriticalFlow) then

                    !Limit to critical flow. Using the critical flow limitation in all cells assumes "slow" flow or
                    !subcritical that is consistent with the formulation used (flow depends on downstream height)
                    !because in supercritical flow it is only dependent on upstream and descritization to describe it would have
                    !to change. Supercritical flow usually exists on hydraulic infraestructures (high drops) and a 
                    !hydraulic jump exists between fast flow and slow flow.
                    
                    !Test Limitation only if free drop exists
!                    if ((level_left .lt. Me%ExtVar%Topography(i,j)) .or. (level_right .lt. Me%ExtVar%Topography(i,j-1))) then

                        !Waterdepth at the center of the face - depending on flow direction since flow
                        !can be in opposite direction of height gradient (AreaU uses the higher water level)              
                        !WaterDepth = Me%AreaU(i,j)/Me%ExtVar%DYY(i,j)
                        if (Me%FaceWaterColumn == WCMaxBottom_) then
                            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))     
                                                            
                            if (Me%lFlowX(i, j) .gt. 0.0) then           
                                WaterDepth = max(Me%MyWaterLevel(i,j-1) - MaxBottom, 0.0)
                            else
                                WaterDepth = max(Me%MyWaterLevel(i,j) - MaxBottom, 0.0)
                            endif
                        elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                            if (Me%lFlowX(i, j) .gt. 0.0) then           
                                WaterDepth = Me%MyWaterColumn(i,j-1)
                            else
                                WaterDepth = Me%MyWaterColumn(i,j)
                            endif                        
                        endif
                        
                        !Critical Flow
                        !CriticalFlow = Me%AreaU(i, j) * sqrt(Gravity * WaterDepth)
                        !m3/s = m * m * m/s
                        CriticalFlow = WaterDepth * Me%ExtVar%DYY(i,j) * sqrt(Gravity * WaterDepth)
                        
                        !only limit if flow higher
                        if (abs(Me%lFlowX(i, j)) > CriticalFlow) then
                            if (Me%lFlowX(i, j) > 0) then
                                Me%lFlowX(i, j) = CriticalFlow
                            else
                                Me%lFlowX(i, j) = -1.0 * CriticalFlow
                            endif
                        endif
 !                   endif
                
                else
                    !Predict water column to avoid negative volumes since 4 fluxes exist and the sum may be more than exists
                    if (Me%lFlowX(i, j) .lt. 0.0) then
                        if (abs(Me%lFlowX(i, j))* LocalDT .gt. Me%myWaterVolumePred(i,j)) then
                            Me%lFlowX(i, j) = - Me%myWaterVolumePred(i,j) / LocalDT
                        endif
                    elseif (Me%lFlowX(i, j) .gt. 0.0) then
                        if (Me%lFlowX(i, j)* LocalDT .gt. Me%myWaterVolumePred(i,j-1)) then
                            Me%lFlowX(i, j) =  Me%myWaterVolumePred(i,j-1) / LocalDT
                        endif
                    endif 
                    
                    !m3 = m3 + (-m3/s * s)
                    Me%myWaterVolumePred(i,j) = Me%myWaterVolumePred(i,j) + (Me%lFlowX(i, j) * LocalDT)
                    Me%myWaterVolumePred(i,j-1) = Me%myWaterVolumePred(i,j-1) - (Me%lFlowX(i, j) * LocalDT) 
                                   
                endif
                
            else
            
                Me%lFlowX(i, j) = 0.0

            endif

        enddo
        enddo        
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "DynamicWaveXX")
        
        
    end subroutine DynamicWaveXX
    
    !--------------------------------------------------------------------------

    subroutine DynamicWaveYY (LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Slope
        real                                        :: level_bottom, level_top
        real                                        :: HydraulicRadius
        real                                        :: Friction
        real                                        :: Pressure
        !real                                        :: upAdv, downAdv, 
        real                                        :: XLeftAdv, XRightAdv, YBottomAdv, YTopAdv
        real                                        :: Advection, Qf, WetPerimeter
        real(8)                                     :: CriticalFlow
        real                                        :: Margin1, Margin2
        integer                                     :: CHUNK, di
        real                                        :: MaxBottom, WaterDepth
        !character(len=StringLength)                 :: Direction

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "DynamicWaveYY")


        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, Slope, level_bottom, level_top, &
        !$OMP HydraulicRadius, Friction, Pressure, XLeftAdv, XRightAdv, YBottomAdv, YTopAdv, Advection, Qf, &
        !$OMP CriticalFlow, Margin1, Margin2, MaxBottom, WaterDepth, di, WetPerimeter)

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ComputeFaceV(i, j) == Compute) then
            
                !Adds to the final level the height of the buidings, if any
                if (Me%Buildings) then
                    level_bottom = Me%myWaterLevel(i-1, j) + Me%BuildingsHeight(i-1, j)
                    level_top    = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                else
                    level_bottom = Me%myWaterLevel(i-1, j)
                    level_top    = Me%myWaterLevel(i, j)
                endif
                
                !!Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_bottom - level_top) / Me%ExtVar%DZY(i-1, j))
                else
                    Slope           = (level_bottom - level_top) / Me%ExtVar%DZY(i-1, j)
                endif
                
                !!Hydraulic Radius
!                Direction = "Y"
!                !This function produced an overhead with openmp so it was abandoned
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_bottom,level_top)
                
                !wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DXX(i, j)
                
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !water Depth consistent with AreaV computed (only water above max bottom)
                    WaterDepth = Me%AreaV(i,j) / Me%ExtVar%DXX(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))

                    !to check wich cell to use since areaV depends on higher water level
                    if (level_bottom .gt. level_top) then
                        di = -1
                    else
                        di = 0
                    endif

                    !bottom Difference to adjacent cells (to check existence of “margins” on the side)
                    Margin1 = Me%ExtVar%Topography(i + di,j+1) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i + di,j-1) - MaxBottom

                    !if positive than there is a “margin” on the side and friction occurs at wet length
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif
                
                !m = m2 / m
                HydraulicRadius = Me%AreaV(i, j) / WetPerimeter
               
                !
                !Sant Venant
                !

                !m3/s             = s  * m/s2    * m2   * m/m
                Pressure          = LocalDT * Gravity * Me%AreaV(i, j) * Slope


                !FRICTION - semi-implicit -----------------------------------------------
                !   -    =  (s * m.s-2  * m3.s-1 * (s.m(-1/3))^2) / (m2 * m(4/3)) = m(10/3) / m(10/3)
                Friction = LocalDT * Gravity * abs(Me%FlowYOld(i, j)) * Me%OverlandCoefficientY(i,j) ** 2. &
                         / ( Me%AreaV(i, j) * HydraulicRadius ** (4./3.) ) 
        

                !Advection
                if ((Me%CalculateAdvection) .and. (Me%myWaterColumn(i,j) .gt. Me%MinimumWaterColumnAdvection)  &
                     .and. (Me%myWaterColumn(i-1,j) .gt. Me%MinimumWaterColumnAdvection)) then
                    
                    !Face YV(i+1,j)
                    if ((Me%ComputeFaceV(i, j) +  Me%ComputeFaceV(i+1, j) == 2)) then 
                        
                        if ((Me%FlowYOld(i, j) * Me%FlowYOld(i+1, j)) .ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i, j) + Me%FlowYOld(i+1, j)) / 2.0

                            if (Qf > 0.0) then
                                YTopAdv = Me%FlowYOld(i, j)   * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                            else
                                YTopAdv = Me%FlowYOld(i+1, j) * Me%FlowYOld(i+1, j) / Me%AreaV(i+1, j)
                            endif
                        else
                            YTopAdv = 0.0
                        endif
                        
                    else
                        YTopAdv = 0.0
                    endif
                    
                    !Face YV(i,j)
                    if ((Me%ComputeFaceV(i-1, j) + Me%ComputeFaceV(i, j) == 2)) then 

                        if ((Me%FlowYOld(i-1, j) * Me%FlowYOld(i, j)) .ge. 0.0) then
                            
                            Qf = (Me%FlowYOld(i-1, j) + Me%FlowYOld(i, j)) / 2.0

                            if (Qf > 0.0) then
                                YBottomAdv = Me%FlowYOld(i-1, j)   * Me%FlowYOld(i-1, j) / Me%AreaV(i-1, j)
                            else
                                YBottomAdv = Me%FlowYOld(i, j) * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                            endif
                        else
                            YBottomAdv = 0.0
                        endif
                        
                    else
                        YBottomAdv = 0.0
                    endif                
        
                    !Faces of V(i,j) that were not being accounted (in 2D need to be accounted)
                    !Face XV(i,j+1)
!                    if (Me%ComputeFaceU(i, j+1) +  Me%ComputeFaceU(i-1, j+1)     &
!                        + Me%ComputeFaceV(i, j) +  Me%ComputeFaceV(i, j+1) == 4) then
                    if (Me%ComputeFaceU(i, j+1) +  Me%ComputeFaceU(i-1, j+1) .ge. 1) then
                        
                       !if flows in same direction, advection is computed                        
                        if ((Me%FlowXOld(i, j+1) * Me%FlowXOld(i-1, j+1)).ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j+1) + Me%FlowXOld(i-1, j+1)) / 2.0

                            if ((Qf > 0.0)) then
                                XRightAdv = Qf   * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                            elseif ((Qf < 0.0) .and. (Me%ComputeFaceV(i,j+1) == Compute)) then
                                XRightAdv = Qf   * Me%FlowYOld(i, j+1) / Me%AreaV(i, j+1)
                            else 
                                XRightAdv = 0.0
                            endif
                        else
                            XRightAdv = 0.0
                        endif
                        
                    else
                        XRightAdv = 0.0
                    endif       

                    !Faces of V(i,j) that were not being accounted (in 2D need to be accounted)
                    !Face XV(i,j)
!                    if (Me%ComputeFaceU(i, j) +  Me%ComputeFaceU(i-1, j)      &
!                        + Me%ComputeFaceV(i, j) +  Me%ComputeFaceV(i, j-1) == 4) then
                    if (Me%ComputeFaceU(i, j) +  Me%ComputeFaceU(i-1, j) .ge. 1) then
                        
                       !if flows in same direction, advection is computed                        
                        if ((Me%FlowXOld(i, j) * Me%FlowXOld(i-1, j)).ge. 0.0) then
                            
                            Qf = (Me%FlowXOld(i, j) + Me%FlowXOld(i-1, j)) / 2.0

                            if ((Qf > 0.0) .and. (Me%ComputeFaceV(i,j-1) == Compute)) then
                                XLeftAdv = Qf   * Me%FlowYOld(i, j-1) / Me%AreaV(i, j-1)
                            elseif ((Qf < 0.0)) then
                                XLeftAdv = Qf   * Me%FlowYOld(i, j) / Me%AreaV(i, j)
                            else
                                XLeftAdv = 0.0
                            endif
                        else
                            XLeftAdv = 0.0
                        endif
                        
                    else
                        XLeftAdv = 0.0
                    endif       
                           
                    !Advection = (upAdv - downAdv) * LocalDT / Me%ExtVar%DVY(i, j)
                    Advection = (YBottomAdv - YTopAdv) * LocalDT / Me%ExtVar%DZY(i-1, j)     &
                                + (XLeftAdv - XRightAdv) * LocalDT / Me%ExtVar%DXX(i, j)
                    
                else
                
                    Advection = 0.0
                    
                endif
                
                Me%lFlowY(i, j) = (Me%FlowYOld(i, j) + Pressure + Advection) / (1.0 + Friction)
                
                
                if (Me%LimitToCriticalFlow) then
                    
!                    if ((level_bottom .lt. Me%ExtVar%Topography(i,j)) .or. (level_top .lt. Me%ExtVar%Topography(i-1,j))) then
                    
                        !Waterdepth at the center of the face - depending on flow direction since flow
                        !can be in opposite direction of height gradient (AreaU uses the higher)
                        !WaterDepth = Me%AreaV(i,j)/Me%ExtVar%DXX(i,j)
                        if (Me%FaceWaterColumn == WCMaxBottom_) then
                            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))
                                                            
                            if (Me%lFlowY(i, j) .gt. 0.0) then           
                                WaterDepth = max(Me%MyWaterLevel(i-1,j) - MaxBottom, 0.0)
                            else
                                WaterDepth = max(Me%MyWaterLevel(i,j) - MaxBottom, 0.0)
                            endif                
                        elseif (Me%FaceWaterColumn == WCAverageBottom_) then
                            if (Me%lFlowY(i, j) .gt. 0.0) then           
                                WaterDepth = Me%MyWaterColumn(i-1,j)
                            else
                                WaterDepth = Me%MyWaterColumn(i,j)
                            endif                                        
                        endif
                        
                        !Critical Flow
                        !CriticalFlow = Me%AreaV(i, j) * sqrt(Gravity * WaterDepth)
                        !m3/s = m * m * m/s
                        CriticalFlow = WaterDepth * Me%ExtVar%DXX(i,j) * sqrt(Gravity * WaterDepth)
                        
                        !only limit if flow higher
                        if (abs(Me%lFlowY(i, j)) > CriticalFlow) then
                            if (Me%lFlowY(i, j) > 0) then
                                Me%lFlowY(i, j) = CriticalFlow
                            else
                                Me%lFlowY(i, j) = -1.0 * CriticalFlow
                            endif
                        endif
 !                   endif
                
                else
                    if (Me%lFlowY(i, j) .lt. 0.0) then
                        if ( abs(Me%lFlowY(i, j))* LocalDT  .gt. Me%myWaterVolumePred(i,j)) then
                            Me%lFlowY(i, j) = - Me%myWaterVolumePred(i,j)  / LocalDT
                        endif
                    elseif (Me%lFlowY(i, j) .gt. 0.0) then
                        if ( Me%lFlowY(i, j)* LocalDT .gt. Me%myWaterVolumePred(i-1,j)) then
                            Me%lFlowY(i, j) =  Me%myWaterVolumePred(i-1,j) / LocalDT
                        endif
                    endif                
                    
                    Me%myWaterVolumePred(i,j) = Me%myWaterVolumePred(i,j) +  (Me%lFlowY(i, j) * LocalDT)                        
                    Me%myWaterVolumePred(i-1,j) = Me%myWaterVolumePred(i-1,j) - (Me%lFlowX(i, j) * LocalDT) 
                    
                endif

            else
            
                Me%lFlowY(i, j) = 0.0
            
            endif
            
        enddo
        enddo         
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleRunOff", "DynamicWaveYY")
        
        
    end subroutine DynamicWaveYY
    
    !-------------------------------------------------------------------------
    
!    real function HydraulicRadius(i,j,Direction, level_before, level_after)
!        
!        !Arguments-------------------------------------------------------------
!        integer                                           :: i,j
!        character(len=StringLength)                       :: Direction
!        real                                              :: level_before, level_after
!        !Local-----------------------------------------------------------------
!        real                                              :: WetPerimeter, WaterDepth, MaxBottom
!        real                                              :: Margin1, Margin2
!        integer                                           :: di, dj
!       
!        
!        if(Direction == "X") then
!        
!            !Hydraulic Radius
!            !Wet perimeter, first is bottom
!            WetPerimeter = Me%ExtVar%DYY(i, j)
!
!            !Water Depth consistent with AreaU computed (only water above max bottom)
!            WaterDepth = Me%AreaU(i,j) / Me%ExtVar%DYY(i, j)
!            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))
!            
!            !to check wich cell to use to use since areaU depends on higher water level
!            if (level_before .gt. level_after) then
!                dj = -1
!            else
!                dj = 0
!            endif
!            
!            !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
!            Margin1 = Me%ExtVar%Topography(i+1, j + dj) - MaxBottom
!            Margin2 = Me%ExtVar%Topography(i-1, j + dj) - MaxBottom
!
!            !if positive than there is a “margin” on the side and friction occurs at wet length
!            !If not basin points than result will be negative.
!            if (Margin1 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
!            endif
!            if (Margin2 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
!            endif
!
!            HydraulicRadius = Me%AreaU(i, j) / WetPerimeter
!                
!        elseif (Direction == "Y") then
!        
!            !Hydraulic Radius
!            !Wet perimeter, first is bottom
!            WetPerimeter = Me%ExtVar%DXX(i, j)
!
!            !Water Depth consistent with AreaV computed (only water above max bottom)
!            WaterDepth = Me%AreaV(i,j) / Me%ExtVar%DXX(i, j)
!            MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))
!
!            !to check wich cell to use since areaV depends on higher water level
!            if (level_before .gt. level_after) then
!                di = -1
!            else
!                di = 0
!            endif
!
!            !Bottom Difference to adjacent cells (to check existence of “margins” on the side)
!            Margin1 = Me%ExtVar%Topography(i + di,j+1) - MaxBottom
!            Margin2 = Me%ExtVar%Topography(i + di,j-1) - MaxBottom
!
!            !if positive than there is a “margin” on the side and friction occurs at wet length
!            if (Margin1 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
!            endif
!            if (Margin2 .gt. 0.0) then
!                WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
!            endif
!            
!            !m = m2 / m
!            HydraulicRadius = Me%AreaV(i, j) / WetPerimeter        
!        endif
!        
!    end function HydraulicRadius
    
    !---------------------------------------------------------------------------
    
    subroutine UpdateWaterLevels(LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: dVol
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


        !X
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaceU(i, j) == BasinPoint) then
                
                !dVol
                dVol                       = Me%lFlowX(i, j) * LocalDT
                
                !Updates Water Volume
                Me%myWaterVolume (i, j-1)  = Me%myWaterVolume (i, j-1) - dVol 
                Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   + dVol 

                !Updates Water Column
                Me%myWaterColumn  (i, j-1) = Me%myWaterVolume (i, j-1) / Me%ExtVar%GridCellArea(i, j-1)
                Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                !Updates Water Level
                Me%myWaterLevel (i, j-1)   = Me%myWaterColumn (i, j-1) + Me%ExtVar%Topography(i, j-1)
                Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)

            endif
        enddo
        enddo

        !Y
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ComputeFaceV(i, j) == BasinPoint) then
                
                !dVol
                dVol                      = Me%lFlowY(i, j) * LocalDT
                
                !Updates Water Volume
                Me%myWaterVolume (i-1, j) = Me%myWaterVolume (i-1, j) - dVol 
                Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   + dVol 

                !Updates Water Column
                Me%myWaterColumn (i-1, j) = Me%myWaterVolume (i-1, j) / Me%ExtVar%GridCellArea(i-1, j)
                Me%myWaterColumn (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                !Updates Water Level
                Me%myWaterLevel (i-1, j)  = Me%myWaterColumn (i-1, j) + Me%ExtVar%Topography(i-1, j)
                Me%myWaterLevel (i, j)    = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)

            endif
        enddo
        enddo
    
    end subroutine UpdateWaterLevels 

    !--------------------------------------------------------------------------
    
    !old routine where flux is not taken into account level difference
    !and had an error where max volume was compared to a flow and not volume (fixed)
    subroutine RouteDFourPoints_v2
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, it, jt
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: dVol, AverageCellLength, FlowMax
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%DFourSinkPoint(i, j) == BasinPoint .and. Me%LowestNeighborI(i, j) /= null_int .and. &
                Me%myWaterColumn(i,  j) > Me%MinimumWaterColumn) then
            

                it = Me%LowestNeighborI(i, j)
                jt = Me%LowestNeighborJ(i, j)
                               
                !Critical Flow                    
                AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0
                !FlowMax = Min(sqrt(Gravity * Me%myWaterColumn(i, j)) *  Me%myWaterColumn(i, j) * AverageCellLength, &
                !              0.1 * Me%myWaterColumn(i, j) * AverageCellLength)
                
                ![m3/s] = [m/s] * [m] * [m]
                FlowMax = sqrt(Gravity * Me%myWaterColumn(i, j)) *  Me%myWaterColumn(i, j) * AverageCellLength
                

                !dVol -> max Critical Flow & Avaliable Volume
                !there was an error in units Flowmax is m3/s and not m3
                !dVol = min(Me%myWaterVolume(i,j), FlowMax)
                ![m3] = [m3/s] * [s]
                dVol = min(Me%myWaterVolume(i,j), FlowMax * Me%ExtVar%DT)
                
                !Updates Water Volume
                Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - dVol 
                Me%myWaterVolume (it, jt)  = Me%myWaterVolume (it, jt) + dVol 

                !Updates Water Column
                Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterColumn  (it, jt) = Me%myWaterVolume (it, jt) / Me%ExtVar%GridCellArea(it, jt)

                !Updates Water Level
                Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                Me%myWaterLevel (it, jt)   = Me%myWaterColumn (it, jt) + Me%ExtVar%Topography(it, jt)

            
            endif
            
        enddo
        enddo
    

    end subroutine RouteDFourPoints_v2

    !--------------------------------------------------------------------------
    
    !new routine where dh is used and only dh may move not all water column. 
    !and water moves in level gradient and not always doenstream
    subroutine RouteDFourPoints_v3
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, it, jt
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: AverageCellLength, Flow, MaxFlow
        real                                        :: WaveHeight, Celerity, dh
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%DFourSinkPoint(i, j) == BasinPoint .and. Me%LowestNeighborI(i, j) /= null_int)  then         

                it = Me%LowestNeighborI(i, j)
                jt = Me%LowestNeighborJ(i, j)
                
                !topography of cell i,j is always higher than it, jt (is the max bottom)
                WaveHeight =  max(Me%myWaterLevel(i, j), Me%myWaterLevel(it,jt)) - Me%ExtVar%Topography(i,j)
                Celerity   = sqrt(Gravity * WaveHeight)

                if (WaveHeight .gt. Me%MinimumWaterColumn) then
                
                    !Critical Flow                    
                    AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0                

                    !dh>0 flow removes water, dh<0 flow brings water
                    dh =  Me%myWaterLevel(i, j) - Me%myWaterLevel(it,jt)
                    
                    !m3/s = m/s * m * m. if dh negative minimum is dh
                    Flow = Celerity *  min(dh, WaveHeight) * AverageCellLength
                    
                    !Max flow is volume given by area * dh
                    !Since it jt has always lower topography if dh negative there is not the
                    !possibility of using an abs(dh) higher than Waveheight (more flux than exists)
                    !if positive dh minimum is positive, if dh negative, negative flux with dh
                    MaxFlow = min(dh, WaveHeight) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                    
                    if (abs(Flow) > abs(MaxFlow)) then
                        Flow = MaxFlow
                    endif
                    
                    Me%iFlowRouteDFour(i,j)    = Flow

                    !Updates Water Volume
                    Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - Flow *  Me%ExtVar%DT
                    Me%myWaterVolume (it, jt)  = Me%myWaterVolume (it, jt) + Flow *  Me%ExtVar%DT 

                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)
                    Me%myWaterColumn  (it, jt) = Me%myWaterVolume (it, jt) / Me%ExtVar%GridCellArea(it, jt)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    Me%myWaterLevel (it, jt)   = Me%myWaterColumn (it, jt) + Me%ExtVar%Topography(it, jt)
                
                else
                    Me%iFlowRouteDFour(i,j)    = 0.0
                endif
                
            endif
            
        enddo
        enddo
    

    end subroutine RouteDFourPoints_v3

    !--------------------------------------------------------------------------

    !new routine where flow is computed from manning. 
    subroutine RouteDFourPoints
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, it, jt, di, dj
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Flow, MaxFlow, dx, dy
        real                                        :: AverageCellLengthSink, AverageCellLengthLower
        real                                        :: WaveHeight, sign !, Celerity
        real                                        :: level_up, level_down, CenterDistance
        real                                        :: Slope, VertArea !, WetPerimeter
        real                                        :: HydraulicRadius, OverlandCoef
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%DFourSinkPoint(i, j) == BasinPoint .and. Me%LowestNeighborI(i, j) /= null_int)  then         

                it = Me%LowestNeighborI(i, j)
                jt = Me%LowestNeighborJ(i, j)
                
                !topography of cell i,j is always higher than it, jt (is the max bottom)
                WaveHeight =  max(Me%myWaterLevel(i, j), Me%myWaterLevel(it,jt)) - Me%ExtVar%Topography(i,j)

                if (WaveHeight .gt. Me%MinimumWaterColumn) then                
                    
                    !applyng manning equation
                    if (Me%Buildings) then
                        level_up   = Me%myWaterLevel(i, j  )  + Me%BuildingsHeight(i, j)
                        level_down = Me%myWaterLevel(it, jt)  + Me%BuildingsHeight(it, jt)
                    else
                        level_up   = Me%myWaterLevel(i, j)
                        level_down = Me%myWaterLevel(it, jt)
                    endif
                    
                    !diagonal is sqrt of squared distances

                    di = it - i
                    !distance to right cell
                    if (di > 0) then
                        dy = Me%ExtVar%DZY(i, j)
                    else
                       !distance to left cell
                        dy = Me%ExtVar%DZY(i-1, j)
                    endif

                    dj = jt - j
                    if (dj > 0) then
                        dx = Me%ExtVar%DZX(i, j)
                    else
                        dx = Me%ExtVar%DZX(i, j-1)
                    endif
                    
                    CenterDistance = sqrt((dx)**2 + (dy)**2)
                    
                    !Slope
                    if (Me%AdjustSlope) then
                        Slope           = AdjustSlope((level_up - level_down) / CenterDistance)
                    else
                        Slope           = (level_up - level_down) / CenterDistance
                    endif

                    if (Slope.LT.0.0) then
                        sign = -1.0
                    else
                        sign = 1.0
                    end if

                    AverageCellLengthSink   = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0  
                    AverageCellLengthLower  = ( Me%ExtVar%DUX (it, jt) + Me%ExtVar%DVY (it, jt) ) / 2.0              
                    VertArea                = ((AverageCellLengthSink + AverageCellLengthLower) / 2.0) * WaveHeight
                    
                    !Wet perimeter approximation to bottom (no walls effect)
                    !WetPerimeter    = (AverageCellLengthSink + AverageCellLengthLower) / 2.0
                    
                    !Same as wave height. short circuit
                    !HydraulicRadius = VertArea / WetPerimeter
                    HydraulicRadius = WaveHeight
                                 
                    OverlandCoef    = (AverageCellLengthSink * Me%OverlandCoefficient(i, j) +      &
                                       AverageCellLengthLower * Me%OverlandCoefficient(it, jt)) /  &
                                       (AverageCellLengthSink + AverageCellLengthLower)               
                    !
                    !MANNING'S EQUATION -  KINEMATIC WAVE
                    !
                    !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                    Flow = sign * VertArea * HydraulicRadius**(2./3.) * sqrt(sign * Slope)          &
                           / OverlandCoef

                    MaxFlow  = sign * VertArea * sqrt(Gravity * WaveHeight)
                    
                    if (abs(Flow) > abs(MaxFlow)) then
                        Flow = MaxFlow
                    endif
                    
                    Me%iFlowRouteDFour(i,j)    = Flow

                    !Updates Water Volume
                    Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - Flow *  Me%ExtVar%DT
                    Me%myWaterVolume (it, jt)  = Me%myWaterVolume (it, jt) + Flow *  Me%ExtVar%DT 

                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)
                    Me%myWaterColumn  (it, jt) = Me%myWaterVolume (it, jt) / Me%ExtVar%GridCellArea(it, jt)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    Me%myWaterLevel (it, jt)   = Me%myWaterColumn (it, jt) + Me%ExtVar%Topography(it, jt)
                
                else
                    Me%iFlowRouteDFour(i,j)    = 0.0
                endif
                
            endif
            
        enddo
        enddo
    

    end subroutine RouteDFourPoints

    !--------------------------------------------------------------------------
    
    subroutine StormWaterDrainage
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: FlowVolume, InfiltrationVolume
        real                                        :: dVol
        integer, dimension(:, :), pointer           :: DrainageDirection
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        integer                                     :: ilowest, jlowest
        real                                        :: dx, dy, dist, flow
        

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Gets Drainage Direction
        call GetDrainageDirection (Me%ObjBasinGeometry, DrainageDirection, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - StormWaterDrainage - ERR01'
        
        !Infiltrates water into the StormWater Drainage
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. &
                Me%myWaterColumn(i, j) > AllmostZero          .and. &
                Me%LowestNeighborI(i, j) /= null_int) then
            
                if (Me%StormWaterDrainageCoef(Me%LowestNeighborI(i, j), Me%LowestNeighborJ(i, j)) > AllmostZero) then

                    !Volume which can infiltrate at avaliable area during time step at max infiltration velocity
                    !m3 
                    FlowVolume                  = Me%StormWaterDrainageCoef(i, j) * Me%ExtVar%GridCellArea(i, j) * &
                                                  Me%StormWaterInfiltrationVelocity * Me%ExtVar%DT

                    !Volume which will be removed from overland flow into the StormWater System
                    !m3
                    InfiltrationVolume          = Min(Me%myWaterVolume(i, j), FlowVolume)

                    !New StormWater Volume at point
                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j) + InfiltrationVolume

                    !New Volume of Overland volume
                    Me%myWaterVolume (i, j)     = Me%myWaterVolume (i, j)   - InfiltrationVolume
                    
                    !New WaterColumn
                    Me%myWaterColumn (i, j)     = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                    !New Level
                    Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j) + Me%ExtVar%Topography(i, j)

                endif
        
            endif

        enddo
        enddo
        
        
        !Routes flows to the lowest Neighbour with constant velocity
        do j = JLB, JUB
        do i = ILB, IUB
            
            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%LowestNeighborI(i, j) /= null_int) then
            
                if (Me%LowestNeighborI(i, j) == i) then
                    dy = 0
                else
                    dy = Me%ExtVar%DZY(i, j)
                endif
                
                if (Me%LowestNeighborJ(i, j) == j) then
                    dx = 0
                else
                    dx = Me%ExtVar%DZX(i, j)
                endif
                
                dist = sqrt(dx**2.0 + dy**2.0)
                flow = Me%StormWaterVolume(i, j) / dist * Me%StormWaterFlowVelocity

                dVol = min(flow * Me%ExtVar%DT, Me%StormWaterVolume(i, j))
                
                flow = dVol/Me%ExtVar%DT
                
                !Output
                Me%StormWaterCenterModulus (i, j) = flow
                if      (Me%LowestNeighborI(i, j) == i-1 .and. Me%LowestNeighborJ(i, j) == j-1) then
                    Me%StormWaterCenterFlowX   (i, j) = -1.0 * sqrt(flow)
                    Me%StormWaterCenterFlowY   (i, j) = -1.0 * sqrt(flow)
                else if (Me%LowestNeighborI(i, j) == i-1 .and. Me%LowestNeighborJ(i, j) == j) then
                    Me%StormWaterCenterFlowX   (i, j) = 0.0
                    Me%StormWaterCenterFlowY   (i, j) = -1.0 * flow
                else if (Me%LowestNeighborI(i, j) == i-1 .and. Me%LowestNeighborJ(i, j) == j+1) then
                    Me%StormWaterCenterFlowX   (i, j) = +1.0 * sqrt(flow)
                    Me%StormWaterCenterFlowY   (i, j) = -1.0 * sqrt(flow)
                else if (Me%LowestNeighborI(i, j) == i .and. Me%LowestNeighborJ(i, j) == j+1) then
                    Me%StormWaterCenterFlowX   (i, j) = flow
                    Me%StormWaterCenterFlowY   (i, j) = 0.0
                else if (Me%LowestNeighborI(i, j) == i+1 .and. Me%LowestNeighborJ(i, j) == j+1) then
                    Me%StormWaterCenterFlowX   (i, j) = sqrt(flow)
                    Me%StormWaterCenterFlowY   (i, j) = sqrt(flow)
                else if (Me%LowestNeighborI(i, j) == i+1 .and. Me%LowestNeighborJ(i, j) == j) then
                    Me%StormWaterCenterFlowX   (i, j) = 0.0
                    Me%StormWaterCenterFlowY   (i, j) = flow
                else if (Me%LowestNeighborI(i, j) == i+1 .and. Me%LowestNeighborJ(i, j) == j-1) then
                    Me%StormWaterCenterFlowX   (i, j) = -1.0 * sqrt(flow)
                    Me%StormWaterCenterFlowY   (i, j) = +1.0 * sqrt(flow)
                else
                    Me%StormWaterCenterFlowX   (i, j) = -1.0 * flow                    
                    Me%StormWaterCenterFlowY   (i, j) = 0.0
                endif                    
                                            
                
                !
                !If the lowest neighbor is a stromwater drainage point, route it there. Otherwise put the water back to the surface
                !                    
                if (Me%StormWaterDrainageCoef(Me%LowestNeighborI(i, j), Me%LowestNeighborJ(i, j)) > AllmostZero) then
                    ilowest = Me%LowestNeighborI(i, j)
                    jlowest = Me%LowestNeighborJ(i, j)
                    Me%StormWaterVolume(ilowest, jlowest) = Me%StormWaterVolume(ilowest, jlowest)   + dVol
                    Me%StormWaterVolume(i,       j      ) = Me%StormWaterVolume(i,       j      )   - dVol
                else
                    
                    !New Volume of Overland volume
                    Me%myWaterVolume   (i, j) = Me%myWaterVolume   (i, j) + dVol
                    
                    !New WaterColumn
                    Me%myWaterColumn   (i, j) = Me%myWaterVolume   (i, j) / Me%ExtVar%GridCellArea(i, j)

                    !New Level
                    Me%myWaterLevel    (i, j) = Me%myWaterColumn   (i, j) + Me%ExtVar%Topography(i, j)
                    
                    Me%StormWaterVolume(i, j) = Me%StormWaterVolume(i, j) - dVol
                endif
                
            else
            
                Me%StormWaterCenterModulus (i, j) = 0.0
                Me%StormWaterCenterFlowX   (i, j) = 0.0
                Me%StormWaterCenterFlowY   (i, j) = 0.0
            
            endif
               
        enddo
        enddo
        
!        
!        !Flow along X
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i, j-1) > AllmostZero) then
!            
!                !Flow from the left to the right
!                if (Me%ExtVar%Topography(i, j-1) > Me%ExtVar%Topography(i, j)) then
!                    Me%StormWaterFlowX(i, j) =      Me%StormWaterVolume(i, j-1) / Me%ExtVar%DZX(i, j-1) * Me%StormWaterFlowVelocity
!                else
!                    Me%StormWaterFlowX(i, j) = -1 * Me%StormWaterVolume(i, j)   / Me%ExtVar%DZX(i, j-1) * Me%StormWaterFlowVelocity
!                endif
!
!            else
!            
!                Me%StormWaterFlowX(i, j) = 0
!            
!            endif
!            
!            
!            
!        enddo
!        enddo
!        
!
!        !Flow along Y
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i-1, j) > AllmostZero) then
!            
!                if (Me%ExtVar%Topography(i-1, j) > Me%ExtVar%Topography(i, j)) then
!                    Me%StormWaterFlowY(i, j) =      Me%StormWaterVolume(i-1, j) / Me%ExtVar%DZY(i-1, j) * Me%StormWaterFlowVelocity
!                else
!                    Me%StormWaterFlowY(i, j) = -1 * Me%StormWaterVolume(i, j)   / Me%ExtVar%DZX(i-1, j) * Me%StormWaterFlowVelocity
!                endif
!            
!            else
!            
!                Me%StormWaterFlowY(i, j) = 0
!            
!            endif
!            
!        enddo
!        enddo
!            
!        !Updates volumes
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i, j-1) > AllmostZero) then
!                
!                if (Me%StormWaterFlowX(i, j) > 0) then
!                
!                    dVol = min(Me%StormWaterFlowX(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i, j-1))
!                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   + dVol
!                    Me%StormWaterVolume(i, j-1) = Me%StormWaterVolume(i, j-1) - dVol
!                    
!                else
!                
!                    dVol = min(-Me%StormWaterFlowX(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i, j))
!                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   - dVol
!                    Me%StormWaterVolume(i, j-1) = Me%StormWaterVolume(i, j-1) + dVol
!                
!                endif
!                
!            endif
!        enddo
!        enddo
!        
!        
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i-1, j) > AllmostZero) then
!                
!                if (Me%StormWaterFlowY(i, j) > 0) then
!                
!                    dVol = min(Me%StormWaterFlowY(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i-1, j))
!                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   + dVol
!                    Me%StormWaterVolume(i-1, j) = Me%StormWaterVolume(i-1, j) - dVol
!                    
!                else
!                
!                    dVol = min(-Me%StormWaterFlowY(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i, j))
!                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   - dVol
!                    Me%StormWaterVolume(i-1, j) = Me%StormWaterVolume(i-1, j) + dVol
!                
!                endif
!            endif
!        enddo
!        enddo
!                    
      
        !Routes water from StormWater Drainage System to river channels
        if (Me%ObjDrainageNetwork /= 0) then 

            call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StormWaterDrainage - ModuleRunOff - ERR01'     

            do j = JLB, JUB
            do i = ILB, IUB
            
                if (Me%ExtVar%RiverPoints(i,j) == BasinPoint) then

                    if (ChannelsWaterLevel (i, j) < Me%myWaterLevel(i, j)) then
                    
                        Me%iFlowToChannels(i, j)  = Me%iFlowToChannels(i, j) + Me%StormWaterVolume(i, j) / Me%ExtVar%DT
                        Me%StormWaterVolume(i, j) = 0.0
                    
                    endif

                endif

            enddo
            enddo
            
            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StormWaterDrainage - ModuleRunOff - ERR05'
            
        endif

        
        call UnGetBasin (Me%ObjBasinGeometry, DrainageDirection, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - StormWaterDrainage - ERR02'
                
        
    end subroutine StormWaterDrainage

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine StreetGutterFlow
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: flow
        real                                        :: AverageCellLength, y0
        integer                                     :: targetI, targetJ
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Calculates the inflow at each street gutter point
        !$OMP PARALLEL PRIVATE(I,J, flow, AverageCellLength)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%myWaterColumn(i, j) > Me%MinimumWaterColumn .and. Me%StreetGutterLength(i, j) > AllmostZero) then 
        
                !Q  = L * K * y0^(3/2) * sqrt(g)
                !L  = Cumprimento Sargeta = 0.5
                !K  = Coef = 0.2
                !y0 = Altura a montante da sargeta
                       
                AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0
                
                !Considering an average side slope of 5% (1/0.05 = 20) of the street
                y0 = sqrt(2.0*Me%myWaterColumn(i, j)*AverageCellLength / 20.0)
                
                !When triangle of street is full, consider new head 
                if (y0 * 20.0 > AverageCellLength) then
                    y0 = AverageCellLength / 40.0 + Me%myWaterColumn(i, j)
                endif
                
                !flow = 0.5 * 0.2 * y0**1.5 * sqrt(Gravity)
                flow = Me%StreetGutterLength(i, j) * 0.2 * y0**1.5 * sqrt(Gravity)
            
                !Flow Rate into street Gutter
                Me%StreetGutterFlow(i, j) = Min(flow, Me%myWaterVolume(i, j) / Me%ExtVar%DT)
                
            else
            
                Me%StreetGutterFlow(i, j) = 0.0
            
            endif

        enddo
        enddo
        !$OMP END DO
        
        !Integrates values from gutter flow at sewer points
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%StormWaterInteraction(i, j) > AllmostZero) then
                Me%SewerInflow(i, j) = 0.0
            endif
        enddo
        enddo    
        !$OMP END DO
        !$OMP END PARALLEL
                    
        !This do loop should not be parallel
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%StreetGutterLength(i, j) > AllmostZero) then
            
                !SEWER interaction point
                targetI = Me%StreetGutterTargetI(i, j)
                targetJ = Me%StreetGutterTargetJ(i, j)
                
                Me%SewerInflow(targetI, targetJ) = Me%SewerInflow(targetI, targetJ) + Me%StreetGutterFlow(i, j)  
            endif
        enddo
        enddo                
        
        
    end subroutine StreetGutterFlow
        
    !--------------------------------------------------------------------------
        
    subroutine AddFlowFromStormWaterModel

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: targetI, targetJ
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J, targetI, targetJ)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        
        !The algorithm below has the following assumptions
        !1. MOHID Land calculates the possible inflow into the sewer system through the street gutters (matrix StreetGutterFlow)
        !2. The values of the StreetGutterFlow are integrated at the nearest StormWaterInteraction points (matrix SewerInflow)
        !3. This matrix (SewerInflow) is provide to SWMM
        !4. Swmm calculates the efective inflow and returns the efective flow at each interaction point (matrix StormWaterModelFlow)
        !5. The algorithm below calculates the efective inflow in each Street Gutter Point
        !5a  - if the flow in the target point is negative (inflow into the sewer system) the flow at each point is the relative one 
        !      StreetGutterFlow/StormWaterModelFlow
        !5b  - if the flow in the target point is positive (outflow from the sewer system), the flow flows out at the target point 
        !      ("saltam as tampas"). 
        !6. The Water Column is reduces due to the final flow


        !Algorithm which calculates the real inflow in each point
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                Me%StormInteractionFlow(i, j) = 0.0

                !If the point is a street gutter point
                !we have to reduce the volume by the total number of associated inlets
                if (Me%StreetGutterLength(i, j) > 0.0) then 
               
                    targetI = Me%StreetGutterTargetI(i, j)
                    targetJ = Me%StreetGutterTargetJ(i, j)

                    if (Me%StormWaterModelFlow(targetI, targetJ) < 0.0 .and. Me%SewerInflow(targetI, targetJ) > AllmostZero) then
                        !Distribute real / potentil
                        Me%StormInteractionFlow(i, j) = -1.0 * Me%StreetGutterFlow(i, j) * &
                                                        Me%StormWaterModelFlow(targetI, targetJ) / Me%SewerInflow(targetI, targetJ)
                                                               
                    endif
                    
                endif
                
                !Overflow of the sewer system
                if (Me%StormWaterInteraction(i, j) > AllmostZero .and. Me%StormWaterModelFlow(i, j) > 0) then
                    Me%StormInteractionFlow(i, j) = -1.0 * Me%StormWaterModelFlow(i, j)
                endif
                
            endif

        enddo
        enddo  
        !$OMP END DO

        !Update water column
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                Me%myWaterColumnOld   (i, j)  = Me%myWaterColumnOld   (i, j)  -         &
                                                Me%StormInteractionFlow(i, j) *         &
                                                Me%ExtVar%DT /                          &
                                                Me%ExtVar%GridCellArea(i, j)
                                             
                if (Me%myWaterColumnOld(i, j) < 0.0) then
                    Me%MassError     (i, j) = Me%MassError(i, j) - Me%myWaterColumnOld(i, j) * &
                                              Me%ExtVar%GridCellArea(i, j)
                    
                    Me%myWaterColumnOld (i, j)  = 0.0
                endif
                
            endif

        enddo
        enddo  
        !$OMP END DO
        !$OMP END PARALLEL
                  
    
    end subroutine AddFlowFromStormWaterModel
    
    !--------------------------------------------------------------------------
    
    subroutine FlowIntoChannels(LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: DifLevel
        real                                        :: Slope, AverageCellLength, dVol
        real                                        :: Area, HydraulicRadius, MaxFlow
        real                                        :: ChannelFreeVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength 
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsVolume
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR02'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR03'        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'   


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !$OMP PARALLEL PRIVATE(I,J, DifLevel, Slope, AverageCellLength, dVol, Area, HydraulicRadius, MaxFlow, ChannelFreeVolume)


        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. ChannelsActiveState(i, j) == BasinPoint) then

                !Checks for Flow from Land -> Channel
                AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0

            
                if (ChannelsWaterLevel (i, j) < Me%myWaterLevel(i, j) .and. Me%myWaterColumn(i, j) > Me%MinimumWaterColumn) then

                    if (ChannelsWaterLevel (i, j) > Me%ExtVar%Topography(i, j)) then
                        DifLevel           = Me%myWaterLevel(i, j) - ChannelsWaterLevel (i, j)
                    else
                        DifLevel           = Me%myWaterColumn(i, j)
                    endif

                    !Volume which can enter the channel
                    ChannelFreeVolume = ChannelsMaxVolume(i, j) - ChannelsVolume (i, j)
                
                    !Channel almost empty... put all water into channel    
!                    if (ChannelFreeVolume / ChannelsMaxVolume(i, j) > 0.01) then

                        !Volume to channel: minimum between free volume and current volume in cell
!                        dVol = min(ChannelFreeVolume, Me%myWaterVolume (i, j))

                        !Flow to channel - positive if enters
!                        Me%lFlowToChannels(i, j) = dVol / LocalDT

!                    else
                
                        Slope                      = AdjustSlope(DifLevel / (AverageCellLength / 4.0))

                        Area                       = DifLevel * ChannelsNodeLength(i, j)
                    
                        HydraulicRadius            = Area / ChannelsNodeLength(i, j)
                
                        !Minium between friction (manning) and critical flow
                        Me%lFlowToChannels(i, j)   = min(Area * HydraulicRadius**(2./3.) * sqrt(Slope) /  &
                                                         Me%OverlandCoefficient(i,j), &
                                                         Area * sqrt(Gravity * DifLevel))
                        
                     
                        !MaxFlow = 0.5 * (DifLevel) * Me%ExtVar%GridCellArea(i, j) / LocalDT

                        MaxFlow = sqrt(Gravity * Me%myWaterColumn(i, j)) * Me%myWaterColumn(i, j) * ChannelsNodeLength(i, j)
                   
                        if (Me%lFlowToChannels(i, j) > MaxFlow) then
                            Me%lFlowToChannels(i, j) = MaxFlow
                        endif
                        
!                    endif
                    
                    
                              
                    !dVol
                    dVol                       = Me%lFlowToChannels(i, j) * LocalDT
                    
                    !Updates Water Volume
                    Me%myWaterVolume (i, j)    = Me%myWaterVolume (i, j)   - dVol 
                    
                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    

                else
                
                    Me%lFlowToChannels(i, j) = 0.0
                
                endif

            
            endif

        enddo
        enddo        
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR10'
        
        
    end subroutine FlowIntoChannels   
    
    !--------------------------------------------------------------------------

    subroutine FlowFromChannels
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: ChannelHeight
        real                                        :: WCR, dVol, VolExcess, NewLevel
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real   , dimension(:, :), pointer           :: ChannelsSurfaceWidth
        real   , dimension(:, :), pointer           :: ChannelsBankSlope
        real   , dimension(:, :), pointer           :: ChannelsBottomLevel
        real                                        :: a0, a1, a2
        real                                        :: x1, x2, MaxFlow, Flow
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR02'

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR03'

        call GetChannelsBankSlope (Me%ObjDrainageNetwork, ChannelsBankSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR04'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR05'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'        


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !$OMP PARALLEL PRIVATE(I,J, ChannelHeight, WCR, dVol, VolExcess, NewLevel, a0, a1, a2, x1, x2, MaxFlow)


        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. ChannelsActiveState(i, j) == BasinPoint) then

                if (ChannelsWaterLevel (i, j) > Me%myWaterLevel(i, j)) then
                
                    ChannelHeight = Me%ExtVar%Topography(i, j) - ChannelsBottomLevel(i, j)                                       
                    !ChannelSlope  = (ChannelsTopWidth(i, j) - ChannelsBottomWidth(i, j)) / ChannelHeight
                    !ChannelSurfaceWidth = ChannelsBottomWidth(i,j) + 2.* ChannelSlope * ChannelHeight
                    
                    !Water Column in River above Topo
                    WCR           = ChannelsWaterLevel (i, j) - Me%ExtVar%Topography(i, j)
                    
                    !Volume above Topography
                    VolExcess    = ChannelsBankSlope(i,j) * WCR * WCR * ChannelsNodeLength(i, j)       &
                                    + WCR * ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) +    &
                                    Me%myWaterVolume(i, j)

                    if (ChannelsBankSlope(i,j) <= AlmostZero) then
                        !Rectangular
                        a1 = ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) + Me%ExtVar%GridCellArea(i, j)
                        NewLevel = VolExcess / a1
                        NewLevel = NewLevel + Me%ExtVar%Topography(i, j)

                    else
                        !Trapezoidal - formula resolvente
                        a0 = ChannelsBankSlope(i,j) * ChannelsNodeLength(i, j)
                        a1 = ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) + Me%ExtVar%GridCellArea(i, j)
                        a2 = -1.0 * VolExcess
                                    
                        !Solves Polynominal
                        x1            = (-a1 + sqrt(a1**2. - 4.*a0*a2)) / (2.*a0)
                        x2            = (-a1 - sqrt(a1**2. - 4.*a0*a2)) / (2.*a0)                        

                        if (x1 > 0. .and. x1 < WCR) then
                            NewLevel  = x1 + Me%ExtVar%Topography(i, j)
                        else
                            NewLevel  = x2 + Me%ExtVar%Topography(i, j)
                        endif
                    endif

                    
                    dVol = (NewLevel - Me%myWaterLevel(i, j)) *  Me%ExtVar%GridCellArea(i, j)
                    
!                    Me%iFlowToChannels(i, j)    = -dVol / Me%ExtVar%DT 
                    !Revision David 10/4/10
                    !Usually for each cell flow has only one direction
                    !But may exist the special case where at the beggining channel level is lower than
                    !runoff level, but with the exchange, the channel level got bigger
                    !and a flow addition (subtraction) is needed    
                    !Me%iFlowToChannels(i, j)    = Me%iFlowToChannels(i, j) -dVol / Me%ExtVar%DT     
                    Flow = -dVol / Me%ExtVar%DT
                    
                    !Limits flow to critical one
                    MaxFlow = -1.0 * sqrt(Gravity * WCR) * WCR * ChannelsNodeLength(i, j)
                    
                    if (Flow > MaxFlow) then
                        Flow = MaxFlow
                    endif
                    
                    Me%iFlowToChannels(i, j)    = Me%iFlowToChannels(i, j) + Flow
                    
                    Me%myWaterVolume (i, j)     = Me%myWaterVolume (i, j) - (Flow *  Me%ExtVar%DT)
                    
                    Me%myWaterColumn  (i, j)    = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)

                
                endif

            
            endif

        enddo
        enddo        
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBankSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        
   
    
    end subroutine FlowFromChannels
    
    !--------------------------------------------------------------------------
    
    subroutine OverLandChannelInteraction
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: dVol, Flow, a1
        real                                        :: TotalVolume, VolExcess, NewLevel
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real   , dimension(:, :), pointer           :: ChannelsSurfaceWidth
        real   , dimension(:, :), pointer           :: ChannelsBankSlope
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR02'

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR03'

        call GetChannelsBankSlope (Me%ObjDrainageNetwork, ChannelsBankSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR04'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'        


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet

                !Total volume in the cell + channel
                TotalVolume = Me%myWaterVolume(i, j) + ChannelsVolume (i, j)
                
                !All water fits into channel?
                if (TotalVolume < ChannelsMaxVolume(i, j)) then
                
                    !Total volume fits into channel. 
                    !Route Flow from overland to channel, using free flow condition
                    !Maximum flow is equal to volume avaliable at watercolumn
                    
                    !Free drop, dh will be the water column
!                    dh      = Me%MyWaterColumn(i,j)
                    !dh      = abs(Me%myWaterLevel(i, j) - ChannelsWaterLevel(i, j))
                    
 !                   if (dh > Me%MinimumWaterColumn) then
!
!                        FlowCel = sqrt(Gravity * dh) * 2.0 * ChannelsNodeLength(i, j) * dh
!                        
!                        MaxFlow = Me%myWaterVolume (i, j) / Me%ExtVar%DT
!                        
!                        Me%iFlowToChannels(i, j) = min(Flow, MaxFlow)
                        
                        !in terms of velocity water already arrived to runoff center cell so it should be
                        !in the river (instantaneously)
                        !!!!Me%iFlowToChannels(i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%DT
                        Flow = Me%myWaterVolume (i, j) / Me%ExtVar%DT
                        
!                        !limit to critical
!                        !m3/s = celerity (m/s) * m2 (Area = (dh * L) * 2)
!                        MaxFlow    = sqrt(Gravity * dh) * 2.0 * ChannelsNodeLength(i, j) * dh
!                                              
!                        if (Me%iFlowToChannels(i, j) > MaxFlow) then
!                            Me%iFlowToChannels(i, j) = MaxFlow
!                        endif
!                    else
!                    
!                        Me%iFlowToChannels(i, j) = 0.0
!                        
!                    endif
                    
                else
                

                    !Total Volume does not fit into the channel
                    !Route flow in a way the the water level becomes horizontal

                    !Limit flow so that volumes to not become negative and critical flow is not exceeded
                    
                    !Volume which does not fit in the channel
                    VolExcess = TotalVolume - ChannelsMaxVolume(i, j)
                    !            
                    !if (ChannelsBankSlope(i,j) <= AlmostZero) then !Rectangular channel
                        
                        !Total area -> channel area + grid cell
                        a1 = ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) + Me%ExtVar%GridCellArea(i, j)
                        
                        !New level = ExcessVolume / area + ground level
                        NewLevel = VolExcess / a1 + Me%ExtVar%Topography(i, j)
                        
                    !else                                            !Trapezoidal - formula resolvente
                    !    
                    !    !VolExcess = newLevel * GridCellAreas + ChannelsNodeLength * (TopWidth + TopWidth * 2 * BankSlope*h) * h) / 2
                    !
                    !    a0 = ChannelsBankSlope(i,j) * ChannelsNodeLength(i, j)
                    !    a1 = ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) + Me%ExtVar%GridCellArea(i, j)
                    !    a2 = -1.0 * VolExcess
                    !
                    !    !Solves Polynominal
                    !    x1            = (-a1 + sqrt(a1**2. - 4.*a0*a2)) / (2.*a0)
                    !    x2            = (-a1 - sqrt(a1**2. - 4.*a0*a2)) / (2.*a0)                        
                    !
                    !    if (x1 > 0. .and. x1 < ChannelsWaterLevel (i, j) - Me%ExtVar%Topography(i, j)) then
                    !        NewLevel  = x1 + Me%ExtVar%Topography(i, j)
                    !    else
                    !        NewLevel  = x2 + Me%ExtVar%Topography(i, j)
                    !    endif
                    !endif
                    
                    !dh will be the maximum height above topography (lateral moving water)
!                    dh = max (Me%myWaterLevel(i, j), ChannelsWaterLevel(i, j)) - Me%ExtVar%Topography(i,j)
                    !dh      = abs(Me%myWaterLevel(i, j) - ChannelsWaterLevel(i, j))

!                    if (dh > Me%MinimumWaterColumn) then
                    
                        !Variation in volume by comparing old level with new level
                        dVol = (Me%myWaterLevel(i, j) - NewLevel ) *  Me%ExtVar%GridCellArea(i, j)
             
                        !Me%iFlowToChannels(i, j)    = dVol / Me%ExtVar%DT
                        Flow      = dVol / Me%ExtVar%DT
                        
                        !Prevent negative volumes
                        if (Flow > 0.0) then
                            Flow = min(Flow, Me%myWaterVolume(i, j) / Me%ExtVar%DT)
                        else
                            Flow = max(Flow, -1.0 * (ChannelsVolume(i,j) - ChannelsMaxVolume(i, j)) / Me%ExtVar%DT)
                        endif


!                        !in terms of velocity water already arrived to runoff/DN center cell so it should be
!                        !in the river or runoff (instantaneously)
!                        !Limits to critical flow to critical one
!                        MaxFlow = sqrt(Gravity * dh) * 2.0 * ChannelsNodeLength(i, j) * dh
!                        if (abs(Flow) > MaxFlow) then
!                            if (Flow > 0.0) then
!                                Flow = MaxFlow
!                            else
!                                Flow = -1.0 * MaxFlow
!                            endif
!                        endif
!                        
!                    else
!                    
!                        Flow = 0.0
!                        
!                    endif                    
                    
                endif

                !!Limits the change to a constant value. Only for test purposes
                !Flow = 0.1 * Flow
                
                !!Important!! flow to channel may have other sources than this, so a sum is needed
                Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Flow

                !Updates Volumes
                !Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - Me%iFlowToChannels    (i, j) * Me%ExtVar%DT
                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                
                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)
                           
            endif

        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBankSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR09'


        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

    
    end subroutine OverLandChannelInteraction
    
    !--------------------------------------------------------------------------
    !Method to instantaneously transport water in river runoff interaction (stability problems)
    !this method should be deleted

    subroutine OverLandChannelInteraction_1
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: dVol
        real                                        :: TotalVolume, VolExcess
        real                                        :: NewH
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real   , dimension(:, :), pointer           :: ChannelsSurfaceWidth
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR05'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR02'

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR03'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR06'        


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet

                !Total volume in the cell + channel
                TotalVolume = Me%myWaterVolume(i, j) + ChannelsVolume (i, j)
                
                !All water fits into channel?
                if (TotalVolume < ChannelsMaxVolume(i, j)) then
                
                    Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Me%myWaterVolume (i, j) / Me%ExtVar%DT
                    
                    Me%myWaterVolume (i, j) = 0.0
                    Me%myWaterColumn (i, j) = 0.0
                    Me%myWaterLevel  (i, j) = Me%ExtVar%Topography  (i, j)

                else               

                    !Total Volume does not fit into the channel
                    !Route flow in a way the the water level becomes horizontal
                    !Limit flow so that volumes to not become negative and critical flow is not exceeded
                    
                    !Volume which does not fit into the channel
                    VolExcess = TotalVolume - ChannelsMaxVolume(i, j)
                    
                    !New Height of water in cell
                    NewH = VolExcess / Me%ExtVar%GridCellArea(i, j)
                    
                    !Flow to or from river is calculated based on the level difference (new to old)
                    !If the new level is higher than the old one, the flow will be positive (flow to channel) 
                    !m3   = (m + m - m) * (m * m)
                    dVol = (NewH + Me%ExtVar%Topography(i, j)    &
                           - ChannelsWaterLevel(i, j)) * (ChannelsNodeLength(i, j) * ChannelsSurfaceWidth(i, j))
                    Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + dVol / Me%ExtVar%DT
                                        
                    !Updates Volumes            
                    Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - dVol                
                    Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                    Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)
                endif
                           
            endif

        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR06'

        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_1 - ModuleRunOff - ERR010'        

    
    end subroutine OverLandChannelInteraction_1
    
    !--------------------------------------------------------------------------
    
    !Method to use celerity as the base for transport water in river runoff interaction
    subroutine OverLandChannelInteraction_2
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: Flow, MaxFlow
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real                                        :: dh, dh_new, WaveHeight, Celerity
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        real  , dimension(:, :), pointer            :: ChannelsSurfaceWidth
        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR02'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'        

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR03'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet


                !dh > 0, flow to channels, dh < 0, flow from channels
                dh         =  Me%myWaterLevel(i, j) - ChannelsWaterLevel(i, j)
                WaveHeight =  max(Me%myWaterLevel(i, j), ChannelsWaterLevel(i, j)) - Me%ExtVar%Topography(i,j)

                Celerity = sqrt(Gravity * WaveHeight)
                
                if (dh > 0) then
                
                    if (Me%myWaterColumn (i, j) > Me%MinimumWaterColumn) then
                    
                        !flux is occuring between dh and with celerity 
                        !m3/s = m/s (celerity) * m2 (Area = (dh * L) * 2)
                        Flow    = Celerity * 2.0 * ChannelsNodeLength(i, j) * min(dh, WaveHeight)
                        
                        !MaxFlow = Me%myWaterVolume (i, j) / Me%ExtVar%DT
                        !if channel level lower than topography - limit is all volume (waveheight is water column)
                        !if channel level higher than topography limit is dh
                        Maxflow = min(dh, WaveHeight) * Me%ExtVar%GridCellArea(i,j) / Me%ExtVar%DT
                    else
                    
                        Flow = 0.0
                        MaxFlow = 0.0
                    
                    endif
                else
                    !Implicit computation of new dh based on celerity dx transport
!                    dh_new = (ChannelsSurfaceWidth(i,j) * dh) /                      &
!                    (ChannelsSurfaceWidth(i,j) + 2 * min (Celerity * Me%ExtVar%DT, 0.5 * Me%ExtVar%DUX(i,j)))
                    
                    !Compute new water height above runoff column based on the distance that water 
                    !will be spread in one dt (surface width + 2 celerity paths - in both ways)
                    ![m] = [m] * [m] / [m] . this is the same as working with volumes where river lenght
                    !would be multiplied in both num and den. dh_new is estimated based on same volume spreading on
                    !wider area
                    dh_new = (ChannelsSurfaceWidth(i,j) * dh) /                      &
                    (ChannelsSurfaceWidth(i,j) + 2 * Celerity * Me%ExtVar%DT)
                    
                    !maximum spread where in one time step all the water above runoff column
                    !will spread along all the cell (DUX)
                    !in case that channel top width is == DUX no flow occurs so this was abandoned
                    !dh_min = (ChannelsSurfaceWidth(i,j) * dh) /                      &
                    !(Me%ExtVar%DUX(i,j))
                    
                    !m3/s = h * L * Length / s
                    Flow    = -1. * (dh_new - dh) * ChannelsSurfaceWidth(i,j) * ChannelsNodeLength(i,j) / Me%ExtVar%DT
                    
                    !MaxFlow = -1. * (dh_min - dh) * ChannelsSurfaceWidth(i,j) * ChannelsNodeLength(i,j) / Me%ExtVar%DT                    
                    !maximum is the channel water above runoff going all to runoff 
                    MaxFlow = dh * ChannelsSurfaceWidth(i,j) * ChannelsNodeLength(i,j) / Me%ExtVar%DT    
                endif

                if (abs(Flow) > abs(MaxFlow)) then
                    Flow = MaxFlow
                endif   
                    
                !!Important!! flow to channel may have other sources than this, so a sum is needed
                Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Flow

                !Updates Volumes
                !Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - Me%iFlowToChannels    (i, j) * Me%ExtVar%DT
                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                
                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)
                           
            endif

        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR020'        

    
    end subroutine OverLandChannelInteraction_2
    
    !--------------------------------------------------------------------------    
    !Method uses complicated method for transport water in river runoff interaction
    !instantaneous if space available in river and a ixture of celerity and instant
    !equalization otherwise. this method should be deleted is unstable as 1 and complex    
    subroutine OverLandChannelInteraction_3
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: Flow, MaxFlow
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real   , dimension(:, :), pointer           :: ChannelsTopArea
        real                                        :: dh, dh_new, WaveHeight, Celerity, dVol
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        real  , dimension(:, :), pointer            :: ChannelsSurfaceWidth
        real                                        :: CellLevel, TotalVolume, VolExcess, NewH
        real                                        :: NewHOnCell, NewHOnRiver
        real                                        :: NewLevel
        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR010'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR020'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR030'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR040'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR050'        

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR060'
        
        call GetChannelsTopArea     (Me%ObjDrainageNetwork, ChannelsTopArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR070'        


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet

                !Total volume in the cell + channel
                TotalVolume = Me%myWaterVolume(i, j) + ChannelsVolume (i, j)
                
                !All water fits into channel?
                if (TotalVolume < ChannelsMaxVolume(i, j)) then
                
                    Flow = Me%myWaterVolume (i, j) / Me%ExtVar%DT
                    
                    Me%myWaterVolume (i, j) = 0.0
                    Me%myWaterColumn (i, j) = 0.0
                    Me%myWaterLevel  (i, j) = Me%ExtVar%Topography  (i, j)                
            
                else
                
                    !Total Volume does not fit into the channel
                    !Route flow in a way the the water level becomes horizontal
                    !Limit flow so that volumes to not become negative and critical flow is not exceeded
                    
                    !Volume which does not fit into the channel
                    VolExcess = TotalVolume - ChannelsMaxVolume(i, j)
                    
                    !New Height of water in cell
                    NewH = VolExcess / Me%ExtVar%GridCellArea(i, j) 
                    NewLevel = NewH + Me%ExtVar%Topography(i, j)
                    
                    if (NewLevel > ChannelsWaterLevel(i,j)) then
                        !There is more water on Runoff than on Channel
                        
                        !Flow to river is calculated based on the level difference (new to old) 
                        !Flow will be positive because there is more water on Runoff than on channel
                        dVol = (NewLevel - ChannelsWaterLevel(i, j)) * ChannelsTopArea(i, j)
                        Flow = dVol / Me%ExtVar%DT
                                        
                        !Updates Volumes            
                        Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - dVol                
                        Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                        Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)  
                        
                    else
                        !Water level on cell is defined by the volume of water on cell divided by the area of cell minus the area of the 
                        !channel plus the topography
                        CellLevel = (Me%myWaterVolume (i, j) / (Me%ExtVar%GridCellArea(i, j) &
                                     - ChannelsTopArea(i, j))) + Me%ExtVar%Topography(i, j)
                
                        !dh > 0, flow to channels, dh < 0, flow from channels
                        dh         =  CellLevel - ChannelsWaterLevel(i, j)
                        WaveHeight =  max(CellLevel, ChannelsWaterLevel(i, j)) - Me%ExtVar%Topography(i,j)

                        Celerity = sqrt(Gravity * WaveHeight)
                
                        !Implicit computation of new dh based on celerity dx transport
                        dh_new = (ChannelsSurfaceWidth(i,j) * dh) /                      &
                                 (ChannelsSurfaceWidth(i,j) + 2 * min (Celerity * Me%ExtVar%DT, 0.5 * Me%ExtVar%DUX(i,j)))
                    
                        !m3/s = h * L * Length / s
                        Flow    = -1. * (dh_new - dh) * ChannelsTopArea(i,j) / Me%ExtVar%DT
                        MaxFlow = -1. * (ChannelsVolume(i,j) - ChannelsMaxVolume(i,j)) / Me%ExtVar%DT
                        if (abs(Flow) > abs(MaxFlow)) then
                            !This is necessary so we can use Height instead of Level on the next pair of equations below
                            Flow = MaxFlow 
                        endif 
                                                
                        NewHOnRiver = (ChannelsVolume(i,j) - ChannelsMaxVolume(i,j)      &
                                       + (Flow * Me%ExtVar%DT)) / ChannelsTopArea(i,j)    
                        NewHOnCell  = (Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT))  &
                                       / (Me%ExtVar%GridCellArea(i, j) - ChannelsTopArea(i,j))
                        
                        if (NewHOnRiver < NewHOnCell) then
                            !If this is true, this means that the celerity will (maybe) make the river and runoff unstable.
                            !So, celerity will not be used and the levels will be harmonized 
                        
                            !Volume which does not fit into the channel
                            VolExcess = TotalVolume - ChannelsMaxVolume(i, j)
                    
                            !New Height of water in cell
                            NewH = VolExcess / Me%ExtVar%GridCellArea(i, j) 
                            NewLevel = NewH + Me%ExtVar%Topography(i, j) 
                            
                            !Flow to or from river is calculated based on the level difference (new to old)
                            !If the new level is higher than the old one, the flow will be positive (flow to channel) 
                            dVol = (NewLevel - ChannelsWaterLevel(i, j)) * ChannelsTopArea(i, j)
                            Flow = dVol / Me%ExtVar%DT
                                        
                            !Updates Volumes            
                            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - dVol                
                            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j) 
                        else
                            !Updates Volumes                            
                            Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)                
                            Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                            Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j) 
                        endif                                               
                    endif                
                endif        
            endif
                    
            !!Important!! flow to channel may have other sources than this, so a sum is needed
            Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Flow            
        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR020'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsTopArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR020'         
    
    end subroutine OverLandChannelInteraction_3
    
    !--------------------------------------------------------------------------     
    !Method to use celerity as the base for transport water in river runoff interaction
    !Difference to method 2 is the max flow definition in case water outside section    
    subroutine OverLandChannelInteraction_4
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: Flow, MaxFlow
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real   , dimension(:, :), pointer           :: ChannelsTopArea
        real                                        :: dh, dh_new, WaveHeight, Celerity, dVol
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        real  , dimension(:, :), pointer            :: ChannelsSurfaceWidth
        real                                        :: TotalVolume, VolExcess, NewH
        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR010'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR020'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR030'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR040'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR050'        

        call GetChannelsSurfaceWidth (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR060'
        
        call GetChannelsTopArea     (Me%ObjDrainageNetwork, ChannelsTopArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR070'        


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet

                !Total volume in the cell + channel
                TotalVolume = Me%myWaterVolume(i, j) + ChannelsVolume (i, j)  
                
                !All water fits into channel?
                if (TotalVolume < ChannelsMaxVolume(i, j)) then                
                    MaxFlow = Me%myWaterVolume (i, j) / Me%ExtVar%DT
                else
                    !Volume which does not fit into the channel
                    VolExcess = TotalVolume - ChannelsMaxVolume(i, j)
                    
                    !New Height of water in cell
                    NewH = VolExcess / Me%ExtVar%GridCellArea(i, j)
                    
                    !MaxFlow to or from river is calculated based on the level difference (new to old)
                    !If the new level is higher than the old level on channel, the max flow will be positive (flow to channel) 
                    dVol = (NewH + Me%ExtVar%Topography(i, j) - ChannelsWaterLevel(i, j)) * ChannelsTopArea(i,j)
                    MaxFlow = dVol / Me%ExtVar%DT                
                endif
                    
                dh         =  Me%myWaterLevel(i, j) - ChannelsWaterLevel(i, j)
                WaveHeight =  max(Me%myWaterLevel(i, j), ChannelsWaterLevel(i, j)) - Me%ExtVar%Topography(i,j)
                Celerity   = sqrt(Gravity * WaveHeight)
                
                if (dh > 0) then !flux is occuring between dh and with celerity                     
                    !m3/s = m/s (celerity) * m2 (Area = (dh * L) * 2)
                    Flow    = Celerity * 2.0 * ChannelsNodeLength(i, j) * min(dh, WaveHeight)
                else
                    !This was updated
                    !dh_new = (ChannelsSurfaceWidth(i,j) * dh) / &
                    !         (ChannelsSurfaceWidth(i,j) + 2 * min (Celerity * Me%ExtVar%DT, 0.5 * Me%ExtVar%DUX(i,j)))
                    dh_new = (ChannelsSurfaceWidth(i,j) * dh) / &
                             (ChannelsSurfaceWidth(i,j) + 2 * Celerity * Me%ExtVar%DT)
                    
                    !m3/s = h * L * Length / s
                    Flow    = -1. * (dh_new - dh) * ChannelsTopArea(i,j) / Me%ExtVar%DT
                endif

                if (abs(Flow) > abs(MaxFlow)) then
                    Flow = MaxFlow
                endif                                  
             
                !Updates Volumes                            
                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)                
                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)
                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)                    
                   
                !!Important!! flow to channel may have other sources than this, so a sum is needed
                Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Flow                  
            endif
                              
        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR080'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR090'
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR100'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR110'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR120'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsSurfaceWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR130'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsTopArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandChannelInteraction_4 - ModuleRunOff - ERR140'         
    
    end subroutine OverLandChannelInteraction_4
    
    !--------------------------------------------------------------------------       
    
    subroutine OverLandChannelInteraction_New
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: Flow, MaxFlow
        real   , dimension(:, :), pointer           :: ChannelsVolume
        real   , dimension(:, :), pointer           :: ChannelsMaxVolume
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength
        real                                        :: dh, WaveHeight
        integer, dimension(:, :), pointer           :: ChannelsActiveState
        

        call GetChannelsVolume      (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'     

        call GetChannelsMaxVolume   (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'   

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR02'

        call GetChannelsActiveState (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'        


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint .and. &   !RiverPoint
                ChannelsActiveState  (i, j) == BasinPoint .and. &   !Active
                ChannelsMaxVolume    (i, j) > 0.0) then             !Not the outlet


                !dh > 0, flow to channels, dh < 0, flow from channels
                dh         =  Me%myWaterLevel(i, j) - ChannelsWaterLevel(i, j)
                WaveHeight =  max(Me%myWaterLevel(i, j), ChannelsWaterLevel(i, j)) - Me%ExtVar%Topography(i,j)

                !flux is occuring between dh and with celerity (with wave height)
                !m3/s = m/s (celerity) * m2 (Area = (dh * L) * 2)
                Flow = sqrt(Gravity * WaveHeight) * 2.0 * ChannelsNodeLength(i, j) * dh

                if (dh > 0) then
                    MaxFlow = Me%myWaterVolume (i, j) / Me%ExtVar%DT
                else
                    MaxFlow = -1. * (ChannelsVolume(i,j) - ChannelsMaxVolume(i,j)) / Me%ExtVar%DT
                endif

                if (abs(Flow) > abs(MaxFlow)) then
                    Flow = MaxFlow
                endif   
                    
                !!Important!! flow to channel may have other sources than this, so a sum is needed
                Me%iFlowToChannels(i, j) = Me%iFlowToChannels(i, j) + Flow

                !Updates Volumes
                !Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - Me%iFlowToChannels    (i, j) * Me%ExtVar%DT
                Me%myWaterVolume (i, j) = Me%myWaterVolume (i, j) - (Flow * Me%ExtVar%DT)
                
                Me%myWaterColumn (i, j) = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                Me%myWaterLevel  (i, j) = Me%myWaterColumn (i, j) + Me%ExtVar%Topography  (i, j)
                           
            endif

        enddo
        enddo        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsMaxVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR07'        

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsActiveState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowFromChannels - ModuleRunOff - ERR010'        

    
    end subroutine OverLandChannelInteraction_New
    
    !--------------------------------------------------------------------------

    subroutine CheckStability (Restart)

        !Arguments-------------------------------------------------------------
        logical                                     :: Restart

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, n_restart
        real                                        :: variation        

        !Begin-----------------------------------------------------------------

        n_restart = 0
        Restart = .false.
        
        !Verifies negative volumes
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                if (Me%myWaterVolume(i, j) < -AllmostZero) then
!                    write(*,*) '-----'
!                    write(*,*) 'OldVolume ', Me%myWaterVolumeOld(i, j)
!                    write(*,*) 'Negative Volume - Me%myWaterVolume (', i, ', ', j, ') =', Me%myWaterVolume (i, j)
!                    write(*,*) '-----'
                    Restart = .true.                 
                    exit do1
                else if (Me%myWaterVolume (i, j) < 0.0) then  
                    Me%myWaterVolume (i, j) = 0.0                 
                endif
            endif
        enddo
        enddo do1
        
        if ((.not. Restart) .and. Me%CV%Stabilize) then
do2:        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                if (Me%StabilityPoints(i, j) == BasinPoint) then
            
                    if ((.not. Me%CV%CheckDecreaseOnly) .or. (Me%myWaterVolumeOld(i, j) > Me%myWaterVolume(i, j))) then
                        if (Me%myWaterVolumeOld(i, j) / Me%ExtVar%GridCellArea(i, j) >= Me%CV%MinimumValueToStabilize) then
                            
                            variation = abs(Me%myWaterVolume(i, j) - Me%myWaterVolumeOld(i, j)) / Me%myWaterVolumeOld(i, j)
                            
                            if (variation > Me%CV%StabilizeFactor) then
                                !Debug routine - may be usefull for using in debug situation
                                !call DebugStability (i,j,variation)                                
                                
                                n_restart = n_restart + 1
                                
                                if (n_restart > Me%CV%MinToRestart) then
                                    Restart = .true.
                                    exit do2
                                endif                                 
                            endif
                        endif
                    endif
                endif
            enddo
            enddo do2              
                        
        endif
        
        if (Restart) then        
            Me%CV%NextNiteration = max(int(Me%CV%NextNiteration * Me%CV%DTSplitFactor), Me%CV%NextNiteration + 1)
                 
            if (Me%CV%NextNiteration >= Me%CV%MaxIterations) then
                 write(*,*)'Number of iterations above maximum: ', Me%CV%NextNiteration
                 stop 'CheckStability - ModuleRunoff - ERR010'
            endif                          
        endif           
        
    end subroutine CheckStability
 
    !--------------------------------------------------------------------------

    subroutine DebugStability(i,j, variation)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: I, J
        real                                        :: variation
        !Local-----------------------------------------------------------------
        character (Len = 5)                         :: str_i, str_j
        character (Len = 15)                        :: str_1, str_2, str_3
        character (len = StringLength)              :: string_to_be_written 
        
        write(str_i, '(i3)') i 
        write(str_j, '(i3)') j
        write(str_1, '(ES10.3)') Me%myWaterVolumeOld(I,J)  
        write(str_2, '(ES10.3)') Me%myWaterVolume(I,J)   
        write(str_3, '(ES10.3)') variation                            
        
        string_to_be_written = ' '//str_i//','//str_j//' '//str_1//' '//str_2//' '//str_3
        
        call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)           
    
    
    end subroutine DebugStability
    
    !--------------------------------------------------------------------------

    subroutine LocalWaterColumn (WaterColumn)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), pointer              :: WaterColumn

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                Me%myWaterVolume(i, j) = WaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)

                Me%myWaterColumn(i, j) = WaterColumn(i, j)

                Me%myWaterLevel (i, j) = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
            endif
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL        
        
    end subroutine LocalWaterColumn            

    !--------------------------------------------------------------------------

    subroutine IntegrateFlow (LocalDT, SumDT)

        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT, SumDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J)

        !Integrates along X Directions
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%iFlowX(i, j) = (Me%iFlowX(i, j) * SumDT + Me%lFlowX(i, j) * LocalDT) / &
                              (SumDT + LocalDT)
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Integrates along Y Directions
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%iFlowY(i, j) = (Me%iFlowY(i, j) * SumDT + Me%lFlowY(i, j) * LocalDT) / &
                              (SumDT + LocalDT)
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Integrates Flow to Channels
        if (Me%ObjDrainageNetwork /= 0) then
           !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%iFlowToChannels(i, j) = (Me%iFlowToChannels(i, j) * SumDT + Me%lFlowToChannels(i, j) * LocalDT) / &
                                           (SumDT + LocalDT)
            enddo
            enddo
            !$OMP END DO NOWAIT
        endif
        
        !Integrates Flow At boundary
!        if (Me%ImposeBoundaryValue) then
!           !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!                Me%iFlowBoundary(i, j) = (Me%iFlowBoundary(i, j) * SumDT + Me%lFlowBoundary(i, j) * LocalDT) / &
!                                         (SumDT + LocalDT)
!            enddo
!            enddo
!            !$OMP END DO NOWAIT
!        endif

        !Integrates Flow Discharges
        if (Me%Discharges) then
           !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%iFlowDischarge(i, j) = (Me%iFlowDischarge(i, j) * SumDT + Me%lFlowDischarge(i, j) * LocalDT) / &
                                          (SumDT + LocalDT)
            enddo
            enddo
            !$OMP END DO NOWAIT
        endif

        !$OMP END PARALLEL        



    end subroutine IntegrateFlow

    !--------------------------------------------------------------------------

    subroutine ImposeBoundaryValue ()
    
        !Arguments-------------------------------------------------------------
        !real                                        :: LocalDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, di, dj
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: dh, dVOl
        !logical                                     :: NearBoundary
        real                                        :: AreaZX, AreaZY !, Width
        real                                        :: WaveHeight, Celerity, MaxFlow

        !Routes water outside the watershed if water is higher then a given treshold values
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !Default is zero
        Me%iFlowBoundary = 0.0
        
        !Sets Boundary values
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j)  == BasinPoint .and. &      !BasinPoint
                Me%BoundaryCells     (i,j)   == BasinPoint .and. &      !BoundaryPoints
!                Me%myWaterColumn     (i, j)  > 0.1         .and. &      !"Flooding"
                Me%myWaterLevel      (i, j)  > Me%BoundaryValue) then   !Level higher then imposed level

!                !Check if near a boundary point
!                NearBoundary = .false.
!                do dj = -1, 1
!                do di = -1, 1
!                    if (Me%ExtVar%BasinPoints(i+di, j+dj) == 0) then
!                        NearBoundary = .true.
!                    endif
!                enddo
!                enddo
!
!                if (NearBoundary) then

                    !Necessary Variation in height - always positive because only evaluates cell as so
                    dh = Me%myWaterLevel (i, j) - Me%BoundaryValue
                    
                    if (dh > Me%MinimumWaterColumn) then
                        
!                        !Cell Width
!                        Width                  = (Me%ExtVar%DYY(i, j) + Me%ExtVar%DXX(i, j)) / 2.0
                        
                        !celerity is limited by water column size and not dh
                        WaveHeight = Me%myWaterColumn(i, j) 

!                        !Area for Flow
!                        Area                   = Width * min(dh, WaveHeight)
                        
                        Celerity = sqrt(Gravity * WaveHeight)

                        !Flow to set cell equal to Boundary Value
                        !m3/s                  = 
!                        Me%lFlowBoundary(i, j) = Min(0.5 * dh * Me%ExtVar%GridCellArea(i, j) / LocalDT,         &
!                                                     0.5 * Area * sqrt(Gravity * dh))

                        !U direction - use middle area because in closed faces does not exist AreaU
                        !flow negative (exiting runoff)
                        AreaZX = Me%ExtVar%DVY(i,j) * Me%myWaterColumn(i,j)
                        do dj = 0, 1
                            if ((Me%ComputeFaceU(i,j+dj) == 0)) then
                                Me%iFlowBoundary(i, j) = Me%iFlowBoundary(i, j) - AreaZX * Celerity
                            endif
                        enddo

                        !V direction - use middle area because in closed faces does not exist AreaV
                        AreaZY = Me%ExtVar%DUX(i,j) * Me%myWaterColumn(i,j)                      
                        do di = 0, 1
                            if ((Me%ComputeFaceV(i+di,j) == 0)) then
                                Me%iFlowBoundary(i, j) = Me%iFlowBoundary(i, j) - AreaZY * Celerity
                            endif
                        enddo
                        
                        !cant remove more than up to boundary or water column if boundary lower than topography 
                        !negative flow 
                        !m3/s = m * m2 / s
                        MaxFlow = - min(dh, WaveHeight) *  Me%ExtVar%GridCellArea(i, j) / Me%ExtVar%DT
                        
                        !m3/s = m2 * m/s 
                        Me%iFlowBoundary(i, j)  = max (Me%iFlowBoundary(i, j), MaxFlow)
                        
                        !dVol
                        dVol = Me%iFlowBoundary(i, j) * Me%ExtVar%DT
                            
                        !Updates Water Volume
                        Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   + dVol 
                            
                        !Updates Water Column
                        Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)  / Me%ExtVar%GridCellArea(i, j)

                        !Updates Water Level
                        Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)  + Me%ExtVar%Topography(i, j)
                        
                    else
                    
                        Me%iFlowBoundary(i, j) = 0.0
                    
                    endif

!                endif
           
            endif
        enddo
        enddo

    
    end subroutine ImposeBoundaryValue
    
    !--------------------------------------------------------------------------

    subroutine ImposeBoundaryValue_v2
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j !, di, dj
        integer                                     :: ILB, IUB, JLB, JUB
!        logical                                     :: NearBoundary
        real                                        :: OldVolume

        !Routes water outside the watershed if water is higher then a given treshold values
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        Me%BoundaryFlowVolume = 0.0
        
        !Sets Boundary values
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j)  == BasinPoint       .and.    &   !BasinPoint
                Me%BoundaryCells     (i,j)   == BasinPoint       .and.    &   !BoundaryPoints
                Me%ExtVar%Topography (i, j)  < Me%MaxDtmForBoundary .and. &   !Low land point where to imposes BC
                Me%myWaterLevel      (i, j)  > Me%BoundaryValue) then   !Level higher then imposed level

!                !Check if near a boundary point
!                NearBoundary = .false.
!                do dj = -1, 1
!                do di = -1, 1
!                    if (Me%ExtVar%BasinPoints(i+di, j+dj) == 0) then
!                        NearBoundary = .true.
!                    endif
!                enddo
!                enddo
!
!                if (NearBoundary) then

                    !Necessary Variation in height
                    Me%myWaterLevel (i, j) = max(Me%BoundaryValue, Me%ExtVar%Topography (i, j))

                    !Updates Water Column
                    Me%myWaterColumn(i, j) = Me%myWaterLevel (i, j) - Me%ExtVar%Topography (i, j) 
                    
                    !Updates Volume and BoundaryFlowVolume
                    OldVolume              = Me%myWaterVolume(i, j)
                    
                    !m3 = m * m2
                    Me%myWaterVolume(i, j) = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                    
                    !m3 = m3 + (m3 - m3)
                    Me%BoundaryFlowVolume  = Me%BoundaryFlowVolume + (OldVolume - Me%myWaterVolume(i, j))
                    
                    !m3/s = m3 / s - always negative exiting runoff
                    Me%iFlowBoundary(i, j) = (Me%myWaterVolume(i, j) - OldVolume) / Me%ExtVar%DT
                    
!                endif
           
            endif
        enddo
        enddo

    
    end subroutine ImposeBoundaryValue_v2       
    !--------------------------------------------------------------------------
    
    subroutine ComputeCenterValues 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: CHUNK

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
       
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
            
        call SetMatrixValue(Me%CenterFlowX, Me%Size, 0.0)
        call SetMatrixValue(Me%CenterFlowY, Me%Size, 0.0)
        call SetMatrixValue(Me%FlowModulus, Me%Size, 0.0)

        call SetMatrixValue(Me%CenterVelocityX, Me%Size, 0.0)
        call SetMatrixValue(Me%CenterVelocityY, Me%Size, 0.0)
        call SetMatrixValue(Me%VelocityModulus, Me%Size, 0.0)

        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
                
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                Me%CenterFlowX(i, j) = (Me%iFlowX(i, j) + Me%iFlowX(i, j+1)) / 2.0
                Me%CenterFlowY(i, j) = (Me%iFlowY(i, j) + Me%iFlowY(i+1, j)) / 2.0
                Me%FlowModulus(i, j) = sqrt (Me%CenterFlowX(i, j)**2. + Me%CenterFlowY(i, j)**2.)
                
                if (Me%myWaterColumn (i,j) > AllmostZero) then
                    Me%CenterVelocityX (i, j) = Me%CenterFlowX (i,j) / ( Me%ExtVar%DYY(i, j) * Me%myWaterColumn (i,j) )
                    Me%CenterVelocityY (i, j) = Me%CenterFlowY (i,j) / ( Me%ExtVar%DXX(i, j) * Me%myWaterColumn (i,j) )
                    Me%VelocityModulus (i, j) = sqrt (Me%CenterVelocityX(i, j)**2.0 + Me%CenterVelocityY(i, j)**2.0)
                else
                    Me%myWaterColumn (i,j) = 0.0
                end if

                if(Me%WriteMaxFlowModulus) then
                    if (Me%FlowModulus(i, j) > Me%MaxFlowModulus(i, j)) then
                        Me%MaxFlowModulus(i, j) = Me%FlowModulus(i, j)
                    end if
                end if

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        
!        if (Me%StormWaterDrainage) then
!        
!            Me%StormWaterCenterFlowX    = 0.0
!            Me%StormWaterCenterFlowY    = 0.0
!            Me%StormWaterCenterModulus  = 0.0
!        
!            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!            do j = JLB, JUB
!            do i = ILB, IUB
!                    
!                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
!
!                    Me%StormWaterCenterFlowX(i, j)   = (Me%StormWaterFlowX(i, j) + Me%StormWaterFlowX(i, j+1)) / 2.0
!                    Me%StormWaterCenterFlowY(i, j)   = (Me%StormWaterFlowY(i, j) + Me%StormWaterFlowY(i+1, j)) / 2.0
!                    Me%StormWaterCenterModulus(i, j) = sqrt (Me%StormWaterCenterFlowX(i, j)**2. + &
!                                                             Me%StormWaterCenterFlowY(i, j)**2.)
!
!                endif
!
!            enddo
!            enddo
!            !$OMP END DO NOWAIT
!        
!        endif
        
        !$OMP END PARALLEL
        
    end subroutine ComputeCenterValues 
    
    !--------------------------------------------------------------------------
    
    subroutine ComputeNextDT (Niter)

        !Arguments-------------------------------------------------------------
        integer                                     :: Niter        
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: nextDTCourant, aux
        real                                        :: nextDTVariation, MaxDT
        logical                                     :: VariableDT
        real                                        :: vel, dist, currentDT
        

        !----------------------------------------------------------------------

        call GetVariableDT(Me%ObjTime, VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModuleRunOff -  ERR010'

        call GetMaxComputeTimeStep(Me%ObjTime, MaxDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModuleRunOff -  ERR020'

        nextDTCourant   = -null_real
        nextDTVariation = -null_real
        
        if (VariableDT) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB
            
            if (Me%CV%LimitDTCourant) then
                        
                do j = JLB, JUB
                do i = ILB, IUB
                        
                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint .and. Me%myWaterColumn (i,j) > Me%MinimumWaterColumn) then

                        vel = sqrt(Gravity * Me%myWaterColumn (i,j))
                        
                        if (vel .gt. 0) then
                        
                            !spatial step, in case of dx = dy, dist = sqrt(2) * dx
                            dist = sqrt ((Me%ExtVar%DZX(i, j)**2) + (Me%ExtVar%DZY(i, j)**2))
                            aux = dist * Me%CV%MaxCourant / vel 
                        
                            nextDTCourant = min(nextDTCourant, aux)
                            
                        endif
                            
                    endif

                enddo
                enddo
                
            endif
            
            if (Niter == 1) then
            
                nextDTVariation = Me%ExtVar%DT * Me%CV%DTFactorUp
                Me%CV%NextNiteration = Niter
                
            elseif (Niter <= Me%CV%MinIterations) then                            
            
                if (Niter > Me%CV%LastGoodNiteration) then

                    nextDTVariation = Me%ExtVar%DT
                    Me%CV%NextNiteration = Niter

                else
                
                    nextDTVariation = Me%ExtVar%DT * Me%CV%DTFactorUp
                    Me%CV%NextNiteration = Niter

                endif
                
            else
            
                if (Niter >= Me%CV%StabilizeHardCutLimit) then
                
                    nextDTVariation = (Me%ExtVar%DT / Niter) * Me%CV%MinIterations
                    Me%CV%NextNiteration = Me%CV%MinIterations
                    
                elseif (Niter > Me%CV%LastGoodNiteration) then
                
                    nextDTVariation = Me%ExtVar%DT / Me%CV%DTFactorDown
                    Me%CV%NextNiteration = max(int(nextDTVariation / Me%CV%CurrentDT), 1)
                    
                else
                
                    nextDTVariation = Me%ExtVar%DT
                    Me%CV%NextNiteration = max(min(int(Niter / Me%CV%DTSplitFactor), Niter - 1), 1)
                    
                endif 
                               
            endif
            
            CurrentDT = nextDTVariation / Me%CV%NextNiteration                                     
                      
            Me%CV%NextDT = min(min(nextDTVariation, nextDTCourant), MaxDT)
            
            if (Me%CV%NextDT < nextDTVariation) then                
                Me%CV%NextNiteration = max(int(Me%CV%NextDT/CurrentDT), 1)
            endif
                       
        else
        
            Me%CV%NextDT = Me%ExtVar%DT
            Me%CV%NextNiteration = Niter            
        
        endif
        
        Me%CV%LastGoodNiteration = Niter
        Me%CV%CurrentDT          = Me%CV%NextDT / Me%CV%NextNiteration
    
    end subroutine ComputeNextDT

    !--------------------------------------------------------------------------
    
    subroutine RunOffOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        real, dimension(6)  , target                :: AuxTime
        real, dimension(:)  , pointer               :: TimePointer       

        if (MonitorPerformance) call StartWatch ("ModuleRunOff", "RunOffOutput")

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


        if (Me%ExtVar%Now >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

            !Writes current time
            call ExtractDate   (Me%ExtVar%Now , AuxTime(1), AuxTime(2),         &
                                                AuxTime(3), AuxTime(4),         &
                                                AuxTime(5), AuxTime(6))
            TimePointer => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR01'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                 "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D      = TimePointer,                    &
                                 OutputNumber = Me%OutPut%NextOutPut,           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR02'

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR03'

            !Writes Flow values
            !Writes the Water Column - should be on runoff
            call HDF5WriteData   (Me%ObjHDF5, "//Results/water column",         &
                                  "water column", "m",                          &
                                  Array2D      = Me%MyWaterColumn,              &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleRunOff - ERR050'

       
            !Writes the Water Level
            call HDF5WriteData   (Me%ObjHDF5, "//Results/water level",          &
                                  "water level", "m",                           &
                                  Array2D      = Me%MyWaterLevel,               &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleRunOff - ERR060'
            


            !Writes Flow X
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/flow X",                                &
                                  "flow X",                                         &   
                                  "m3/s",                                           &
                                  Array2D      = Me%CenterFlowX,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'

            
            !Writes Flow Y
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/flow Y",                                &
                                  "flow Y",                                         &   
                                  "m3/s",                                           &
                                  Array2D      = Me%CenterFlowY,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

             !Writes Flow Modulus
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/"//trim(GetPropertyName (FlowModulus_)),&
                                  trim(GetPropertyName (FlowModulus_)),             &   
                                  "m3/s",                                           &
                                  Array2D      = Me%FlowModulus,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

             !Writes Velocity X 
            call HDF5WriteData   (Me%ObjHDF5,                                          &
                                  "/Results/"//trim(GetPropertyName (VelocityU_)),     &
                                  trim(GetPropertyName (VelocityU_)),                  &
                                  "m/s",                                               &
                                  Array2D      = Me%CenterVelocityX,                   &
                                  OutputNumber = Me%OutPut%NextOutPut,                 &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

             !Writes Velocity Y 
            call HDF5WriteData   (Me%ObjHDF5,                                          &
                                  "/Results/"//trim(GetPropertyName (VelocityV_)),     &
                                  trim(GetPropertyName (VelocityV_)),                  &
                                  "m/s",                                               &
                                  Array2D      = Me%CenterVelocityY,                   &
                                  OutputNumber = Me%OutPut%NextOutPut,                 &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

            !Writes Velocity Modulus
            call HDF5WriteData   (Me%ObjHDF5,                                                &
                                  "/Results/"//trim(GetPropertyName (VelocityModulus_)),     &
                                  trim(GetPropertyName (VelocityModulus_)),                  &
                                  "m/s",                                                     &
                                  Array2D      = Me%VelocityModulus,                         &
                                  OutputNumber = Me%OutPut%NextOutPut,                       &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

            !Writes Storm Water Volume of each Cell
            if (Me%StormWaterDrainage) then
                call HDF5WriteData   (Me%ObjHDF5, "//Results/storm water volume",      &
                                      "storm water volume", "m3",                      &
                                      Array2D      = Me%StormWaterVolume,              &
                                      OutputNumber = Me%OutPut%NextOutPut,             &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleRunOff - ERR080'
                
                
                !Writes Flow X
                call HDF5WriteData   (Me%ObjHDF5,                                       &
                                      "//Results/storm water flow X",                   &
                                      "storm water flow X",                             &   
                                      "m3/s",                                           &
                                      Array2D      = Me%StormWaterCenterFlowX,          &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'

                
                !Writes SW Flow Y
                call HDF5WriteData   (Me%ObjHDF5,                                       &
                                      "//Results/storm water flow Y",                   &
                                      "storm water flow Y",                             &   
                                      "m3/s",                                           &
                                      Array2D      = Me%StormWaterCenterFlowY,          &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'
                

               
                !Writes SW Modulus
                call HDF5WriteData   (Me%ObjHDF5,                                       &
                                      "//Results/storm water flow modulus",             &
                                      "storm water flow modulus",                       &   
                                      "m3/s",                                           &
                                      Array2D      = Me%StormWaterCenterModulus,        &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'
                
                
            endif
            
            if (Me%StormWaterModel) then
            
                call HDF5WriteData   (Me%ObjHDF5, "//Results/storm water flow",         &
                                      "storm water flow", "m3/s",                       &
                                      Array2D      = Me%StormWaterModelFlow,            &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleRunOff - ERR085'
                
                call HDF5WriteData   (Me%ObjHDF5, "//Results/street gutter flow",       &
                                      "street gutter flow", "m3/s",                     &
                                      Array2D      = Me%StreetGutterFlow,               &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleRunOff - ERR085'
            
                call HDF5WriteData   (Me%ObjHDF5, "//Results/storm water real flow",    &
                                      "storm water real flow", "m3/s",                  &
                                      Array2D      = Me%StormInteractionFlow,           &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleRunOff - ERR085'
           
            endif

           
            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR99'

            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

        endif

         if (MonitorPerformance) call StopWatch ("ModuleRunOff", "RunOffOutput")
        
    end subroutine RunOffOutput

    !--------------------------------------------------------------------------


    subroutine OutputMaxWaterColumn

        !Locals----------------------------------------------------------------
        integer                                 :: ILB,IUB, JLB, JUB, i, j
        integer                                 :: STAT_CALL
        real, dimension(:,:), pointer           :: ChannelsWaterLevel
        

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
   
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                !Water Column of overland flow
                if (Me%myWaterColumn(i, j) > Me%MaxWaterColumn(i, j)) then
                    Me%MaxWaterColumn(i, j) = Me%myWaterColumn(i, j)
                endif
            endif

        enddo
        enddo

        if (Me%ObjDrainageNetwork /= 0) then

            call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputMaxWaterColumn - ModuleRunOff - ERR01'     

            do j = JLB, JUB
            do i = ILB, IUB
       
                !Water Column of River Network
                if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) then
                    if (ChannelsWaterLevel(i, j) - Me%ExtVar%Topography(i, j) > Me%MaxWaterColumn(i, j)) then
                        Me%MaxWaterColumn(i, j) = ChannelsWaterLevel(i, j) - Me%ExtVar%Topography(i, j)
                    endif
                endif

            enddo
            enddo


            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutputMaxWaterColumn - ModuleRunOff - ERR05'

        endif


    end subroutine OutputMaxWaterColumn

    !---------------------------------------------------------------------------

!    subroutine  WriteChannelsLevelData
!
!        !Local-------------------------------------------------------------------
!        integer                                                 :: ILB,IUB, JLB, JUB
!        integer                                                 :: STAT_CALL,i,j
!        integer, dimension (:,:), pointer                       :: ChannelsID
!        character(len=StringLength), dimension (:,:), pointer   :: ChannelsStationName
!
!        !------------------------------------------------------------------------
!
!        call GetChannelsID  (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR01'
!
!        call GetChannelsStationName  (Me%ObjDrainageNetwork, ChannelsStationName, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR02'
!
!        call GetRiverPoints (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR02a'
!
!
!        !GetNodeID
!        !GetNodeStationName
!
!        open(UNIT=UnitMax, FILE=Me%MaxWaterColumnFile, ACTION='WRITE', STATUS='REPLACE', IOSTAT=STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR03'
!
!
!
!        write(UnitMax,*) 'NodeID     MaxWaterColumn DateTime            StationName'
!
!        ILB = Me%WorkSize%ILB
!        IUB = Me%WorkSize%IUB
!        JLB = Me%WorkSize%JLB
!        JUB = Me%WorkSize%JUB
!        
!        do j = JLB, JUB
!        do i = ILB, IUB
!
!            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) &
!                write(UnitMax,100) ChannelsID(i,j), Me%MaxWaterColumn(i,j), Me%MaxWaterColumnTime(i,j), &
!                trim(adjustl(ChannelsStationName(i,j)))
!
!        enddo
!        enddo
!       
!        close(UnitMax)
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR04'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsStationName, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR05'
!
!        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteChannelsLevelData - ModuleRunOff - ERR05a'
!
!
!        100 format(I10,1x, f16.3, 1x, A19, 1x, A)   
!
!    end subroutine  WriteChannelsLevelData


    !----------------------------------------------------------------------------
   
    real function AdjustSlope (Slope)
    
        !Arguments--------------------------------------------------------------
        real                                    :: Slope
        real                                    :: sign

        !Slope correction given by City of Albuquerque, 1997, p.22-26
        !http://www.hkh-friend.net.np/rhdc/training/lectures/HEGGEN/Tc_3.pdf


        if (Slope.LT.0.0) then
            sign = -1.0
        else
            sign = 1.0
        end if

        Slope = abs (Slope)
        
        if (Slope.GE.0.04 .and. Me%AdjustSlope) then
            Slope = 0.05247 + 0.06363 * Slope - 0.182 * exp (-62.38 * Slope)
        end if
        
        AdjustSlope = sign * Slope
        

    end function AdjustSlope

    !----------------------------------------------------------------------------

    subroutine CalculateTotalStoredVolume

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        
        Me%TotalStoredVolume = 0.0
        
        Me%VolumeStoredInSurface     = 0.0
        Me%VolumeStoredInStormSystem = 0.0

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
                !m3 = m3 + m3
                Me%TotalStoredVolume = Me%TotalStoredVolume + Me%MyWaterVolume(i, j)
                
                !m3 = m3 + m3
                Me%VolumeStoredInSurface = Me%VolumeStoredInSurface + Me%MyWaterVolume(i, j)
            endif

        enddo
        enddo


        if (Me%StormWaterDrainage) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                if (Me%ExtVar%BasinPoints(i, j) == 1) then                    
                    !m3 = m3 + m3
                    Me%TotalStoredVolume = Me%TotalStoredVolume + Me%StormWaterVolume(i, j)
                         
                    !m3 = m3 + m3
                    Me%VolumeStoredInStormSystem = Me%VolumeStoredInStormSystem + Me%StormWaterVolume(i, j)                    
                endif

            enddo
            enddo
                
        endif
        

    end subroutine CalculateTotalStoredVolume

    !--------------------------------------------------------------------------
    
    subroutine WriteFinalFile
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: FinalFile
        integer                                     :: STAT_CALL
        character(LEN = PathLength)                 :: FileName
        
        !----------------------------------------------------------------------

        !Gets Date
        call ExtractDate(Me%ExtVar%Now, Year_File, Month_File, Day_File,               &
                         Hour_File, Minute_File, Second_File)
        
        
        if (Me%ExtVar%Now == Me%EndTime) then
            FileName = Me%Files%FinalFile
        else
            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%ExtVar%Now))//".fin")
        endif            
        
        call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR01'

        open(Unit = FinalFile, File = FileName, Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR02'

        !Writes Date
        write(FinalFile) Year_File, Month_File, Day_File, Hour_File, Minute_File,       &
                         Second_File

        write(FinalFile)Me%myWaterColumn
        
        if (Me%StormWaterDrainage) then
            write(FinalFile)Me%StormWaterVolume
        endif


        call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoff - ERR03'

    end subroutine WriteFinalFile

    !------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillRunOff(RunOffID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: RunOffID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, STAT_CALL    
        character(len=StringLength)         :: MassErrorFile

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunOffID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then


            nUsers = DeassociateInstance(mRUNOFF_,  Me%InstanceID)

            if (nUsers == 0) then

                !Writes file with final condition
                call WriteFinalFile

                !Writes Mass Error
                call ReadFileName("ROOT_SRT", MassErrorFile, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR02a'
                MassErrorFile = trim(adjustl(MassErrorFile))//"MassError.dat"
                
                call WriteGridData  (MassErrorFile,                            &
                     COMENT1          = "MassErrorFile",                       &
                     COMENT2          = "MassErrorFile",                       &
                     HorizontalGridID = Me%ObjHorizontalGrid,                  &
                     FillValue        = -99.0,                                 &
                     OverWrite        = .true.,                                &
                     GridData2D_Real  = Me%MassError,                          &
                     STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR00'
        
                if(Me%WriteMaxFlowModulus) then
                    call WriteGridData  (Me%MaxFlowModulusFile,                &
                         COMENT1          = "MaxFlowModulusFile",              &
                         COMENT2          = "MaxFlowModulusFile",              &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%MaxFlowModulus,                 &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR00'
                endif
                
                if (Me%WriteMaxWaterColumn) then
                    call WriteGridData  (Me%MaxWaterColumnFile,                &
                         COMENT1          = "MaxWaterColumnFile",              &
                         COMENT2          = "MaxWaterColumnFile",              &
                         HorizontalGridID = Me%ObjHorizontalGrid,              &
                         FillValue        = -99.0,                             &
                         OverWrite        = .true.,                            &
                         GridData2D_Real  = Me%MaxWaterColumn,                 &
                         STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - RunOff - ERR00'
                endif


                if (Me%ObjDrainageNetwork /= 0) then
 
!                    if(Me%WriteMaxWaterColumn) call WriteChannelsLevelData

                    nUsers = DeassociateInstance (mDRAINAGENETWORK_, Me%ObjDrainageNetwork)
                    if (nUsers == 0) stop 'KillRunOff - RunOff - ERR01'
                endif

                if (Me%OutPut%Yes) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunOff - ERR05'
                endif
                
                if (Me%Discharges) then
                    call Kill_Discharges(Me%ObjDischarges, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunOff - ERR05a'
                endif

                !Deassociates External Instances
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR05'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR06'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjGridData)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR07'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR08'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillRunOff - RunOff - ERR09'
                
                deallocate(Me%myWaterColumnOld)
                
                deallocate (Me%iFlowX)
                deallocate (Me%iFlowY)
                deallocate (Me%lFlowX)
                deallocate (Me%lFlowY)
                deallocate (Me%iFlowToChannels)
                deallocate (Me%lFlowToChannels)
                deallocate (Me%lFlowBoundary)
                deallocate (Me%iFlowBoundary)
                deallocate (Me%iFlowRouteDFour)

                nullify    (Me%iFlowX)
                nullify    (Me%iFlowY)
                nullify    (Me%lFlowX)
                nullify    (Me%lFlowY)
                nullify    (Me%iFlowToChannels)
                nullify    (Me%lFlowToChannels)
                nullify    (Me%lFlowBoundary)
                nullify    (Me%iFlowBoundary)
                nullify    (Me%iFlowRouteDFour)


                !Deallocates Instance
                call DeallocateInstance ()

                RunOffID   = 0
                STAT_      = SUCCESS_

            end if


        end if cd1


        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillRunOff

    !------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_RunOff), pointer                    :: AuxObjRunOff
        type (T_RunOff), pointer                    :: PreviousObjRunOff

        !Updates pointers
        if (Me%InstanceID == FirstObjRunOff%InstanceID) then
            FirstObjRunOff => FirstObjRunOff%Next
        else
            PreviousObjRunOff => FirstObjRunOff
            AuxObjRunOff      => FirstObjRunOff%Next
            do while (AuxObjRunOff%InstanceID /= Me%InstanceID)
                PreviousObjRunOff => AuxObjRunOff
                AuxObjRunOff      => AuxObjRunOff%Next
            enddo

            !Now update linked list
            PreviousObjRunOff%Next => AuxObjRunOff%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready (RunOffID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (RunOffID > 0) then
            call LocateObjRunOff (RunOffID)
            ready_ = VerifyReadLock (mRUNOFF_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjRunOff (ObjRunOffID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunOffID

        !Local-----------------------------------------------------------------

        Me => FirstObjRunOff
        do while (associated (Me))
            if (Me%InstanceID == ObjRunOffID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleRunOff - LocateObjRunOff - ERR01'

    end subroutine LocateObjRunOff

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar (StaticOnly)
        
        !Arguments-------------------------------------------------------------
        logical                                     :: StaticOnly

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Time Stuff
        call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR01'

        call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR02'

        !Gets Basin Points
        call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR03'
        
        !Gets cell slope
        call GetCellSlope   (Me%ObjBasinGeometry, Me%ExtVar%CellSlope, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR04'

        !Gets River Points
        call GetRiverPoints (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR05'

        !Gets Horizontal Grid
        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                     &
                               DUX    = Me%ExtVar%DUX,    DVY    = Me%ExtVar%DVY,        &
                               DXX    = Me%ExtVar%DXX,    DYY    = Me%ExtVar%DYY,        &
                               DZX    = Me%ExtVar%DZX,    DZY    = Me%ExtVar%DZY,        &
                               XX2D_Z = Me%ExtVar%XX2D_Z, YY2D_Z = Me%ExtVar%YY2D_Z,     &
                               STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR06'

        call GetGridCellArea  (Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea,             &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR06a'

        !Gets a pointer to Topography
        call GetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR07'

        if (.not. StaticOnly) then

            !Gets Boundary Points
            call GetBoundaries    (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunOff - ERR10'
        
        endif

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar(StaticOnly)
        
        !Arguments-------------------------------------------------------------
        logical                                     :: StaticOnly
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Unget Basin Points
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR01'

        !Unget River Points
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR02'

        !Unget Cell Slope
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%CellSlope, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR02a'

        !Unget Horizontal Grid
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR03'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR04'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DXX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR06'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR06'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%XX2D_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR07'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%YY2D_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR08'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR09'

        !Ungets the Topography
        call UngetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR10'

        if (.not. StaticOnly) then

            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunOff - ERR11'

        endif 
        
    end subroutine ReadUnLockExternalVar
    
#ifdef _OPENMI_


    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::IsUrbanDrainagePoint
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_ISURBANDRAINAGEPOINT"::IsUrbanDrainagePoint
    !DEC$ ENDIF
    logical function IsUrbanDrainagePoint(RunOffID, i, j)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer                                     :: i, j
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(RunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            if (Me%StormWaterInteraction(i, j) > AllmostZero) then        
                IsUrbanDrainagePoint = .true.
            else
                IsUrbanDrainagePoint = .false.
            endif
        else 
            IsUrbanDrainagePoint = .false.
        end if
           
        return

    end function IsUrbanDrainagePoint



    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetPondedWaterColumn
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETPONDEDWATERCOLUMN"::GetPondedWaterColumn
    !DEC$ ENDIF
    logical function GetPondedWaterColumn(RunOffID, nComputePoints, waterColumn)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer                                     :: nComputePoints
        real(8), dimension(nComputePoints)          :: waterColumn
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        integer                                     :: i, j, idx

        call Ready(RunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            idx = 1
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%StormWaterInteraction(i, j) > AllmostZero) then
                    waterColumn(idx) = Max(Me%MyWaterColumn(i, j) - Me%MinimumWaterColumn, 0.0)
                    idx = idx + 1
                endif
            enddo
            enddo

            GetPondedWaterColumn = .true.
        else 
            call PlaceErrorMessageOnStack("Runoff not ready")
            GetPondedWaterColumn = .false.
        end if

    end function GetPondedWaterColumn

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetInletInFlow
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETINLETINFLOW"::GetInletInFlow
    !DEC$ ENDIF
    logical function GetInletInFlow(RunOffID, nComputePoints, inletInflow)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer                                     :: nComputePoints
        real(8), dimension(nComputePoints)          :: inletInflow
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        integer                                     :: i, j, idx
        integer                                     :: targetI
        integer                                     :: targetJ

        call Ready(RunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            !Puts values into 1D OpenMI matrix
            idx = 1
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%StormWaterInteraction (i, j) > AllmostZero) then 
                
                    !inlet Flow rate min between 
                    inletInflow(idx) = Me%SewerInflow(i, j)
                    idx = idx + 1
                endif
                    
            enddo
            enddo


            GetInletInFlow = .true.
        else 
            call PlaceErrorMessageOnStack("Runoff not ready")
            GetInletInFlow = .false.
        end if

    end function GetInletInFlow

    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::SetStormWaterModelFlow
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_SETSTORMWATERMODELFLOW"::SetStormWaterModelFlow
    !DEC$ ENDIF
    logical function SetStormWaterModelFlow(RunOffID, nComputePoints, overlandToSewerFlow)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: RunOffID
        integer                                     :: nComputePoints
        real(8), dimension(nComputePoints)          :: overlandToSewerFlow
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        integer                                     :: i, j, idx

        call Ready(RunOffID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            idx = 1
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%StormWaterInteraction(i, j) > AllmostZero) then
                    Me%StormWaterModelFlow(i, j) = overlandToSewerFlow(idx)
                    idx = idx + 1
                endif
            enddo
            enddo

            SetStormWaterModelFlow = .true.
        else 
            call PlaceErrorMessageOnStack("Runoff not ready")
            SetStormWaterModelFlow = .false.
        end if
           

    end function SetStormWaterModelFlow
        

#endif

end module ModuleRunOff
