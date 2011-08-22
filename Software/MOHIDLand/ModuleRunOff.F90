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
    !use ModuleTimeSerie
    use ModuleHDF5
    use ModuleFunctions         ,only : TimeToString, SetMatrixValue, ChangeSuffix
    use ModuleHorizontalGrid    ,only : GetHorizontalGridSize, GetHorizontalGrid,        &
                                        UnGetHorizontalGrid, WriteHorizontalGrid,        &
                                        GetGridCellArea
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
                                        GetChannelsID, GetChannelsStationName
    use ModuleDischarges        ,only : Construct_Discharges, GetDischargesNumber,       &
                                        GetDischargesGridLocalization,                   &
                                        GetDischargeWaterFlow, Kill_Discharges
                                        
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
    public  ::  GetFlowAtBoundary
    public  ::  GetFlowDischarge
    public  ::  GetRunoffWaterLevel
    public  ::  GetRunoffWaterColumn
    public  ::  GetRunoffWaterColumnOld
    public  ::  GetRunoffCenterVelocity
    public  ::  GetRunoffTotalStoredVolume
    public  ::  GetNextRunOffDT
    public  ::  SetBasinColumnToRunoff
    public  ::  UnGetRunOff
    

    !Modifier
    public  ::  ModifyRunOff
    private ::      LocalWaterColumn
    private ::      IntegrateFlow
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
    integer, parameter                              :: Manning_         = 1
    integer, parameter                              :: Chezy_           = 2

    integer, parameter                              :: UnitMax          = 80

    !Types---------------------------------------------------------------------
    type T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime                  => null()
         integer                                    :: NextOutPut
         logical                                    :: Yes = .false.
         type (T_Time), dimension(:), pointer       :: RestartOutTime
        logical                                     :: WriteRestartFile     = .false.
        logical                                     :: RestartOverwrite     = .false.
        integer                                     :: NextRestartOutput    = 1         
    end type T_OutPut


    type T_Files
        character(PathLength)                       :: DataFile
        character(PathLength)                       :: InitialFile
        character(PathLength)                       :: FinalFile
        character(PathLength)                       :: TransientHDF
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
        real                                        :: DT
    end type T_ExtVar


  
    type  T_RunOff
        integer                                     :: InstanceID
        character(len=StringLength)                 :: ModelName
        integer                                     :: ObjBasinGeometry         = 0
        integer                                     :: ObjTime                  = 0
        integer                                     :: ObjHorizontalGrid        = 0
        integer                                     :: ObjHorizontalMap         = 0
        integer                                     :: ObjGridData              = 0
        integer                                     :: ObjHDF5                  = 0
        integer                                     :: ObjDrainageNetwork       = 0
        integer                                     :: ObjDischarges            = 0
        type (T_OutPut   )                          :: OutPut
        type (T_ExtVar)                             :: ExtVar
        type (T_Files)                              :: Files
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime
        real(8), dimension(:,:), pointer            :: myWaterLevel             => null()
        real(8), dimension(:,:), pointer            :: myWaterColumn            => null()
        real(8), dimension(:,:), pointer            :: myWaterVolume            => null() 
        real(8), dimension(:,:), pointer            :: myWaterColumnOld         => null() !OldColumn from Basin
        real(8), dimension(:,:), pointer            :: myWaterVolumeOld         => null()
        real,    dimension(:,:), pointer            :: lFlowToChannels          => null() !Instantaneous Flow To Channels
        real,    dimension(:,:), pointer            :: iFlowToChannels          => null() !Integrated    Flow
        real,    dimension(:,:), pointer            :: lFlowBoundary            => null() !Instantaneous Flow to impose BC
        real,    dimension(:,:), pointer            :: iFlowBoundary            => null() !Integrated    Flow to impose BC
        real,    dimension(:,:), pointer            :: lFlowDischarge           => null() !Instantaneous Flow of discharges
        real,    dimension(:,:), pointer            :: iFlowDischarge           => null() !Integrated    Flow of discharges
        real,    dimension(:,:), pointer            :: lFlowToSewerSystem       => null() !Instantaneous Flow to sewer system
        real,    dimension(:,:), pointer            :: iFlowToSewerSystem       => null() !Integrated    Flow to sewer system
        real(8), dimension(:,:), pointer            :: lFlowX, lFlowY           => null() !Instantaneous OverLandFlow (LocalDT   )
        real(8), dimension(:,:), pointer            :: iFlowX, iFlowY           => null() !Integrated    OverLandFlow (AfterSumDT)
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
        real,    dimension(:,:), pointer            :: BuildingsHeight          => null() !Height of Buildings
        real, dimension(:,:), pointer               :: CenterFlowX, CenterFlowY
        real, dimension(:,:), pointer               :: CenterVelocityX, CenterVelocityY
        real, dimension(:,:), pointer               :: FlowModulus, VelocityModulus
        type(T_PropertyID)                          :: OverLandCoefficientID
        logical                                     :: StormWaterModel = .false.          !If connected to SWMM
        real,    dimension(:,:), pointer            :: StormWaterModelFlow      => null() !Flow from SWMM
        real                                        :: MinSlope
        logical                                     :: AdjustSlope
        logical                                     :: Stabilize
        logical                                     :: Discharges           = .false.
        logical                                     :: StormWaterDrainage   = .false.
        real                                        :: StormWaterInfiltrationVelocity  = 1.4e-5  !~50mm/h
        real                                        :: StormWaterFlowVelocity          = 0.2     !velocity in pipes
        logical                                     :: Buildings            = .false.
        real                                        :: StabilizeFactor
        integer                                     :: Routing              = Manning_
        logical                                     :: ImposeMaxVelocity    = .false.
        real                                        :: ImposedMaxVelocity   = 0.1
        logical                                     :: DynamicAdjustManning = .false.
        integer                                     :: LastGoodNiter        = 1
        real                                        :: MinimumWaterColumn
        real                                        :: NextDT               = null_real
        real                                        :: DTFactor
        logical                                     :: ImposeBoundaryValue  = .false.
        real                                        :: BoundaryValue
        real(8)                                     :: FlowAtBoundary       = 0.0
        integer                                     :: MaxIterations        = 10

        logical                                     :: WriteMaxFlowModulus  = .false.
        character(Pathlength)                       :: MaxFlowModulusFile
        real, dimension(:,:), pointer               :: MaxFlowModulus

        logical                                     :: WriteMaxWaterColumn  = .false.        
        character(Pathlength)                       :: MaxWaterColumnFile
        real, dimension(:,:), pointer               :: MaxWaterColumn
        
        logical                                     :: Continuous          = .false.
        logical                                     :: StopOnWrongDate     = .true.
        
        real(8)                                     :: TotalStoredVolume  = 0.
        real                                        :: InitialWaterColumn
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
                               InitialWaterColumn,                              &
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
        real, intent(OUT)                               :: InitialWaterColumn
        integer, optional, intent(OUT)                  :: STAT     

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
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR00'

            call ReadLockExternalVar (StaticOnly = .false.)


            !Gets the size of the grid
            call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                        Size     = Me%Size,                              &
                                        WorkSize = Me%WorkSize,                          &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunOff - ModuleRunOff - ERR01'
            
            call AllocateVariables
            
            call ReadDataFile

            call InitializeVariables

            InitialWaterColumn = Me%InitialWaterColumn

            call ConstructOverLandCoefficient

            !Checks if River Network is consistent with the one previously constructed
            if (DrainageNetworkID /= 0) then
                call CheckRiverNetWorkConsistency
            endif

            !Constructs Discharges
            if (Me%Discharges) then
                call Construct_Discharges(Me%ObjDischarges,                              &
                                          Me%ObjTime,                                    &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ConstructRunOff - ERR02' 
            endif
            
            !Constructs StormWaterDrainage
            if (Me%StormWaterDrainage .or. Me%StormWaterModel) then
                call ConstructStormWaterDrainage
            endif
            
            if (Me%OutPut%Yes) then
                call ConstructHDF5Output
            endif

            !Reads conditions from previous run
            if (Me%Continuous) call ReadInitialFile

            call CalculateTotalStoredVolume

            !Output Results
            if (Me%OutPut%Yes) then
                call ComputeCenterValues                   
                call RunOffOutput
            endif


            call ReadUnLockExternalVar (StaticOnly = .false.)

            !Returns ID
            RunOffID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleRunOff - ConstructRunOff - ERR99' 

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
        integer                                     :: ObjEnterData = 0
        integer                                     :: STAT_CALL
        type(T_PropertyID)                          :: OverLandCoefficientDeltaID
        type(T_PropertyID)                          :: StormWaterDrainageID
        type(T_PropertyID)                          :: BuildingsHeightID
        integer                                     :: iflag, ClientNumber
        logical                                     :: BlockFound
        integer                                     :: i, j

        !Reads the name of the data file from nomfich
        call ReadFileName ('RUNOFF_DATA', Me%Files%DataFile, "RunOff Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR01'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('RUNOFF_HDF', Me%Files%TransientHDF, "RunOff HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR02'

        call ReadFileName('RUNOFF_FIN', Me%Files%FinalFile,                              &
                           Message = "RunOff Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR05'

        !Constructs the DataFile
        call ConstructEnterData (ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR03'

        !Initial Water Column
        call GetData(Me%InitialWaterColumn,                                              &
                     ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'INITIAL_WATER_COLUMN',                              &
                     default      = 0.0,                                                 & 
                     ClientModule = 'ModuleRunOff',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR03.5'

         !Gets Minimum Slope 
        call GetData(Me%MinSlope,                                               &
                     ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'MIN_SLOPE',                                &
                     default      = 0.0,                                        &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04'

        if (Me%MinSlope < 0.0 .or. Me%MinSlope >= 1.) then
            write (*,*) 'Invalid Minimum Slope [MIN_SLOPE]'
            stop 'ReadDataFile - ModuleRunOff - ERR07'
        end if

        !Adjusts Slope according to
        !http://www.hkh-friend.net.np/rhdc/training/lectures/HEGGEN/Tc_3.pdf
        call GetData(Me%AdjustSlope,                                            &
                     ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'ADJUST_SLOPE',                             &
                     default      = .true.,                                     &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04'


        !Gets Routing method
        call GetData(Me%Routing,                                                &
                     ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'ROUTING',                                  &
                     default      = Manning_,                                   &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04'

        !if (Me%Routing /= Manning_ .and. Me%Routing /= Chezy_) then
        if (Me%Routing /= Manning_) then
            write (*,*) 'Invalid Routing Method [ROUTING]'
            stop 'ReadDataFile - ModuleRunOff - ERR07'
        end if
        
        !Gets if solution is limited by an maximum velocity
        call GetData(Me%ImposeMaxVelocity,                                      &
                     ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'IMPOSE_MAX_VELOCITY',                      &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04b'

        if (Me%ImposeMaxVelocity) then
        
            !Gets if solution is limited by an maximum velocity
            call GetData(Me%ImposedMaxVelocity,                                     &
                         ObjEnterData, iflag,                                       &
                         SearchType   = FromFile,                                   &
                         keyword      = 'MAX_VELOCITY',                             &
                         default      = 0.1,                                        &
                         ClientModule = 'ModuleRunOff',                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04c'
        
        endif


        !Gets if Manning Coeficient is increased with water depth
        call GetData(Me%DynamicAdjustManning,                                   &
                     ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'DYNAMIC_ADJUST_MANNING',                   &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04d'


        !Minimum Water Column for overland flow
        call GetData(Me%MinimumWaterColumn,                                              &
                     ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MIN_WATER_COLUMN',                                  &
!                     default      = 0.001,                                              &
                     ClientModule = 'ModuleRunOff',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04e'
        if (iflag == 0) then
            write(*,*)'MIN_WATER_COLUMN must be defined in module Runoff'
            stop 'ReadDataFile - ModuleRunOff - ERR07a'
        endif

        !Continuous Computation
        call GetData(Me%Continuous,                                                      &
                     ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONTINUOUS',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleRunoff',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04f'

        if (Me%Continuous) then
            call ReadFileName('RUNOFF_INI', Me%Files%InitialFile,                         &
                               Message = "Runoff Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04g'
        endif

        call GetData(Me%StopOnWrongDate,                                                 &
                     ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'STOP_ON_WRONG_DATE',                                &
                     default      = .true.,                                              &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04h'

        !Factor for DT Prediction
        call GetData(Me%DTFactor,                                           &
                     ObjEnterData, iflag,                                   &  
                     keyword      = 'DT_FACTOR',                            &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = 0.25,                                   &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR19'        

        !Stabilize Solution
        call GetData(Me%Stabilize,                                          &
                     ObjEnterData, iflag,                                   &  
                     keyword      = 'STABILIZE',                            &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .true.,                                 &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR19'        
        
        if (Me%Stabilize) then
            call GetData(Me%StabilizeFactor,                                    &
                         ObjEnterData, iflag,                                   &  
                         keyword      = 'STABILIZE_FACTOR',                     &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = 0.1,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR19'        
        endif
        
        !Impose Boundary Value
        call GetData(Me%ImposeBoundaryValue,                                    &
                     ObjEnterData, iflag,                                       &  
                     keyword      = 'IMPOSE_BOUNDARY_VALUE',                    &
                     ClientModule = 'ModuleRunOff',                             &
                     SearchType   = FromFile,                                   &
                     Default      = .false.,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR31'        
        
        if (Me%ImposeBoundaryValue) then
            call GetData(Me%BoundaryValue,                                      &
                         ObjEnterData, iflag,                                   &  
                         keyword      = 'BOUNDARY_VALUE',                       &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = 0.0,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR32'        
        endif
        
        !Discharges
        call GetData(Me%Discharges,                                         &
                     ObjEnterData, iflag,                                   &  
                     keyword      = 'DISCHARGES',                           &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR32a'        

        !Storm Water Drainage
        call GetData(Me%StormWaterDrainage,                                 &
                     ObjEnterData, iflag,                                   &  
                     keyword      = 'STORM_WATER',                          &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR33'

        if (Me%StormWaterDrainage) then
        
            !Storm Water Infiltration Velocity
            call GetData(Me%StormWaterInfiltrationVelocity,                     &
                         ObjEnterData, iflag,                                   &  
                         keyword      = 'STORM_WATER_INF_VEL',                  &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         default      = 1.4e-5,                                 & !~50mm/h
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR34'

            !Storm Water Transfer Coeficient
            call GetData(Me%StormWaterFlowVelocity,                             &
                         ObjEnterData, iflag,                                   &  
                         keyword      = 'STORM_WATER_FLOW_VEL',                 &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         default      = 0.2,                                    & !0.2m/s
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR34a'

        endif

        !If Buildings are to be simulated (flow ocuation in urban areas)
        call GetData(Me%Buildings,                                          &
                     ObjEnterData, iflag,                                   &  
                     keyword      = 'BUILDINGS',                            &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR35'
        
        !If Connected to a StormWater model
        call GetData(Me%StormWaterModel,                                    &
                     ObjEnterData, iflag,                                   &  
                     keyword      = 'STORM_WATER_MODEL_LINK',               &
                     ClientModule = 'ModuleRunOff',                         &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - ReadDataFile - ERR36'
        
                   

        !Gets Output Time 
        call GetOutPutTime(ObjEnterData,                                                 &
                           CurrentTime = Me%ExtVar%Now,                                  &
                           EndTime     = Me%EndTime,                                     &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%Yes,                                  &
                           STAT        = STAT_CALL)
        Me%OutPut%NextOutPut = 1    


        !Output for restart
        call GetOutPutTime(ObjEnterData,                                                &
                           CurrentTime  = Me%ExtVar%Now,                                &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoff - ERR40'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     ObjEnterData,                                                      &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleBasin',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleRunoff - ERR50'



        call RewindBuffer (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR05'

        !Gets Block for OverLand Coef
        call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                 &
                                    '<BeginOverLandCoefficient>',               &
                                    '<EndOverLandCoefficient>', BlockFound,     &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR05'
        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = Me%OverLandCoefficientID,     &
                                        EnterDataID      = ObjEnterData,                 &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%OverLandCoefficient,       &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR06'

            call KillFillMatrix(Me%OverLandCoefficientID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR07'


            !Check that manning values entered are not zero or negative
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    if (.not. Me%OverLandCoefficient(i,j) .gt. 0.0) then
                        write(*,*) 'Found Manning Overland coefficient zero or negative in input'
                        write(*,*) 'in cell', i, j
                        stop 'ReadDataFile - Module Runoff - ERR08'
                    endif
                
                
                endif
                
            enddo
            enddo

        else
            write(*,*)'Missing Block <BeginOverLandCoefficient> / <EndOverLandCoefficient>' 
            stop      'ReadDataFile - ModuleRunOff - ERR08'
        endif
        

        
        !Gets Block for OverLand Coef Difference 
        !To compute overland resistance in bottom for shear computation (erosion/deposition).
        !This process was created to remove from manning the resistance given by 
        !aerial vegetation parts that affect flow but do not affect bottom shear. Without that,
        !a manning increase (e.g. forestation) in one cell increases water depth (and reduces velocity)
        !but may increase shear stress (because water height increase is transformed in bottom resistance 
        !using manning - chezy see module runoff properties)
        call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                  &
                                    '<BeginOverLandCoefficientDelta>',           &
                                    '<EndOverLandCoefficientDelta>', BlockFound, &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR015'
        if (BlockFound) then
            call ConstructFillMatrix  ( PropertyID       = OverLandCoefficientDeltaID,   &
                                        EnterDataID      = ObjEnterData,                 &
                                        TimeID           = Me%ObjTime,                   &
                                        HorizontalGridID = Me%ObjHorizontalGrid,         &
                                        ExtractType      = FromBlock,                    &
                                        PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                        Matrix2D         = Me%OverLandCoefficientDelta,  &
                                        TypeZUV          = TypeZ_,                       &
                                        STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR020'

            call KillFillMatrix(OverLandCoefficientDeltaID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR030'

            !Check that final manning values are not zero or negative
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    
                    if (.not. (Me%OverLandCoefficient(i,j) - Me%OverLandCoefficientDelta(i,j)) .gt. 0.0) then
                        write(*,*) 'Manning Overland coefficient delta found zero or negative in input'
                        write(*,*) 'in cell', i, j
                        stop 'ReadDataFile - Module Runoff - ERR09'
                    endif
                
                
                endif
                
            enddo
            enddo

        else
            !Do not remove aerial vegetation effect from manning 
            Me%OverLandCoefficientDelta(:,:) = 0.0
        endif
        
        
        !Looks for StormWater DrainageCoef
        if (Me%StormWaterDrainage .or. Me%StormWaterModel) then

            allocate(Me%StormWaterDrainageCoef (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%StormWaterDrainageCoef  = null_real
            

            call RewindBuffer (ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR20'

            !Gets Flag with Sewer Points
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        '<BeginStormWaterDrainage>',                         &
                                        '<EndStormWaterDrainage>', BlockFound,               &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR21'

            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = StormWaterDrainageID,      &
                                            EnterDataID      = ObjEnterData,                 &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%StormWaterDrainageCoef,    &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR22'

                call KillFillMatrix(StormWaterDrainageID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR23'

            else
                write(*,*)'Missing Block <BeginStormWaterDrainage> / <EndStormWaterDrainage>' 
                stop      'ReadDataFile - ModuleRunOff - ERR08'
            endif
            
        endif

        if (Me%Buildings) then
        
            allocate(Me%BuildingsHeight(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            call RewindBuffer (ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR20'

            !Gets Flag with Sewer Points
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                        '<BeginBuildingsHeight>',                            &
                                        '<EndBuildingsHeight>', BlockFound,                  &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR21'

            if (BlockFound) then
                call ConstructFillMatrix  ( PropertyID       = BuildingsHeightID,            &
                                            EnterDataID      = ObjEnterData,                 &
                                            TimeID           = Me%ObjTime,                   &
                                            HorizontalGridID = Me%ObjHorizontalGrid,         &
                                            ExtractType      = FromBlock,                    &
                                            PointsToFill2D   = Me%ExtVar%BasinPoints,        &
                                            Matrix2D         = Me%BuildingsHeight,           &
                                            TypeZUV          = TypeZ_,                       &
                                            STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR22'

                call KillFillMatrix(BuildingsHeightID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR23'

            else
                write(*,*)'Missing Block <BeginBuildingsHeight> / <EndBuildingsHeight>' 
                stop      'ReadDataFile - ModuleRunOff - ERR08'
            endif
        
        endif
        

         !Write Max Flow Modulus File 
        call GetData(Me%WriteMaxFlowModulus,                                    &
                     ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'WRITE_MAX_FLOW_FILE',                      &
                     default      = .false.,                                    &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04'

        if(Me%WriteMaxFlowModulus) then
            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%MaxFlowModulusFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR02a'
            Me%MaxFlowModulusFile = trim(adjustl(Me%MaxFlowModulusFile))//"MaxRunOff.dat"
        end if

        !Write Max Channels Level  
        call GetData(Me%WriteMaxWaterColumn,                                    &
                     ObjEnterData, iflag,                                       &
                     SearchType   = FromFile,                                   &
                     keyword      = 'WRITE_MAX_WATER_COLUMN',                   &
                     default      = .true.,                                     &
                     ClientModule = 'ModuleRunOff',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR04'

        if(Me%WriteMaxWaterColumn) then
            !Gets the root path from the file nomfich.dat
            call ReadFileName("ROOT_SRT", Me%MaxWaterColumnFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR02a'
            Me%MaxWaterColumnFile = trim(adjustl(Me%MaxWaterColumnFile))//"MaxWaterColumn.dat"
        end if


        !Closes Data File
        call KillEnterData      (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR09'


    end subroutine ReadDataFile

    !--------------------------------------------------------------------------

    subroutine InitializeVariables

        !Arguments-------------------------------------------------------------
        
        !Local----------------------------------------------------------------- 
        integer                                           :: i,j
        !Begin-----------------------------------------------------------------       
        
        
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                Me%myWaterLevel(i, j)    = Me%InitialWaterColumn + Me%ExtVar%Topography(i, j)
                Me%MyWaterColumn(i,j)    = Me%InitialWaterColumn
                Me%MyWaterColumnOld(i,j) = Me%MyWaterColumn(i,j)
            endif

        enddo
        enddo
        
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
        

        allocate(Me%myWaterLevel         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterColumn        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolume        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterColumnOld     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%myWaterVolumeOld     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        Me%myWaterLevel            = null_real
        Me%myWaterColumn           = null_real
        Me%myWaterVolume           = null_real
        Me%myWaterColumnOld        = null_real
        Me%myWaterVolumeOld        = null_real

        allocate(Me%iFlowX               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%iFlowY               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowX               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%lFlowY               (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        
        Me%iFlowX               = 0.0
        Me%iFlowY               = 0.0
        Me%lFlowX               = 0.0
        Me%lFlowY               = 0.0

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
            Me%StormWaterCenterModulus = null_real

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
            Me%StormWaterModelFlow    = 0
            
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
        integer                                     :: STAT_CALL

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

        read(InitialFile)Me%MyWaterLevel
        
        if (Me%StormWaterDrainage) then
            read(InitialFile)Me%StormWaterVolume
        endif
        

   10   continue


        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleRunoff - ERR05'  
      
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

    subroutine GetFlowAtBoundary (ObjRunOffID, FlowAtBoundary, STAT)

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

    end subroutine GetFlowAtBoundary
    
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

            DT        = Me%NextDT

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
                    Me%myWaterVolume(i, j) = WaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
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
        real                                        :: LocalDT, SumDT
        logical                                     :: Restart
        integer                                     :: Niter, iter
        logical, save                               :: XFirst = .true.

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunOffID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (MonitorPerformance) call StartWatch ("ModuleRunOff", "ModifyRunOff")

            !Time Stuff
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunOff - ModuleRunOff - ERR01'

            !Stores initial values = from basin
            Me%myWaterColumnOld = Me%myWaterColumn

            !Adds Flow from SEWER OverFlow to Water Column OLD
            if (Me%StormWaterModel) then
                call ReadLockExternalVar   (StaticOnly = .true.)
                call AddFlowFromStormWaterModel
                call ReadUnLockExternalVar (StaticOnly = .true.)
            endif

            
            Restart     = .true.
            do while (Restart)

                !Calculates local Watercolumn
                call ReadLockExternalVar   (StaticOnly = .true.)
                call LocalWaterColumn      (Me%myWaterColumnOld)
                call ReadUnLockExternalVar (StaticOnly = .true.)


                SumDT       = 0.0
                Restart     = .false.
                Niter       = max(Me%LastGoodNiter - 1, 1)
                iter        = 1
                LocalDT     = Me%ExtVar%DT / Niter
                Me%iFlowX          = 0.0
                Me%iFlowY          = 0.0
                Me%iFlowToChannels = 0.0
                Me%iFlowBoundary   = 0.0
doIter:         do while (iter <= Niter)

                    !Gets ExternalVars
                    call ReadLockExternalVar (StaticOnly = .false.)

                    !Stores WaterVolume for convergence test
                    Me%myWaterVolumeOld = Me%myWaterVolume

                    !Inputs Water from discharges
                    if (Me%Discharges) then
                        call ModifyWaterDischarges  (LocalDT)                
                    endif

                    !Calculates Flow Direction
                    XFirst = .not. XFirst
                    if (XFirst) then
                        call DirectX(LocalDT)
                        call DirectY(LocalDT)
                    else
                        call DirectY(LocalDT)
                        call DirectX(LocalDT)
                    endif
                    
                    !Interaction with channels
                    if (Me%ObjDrainageNetwork /= 0) then
                        call FlowIntoChannels       (LocalDT)
                    endif

                    !Boundary Condition
                    if (Me%ImposeBoundaryValue) then
                        call ImposeBoundaryValue    (LocalDT)
                    endif

                    
                    call CheckStability(Restart, Niter)
                    
                    call ReadUnLockExternalVar (StaticOnly = .false.)
                    
                    if (Restart) then
                        Me%LastGoodNiter   = Me%LastGoodNiter + 2
                        exit doIter
                    endif

                    call IntegrateFlow     (LocalDT, SumDT)  
                        
                    SumDT = SumDT + LocalDT
                    iter  = iter  + 1
                    
                    
                    
                enddo doIter
            
                
                
            enddo
            

            Me%LastGoodNiter = Niter

            !Gets ExternalVars
            call ReadLockExternalVar (StaticOnly = .false.)

                    
            !StormWater Drainage
            if (Me%StormWaterDrainage) then
                call StormWaterDrainage
            endif


            !Calculates flow from channels to land
            if (Me%ObjDrainageNetwork /= 0) then
                call FlowFromChannels 
            endif


            !Calculates center flow and velocities (for output and next DT)
            call ComputeCenterValues
            
            !Sets Next DT
            !Me%NextDT   = LocalDT * 1.5 !/ max(Me%LastGoodNiter - 1, 1)

            if (Niter <= 5) then
                Me%NextDT = Me%ExtVar%DT * 1.50
            else if (Niter <= 10) then
                Me%NextDT = Me%ExtVar%DT
            else
                Me%NextDT = Me%ExtVar%DT / 2.0
            endif

            !Output Results
            if (Me%OutPut%Yes) then                   
                call RunOffOutput
            endif

            if (Me%ObjDrainageNetwork /= 0 .and. Me%WriteMaxWaterColumn) &
                call OutputOutputMaxWaterColumn

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
        real                                    :: Flow    
        integer                                 :: STAT_CALL

        !Sets to 0
        Me%lFlowDischarge = 0.0

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

                Me%lFlowDischarge(i, j)     = Me%lFlowDischarge(i, j) + Flow

                !Updates Water Volume
                Me%myWaterVolume(i, j)      = Me%myWaterVolume(i, j) + Flow * LocalDT

                !Updates Water Column
                Me%myWaterColumn  (i, j)    = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                !Updates Water Level
                Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j) + Me%ExtVar%Topography(i, j)

                !if (Me%CheckMass) Me%TotalInputVolume = Me%TotalInputVolume + Me%DischargesFlow(iDis) * LocalDT
 
            endif
           
        enddo

    end subroutine ModifyWaterDischarges    
    
    !--------------------------------------------------------------------------
    
    subroutine DirectX(LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: MaxBottom
        real(8)                                     :: WCL, WCR, WCA
        real                                        :: Slope, dVol, MaxFlow
        real                                        :: level_left, level_right

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint .and. Me%ExtVar%BasinPoints(i, j-1) == BasinPoint) then
            
                !Maximum Bottom Level
                MaxBottom = max(Me%ExtVar%Topography(i, j-1), Me%ExtVar%Topography(i, j))
                
                !Water Column Left
                WCL       = max(Me%myWaterLevel(i, j-1) - MaxBottom, dble(0.0))
            
                !Water Column Right
                WCR       = max(Me%myWaterLevel(i, j  ) - MaxBottom, dble(0.0))
                
                !Average Water Column
                WCA       = (WCL + WCR) / 2.0
                
                if (WCA > Me%MinimumWaterColumn) then
                
                    !Adds to the final level the height of the buidings, if any
                    if (Me%Buildings) then
                        level_left  = Me%myWaterLevel(i, j-1) + Me%BuildingsHeight(i, j-1)
                        level_right = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_left  = Me%myWaterLevel(i, j-1)
                        level_right = Me%myWaterLevel(i, j)
                    endif
                    
                
                    if (level_left > level_right) then
                
                        !Slope
                        Slope           = AdjustSlope((level_left - level_right) / Me%ExtVar%DZX(i, j-1))

                        !Flow
                        !Me%lFlowX(i, j) = FlowRouting(Me%Routing, WCA, Me%ExtVar%DYY(i, j), Slope, Me%OverlandCoefficientX(i,j))
                        Me%lFlowX(i, j) = FlowRouting(Me%Routing, WCL, Me%ExtVar%DYY(i, j), Slope, Me%OverlandCoefficientX(i,j))
                        
                        !Estimate max flow (considering both cell the same area)
                        MaxFlow = 0.5 * (level_left - level_right) * Me%ExtVar%GridCellArea(i, j) / LocalDT
                        
                        !Test for to high flow
                        if (Me%lFlowX(i, j) > MaxFlow) then
                            Me%lFlowX(i, j) = MaxFlow
                        endif
                    
                    else
                    
                        !Slope
                        Slope           = AdjustSlope((level_right - level_left) / Me%ExtVar%DZX(i, j-1))
                    
                        !Flow
                        !Me%lFlowX(i, j) = -1.0 * FlowRouting(Me%Routing, WCA, Me%ExtVar%DYY(i, j), 
                        !Slope, Me%OverlandCoefficientX(i,j))

                        Me%lFlowX(i, j) = -1.0 * FlowRouting(Me%Routing, WCR, Me%ExtVar%DYY(i, j), Slope,        &
                                                             Me%OverlandCoefficientX(i,j))
                    
                        !Estimate max flow (considering both cell the same area)
                        MaxFlow = - 0.5 * (level_right - level_left) * Me%ExtVar%GridCellArea(i, j) / LocalDT

                        if (Me%lFlowX(i, j) < MaxFlow) then
                            Me%lFlowX(i, j) = MaxFlow
                        endif
                   
                    endif
                    
                    
                    !dVol
                    dVol = Me%lFlowX(i, j) * LocalDT
                    
                    !Updates Water Volume
                    Me%myWaterVolume (i, j-1) = Me%myWaterVolume (i, j-1) - dVol 
                    Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   + dVol 
                    
                   
                    !Updates Water Column
                    Me%myWaterColumn  (i, j-1) = Me%myWaterVolume (i, j-1) / Me%ExtVar%GridCellArea(i, j-1)
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    !Updates Water Level
                    Me%myWaterLevel (i, j-1)   = Me%myWaterColumn (i, j-1) + Me%ExtVar%Topography(i, j-1)
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    
                else
                
                    Me%lFlowX(i, j) = 0.0

                endif

            
            endif
        enddo
        enddo        
        
    end subroutine DirectX
    
    !--------------------------------------------------------------------------
    
    subroutine DirectY(LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: MaxBottom
        real(8)                                     :: WCL, WCR, WCA
        real                                        :: Slope, dVol, MaxFlow
        real                                        :: level_bottom, level_top

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        


        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint .and. Me%ExtVar%BasinPoints(i-1, j) == BasinPoint) then
            
                !Maximum Bottom Level
                MaxBottom = max(Me%ExtVar%Topography(i-1, j), Me%ExtVar%Topography(i, j))
                
                !Water Column Bottom
                WCL       = max(Me%myWaterLevel(i-1, j) - MaxBottom, dble(0.0))
            
                !Water Column Right
                WCR       = max(Me%myWaterLevel(i, j  ) - MaxBottom, dble(0.0))
                
                !Average Water Column
                WCA       = (WCL + WCR) / 2.0
                
                if (WCA > Me%MinimumWaterColumn) then
                
                    !Adds to the final level the height of the buidings, if any
                    if (Me%Buildings) then
                        level_bottom = Me%myWaterLevel(i-1, j) + Me%BuildingsHeight(i-1, j)
                        level_top    = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_bottom = Me%myWaterLevel(i-1, j)
                        level_top    = Me%myWaterLevel(i, j)
                    endif
                
                    if (level_bottom > level_top) then
                
                        !Slope
                        Slope           = AdjustSlope((level_bottom - level_top) / Me%ExtVar%DZY(i-1, j))

                        !Flow
                        !Me%lFlowY(i, j) = FlowRouting(Me%Routing, WCA, Me%ExtVar%DXX(i, j), Slope, Me%OverlandCoefficientY(i,j))
                        Me%lFlowY(i, j) = FlowRouting(Me%Routing, WCL, Me%ExtVar%DXX(i, j), Slope, Me%OverlandCoefficientY(i,j))
                 
                        !Estimate max flow (considering both cell the same area)
                        MaxFlow = 0.5 * (level_bottom - level_top) * Me%ExtVar%GridCellArea(i, j) / LocalDT

                        if (Me%lFlowY(i, j) > MaxFlow) then
                            Me%lFlowY(i, j) = MaxFlow
                        endif
                        
                     else
                    
                        !Slope
                        Slope           = AdjustSlope((level_top - level_bottom) / Me%ExtVar%DZY(i-1, j))
                    
                        !Flow
                        !Me%lFlowY(i, j) = -1.0 * FlowRouting(Me%Routing, WCA, Me%ExtVar%DXX(i, j), 
                        !Slope, Me%OverlandCoefficientY(i,j))

                        Me%lFlowY(i, j) = -1.0 * FlowRouting(Me%Routing, WCR, Me%ExtVar%DXX(i, j), Slope,      &
                                                             Me%OverlandCoefficientY(i,j))
                    
                        !Estimate max flow (considering both cell the same area)
                        MaxFlow = - 0.5 * (level_top - level_bottom) * Me%ExtVar%GridCellArea(i, j) / LocalDT

                        if (Me%lFlowY(i, j) < MaxFlow) then
                            Me%lFlowY(i, j) = MaxFlow
                        endif

                    endif
                    
                    
                    !dVol
                    dVol = Me%lFlowY(i, j) * LocalDT
                    
                    !Updates Water Volume
                    Me%myWaterVolume (i-1, j) = Me%myWaterVolume (i-1, j) - dVol 
                    Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   + dVol 

                    !Updates Water Column
                    Me%myWaterColumn  (i-1, j) = Me%myWaterVolume (i-1, j) / Me%ExtVar%GridCellArea(i-1, j)

                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    !Updates Water Level
                    Me%myWaterLevel (i-1, j)   = Me%myWaterColumn (i-1, j) + Me%ExtVar%Topography(i-1, j)
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    
                else
                
                    Me%lFlowY(i, j) = 0.0
                
                endif

            
            endif
        enddo
        enddo        
        
    end subroutine DirectY    

    !--------------------------------------------------------------------------
    
    subroutine StormWaterDrainage
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: AvaliableVolume, FlowVolume, InfiltrationVolume
        real                                        :: dVol
        integer, dimension(:, :), pointer           :: DrainageDirection
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 

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
        
            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%myWaterColumn(i, j) > AllmostZero) then

                !Volume which is avaliable on the cell and might infiltrate
                !m3                  
                AvaliableVolume             = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
                
                !Volume which can infiltrate at avaliable area during time step at max infiltration velocity
                !m3 
                FlowVolume                  = Me%StormWaterDrainageCoef(i, j) * Me%ExtVar%GridCellArea(i, j) * &
                                              Me%StormWaterInfiltrationVelocity * Me%ExtVar%DT

                !Volume which will be removed from overland flow into the StormWater System
                !m3
                InfiltrationVolume          = Min(AvaliableVolume, FlowVolume)

                !New StormWater Volume at point
                Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j) + InfiltrationVolume

                !New Volume of Overland volume
                Me%myWaterVolume (i, j)     = Me%myWaterVolume (i, j)   - InfiltrationVolume
                
                !New WaterColumn
                Me%myWaterColumn (i, j)     = Me%myWaterVolume (i, j) / Me%ExtVar%GridCellArea(i, j)

                !New Level
                Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j) + Me%ExtVar%Topography(i, j)

            endif

        enddo
        enddo
        
        
        !Flow along X
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i, j-1) > AllmostZero) then
            
                !Flow from the left to the right
                if (Me%ExtVar%Topography(i, j-1) > Me%ExtVar%Topography(i, j)) then
                    Me%StormWaterFlowX(i, j) =      Me%StormWaterVolume(i, j-1) / Me%ExtVar%DZX(i, j-1) * Me%StormWaterFlowVelocity
                else
                    Me%StormWaterFlowX(i, j) = -1 * Me%StormWaterVolume(i, j)   / Me%ExtVar%DZX(i, j-1) * Me%StormWaterFlowVelocity
                endif

            else
            
                Me%StormWaterFlowX(i, j) = 0
            
            endif
            
            
            
        enddo
        enddo
        

        !Flow along Y
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i-1, j) > AllmostZero) then
            
                if (Me%ExtVar%Topography(i-1, j) > Me%ExtVar%Topography(i, j)) then
                    Me%StormWaterFlowY(i, j) =      Me%StormWaterVolume(i-1, j) / Me%ExtVar%DZY(i-1, j) * Me%StormWaterFlowVelocity
                else
                    Me%StormWaterFlowY(i, j) = -1 * Me%StormWaterVolume(i, j)   / Me%ExtVar%DZX(i-1, j) * Me%StormWaterFlowVelocity
                endif
            
            else
            
                Me%StormWaterFlowY(i, j) = 0
            
            endif
            
        enddo
        enddo
            
        !Updates volumes
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i, j-1) > AllmostZero) then
                
                if (Me%StormWaterFlowX(i, j) > 0) then
                
                    dVol = min(Me%StormWaterFlowX(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i, j-1))
                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   + dVol
                    Me%StormWaterVolume(i, j-1) = Me%StormWaterVolume(i, j-1) - dVol
                    
                else
                
                    dVol = min(-Me%StormWaterFlowX(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i, j))
                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   - dVol
                    Me%StormWaterVolume(i, j-1) = Me%StormWaterVolume(i, j-1) + dVol
                
                endif
                
            endif
        enddo
        enddo
        
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%StormWaterDrainageCoef(i, j) > AllmostZero .and. Me%StormWaterDrainageCoef(i-1, j) > AllmostZero) then
                
                if (Me%StormWaterFlowY(i, j) > 0) then
                
                    dVol = min(Me%StormWaterFlowY(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i-1, j))
                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   + dVol
                    Me%StormWaterVolume(i-1, j) = Me%StormWaterVolume(i-1, j) - dVol
                    
                else
                
                    dVol = min(-Me%StormWaterFlowY(i, j) * Me%ExtVar%DT, Me%StormWaterVolume(i, j))
                    Me%StormWaterVolume(i, j)   = Me%StormWaterVolume(i, j)   - dVol
                    Me%StormWaterVolume(i-1, j) = Me%StormWaterVolume(i-1, j) + dVol
                
                endif
            endif
        enddo
        enddo
                    
      
        !Routes water from StormWater Drainage System to river channels
        if (Me%ObjDrainageNetwork /= 0) then 

            call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR01'     

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
            if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'
            
        endif

        
        call UnGetBasin (Me%ObjBasinGeometry, DrainageDirection, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunOff - StormWaterDrainage - ERR02'
                
        
    end subroutine StormWaterDrainage

    !--------------------------------------------------------------------------
    
    subroutine AddFlowFromStormWaterModel

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                Me%myWaterColumnOld   (i, j)  = Me%myWaterColumnOld   (i, j) +          &
                                                Me%StormWaterModelFlow(i, j) *          &
                                                Me%ExtVar%DT /                          &
                                                Me%ExtVar%GridCellArea(i, j)
                                             
                if (Me%myWaterColumnOld(i, j) < 0.0) then
                    !Write(*,*)'Negative WC: ', i, j, Me%myWaterColumnOld(i, j)
                    Me%myWaterColumnOld   (i, j)  = 0.0
                endif
                
            endif

        enddo
        enddo  
                  
    
    end subroutine AddFlowFromStormWaterModel
    
    !--------------------------------------------------------------------------
    
    subroutine FlowIntoChannels(LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB, STAT_CALL
        real                                        :: DifLevel
        real                                        :: Slope, AverageCellLength, dVol, MaxFlow
        real   , dimension(:, :), pointer           :: ChannelsWaterLevel 
        real   , dimension(:, :), pointer           :: ChannelsNodeLength 


        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR01'     

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR04'


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) then

                !Checks for Flow from Land -> Channel
                AverageCellLength  = ( Me%ExtVar%DUX (i, j) + Me%ExtVar%DVY (i, j) ) / 2.0

            
                if (ChannelsWaterLevel (i, j) < Me%myWaterLevel(i, j) .and. Me%myWaterColumn(i, j) > Me%MinimumWaterColumn) then

                    if (ChannelsWaterLevel (i, j) > Me%ExtVar%Topography(i, j)) then
                        DifLevel           = Me%myWaterLevel(i, j) - ChannelsWaterLevel (i, j)
                    else
                        DifLevel           = Me%myWaterColumn(i, j)
                    endif

                    Slope              = AdjustSlope(DifLevel / (AverageCellLength / 4.0))                

                    Me%lFlowToChannels(i, j) = FlowRouting (Me%Routing,                     &
                                                            real(DifLevel, 8),              &
                                                            ChannelsNodeLength(i, j),       &
                                                            Slope,                          &
                                                            Me%OverlandCoefficient(i,j))
                
                    MaxFlow = 0.5 * (DifLevel) * Me%ExtVar%GridCellArea(i, j) / LocalDT
                
                    if (Me%lFlowToChannels(i, j) > MaxFlow) then
                        Me%lFlowToChannels(i, j) = MaxFlow
                    endif
                
                
                    !dVol
                    dVol = Me%lFlowToChannels(i, j) * LocalDT
                    
                    !Updates Water Volume
                    Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   - dVol 
                    
                    !Updates Water Column
                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    !Updates Water Level
                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
                    
!                else if (ChannelsWaterLevel (i, j) > Me%ExtVar%Topography(i, j)) then
!                
!
!                    DifLevel           = ChannelsWaterLevel (i, j) - Me%myWaterLevel (i, j)
!                    
!                    !Checks for Flow Channel - Land
!                    Slope              = AdjustSlope(DifLevel / (AverageCellLength / 4.0))
!                    
!
!                    
!                    Me%lFlowToChannels(i, j) = -1.0* FlowRouting (Me%Routing,                     &
!                                                                  DifLevel,                       &
!                                                                  ChannelsNodeLength(i, j),       &
!                                                                  Slope,                          &
!                                                                  Me%OverlandCoefficient(i,j))
!                    
!                    !dVol
!                    dVol = Me%lFlowToChannels(i, j) * LocalDT
!                    
!                    !Updates Water Volume
!                    Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   - dVol 
!                    
!                    !Updates Water Column
!                    Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)
!                    
!                    !Updates Water Level
!                    Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
!                   
!                
                else
                
                    Me%lFlowToChannels(i, j) = 0.0
                
                endif

            
            endif

        enddo
        enddo        
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR05'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FlowIntoChannels - ModuleRunOff - ERR08'

        
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
        real                                        :: x1, x2


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


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) then

                if (ChannelsWaterLevel (i, j) > Me%myWaterLevel(i, j)) then
                
                    ChannelHeight = Me%ExtVar%Topography(i, j) - ChannelsBottomLevel(i, j)                                       
                    !ChannelSlope  = (ChannelsTopWidth(i, j) - ChannelsBottomWidth(i, j)) / ChannelHeight
                    !ChannelSurfaceWidth = ChannelsBottomWidth(i,j) + 2.* ChannelSlope * ChannelHeight
                    
                    !Water Column in River above Topo
                    WCR           = ChannelsWaterLevel (i, j) - Me%ExtVar%Topography(i, j)
                    
                    !Volume above Topography
                    VolExcess    = ChannelsBankSlope(i,j) * WCR * WCR * ChannelsNodeLength(i, j)       &
                                    + WCR * ChannelsSurfaceWidth(i, j) * ChannelsNodeLength(i, j) +    &
                                    Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)

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
                    Me%iFlowToChannels(i, j)    = Me%iFlowToChannels(i, j) -dVol / Me%ExtVar%DT     
            
                    Me%myWaterVolume (i, j)     = Me%myWaterVolume (i, j) + dVol 
                    
                    Me%myWaterColumn  (i, j)    = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                    Me%myWaterLevel (i, j)      = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)

                
                endif

            
            endif

        enddo
        enddo        
        
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

        
    
    end subroutine FlowFromChannels
    
    !--------------------------------------------------------------------------

    subroutine CheckStability (Restart, Niter)

        !Arguments-------------------------------------------------------------
        logical                                     :: Restart
        integer                                     :: Niter

        !Local-----------------------------------------------------------------
        integer                                     :: i, j

        !Begin-----------------------------------------------------------------
        
        !Verifies negative volumes
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                if (Me%myWaterVolume (i, j) < -1.0 * AllmostZero) then
                    Restart = .true.
                    return
                else if (Me%myWaterVolume (i, j) < 0.0) then
                    Me%myWaterVolume (i, j) = 0.0
                endif
            endif
        enddo
        enddo
        
        
        !Verifies stabilize criteria
        if (Me%Stabilize .and. Niter < Me%MaxIterations) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%myWaterVolumeOld(i, j) / Me%ExtVar%GridCellArea(i, j) > 0.01) then
                    if (abs(Me%myWaterVolume(i, j) - Me%myWaterVolumeOld(i, j)) / Me%myWaterVolumeOld(i, j)     &
                         > Me%StabilizeFactor) then
                        Restart = .true.
                        return
                    endif
                endif
            enddo
            enddo
        endif
        
        
    end subroutine CheckStability
 
    !--------------------------------------------------------------------------

    subroutine LocalWaterColumn (WaterColumn)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), pointer              :: WaterColumn

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Estimates Flow along X Direction
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                Me%myWaterColumn(i, j) = WaterColumn(i, j)

                Me%myWaterLevel (i, j) = Me%myWaterColumn(i, j) + Me%ExtVar%Topography(i, j)
                Me%myWaterVolume(i, j) = Me%myWaterColumn(i, j) * Me%ExtVar%GridCellArea(i, j)
            endif
        enddo
        enddo
        
    end subroutine LocalWaterColumn            

    !--------------------------------------------------------------------------

    subroutine IntegrateFlow (LocalDT, SumDT)

        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT, SumDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j

        !Integrates along X Directions
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%iFlowX(i, j) = (Me%iFlowX(i, j) * SumDT + Me%lFlowX(i, j) * LocalDT) / &
                              (SumDT + LocalDT)
        enddo
        enddo

        !Integrates along Y Directions
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%iFlowY(i, j) = (Me%iFlowY(i, j) * SumDT + Me%lFlowY(i, j) * LocalDT) / &
                              (SumDT + LocalDT)
        enddo
        enddo

        !Integrates Flow to Channels
        if (Me%ObjDrainageNetwork /= 0) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%iFlowToChannels(i, j) = (Me%iFlowToChannels(i, j) * SumDT + Me%lFlowToChannels(i, j) * LocalDT) / &
                                           (SumDT + LocalDT)
            enddo
            enddo
        endif
        
        !Integrates Flow At boundary
        if (Me%ImposeBoundaryValue) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%iFlowBoundary(i, j) = (Me%iFlowBoundary(i, j) * SumDT + Me%lFlowBoundary(i, j) * LocalDT) / &
                                         (SumDT + LocalDT)
            enddo
            enddo
        endif

        !Integrates Flow Discharges
        if (Me%Discharges) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%iFlowDischarge(i, j) = (Me%iFlowDischarge(i, j) * SumDT + Me%lFlowDischarge(i, j) * LocalDT) / &
                                          (SumDT + LocalDT)
            enddo
            enddo
        endif


    end subroutine IntegrateFlow

    !--------------------------------------------------------------------------

    subroutine ImposeBoundaryValue (LocalDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: dh, dVOl

        !Empiric estimation of the water flow into the channels
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        !Default is zero
        Me%lFlowBoundary = 0.0
        
        !Sets Boundary values
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BoundaryPoints2D(i, j) == 1 .and. (j == Me%WorkSize%JLB .or. j == Me%WorkSize%JUB    &
                 .or. i == Me%WorkSize%ILB .or. i == Me%WorkSize%IUB)) then

                !Necessary Variation in height            
                dh = Me%myWaterColumn (i, j) - Me%BoundaryValue
                
                !Flow to set cell equal to Boundary Value
                !m3/s                         = m  * 
                Me%lFlowBoundary(i, j) = dh * Me%ExtVar%GridCellArea(i, j) / LocalDT
           
                !dVol
                dVol = Me%lFlowBoundary(i, j) * LocalDT
                    
                !Updates Water Volume
                Me%myWaterVolume (i, j)   = Me%myWaterVolume (i, j)   - dVol 
                    
                !Updates Water Column
                Me%myWaterColumn  (i, j)   = Me%myWaterVolume (i, j)   / Me%ExtVar%GridCellArea(i, j)

                !Updates Water Level
                Me%myWaterLevel (i, j)     = Me%myWaterColumn (i, j)   + Me%ExtVar%Topography(i, j)
           
           
            endif
        enddo
        enddo
    
    end subroutine
    
    !--------------------------------------------------------------------------
    
    subroutine ComputeCenterValues 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
       
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
            
           
        Me%CenterFlowX     = 0.0
        Me%CenterFlowY     = 0.0
        Me%FlowModulus     = 0.0
        Me%CenterVelocityX = 0.0
        Me%CenterVelocityY = 0.0
        Me%VelocityModulus = 0.0

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
                end if

                if(Me%WriteMaxFlowModulus) then
                    if (Me%FlowModulus(i, j) > Me%MaxFlowModulus(i, j)) then
                        Me%MaxFlowModulus(i, j) = Me%FlowModulus(i, j)
                    end if
                end if

            endif

        enddo
        enddo
        
        if (Me%StormWaterDrainage) then
        
            Me%StormWaterCenterFlowX    = 0.0
            Me%StormWaterCenterFlowY    = 0.0
            Me%StormWaterCenterModulus  = 0.0
        
            do j = JLB, JUB
            do i = ILB, IUB
                    
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                    Me%StormWaterCenterFlowX(i, j)   = (Me%StormWaterFlowX(i, j) + Me%StormWaterFlowX(i, j+1)) / 2.0
                    Me%StormWaterCenterFlowY(i, j)   = (Me%StormWaterFlowY(i, j) + Me%StormWaterFlowY(i+1, j)) / 2.0
                    Me%StormWaterCenterModulus(i, j) = sqrt (Me%StormWaterCenterFlowX(i, j)**2. + &
                                                             Me%StormWaterCenterFlowY(i, j)**2.)

                endif

            enddo
            enddo
        
        endif
        
        
    end subroutine ComputeCenterValues 

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
                                  "/Results/flow/flow X",                           &
                                  "flow X",                                         &   
                                  "m3/s",                                           &
                                  Array2D      = Me%CenterFlowX,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'

            
            !Writes Flow Y
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/flow/flow Y",                           &
                                  "flow Y",                                         &   
                                  "m3/s",                                           &
                                  Array2D      = Me%CenterFlowY,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

             !Writes Flow Modulus
            call HDF5WriteData   (Me%ObjHDF5,                                       &
                                  "/Results/flow/"//trim(GetPropertyName (FlowModulus_)),&
                                  trim(GetPropertyName (FlowModulus_)),             &   
                                  "m3/s",                                           &
                                  Array2D      = Me%FlowModulus,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,              &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

             !Writes Velocity X 
            call HDF5WriteData   (Me%ObjHDF5,                                          &
                                  "/Results/velocity/"//trim(GetPropertyName (VelocityU_)),     &
                                  trim(GetPropertyName (VelocityU_)),                  &
                                  "m/s",                                               &
                                  Array2D      = Me%CenterVelocityX,                   &
                                  OutputNumber = Me%OutPut%NextOutPut,                 &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

             !Writes Velocity Y 
            call HDF5WriteData   (Me%ObjHDF5,                                          &
                                  "/Results/velocity/"//trim(GetPropertyName (VelocityV_)),     &
                                  trim(GetPropertyName (VelocityV_)),                  &
                                  "m/s",                                               &
                                  Array2D      = Me%CenterVelocityY,                   &
                                  OutputNumber = Me%OutPut%NextOutPut,                 &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

            !Writes Velocity Modulus
            call HDF5WriteData   (Me%ObjHDF5,                                                &
                                  "/Results/velocity/"//trim(GetPropertyName (VelocityModulus_)),     &
                                  trim(GetPropertyName (VelocityModulus_)),                  &
                                  "m/s",                                                     &
                                  Array2D      = Me%VelocityModulus,                         &
                                  OutputNumber = Me%OutPut%NextOutPut,                       &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'

            !Writes Storm Water Volume of each Cell
            if (Me%StormWaterDrainage) then
                call HDF5WriteData   (Me%ObjHDF5, "//Results/storm water/water volume",     &
                                      "storm water volume", "m3",                           &
                                      Array2D      = Me%StormWaterVolume,                   &
                                      OutputNumber = Me%OutPut%NextOutPut,                  &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleRunOff - ERR080'
                
                
                !Writes Flow X
                call HDF5WriteData   (Me%ObjHDF5,                                       &
                                      "//Results/storm water/storm flow X",             &
                                      "storm water flow X",                             &   
                                      "m3/s",                                           &
                                      Array2D      = Me%StormWaterCenterFlowX,          &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'

                
                !Writes SW Flow Y
                call HDF5WriteData   (Me%ObjHDF5,                                       &
                                      "//Results/storm water/storm flow Y",             &
                                      "storm water flow Y",                             &   
                                      "m3/s",                                           &
                                      Array2D      = Me%StormWaterCenterFlowY,          &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'
                

               
                !Writes SW Modulus
                call HDF5WriteData   (Me%ObjHDF5,                                       &
                                      "//Results/storm water/storm flow modulus",       &
                                      "storm water flow modulus",                       &   
                                      "m3/s",                                           &
                                      Array2D      = Me%StormWaterCenterModulus,        &
                                      OutputNumber = Me%OutPut%NextOutPut,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'
                
                
            endif
            
            if (Me%StormWaterModel) then
            
                call HDF5WriteData   (Me%ObjHDF5, "//Results/storm water/water flow",       &
                                      "storm water flow", "m3",                             &
                                      Array2D      = Me%StormWaterModelFlow,                &
                                      OutputNumber = Me%OutPut%NextOutPut,                  &
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


    subroutine OutputOutputMaxWaterColumn

        !Locals----------------------------------------------------------------
        integer                                 :: ILB,IUB, JLB, JUB, i, j
        integer                                 :: STAT_CALL
        real, dimension(:,:), pointer           :: ChannelsWaterLevel
        
        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputOutputMaxWaterColumn - ModuleRunOff - ERR01'     

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
                
                !Water Column of River Network
                if (Me%ExtVar%RiverPoints(i, j) == BasinPoint) then
                    if (ChannelsWaterLevel(i, j) - Me%ExtVar%Topography(i, j) > Me%MaxWaterColumn(i, j)) then
                        Me%MaxWaterColumn(i, j) = ChannelsWaterLevel(i, j) - Me%ExtVar%Topography(i, j)
                    endif
                endif

            endif

        enddo
        enddo

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutputOutputMaxWaterColumn - ModuleRunOff - ERR05'

    end subroutine OutputOutputMaxWaterColumn

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

    !--------------------------------------------------------------------------

    real function FlowRouting (RoutingMethod, WaterColumn, Width, Slope, OverLandCoefficient)

        !Arguments--------------------------------------------------------------
        integer                                 :: RoutingMethod
        real(8)                                 :: WaterColumn
        real                                    :: Width, Slope
        real                                    :: OverLandCoefficient
        real                                    :: VelAux, Area, Coef

        if (WaterColumn > Me%MinimumWaterColumn) then

            if (RoutingMethod == Manning_) then
            
                if (Me%DynamicAdjustManning) then
                    !New Over Land Coef
                    Coef = OverLandCoefficient * (1.0 + 100.0 * (WaterColumn - Me%MinimumWaterColumn))
                else
                    Coef = OverLandCoefficient 
                endif
                
                !m3.s-1 = m * m(5/3) / (s.m(-1/3)) = m * m2 / s
                FlowRouting = Width * (WaterColumn - Me%MinimumWaterColumn)**(5./3.) * sqrt(Slope) / Coef

!            !Chezy units are wrong so it was removed
!            else if (RoutingMethod == Chezy_) then
!                FlowRouting = Width * (WaterColumn - Me%MinimumWaterColumn) * sqrt(Slope) * OverLandCoefficient
            end if
            
            if (Me%ImposeMaxVelocity) then
            
                Area   = (Width  * (WaterColumn - Me%MinimumWaterColumn))
                VelAux = FlowRouting / Area
                if (VelAux > Me%ImposedMaxVelocity) then
                    FlowRouting = Me%ImposedMaxVelocity * Area
                endif
            
            endif


        else
            
            FlowRouting = 0.0
        
        endif
        
        
        

    end function FlowRouting

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

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
                !m3 = m3 + (m * m2)
                Me%TotalStoredVolume = Me%TotalStoredVolume + (Me%MyWaterColumn(i,j)             * &
                                       Me%ExtVar%GridCellArea(i,j))
                                       
            endif

        enddo
        enddo


        if (Me%StormWaterDrainage) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                if (Me%ExtVar%BasinPoints(i, j) == 1) then
                    
                    !m3 = m3 + m3
                    Me%TotalStoredVolume = Me%TotalStoredVolume + Me%StormWaterVolume(i, j)
                                           
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

        write(FinalFile)Me%MyWaterLevel
        
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


        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunOffID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then


            nUsers = DeassociateInstance(mRUNOFF_,  Me%InstanceID)

            if (nUsers == 0) then

                !Writes file with final condition
                call WriteFinalFile

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


                nullify    (Me%iFlowX)
                nullify    (Me%iFlowY)
                nullify    (Me%lFlowX)
                nullify    (Me%lFlowY)
                nullify    (Me%iFlowToChannels)
                nullify    (Me%lFlowToChannels)
                nullify    (Me%lFlowBoundary)
                nullify    (Me%iFlowBoundary)


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
        
            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPondedWaterColumn - ModuleRunOff - ERR01'
        
            !Gets a pointer to Topography
            call GetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPondedWaterColumn - ModuleRunOff - ERR02'
        
            idx = 1
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    if (Me%StormWaterDrainageCoef(i, j) > AllmostZero) then
                        waterColumn(idx) = Me%MyWaterColumn(i, j)

                        !Removes Water from this module                        
                        Me%myWaterColumnOld(i, j) = 0.0
                        Me%myWaterColumn(i, j)    = 0.0
                        Me%myWaterLevel (i, j) = Me%ExtVar%Topography(i, j)
                        Me%myWaterVolume(i, j) = 0.0
                    else
                        waterColumn(idx) = 0.0
                    endif
                    idx = idx + 1 
                endif
            enddo
            enddo

            !Ungets the Topography
            call UngetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPondedWaterColumn - ModuleRunOff - ERR10'

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPondedWaterColumn - ModuleRunOff - ERR11'

            GetPondedWaterColumn = .true.
        else 
            call PlaceErrorMessageOnStack("Runoff not ready")
            GetPondedWaterColumn = .false.
        end if
           

    end function GetPondedWaterColumn
    
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
        
            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPondedWaterColumn - ModuleRunOff - ERR01'
        
            idx = 1
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    Me%StormWaterModelFlow(i, j) = overlandToSewerFlow(idx)
                    idx = idx + 1 
                endif
            enddo
            enddo

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPondedWaterColumn - ModuleRunOff - ERR02'

            SetStormWaterModelFlow = .true.
        else 
            call PlaceErrorMessageOnStack("Runoff not ready")
            SetStormWaterModelFlow = .false.
        end if
           

    end function SetStormWaterModelFlow
        

#endif

end module ModuleRunOff
