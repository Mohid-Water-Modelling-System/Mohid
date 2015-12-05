!--------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!--------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land
! MODULE        : Irrigation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jul 2015
! REVISION      : Eduardo Jauch - v4.0
! DESCRIPTION   : Module which deals with irrigation
!
!--------------------------------------------------------------------------------
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
!--------------------------------------------------------------------------------
   
module ModuleIrrigation

    use ModuleGlobalData
    use ModuleTime
    use ModuleTimeSerie         ,only : StartTimeSerie, StartTimeSerieInput,             &
                                        KillTimeSerie, GetNumberOfTimeSeries,            &
                                        GetTimeSerieInitialData, GetTimeSerieValue,      &
                                        GetTimeSerieLocation, GetTimeSerieName,          &
                                        TryIgnoreTimeSerie, CorrectsCellsTimeSerie,      &
                                        WriteTimeSerieLineNow
    
    use ModuleEnterData
    use ModuleHDF5
    use ModuleFunctions         ,only : TimeToString, SetMatrixValue, ChangeSuffix,      &
                                        CHUNK_J, LinearInterpolation,                    &
                                        InterpolateValueInTime, ConstructPropertyID
    use ModuleHorizontalGrid    ,only : GetHorizontalGridSize, GetHorizontalGrid,        &
                                        UnGetHorizontalGrid, WriteHorizontalGrid,        &
                                        GetGridCellArea, GetXYCellZ,                     &
                                        GetCellZInterceptByLine,                         &
                                        GetCellZInterceptByPolygon
    use ModuleHorizontalMap     ,only : GetBoundaries, UngetHorizontalMap
    use ModuleGridData          ,only : GetGridData, UngetGridData, WriteGridData
    use ModuleGeometry
    use ModuleBasinGeometry     ,only : GetBasinPoints, GetRiverPoints, GetCellSlope,    &
                                        GetDrainageDirection, TargetPoint,               &
                                        UnGetBasin
    use ModuleStopWatch         ,only : StartWatch, StopWatch
    use ModuleFillMatrix        ,only : ConstructFillMatrix, ModifyFillMatrix,           &
                                        KillFillMatrix
   
    implicit none

    private 
     
    !Subroutines-----------------------------------------------------------------

    !Constructor
    public  ::  ConstructIrrigation
    private ::      AllocateInstance
    private ::      ReadDataFile
    private ::      AllocateVariables
    private ::      InitializeVariables
     
    !Selector
    public  ::  GetIrrigationDTForNextEvent
    public  ::  GetIrrigationPMIsRequired
    public  ::  GetIrrigationThresholds
    public  ::  GetIrrigationRequirements
    public  ::  GetIrrigationFlux
    public  ::  SetIrrigationThresholds
    public  ::  UnGetIrrigation
   
    !Modifier
    public  ::  ModifyIrrigation
   
    !Destructor
    public  ::  KillIrrigation     
   
    !Management
    private ::  ReadLockExternalVar
    private ::  ReadUnLockExternalVar
    private ::  Ready
    private ::      LocateObjIrrigation    
   
    !Interfaces------------------------------------------------------------------
    private :: UnGetIrrigation2D_R4
    private :: UnGetIrrigation2D_R8
    interface  UnGetIrrigation
        module procedure UnGetIrrigation2D_R4
        module procedure UnGetIrrigation2D_R8
    end interface  UnGetIrrigation
   
    !Parameters------------------------------------------------------------------
    integer, parameter                              :: ApplicationArea_                 = 1
    integer, parameter                              :: StartHeadThreshold_              = 2
    integer, parameter                              :: StartHeadLimitThreshold_         = 3
    integer, parameter                              :: EndHeadThreshold_                = 4
    integer, parameter                              :: IrrigationProperty_              = 5
    integer, parameter                              :: MinimumIntervalBetweenEvents_    = 6
    integer, parameter                              :: StartInstantThreshold_           = 7
    integer, parameter                              :: GearType_                        = 8
    integer, parameter                              :: GearEfficiency_                  = 9
    integer, parameter                              :: Debit_                           = 10
    integer, parameter                              :: Position_                        = 11
    integer, parameter                              :: AccIrrigation_                   = 12
    integer, parameter                              :: SecondsToEventEnd_               = 13
    
    character(20), parameter                        :: BeginSchedule    = "<BeginSchedule>"
    character(20), parameter                        :: EndSchedule      = "<EndSchedule>"

    character(20), parameter                        :: BeginIrriProp    = "<BeginProperty>"
    character(20), parameter                        :: EndIrriProp      = "<EndProperty>"
    
    integer, parameter                              :: FixedIrrigation      = 1
    integer, parameter                              :: IrrigationBySteps    = 2
    integer, parameter                              :: ContinuousIrrigation = 3
    
    !Gear Type
    integer, parameter                              :: CenterPivot      = 1
    integer, parameter                              :: LinearPivot      = 2
    integer, parameter                              :: Sprinkler        = 3
    integer, parameter                              :: GravitySystem    = 4
    integer, parameter                              :: DripIrrigation   = 5
    
    !Types-----------------------------------------------------------------------
    
    !============================================================================
    !T_OutPut
    !
    !Contains the information required to set the output (HDF), the FIN file and
    !RESTART files
    !============================================================================
    type T_OutPut
        
        type (T_Time), pointer, dimension(:)        :: OutTime => null()
        
        integer                                     :: NextOutPut = 1
        logical                                     :: Yes = .false.
        
        type (T_Time), dimension(:), pointer        :: RestartOutTime => null()
        logical                                     :: WriteRestartFile = .false.
        logical                                     :: RestartOverwrite = .false.
        integer                                     :: NextRestartOutput = 1 
        
        logical                                     :: BoxFluxes = .false.        
        logical                                     :: TimeSerie_On = .false.        
        logical                                     :: HDF_On = .false.
        
        integer                                     :: Number = 0
        
    end type T_OutPut
    
    type T_AccOutput
        
        logical                                     :: Yes = .false.
        
        type (T_Time)                               :: NextOutputTime
        real                                        :: OutputInterval = 86400.0
        
        real                                        :: Value = 0.0
        
    end type T_AccOutput
    
    !============================================================================
    !T_IrriProperty
    !
    !Contains the information required to deal with Properties (user defined) and
    !INTERNAL Properties (module defined)
    !============================================================================
    type T_IrriProperty
        
        type(T_PropertyID)                          :: ID
               
        real                                        :: DefaultValue = 0.0
        real, pointer, dimension (:,:)              :: Field => null()
        real, pointer, dimension (:,:,:)            :: Field3D => null()
        logical, pointer, dimension (:,:)           :: LogicalField => null()
        logical, pointer, dimension (:,:,:)         :: LogicalField3D => null()
        
        logical                                     :: IsLogical = .false.
        
        type(T_IrriProperty), pointer               :: Next => null()
        type(T_IrriProperty), pointer               :: Prev => null()
        
        logical                                     :: Old = .false.
        
        logical                                     :: OutputToHDF = .false.
        logical                                     :: OutputToTimeseries = .false.
        logical                                     :: SaveToFinalFile = .false.
        
    end type T_IrriProperty
    
    !============================================================================
    !T_DailySchedule
    !
    !Contains the information required to set the daily irrigation schedules.
    !============================================================================
    type T_DailySchedule
        
        type(T_Time)                                :: StartInstant
        type(T_Time)                                :: EndInstant

        type(T_Time)                                :: StartDay
        type(T_Time)                                :: EndDay
        
        type(T_STime)                               :: SStartInstant
        type(T_STime)                               :: SEndInstant
        
        real                                        :: Irrigation = 0.0 !in 'mm/s'
        real                                        :: TotalIrrigation = 0.0 !in 'mm'
        
        type(T_DailySchedule), pointer              :: Next => null()
        type(T_DailySchedule), pointer              :: Prior => null()
        
    end type T_DailySchedule
    
    !============================================================================
    !T_IrriSchedule
    !
    !Contains the information required to deal with every SCHEDULE provided by 
    !the user in the MODULE INPUT FILE.
    !============================================================================
    type T_IrriSchedule
        
        type(T_PropertyID)                          :: ID
        
        integer                                     :: ObjTimeSeries = 0
        
        integer                                     :: IrrigationMethod = FixedIrrigation
        
        logical                                     :: IntegratedSchedule = .true.
        logical                                     :: AutoIrrigation = .false.
        
        logical                                     :: IsIrrigating = .false.
        type(T_DailySchedule), pointer              :: FirstDailySchedule => null()
        type(T_DailySchedule), pointer              :: LastDailySchedule => null()
        integer                                     :: DailySchedulesNumber = 0
        
        
        logical                                     :: SingleSystem = .false. !If set to .true. will integrate the area and give a single value
                                                                              !For detailed simulations where the system area is larger than the 
                                                                              !cell area
        
        integer, pointer, dimension(:,:)            :: RootsKLB
        
        !ApplicationAreaMap(logical) - Application area map (mask)
        !type(T_IrriProperty), pointer               :: ApplicationAreaMap => null() 
        type(T_IrriProperty)               :: ApplicationAreaMap 
        
        !ApplicationArea(m2) - Application area
        real                                        :: ApplicationArea = 0.0 !In 'm2'
        
        !MinimumAreaToStartIrrigation(m2) - Minimum area to start irrigation
        real                                        :: MinimumAreaToStartIrrigation = 0.0 !In 'm2'
        
        !Irrigation(mm) - Property Irrigation        
        type(T_IrriProperty), pointer               :: Irrigation => null() !in 'mm'
        real                                        :: ActualIrrigation = 0.0
                
        !HeadWiltingPoint(m) - Head value to start irrigating
        real                                        :: HeadWiltingPoint = 0 
                
        !HeadTarget(m) - Head value to be achieved by irrigation
        real                                        :: HeadTarget = 0
        
        real                                        :: TimeSinceLastEvent = 9.9E+15 !in 's'
        real                                        :: MinimumIntervalBetweenEvents = 0.0 !In 's'
        integer                                     :: MaxConsecutiveDays = 2
        real                                        :: MaxDailyIrrigationTime = 86400.0
        
        real                                        :: StartInstantThreshold = 0.0 !in 'hours'
        real                                        :: EndInstantThreshold = 0.0 !in 'hours'
               
        real                                        :: MaxDepthToCheckSaturation = 0.2 !In 'm'
        real                                        :: SaturationThreshold = 0.9 !In 'decimal %' (0 to 1)
        real                                        :: MaxSaturatedFraction = 0.3 !In 'decimal %' (0 to 1)
    
        integer                                     :: GearType = 1 !1 to 5        
        real                                        :: GearEfficiency = 1 !in '%': 0.01 - 1.00        
        real                                        :: GearDebit = 10.0 !in 'mm/hour'
        real                                        :: GearMinimumVelocity = 5.0 !in 'km/hour'
        real                                        :: GearMaximumVelocity = 5.0 !in 'km/hour'
        integer                                     :: GearPosition = 0 !0 for Above Plant, 1 for Below Plant
        
        real                                        :: MinimumToIrrigate = 5.0 !in 'mm'
        
        type(T_AccOutput)                           :: AccOutput
        
        !type(T_IrriProperty), pointer               :: AccIrrigation => null() !in 'mm'
        !type(T_IrriProperty), pointer               :: ToIrrigate => null() !In 'mm'
        
        type(T_IrriProperty), pointer               :: WaterContentTarget => null() !In Theta % (0 - 100). Must be higher than the Theta Residual and lower than the Theta Saturation
        type(T_IrriProperty), pointer               :: WaterContentEasy => null() !
        type(T_IrriProperty), pointer               :: WaterContentWiltingPoint => null() !In Theta % (0 - 100). Must be higher than the Theta Residual and lower than the Theta Saturation
        
        integer                                     :: PropertiesNumber = 0
        type(T_IrriProperty), pointer               :: FirstProperty => null()
        type(T_IrriProperty), pointer               :: LastProperty => null()
        
        logical                                     :: IsIrrigationScheduledForToday = .false.
        
    end type T_IrriSchedule
    
    !============================================================================
    !T_Files
    !
    !Stores the names of the files used by the module
    !============================================================================
    type T_Files
       
        character(PathLength)                       :: DataFile = null_str
        character(PathLength)                       :: OutputHDFFile = null_str
        character(PathLength)                       :: InitialFile = null_str
        character(PathLength)                       :: FinalFile = null_str
      
    end type T_Files     
     
    !============================================================================
    !T_IrrigationData
    !
    !Stores the pointers to the arrays that are passed to the MODIFY subroutine 
    !from the BASIN module, with information from other modules used to compute
    !the irrigation needs
    !============================================================================
    public :: T_IrrigationData
    type T_IrrigationData
       
        integer, dimension(:,:), pointer            :: BasinPoints => null()
        real, dimension(:,:), pointer               :: Topography => null()
        real, dimension(:,:), pointer               :: RootsDepth => null()
        real, dimension(:,:), pointer               :: LAISenescence => null()
        real, dimension(:,:,:), pointer             :: SoilWaterContent => null()
        real, dimension(:,:,:), pointer             :: SoilRelativeWaterContent => null()
        real, dimension(:,:,:), pointer             :: SoilHead => null()
        real, dimension(:,:,:), pointer             :: DWZ => null()
        real, dimension(:,:), pointer               :: Areas => null()
      
    end type T_IrrigationData   
   
    !============================================================================
    !T_Irrigation
    !
    !Main 'object' of the module
    !============================================================================
    type  T_Irrigation
       
        integer                                     :: InstanceID = 0
        character(len=StringLength)                 :: ModelName = null_str
      
        integer                                     :: ObjBasinGeometry = 0
        integer                                     :: ObjTime = 0
        integer                                     :: ObjHorizontalGrid = 0
        integer                                     :: ObjHorizontalMap = 0
        integer                                     :: ObjGridData = 0
        integer                                     :: ObjHDF5 = 0
        integer                                     :: ObjIniHDF5 = 0
        integer                                     :: ObjEnterData = 0  
        integer                                     :: ObjTimeSerie = 0
      
        type (T_IrrigationData), pointer            :: Data
        type (T_Files)                              :: Files
        type (T_Output)                             :: Output
      
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime      
      
        !Grid size
        type (T_Size3D)                             :: Size
        type (T_Size3D)                             :: WorkSize
      
        logical                                     :: Continuous = .false.
        logical                                     :: StopOnWrongDate = .true.      
     
        real, dimension(:,:), pointer               :: Irrigation => null()
        real, dimension(:,:), pointer               :: IrrigationFlux => null()
        
        integer                                     :: NumberOfSchedules = 0
        type(T_IrriSchedule), dimension(:), pointer :: Schedules => null()
        
        logical, dimension(:,:), pointer            :: ComputePoints => null()
        
        type (T_Time)                               :: Now
        type (T_Time)                               :: Day
        type (T_STime)                              :: SNow
        type (T_STime)                              :: SDay
        real                                        :: DT = null_real
        real                                        :: DTForNextEvent = null_real
        
        real, pointer, dimension(:)                 :: AccIrrigation => null()
        
        character(len=PathLength)                   :: TimeSeriesLocation = "*****"

        type(T_Irrigation), pointer                 :: Next => null()
        
        integer                                     :: GlobalCounter = 0

    end type T_Irrigation
   
    !Global Module Variables
    type (T_Irrigation), pointer                    :: FirstObjIrrigation => null()
    type (T_Irrigation), pointer                    :: Me => null()
 
    !----------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    subroutine ConstructIrrigation (ModelName,          &            
                                    id,                 &
                                    ComputeTimeID,      &
                                    HorizontalGridID,   &
                                    HorizontalMapID,    &
                                    GridDataID,         &
                                    BasinGeometryID,    &
                                    GeometryID,         &                                    
                                    STAT)

        !Arguments---------------------------------------------------------------
        character(len=*)                                :: ModelName
        integer                                         :: id
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: GridDataID
        integer                                         :: BasinGeometryID
        integer                                         :: GeometryID
        integer, optional, intent(OUT)                  :: STAT           

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, stat_call

        !------------------------------------------------------------------------
        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.NOT. ModuleIsRegistered(mIrrigation_)) then
            
            nullify (FirstObjIrrigation)
            call RegisterModule (mIrrigation_) 
            
        endif

        ready_ = Ready (id)
        
cd0:    if (ready_ == OFF_ERR_) then

            call AllocateInstance
            
            Me%ModelName = ModelName
            
            !Associates External Instances
            Me%ObjTime            = AssociateInstance (mTIME_           , ComputeTimeID     )
            Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_ , HorizontalGridID  )
            Me%ObjHorizontalMap   = AssociateInstance (mHORIZONTALMAP_  , HorizontalMapID   )
            Me%ObjGridData        = AssociateInstance (mGRIDDATA_       , GridDataID        )
            Me%ObjBasinGeometry   = AssociateInstance (mBASINGEOMETRY_  , BasinGeometryID   )

            call ReadLockExternalVar ()

            !Gets Compute Time Limits
            call GetComputeTimeLimits (Me%ObjTime, BeginTime = Me%BeginTime, &
                                       EndTime = Me%EndTime, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ConstructIrrigation - ModuleIrrigation - ERR020'            
            
            !Geometry Size
            call GetGeometrySize(GeometryID,                 &    
                                 Size     = Me%Size,         &
                                 WorkSize = Me%WorkSize,     &
                                 STAT = stat_call)
            if (stat_call /= SUCCESS_) stop 'ConstructIrrigation - ModuleIrrigation - ERR030'
                      
            call AllocateVariables
            call ReadDataFile                     
            
            if (Me%Continuous) call OpenInitialFile                     
            call ConstructSchedules
            call InitializeVariables
            if (Me%Continuous) call CloseInitialFile
            
            call CheckConfiguration
            
            call ConstructTimeSeriesOutput ()
         
            call ReadUnLockExternalVar ()

            !Returns ID
            id = Me%InstanceID
      
            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructIrrigation - ModuleIrrigation - ERR040' 

        end if cd0

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine ConstructIrrigation
 
    !----------------------------------------------------------------------------
                              
    subroutine AllocateInstance

        !Arguments---------------------------------------------------------------
                                                    
        !Local-------------------------------------------------------------------
        type (T_Irrigation), pointer :: obj_new
        type (T_Irrigation), pointer :: obj_prev

        !------------------------------------------------------------------------
        !Allocates new instance
        allocate (obj_new)
        nullify  (obj_new%Next)

        !Insert New Instance into list and makes Current point to it
        if (.NOT. associated(FirstObjIrrigation)) then
            
            FirstObjIrrigation => obj_new
            Me => obj_new
            
        else
            
            obj_prev => FirstObjIrrigation
            Me => FirstObjIrrigation%Next
            
            do while (associated(Me))
                
                obj_prev => Me
                Me => Me%Next
                
            enddo
            
            Me => obj_new
            obj_prev%Next => obj_new
            
        endif

        Me%InstanceID = RegisterNewInstance (mIrrigation_)

    end subroutine AllocateInstance

    !----------------------------------------------------------------------------

    subroutine ReadDataFile

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: stat_call
        integer                                     :: iflag

        !Reads the name of the data file from nomfich
        call ReadFileName ('IRRIGATION_DATA', Me%Files%DataFile, "Irrigation Data File", STAT = stat_call)
        if (stat_call /= SUCCESS_) stop 'ReadDataFile - ModuleIrrigation - ERR010'

        call ReadFileName ('IRRIGATION_FIN', Me%Files%FinalFile, Message = "Irrigation Final File", STAT = stat_call)
        if (stat_call /= SUCCESS_) stop 'ReadDataFile - ModuleIrrigation - ERR020'                
        
        !Constructs the DataFile
        call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = stat_call)
        if (stat_call /= SUCCESS_) stop 'ReadDataFile - ModuleIrrigation - ERR030'

        !Continuous Computation
        call GetData (Me%Continuous,                        &
                      Me%ObjEnterData, iflag,               &
                      SearchType   = FromFile,              &
                      keyword      = 'CONTINUOUS',          &
                      default      = .false.,               &
                      ClientModule = 'ModuleIrrigation',    &
                      STAT         = stat_call)
        if (stat_call /= SUCCESS_) stop 'ReadDataFile - ModuleIrrigation - ERR040'
        
        call GetData(Me%TimeSeriesLocation,                 &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'TIME_SERIE_LOCATION',  &
                     ClientModule = 'ModuleIrrigation',     &
                     Default      = Me%Files%DataFile,      &
                     STAT         = stat_call)
        if (stat_call .NE. SUCCESS_) stop 'ReadDataFile - ModuleIrrigation - ERR050'
        
      
    end subroutine ReadDataFile
   
    !----------------------------------------------------------------------------

    subroutine AllocateVariables

        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        !Begin-------------------------------------------------------------------
        allocate (Me%Irrigation (Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
        Me%Irrigation = 0.0
        
        allocate (Me%IrrigationFlux (Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
        Me%IrrigationFlux = 0.0
    
        allocate (Me%AccIrrigation (1))
        Me%AccIrrigation = 0.0
        
        allocate (Me%ComputePoints (Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
        Me%ComputePoints = 0
                  
    end subroutine AllocateVariables
       
    !----------------------------------------------------------------------------
   
    subroutine InitializeVariables

        !Arguments---------------------------------------------------------------
        !Local-------------------------------------------------------------------
        integer                             :: i
        type(T_IrriSchedule), pointer       :: schedule       
        integer                             :: ILB,IUB
        integer                             :: JLB,JUB
        integer                             :: KLB,KUB
      
        !Begin-------------------------------------------------------------------
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
        
        do i = 1, Me%NumberOfSchedules
            allocate (Me%Schedules(i)%WaterContentTarget)
            allocate (Me%Schedules(i)%WaterContentTarget%Field3D (ILB:IUB,JLB:JUB,KLB:KUB))
            Me%Schedules(i)%WaterContentTarget%Field3D = 0.0
            
            allocate (Me%Schedules(i)%WaterContentWiltingPoint)
            allocate (Me%Schedules(i)%WaterContentWiltingPoint%Field3D (ILB:IUB,JLB:JUB,KLB:KUB))
            Me%Schedules(i)%WaterContentWiltingPoint%Field3D = 0.0
            
            allocate (Me%Schedules(i)%WaterContentEasy)
            allocate (Me%Schedules(i)%WaterContentEasy%Field3D (ILB:IUB,JLB:JUB,KLB:KUB))
            Me%Schedules(i)%WaterContentEasy%Field3D = 0.0
            
            allocate (Me%Schedules(i)%RootsKLB (ILB:IUB,JLB:JUB))
            Me%Schedules(i)%RootsKLB = 0.0
        enddo
        
    end subroutine InitializeVariables
       
    !----------------------------------------------------------------------------
    
    subroutine CheckConfiguration ()
   
        !Arguments---------------------------------------------------------------
   
        !Local-------------------------------------------------------------------

        !Begin-------------------------------------------------------------------
      
      
    end subroutine CheckConfiguration
    
    !----------------------------------------------------------------------------  
      
    subroutine ConstructSchedules
      
        !Arguments---------------------------------------------------------------
   
        !Local-------------------------------------------------------------------
        integer                                 :: client_number
        integer                                 :: stat_call
        logical                                 :: block_found
        type (T_IrriSchedule), pointer          :: new_schedule
        integer                                 :: n

        !Begin-------------------------------------------------------------------
        call GetNumberOfBlocks (Me%ObjEnterData, BeginSchedule, EndSchedule,    &
                                FromFile_, Me%NumberOfSchedules,                &
                                client_number, stat_call)
        if (stat_call /= SUCCESS_) stop 'ConstructSchedules - ModuleIrrigation - ERR010'
        if (Me%NumberOfSchedules <= 0) stop 'ConstructSchedules - ModuleIrrigation - ERR020'
        
        allocate (Me%Schedules (Me%NumberOfSchedules))
        
        do n = 1, Me%NumberOfSchedules
            call ExtractBlockFromBuffer (Me%ObjEnterData,                   &
                                         ClientNumber    = client_number,   &
                                         block_begin     = BeginSchedule,   &
                                         block_end       = EndSchedule,     &
                                         BlockFound      = block_found,     &
                                         STAT            = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ConstructSchedules - ModuleIrrigation - ERR030'
            
            if (.NOT. block_found) &
                stop 'ConstructSchedules - ModuleIrrigation - ERR040'

            !Construct the new schedule
            new_schedule => Me%Schedules(n)
            Call ConstructSchedule (new_schedule, client_number)

        enddo
                
        if (Me%NumberOfSchedules > 0) then
            call RewindBuffer(Me%ObjEnterData, STAT = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ConstructSchedules - ModuleIrrigation - ERR050'
        endif
        
    end subroutine ConstructSchedules
   
    !----------------------------------------------------------------------------
    
    subroutine ConstructSchedule (new_schedule, client_number)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: new_schedule
        integer                             :: client_number

        !External----------------------------------------------------------------
        integer                             :: stat_call
        integer                             :: iflag

        !Begin-------------------------------------------------------------------
        call GetData (new_schedule%ID%Name,                 &
                      Me%ObjEnterData, iflag,               &
                      Keyword      = 'NAME',                &
                      ClientModule = 'ModuleIrrigation',    &
                      SearchType   = FromBlock,             &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR010'
        
        call GetData (new_schedule%IrrigationMethod,        &
                      Me%ObjEnterData, iflag,               &
                      Keyword      = 'METHOD',              &
                      ClientModule = 'ModuleIrrigation',    &
                      Default      = FixedIrrigation,       &
                      SearchType   = FromBlock,             &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR020'
        
        call GetData (new_schedule%SingleSystem,            &
                      Me%ObjEnterData, iflag,               &
                      Keyword      = 'SINGLE_SYSTEM',       &
                      ClientModule = 'ModuleIrrigation',    &
                      Default      = .false.,               &
                      SearchType   = FromBlock,             &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR030'
        
        call GetData (new_schedule%HeadWiltingPoint,        &
                      Me%ObjEnterData, iflag,               &
                      Keyword      = 'HEAD_WILTING_POINT',  &
                      ClientModule = 'ModuleIrrigation',    &
                      Default      = -150.0,                &
                      SearchType   = FromBlock,             &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR040'
        
        call GetData (new_schedule%HeadTarget,              &
                      Me%ObjEnterData, iflag,               &
                      Keyword      = 'HEAD_TARGET',         &
                      ClientModule = 'ModuleIrrigation',    &
                      Default      = -10.0,                 &
                      SearchType   = FromBlock,             &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR050'
        
        call GetData (new_schedule%MinimumIntervalBetweenEvents,    &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'MIN_INTERVAL_BETWEEN_EVENTS', &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 86400.0,                       &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR060'
        
        new_schedule%TimeSinceLastEvent = new_schedule%MinimumIntervalBetweenEvents + 1
        
        call GetData (new_schedule%MaxConsecutiveDays,              &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'MAX_CONSECUTIVE_DAYS',        &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 2,                             &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR070'
        
        call GetData (new_schedule%MaxDailyIrrigationTime,          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'MAX_DAILY_IRRIGATION_TIME',   &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 86400.0,                       &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR080'
        
        call GetData (new_schedule%StartInstantThreshold,           &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'START_INSTANT_THRESHOLD',     &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 0.0,                           &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR090'    
        
        call GetData (new_schedule%EndInstantThreshold,             &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'END_INSTANT_THRESHOLD',       &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 23.0,                          &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR100'   
        
        call GetData (new_schedule%MaxDepthToCheckSaturation,       &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'MAX_DEPTH_TO_CHECK',          &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 0.2,                           &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR110'
        
        call GetData (new_schedule%SaturationThreshold,             &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'SATURATION_THRESHOLD',        &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 0.9,                           &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR120'
        
        call GetData (new_schedule%MaxSaturatedFraction,            &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'MAX_SATURATED_FRACTION',      &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 0.3,                           &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR130'
        
        call GetData (new_schedule%GearType,                        &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'GEAR_TYPE',                   &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 1,                             &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR140'

        call GetData (new_schedule%GearEfficiency,                  &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'GEAR_EFFICIENCY',             &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 0.85,                          &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR150'
        
        call GetData (new_schedule%GearDebit,                       &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'GEAR_DEBIT',                  &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 10.0,                          &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR160'
        
        call GetData (new_schedule%GearMinimumVelocity,             &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'GEAR_MIN_VEL',                &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 0.2,                           &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR170'
        
        call GetData (new_schedule%GearMaximumVelocity,             &
                      Me%ObjEnterData, iflag,                       &
                      Keyword      = 'GEAR_MAX_VEL',                &
                      ClientModule = 'ModuleIrrigation',            &
                      Default      = 4.0,                           &
                      SearchType   = FromBlock,                     &
                      STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructSchedule - ModuleIrrigation - ERR180'
        
        call ConstructProperties (new_schedule, client_number)
        call ConstructInternalProperties (new_schedule)
        
        if (Me%Continuous) then
            call ConstructDailySchedules (new_schedule, Me%ObjIniHDF5)
        endif
    
    end subroutine ConstructSchedule
    
    !----------------------------------------------------------------------------

    subroutine ConstructProperties (new_schedule, client_number)
   
        !Arguments---------------------------------------------------------------
        type (T_IrriSchedule), pointer      :: new_schedule
        integer                             :: client_number
   
        !Local-------------------------------------------------------------------
        integer                             :: stat_call
        logical                             :: block_found
        type (T_IrriProperty), pointer      :: new_property

        !Begin-------------------------------------------------------------------
      
do1 :   do      
            call ExtractBlockFromBlock (Me%ObjEnterData,                            &
                                        ClientNumber        = client_number,        &
                                        block_begin         = BeginIrriProp,        &
                                        block_end           = EndIrriProp,          &
                                        BlockInBlockFound   = block_found,          &
                                        STAT                = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ConstructProperties - ModuleIrrigation - ERR010'
            
cd1 :       if (block_found) then

                !Construct a New Property 
                Call ConstructProperty (new_schedule, new_property)

                !Add new Property to the SoilProperties List 
                Call AddProperty (new_schedule, new_property)   
   
            else cd1

                call Block_Unlock(Me%ObjEnterData, client_number, STAT = stat_call) 
                if (stat_call /= SUCCESS_) &
                    stop 'ConstructProperties - ModuleIrrigation - ERR020'
               
                exit do1    !No more blocks            
                
            end if cd1

        enddo do1
      
    end subroutine ConstructProperties
   
    !----------------------------------------------------------------------------
    
    subroutine ConstructInternalProperties (new_schedule)
   
        !Arguments---------------------------------------------------------------
        type (T_IrriSchedule), pointer      :: new_schedule
   
        !Local-------------------------------------------------------------------
        type (T_IrriProperty), pointer      :: new_property
        integer                             :: stat_call
      
        !Begin-------------------------------------------------------------------
      
        !Accumulated Irrigation
        allocate (new_property, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'AllocateInternalProperties - ModuleIrrigation - ERR010'        
        call ConstructInternalProperty (new_schedule, new_property, 'acc. irrigation', Me%Continuous, .true.)

        !allocate (new_property, STAT = stat_call)
        !if (stat_call /= SUCCESS_) &
        !    stop 'AllocateInternalProperties - ModuleIrrigation - ERR010'        
        !call ConstructInternalProperty (new_schedule, new_property, 'acc. irrigation', Me%Continuous, .true.)

    end subroutine ConstructInternalProperties
   
    !----------------------------------------------------------------------------
    
    subroutine ConstructInternalProperty (new_schedule, new_property, name, old, save_to_final_file)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: new_schedule
        type(T_IrriProperty), pointer       :: new_property
        character(LEN=*)                    :: name
        logical                             :: old
        logical                             :: save_to_final_file

        !------------------------------------------------------------------------
        new_property%ID%Name = name
        new_property%Old = old
        new_property%SaveToFinalFile = save_to_final_file
        
        call CheckIrriPropertyName (new_schedule, new_property)
        call ConstructPropertyValues (new_schedule, new_property, .true.)
        call AddProperty (new_schedule, new_property)
        
    end subroutine ConstructInternalProperty
    
    !----------------------------------------------------------------------------
     
    subroutine ConstructProperty (new_schedule, new_property)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: new_schedule
        type(T_IrriProperty), pointer       :: new_property

        !External----------------------------------------------------------------
        integer                             :: stat_call

        !------------------------------------------------------------------------
             
        allocate (new_property, STAT = stat_call)
        if(stat_call /= SUCCESS_) then
            print *, '[ModuleIrrigation] Was not possible to allocate memory to new property'
            print *, '[ModuleIrrigation] during module initialization.'
            stop 'ConstructProperty - ModuleIrrigation - ERR010'
        endif
              
        call ConstructPropertyID        (new_property%ID, Me%ObjEnterData, FromBlockInBlock, .false.)
        call CheckIrriPropertyName      (new_schedule, new_property)
        call ConstructPropertyValues    (new_schedule, new_property, .false.)        

    end subroutine ConstructProperty
    
    !----------------------------------------------------------------------------
    
    subroutine CheckIrriPropertyName (new_schedule, new_property)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: new_schedule
        type(T_IrriProperty), pointer       :: new_property

        !External----------------------------------------------------------------

        !------------------------------------------------------------------------
        if ((new_property%ID%Name) == "application area") then
            new_property%ID%IDNumber = ApplicationArea_
            !new_schedule%ApplicationAreaMap => new_property
            allocate (new_schedule%ApplicationAreaMap%LogicalField(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))            
            new_property%IsLogical = .true.
            
        !elseif (trim(new_property%ID%Name) == "start head threshold") then
        !    new_property%ID%IDNumber = StartHeadThreshold_
        !    new_schedule%HeadWiltingPoint => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "start head limit threshold") then
        !    new_property%ID%IDNumber = StartHeadLimitThreshold_
        !    new_schedule%HeadTarget => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "end head threshold") then
        !    new_property%ID%IDNumber = EndHeadThreshold_
        !    new_schedule%EndHeadTreshold => new_property
        !    
        elseif (trim(new_property%ID%Name) == "fixed irrigation") then
            
            new_property%ID%IDNumber = IrrigationProperty_
            new_schedule%Irrigation => new_property
            
        !elseif (trim(new_property%ID%Name) == "minimum interval between events") then
        !    new_property%ID%IDNumber = MinimumIntervalBetweenEvents_
        !    new_schedule%MinimumIntervalBetweenEvents => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "start instant threshold") then
        !    new_property%ID%IDNumber = StartInstantThreshold_
        !    new_schedule%StartInstantThreshold => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "gear type") then
        !    new_property%ID%IDNumber = GearType_
        !    new_schedule%GearType => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "gear efficiency") then
        !    new_property%ID%IDNumber = GearEfficiency_
        !    new_schedule%GearEfficiency => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "debit") then
        !    new_property%ID%IDNumber = Debit_
        !    new_schedule%Debit => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "position") then
        !    new_property%ID%IDNumber = Position_
        !    new_schedule%Position => new_property
        !    
        !elseif (trim(new_property%ID%Name) == "acc. irrigation") then
        !    new_property%ID%IDNumber = AccIrrigation_
        !    
        elseif (trim(new_property%ID%Name) == "acc. irrigation") then
            new_property%ID%IDNumber = AccIrrigation_
            
        else
            
            print *, '[ModuleIrrigation] Found unknown property in the data file: '
            print *, '[ModuleIrrigation] ', trim(new_property%ID%Name)
            stop 'CheckIrriPropertyName - ModuleIrrigation - ERR010'
            
        endif
        
    end subroutine CheckIrriPropertyName
    
    !----------------------------------------------------------------------------
   
    subroutine ConstructPropertyValues (new_schedule, new_property, is_internal)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: new_schedule
        type(T_IrriProperty), pointer       :: new_property
        logical                             :: is_internal

        !External----------------------------------------------------------------
        integer                             :: stat_call

        !Local-------------------------------------------------------------------
        integer                             :: iflag
        integer                             :: i,j
        integer                             :: ILB,IUB
        integer                             :: JLB,JUB
        integer                             :: WorkSizeILB, WorkSizeIUB
        integer                             :: WorkSizeJLB, WorkSizeJUB
        integer, dimension(:,:), pointer    :: basin_points
      
        !Begin-------------------------------------------------------------------
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        WorkSizeILB = Me%WorkSize%ILB
        WorkSizeIUB = Me%WorkSize%IUB
        WorkSizeJLB = Me%WorkSize%JLB
        WorkSizeJUB = Me%WorkSize%JUB

        !Get water points
        call GetBasinPoints (Me%ObjBasinGeometry, basin_points, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructPropertyValues - ModuleIrrigation - ERR001'
        
        allocate (new_property%Field (ILB:IUB, JLB:JUB), STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructPropertyValues - ModuleIrrigation - ERR010'
      
        new_property%Field(:,:) = 0.0
        
        if (is_internal) then

            new_property%Old = Me%Continuous
        
        else
            
            call GetData (new_property%Old,                         &
                          Me%ObjEnterData, iflag,                   &
                          keyword      = 'OLD',                     &
                          Default      = Me%Continuous,             &                        
                          SearchType   = FromBlockInBlock_,         &
                          ClientModule = 'ModuleIrrigation',        &
                          STAT         = stat_call)              
            if (stat_call /= SUCCESS_) &
                stop 'ConstructPropertyValues - ModuleIrrigation - ERR020'
        
        endif
        
        if ((.NOT. Me%Continuous) .AND. new_property%Old) then
            write (*,*) '[ModuleIrrigation] Property ', trim(new_property%ID%Name), &
                        '[ModuleIrrigation] has OLD set to TRUE, but the CONTINOUS file keyword is missing or set to FALSE'
            stop 'ConstructPropertyValues - ModuleIrrigation - ERR030'
        endif

        if (.NOT. new_property%old .and. .not. is_internal) then
            
            call ConstructFillMatrix (PropertyID           = new_property%ID,           &
                                      EnterDataID          = Me%ObjEnterData,           &
                                      TimeID               = Me%ObjTime,                &
                                      HorizontalGridID     = Me%ObjHorizontalGrid,      &
                                      ExtractType          = FromBlockInBlock,          &
                                      PointsToFill2D       = basin_points,              &
                                      Matrix2D             = new_property%Field,        &
                                      TypeZUV              = TypeZ_,                    &
                                      STAT                 = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ConstructPropertyValues - ModuleIrrigation - ERR040'

            if(.NOT. new_property%ID%SolutionFromFile )then

                call KillFillMatrix (new_property%ID%ObjFillMatrix, STAT = stat_call)
                if (stat_call /= SUCCESS_) &
                    stop 'ConstructPropertyValues - ModuleIrrigation - ERR050'
            
            end if
                       
        elseif (new_property%old) then !It is a property with OLD set to 1 (run continuation)
            
            ! If the property is old then the program is going to try to find a property
            ! with the same name in the Irrigation initial file written in HDF format
            call ReadOldValueFromHDF (new_schedule, new_property)
            
        end if
        
        if (new_property%IsLogical) then
            
            allocate (new_property%LogicalField (ILB:IUB, JLB:JUB), STAT = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ConstructPropertyValues - ModuleIrrigation - ERR010'
        
            !print *, "Get to here..."
            !print *, basin_points
            
            do i = ILB, IUB
            do j = JLB, JUB
                
                if (basin_points(i,j)) then
                    new_schedule%ApplicationAreaMap%LogicalField(i,j) = .true.
                    if (new_property%Field(i,j) > 0.5) then
                        new_property%LogicalField(i,j) = .true.
                    else
                        new_property%LogicalField(i,j) = .false.
                    endif
                else
                    new_schedule%ApplicationAreaMap%LogicalField(i,j) = .false.
                    new_property%LogicalField(i,j) = .false.
                endif
                
            enddo
            enddo
            
        endif
        
        call UnGetBasin(Me%ObjBasinGeometry, basin_points, STAT = stat_call) 
        if (stat_call /= SUCCESS_) &
            stop 'ConstructPropertyValues - ModuleIrrigation - ERR060'

    end subroutine ConstructPropertyValues

    !----------------------------------------------------------------------------
        
    subroutine AddProperty (new_schedule, new_property)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: new_schedule
        type(T_IrriProperty), pointer       :: new_property

        !------------------------------------------------------------------------

        if (.NOT. associated(new_schedule%FirstProperty)) then
            
            new_schedule%PropertiesNumber = 1
            new_schedule%FirstProperty => new_property
            new_schedule%LastProperty => new_property
            
        else
            
            new_property%Prev => new_schedule%LastProperty
            new_schedule%LastProperty%Next => new_property
            new_schedule%LastProperty => new_property
            new_schedule%PropertiesNumber = new_schedule%PropertiesNumber + 1
            
        end if 

    end subroutine AddProperty 

    !----------------------------------------------------------------------------        
    
    subroutine ReadOldValueFromHDF (new_schedule, new_property)

        !Arguments---------------------------------------------------------------
        type (T_IrriSchedule), pointer      :: new_schedule
        type (T_IrriProperty), pointer      :: new_property

        !Local-------------------------------------------------------------------
        integer                             :: stat_call

        !------------------------------------------------------------------------

        if (associated(new_property%field)) then
                        
            call HDF5SetLimits (Me%ObjIniHDF5,                     &
                                Me%WorkSize%ILB, Me%WorkSize%IUB,  &
                                Me%WorkSize%JLB, Me%WorkSize%JUB,  &
                                STAT = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ReadOldValueFromHDF - ModuleIrrigation - ERR010'
            
            call HDF5ReadData (Me%ObjIniHDF5,                                                       &
                               "/Results/"//trim(adjustl(new_schedule%ID%Name))//"/properties",     &
                               trim(adjustl(new_property%ID%Name)),                                 &
                               Array2D = new_property%field,                                        &
                                STAT    = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ReadOldValueFromHDF - ReadOldValueFromHDF - ERR020'
            
        elseif (associated(new_property%field3D)) then
            
            call HDF5SetLimits (Me%ObjIniHDF5,                     &
                                Me%WorkSize%ILB, Me%WorkSize%IUB,  &
                                Me%WorkSize%JLB, Me%WorkSize%JUB,  &
                                Me%WorkSize%KLB, Me%WorkSize%KUB,  &
                                STAT = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ReadOldValueFromHDF - ModuleIrrigation - ERR030'
            
            call HDF5ReadData (Me%ObjIniHDF5,                                                       &
                               "/Results/"//trim(adjustl(new_schedule%ID%Name))//"/properties",     &
                               trim(adjustl(new_property%ID%Name)),                                 &
                               Array3D = new_property%field3D,                                      &
                                STAT    = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ReadOldValueFromHDF - ReadOldValueFromHDF - ERR040'

        else

            stop 'ReadOldValueFromHDF - ReadOldValueFromHDF - ERR050'
            
        endif
        
    end subroutine ReadOldValueFromHDF

    !----------------------------------------------------------------------------
    
    subroutine ConstructDailySchedules (schedule, obj_hdf)

        !Arguments---------------------------------------------------------------
        type (T_IrriSchedule), pointer      :: schedule
        integer                             :: obj_hdf

        !Local-------------------------------------------------------------------
        integer                             :: stat_call
        integer, dimension(2), target       :: schedule_array
        real, dimension(15), target         :: daily_sch_array
        integer, dimension(:), pointer      :: array_integer_ptr
        real, dimension(:), pointer         :: array_real_ptr
        integer                             :: index
        type (T_DailySchedule), pointer     :: new_day
        
        !------------------------------------------------------------------------
        array_integer_ptr => schedule_array
        array_real_ptr => daily_sch_array
        
        call HDF5SetLimits (Me%ObjIniHDF5, 1, 2, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructDailySchedules - ModuleIrrigation - ERR010'
        
        call HDF5ReadData (obj_hdf, "/Schedules/"//trim(adjustl(schedule%ID%Name)), &
                           "info",                                                  &
                           Array1D = array_integer_ptr,                             &
                           STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'ConstructDailySchedules - ModuleIrrigation - ERR020'
        
        if (array_integer_ptr(2) > 0) then
            
            call HDF5SetLimits (Me%ObjIniHDF5, 1, 15, STAT = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'ConstructDailySchedules - ModuleIrrigation - ERR030'
            
            do index = 1, array_integer_ptr(2)                            
                
                call HDF5ReadData (obj_hdf, "/Schedules/"//trim(adjustl(schedule%ID%Name))//"/daily schedules", &
                                   "daily_schedule",                                                            &
                                   Array1D = array_real_ptr,                                                    &
                                   OutputNumber = index,                                                        &
                                   STAT = stat_call)
                if (stat_call /= SUCCESS_) &
                    stop 'ConstructDailySchedules - ModuleIrrigation - ERR040'
                
                call AddDailyScheduleToList (schedule, new_day)
                
                new_day%SStartInstant%Year = array_real_ptr(2)
                new_day%SStartInstant%Month = array_real_ptr(3)
                new_day%SStartInstant%Day = array_real_ptr(4)
                new_day%SStartInstant%Hour = array_real_ptr(5)
                new_day%SStartInstant%Minute = array_real_ptr(6)
                new_day%SStartInstant%Second = array_real_ptr(7)
                
                new_day%SEndInstant%Year = array_real_ptr(8)
                new_day%SEndInstant%Month = array_real_ptr(9)
                new_day%SEndInstant%Day = array_real_ptr(10)
                new_day%SEndInstant%Hour = array_real_ptr(11)
                new_day%SEndInstant%Minute = array_real_ptr(12)
                new_day%SEndInstant%Second = array_real_ptr(13)
                
                new_day%Irrigation = array_real_ptr(14)
                new_day%TotalIrrigation = array_real_ptr(15)
                
                call SetDate (new_day%StartInstant,         &
                              new_day%SStartInstant%Year,   &
                              new_day%SStartInstant%Month,  &
                              new_day%SStartInstant%Day,    &
                              new_day%SStartInstant%Hour,   &
                              new_day%SStartInstant%Minute, &
                              new_day%SStartInstant%Second)

                call SetDate (new_day%EndInstant,           &
                              new_day%SEndInstant%Year,     &
                              new_day%SEndInstant%Month,    &
                              new_day%SEndInstant%Day,      &
                              new_day%SEndInstant%Hour,     &
                              new_day%SEndInstant%Minute,   &
                              new_day%SEndInstant%Second)
                              
                call SetSDate (new_day%StartDay, new_day%SStartInstant, ignore_time = .true.)
                call SetSDate (new_day%EndDay, new_day%SEndInstant, ignore_time = .true.)
            enddo
            
        endif
        
    end subroutine ConstructDailySchedules
      
    !----------------------------------------------------------------------------

    subroutine OpenInitialFile

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                             :: stat_call
        integer                             :: WorkILB, WorkIUB
        integer                             :: WorkJLB, WorkJUB
        integer                             :: HDF5_READ

        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB       
      
        call ReadFileName ('IRRIGATION_INI', Me%Files%InitialFile, "Irrigation Initial File", STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'OpenInitialFile - ModuleIrrigation - ERR010'      

        !Gets File Access Code
        call GetHDF5FileAccess (HDF5_READ = HDF5_READ)

        Me%ObjIniHDF5 = 0

        !Opens HDF5 File
        call ConstructHDF5 (Me%ObjIniHDF5,              &
                            trim(Me%Files%InitialFile), &
                            HDF5_READ, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'OpenInitialFile - ModuleIrrigation - ERR020'      
            
    end subroutine OpenInitialFile
   
    !----------------------------------------------------------------------------
   
    subroutine CloseInitialFile
   
        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: stat_call   
   
        !------------------------------------------------------------------------
      
        call KillHDF5 (Me%ObjIniHDF5, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'CloseInitialFile - ModuleIrrigation - ERR010'
   
    end subroutine CloseInitialFile
    
    
    !============================================================================
    !============================================================================
    !The objective of the ConstructTimeSeriesOutput subroutine is to pass 
    !construct the TimeSerie object that will deal with the timeseries output
    !============================================================================   
    subroutine ConstructTimeSeriesOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: list
        integer                                             :: prop_count
        integer                                             :: stat_call
        integer                                             :: schedule_i
        type(T_IrriSchedule), pointer                       :: schedule
        
        !Begin-----------------------------------------------------------------
        
        prop_count = 1
        allocate(list(prop_count))
        
        list(1)  = 'Irrigation[mm]'

        do schedule_i = 1, Me%NumberOfSchedules

            schedule => Me%Schedules(schedule_i)
            
            call StartTimeSerie (schedule%ObjTimeSeries, Me%ObjTime,                        &
                                 TimeSerieDataFile = trim(Me%TimeSeriesLocation),           &
                                 PropertyList      = list,                                  &
                                 Extension         = "sri",                                 &
                                 ResultFileName    = 'Irrigation-'//trim(schedule%ID%Name), &
                                 ModelName         = Me%ModelName,                          &
                                 STAT              = stat_call)
            if (stat_call .NE. SUCCESS_) then 
                write (*,*) '[ModuleIrrigation] Error during Timeseries creation for '//trim(schedule%ID%Name)//' schedule'
                stop 'ConstructTimeSeriesOutput - ModuleIrrigation - ERR010'
            endif
        
            schedule%AccOutput%NextOutputTime = Me%BeginTime + schedule%AccOutput%OutputInterval
            
        enddo

        !Deallocates PropertyList
        deallocate(list)
       
    end subroutine ConstructTimeSeriesOutput
   
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !MODULE API MODULE API MODULE API MODULE API MODULE API MODULE API MODULE API
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !============================================================================
    !============================================================================
    !The objective of the GetIrrigationDTForNextEvent function is to return the
    !DT for the next irrigation event. THis DT is found by the ComputeNextDT
    !subroutine at the end of the MODIFY step. THis routine looks for the nearest
    !event, looking in all the user defined schedules
    !
    !Arguments:
    !
    !id          : The ModuleIrrigation id
    !stat        : Indicates if the operation was successfull or not
    !
    !Result:
    !
    !A real number indicating the seconds (DT) to the next irrigation event
    !============================================================================  
    function GetIrrigationDTForNextEvent (id, stat) result(res)
    
        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: id
        integer, intent(OUT), optional      :: stat
        real                                :: res

        !Local-------------------------------------------------------------------
        integer                             :: stat_, ready_
        integer, save                       :: schedule_
        
        !------------------------------------------------------------------------
        
        res = 9.9E+15
        stat_ = UNKNOWN_
        
        ready_ = Ready (id) 
        
        if ((ready_ == IDLE_ERR_     ) .OR. &
            (ready_ == READ_LOCK_ERR_)) then
            
            res = Me%DTForNextEvent
                
            stat_ = SUCCESS_            
                
        else
                
            stat_ = ready_
            
        end if

        if (present(stat)) stat = stat_
        
    end function GetIrrigationDTForNextEvent
    
    !============================================================================
    !============================================================================
    !The objective of the GetIrrigationUserHeads subroutine is to pass the user
    !defined HEADS (limits) defined for each user schedule. This is used in Basin
    !to ask to PorousMedia for the Water Content for each head, that will be 
    !passed again to the Irrigation module.
    !
    !This routine must be called by basin in a looping, until the 'finish' 
    !argument returns .TRUE.
    !
    !Arguments:
    !
    !id          : The ModuleIrrigation id
    !start_head  : HEAD value to start irrigating
    !target_head : HEAD value to achieve after irrigation
    !schedule_id : ID of the user schedule, to be used when setting the water
    !            : content threshold values
    !mask        : The mask (2D) based on the application area map
    !finish      : if .TRUE. indicates that there is no more user schedules
    !stat        : Indicates if the operation was successfull or not
    !============================================================================      
    function GetIrrigationPMIsRequired (id, stat) result(res)
    
        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: id
        integer, intent(OUT), optional      :: stat
        logical                             :: res

        !Local-------------------------------------------------------------------
        integer                             :: stat_, ready_
        integer, save                       :: schedule_
        
        !------------------------------------------------------------------------
        
        res = .FALSE.
        stat_ = UNKNOWN_
        
        ready_ = Ready (id) 
        
        if ((ready_ == IDLE_ERR_     ) .OR. &
            (ready_ == READ_LOCK_ERR_)) then
            
do1:        do schedule_ = 1, Me%NumberOfSchedules
                
                if (Me%Schedules(schedule_)%IrrigationMethod /= FixedIrrigation) then
                    
                    res = .TRUE. 
                    exit do1
                    
                endif

            enddo do1
                
            stat_ = SUCCESS_            
                
        else
                
            stat_ = ready_
            
        end if

        if (present(stat)) stat = stat_
    
    end function GetIrrigationPMIsRequired
    
    !============================================================================
    !============================================================================
    !The objective of the GetIrrigationUserHeads subroutine is to pass the user
    !defined HEADS (limits) defined for each user schedule. This is used in Basin
    !to ask to PorousMedia for the Water Content for each head, that will be 
    !passed again to the Irrigation module.
    !
    !This routine must be called by basin in a looping, until the 'finish' 
    !argument returns .TRUE.
    !
    !Arguments:
    !
    !id          : The ModuleIrrigation id
    !start_head  : HEAD value to start irrigating
    !target_head : HEAD value to achieve after irrigation
    !schedule_id : ID of the user schedule, to be used when setting the water
    !            : content threshold values
    !mask        : The mask (2D) based on the application area map
    !finish      : if .TRUE. indicates that there is no more user schedules
    !stat        : Indicates if the operation was successfull or not
    !============================================================================    
    subroutine GetIrrigationThresholds (id, start_head, target_head, schedule_id, mask, finish, stat)
    
        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: id
        real, intent(OUT)                   :: start_head
        real, intent(OUT)                   :: target_head
        integer, intent(OUT)                :: schedule_id
        logical, pointer, intent(OUT)       :: mask(:,:)
        logical, intent(OUT)                :: finish
        integer, intent(OUT), optional      :: stat

        !Local-------------------------------------------------------------------
        integer                             :: stat_, ready_
        integer, save                       :: schedule_ = 0
        type(T_IrriSchedule), pointer       :: schedule
        logical                             :: must_continue
        
        !------------------------------------------------------------------------
        
        start_head = 0.0
        target_head = 0.0
        schedule_id = -1
        finish = .true.
        stat_ = UNKNOWN_
        
        must_continue = .true.
        
        ready_ = Ready (id) 
        
        if ((ready_ == IDLE_ERR_     ) .OR. &
            (ready_ == READ_LOCK_ERR_)) then
            
            do while (must_continue)
                schedule_ = schedule_ + 1
                
                if (schedule_ > Me%NumberOfSchedules) then
                    
                    finish = .true.
                    schedule_ = 0
                    must_continue = .false.
                    
                else
                    
                    schedule => Me%Schedules(schedule_)
                    
                    if (schedule%IrrigationMethod == FixedIrrigation) then
                        
                        must_continue = .true.
                        
                    else
                        
                        schedule_id = schedule_
                        start_head = schedule%HeadWiltingPoint
                        target_head = schedule%HeadTarget
                        schedule%ApplicationAreaMap%LogicalField = .false.
                        schedule%ApplicationAreaMap%LogicalField(2,2) = .true.
                        mask => schedule%ApplicationAreaMap%LogicalField
                        finish = .false.
                        must_continue = .false.
                        
                    endif
                endif
            enddo
                
            stat_ = SUCCESS_            
                
        else
                
            stat_ = ready_
            
        end if

        if (present(stat)) stat = stat_
    
    end subroutine GetIrrigationThresholds
        
    !============================================================================
    !============================================================================
    !The objective of the GetIrrigationRequirements subroutine is to pass the 
    !requirements of ModuleIrrigation to run properly, as this is used by Basin 
    !to know which data must be passed to the Irrigation module on Modify
    !
    !Arguments:
    !
    !id                 : The ModuleIrrigation id
    !roots_depth        : .TRUE. if this module needs the root depth
    !roots_klb          : .TRUE. if this module needs the root KLB
    !soil_water_content : .TRUE. if this module needs Soil Water Content
    !soil_head          : .TRUE. if this module needs Soil HEAD
    !stat               : Indicates if the operation was successfull or not
    !============================================================================ 
    subroutine GetIrrigationRequirements (id,                   &
                                          roots_depth,          &
                                          roots_klb,            &
                                          soil_water_content,   &
                                          soil_head,            &
                                          stat)
        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: id
        logical, intent(OUT)                :: roots_depth
        logical, intent(OUT)                :: roots_klb
        logical, intent(OUT)                :: soil_water_content
        logical, intent(OUT)                :: soil_head
        integer, intent(OUT), optional      :: stat

        !Local-------------------------------------------------------------------
        integer                             :: stat_, ready_

        !------------------------------------------------------------------------
        stat_ = UNKNOWN_
        
        roots_depth = .false.
        roots_klb = .false.
        soil_water_content = .false.
        soil_head = .false.

        ready_ = Ready (id) 
        
        if ((ready_ == IDLE_ERR_     ) .OR. &
            (ready_ == READ_LOCK_ERR_)) then
            
            call Read_Lock(mIrrigation_, Me%InstanceID) 

            !ToDo: Add the code here

            stat_ = SUCCESS_
            
        else
                
            stat_ = ready_
            
        end if

        if (present(STAT))STAT = STAT_
        
    end subroutine GetIrrigationRequirements
 
    !============================================================================
    !============================================================================
    !The objective of the GetIrrigationFlux subroutine is to return the Irrigation
    !matrix to be used in Basin
    !
    !Arguments:
    !
    !id                 : The ModuleIrrigation id
    !irrigation         : Irrigation matrix
    !stat               : Indicates if the operation was successfull or not
    !============================================================================ 
    subroutine GetIrrigationFlux (id,           &
                                  irrigation,   &
                                  stat)
   
        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: id
        real, intent(INOUT), pointer        :: irrigation(:,:)
        integer, intent(OUT), optional      :: stat

        !Local-------------------------------------------------------------------
        integer                             :: stat_, ready_

        !------------------------------------------------------------------------
      
        stat_ = UNKNOWN_
        irrigation => null()

        ready_ = Ready (id) 
        
        if ((ready_ == IDLE_ERR_     ) .OR. &
            (ready_ == READ_LOCK_ERR_)) then
            
            call Read_Lock(mIrrigation_, Me%InstanceID) 

            !call UpdateIrrigation
            irrigation => Me%IrrigationFlux

            STAT_ = SUCCESS_
            
        else
                
            STAT_ = ready_
            
        end if

        if (present(stat)) stat = stat_
      
    end subroutine GetIrrigationFlux
   
    !============================================================================
    !============================================================================
    !The objective of the SetIrrigationThresholds subroutine is to set in the 
    !Irrigation module the water content thresholds for each schedule
    !matrix to be used in Basin
    !
    !Arguments:
    !
    !id                 : The ModuleIrrigation id
    !schedule           : ID of the user schedule to be set
    !start_wc           : Water content to start irrigation
    !target_wc          : Target water content for irrigation
    !stat               : Indicates if the operation was successfull or not
    !============================================================================ 
    subroutine SetIrrigationThresholds (id,         &
                                        schedule,   &
                                        start_wc,   &
                                        target_wc,  &
                                        mapping,    &
                                        stat)
        
        !Arguments---------------------------------------------------------------
        integer, intent(IN)                 :: id
        integer, intent(IN)                 :: schedule
        real, intent(IN), pointer           :: start_wc(:,:,:)
        real, intent(IN), pointer           :: target_wc(:,:,:)
        integer, pointer                    :: mapping(:,:,:)
        integer, intent(OUT), optional      :: stat

        !Local-------------------------------------------------------------------
        integer                             :: stat_, ready_
        integer                             :: i, j, k
        
        !------------------------------------------------------------------------
      
        stat_ = UNKNOWN_        

        ready_ = Ready (id) 
        
        if ((ready_ == IDLE_ERR_     ) .OR. &
            (ready_ == READ_LOCK_ERR_)) then
            
            !call Read_Lock(mIrrigation_, Me%InstanceID) 

            call SetMatrixValue (Me%Schedules(schedule)%WaterContentTarget%Field3D, Me%WorkSize, target_wc)
            call SetMatrixValue (Me%Schedules(schedule)%WaterContentWiltingPoint%Field3D, Me%WorkSize, start_wc)
            
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                
                if (mapping(i,j,k) == 1) then
                    Me%Schedules(schedule)%WaterContentEasy%Field3D(i,j,k) =                    &
                             Me%Schedules(schedule)%WaterContentWiltingPoint%Field3D(i,j,k) +   &
                             (Me%Schedules(schedule)%WaterContentTarget%Field3D(i,j,k) -        &
                             Me%Schedules(schedule)%WaterContentWiltingPoint%Field3D(i,j,k)) / 2
                else
                    Me%Schedules(schedule)%WaterContentEasy%Field3D(i,j,k) = 0.0
                endif
            enddo
            enddo
            enddo
                

            STAT_ = SUCCESS_
            
        else
                
            STAT_ = ready_
            
        end if

        if (present(stat)) stat = stat_
        
    end subroutine SetIrrigationThresholds

    !----------------------------------------------------------------------------
                                    
    subroutine UnGetIrrigation2D_R4 (id, array, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: id
        real(4), dimension(:,:), pointer    :: array
        integer, intent(OUT), optional      :: STAT

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        ready_ = Ready (id)

        if (ready_ == READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mIrrigation_, Me%InstanceID, "UnGetIrrigationMelted2D_R4")

            STAT_ = SUCCESS_
      
        else               
            STAT_ = ready_
      
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetIrrigation2D_R4

    !--------------------------------------------------------------------------

    subroutine UnGetIrrigation2D_R8 (id, Array, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: id
        real(8), dimension(:,:), pointer    :: array
        integer, intent(OUT), optional      :: STAT

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, ready_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        ready_ = Ready (id)

        if (ready_ == READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mIrrigation_, Me%InstanceID, "UnGetIrrigationMelted2D_R4")

            STAT_ = SUCCESS_
      
        else               
            STAT_ = ready_
      
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetIrrigation2D_R8

    !----------------------------------------------------------------------------
        
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
    subroutine ModifyIrrigation (id, data, stat)

        !Arguments---------------------------------------------------------------
        integer                                     :: id
        type(T_IrrigationData), pointer             :: data
        integer, intent(OUT), optional              :: stat

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: schedule_i
        type(T_IrriSchedule), pointer               :: schedule_x
        
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        ready_ = Ready (id)

        if (ready_ == IDLE_ERR_) then

            if (MonitorPerformance) &
                call StartWatch ("ModuleIrrigation", "ModifyIrrigation")
                  
            Me%Data => data
            
            Me%Irrigation = 0.0
            
            call ReadLockExternalVar()   
            
            do schedule_i = 1, Me%NumberOfSchedules
            
                schedule_x => me%Schedules(schedule_i)                
                call ReadPropertiesFromFile (schedule_x)                                
                
                !First, compute irrigation needs
                select case (schedule_x%IrrigationMethod)
                case (FixedIrrigation)
                    
                    !Fixed irrigation is provided by the user.
                    !The Irrigation was already set at 'ReadPropertiesFromFile'
                    
                case (IrrigationBySteps)
                    
                    call ComputeIrrigationBySteps (schedule_x)
                    
                case (ContinuousIrrigation)
                    
                    call ComputeContinuousIrrigation (schedule_x)
                    
                case default
                
                    stop 'ModifyIrrigation - ModuleIrrigation - ERR030'
                    
                end select                            
                               
                !Remove any old (before today) schedule 
                call RemoveOldDailySchedules (schedule_x)
                
                !Update the irrigation property used by Basin
                call UpdateIrrigation (schedule_x)
                
                !Update the output (accumulated) timeseries
                call UpdateOutputTimeSeries (schedule_x)
            enddo
            
            !Restart Output File
            if (Me%Output%WriteRestartFile .AND. .NOT. (Me%Now == Me%EndTime)) then
                if(Me%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    
                    call WriteRestartFile(.false.)
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                    
                endif
            endif
            
            call ComputeNextDT

            call ReadUnLockExternalVar()
            
            STAT_ = SUCCESS_                        
            
            if (MonitorPerformance) &
                call StopWatch ("ModuleIrrigation", "ModifyIrrigation")

        else               
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyIrrigation
    
    !----------------------------------------------------------------------------
    
    subroutine ComputeNextDT
    
        !Arguments---------------------------------------------------------------
    
        !Local-------------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        type(T_DailySchedule), pointer      :: sch_day
        integer                             :: i
        real                                :: dt_for_next_event
        
        !------------------------------------------------------------------------
        Me%DTForNextEvent = 9.9E+15
        
do1:    do i = 1, Me%NumberOfSchedules
        
            schedule => Me%Schedules(i)
            sch_day => schedule%FirstDailySchedule
        
            if (associated (sch_day)) then
        
                dt_for_next_event = min(0.0, (sch_day%StartInstant - Me%Now))
                
                if (dt_for_next_event < Me%DTForNextEvent) &
                    Me%DTForNextEvent = dt_for_next_event
        
            endif
            
            if (Me%DTForNextEvent == 0.0) exit do1

        enddo do1
            
    end subroutine ComputeNextDT

    !----------------------------------------------------------------------------
    
    subroutine FindRootsKLB (schedule, roots_depth, lai_senescence)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        real, pointer, dimension(:,:)       :: roots_depth
        real, pointer, dimension(:,:)       :: lai_senescence
        
        !Local-------------------------------------------------------------------
        integer                             :: i, j, k
        real                                :: acc_depth
        
        !Begin-------------------------------------------------------------------
        
        schedule%RootsKLB = 9999999
        
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            
            if (Me%Data%BasinPoints(i,j)) then
            if (lai_senescence(i,j) < 1.0) then
                
                acc_depth = 0.0            
            
do1:            do k = Me%WorkSize%KUB, Me%WorkSize%KLB, -1
            
                    acc_depth = acc_depth + Me%Data%DWZ(i,j,k)
                                
                    if (acc_depth >= roots_depth(i,j)) then
                    
                        schedule%RootsKLB(i,j) = k
                        !print *, "Root KLB = ", k
                        exit do1
                    
                    endif
                            
                enddo do1

            endif
            endif
        
        enddo
        enddo
        
        
    end subroutine
    
    !----------------------------------------------------------------------------

    subroutine ComputeIrrigationBySteps (schedule)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        
        !Local-------------------------------------------------------------------
        real                                :: irrigation
      
        !Begin-------------------------------------------------------------------
        
        !First, check if can compute irrigation and if can, compute the 
        !irrigation needs
        if (CanComputeIrrigationNeeds(schedule)) then
        
            !print *, "Can compute irrigation needs"
            
            call FindRootsKLB (schedule, Me%Data%RootsDepth, Me%Data%LAISenescence)
            
            if (SetComputePoints (schedule, Me%ComputePoints)) then
                
                !print *, "Will check irrigation needs"
                
                irrigation = ComputeIrrigationNeed (schedule, Me%ComputePoints)
                
                if (irrigation > 0.0) then

                    !print *, "Irrigation needs is = ", irrigation
                    Me%GlobalCounter = Me%GlobalCounter + 1
                    !print *, Me%GlobalCounter
                    call SetupDailySchedule (schedule, irrigation)
            
                endif
                           
            endif
            
        endif
        
        if (.not. schedule%IsIrrigationScheduledForToday) then
            schedule%TimeSinceLastEvent = schedule%TimeSinceLastEvent + Me%DT
        endif
        
    end subroutine ComputeIrrigationBySteps
    
    !----------------------------------------------------------------------------
    
    subroutine ComputeContinuousIrrigation (schedule)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        
        !Local-------------------------------------------------------------------
      
        !Begin-------------------------------------------------------------------
        
    
    end subroutine ComputeContinuousIrrigation
    
    !============================================================================
    !============================================================================
    !The objective of the UpdateIrrigation subroutine is to update the Irrigation
    !matrix that will be used by the other modules
    !============================================================================
    subroutine UpdateIrrigation (schedule)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        
        !Local-------------------------------------------------------------------
        type(T_DailySchedule), pointer      :: sch_day
        integer                             :: i, j, index
        
        !Begin-------------------------------------------------------------------        
        
        Me%Irrigation = 0.0
        
        do index = 1, Me%NumberOfSchedules
        
            schedule => Me%Schedules(index)
            sch_day => schedule%FirstDailySchedule

            schedule%ActualIrrigation = 0.0
            
            if (associated (sch_day)) then
            if (Me%Now > sch_day%StartInstant) then
            if (Me%Now <= sch_day%EndInstant) then
        
                ![mm]                     = (     [mm/s]        *  [s] )
                schedule%ActualIrrigation = (sch_day%Irrigation * Me%DT)
                if (schedule%ActualIrrigation > 0.0) then
                    schedule%TimeSinceLastEvent = 0.0
                endif
                
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            
                    if (Me%Data%BasinPoints(i,j) == 1) then ! .AND. schedule%ApplicationAreaMap%LogicalField(i,j)) then                        
                        ![mm]              =        [mm]        + (      [mm/s]       *  [s] )
                        Me%Irrigation(i,j) = Me%Irrigation(i,j) + (sch_day%Irrigation * Me%DT)
                    else
                        Me%Irrigation(i,j) = 0.0
                    endif
                
                enddo
                enddo                                
                
            endif
            endif
            endif
        
        enddo

        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            
            if (Me%Data%BasinPoints(i,j) == 1) then ! .AND. schedule%ApplicationAreaMap%LogicalField(i,j)) then
                ![m3/s]                =       [mm]         / [mm/m] *      [m2]          /  [s]
                Me%IrrigationFlux(i,j) = Me%Irrigation(i,j) / 1000.0 * Me%Data%Areas(i,j) / Me%DT
            endif
                
        enddo
        enddo
                
    end subroutine UpdateIrrigation
    
    !============================================================================
    !============================================================================
    !The objective of the RemoveOldDailySchedules subroutine is to remove
    !daily schedules prior to the actual day.
    !
    !Arguments:
    !
    !schedule   : The user defined (via input file) schedule
    !============================================================================
    subroutine RemoveOldDailySchedules (schedule)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        
        !Local-------------------------------------------------------------------
        type(T_DailySchedule), pointer      :: sch_day
        logical                             :: quit
        logical                             :: can_delete
        
        !Begin-------------------------------------------------------------------
        quit = .false.
        
        do while (.NOT. quit)
            
            can_delete = .false.
            
            sch_day => schedule%FirstDailySchedule
do1:        do while (associated (sch_day))
                    
                if (Me%Day > sch_day%EndDay) then
                
                    can_delete = .true.
                    exit do1
                
                endif
            
                sch_day => sch_day%Next
            
            end do do1

            if (.NOT. associated (sch_day)) then
                
                quit = .true.
                
            else
                
                if (can_delete) then
                    
                    call RemoveDailyScheduleFromList (schedule, sch_day)
                    deallocate (sch_day)
                    nullify (sch_day)
                    
                endif
                
            endif
        
        end do
        
    end subroutine RemoveOldDailySchedules
    
    !============================================================================
    !============================================================================
    !The objective of the CanComputeIrrigationNeeds function is to verify if 
    !the model should compute irrigation needs. 
    !
    !Arguments:
    !
    !schedule   : The user defined (via input file) schedule
    !now        : Actual time instant
    !
    !Function results:
    !
    !res : .TRUE. if the model have to compute needs, .FALSE. otherwise
    !============================================================================
    function CanComputeIrrigationNeeds (schedule) result (res)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        logical                             :: res
        
        !Local-------------------------------------------------------------------
        
        !------------------------------------------------------------------------
        res = .false.
        schedule%IsIrrigationScheduledForToday = IsIrrigationScheduledForToday(schedule)
        
        if (.NOT. schedule%IsIrrigating) then    
        if (.NOT. schedule%IsIrrigationScheduledForToday) then
        if (schedule%TimeSinceLastEvent >= schedule%MinimumIntervalBetweenEvents) then
            
            res = .true.
            
        endif
        endif
        endif
    
    end function CanComputeIrrigationNeeds
    
    !============================================================================
    !============================================================================
    !The objective of the SetComputePoints subroutine is to set the 
    !compute_points matrix, that will be used in the ComputeIrrigationNeeds 
    !subroutine
    !
    !Arguments:
    !
    !schedule       : The user defined (via input file) schedule
    !roots_depth    : Roots depth matrix
    !compute_points : Compute Points matrix to be set
    !============================================================================
    function SetComputePoints (schedule, compute_points) result (res)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule        
        logical, dimension(:,:), pointer    :: compute_points
        logical                             :: res
        
        !Local-------------------------------------------------------------------
        integer                             :: i, j
        logical                             :: surface_is_saturated
        logical                             :: there_are_roots
        logical                             :: must_apply_at_this_location
        
        !------------------------------------------------------------------------
       
        res = .false.
        
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            
            if (Me%Data%BasinPoints(i,j) == 1) then
                
                schedule%ApplicationAreaMap%LogicalField(i,j) = .true.
                must_apply_at_this_location = schedule%ApplicationAreaMap%LogicalField(i,j)
                
                if (Me%Data%RootsDepth(i,j) > 0.0) then
                    there_are_roots = .true.
                else
                    there_are_roots = .false.
                endif
                
                surface_is_saturated = IsSurfaceSaturated(schedule, i, j)
                
                !print *, "Before: ", schedule%ApplicationAreaMap%LogicalField(i,j), Me%Data%RootsDepth(i,j)
                !print *, "Compute points: ", must_apply_at_this_location, there_are_roots, surface_is_saturated
                              
                if (must_apply_at_this_location .AND. there_are_roots .AND. (.NOT. surface_is_saturated)) then
                    compute_points(i,j) = .true.
                    res = .true.
                else
                    compute_points(i,j) = .false.
                endif               

            else
                compute_points(i,j) = .false.
            endif
        enddo
        enddo
        
    end function
    
    !============================================================================
    !============================================================================
    !The objective of the IsIrrigationScheduledForToday function is to verify if 
    !there is an irrigation already scheduled to the actual day. 
    !
    !Arguments:
    !
    !schedule   : The user defined (via input file) schedule
    !now        : Actual time instant
    !
    !Function results:
    !
    !res : .TRUE. if there is a scheduled irrigation, .FALSE. otherwise
    !============================================================================
    function IsIrrigationScheduledForToday (schedule) result (res)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        logical                             :: res
        
        !Local-------------------------------------------------------------------
        type(T_DailySchedule), pointer      :: sch_day
        
        !------------------------------------------------------------------------
        res = .false.
        
        sch_day => schedule%FirstDailySchedule
do1:    do while (associated (sch_day))
                    
            if (Me%Day == sch_day%StartDay) then
                
                res = .true.
                exit do1
                
            endif
            
            if (Me%Day == sch_day%EndDay) then

                res = .true.
                exit do1
                
            endif            
            
            sch_day => sch_day%Next
            
        end do do1
    
    end function IsIrrigationScheduledForToday
    
    !============================================================================
    !============================================================================
    !The objective of the IsSurfaceSaturated routine is to verify if the surface
    !is with too much water, to avoid start irrigation that will cause runoff
    !
    !Arguments:
    !
    !schedule       : The user defined (via input file) schedule
    !i, j           : Cell location
    !k_bottom, kub  : Layers limits for this cell
    !root_depth     : Roots depth
    !dwz            : Cell height (center) matrix
    !rwc            : Relative Water Content matrix
    !
    !Result:
    !
    !res            : .TRUE. if it is saturated, .FALSE. otherwise
    !============================================================================
    function IsSurfaceSaturated (schedule, i, j) result(res)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        integer                             :: i, j
        logical                             :: res
        
        !Local-------------------------------------------------------------------
        integer                             :: k
        real                                :: saturated_height
        real                                :: acc_depth
        real                                :: saturated_fraction
        
        !------------------------------------------------------------------------
        
        saturated_height = 0.0
        acc_depth = 0.0
        res = .false.
        
        if (Me%Data%RootsDepth(i,j) > 0.0) then
        
do1:        do k = Me%WorkSize%KUB, schedule%RootsKLB(i,j), -1            
    
                if (acc_depth > schedule%MaxDepthToCheckSaturation) &
                    exit do1
            
                if (acc_depth > Me%Data%RootsDepth(i,j)) &
                    exit do1
            
                if (Me%Data%SoilRelativeWaterContent(i,j,k) >= schedule%SaturationThreshold) &
                    saturated_height = saturated_height + Me%Data%DWZ(i,j,k)
                
                acc_depth = acc_depth + Me%Data%DWZ(i,j,k)
            
            end do do1
        
            if (saturated_height > 0.0 .AND. acc_depth > 0.0) then
            
                saturated_fraction = saturated_height / acc_depth
            
                if (saturated_fraction > schedule%MaxSaturatedFraction) &
                    res = .true.
        
            endif
            
        endif
        
    end function IsSurfaceSaturated
    
    !============================================================================
    !============================================================================
    !The objective of the SetupDailySchedule routine is to setup the 
    !calculated irrigation needs, distributing it through the current and next 
    !days accordinly the irrigation gear characteristics and user options
    !
    !Arguments:
    !
    !schedule   : The user defined (via input file) schedule
    !irrigation : The ammount of irrigation calculated
    !now        : Actual time instant
    !============================================================================
    subroutine SetupDailySchedule (schedule, irrigation)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        real                                :: irrigation
        
        !Local-------------------------------------------------------------------
        integer                             :: days_count
        logical                             :: finished
        real                                :: rest
        type(T_DailySchedule), pointer      :: new_day
        type(T_Time)                        :: time_aux
        type(T_Time)                        :: start_instant
        type(T_STime)                       :: sstart_instant
        real                                :: duration
        real                                :: irrigation_
        real                                :: hours
        !------------------------------------------------------------------------
        !Initialize local variables
        days_count = 0
        finished = .false.
        rest = irrigation
        start_instant = Me%Now - Me%DT
        call SetSDate (sstart_instant, start_instant) 
        
        !Daily Schedule Looping
        do while (.NOT. finished)
            
            if (rest > schedule%MinimumToIrrigate) then
                                
                days_count = days_count + 1
                
                if (days_count <= schedule%MaxConsecutiveDays) then
                    
                    !Creates a new schedule day
                    call AddDailyScheduleToList (schedule, new_day)
                    
                    !Compute the irrigation start instant
                    if (sstart_instant%Hour >= schedule%StartInstantThreshold) then
                        if (sstart_instant%Hour <= schedule%EndInstantThreshold) then
                            
                            new_day%StartInstant = start_instant
                            
                        else
                                                        
                            call SetDate(start_instant, sstart_instant%Year, sstart_instant%Month, sstart_instant%Day, schedule%StartInstantThreshold, 0.0, 0.0)
                            new_day%StartInstant = start_instant + 86400.0
                            
                        endif
                    else
                                                
                        call SetDate(start_instant, Me%SNow%Year, Me%SNow%Month, Me%SNow%Day, schedule%StartInstantThreshold, 0.0, 0.0)
                        new_day%StartInstant = start_instant
                        
                    endif
                    
                    call SetSDate (new_day%SStartInstant, new_day%StartInstant)
                    call SetSDate (new_day%StartDay, new_day%SStartInstant, ignore_time = .true.)
                    
                    ![h]  = [mm] /      [mm/h]
                    hours = rest / schedule%GearDebit

                    ![s] = min([h] * [s/h], [s])
                    duration = min (hours * 3600.0, schedule%MaxDailyIrrigationTime)

                    !Corrects the hours to the limit if necessary
                    ![h]  =   [s]    / [s/h]
                    hours = duration / 3600.0 
                    
                    !![mm]      =       [mm/h]       *  [h]
                    irrigation_ = schedule%GearDebit * hours
                    
                    ![mm/s]            =     [mm]    /   [s]
                    new_day%Irrigation = irrigation_ / duration
                    
                    !Compute the instant to stop irrigation ofr this day
                    new_day%EndInstant = new_day%StartInstant + duration
                    
                    call SetSDate (new_day%SEndInstant, new_day%EndInstant)
                    call SetSDate (new_day%EndDay, new_day%SEndInstant, ignore_time = .true.)
                    
                    start_instant = start_instant + 86400.0
                    call SetSDate (sstart_instant, start_instant)
                                                          
                    !write (*, 100) int(new_day%SStartInstant%Year), int(new_day%SStartInstant%Month), int(new_day%SStartInstant%Day), int(new_day%SStartInstant%Hour), int(new_day%SStartInstant%Minute), int(new_day%SStartInstant%Second)
                    !write (*, 110) int(new_day%SEndInstant%Year), int(new_day%SEndInstant%Month), int(new_day%SEndInstant%Day), int(new_day%SEndInstant%Hour), int(new_day%SEndInstant%Minute), int(new_day%SEndInstant%Second)
                    !print *, 'Initial           (mm)  : ', rest
                    !print *, 'Irrigation Volume (mm)  : ', irrigation_

                    !Compute the ammount of irrigation to apply next day
                    ![mm]= [mm] - [mm]
                    rest = rest - irrigation_                    
                    
                    !print *, 'Rest                    : ', rest
                    
                    !100 format(1X, "Irrigation Start Instant: ",(i4,":"),4(i2, ":"), i2)
                    !110 format(1X, "Irrigation End Instant  : ",(i4,":"),4(i2, ":"), i2)
                    
                    if (rest <= 0.0) then
                        
                        finished = .true.
                        
                    endif
                    
                else
                    
                    finished = .true.
                    
                endif
                
            else
                
                finished = .true.
                
            endif
            
        enddo
        
    end subroutine SetupDailySchedule
    
    !============================================================================
    !============================================================================
    !The objective of the ComputeIrrigationNeed routine computes the ammount of
    !irrigation required for the selected area defined in the Schedule
    !
    !Arguments:
    !
    !schedule       : The user defined (via input file) schedule
    !compute_points : Points to compute
    !
    !Result:
    !
    !res            : Ammount of water to irrigate
    !============================================================================
    function ComputeIrrigationNeed (schedule, compute_points) result(res)
    
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        logical, dimension(:,:), pointer    :: compute_points
        real                                :: res
        
        !Local-------------------------------------------------------------------
        integer                             :: i, j, k
        real                                :: acc_to_irrigate
        real                                :: acc_area
        real                                :: defice
        real                                :: iv
        real                                :: ewc, wc, twc
        
        !------------------------------------------------------------------------
        
        acc_to_irrigate = 0.0
        acc_area = 0.0
        res = 0.0

        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            
            if (compute_points(i,j)) then
            
                defice = 0.0
                
                do k = schedule%RootsKLB(i,j), Me%WorkSize%KUB
                    
                    ewc = schedule%WaterContentEasy%Field3D(i,j,k)
                    wc = Me%Data%SoilWaterContent(i,j,k)
                    twc = schedule%WaterContentTarget%Field3D(i,j,k)
                    
                    !print *, "For k = ", k, ewc, wc, twc
                    
                    if (ewc >= wc) then
                        
                        defice = defice + (twc - wc)
                        !print *, "Defice (acc) = ", defice
                        
                    endif
            
                end do
        
                if (defice > 0.0) then
                    
                    acc_to_irrigate = acc_to_irrigate + max(0.0, defice * Me%Data%RootsDepth(i,j) * 1000)
                    acc_area = acc_area + Me%Data%Areas(i,j)
                    
                    !print *, "Acc to irrigate = ", acc_to_irrigate
                    !print *, "Acc Area = ", acc_area
                    
                endif
                
            endif
            
        enddo
        enddo
        
        if (acc_area > 0.0 .and. acc_area >= schedule%MinimumAreaToStartIrrigation) then
            
            iv = acc_to_irrigate / acc_area
            if (iv >= schedule%MinimumToIrrigate) then
                !print *, "Irrigation ", iv, " >= minimum of ", schedule%MinimumToIrrigate
                res = iv
            endif
        endif
        
    end function ComputeIrrigationNeed
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !DATA INPUT DATA INPUT DATA INPUT DATA INPUT DATA INPUT DATA INPUT DATA INPUT
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    
    !============================================================================
    !============================================================================
    !The objective of the ReadPropertiesFromFile routine is to update any property set
    !to be read from FILE (HDF, TIMESERIES, etc)
    !
    !Arguments:
    !
    !schedule       : The user defined (via input file) schedule
    !============================================================================
    subroutine ReadPropertiesFromFile (schedule)
 
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        
        !Local-------------------------------------------------------------------
        type (T_IrriProperty), pointer      :: property_x
        integer                             :: stat_call
      
        !Begin-------------------------------------------------------------------
        
        property_x => schedule%FirstProperty
        do while (associated(property_x))

            if (property_x%ID%SolutionFromFile) then

                call ModifyFillMatrix (FillMatrixID   = property_x%ID%ObjFillMatrix,    &
                                       Matrix2D       = property_x%Field,               &
                                       PointsToFill2D = Me%Data%BasinPoints,          &
                                       STAT           = stat_call)
                if (stat_call /= SUCCESS_) then
                    write (*,*) "[ModuleIrrigation] ATTENTION"
                    write (*,*) "[ModuleIrrigation] Was not possible to read property '", trim(property_x%ID%Name), "' from file."  
                    stop 'ReadPropertiesFromFile - ModuleIrrigation - ERR010'
                endif
               
            endif
                      
            property_x => property_x%Next
        
        enddo

    end subroutine ReadPropertiesFromFile

    !============================================================================
    !============================================================================
    !The objective of the ReadLockExternalVar is to get a reference to a set of
    !arrays and other variables from other modules to be used during the 
    !ModuleIrrigation MODIFY step. In the case of the arrays, usually the other
    !module GET routine will execute a LOCK that must be released at the end of
    !the MODIFY step (through the ReadUnLockExternalVar subroutine)
    !============================================================================    
    subroutine ReadLockExternalVar ()
        
        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                             :: stat_call

        !------------------------------------------------------------------------
        !Time Stuff
        call GetComputeCurrentTime  (Me%ObjTime, Me%Now, STAT = stat_call)
        if (stat_call /= SUCCESS_) stop 'ReadLockExternalVar - ModuleIrrigation - ERR01'
      
        call SetSDate (Me%SNow, Me%Now)      
        call SetSDate (Me%SDay, Me%Now, ignore_time = .true.)
        call SetSDate (Me%Day, Me%SDay)

        call GetComputeTimeStep     (Me%ObjTime, Me%DT, STAT = stat_call)
        if (stat_call /= SUCCESS_) stop 'ReadLockExternalVar - ModuleIrrigation - ERR02'

        !Gets Basin Points
        !call GetBasinPoints (Me%ObjBasinGeometry, Me%Data%BasinPoints, STAT = stat_call)
        !if (stat_call /= SUCCESS_) stop 'ReadLockExternalVar - ModuleIrrigation - ERR03'

    end subroutine ReadLockExternalVar
   
    !============================================================================
    !============================================================================
    !The objective of the ReadUnLockExternalVar is to release the locks generated
    !during the previous call to ReadLockExternalVar
    !============================================================================
    subroutine ReadUnLockExternalVar()
        
        !Arguments---------------------------------------------------------------
        
        !Local-------------------------------------------------------------------
        integer                                     :: stat_call

        !------------------------------------------------------------------------
        !Unget Basin Points
        !call UnGetBasin (Me%ObjBasinGeometry, Me%Data%BasinPoints, STAT = stat_call)
        !if (stat_call /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleIrrigation - ERR01'
        
    end subroutine ReadUnLockExternalVar
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT OUTPUT
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    subroutine UpdateOutputTimeSeries (schedule)
        
        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule

        !Local-------------------------------------------------------------------
        integer                             :: stat_call

        !------------------------------------------------------------------------

        schedule%AccOutput%Value = schedule%AccOutput%Value + schedule%ActualIrrigation
        
        if (Me%Now .ge. schedule%AccOutput%NextOutputTime .or. Me%Now .ge. Me%EndTime) then

            Me%AccIrrigation(1) = schedule%AccOutput%Value
            
            call WriteTimeSerieLineNow(schedule%ObjTimeSeries,      &
                                       DataLine = Me%AccIrrigation, &
                                       STAT = stat_call)
            if (stat_call /= SUCCESS_) stop 'ModuleDrainageNetwork - OutputIntFlow - ERR01'                              
            
            schedule%AccOutput%NextOutputTime = schedule%AccOutput%NextOutputTime + &
                                                schedule%AccOutput%OutputInterval
        
            schedule%AccOutput%Value = 0.0
            
        endif
        
    end subroutine UpdateOutputTimeSeries
    
    !============================================================================
    !============================================================================
    !The objective of the WriteRestartFile routine is to write restart files, both
    !the last restart file (final file) as well as intermediate restart files
    !
    !Arguments:
    !
    !is_last        : If .TRUE. tells the routine this is the last restart file.
    !============================================================================
    subroutine WriteRestartFile (is_last)
              
        !Arguments
        logical                             :: is_last
        
        !Local-------------------------------------------------------------------
        type (T_IrriProperty), pointer      :: property_x
        type (T_IrriSchedule), pointer      :: schedule_x
        integer                             :: schedule_i, daily_sch_i
        integer                             :: stat_call      
        integer                             :: HDF5_CREATE
        character(LEN = PathLength)         :: fileName
        integer                             :: ObjHDF5
        real, dimension(6), target          :: AuxTime
        real, dimension(:), pointer         :: TimePtr
        type (T_Time)                       :: Actual
        real, dimension(2), target          :: schedule_array
        real, dimension(15), target         :: daily_schedule_array
        real, dimension(:), pointer         :: aux
        type(T_DailySchedule), pointer      :: day_schedule
        
        !Begin-------------------------------------------------------------------
        !Gets a pointer to Topography
        call GetGridData (Me%ObjGridData, Me%Data%Topography, STAT = stat_call)
        if (stat_call /= SUCCESS_) stop 'WriteRestartFile - ModuleIrrigation - ERR010'

        call GetBasinPoints   (Me%ObjBasinGeometry, Me%Data%BasinPoints, STAT = stat_call)
        if (stat_call /= SUCCESS_) stop 'WriteRestartFile - ModuleIrrigation - ERR020'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        if (is_last .or. Me%Output%RestartOverwrite) then
            filename = trim(Me%Files%FinalFile)
        else
            fileName = ChangeSuffix(Me%Files%FinalFile, &
                       "_"//trim(TimeToString(Me%Now))//".fin")
        endif

        ObjHDF5 = 0
        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5,        &
                            trim(filename), &
                            HDF5_CREATE, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR030'

        Actual = Me%Now
         
        call ExtractDate (Actual, AuxTime(1), AuxTime(2), AuxTime(3), &
                          AuxTime(4), AuxTime(5), AuxTime(6))
        !Writes Time
        TimePtr => AuxTime
        call HDF5SetLimits (ObjHDF5, 1, 6, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR040'

        call HDF5WriteData (ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                            Array1D = TimePtr, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR050'

        !Sets limits for next write operations
        call HDF5SetLimits (ObjHDF5,           &
                            Me%WorkSize%ILB,   &
                            Me%WorkSize%IUB,   &
                            Me%WorkSize%JLB,   &
                            Me%WorkSize%JUB,   &
                            STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR060'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, ObjHDF5, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR070'

        !Writes the Grid
        call HDF5WriteData (ObjHDF5, "//Grid", "Topography", "m", &
                            Array2D = Me%Data%Topography, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR080'

        !WriteBasinPoints
        call HDF5WriteData (ObjHDF5, "//Grid", "BasinPoints", "-", &
                            Array2D = Me%Data%BasinPoints, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR090'
        
        !Write each schedule to the file
        do schedule_i = 1, Me%NumberOfSchedules
            
            schedule_x => Me%Schedules(schedule_i)
            
            call HDF5SetLimits  (ObjHDF5, 1, 2, STAT = stat_call)
            if (stat_call /= SUCCESS_) stop 'WriteRestartFile - ModuleIrrigation - ERR090'
            
            aux => schedule_array
            
            aux(1) = schedule_x%ID%IDNumber
            aux(2) = schedule_x%DailySchedulesNumber

            call HDF5WriteData (ObjHDF5,                                    &
                                "/Schedules/"//trim(schedule_x%ID%Name),    &
                                "info",                                     &
                                "-",                                        &
                                Array1D = aux,                              &
                                STAT = stat_call)
            if (stat_call /= SUCCESS_) &
                stop 'WriteRestartFile - ModuleIrrigation - ERR110'

            !if there are active daily schedules, save them to the file
            if (schedule_x%DailySchedulesNumber > 0) then
                
                call HDF5SetLimits  (ObjHDF5, 1, 15, STAT = stat_call)
                if (stat_call /= SUCCESS_) stop 'WriteRestartFile - ModuleIrrigation - ERR090'
                
                aux => daily_schedule_array
                
                daily_sch_i = 1
                day_schedule => schedule_x%FirstDailySchedule
                do while (associated(day_schedule))
                
                    aux(1) = real(daily_sch_i)
                    
                    aux(2) = day_schedule%SStartInstant%Year
                    aux(3) = day_schedule%SStartInstant%Month
                    aux(4) = day_schedule%SStartInstant%Day
                    aux(5) = day_schedule%SStartInstant%Hour
                    aux(6) = day_schedule%SStartInstant%Minute
                    aux(7) = day_schedule%SStartInstant%Second
                    
                    aux(8)  = day_schedule%SEndInstant%Year
                    aux(9)  = day_schedule%SEndInstant%Month
                    aux(10) = day_schedule%SEndInstant%Day
                    aux(11) = day_schedule%SEndInstant%Hour
                    aux(12) = day_schedule%SEndInstant%Minute
                    aux(13) = day_schedule%SEndInstant%Second
            
                    aux(14) = day_schedule%Irrigation
                    aux(15) = day_schedule%TotalIrrigation
                    
                    call HDF5WriteData (ObjHDF5,                                                        &
                                        "/Schedules/"//trim(schedule_x%ID%Name)//"/daily schedules",    &
                                        "daily schedule",                                               &
                                        "-",                                                            &
                                        Array1D = aux,                                                  &
                                        OutputNumber = daily_sch_i,                                     &
                                        STAT = stat_call)
                    if (stat_call /= SUCCESS_) &
                        stop 'WriteRestartFile - ModuleIrrigation - ERR110'
                    
                    daily_sch_i = daily_sch_i + 1
                    day_schedule => day_schedule%Next
                    
                enddo
            
            endif
            
            if (schedule_x%IrrigationMethod == FixedIrrigation) cycle
            
            !Loop through the properties to find the ones to save to the final file
            property_x => schedule_x%FirstProperty
            do while (associated(property_x))

                !If set to be saved in the final file, write the property to the file
                if (property_x%SaveToFinalFile) then
                    call HDF5SetLimits (ObjHDF5,         &
                                        Me%WorkSize%ILB, &
                                        Me%WorkSize%IUB, &
                                        Me%WorkSize%JLB, &
                                        Me%WorkSize%JUB, &
                                        STAT = stat_call)
                    if (stat_call /= SUCCESS_) &
                        stop 'WriteRestartFile - ModuleIrrigation - ERR100'

                    call HDF5WriteData (ObjHDF5,                                            &
                                    "/Results/"//trim(schedule_x%ID%Name)//"/properties",   &
                                    trim(property_x%ID%Name),                               &
                                    trim(property_x%ID%Units),                              &
                                    Array2D = property_x%Field,                             &
                                    STAT = stat_call)
                    if (stat_call /= SUCCESS_) &
                        stop 'WriteRestartFile - ModuleIrrigation - ERR110'
                
                endif
                
                property_x => property_x%Next

            enddo
            
        enddo

        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR030'

        !Unget
        call UnGetBasin (Me%ObjBasinGeometry, Me%Data%BasinPoints, stat_call)
        if (stat_call /= SUCCESS_) stop 'WriteRestartFile - ModuleIrrigation - ERR90'  

        !UnGets Topography
        call UnGetGridData (Me%ObjGridData, Me%Data%Topography, STAT = stat_call)
        if (stat_call /= SUCCESS_) &
            stop 'WriteRestartFile - ModuleIrrigation - ERR100'

    end subroutine WriteRestartFile

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !LISTS MANAGEMENT LISTS MANAGEMENT LISTS MANAGEMENT LISTS MANAGEMENT LISTS MA
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !============================================================================
    !============================================================================
    !The objective of the AddDailyScheduleToList routine is to add a daily 
    !schedule previously created to the daily schedule list. It will add it to 
    !the end of the list.
    !
    !Arguments:
    !
    !schedule       : The user defined (via input file) schedule
    !daily_schedule : The daily schedule to be add to the list
    !============================================================================
    subroutine AddDailyScheduleToList (schedule, daily_schedule)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer               :: schedule
        type(T_DailySchedule), pointer, intent(out) :: daily_schedule
        
        !Local-------------------------------------------------------------------
        integer                                     :: stat_call

        !------------------------------------------------------------------------

        allocate (daily_schedule, STAT=stat_call)
        if (stat_call /= SUCCESS_) then                        
            print *, '[ModuleIrrigation] Memory allocation error.'
            stop 'SetupIrrigationSchedule - ModuleIrrigation - ERR010'
        endif
        
        if (.NOT. associated(schedule%FirstDailySchedule)) then
            
            schedule%DailySchedulesNumber = 1
            schedule%FirstDailySchedule => daily_schedule
            schedule%LastDailySchedule => daily_schedule
            
        else
            
            daily_schedule%Prior => schedule%LastDailySchedule
            schedule%LastDailySchedule%Next => daily_schedule
            schedule%LastDailySchedule => daily_schedule
            schedule%DailySchedulesNumber = schedule%DailySchedulesNumber + 1
            
        end if 

    end subroutine AddDailyScheduleToList
    
    !============================================================================
    !============================================================================
    !The objective of the RemoveDailyScheduleFromList routine is to remove a  
    !daily schedule from the daily schedule list. 
    !
    !Arguments:
    !
    !schedule       : The user defined (via input file) schedule
    !daily_schedule : The daily schedule to be removed from the list
    !============================================================================
    subroutine RemoveDailyScheduleFromList (schedule, daily_schedule)

        !Arguments---------------------------------------------------------------
        type(T_IrriSchedule), pointer       :: schedule
        type(T_DailySchedule), pointer      :: daily_schedule

        !------------------------------------------------------------------------
    
        if (.NOT. associated(schedule%FirstDailySchedule)) then
            
            print *, '[ModuleIrrigation] Trying to remove a Daily Schedule from an EMPTY list'
            print *, '[ModuleIrrigation] Please, contact the module developer team.'
            stop 'RemoveDailyScheduleFromList - ModuleIrrigation - ERR010'
            
        else
            if (associated (daily_schedule%Next)) then
                daily_schedule%Next%Prior => daily_schedule%Prior
            else
                schedule%LastDailySchedule => daily_schedule%Prior
            endif
            
            if (associated (daily_schedule%Prior)) then
                daily_schedule%Prior%Next => daily_schedule%Next
            else
                schedule%FirstDailySchedule => daily_schedule%Next
            endif
            
            schedule%DailySchedulesNumber = schedule%DailySchedulesNumber - 1
            
            if (.NOT. associated(schedule%FirstDailySchedule)) then
                schedule%IsIrrigating = .false.
            endif
        end if 
        
    end subroutine RemoveDailyScheduleFromList

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillIrrigation (id, stat)

        !Arguments---------------------------------------------------------------
        integer                             :: id
        integer, optional, intent(OUT)      :: stat

        !Local-------------------------------------------------------------------
        integer                             :: ready_              
        integer                             :: stat_, nUsers
        integer                             :: i
        integer                             :: STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        ready_ = Ready (id)
        
        if (ready_ /= OFF_ERR_) then

            nUsers = DeassociateInstance(mIrrigation_,  Me%InstanceID)

            if (nUsers == 0) then

                !Writes file with final condition
                call WriteRestartFile(.TRUE.)

                do i = 1, Me%NumberOfSchedules
                    
                    if (Me%Schedules(i)%ObjTimeSeries /= 0) then
                        call KillTimeSerie(Me%Schedules(i)%ObjTimeSeries, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - ModuleIrrigation - ERR010'
                    endif
                    
                enddo
                
                if (Me%OutPut%Yes) then
                
                    call KillHDF5 (Me%ObjHDF5, STAT = stat_)
                    if (stat_ /= SUCCESS_) &
                        stop 'KillIrrigation - ModuleIrrigation - ERR020'
               
                endif
                
                !Deassociates External Instances
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) &
                    stop 'KillIrrigation - Irrigation - ERR030'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) &
                    stop 'KillIrrigation - Irrigation - ERR040'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjGridData)
                if (nUsers == 0) &
                    stop 'KillIrrigation - Irrigation - ERR050'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) &
                    stop 'KillIrrigation - Irrigation - ERR060'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) &
                    stop 'KillIrrigation - Irrigation - ERR070'

                !Deallocates Instance
                call DeallocateInstance ()

                id = 0
                stat_ = SUCCESS_

            endif
            
        end if

        if (present(stat)) stat = stat_

    end subroutine KillIrrigation

    !----------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        type (T_Irrigation), pointer        :: obj_aux
        type (T_Irrigation), pointer        :: obj_prev

        !------------------------------------------------------------------------        
        !Updates pointers
        if (Me%InstanceID == FirstObjIrrigation%InstanceID) then
            
            FirstObjIrrigation => FirstObjIrrigation%Next
            
        else
            
            obj_prev => FirstObjIrrigation
            obj_aux => FirstObjIrrigation%Next
        
            do while (obj_aux%InstanceID /= Me%InstanceID)
                
                obj_prev => obj_aux
                obj_aux => obj_aux%Next
                
            enddo

            !Now update linked list
            obj_prev%Next => obj_aux%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify (Me) 

    end subroutine DeallocateInstance   
   
    !----------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !MODULE MANAGEMENT MODULE MANAGEMENT MODULE MANAGEMENT MODULE MANAGEMENT MODU
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !----------------------------------------------------------------------------

    function Ready (id) result (res)

        !Arguments---------------------------------------------------------------
        integer                             :: id
        integer                             :: res

        !------------------------------------------------------------------------

        nullify (Me)

        if (id > 0) then
            
            call LocateObjIrrigation (id)
            res = VerifyReadLock (mIrrigation_, Me%InstanceID)
            
        else
            
            res = OFF_ERR_
      
        end if

        !------------------------------------------------------------------------

    end function Ready

    !----------------------------------------------------------------------------

    subroutine LocateObjIrrigation (id)

        !Arguments---------------------------------------------------------------
        integer                             :: id

        !Local-------------------------------------------------------------------

        Me => FirstObjIrrigation
        do while (associated (Me))
            
            if (Me%InstanceID == id) exit
            Me => Me%Next
            
        enddo

        if (.NOT. associated(Me)) &
            stop 'ModuleIrrigation - LocateObjIrrigation - ERR010'

    end subroutine LocateObjIrrigation
   
end module ModuleIrrigation   