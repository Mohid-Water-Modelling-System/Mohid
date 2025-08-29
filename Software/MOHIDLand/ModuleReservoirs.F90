!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Module Reservoirs
! PROJECT       : Mohid Land
! MODULE        : Reservoirs
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2015
! REVISION      : David
! DESCRIPTION   : Module to compute reservoir dynamics (0D and linked to DN)
!
!------------------------------------------------------------------------------


Module ModuleReservoirs

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleDischarges
    use ModuleFunctions,    only: LinearInterpolation, ConstructPropertyID,      &
                                  SetMatrixValue, TimeToString, ChangeSuffix,    &
                                  InterpolateValueInTime
    use ModuleTimeSerie,    only: StartTimeSerie, WriteTimeSerieLine, KillTimeSerie, &
                                  StartTimeSerieInput, GetTimeSerieValue
    use ModuleHDF5
    
    use ModuleHorizontalGrid, only: ConstructHorizontalGrid, KillHorizontalGrid, &
                                    GetHorizontalGridSize, UngetHorizontalGrid,  &
                                    WriteHorizontalGrid
    use ModuleGridData,       only: ConstructGridData, KillGridData, GetGridData, &
                                    UngetGridData
    use ModuleBasinGeometry,  only: ConstructBasinGeometry, KillBasinGeometry,   &
                                    GetBasinPoints, UngetBasin
    
    use ModuleFillMatrix,     only: ConstructFillMatrix, KillFillMatrix
    
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructReservoirs
    private ::      AllocateInstance

    !Selector
    public  :: UnGetReservoirs
    public  :: GetReservoirsnProperties
    public  :: GetReservoirsPropertiesIDByIdx
    public  :: GetReservoirsConcentration
    public  :: GetNumberOfReservoirs
    public  :: GetReservoirsNodeIDs    
    public  :: SetReservoirsInflow           !Reservoirs gets inflow from DN (DN outflow)
    public  :: GetReservoirsOutflow         !DrainageNetwork gets outflow from Reservoirs (DN inflow)
    public  :: SetDNConcReservoirs          !Reservoirs gets the drainage network concentrations
    public  :: CheckReservoirProperty
    
    !Modifier
    public  :: ModifyReservoirs

    !Destructor
    public  :: KillReservoirs                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjReservoirs 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetReservoirs1D_I
    private :: UnGetReservoirs1D_R8
    interface  UnGetReservoirs
        module procedure UnGetReservoirs1D_I
        module procedure UnGetReservoirs1D_R8
    end interface  UnGetReservoirs

    !Parameter-----------------------------------------------------------------
    !Operation - State variable - Outflow Estimate
    integer,    parameter                             :: Operation_Level_Outflow_        = 1
    integer,    parameter                             :: Operation_Level_PercInflow_     = 2
    integer,    parameter                             :: Operation_PercVol_Outflow_      = 3    
    integer,    parameter                             :: Operation_PercVol_PercInflow_   = 4
    integer,    parameter                             :: Operation_PercVol_PercMaxOutflow_   = 5
    
	!Options for maximum volume definition
    integer,    parameter                             :: UniqueValue_                    = 1
    integer,    parameter                             :: MonthlyValues_                  = 2   
    !Initial Conditions
    integer,    parameter                             :: StartPercentageFull_            = 1
    !integer,    parameter                             :: StartMinimumVolume_             = 2
    !integer,    parameter                             :: StartMaximumVolume_             = 3     
    
    !method to compute concentration
    integer,    parameter                             :: InstantMixing_                  = 1
    integer,    parameter                             :: RetentionTimeMixing_            = 2
    
    !Types---------------------------------------------------------------------


    type       T_ExtVar
        real                                        :: DT         = null_real
        type(T_Time     )                           :: BeginTime
        type(T_Time     )                           :: EndTime
        type(T_Time     )                           :: ActualTime        
        character(len=Pathlength)                   :: TopographicFile = null_str        
        type (T_Size2D)                             :: Size2D, WorkSize2D
    end type T_ExtVar    
    
    type T_TimeSerie
        integer                                                 :: ObjEnterData = 0
        character(PathLength)                                   :: Location     = null_str   
        integer                                                 :: nNodes       = 0
        integer                                                 :: nProp        = 0
        integer                 , dimension (:), pointer        :: ObjTimeSerie => null()
        integer                 , dimension (:), pointer        :: ReservoirIDs => null() 
        character(StringLength) , dimension (:), pointer        :: Name         => null()
        real                    , dimension (:), pointer        :: X            => null()
        real                    , dimension (:), pointer        :: Y            => null()
        real, dimension(:), pointer                             :: DataLine     => null()
        
    end type T_TimeSerie    
    
    type       T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime            => null()
         type (T_Time), dimension(:), pointer       :: RestartOutTime     => null()
         logical                                    :: HDF                = .false.
         logical                                    :: HDFActive          = .false.
         integer                                    :: Number             = null_int 
         logical                                    :: TimeSerie          = .false.  
         integer                                    :: NextOutPut         = null_int
         logical                                    :: WriteRestartFile     = .false.       
         logical                                    :: RestartOverwrite     = .false.
         integer                                    :: NextRestartOutput    = 1           
    end type T_OutPut    
    
    type T_ComputeOptions
        logical                                     :: Discharges             = .false.
        logical                                     :: T90_Decay                = .false.
        logical                                     :: Generic_Decay            = .false.
        logical                                     :: SurfaceFluxes            = .false.
        logical                                     :: BottomFluxes             = .false.
        logical                                     :: Erosion                  = .false.
        logical                                     :: Deposition               = .false.
        logical                                     :: WaterQuality             = .false.
        logical                                     :: Benthos                  = .false.
        logical                                     :: CeQualW2                 = .false.
        logical                                     :: Life                     = .false.        
        logical                                     :: MinConcentration       = .false.   
        logical                                     :: WarnOnNegativeValues     = .false.        
        logical                                     :: SumTotalConc             = .false.
        logical                                     :: TimeSerie                = .false.
        logical                                     :: HDF                      = .false.
        logical                                     :: DTIntervalAssociated     = .false.
        integer                                     :: PropertyComputeMethod    = InstantMixing_
    end type T_ComputeOptions    
    
    type        T_Property
        type (T_PropertyID)                         :: ID
        type (T_ComputeOptions)                     :: ComputeOptions

        !Concentrations
        real, dimension (:), pointer                :: Concentration            => null()
        real, dimension (:), pointer                :: ConcentrationOld         => null()
        real, dimension (:), pointer                :: MassCreated              => null()   !kg
        real, dimension (:), pointer                :: DNInflowConc             => null()
        real, dimension (:), pointer                :: MassInKg                 => null()   !kg (run with Benthos)
        real, dimension (:), pointer                :: TotalConc                => null()   !**WASSIM 16/11/2005 
        real, dimension (:), pointer                :: BottomConc               => null()   !kg m-2
        real, dimension (:), pointer                :: ErosionRate              => null()   !kg m-2 s-1
        real, dimension (:), pointer                :: DepositionRate           => null()   !kg m-3 s-1
        real, dimension (:), pointer                :: Ws                       => null()   !m s-1 (vertical velocity)
                                                                                            !positive direction is downswards                                                                                          
        real, dimension (:), pointer                :: OutputTime               => null()    !s
        real, dimension (:), pointer                :: InitialOutputTime        => null()    !s
        
        real                                        :: MinValue                 = null_real
        logical                                     :: WarnOnNegativeValues     = .false.
        real                                        :: InitialValue             = null_real
        real                                        :: BottomMinConc            = null_real  !kg m-2
        real                                        :: BoundaryConcentration    = null_real
        
        logical                                     :: Old                      = .false.
        
        !Advection Diffusion
        real                                        :: Diffusivity              = null_real
        !integer                                     :: Advection_Scheme         = null_int
        !integer                                     :: Diffusion_Scheme         = null_int
        
        !Decay
        real                                        :: DecayRate                = null_real    

        real                                        :: IScoefficient            = null_real
        real                                        :: ExtinctionCoefficient    = null_real
        real                                        :: ErosionCriticalShear     = null_real   
        real                                        :: DepositionCriticalShear  = null_real
        real                                        :: ErosionCoefficient       = null_real
        real                                        :: CHS                      = null_real
        integer                                     :: Ws_Type                  = null_int
        real                                        :: Ws_Value                 = null_real
        real                                        :: KL                       = null_real
        real                                        :: KL1                      = null_real
        real                                        :: ML                       = null_real
        real                                        :: M                        = null_real

        character(PathLength)                       :: OutputName               = null_str
        type (T_Property), pointer                  :: Next                     => null()
        type (T_Property), pointer                  :: Prev                     => null()
       
       !property dt in quality modules
        real                                        :: DTInterval               = null_real
        type(T_Time)                                :: LastCompute 
        type(T_Time)                                :: NextCompute        
    end type    T_Property        
    
    type T_Files
         character(len=Pathlength)                    :: ConstructData    = null_str 
         character(len=Pathlength)                    :: ReservoirFile    = null_str 
         character(len=Pathlength)                    :: InitialFile      = null_str 
         character(len=Pathlength)                    :: ResultsFile      = null_str 
         character(len=Pathlength)                    :: FinalFile        = null_str 
    end type T_Files    

    type T_TimeSerieImposed
        integer                                      :: ID                   = null_int 
        character(len=Pathlength)                    :: FileName             = null_str 
        integer                                      :: Column               = null_int 
    end  type T_TimeSerieImposed
    
    type T_FlowOver
        real                                         :: WeirLength           = null_real 
        real                                         :: DischargeCoeficient  = null_real 
        real                                         :: CrestLevel           = null_real 
    end  type T_FlowOver
    
    type T_Management
        logical                                       :: ON                   = .false.
        
        logical                                       :: HasDischarge         = .false.  
        
        !how to compute outflow
        logical                                       :: ImposedOutflow       = .false.    
        logical                                       :: ImposedLevel         = .false.    
        logical                                       :: ImposedOperation     = .false.    !operation rules        
        
        !Operation
        integer                                       :: OperationType        = null_int   
        real, dimension(:,:), pointer                 :: OperationCurve       => null()    
        integer                                       :: OperationCurvePoints = 0
        
        real                                          :: MinOutflow          = null_real  !Environmental flow
        real                                          :: MaxOutflow          = null_real  !Projected
        
        real, dimension(:,:), pointer                 :: AccVolumeCurve       => null()
        integer                                       :: AccVolumeCurvePoints = 0
        
    end type T_Management
    
    
    !Each reservoir properties
    type T_Reservoir
        integer                                       :: ID                   = null_int
        integer                                       :: Position             = null_int
        character(len = StringLength)                 :: Name                 = null_str
        real                                          :: CoordinateX          = null_real
        real                                          :: CoordinateY          = null_real

        integer                                       :: GridI                = null_int
        integer                                       :: GridJ                = null_int        
        
        !DN
        integer                                       :: DNNodeID             = null_int
        
        !State Variables
        real                                          :: VolumeOld            = null_real  !volume at start of timestep
        real                                          :: VolumeNew            = null_real  !volume conti being updated with fluxes
        real                                          :: VolumeTarget         = null_real  !volume final used if imposed level
        real                                          :: PercFull             = null_real
        
        !Derived State variables
        real                                          :: WaterLevel           = null_real
        
        !Fluxes
        real                                          :: DNInflow             = null_real
        real                                          :: Outflow              = null_real
        real                                          :: Discharges           = null_real
        real                                          :: SurfaceFluxes        = null_real    
		real                                          :: ExtraOutflow         = null_real
        
        !Management
        type(T_Management)                            :: Management
        type(T_FlowOver)                              :: FlowOver
        
        logical                                       :: IsWeir              = .false.
        
        !volumes
		integer                                       :: MaxVolumeType        = null_int    !Added by Ana Oliveira
        real                                          :: MinVolume            = null_real
        real                                          :: MaxVolume            = null_real
		real, dimension (12)                          :: MaxVolume_Monthly    = null_real   !Added by Ana Oliveira
        real                                          :: InitialVolume        = null_real
        logical                                       :: InitialVolumeDefined = .false.
        !Geometry
        real                                          :: WallHeight           = null_real
        real                                          :: SurfaceArea          = null_real
        real                                          :: BottomArea           = null_real
        
        !Other Info
        character(len = StringLength)                 :: WaterUse             = null_str
        character(len = StringLength)                 :: Owner                = null_str        
        
        logical                                       :: TimeSerie            = .false.
        character(len = StringLength)                 :: TimeSeriesName       = null_str
        
        type (T_Reservoir), pointer                   :: Next                 => null()
        type (T_Reservoir), pointer                   :: Prev                 => null()       
        type (T_TimeSerieImposed)                     :: TimeSeries           
        
        integer                                       :: ConstructionYear     = null_int    !Added by Ana Oliveira
    end type T_Reservoir    
    
    private :: T_Reservoirs
    type       T_Reservoirs
        integer                                     :: InstanceID
        character(len=StringLength)                 :: ModelName             = null_str
        type (T_Size1D)                             :: Size1D, WorkSize1D
        
        !Matrix for interaction with other modules
        real, dimension(:),  pointer                :: ReservoirVolumes      => null()
        real, dimension(:),  pointer                :: ReservoirPercFull     => null()
        real, dimension(:),  pointer                :: ReservoirInflows      => null()    !input from DN
        real, dimension(:),  pointer                :: ReservoirOutflows     => null()    !output to DN
        !real, dimension(:, :)),  pointer            :: ReservoirConc         => null()    !output to DN
        
        integer, dimension(:),  pointer             :: ReservoirsNodeIDs                  !DN Id's    
        
        integer, dimension(:),  pointer             :: ReservoirDischargeLink   => null()
        real, dimension(:),  pointer                :: ReservoirDischargeFlow   => null()
        real, dimension(:, :),  pointer             :: ReservoirDischargeConc   => null()
        logical, dimension(:),  pointer             :: DischargesActive         => null()
        
        !Initial condition
        integer                                     :: InitialVolumeDefaultMethod  = null_int
        real                                        :: StartPercentageFull   = null_real
        
        !Continuous
        logical                                     :: Continuous            = .false.
        
        type (T_Files)                              :: Files   
        type (T_Output)                             :: Output   
        type (T_ExtVar)                             :: ExtVar
        type (T_ComputeOptions)                     :: ComputeOptions
        type (T_TimeSerie)                          :: TimeSerie
        
       !Instance of Module_EnterData
        integer                                     :: ObjEnterData              = 0    !Data File - ConstructData        
        integer                                     :: ObjEnterDataReservoirFile = 0 
        integer                                     :: ObjDischarges             = 0
        integer                                     :: ObjTime                   = 0
        integer                                     :: ObjHDF5                   = 0
        
        type(T_Reservoirs), pointer                 :: Next
        type(T_Reservoir),  pointer                 :: FirstReservoir
        type(T_Reservoir),  pointer                 :: LastReservoir
        
        type(T_Property),   pointer                 :: FirstProperty
        type(T_Property),   pointer                 :: LastProperty
        
        logical                                     :: HasProperties         = .false.
        
        integer                                     :: nReservoirs           = 0
        integer                                     :: nProperties           = 0
        integer                                     :: nDischarges           = 0
        integer                                     :: nPropWithDischarges   = 0
    end type  T_Reservoirs

    !Global Module Variables
    type (T_Reservoirs), pointer                         :: FirstObjReservoirs
    type (T_Reservoirs), pointer                         :: Me

    !integer                                         :: mReservoirs_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructReservoirs(ModelName, ObjReservoirsID, TimeID, TopographicFile, STAT)

        !Arguments---------------------------------------------------------------
        character(len=*)                                :: ModelName
        integer                                         :: ObjReservoirsID 
        integer                                         :: TimeID 
        integer, optional, intent(OUT)                  :: STAT     
        character(len=Pathlength), optional             :: TopographicFile
        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mReservoirs_)) then
            nullify (FirstObjReservoirs)
            call RegisterModule (mReservoirs_) 
        endif

        call Ready(ObjReservoirsID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%ModelName        = ModelName
            
            !Associates module Time
            Me%ObjTime = AssociateInstance   (mTIME_, TimeID)            
            
            !Returns ID
            ObjReservoirsID          = Me%InstanceID

            if (present (TopographicFile)) then
                Me%ExtVar%TopographicFile = TopographicFile
                Me%Output%HDFActive = .true.
            endif            

            !Gets Current Compute Time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExTVar%ActualTime, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ConstructReservoirs - ModuleReservoirs - ERR00'
            
            !Gets Compute Time Limits
            call GetComputeTimeLimits (Me%ObjTime, BeginTime = Me%ExtVar%BeginTime,           &
                                       EndTime = Me%ExtVar%EndTime, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ConstructReservoirs - ModuleReservoirs - ERR01a'            
            
            
            call ReadFileNames
          
            
            !Open data file
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirs - ModuleReservoirs - ERR010' 
            
            call ConstructGlobalVariables            
                                 
            call ConstructReservoirList
            
            call ConstructPropertyList                        
            
            
            if (Me%ComputeOptions%Discharges) then
                call ConstructDischarges
            endif
            
            call VerifyOptions
            
            !Intial volume condition
            if (Me%Continuous) then
                call ReadInitialVolume
            else
                
                !Initialize all reservoirs volume
                call InitializeReservoirsVolume
            endif
            
            !Allocate and reset Variables
            call InitializeVariables(.true.)            
            
            if (Me%Output%HDF .or. Me%Output%TimeSerie) then
                call ConstructOutput
                !First Output
                call ModifyOutput
            endif
            
            !Close data file
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirs - ModuleReservoirs - ERR020'            
            
            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleReservoirs - ConstructReservoirs - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructReservoirs
 
    !--------------------------------------------------------------------------
    
    subroutine ReadFileNames
        
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------    
        integer                                          :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        !Read data file name
        call ReadFileName('RESERVOIRS_DAT', Me%Files%ConstructData, Message = "Reservoirs Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirs - ModuleReservoirs - ERR01b'      
    
        !Read the file name to place results
        call ReadFileName('RESERVOIRS_HDF', Me%Files%ResultsFile, Message = "Reservoirs HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirs - ModuleReservoirs - ERR01'        
            
        !Reads the name of the file where to write the final data
        call ReadFileName ('RESERVOIRS_FIN', Me%Files%FinalFile, 'Reservors Initial File', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialConc - ModuleReservoirs - ERR20'         
        
    
    end subroutine ReadFileNames
    
    !--------------------------------------------------------------------------    
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Reservoirs), pointer                         :: NewObjReservoirs
        type (T_Reservoirs), pointer                         :: PreviousObjReservoirs


        !Allocates new instance
        allocate (NewObjReservoirs)
        nullify  (NewObjReservoirs%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjReservoirs)) then
            FirstObjReservoirs         => NewObjReservoirs
            Me                    => NewObjReservoirs
        else
            PreviousObjReservoirs      => FirstObjReservoirs
            Me                    => FirstObjReservoirs%Next
            do while (associated(Me))
                PreviousObjReservoirs  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjReservoirs
            PreviousObjReservoirs%Next => NewObjReservoirs
        endif

        Me%InstanceID = RegisterNewInstance (mReservoirs_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------


    subroutine ConstructReservoirList

        !Local-----------------------------------------------------------------
        type (T_Reservoir), pointer                 :: NewReservoir      => null()
        integer                                     :: ClientNumber, numberReservoirs
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound
        character(LEN = StringLength)               :: block_begin          = '<beginreservoir>'
        character(LEN = StringLength)               :: block_end            = '<endreservoir>'        
        !----------------------------------------------------------------------   
        
        !Open data file
        call ConstructEnterData(Me%ObjEnterDataReservoirFile, Me%Files%ReservoirFile, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirList - ModuleReservoirs - ERR05'             
        
        !list with DN nodes where reservoirs are
        numberReservoirs = CountTotalReservoirs()
        allocate(Me%ReservoirsNodeIDs(numberReservoirs))
        
        ! Initialize the Reservoirs list
        Me%nReservoirs = 0
        nullify (Me%FirstReservoir)
        nullify (Me%LastReservoir)

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterDataReservoirFile, ClientNumber,          &
                                        block_begin, block_end, BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                
                    
                    ! Construct a New Reservoir 
                    call ConstructReservoir(NewReservoir, ClientNumber)

                    ! Add new Reservoir to the Reservoirs List 
                    call AddReservoir(NewReservoir)
                    
                else
                    
                    call Block_Unlock(Me%ObjEnterDataReservoirFile, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructReservoirList - ModuleReservoirs - ERR010'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructReservoirList - ModuleReservoirs - ERR020'
            end if cd1
        end do do1

        call KillEnterData(Me%ObjEnterDataReservoirFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoirList - ModuleReservoirs - ERR030'        
        

    end subroutine ConstructReservoirList

    !--------------------------------------------------------------------------    
    
    real function CountTotalReservoirs()

        !This subroutine counts the total number of reservoirs and checks the
        !existence of valid and repeated reservoirsIDs
        !Local------------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine, Reservoirs           
        integer                                     :: STAT_CALL
        integer                                     :: ReservoirID, OldReservoirID
        integer                                     :: flag
        character(LEN = StringLength)               :: block_begin          = '<beginreservoir>'
        character(LEN = StringLength)               :: block_end            = '<endreservoir>' 
        
        Reservoirs = 0
        OldReservoirID = null_int

do1:    do 

            call ExtractBlockFromBuffer(Me%ObjEnterDataReservoirFile, ClientNumber,     &
                                        block_begin, block_end, BlockFound,                 &
                                        FirstLine, LastLine, STAT_CALL) 

if1:        if (STAT_CALL .EQ. SUCCESS_) then    

if2:            if (BlockFound) then                 
                    
                    !Gets ID
                    call GetData(ReservoirID,                                   &
                                 Me%ObjEnterDataReservoirFile, flag,      & 
                                 keyword      = 'ID',                           &
                                 ClientModule = 'ModuleReservoirs',             &
                                 SearchType   = FromBlock,                      &
                                 STAT         = STAT_CALL)                                  
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - CountTotalReservoirs - ERR01'

                    if (flag /= 1) then
                        write (*,*)'Invalid Reservoir ID [ID]'
                        stop 'ModuleReservoirs - CountTotalReservoirs - ERR02'
                    endif                                  
                   
                    if (ReservoirID .EQ. OldReservoirID ) then
                        write (*,*) 'Repeated Reservoir ID = ', ReservoirID
                        stop 'ModuleReservoirs - CountTotalReservoirs - ERR03'
                    else
                        OldReservoirID = ReservoirID
                    end if
                    
                    Reservoirs = Reservoirs + 1  
                    
                else if2                    

                    call Block_Unlock(Me%ObjEnterDataReservoirFile, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - CountTotalReservoirs - ERR01'

                    call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - CountTotalReservoirs - ERR02'
                    
                    exit do1    !No more blocks

                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1

                stop 'ModuleReservoirs - CountTotalReservoirs - ERR03.'

            end if if1

        end do do1

        CountTotalReservoirs = Reservoirs
        
    end function CountTotalReservoirs
       
    !---------------------------------------------------------------------------    
    
    !This subroutine reads all the information needed to construct a new property.           
    subroutine ConstructReservoir(NewReservoir, ClientID)

        !Arguments-------------------------------------------------------------
        type(T_Reservoir), pointer                  :: NewReservoir
        integer                                     :: ClientID

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        real, dimension (2)                         :: AuxCoord 
        !----------------------------------------------------------------------
             

        !Allocates new property
        allocate (NewReservoir, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProperty - ModuleReservoirs - ERR01'

        nullify(NewReservoir%Next     )
        nullify(NewReservoir%Prev     )
        
        call GetData(NewReservoir%ID,                                                   &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='ID',                                                &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR010'        

        if (iflag /= 1) then
            write (*,*)'Invalid Node ID [ID]'
            stop 'ConstructReservoir - ModuleReservoirs - ERR011'
        endif        
        
        call GetData(NewReservoir%Name,                                                 &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='NAME',                                              &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR020'
        if (iflag == 0) then
            write(*,*) 'Not Found Reservoir Name [NAME]'
            stop 'ConstructReservoirValues - ModuleReservoirs - ERR021'
        endif        
        
        call GetData(NewReservoir%DNNodeID,                                             &
                     Me%ObjEnterDataReservoirFile, iflag,                               &
                     SearchType   = FromBlock,                                          &
                     keyword      ='DN_NODE_ID',                                        &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR025'
        if (iflag == 0) then
            write(*,*) 'Not Found Reservoir Drainage Network Node ID where reservoir is located [DN_NODE_ID]'
            stop 'ConstructReservoirValues - ModuleReservoirs - ERR026'
        endif
        
        !add to global matrix
        Me%ReservoirsNodeIDs(Me%nReservoirs + 1) = NewReservoir%DNNodeID
        
        call GetData(AuxCoord,                                                          &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='COORDINATES' ,                                      &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR030'
        if (iflag == 2) then
            NewReservoir%CoordinateX = AuxCoord(1)
            NewReservoir%CoordinateY = AuxCoord(2)                       
        else
            write(*,*) 'Invalid Reservoir Coordenates [COORDINATES]'
            stop 'ConstructReservoir - ModuleReservoirs - ERR31'
        end if 
            
        call GetData(NewReservoir%GridI,                                                &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='GRID_I',                                            &
                     default      = null_int,                                           &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR040'

        call GetData(NewReservoir%GridJ,                                                &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='GRID_J',                                            &
                     default      = null_int,                                           &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR050'   
        
        !not accounting dead storage
        call GetData(NewReservoir%MinVolume,                                            &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='MIN_VOLUME',                                        &
                     default      = 0.0,                                                &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR060'         
		
		!Choose if the maximum value of reservoir varies every month or if it is constant - by Ana Oliveira
        call GetData(NewReservoir%MaxVolumeType,                                        &
                     Me%ObjEnterDataReservoirFile, iflag,                               &
                     SearchType     = FromBlock,                                        &
                     Keyword        = 'MAX_VOLUME_TYPE',                                &
                     default        = 1,                                                &
                     ClientModule   = 'ModuleReservoirs',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR70'  
        
        !Total Capacity (at full supply level)
        if (NewReservoir%MaxVolumeType == UniqueValue_) then
            call GetData(NewReservoir%MaxVolume,                                            &
                         Me%ObjEnterDataReservoirFile, iflag,                                            &
                         SearchType   = FromBlock,                                          &
                         keyword      ='MAX_VOLUME',                                        &
                         ClientModule = 'ModuleReservoirs',                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR070.1'
        endif
        
        if (NewReservoir%MaxVolumeType == MonthlyValues_) then
            call GetData(NewReservoir%MaxVolume_Monthly,                                    &
                         Me%ObjEnterDataReservoirFile, iflag,                               &
                         SearchType   = FromBlock,                                          &
                         keyword      ='MAX_VOLUME',                                        &
                         ClientModule = 'ModuleReservoirs',                                 &
                         STAT         = STAT_CALL)
            NewReservoir%MaxVolume = maxval(NewReservoir%MaxVolume_Monthly)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR070.2'
        endif
        if (iflag == 0) then
            write(*,*) 'Not Found Reservoir maximum volume [MAX_VOLUME]'
            stop 'ConstructReservoirValues - ModuleReservoirs - ERR061'
        endif        
        
        !Initial value
        call GetData(NewReservoir%InitialVolume,                                        &
                     Me%ObjEnterDataReservoirFile, iflag,                               &
                     SearchType   = FromBlock,                                          &
                     keyword      ='INITIAL_VOLUME',                                    &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR065'       
        if (iflag /= 0 .and. NewReservoir%InitialVolume >= 0.0 .and. NewReservoir%InitialVolume <= NewReservoir%MaxVolume) then
            NewReservoir%InitialVolumeDefined = .true.            
        endif
        
        call GetData(NewReservoir%WallHeight,                                           &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='WALL_HEIGHT',                                       &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR070'            

        call GetData(NewReservoir%SurfaceArea,                                          &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='SURFACE_AREA',                                      &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR080'           
        
        if (Me%ComputeOptions%BottomFluxes) then
            
            call GetData(NewReservoir%BottomArea,                                          &
                         Me%ObjEnterDataReservoirFile, iflag,                                            &
                         SearchType   = FromBlock,                                          &
                         keyword      ='BOTTOM_AREA',                                       &
                         ClientModule = 'ModuleReservoirs',                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR085'  
            if (iflag == 0) then
                write(*,*) 'Not Found Reservoir bottom area [BOTTOM_AREA]'
                write(*,*) 'It is mandatory when computing bottom fluxes'
                stop 'ConstructReservoirValues - ModuleReservoirs - ERR086'                
            endif
        endif
        
        call GetData(NewReservoir%WaterUSe,                                             &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='WATER_USE',                                         &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR090'           

        call GetData(NewReservoir%Owner,                                                &
                     Me%ObjEnterDataReservoirFile, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      ='OWNER',                                             &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR100'           

        !reset
        NewReservoir%Management%ON                    = .false.
        NewReservoir%Management%ImposedOperation      = .false.
        NewReservoir%Management%ImposedLevel          = .false.
        NewReservoir%Management%ImposedOutflow        = .false.               
        
        
        !Get curve with acc volume for level
        !Always fetched is used for transforming volume to level
        !In case of weir, imposed level or operations based on level is mandatory
        call GetVolumeAccCurve(NewReservoir, ClientID)        
        
        
        !Get if reservoir has imposed level
        call GetImposedReservoirLevel(NewReservoir, ClientID, .true.)        
        
        !Get management options (operation curves, environmental flow)
        !or unmanaged weir
        call GetManagementOptions(NewReservoir, ClientID)
        
      
        
        !max ouflow trough all discharges (projected). if not defined almost infinite
        call GetData(NewReservoir%Management%MaxOutflow,                                &
                     Me%ObjEnterDataReservoirFile, iflag,                               &
                     SearchType   = FromBlock,                                          &
                     keyword      ='MAX_OUTFLOW',                                       &
                     default      = -null_real,                                         &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR130'         
        
        !In case that is forcing with operation pec max outflow this property is mandatory and not dummy
        if (NewReservoir%Management%OperationType == Operation_PercVol_PercMaxOutflow_) then
            if (NewReservoir%Management%MaxOutflow .gt. -null_real / 2.0) then
                write(*,*) 'Using operation curve of type percentage max outflow and'
                write(*,*) 'MAX_OUTFLOW is not defined '
                write(*,*) 'in reservoir ID : ', NewReservoir%ID                
                stop 'ConstructReservoir - ModuleReservoirs - ERR0140'                  
            endif
        endif
        
        ! Added by Ana Oliveira
        call GetData(NewReservoir%ConstructionYear,                                     &
                     Me%ObjEnterDataReservoirFile, iflag,                               &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'CONSTRUCTION_YEAR',                                &
                     default      = 0,                                                  &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR150'
        

    end subroutine ConstructReservoir

    !--------------------------------------------------------------------------
    
    subroutine GetManagementOptions(NewReservoir, ClientID)
    
        !Arguments-------------------------------------------------------------
        type(T_Reservoir), pointer                  :: NewReservoir
        integer                                     :: ClientID

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine, NLayers, l, line
        real, dimension(:), allocatable             :: Aux 
        !----------------------------------------------------------------------
        
        call RewindBlock(Me%ObjEnterDataReservoirFile, ClientId, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'GetManagementOptions - ModuleReservoirs - ERR01'        
        
        !Get Management Options. Start from more complex managed to unmanaged        
        call GetData(NewReservoir%Management%OperationType,                             &
                     Me%ObjEnterDataReservoirFile, iflag,                               &
                     SearchType   = FromBlock,                                          &
                     keyword      ='OPERATION_TYPE',                                    &
                     default      = null_int,                                           &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR110'         

        if (NewReservoir%Management%OperationType /= null_int) then
            
            if (NewReservoir%Management%OperationType /= Operation_Level_Outflow_ .and.  &
                NewReservoir%Management%OperationType /= Operation_Level_PercInflow_ .and.  &
                NewReservoir%Management%OperationType /= Operation_PercVol_Outflow_ .and.  &
                NewReservoir%Management%OperationType /= Operation_PercVol_PercInflow_ .and. &
                NewReservoir%Management%OperationType /= Operation_PercVol_PercMaxOutflow_) then
                write(*,*) 'Unknown OPERATION_TYPE'
                write(*,*) 'in reservoir ID : ', NewReservoir%ID
                stop 'ConstructReservoir - ModuleReservoirs - ERR0110a'                     
            endif
                      
            call ExtractBlockFromBlock(Me%ObjEnterDataReservoirFile, ClientID,     &
                                        '<<beginoperation>>', '<<endoperation>>', BlockFound,     &
                                        FirstLine = FirstLine, LastLine = LastLine,               &            
                                        STAT = STAT_CALL)

            if(STAT_CALL .EQ. SUCCESS_ .and. BlockFound) then    
                                                
                NLayers =  LastLine - FirstLine - 1
                
                if (NLayers > 0) then
                    
                    NewReservoir%Management%ON                    = .true.
                    NewReservoir%Management%ImposedOperation      = .true.
                    NewReservoir%Management%OperationCurvePoints  = NLayers
                    
                    allocate (NewReservoir%Management%OperationCurve(NLayers, 2))
 
                    !Allocates auxiliar variables
                    allocate (Aux(2))
            
                    l = 1
                    do line = FirstLine + 1, LastLine - 1

                        call GetData(Aux, EnterDataID = Me%ObjEnterDataReservoirFile, flag = iflag, &
                                        SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoir - ModuleReservoirs - ERR111'

                        NewReservoir%Management%OperationCurve(l, 1) = Aux(1)
                        NewReservoir%Management%OperationCurve(l, 2) = Aux(2)
                        
                        !cant have equal or lower level, or perentage volumes. it is ascending order
                        if (l > 1) then
                            if (NewReservoir%Management%OperationCurve(l, 1) <=     &
                                 NewReservoir%Management%OperationCurve(l - 1, 1)) then
                                write(*,*) 'operation curve needs to be in ascending order of level or '
                                write(*,*) 'percentage volume (1st column)'
                                write(*,*) 'in reservoir ID : ', NewReservoir%ID
                                stop 'ConstructReservoir - ModuleReservoirs - ERR0112'  
                            endif
                        endif
                        
                        if (NewReservoir%Management%OperationType /= Operation_Level_PercInflow_ .or.    &
                            NewReservoir%Management%OperationType /= Operation_PercVol_PercInflow_) then
                            if (NewReservoir%Management%OperationCurve(l, 2) < 0.0 .or.                  &
                                NewReservoir%Management%OperationCurve(l, 2) > 1.0) then
                                write(*,*) 'operation curve for otflow from percentage of inflow cant have values < 0 or > 1.0'
                                write(*,*) 'in reservoir ID : ', NewReservoir%ID
                                stop 'ConstructReservoir - ModuleReservoirs - ERR0113'                                
                            endif
                        endif
                        
                        l = l + 1

                    enddo

                    deallocate(Aux)
                    
                endif
            end if 

        endif
        
        !if no operation found, search for minimum flow
        if (.not. NewReservoir%Management%ON) then
            
            !Environmental flow. if not defined zero
            call GetData(NewReservoir%Management%MinOutflow,                                &
                         Me%ObjEnterDataReservoirFile, iflag,                               &
                         SearchType   = FromBlock,                                          &
                         keyword      ='MIN_OUTFLOW',                                       &
                         default      = 0.0,                                                &
                         ClientModule = 'ModuleReservoirs',                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR113'  
            
            !simple amanagement
            if (NewReservoir%Management%MinOutflow > 0.0) then
                NewReservoir%Management%ON                    = .true.             
            
            
            else !no management - type ditch or weir                          
                
                call GetData(NewReservoir%IsWeir,                                               &
                             Me%ObjEnterDataReservoirFile, iflag,                               &
                             SearchType   = FromBlock,                                          &
                             keyword      ='IS_WEIR',                                           &
                             default      = .false.,                                            &
                             ClientModule = 'ModuleReservoirs',                                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ConstructReservoir - ModuleReservoirs - ERR114'                  
                
                if (NewReservoir%IsWeir) then                    
                    
                  call GetData(NewReservoir%FlowOver%WeirLength,                              &
                                 Me%ObjEnterDataReservoirFile,                                  &
                                 iflag,                                                         &
                                 FromBlock,                                                     &
                                 keyword      ='WEIR_LENGTH',                                   &
                                 ClientModule = 'ModuleDischarges',                             &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoir - ModuleReservoirs - ERR115'
                    if (iflag /= 1) then
                        write(*,*) 'Weir Length Missing'
                        stop ' Construct_FlowValues - ModuleDischarges - ERR60'
                    endif

                    call GetData(NewReservoir%FlowOver%DischargeCoeficient,                     &
                                 Me%ObjEnterDataReservoirFile,                                  &
                                 iflag,                                                         &
                                 FromBlock,                                                     &
                                 keyword      ='WEIR_COEF',                                     &
                                 ClientModule = 'ModuleDischarges',                             &
                                 default      = 0.4,                                            &
                                STAT         = STAT_CALL)
            
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoir - ModuleReservoirs - ERR116'

                    call GetData(NewReservoir%FlowOver%CrestLevel,                              &
                                 Me%ObjEnterDataReservoirFile,                                  &
                                 iflag,                                                         &
                                 FromBlock,                                                     &
                                 keyword      ='CREST_LEVEL',                                   &
                                 ClientModule = 'ModuleDischarges',                             &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoir - ModuleReservoirs - ERR117'
            
                    if (iflag /= 1) then
                        write(*,*) 'Crest Height Missing'
                        stop ' ConstructReservoir - ModuleReservoirs - ERR90'
                    endif                  
                endif
            endif
        
        endif
        
        
        
    
    end subroutine GetManagementOptions
    
    !--------------------------------------------------------------------------
    
    subroutine GetVolumeAccCurve(NewReservoir, ClientID)
    
        !Arguments-------------------------------------------------------------
        type(T_Reservoir), pointer                  :: NewReservoir
        integer                                     :: ClientID

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine, NLayers, l, line
        real, dimension(:), allocatable             :: Aux 
        !----------------------------------------------------------------------
    
        call RewindBlock(Me%ObjEnterDataReservoirFile, ClientId, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'GetVolumeAccCurve - ModuleReservoirs - ERR01'        
        
        
        !Need volume - level curve
        !if (NewReservoir%IsWeir .or. NewReservoir%Management%OperationType == Operation_Level_Outflow_                      &
        !    .or. NewReservoir%Management%OperationType == Operation_Level_PercInflow_) then
        !Always get
            call ExtractBlockFromBlock(Me%ObjEnterDataReservoirFile, ClientID,               &
                                        '<<beginaccvolumecurve>>', '<<endaccvolumecurve>>', BlockFound,     &        
                                            FirstLine = FirstLine, LastLine = LastLine,                        & 
                                        STAT = STAT_CALL)

            if(STAT_CALL .EQ. SUCCESS_ .and. BlockFound) then    
                                                
                NLayers =  LastLine - FirstLine - 1
                
                if (NLayers > 0) then
                    
                    NewReservoir%Management%AccVolumeCurvePoints        = NLayers
                    
                    allocate (NewReservoir%Management%AccVolumeCurve(NLayers, 2))
 
                    !Allocates auxiliar variables
                    allocate (Aux(2))
            
                    l = 1
                    do line = FirstLine + 1, LastLine - 1

                        call GetData(Aux, EnterDataID = Me%ObjEnterDataReservoirFile, flag = iflag, &
                                        SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructReservoir - ModuleReservoirs - ERR118'

                        NewReservoir%Management%AccVolumeCurve(l, 1) = Aux(1)
                        NewReservoir%Management%AccVolumeCurve(l, 2) = Aux(2)
                        l = l + 1

                    enddo
                        
                    !Verify if volumes are comprised in limits
                    if (NewReservoir%Management%AccVolumeCurve(1,1) > NewReservoir%MinVolume) then
                        write(*,*) 'First line of accumulated volume curve should be at least'
                        write(*,*) 'lower or equal the reservoir minimum volume in reservoir ID : ', NewReservoir%ID
                        stop 'ConstructReservoir - ModuleReservoirs - ERR0119'                                
                    endif
                    if (NewReservoir%Management%AccVolumeCurve(NLayers, 1) < NewReservoir%MaxVolume) then
                        write(*,*) 'Last line of accumulated volume curve should be at least'
                        write(*,*) 'higher or equal the reservoir maximum volume in reservoir ID : ', NewReservoir%ID
                        stop 'ConstructReservoir - ModuleReservoirs - ERR0120'                                
                    endif    
                        
                    deallocate(Aux)
                    
                endif
            !endif
                             
            end if 
        !endif        
        
        
    end subroutine GetVolumeAccCurve    
    
    !---------------------------------------------------------------------------
    
   subroutine ConstructDischarges

        !Arguments--------------------------------------------------------------
                                                    
        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iDis, ReservoirID
        logical                                         :: Found
        type (T_Reservoir), pointer                     :: Reservoir
        logical                                         :: IsImposedOutflow


        call Construct_Discharges(Me%ObjDischarges,                              &
                                  Me%ObjTime,                                    &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleReservoirs - ConstructDischarges - ERR02'  
        
        !Build Discharge ReservoirID link
        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, Me%nDischarges, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleReservoirs - ConstructDischarges - ERR03'

        allocate(Me%ReservoirDischargeLink(Me%nDischarges))
        allocate(Me%ReservoirDischargeFlow(Me%nDischarges))
        allocate(Me%ReservoirDischargeConc(Me%nDischarges, Me%nPropWithDischarges))
        allocate(Me%DischargesActive(Me%nDischarges))
        
        do iDis = 1, Me%nDischarges

            call GetDischargesReservoirID  (Me%ObjDischarges, iDis, ReservoirID, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModuleReservoirs - ConstructDischarges - ERR04'

            !ignore discharges from DN (that can cohesist)
            if (ReservoirID > 0) then
                
                Me%DischargesActive(iDis) = .true.            
            
                call FindReservoir   (ReservoirID, Reservoir, Found)
                            
                if (Found) then
                    Reservoir%Management%HasDischarge = .true.                
                    Me%ReservoirDischargeLink(iDis) = Reservoir%ID
                else
                    write (*,*) 'Discharge Reservoir not found'
                    write (*,*) 'Reservoir ID = ', ReservoirID            
                    stop 'ModuleReservoirs - ConstructDischarges - ERR05'
                end if
            
                !has imposed outflow discharge?
                call GetIsReservoirOutflow(Me%ObjDischarges, iDis, IsImposedOutflow, STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_) stop 'ModuleReservoir - ConstructDischarges - ERR06'         
            
                !if any of the discharges are outflow, mark it (the reservoir will not compute outflow, it will be imposed)
                !Imposed outflow imposes itself to the other. Sometimes when there is data the reservoir can have
                !imposed flow and when removed the latter conditons will prevail
                if (IsImposedOutflow) then
                    Reservoir%Management%ON                 = .true.
                    Reservoir%Management%ImposedOutflow     = .true.
                    

                endif
            else
                Me%DischargesActive(iDis) = .false.  
            endif
        end do
    
    
    end subroutine ConstructDischarges            
    
    !--------------------------------------------------------------------------    
    
    subroutine FindReservoir (ReservoirID, Reservoir, Found)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                             :: ReservoirID
        type (T_Reservoir), pointer, intent(OUT)        :: Reservoir
        logical, intent(OUT)                            :: Found
        !Local------------------------------------------------------------------


        Found = .FALSE.
        
        nullify(Reservoir)
        Reservoir => Me%FirstReservoir
        
        do while (associated(Reservoir))
        
            if (Reservoir%ID == ReservoirID) then
                Found = .TRUE.
                exit
            end if
            Reservoir => Reservoir%Next
        end do

    end subroutine FindReservoir

    !---------------------------------------------------------------------------   
   
   subroutine ConstructPropertyList

        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: NewProperty
        integer                                     :: ClientNumber
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound
        character(LEN = StringLength)               :: block_begin          = '<beginproperty>'
        character(LEN = StringLength)               :: block_end            = '<endproperty>' 

        ! Initialize the properties number   
        Me%nProperties = 0

        ! Initialize the properties list   
        nullify (Me%FirstProperty)        

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                        block_begin, block_end, BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                                        
                    call ConstructProperty(NewProperty)
                    
                    call AddProperty(NewProperty)

                else
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructPropertyList - ModuleReservoirs - ERR01'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructPropertyList - ModuleReservoirs - ERR02'
            end if cd1
        end do do1

        if (Me%nProperties .GT. 0) then
            Me%HasProperties = .true.
        end if

    end subroutine ConstructPropertyList

    !---------------------------------------------------------------------------
        
    subroutine ConstructProperty(NewProperty)

        !Arguments--------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External---------------------------------------------------------------
        integer                         :: STAT_CALL

        !-----------------------------------------------------------------------
             

        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProperty - ModuleReservoirs - ERR01'

        allocate (NewProperty%Concentration            (1:Me%nReservoirs))
        allocate (NewProperty%ConcentrationOld         (1:Me%nReservoirs))
        allocate (NewProperty%MassCreated              (1:Me%nReservoirs))
        allocate (NewProperty%DNInflowConc             (1:Me%nReservoirs))
        allocate (NewProperty%TotalConc                (1:Me%nReservoirs))
        allocate (NewProperty%MassInKg                 (1:Me%nReservoirs))
        allocate (NewProperty%OutputTime               (1:Me%nReservoirs))
        allocate (NewProperty%InitialOutputTime        (1:Me%nReservoirs))            
       
        NewProperty%Concentration           = 0.0
        NewProperty%ConcentrationOld        = 0.0
        NewProperty%MassCreated             = 0.0
        NewProperty%DNInflowConc            = 0.0
        NewProperty%TotalConc               = 0.0
        NewProperty%MassInKg                = 0.0
        NewProperty%InitialOutputTime       = 0.0
        NewProperty%OutputTime              = 0.0

        call ConstructPropertyID     (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call ConstructPropertyValues (NewProperty)
        

    end subroutine ConstructProperty

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyValues (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),   pointer                 :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        real                                        :: BottomInitialConc
        real                                        :: ModelDT, auxFactor, errorAux, DTaux
        integer                                     :: ObjHorizontalGrid = 0
        integer                                     :: ObjBasinGeometry  = 0      
        integer                                     :: ObjGridData = 0
        integer, dimension(:,:), pointer            :: BasinPoints       => null()
        real, dimension(:,:), pointer               :: ReservoirConc
        type (T_Reservoir), pointer                 :: CurrReservoir
        !Begin-----------------------------------------------------------------
        
         !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModuleReservoirs - ERR120'
          
        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized (user can provide a grid data or constant value)
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not. NewProperty%Old) then

            call ConstructHorizontalGrid(ObjHorizontalGrid, Me%ExtVar%TopographicFile, &
                                            STAT = STAT_CALL)           
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR010'
            
            !Constructs GridData
            call ConstructGridData      (ObjGridData, ObjHorizontalGrid,           &
                                            FileName = Me%ExtVar%TopographicFile,  &
                                            STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR020'

            !Constructs BasinGeometry
            call ConstructBasinGeometry (ObjBasinGeometry, ObjGridData,            &
                                            ObjHorizontalGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR030'
            
            
            !Gets BasinPoints
            call GetBasinPoints         (ObjBasinGeometry, BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR040'     

            !Gets the size of the grid
            call GetHorizontalGridSize (ObjHorizontalGrid,                               &
                                        Size     = Me%ExtVar%Size2D,                     &
                                        WorkSize = Me%ExtVar%WorkSize2D,                 &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR050'               
            
            allocate(ReservoirConc           (Me%ExtVar%Size2D%ILB:Me%ExtVar%Size2D%IUB,   &
                                              Me%ExtVar%Size2D%JLB:Me%ExtVar%Size2D%JUB))
            call SetMatrixValue (ReservoirConc,          Me%ExtVar%Size2D, 0.0)            
            
            
            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = ObjHorizontalGrid,                &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = BasinPoints,                      &
                                       Matrix2D             = ReservoirConc,                    &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'ConstructPropertyValues - ModuleReservoirs - ERR140'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'ConstructPropertyValues - ModuleReservoirs - ERR0150'
                
            else
                write(*,*)
                write(*,*)'Cant impose reservoir properties in time'
                stop 'ConstructPropertyValues - ModuleReservoirs - ERR0151'
            end if

            
            CurrReservoir => Me%FirstReservoir
            do while (associated(CurrReservoir))
                
                NewProperty%Concentration(CurrReservoir%Position) = ReservoirConc(CurrReservoir%GridI, CurrReservoir%GridJ)
                                                
                CurrReservoir => CurrReservoir%Next
            enddo                
            
            deallocate (ReservoirConc)
            
            !Ungets BasinPoints
            call UngetBasin             (ObjBasinGeometry, BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR152' 
            
            call KillBasinGeometry      (ObjBasinGeometry,   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR160'

            call KillGridData           (ObjGridData,        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR170'            
            
            call KillHorizontalGrid     (ObjHorizontalGrid,  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR180'             

        else

            if (Me%Continuous) then
                call ConstructHorizontalGrid(ObjHorizontalGrid, Me%ExtVar%TopographicFile, &
                                                STAT = STAT_CALL)           
                if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR010'
                ! If the property is old then the program is going to try to find a property
                ! with the same name in the Water properties initial file written in HDF format
                call GetHorizontalGridSize (ObjHorizontalGrid,                               &
                                            Size     = Me%ExtVar%Size2D,                     &
                                            WorkSize = Me%ExtVar%WorkSize2D,                 &
                                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR050'
                
                call ReadInitialConc(NewProperty)
                
                call KillHorizontalGrid     (ObjHorizontalGrid,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleReservoirs - ERR180' 
            else
                write(*,*)
                write(*,*)'Cant use OLD keyword on properties if CONTINUOUS is OFF'
                write(*,*)'Property : ', trim(NewProperty%ID%Name)
                stop 'ConstructPropertyValues - ModuleReservoirs - ERR190' 
            endif

        end if  


        call GetData(NewProperty%MinValue,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'MIN_VALUE',                              &
                     ClientModule   = 'ModuleReservoirs',                  &
                     SearchType     = FromBlock,                                &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR03'
        
        if (iflag==1)  then
            NewProperty%ComputeOptions%MinConcentration = .true.
            Me%ComputeOptions%MinConcentration          = .true.
        else
            NewProperty%ComputeOptions%MinConcentration = .false.
        endif

        call GetData(NewProperty%ComputeOptions%WarnOnNegativeValues,                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType     = FromBlock,                                         &
                     keyword        = 'WARN_ON_NEGATIVE_VALUES',                         &
                     ClientModule   = 'ModuleReservoirs',                           &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR3.5' 
        
        if (NewProperty%ComputeOptions%WarnOnNegativeValues) Me%ComputeOptions%WarnOnNegativeValues = .true.
        
!        call GetData(NewProperty%ComputeOptions%AdvectionDiffusion,                 &
!                     Me%ObjEnterData, iflag,                                        &
!                     Keyword        = 'ADVECTION_DIFUSION',                         &
!                     ClientModule   = 'ModuleReservoirs',                      &
!                     SearchType     = FromBlock,                                    &
!                     Default        = ON,                                           &
!                     STAT           = STAT_CALL)              
!        if (STAT_CALL .NE. SUCCESS_)                                                &
!            stop 'ModuleReservoirs - ConstructPropertyValues - ERR04'
!
!if1:    if (NewProperty%ComputeOptions%AdvectionDiffusion) then
!            
!            Me%ComputeOptions%AdvectionDiffusion = .true.
!                
!            !Numerical Discretization of Advection
!            call GetData(NewProperty%Advection_Scheme,                             &
!                         Me%ObjEnterData, iflag,                                   &
!                         Keyword        = 'ADVECTION_SCHEME',                      &
!                         ClientModule   = 'ModuleReservoirs',                 &
!                         SearchType     = FromBlock,                               &
!                         Default        = UpwindOrder1,                            &
!                         STAT           = STAT_CALL)              
!            if (STAT_CALL .NE. SUCCESS_)                                           &
!                stop 'ModuleReservoirs - ConstructPropertyValues - ERR05'
!
!            if (NewProperty%Advection_Scheme /= UpwindOrder1) then
!                write (*,*) 'Invalid option for keyword [ADVECTION_SCHEME]'
!                stop 'ConstructPropertyValues - ModuleReservoirs - ERR06'
!            end if
!!            if (NewProperty%Advection_Scheme /= UpwindOrder1 .AND.                 &
!!                NewProperty%Advection_Scheme /= CentralDif) then
!!                write (*,*) 'Invalid keyword [ADVECTION_SCHEME]'
!!                stop 'ConstructPropertyValues - ModuleReservoirs - ERR06'
!!            end if
!
!
!             call GetData(NewProperty%Diffusion_Scheme,                             &
!                          Me%ObjEnterData, iflag,                                   &
!                          Keyword        = 'DIFFUSION_SCHEME',                      &
!                          ClientModule   = 'ModuleReservoirs',                 &
!                          SearchType     = FromBlock,                               &
!                          Default        = CentralDif,                              &
!                          STAT           = STAT_CALL)              
!            if (STAT_CALL .NE. SUCCESS_)                                            &
!                stop 'ModuleReservoirs - ConstructPropertyValues - ERR07'
! 
!             
!            if (NewProperty%Diffusion_Scheme /= CentralDif) then
!                write (*,*) 'Invalid keyword [DIFFUSION_SCHEME]'
!                stop 'ConstructPropertyValues - ModuleReservoirs - ERR08'
!            end if
! 
!            !Molecular diffusivity of property in m2/s
!            call GetData(NewProperty%Diffusivity,                                   &
!                         Me%ObjEnterData, iflag,                                    &
!                         Keyword        = 'DIFFUSIVITY',                            &
!                         ClientModule   = 'ModuleReservoirs',                  &
!                         SearchType     = FromBlock,                                &
!                         Default        = 1e-8,                                     &
!                         STAT           = STAT_CALL)              
!            if (STAT_CALL .NE. SUCCESS_)                                            &
!                stop 'ModuleReservoirs - ConstructPropertyValues - ERR09' 
!        
!        
!        end if if1

        call GetData(NewProperty%ComputeOptions%Discharges,                         &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'DISCHARGES',                                 &
                     ClientModule   = 'ModuleReservoirs',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR10' 

        if (NewProperty%ComputeOptions%Discharges .and. .not. Me%ComputeOptions%Discharges) then
            call SetError(WARNING_, INTERNAL_, 'Missing keyword [DISCHARGES]', ON)            
            Me%ComputeOptions%Discharges = ON
        end if
        
        if (NewProperty%ComputeOptions%Discharges) Me%nPropWithDischarges = Me%nPropWithDischarges + 1

        
        !Dacay Time
        call GetData(NewProperty%ComputeOptions%T90_Decay,                          &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'DECAY_T90',                                  &
                     ClientModule   = 'ModuleReservoirs',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR11' 
        
        if (NewProperty%ComputeOptions%T90_Decay .and. .not. NewProperty%ComputeOptions%Discharges) then
            write (*,*) 'Decaying properties must be discharged',                   &
                         trim(adjustl(adjustr(NewProperty%ID%Name)))
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR12' 
        end if
         
        if (NewProperty%ComputeOptions%T90_Decay ) then
            Me%ComputeOptions%T90_Decay = ON
        endif

        !Generic decay
        call GetData(NewProperty%ComputeOptions%Generic_Decay,                           &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DECAY_GENERIC',                                     &
                     ClientModule = 'ModuleReservoirs',                             &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ConstructPropertyValues - ModuleReservoirs - ERR50'
        if (NewProperty%ComputeOptions%Generic_Decay) then
            Me%ComputeOptions%Generic_Decay       = .true.
        endif

        if (NewProperty%ComputeOptions%Generic_Decay) then
            
            !Decay rate k (s-1) in P = Po*exp(-kt)
            call GetData(NewProperty%DecayRate,                                              &
                         Me%ObjEnterData,iflag,                                              &
                         SearchType   = FromBlock,                                           &
                         keyword      = 'DECAY_RATE',                                        &
                         ClientModule = 'ModuleReservoirs',                             &
                         default      = 0.,                                                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                     &
                stop 'ConstructPropertyValues - ModuleReservoirs - ERR70'            
            
        endif

        !Checks for Surface Fluxes
        call GetData(NewProperty%ComputeOptions%SurfaceFluxes,                      &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'SURFACE_FLUXES',                               &
                     Default      = .false.,                                        & 
                     ClientModule = 'ModuleReservoirs',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR80' 

        if (NewProperty%ComputeOptions%SurfaceFluxes .and. .not. Me%ComputeOptions%SurfaceFluxes) then
            call SetError(WARNING_, INTERNAL_, 'Missing keyword [SURFACE_FLUXES]', ON)            
            Me%ComputeOptions%SurfaceFluxes = ON
        end if
        
        !Checks for Bottom Fluxes
        call GetData(NewProperty%ComputeOptions%BottomFluxes,                       &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'BOTTOM_FLUXES',                                &
                     Default      = .false.,                                        & 
                     ClientModule = 'ModuleReservoirs',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR90' 

        if (NewProperty%ComputeOptions%BottomFluxes .and. .not. Me%ComputeOptions%BottomFluxes) then
            call SetError(WARNING_, INTERNAL_, 'Missing keyword [BottomFluxes]', ON)            
            Me%ComputeOptions%BottomFluxes = ON
        end if
        
        if (NewProperty%ComputeOptions%BottomFluxes) then
!~             if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
            if(.not. NewProperty%ID%IsParticulate) then
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'ModuleReservoirs - ConstructPropertyValues - ERR100' 
            end if            
        endif

       !in Reservoirs all properties recognized by the model as particulate need to
       !have Bottom Fluxes because if all water exits reservoir the mass needs to go somewhere
       !and so needs the bottom concentration 
!~         if(Check_Particulate_Property(NewProperty%ID%IDNumber) .and.  &
       if (NewProperty%ID%IsParticulate .and. (.not. NewProperty%ComputeOptions%BottomFluxes)) then 
            write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' has not BOTTOM_FLUXES ON'
            write(*,*) 'but is recognised by the model as particulate.'
            write(*,*) 'Particulated recognized properties can accumulate in bottom and'
            write(*,*) 'need BOTTOM_FLUXES to be active for the propery in Reservoirs'
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR110'
        end if      

ifB:    if (NewProperty%ComputeOptions%BottomFluxes) then
            Me%ComputeOptions%BottomFluxes = .true.

            allocate (NewProperty%BottomConc     (1:Me%nReservoirs))
            allocate (NewProperty%ErosionRate    (1:Me%nReservoirs))
            allocate (NewProperty%DepositionRate (1:Me%nReservoirs))
            allocate (NewProperty%Ws             (1:Me%nReservoirs))

            NewProperty%BottomConc              = 0.0
            NewProperty%ErosionRate             = 0.0
            NewProperty%DepositionRate          = 0.0
            NewProperty%Ws                      = 0.0


            !Bottom Initial Concentration
            call GetData(BottomInitialConc,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOTTOM_CONC',                                  &
                         Default      =  0.0,                                           & 
                         ClientModule = 'ModuleReservoirs',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR120' 
            
            NewProperty%BottomConc = BottomInitialConc

            !Bottom Initial Concentration
            call GetData(NewProperty%BottomMinConc,                                     &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOTTOM_MIN_CONC',                              &
                         Default      =  0.0,                                           & 
                         ClientModule = 'ModuleReservoirs',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR130' 

            !Compute erosion fluxes
            call GetData(NewProperty%ComputeOptions%Erosion,                            &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'EROSION',                                      &
                         Default      =  .false.,                                       & 
                         ClientModule = 'ModuleReservoirs',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR140' 

            if (NewProperty%ComputeOptions%Erosion) then
                !Critial Erosion Shear Stress [Pa]
                call GetData(NewProperty%ErosionCriticalShear,                              &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'CRIT_SS_EROSION',                              &
                             Default      =  0.2,                                           & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR150' 

                !Erosion Coefficient [kg m-2 s-1]
                call GetData(NewProperty%ErosionCoefficient,                                &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'EROSION_COEF',                                 &
                             Default      =  5.0E-4,                                        & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR160' 
            
            end if


            !Compute deposition fluxes
            call GetData(NewProperty%ComputeOptions%Deposition,                         &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'DEPOSITION',                                   &
                         Default      =  .false.,                                       & 
                         ClientModule = 'ModuleReservoirs',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR170' 

            if (NewProperty%ComputeOptions%Deposition) then
                !Critial Deposition Shear Stress [Pa]
                call GetData(NewProperty%DepositionCriticalShear,                           &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'CRIT_SS_DEPOSITION',                           &
                             Default      =  0.1,                                           & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR180' 

                if (NewProperty%ComputeOptions%Erosion .and. &
                    NewProperty%DepositionCriticalShear >= NewProperty%ErosionCriticalShear) then
                    write (*,*) '[CRIT_SS_EROSION] must be higher than [CRIT_SS_DEPOSITION]'
                    stop 'ModuleReservoirs - ConstructPropertyValues - ERR14' 
                end if
    

                !See ModuleFreeVerticalMovement - Hindered settling  - CHS - kg m-3
                call GetData(NewProperty%CHS,                                               &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'CHS',                                          &
                             Default      =  4.0,                                           & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR190' 

                !See ModuleFreeVerticalMovement - Hindered settling  - CHS
                !Settling type: WSConstant = 1, SPMFunction = 2
                call GetData(NewProperty%Ws_Type,                                           &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'WS_TYPE',                                      &
                             Default      =  WSConstant,                                    & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR200' 

                !See ModuleFreeVerticalMovement - Constant settling velocity [m s-1]
                call GetData(NewProperty%Ws_Value,                                          &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'WS_VALUE',                                     &
                             Default      =  0.0001,                                        & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR210' 


                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%KL,                                                &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'KL',                                           &
                             Default      =  0.1,                                           & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR220' 

                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%KL1,                                               &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'KL1',                                          &
                             Default      =  0.1,                                           & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR230' 

                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%ML,                                                &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'ML',                                           &
                             Default      =  4.62,                                          & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR240' 


                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%M,                                                 &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'M',                                            &
                             Default      =  1.0,                                           & 
                             ClientModule = 'ModuleReservoirs',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR250' 

            end if

        end if ifB
        
        !Checks for Water Quality model
        call GetData(NewProperty%ComputeOptions%WaterQuality,                       &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'WATER_QUALITY',                              &
                     ClientModule   = 'ModuleReservoirs',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR260' 
        
        if (NewProperty%ComputeOptions%WaterQuality) then
            Me%ComputeOptions%WaterQuality = .true.
        end if


        !Checks for Benthos model
        call GetData(NewProperty%ComputeOptions%Benthos,                            &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'BENTHOS',                                    &
                     ClientModule   = 'ModuleReservoirs',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR270' 
        
        if (NewProperty%ComputeOptions%Benthos) then
            
            Me%ComputeOptions%Benthos = .true.

        end if


        !Checks for CEQUAL_W2 model
        call GetData(NewProperty%ComputeOptions%CeQualW2,                           &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'CEQUALW2',                                   &
                     ClientModule   = 'ModuleReservoirs',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR280' 
        
        if (NewProperty%ComputeOptions%CeQualW2) then
            Me%ComputeOptions%CeQualW2 = .true.
        end if

        !Checks for Life model
        call GetData(NewProperty%ComputeOptions%Life,                               &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'LIFE',                                       &
                     ClientModule   = 'ModuleReservoirs',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR290' 
        
        if (NewProperty%ComputeOptions%Life) then
            Me%ComputeOptions%Life = .true.
        end if

        !Checks if user wants to calculate total Concentration (Column + Bottom)
        call GetData(NewProperty%ComputeOptions%SumTotalConc,                       &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'SUMTOTALCONC',                                 &
                     Default      =  .FALSE.,                                       & 
                     ClientModule = 'ModuleReservoirs',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleReservoirs - ConstructPropertyValues - ERR300'

        if (NewProperty%ComputeOptions%SumTotalConc) then
!~             if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
            if (.not. NewProperty%ID%IsParticulate) then
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' is not'
                write(*,*) 'recognised as PARTICULATE and does not have Bottom_ or total_Conc'
                stop 'ModuleReservoirs - ConstructPropertyValues - ERR16b' 
            end if            
            
            Me%ComputeOptions%SumTotalConc = .true.
        endif


        !IS Coeficient
        call GetData(NewProperty%IScoefficient,                                     &
                     Me%ObjEnterData, iflag,                                        &
                     KeyWord        = 'IS_COEF',                                    &
                     Default        = 1.e-3,                                        &      
                     SearchType     = FromBlock,                                    &
                     ClientModule   = 'ModuleReservoirs',                      &
                     STAT           = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR310' 

        !IS ExtinctionCoef
        call GetData(NewProperty%ExtinctionCoefficient,                             &
                     Me%ObjEnterData, iflag,                                        &
                     KeyWord        = 'EXTINCTION_PARAMETER',                       &
                     Default        = 1.0,                                          &      
                     SearchType     = FromBlock,                                    &
                     ClientModule   = 'ModuleReservoirs',                      &
                     STAT           = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR320' 

        call GetData(NewProperty%ComputeOptions%TimeSerie,                          &
             Me%ObjEnterData, iflag,                                                &
             Keyword        = 'TIME_SERIE',                                         &
             ClientModule   = 'ModuleWaterProperties',                              &
             SearchType     = FromBlock,                                            &
             Default        = .false.,                                              &
             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR330'

        call GetData(NewProperty%ComputeOptions%HDF,                                &
             Me%ObjEnterData, iflag,                                                &
             Keyword        = 'OUTPUT_HDF',                                         &
             ClientModule   = 'ModuleWaterProperties',                              &
             SearchType     = FromBlock,                                            &
             Default        = .false.,                                              &
             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR335'        
        
        call GetData(NewProperty%OutputName,                                        &
             Me%ObjEnterData, iflag,                                                &
             Keyword        = 'OUTPUT_NAME',                                        &
             ClientModule   = 'ModuleWaterProperties',                              &
             SearchType     = FromBlock,                                            &
             Default        = 'NAME',                                               &
             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR340'

        if (NewProperty%OutputName /= 'NAME' .and. NewProperty%OutputName /= 'DESCRIPTION') &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR350'

        ModelDT = Me%ExtVar%DT

        call GetData(NewProperty%DTInterval,                                         &
                     Me%ObjEnterData, iflag,                                         &
                     SearchType   = FromBlock,                                       &
                     keyword      = 'DTINTERVAL',                                    &
                     Default      = ModelDT,                                         &
                     ClientModule = 'ModuleReservoirs',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                 &
            stop 'ModuleReservoirs - ConstructPropertyValues - ERR390'
        
        if (iflag == 1) then
            NewProperty%ComputeOptions%DTIntervalAssociated = .true.
            Me%ComputeOptions%DTIntervalAssociated          = .true.                           
        
            if (NewProperty%DTInterval < ModelDT) then
                write(*,*) 
                write(*,*) 'Property time step is smaller then model time step'
                stop 'ModuleReservoirs - ConstructPropertyValues - ERR400'

            elseif (NewProperty%DTInterval > ModelDT) then 

                !Property time step must be a multiple of the model time step
                auxFactor = NewProperty%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) 'Property time step must be a multiple of model time step.'
                    write(*,*) 'Please review your input data.'
                    stop 'ModuleReservoirs - ConstructPropertyValues - ERR410'
                endif

                !Run period in seconds
                DTaux = Me%ExtVar%EndTime - Me%ExtVar%ActualTime

                !The run period   must be a multiple of the Property DT
                auxFactor = DTaux / NewProperty%DTInterval

                ErrorAux = auxFactor - int(auxFactor)
                if (ErrorAux /= 0) then

                    write(*,*) 
                    write(*,*) 'Property time step is not a multiple of model time step.'
                    stop 'ModuleReservoirs - ConstructPropertyValues - ERR420'
                end if
            endif

            NewProperty%NextCompute = Me%ExtVar%ActualTime + NewProperty%DTInterval
        endif
        
    end subroutine ConstructPropertyValues    
    
    !---------------------------------------------------------------------------    
    
    subroutine ReadInitialConc(NewProperty)

        !Arguments-------------------------------------------------------------
        type (T_Property)                           :: NewProperty
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: EXIST
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ                      
        real, dimension(:,:), pointer               :: ReservoirConc, ReservoirBottomConc      
        type (T_Reservoir), pointer                 :: CurrReservoir   
        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%ExtVar%WorkSize2D%ILB 
        WorkIUB = Me%ExtVar%WorkSize2D%IUB 

        WorkJLB = Me%ExtVar%WorkSize2D%JLB 
        WorkJUB = Me%ExtVar%WorkSize2D%JUB 

        !----------------------------------------------------------------------

        
        allocate(ReservoirConc           (WorkILB:WorkIUB, WorkJLB:WorkJUB))
        allocate(ReservoirBottomConc     (WorkILB:WorkIUB, WorkJLB:WorkJUB))
        call SetMatrixValue (ReservoirConc,          Me%ExtVar%WorkSize2D, 0.0)
        call SetMatrixValue (ReservoirBottomConc,    Me%ExtVar%WorkSize2D, 0.0)                


        inquire (FILE=trim(Me%Files%InitialFile), EXIST = Exist)

cd0:    if (Exist) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%InitialFile),                              &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialConc - ModuleReservoirs - ERR01'


            ! Reads from HDF file the reservoir volume
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialConc - ModuleReservoirs - ERR02'


            call HDF5ReadData   (ObjHDF5, "/Results/"//trim(adjustl(NewProperty%ID%Name)),        &
                                    trim(adjustl(NewProperty%ID%Name)),                              &
                                    Array2D = ReservoirConc,                                         &
                                    STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'ReadInitialConc - ModuleReservoirs - ERR020'                
                
            if (NewProperty%ComputeOptions%BottomFluxes) then  
                call HDF5ReadData   (ObjHDF5, "/Results/"//trim(adjustl(NewProperty%ID%Name))//"_Bottom",        &
                                        trim(adjustl(NewProperty%ID%Name))//"_Bottom",                              &
                                        Array2D = ReservoirBottomConc,                                 &
                                        STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                          &
                    stop 'ReadInitialConc - ModuleReservoirs - ERR030'                      
            endif
                
            CurrReservoir => Me%FirstReservoir
            do while (associated(CurrReservoir))
                
                NewProperty%Concentration(CurrReservoir%Position) = ReservoirConc(CurrReservoir%GridI, CurrReservoir%GridJ)
                                                                           
                if(NewProperty%ComputeOptions%BottomFluxes) then 
                        
                    NewProperty%BottomConc(CurrReservoir%Position) = ReservoirBottomConc(CurrReservoir%GridI,CurrReservoir%GridJ)

                endif                                            
                        
                        
                CurrReservoir => CurrReservoir%Next
            enddo                                                                                               

            deallocate(ReservoirConc           )   
            deallocate(ReservoirBottomConc     )  
            
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialConc - ModuleReservoirs - ERR040'

        else
            
            write(*,*)
            stop 'ReadInitialConc - ModuleReservoirs - ERR050'

        end if cd0
    
    end subroutine ReadInitialConc
    
    !--------------------------------------------------------------------------    
    
    subroutine ConstructGlobalVariables
        
        !External-----------------------------------------------------------------
        integer                                :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------      
        
        !File wth reservoirs info
        call GetData(Me%Files%ReservoirFile,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      ='RESERVOIR_FILE',                                    &
                     ClientModule = 'ModuleReservoirs',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructGlobalVariables - ModuleReservoirs - ERR05'        

        if (iflag /= 1) then
            write (*,*)'Not Found reservoir file [RESERVOIR_FILE]'
            stop 'ConstructGlobalVariables - ModuleReservoirs - ERR01'
        endif               
        
        !allow discharges
        call GetData(Me%ComputeOptions%Discharges,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'DISCHARGES',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleReservoirs',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR07'         
        
        !Checks for Surface Fluxes
        call GetData(Me%ComputeOptions%SurfaceFluxes,                               &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'SURFACE_FLUXES',                               &
                     Default      = .false.,                                        & 
                     ClientModule = 'ModuleReservoirs',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ConstructPropertyValues - ERR08'         
        
        !Checks for Bottom Fluxes
        call GetData(Me%ComputeOptions%BottomFluxes,                                &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'BOTTOM_FLUXES',                                &
                     Default      = .false.,                                        & 
                     ClientModule = 'ModuleReservoirs',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ConstructPropertyValues - ERR09'            
        
        !allow discharges
        call GetData(Me%ComputeOptions%PropertyComputeMethod,                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'PROP_COMPUTE_METHOD',                            &
                     Default        = InstantMixing_,                                   &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleReservoirs',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR09a' 
        
        if (Me%ComputeOptions%PropertyComputeMethod /= InstantMixing_ .and.   &
            Me%ComputeOptions%PropertyComputeMethod /= RetentionTimeMixing_) then
            write (*,*)'Unknown [PROP_COMPUTE_METHOD]'
            stop 'ConstructGlobalVariables - ModuleReservoirs - ERR09b'            
        endif
        
        !Is a continuous simulation
        call GetData(Me%Continuous,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'CONTINUOUS',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleReservoirs',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR10'           
        
        if (Me%Continuous) then   

            !Reads the name of the file where to read initial data
            call ReadFileName ('RESERVOIRS_INI', Me%Files%InitialFile, 'Reservors Initial File', STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialConc - ModuleReservoirs - ERR20'             
            
        else
            
            !Initial conditions when not defined already (in reservoir file or in timeseries)
            call GetData(Me%InitialVolumeDefaultMethod,                                     &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'INITIAL_VOLUME_DEFAULT_METHOD',                  &
                         Default        = StartPercentageFull_,                             &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleReservoirs',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR30'                   
            
            if (Me%InitialVolumeDefaultMethod /= StartPercentageFull_) then
                write (*,*)'Unknown [INITIAL_CONDITION]'
                stop 'ConstructGlobalVariables - ModuleReservoirs - ERR035'                
            endif
            
            if (Me%InitialVolumeDefaultMethod == StartPercentageFull_) then
                call GetData(Me%StartPercentageFull,                                            &
                             Me%ObjEnterData, iflag,                                            &
                             Keyword        = 'START_PERCENTAGE_FULL',                          &
                             Default        = 50.,                                              &
                             SearchType     = FromFile,                                         &
                             ClientModule   = 'ModuleReservoirs',                               &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR40'                    
            endif
            
        endif
        

        !Output
        call GetOutPutTime(Me%ObjEnterData,                                     &
                           CurrentTime   = Me%ExtVar%ActualTime,                &
                           EndTime       = Me%ExtVar%EndTime,                   &
                           keyword       = 'OUTPUT_TIME',                       &
                           SearchType    = FromFile,                            &
                           OutPutsTime   = Me%OutPut%OutTime,                   &
                           OutPutsOn     = Me%OutPut%HDF,                       &
                           OutPutsNumber = Me%OutPut%Number,                    &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR50' 

        if (Me%OutPut%HDF) Me%OutPut%NextOutPut = 1
           

        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%ExtVar%ActualTime,                         &
                           EndTime      = Me%ExtVar%EndTime,                             &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR51'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleDrainageNetwork',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ConstructGlobalVariables - ModuleReservoirs - ERR52' 
        
        
        call GetData(Me%OutPut%TimeSerie,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TIME_SERIE',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleReservoirs',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleReservoirs - ERR60'        
        
    end subroutine ConstructGlobalVariables
    
    !-------------------------------------------------------------------------
    
    subroutine GetImposedReservoirlevel(NewReservoir, ClientID, Construct)
    

        !Arguments-------------------------------------------------------------
        type(T_Reservoir), pointer                  :: NewReservoir
        integer, optional                           :: ClientID
        logical                                     :: Construct
        
        !Local------------------------------------------------------------------
        integer                                         :: NumberOfSources, index, iflag
        integer                                         :: STAT_CALL
        logical                                         :: FoundBlock
        real                                            :: ReservoirVolume
 
        !Begin------------------------------------------------------------------
        
           
        
        if (Construct) then
            
            
            call RewindBlock(Me%ObjEnterDataReservoirFile, ClientId, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR01'             
            
           !Gets the number of Source blocks
            call GetNumberOfBlocks(Me%ObjEnterDataReservoirFile, "<<beginleveltimeseries>>", "<<endleveltimeseries>>",   &
                                   FromBlock, NumberOfSources,    &
                                   ClientID, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR010'        
        
            if (NumberOfSources > 0) then
          
                do index = 1, NumberOfSources

                    call ExtractBlockFromBlock(Me%ObjEnterDataReservoirFile,                    &
                                               ClientNumber      = ClientID,                    &
                                               block_begin       = "<<beginleveltimeseries>>",    &
                                               block_end         = "<<endleveltimeseries>>",      &
                                               BlockInBlockFound = FoundBlock,                  &
                                               STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR020'

                    if (FoundBlock) then
            

                        call GetData(NewReservoir%TimeSeries%Filename,                &
                                        Me%ObjEnterDataReservoirFile, iflag,          &
                                        SearchType   = FromBlockInBlock,              &
                                        keyword      = 'FILENAME',                    &                             
                                        ClientModule = 'ModuleReservoirs',            &
                                        STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR025'
                        
                        call GetData(NewReservoir%TimeSeries%Column,                  &
                                        Me%ObjEnterDataReservoirFile , iflag,         &
                                        SearchType   = FromBlockInBlock,              &
                                        keyword      = 'DATA_COLUMN',                 &
                                        ClientModule = 'ModuleReservoirs',            &
                                        STAT         = STAT_CALL)                                      
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR026'
                        if (iflag /= 1) &
                            stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR027'                

                        !Starts Time Serie
                        NewReservoir%TimeSeries%ID = 0
                        call StartTimeSerieInput(NewReservoir%TimeSeries%ID,       &
                                                    NewReservoir%TimeSeries%FileName, &
                                                    Me%ObjTime,                        &
                                                    STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR028'
                    
                    
                        !Get value for start instant
                        ReservoirVolume = GetReservoirVolume(NewReservoir, Me%ExtVar%ActualTime)
                    
                        !apply reservoir value. it will overwrite value defined in reservoir file (if defined)
                        if (ReservoirVolume <= NewReservoir%MaxVolume .and. ReservoirVolume >= 0.0) then
                            NewReservoir%InitialVolume = ReservoirVolume
                            NewReservoir%InitialVolumeDefined = .true.      
                            
                            NewReservoir%Management%ON                    = .true.
                            NewReservoir%Management%ImposedLevel          = .true.                              
                        endif
                        
                      
                            
                        !call KillTimeSerie (ObjTimeSerie, STAT = STAT_CALL)
                        !if (STAT_CALL /= SUCCESS_) stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR30'
                        

                    
                    else

                        stop 'GetImposedReservoirLevel - ModuleReservoirs - ERR040'

                    endif

                enddo          
        
            endif
            
        else
            
            if (NewReservoir%Management%ImposedLevel) then
                
                !Get value for timestep start
                ReservoirVolume = GetReservoirVolume(NewReservoir, Me%ExtVar%ActualTime - Me%ExtVar%DT)
                   
                !apply reservoir value. 
                if (ReservoirVolume <= NewReservoir%MaxVolume .and. ReservoirVolume >= 0.0) then
                    NewReservoir%VolumeOld = ReservoirVolume   
                    NewReservoir%VolumeNew = NewReservoir%VolumeOld
                endif            
            
                !Get value for timestep end
                !Get reservoir volume that will be the target after all fluxes to compute outflow
                NewReservoir%VolumeTarget = GetReservoirVolume(NewReservoir, Me%ExtVar%ActualTime)
            
            endif
            
        endif
        
    end subroutine GetImposedReservoirLevel
    
    !-------------------------------------------------------------------------
    
    real function GetReservoirVolume(CurrReservoir, Time)
    
        !Arguments----------------------------------------------------------------
        type(T_Reservoir),           pointer     :: CurrReservoir
        type(T_Time)                             :: Time    
        
        !Local--------------------------------------------------------------------
        type (T_Time)                                   :: Time1, Time2
        real                                            :: Value1, Value2
        logical                                         :: TimeCycle        
        real                                            :: TimeSerieValue
        integer                                         :: STAT_CALL
        !Begin--------------------------------------------------------------------

        call GetTimeSerieValue (CurrReservoir%TimeSeries%ID, Time,          &
                                CurrReservoir%TimeSeries%Column,             &
                                Time1, Value1, Time2, Value2, TimeCycle,     &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetReservoirVolume - ModuleReservoirs - ERR29'
     
        if (TimeCycle) then
            TimeSerieValue = Value1

        else

            !Interpolates Value for current instant
            call InterpolateValueInTime(Time, Time1, Value1, Time2, Value2, TimeSerieValue)

        endif     
                        
        !Compute Volume from level
        GetReservoirVolume = ComputeReservoirVolume(CurrReservoir, TimeSerieValue)    
    
    end function 
    
    !-------------------------------------------------------------------------
    
    subroutine VerifyOptions
        
        !Local----------------------------------------------------------------
        type(T_Reservoir),           pointer     :: CurrentReservoir
        integer                                  :: nImposed
    
       CurrentReservoir => Me%FirstReservoir
        do while (associated(CurrentReservoir))
            
            nImposed = 0
            if (CurrentReservoir%Management%ImposedLevel) then
                nImposed = nImposed + 1
            endif
            
            if (CurrentReservoir%Management%ImposedOutflow) then
                nImposed = nImposed + 1
            endif
            
            if (CurrentReservoir%Management%ImposedOperation) then
                nImposed = nImposed + 1
            endif       
            
            !inconsistent options  can only impose one
            if (nImposed > 1) then
                write (*,*) 'Can only impose outflow OR level OR operation curves in'
                write (*,*) 'Reservoir ID = ', CurrentReservoir%ID            
                stop 'VerifyOptions - ModuleReservoirs - ERR01'                        
            endif
                    
            
            !need acc curve?
            if (CurrentReservoir%IsWeir .or. CurrentReservoir%Management%OperationType == Operation_Level_Outflow_  &
                .or. CurrentReservoir%Management%OperationType == Operation_Level_PercInflow_                       &
                .or. CurrentReservoir%Management%ImposedLevel) then
                
                if (.not. associated(CurrentReservoir%Management%AccVolumeCurve)) then
                    write(*,*) 'Not Found Reservoir accumulated volumes curve'
                    write(*,*) 'It is mandatory when reservoir is weir, level is imposed or operation curve depends on level'
                    write(*,*)  'in reservoir ID : ', CurrentReservoir%ID
                    stop 'VerifyOptions - ModuleReservoirs - ERR010'   
                endif
            endif
            CurrentReservoir => CurrentReservoir%Next
        enddo    
    
    end subroutine VerifyOptions
    
    !-------------------------------------------------------------------------
    
    subroutine AddReservoir(NewReservoir)

        !Arguments--------------------------------------------------------------
        type(T_Reservoir),           pointer     :: NewReservoir


        if (.not.associated(Me%FirstReservoir)) then
            Me%nReservoirs  = 1
            Me%FirstReservoir      => NewReservoir
            Me%LastReservoir      => NewReservoir
        else
            NewReservoir%Prev     => Me%LastReservoir
            Me%LastReservoir%Next => NewReservoir
            Me%LastReservoir      => NewReservoir
            Me%nReservoirs        = Me%nReservoirs + 1            
        end if 

        NewReservoir%Position = Me%nReservoirs
        
    end subroutine AddReservoir     
    
    !-------------------------------------------------------------------------      
    
    subroutine AddProperty(NewProperty)

        !Arguments--------------------------------------------------------------
        type(T_Property),           pointer     :: NewProperty


        if (.not.associated(Me%FirstProperty)) then
            Me%nProperties  = 1
            Me%FirstProperty     => NewProperty
            Me%LastProperty      => NewProperty
        else
            NewProperty%Prev     => Me%LastProperty
            Me%LastProperty%Next => NewProperty
            Me%LastProperty      => NewProperty
            Me%nProperties       = Me%nProperties + 1
        end if 

    end subroutine AddProperty     
    
    !-------------------------------------------------------------------------    
    
    subroutine InitializeVariables(Constructing)

        !Arguments------------------------------------------------------------
        logical                                 :: Constructing
        !Local-----------------------------------------------------------------
        type(T_Reservoir),          pointer     :: CurrentReservoir
        type(T_Property),           pointer     :: CurrentProperty
        
        !Begin----------------------------------------------------------------
        
        if (Constructing) then
            
            allocate(Me%ReservoirVolumes(Me%nReservoirs))
            allocate(Me%ReservoirInflows(Me%nReservoirs))    
            allocate(Me%ReservoirOutflows(Me%nReservoirs))
            allocate(Me%ReservoirPercFull(Me%nReservoirs))
            Me%Size1D%ILB = 1
            Me%Size1D%IUB = Me%nReservoirs
            
            call SetMatrixValue(Me%ReservoirInflows, Me%Size1D, 0.0)
        
        endif
        
        call SetMatrixValue(Me%ReservoirVolumes, Me%Size1D, 0.0)        
        !Cant reset since it is imposed
        !call SetMatrixValue(Me%ReservoirInflows, Me%Size1D, 0.0)        
        call SetMatrixValue(Me%ReservoirOutflows, Me%Size1D, 0.0)
        call SetMatrixValue(Me%ReservoirPercFull, Me%Size1D, 0.0)
        
        CurrentReservoir => Me%FirstReservoir
        do while (associated(CurrentReservoir))
            
            if (Constructing) then
                CurrentReservoir%DNInflow    = 0.0
            endif
            
            CurrentReservoir%VolumeOld      = CurrentReservoir%VolumeNew
            CurrentReservoir%PercFull       = CurrentReservoir%VolumeOld / CurrentReservoir%MaxVolume  * 100
            !CurrentReservoir%DNInflow       = 0.0  ! Cant reset since it is imposed
            CurrentReservoir%Outflow        = 0.0
            CurrentReservoir%Discharges     = 0.0
            CurrentReservoir%SurfaceFluxes  = 0.0
            
            if (Me%HasProperties) then
                CurrentProperty => Me%FirstProperty
        
                do while (associated(CurrentProperty))
                    CurrentProperty%ConcentrationOld(CurrentReservoir%Position) =   &
                             CurrentProperty%Concentration(CurrentReservoir%Position)
                    CurrentProperty%MassCreated(CurrentReservoir%Position)    = 0.0
                    
                    if (CurrentProperty%ComputeOptions%BottomFluxes) then
                        CurrentProperty%ErosionRate(CurrentReservoir%Position)    = 0.0
                        CurrentProperty%DepositionRate(CurrentReservoir%Position) = 0.0
                    endif
                    
                    CurrentProperty => CurrentProperty%Next
                enddo
            endif
                
            CurrentReservoir => CurrentReservoir%Next
        enddo

        

    end subroutine InitializeVariables

    !--------------------------------------------------------------------------      
    
    subroutine ReadInitialVolume()

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: EXIST
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ                      
        real, dimension(:,:), pointer               :: ReservoirVolume   
        type (T_Reservoir), pointer                 :: CurrReservoir   
        integer                                     :: ObjHorizontalGrid = 0
        !----------------------------------------------------------------------

        
        !Constructs Horizontal Grid
        call ConstructHorizontalGrid(ObjHorizontalGrid, Me%ExtVar%TopographicFile, &
                                        STAT = STAT_CALL)           
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialVolume - ModuleReservoirs - ERR000'   
        
        !Gets the size of the grid
        call GetHorizontalGridSize (ObjHorizontalGrid,                               &
                                    Size     = Me%ExtVar%Size2D,                     &
                                    WorkSize = Me%ExtVar%WorkSize2D,                 &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialVolume - ModuleReservoirs - ERR001'   
                    
        
        
        !Bounds
        WorkILB = Me%ExtVar%WorkSize2D%ILB 
        WorkIUB = Me%ExtVar%WorkSize2D%IUB 

        WorkJLB = Me%ExtVar%WorkSize2D%JLB 
        WorkJUB = Me%ExtVar%WorkSize2D%JUB 

        !----------------------------------------------------------------------

        
        allocate(ReservoirVolume         (WorkILB:WorkIUB, WorkJLB:WorkJUB))
        call SetMatrixValue (ReservoirVolume,        Me%ExtVar%WorkSize2D, 0.0)        


        inquire (FILE=trim(Me%Files%InitialFile), EXIST = Exist)

cd0:    if (Exist) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%InitialFile),                              &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialVolume - ModuleReservoirs - ERR01'


            ! Reads from HDF file the reservoir volume
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialVolume - ModuleReservoirs - ERR02'

            call HDF5ReadData   (ObjHDF5, "/Results/volume",        &
                                 "volume",                          &
                                 Array2D = ReservoirVolume,         &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                              &
                stop 'ReadInitialVolume - ModuleReservoirs - ERR03'     
            
            
            CurrReservoir => Me%FirstReservoir
            do while (associated(CurrReservoir))
                
                CurrReservoir%VolumeOld = ReservoirVolume(CurrReservoir%GridI, CurrReservoir%GridJ)
                CurrReservoir%VolumeNew = CurrReservoir%VolumeOld
                
                CurrReservoir => CurrReservoir%Next
            enddo                                                        

            deallocate(ReservoirVolume         )
            
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialVolume - ModuleReservoirs - ERR040'

        else
            
            write(*,*)
            stop 'ReadInitialVolume - ModuleReservoirs - ERR050'

        end if cd0
        

        call KillHorizontalGrid     (ObjHorizontalGrid,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialVolume - ModuleReservoirs - ERR180'            
    
    end subroutine ReadInitialVolume
    
    !--------------------------------------------------------------------------

    subroutine InitializeReservoirsVolume()

        !Arguments------------------------------------------------------------
        !Local-----------------------------------------------------------------        
        type(T_Reservoir),           pointer    :: CurrentReservoir
        real, dimension(6), target              :: AuxTime

        !Begin----------------------------------------------------------------

        CurrentReservoir => Me%FirstReservoir
            
        !To get the first year of the simulation - Ana Oliveira
        call ExtractDate   (Me%ExtVar%BeginTime , AuxTime(1), AuxTime(2),         &
                            AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))
        
        do while (associated(CurrentReservoir))
                
            !if value defined in reservoir file or in timeseries, aplly it
            if (CurrentReservoir%InitialVolumeDefined) then
                CurrentReservoir%VolumeOld = CurrentReservoir%InitialVolume
            
            !if value was not defined, apply a default precentage of max volume
            else
                if (Me%InitialVolumeDefaultMethod == StartPercentageFull_) then
                    CurrentReservoir%VolumeOld = CurrentReservoir%MaxVolume * Me%StartPercentageFull / 100.0
                else
                    write(*,*)'Not recognized Initial Condition Method'
                    stop 'InitializeReservoirs - Module Reservoirs - ERR01'
                endif 
            endif
            
            !If the year of construction is bigger than the current year the percentage full must be 0 - Ana Oliveira
            if (CurrentReservoir%ConstructionYear <= AuxTime(1)) then
                CurrentReservoir%VolumeNew = CurrentReservoir%VolumeOld
            else
                CurrentReservoir%VolumeNew = 0.0
            endif
            
            if (associated(CurrentReservoir%Management%AccVolumeCurve)) then
                CurrentReservoir%WaterLevel = ComputeReservoirLevel(CurrentReservoir)
            endif
                        
            CurrentReservoir => CurrentReservoir%Next
        enddo                

    end subroutine InitializeReservoirsVolume

    !--------------------------------------------------------------------------    
    
    real function ComputeReservoirLevel(CurrentReservoir)

        !Arguments------------------------------------------------------------
        type(T_Reservoir),           pointer     :: CurrentReservoir

        !Local-----------------------------------------------------------------
        integer                                 :: i
        real                                    :: WaterLevel, ReservoirVolume, PreviousCurveLevel
        real                                    :: PreviousCurveVolume, NextCurveLevel, NextCurveVolume
        !Begin----------------------------------------------------------------
        WaterLevel = null_real
        if (associated(CurrentReservoir%Management%AccVolumeCurve)) then        
            
            ReservoirVolume              = CurrentReservoir%VolumeNew
            PreviousCurveLevel           = CurrentReservoir%Management%AccVolumeCurve(1, 2)
            PreviousCurveVolume          = CurrentReservoir%Management%AccVolumeCurve(1, 1)
            
            !go trough all points to find where belongs
            do i = 2, CurrentReservoir%Management%AccVolumeCurvePoints
                
                NextCurveLevel           = CurrentReservoir%Management%AccVolumeCurve(i, 2)
                NextCurveVolume          = CurrentReservoir%Management%AccVolumeCurve(i, 1)
                
                if (ReservoirVolume >= PreviousCurveVolume &
                    .and. ReservoirVolume <= NextCurveVolume) then
                    
                    WaterLevel = LinearInterpolation(PreviousCurveVolume, PreviousCurveLevel,       &
                                        NextCurveVolume, NextCurveLevel, ReservoirVolume)                           
                    
                endif
                
                PreviousCurveLevel  = NextCurveLevel
                PreviousCurveVolume = NextCurveVolume
            
            enddo
               
        else
            !not all reservoitrs need to have it defined
            WaterLevel = null_real
        endif
            
        ComputeReservoirLevel = WaterLevel

    end function ComputeReservoirLevel
    
    !---------------------------------------------------------------------------
    
    real function ComputeReservoirVolume(CurrentReservoir, ReservoirLevel)

        !Arguments------------------------------------------------------------
        type(T_Reservoir),           pointer     :: CurrentReservoir
        real                                     :: ReservoirLevel

        !Local-----------------------------------------------------------------
        integer                                 :: i
        real                                    :: WaterVolume, PreviousCurveLevel
        real                                    :: PreviousCurveVolume, NextCurveLevel, NextCurveVolume
        !Begin----------------------------------------------------------------

        if (associated(CurrentReservoir%Management%AccVolumeCurve)) then        
            
            PreviousCurveLevel           = CurrentReservoir%Management%AccVolumeCurve(1, 2)
            PreviousCurveVolume          = CurrentReservoir%Management%AccVolumeCurve(1, 1)
            
            !go trough all points to find where belongs
            do i = 2, CurrentReservoir%Management%AccVolumeCurvePoints
                
                NextCurveLevel           = CurrentReservoir%Management%AccVolumeCurve(i, 2)
                NextCurveVolume          = CurrentReservoir%Management%AccVolumeCurve(i, 1)
                
                if (ReservoirLevel >= PreviousCurveLevel &
                    .and. ReservoirLevel <= NextCurveLevel) then
                    
                    WaterVolume = LinearInterpolation(PreviousCurveLevel, PreviousCurveVolume,       &
                                        NextCurveLevel, NextCurveVolume, ReservoirLevel)                           
                    
                endif
                
                PreviousCurveLevel  = NextCurveLevel
                PreviousCurveVolume = NextCurveVolume
            
            enddo
               
        else
            !not all reservoitrs need to have it defined
            WaterVolume = null_real
        endif
            
        ComputeReservoirVolume = WaterVolume

    end function ComputeReservoirVolume
    
    !---------------------------------------------------------------------------    
    
    subroutine ConstructOutput

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------

        if (Me%Output%TimeSerie) then
            !Verifies in whichs nodes to write time series 
            call ReadTimeSerieReservoirList
       
            !Basic Time Series Data (Hydrodynamic Properties) & Transported Properties
            call ConstructTimeSerieList

            !Opens all Time Series Data Files
            call ConstructTimeSeries
        endif
        
        !Opens HDF5 File
        if (Me%Output%HDF .and. Me%Output%HDFActive) then
            call ConstructHDF5Output 
        endif
        
       
    end subroutine ConstructOutput

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ReadTimeSerieReservoirList

        !This subroutine reads selected NodeID and selects also the downstream reach.
        !Local------------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound, Found
        integer                                     :: FirstLine, LastLine            
        integer                                     :: STAT_CALL, flag
        integer                                     :: ReservoirID
        type (T_Reservoir), pointer                 :: CurrReservoir              
        character(LEN = StringLength)               :: block_begin          = '<beginreservoirtimeserie>'
        character(LEN = StringLength)               :: block_end            = '<endreservoirtimeserie>'  
 
        !-----------------------------------------------------------------------

        
        call GetData(Me%TimeSerie%Location,                                 &
                     Me%ObjEnterData, flag,                                 &
                     SearchType   = FromFile,                               &
                     keyword      = 'TIME_SERIE_LOCATION',                  &
                     ClientModule = 'ModuleReservoirs',                &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'ReadTimeSerieReservoirList - ModuleReservoirs - ERR02'
            
if1:    if (flag==1) then   
     
            Me%TimeSerie%nNodes =  0

            call ConstructEnterData (Me%TimeSerie%ObjEnterData,             &
                                     Me%TimeSerie%Location, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadTimeSerieReservoirList - ModuleReservoirs - ERR03'
        

do2:        do
            
                !Nodes--------------------------------------------------------------
                call ExtractBlockFromBuffer(Me%TimeSerie%ObjEnterData, ClientNumber,    &
                                            block_begin, block_end,       &
                                            BlockFound, FirstLine, LastLine, STAT_CALL) 

                if (STAT_CALL .EQ. SUCCESS_) then    

                    if (BlockFound) then                 
                        
                        call GetData(ReservoirID, Me%TimeSerie%ObjEnterData,             &
                                     flag,                                          &
                                     keyword      = 'RESERVOIR_ID',                 &
                                     ClientModule = 'ModuleReservoirs',             &
                                     SearchType   = FromBlock,                      &
                                     STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                  &
                            stop 'ReadTimeSerieReservoirList - ModuleReservoirs - ERR04'
                        
                        call FindReservoir (ReservoirID, CurrReservoir, Found)

                        if (.NOT.Found) then
                            write (*,*) 'Reservoir not found'
                            write (*,*) 'Reservoir ID = ', ReservoirID                           
                            stop 'ReadTimeSerieReservoirList - ModuleDrainageNetwork - ERR06'
                        end if                
                        
                        if (CurrReservoir%TimeSerie) then
                            write (*,*) 'Repeated reservoir in time series: ', CurrReservoir%ID
                            stop 'ReadTimeSerieReservoirList - ModuleReservoirs - ERR07'
                        end if
                            
                        CurrReservoir%TimeSerie = .TRUE.
                        Me%TimeSerie%nNodes = Me%TimeSerie%nNodes + 1
                        
                        call GetData(CurrReservoir%TimeSeriesName,                  &
                                     Me%TimeSerie%ObjEnterData,                     &
                                     flag,                                          &
                                     keyword      = 'NAME',                         &
                                     ClientModule = 'ModuleReservoirs',             &
                                     SearchType   = FromBlock,                      &
                                     STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                  &
                            stop 'ReadTimeSerieReservoirList - ModuleReservoirs - ERR08'
                                                                                          
                    else 

                        call Block_Unlock(Me%TimeSerie%ObjEnterData, ClientNumber)
                        exit do2
                        
                    end if

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then 

                    stop 'ReadTimeSerieReservoirList - ModuleReservoirs - ERR11'

                end if
            enddo do2
        end if if1


    end subroutine ReadTimeSerieReservoirList

    !---------------------------------------------------------------------------
    
   subroutine ConstructTimeSerieList

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer                          :: Property
        integer                                             :: i
        character(LEN = StringLength)                       :: aux, str_Length, TextFormat
        type(T_Reservoir  ), pointer                        :: CurrReservoir     

        
        !Basic output: reservoir volume, percent full, level, inflow, discharges, surface fluxes, outflow
        Me%TimeSerie%nProp = 7

        if (Me%HasProperties) then
            Property => Me%FirstProperty
            do while (associated (Property))            
                if (Property%ComputeOptions%TimeSerie   ) then
                    Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1
                end if
                Property => Property%Next    
            end do                   
            
        end if


        allocate (Me%TimeSerie%ObjTimeSerie (1:Me%TimeSerie%nNodes))
        allocate (Me%TimeSerie%Name         (1:Me%TimeSerie%nNodes))
        allocate (Me%TimeSerie%X            (1:Me%TimeSerie%nNodes))
        allocate (Me%TimeSerie%Y            (1:Me%TimeSerie%nNodes))

        Me%TimeSerie%ObjTimeSerie = 0
        i = 0
        
        CurrReservoir => Me%FirstReservoir        
        do while(associated(CurrReservoir))
            
            
            if (CurrReservoir%TimeSerie) then
                    
                i = i + 1                    
                aux = ''
                write (str_Length, '(i10)') StringLength
                TextFormat   = '(a'//trim(adjustl(adjustr(str_Length)))//')'
                    
                write(aux,TextFormat) CurrReservoir%TimeSeriesName                    
                Me%TimeSerie%Name(i)  = 'Reservoir_'//trim(adjustl(adjustr(aux)))
                Me%TimeSerie%X(i)     = CurrReservoir%CoordinateX
                Me%TimeSerie%Y(i)     = CurrReservoir%CoordinateY
            end if
            
            CurrReservoir => CurrReservoir%Next
        end do


    end subroutine ConstructTimeSerieList

    !---------------------------------------------------------------------------    
    
    
    subroutine ConstructTimeSeries

        !Local------------------------------------------------------------------
        integer                                              :: STAT_CALL      
        integer                                              :: i
        character(LEN = StringLength), dimension(:), pointer :: PropHeaderList 
       
        
if1:    if (Me%TimeSerie%nNodes > 0) then
            
            allocate (PropHeaderList           (1:Me%TimeSerie%nProp ))
            allocate (Me%TimeSerie%DataLine    (1:Me%TimeSerie%nProp ))            

            call FillPropNameVector (PropHeaderList)

            do i = 1, Me%TimeSerie%nNodes

                call StartTimeSerie(Me%TimeSerie%ObjTimeSerie(i), Me%ObjTime,       &
                                    TimeSerieDataFile = trim(Me%TimeSerie%Location),&
                                    PropertyList      = PropHeaderList,             &
                                    Extension         = "sra",                      &
                                    ResultFileName    = Me%TimeSerie%Name (i),      &
                                    ModelName         = Me%ModelName,               &
                                    CoordX            = Me%TimeSerie%X    (i),      &
                                    CoordY            = Me%TimeSerie%Y    (i),      &
                                    STAT              = STAT_CALL)
                if (STAT_CALL /= 0) stop 'ConstructTimeSeries - ModuleReservoirs- ERR01'
                
            end do

            deallocate (PropHeaderList)
            
            
        end if if1


    end subroutine ConstructTimeSeries
    
    !---------------------------------------------------------------------------
   
    subroutine FillPropNameVector (PropVector)

        !Arguments--------------------------------------------------------------
        character(StringLength) , dimension (:), pointer    :: PropVector  

        !Local------------------------------------------------------------------
        type (T_Property), pointer                          :: Property
        integer                                             :: i   


        PropVector (1) = 'Volume'
        PropVector (2) = 'Percent Full'
        PropVector (3) = 'Water Level'
        PropVector (4) = 'Inflow'
        PropVector (5) = 'Discharges'
        PropVector (6) = 'Surface Fluxes'
        PropVector (7) = 'Outflow'
        i = 7

if0:    if (Me%HasProperties) then

            Property => Me%FirstProperty
            i = i + 1
            
            do while (associated (Property))

                if (Property%ComputeOptions%TimeSerie) then                                    
                    
                    if (Property%OutputName == 'NAME') then
                        PropVector (i) = trim(adjustl(adjustr(Property%ID%Name)))
                    else if (Property%OutputName == 'DESCRIPTION') then !if2
                        PropVector (i) = trim(adjustl(adjustr(Property%ID%Description)))
                    end if
                    i = i + 1              
                    
                end if                        
                    
                Property => Property%Next
            end do

        end if if0    
        
        nullify (Property) 
        

    end subroutine FillPropNameVector

    !--------------------------------------------------------------------------    
    
    subroutine ConstructHdf5Output

        !Arguments--------------------------------------------------------------     

        !Local-----------------------------------------------------------------
        real, pointer, dimension(:, :)              :: GridData
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: HDF5_CREATE
        integer                                     :: ObjHorizontalGrid = 0
        integer                                     :: ObjGridData       = 0
        integer                                     :: ObjBasinGeometry  = 0
        integer, dimension(:,:), pointer            :: BasinPoints       => null()
        !----------------------------------------------------------------------
       
        
        
        !Build all Topography stuff
        
        !Constructs Horizontal Grid
        call ConstructHorizontalGrid(ObjHorizontalGrid, Me%ExtVar%TopographicFile, &
                                        STAT = STAT_CALL)           
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR010'

        !Constructs GridData
        call ConstructGridData      (ObjGridData, ObjHorizontalGrid,           &
                                        FileName = Me%ExtVar%TopographicFile,  &
                                        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR020'

        !Constructs BasinGeometry
        call ConstructBasinGeometry (ObjBasinGeometry, ObjGridData,            &
                                        ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR030'

        !Gets BasinPoints
        call GetBasinPoints         (ObjBasinGeometry, BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR040'        
        
        !Gets the size of the grid
        call GetHorizontalGridSize (ObjHorizontalGrid,                               &
                                    Size     = Me%ExtVar%Size2D,                     &
                                    WorkSize = Me%ExtVar%WorkSize2D,                 &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR050'   
            
            
            

        !Bounds
        WorkILB = Me%ExtVar%WorkSize2D%ILB 
        WorkIUB = Me%ExtVar%WorkSize2D%IUB 

        WorkJLB = Me%ExtVar%WorkSize2D%JLB 
        WorkJUB = Me%ExtVar%WorkSize2D%JUB             
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%ResultsFile)//"5", &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR060'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(ObjHorizontalGrid, Me%ObjHDF5,       &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR070'

        !Gets a pointer to GridData
        call GetGridData      (ObjGridData, GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR080'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR090'

        !Writes the GridData
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",            &
                              Array2D = GridData,                                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR100'

        !Writes the WaterPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",         &
                              Array2D = BasinPoints,            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR110'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR120'

        
        
        
        !Unget all stuff
        !Ungets the GridData
        call UngetGridData (ObjGridData, GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR130'        
                    
        !Ungets BasinPoints
        call UngetBasin             (ObjBasinGeometry, BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR150'                     
        
        
        call KillBasinGeometry      (ObjBasinGeometry,   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR160'

        call KillGridData           (ObjGridData,        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR170'

        call KillHorizontalGrid     (ObjHorizontalGrid,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR180'        
        
        
    end subroutine ConstructHdf5Output

    !---------------------------------------------------------------------------    
    
    !---------------------------------------------------------------------------    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
   subroutine GetNumberOfReservoirs (ObjReservoirsID, nReservoirs, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        integer                                         :: nReservoirs
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            nReservoirs = Me%nReservoirs                             

            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

   end subroutine GetNumberOfReservoirs
    
   !---------------------------------------------------------------------------
   
    subroutine GetReservoirsNodeIDs (ObjReservoirsID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        integer, dimension(:),  pointer                 :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mReservoirs_, Me%InstanceID)

            Matrix => Me%ReservoirsNodeIDs

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetReservoirsNodeIDs
      
    !--------------------------------------------------------------------------
    
   subroutine SetReservoirsInflow   (ObjReservoirsID, ReservoirInflow, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        real, dimension(:), pointer                     :: ReservoirInflow
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: ReservoirPos
        integer                                         :: STAT_, ready_
        type(T_Reservoir), pointer                      :: CurrReservoir

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
            ReservoirPos = 0
            CurrReservoir => Me%FirstReservoir
            do while(associated(CurrReservoir))
                
                ReservoirPos = ReservoirPos + 1
                !fill global matrix
                Me%ReservoirInflows(ReservoirPos) = ReservoirInflow(ReservoirPos)
                
                !fill also individual reservoir values
                CurrReservoir%DNInflow = ReservoirInflow(ReservoirPos)
                CurrReservoir => CurrReservoir%Next
            enddo
            
            
            STAT_ = SUCCESS_ 
        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetReservoirsInflow

    !---------------------------------------------------------------------------          

    subroutine GetReservoirsOutflow (ObjReservoirsID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        real, dimension(:),  pointer                    :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mReservoirs_, Me%InstanceID)

            Matrix => Me%ReservoirOutflows

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetReservoirsOutflow
    
    !--------------------------------------------------------------------------

    subroutine SetDNConcReservoirs   (ObjReservoirsID, ConcentrationX,      &
                                                PropertyXIDNumber, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        real, dimension(:), pointer                     :: ConcentrationX
        integer                                         :: PropertyXIDNumber
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: ReservoirPos
        integer                                         :: STAT_, ready_
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_)

            if (STAT_ == SUCCESS_) then
            
                do ReservoirPos = 1, Me%nReservoirs            
            
                    PropertyX%DNInflowConc(ReservoirPos) = ConcentrationX(ReservoirPos)

                enddo           

            else
                write(*,*) 'Looking for Drainage Network Property in Reservoirs', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetDNConcReservoirs - ModuleReservoirs - ERR010'
            end if

        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetDNConcReservoirs

    !---------------------------------------------------------------------------       

                                                
    subroutine GetReservoirsConcentration(ObjReservoirsID, ConcentrationX, PropertyXIDNumber, &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: ObjReservoirsID
        real, pointer, dimension(:)                 :: ConcentrationX
        character(LEN = *), optional, intent(OUT)   :: PropertyXUnits
        integer,                      intent(IN )   :: PropertyXIDNumber
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_CALL              
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: UnitsSize
        integer                                     :: STAT_    

        !------------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mReservoirs_, Me%InstanceID) 

            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_CALL)

            if (STAT_CALL == SUCCESS_) then
                ConcentrationX => PropertyX%concentration

                if (present(PropertyXUnits)) then 
                   UnitsSize      = LEN (PropertyXUnits)
                   PropertyXUnits = PropertyX%ID%Units(1:UnitsSize)
                end if

                STAT_ = SUCCESS_
            else
                write(*,*) 'Looking for Property in Reservois', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'GetReservoirsConcentration - ModuleReservoirs - ERR010'
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetReservoirsConcentration

    !---------------------------------------------------------------------------

     subroutine CheckReservoirProperty (ObjReservoirsID,                        &
                                  PropertyID,                               &
                                  STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        integer                                         :: PropertyID
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        type (T_Property), pointer                      :: PropertyX
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then 
                write(*,*)
                write(*,*)'Could not find property', GetPropertyName(PropertyID)
                write(*,*)'in Reservoirs Module'
                stop 'CheckReservoirProperty - ModuleReservoirs - ERR010'
            endif
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine CheckReservoirProperty

    !---------------------------------------------------------------------------                                    
                                
   subroutine GetReservoirsnProperties (ObjReservoirsID, nProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        integer                                         :: nProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            nProperties = Me%nProperties                             

            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

   end subroutine GetReservoirsnProperties
    
   !---------------------------------------------------------------------------     

    subroutine GetReservoirsPropertiesIDByIdx (ObjReservoirsID, Idx, ID, Particulate, OutputName, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        integer, intent(IN)                             :: Idx
        integer, intent(OUT)                            :: ID
        logical, intent(OUT), optional                  :: Particulate
        integer, intent(OUT), optional                  :: STAT
        character (Len = StringLength), optional        :: OutputName

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, i
        type (T_Property), pointer                      :: CurrProp

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrProp => Me%FirstProperty
            do i = 1, idx - 1
                CurrProp => CurrProp%Next
            enddo

            ID          = CurrProp%ID%IDNumber
            
            if (present (Particulate)) then
!~                 Particulate = Check_Particulate_Property(CurrProp%ID%IDNumber)
                Particulate = CurrProp%ID%IsParticulate
            endif
            
            if (present (OutputName)) then
                if (CurrProp%OutputName == 'NAME') then
                    OutputName = CurrProp%ID%Name
                else if (CurrProp%OutputName == 'DESCRIPTION') then
                    OutputName = CurrProp%ID%Description
                end if
            end if

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetReservoirsPropertiesIDByIdx

    !---------------------------------------------------------------------------
                                                
    subroutine SearchProperty(PropertyX, PropertyXIDNumber, PrintWarning, STAT)

        !Arguments--------------------------------------------------------------
        type(T_Property), optional, pointer         :: PropertyX
        integer         , optional, intent (IN)     :: PropertyXIDNumber
        logical,          optional, intent (IN)     :: PrintWarning
        integer         , optional, intent (OUT)    :: STAT

        !Local------------------------------------------------------------------
        integer                                     :: STAT_ 
        
        !-----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstProperty

do2 :   do while (associated(PropertyX)) 
if5 :       if (PropertyX%ID%IDNumber==PropertyXIDNumber) then
                exit do2 
            else
                PropertyX => PropertyX%Next                 
            end if if5
        end do do2

        !A PropertyX was found
        if (associated(PropertyX)) then
            STAT_ = SUCCESS_  
        else
            if (present(PrintWarning)) then
                if (PrintWarning) write (*,*)'Property Not Found in Module Reservoirs ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if

        if (present(STAT)) STAT = STAT_

        !-----------------------------------------------------------------------

    end subroutine SearchProperty

   !----------------------------------------------------------------------------                                                
                                                
    subroutine UnGetReservoirs1D_I(ObjReservoirsID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mReservoirs_, Me%InstanceID, "UnGetReservoirs3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetReservoirs1D_I

    !--------------------------------------------------------------------------

    subroutine UnGetReservoirs1D_R8(ObjReservoirsID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjReservoirsID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mReservoirs_, Me%InstanceID,  "UnGetReservoirs3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetReservoirs1D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyReservoirs(ObjReservoirsID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjReservoirsID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, STAT_CALL, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !Gets Current Time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%ActualTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyReservoirs - ModuleReservoirs - ERR01'
            
            call GetComputeTimeStep (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyDrainageNetLocal - ModuleReservoirs - ERR02'            
            
            call InitializeVariables(.false.)
            
            !verify if level (-> volume) is imposed and update it
            call GetImposedLevel
            
            !Get discharge flow and concentration. include possible imposed outflow 
            if (Me%ComputeOptions%Discharges) then
                call ComputeReservoirDischarges
            endif
            
            if (Me%ComputeOptions%SurfaceFluxes) then
                call ComputeReservoirSurfaceFluxes
            endif                    
            
            if (Me%ComputeOptions%BottomFluxes) then
                call ComputeReservoirBottomFluxes
            endif              
            
            !Update volumes and concentrations based on computed fluxes (except outflow)
            call UpdateReservoirStateVariables
                                    
            !Compute outflows from actual volume condition
            call ComputeReservoirOutflows
            
            !Update finally the new volumes, and derived matrixes
            call UpdateVolumes
        

            if (Me%Output%HDF .or. Me%Output%TimeSerie) then
                call ModifyOutput
            endif            

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyReservoirs

    !-------------------------------------------------------------------------
    
    subroutine GetImposedLevel()
    
        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Reservoir), pointer             :: CurrReservoir
        logical                                 :: Construct            
    
        CurrReservoir => Me%FirstReservoir
        do while(associated(CurrReservoir))
            
            Construct = .false.
            call GetImposedReservoirLevel(NewReservoir = CurrReservoir, Construct = Construct)

            CurrReservoir => CurrReservoir%Next
        enddo
    
    end subroutine GetImposedLevel
    
    !-------------------------------------------------------------------------
    
    subroutine ComputeReservoirDischarges ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer              :: Property
        type (T_Reservoir), pointer             :: CurrReservoir
        integer                                 :: iDis, ReservoirID        
        integer                                 :: STAT_CALL
        integer                                 :: iProp
        logical                                 :: IsImposedOutflow, Found

        !Actualize Volumes
        

        do iDis = 1, Me%nDischarges                  

            if (Me%DischargesActive(iDis)) then            
            
                ReservoirID = Me%ReservoirDischargeLink(iDis)

                call FindReservoir(ReservoirID, CurrReservoir, Found)
            
                call GetDischargeWaterFlow(Me%ObjDischarges,                                &
                                        Me%ExtVar%ActualTIme, iDis,                         &
                                        CurrReservoir%WaterLevel,                           &
                                        Me%ReservoirDischargeFlow(iDis), STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_) stop 'ModuleReservoir - ComputeReservoirDischarges - ERR04'
            
                !if imposed outflow, can only be negative value
                if (CurrReservoir%Management%ImposedOutflow) then
            
                    !this discharge is one imposed outflow discharge?
                    call GetIsReservoirOutflow(Me%ObjDischarges, iDis, IsImposedOutflow, STAT = STAT_CALL)
                    if (STAT_CALL/=SUCCESS_) stop 'ModuleReservoir - ConstructDischarges - ERR06'
                
                    if (IsImposedOutflow .and. Me%ReservoirDischargeFlow(iDis) > 0.0) then   
                        write(*,*) 'Found positive discharge value '
                        write(*,*) 'in Reservoir ID: ', CurrReservoir%ID
                        write(*,*) 'while imposing outflow discharge trough keyword IS_OUTFLOW'
                        write(*,*) 'in Discharge data file. Can onl have zero or negative values.'
                        stop 'ModuleReservoir - ConstructDischarges - ERR010'
                    endif
                endif
            
                nullify (Property)
                Property => Me%FirstProperty
                iProp = 0
                do while (associated (Property))
                
                    if (Property%ComputeOptions%Discharges) then 
                    
                        iProp = iProp + 1

                        call GetDischargeConcentration (Me%ObjDischarges,                              &
                                                        Me%ExtVar%ActualTime,                          &
                                                        iDis, Me%ReservoirDischargeConc(iDis, iProp),  &
                                                        Property%ID%IDNumber,                          &
                                                        STAT = STAT_CALL)
                        if (STAT_CALL/=SUCCESS_) then
                            if (STAT_CALL == NOT_FOUND_ERR_) then 
                                !When a property is not found associated to a discharge
                                !by default is consider that the concentration is zero
                                Me%ReservoirDischargeConc(iDis, iProp) = 0.
                            else
                                stop 'ModuleReservoir - ComputeReservoirDischarges - ERR020'
                            endif
                        endif                
                    
                    end if
                                
                    Property => Property%Next

                enddo                
            endif
        enddo

    end subroutine ComputeReservoirDischarges
    
    !---------------------------------------------------------------------------    
    
   subroutine ComputeReservoirSurfaceFluxes ()


        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------

        !Begin------------------------------------------------------------------
        
        

   end subroutine ComputeReservoirSurfaceFluxes
    
    !---------------------------------------------------------------------------    
   
   subroutine ComputeReservoirBottomFluxes ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------

        !Begin------------------------------------------------------------------
  



    end subroutine ComputeReservoirBottomFluxes
    
    !---------------------------------------------------------------------------       
    
    subroutine UpdateReservoirStateVariables ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer              :: Property
        type (T_Reservoir), pointer             :: CurrReservoir
        integer                                 :: iDis, ReservoirID        
        real                                    :: VolumeNew    
        integer                                 :: STAT_CALL
        integer                                 :: iProp
        logical                                 :: IsImposedOutflow, Found
        real                                    :: Concentration, ConcDif, RetentionTime
        logical                                 :: ProcessDischargeHere     = .true.
        
       !Actualize Volumes and concentration due to DN inflow, computed outflow and surface fluxes    
        CurrReservoir => Me%FirstReservoir
        do while (associated (CurrReservoir))        
            
            !DNInflow Conc
            nullify (Property)
            Property => Me%FirstProperty
            do while (associated (Property))
                                
                Concentration = DischargeProperty (Me%ReservoirInflows(CurrReservoir%Position),       &
                                                     Property%DNInflowConc(CurrReservoir%Position),   &
                                                    CurrReservoir,  Property, Property%IScoefficient, .false.)                
                    
                !Instant mixing
                if (Me%ReservoirInflows(CurrReservoir%Position) .le. 0.0     &
                    .or. Me%ComputeOptions%PropertyComputeMethod == InstantMixing_) then
                    
                    Property%Concentration(CurrReservoir%Position) = Concentration
                
                !complete mixing only ocurs at retention time
                else if (Me%ComputeOptions%PropertyComputeMethod == RetentionTimeMixing_) then                       
                        
                    ConcDif = Concentration - Property%Concentration(CurrReservoir%Position)
                        
                    ![s] = [m3] / [m3/s]
                    RetentionTime = CurrReservoir%VolumeOld / Me%ReservoirInflows(CurrReservoir%Position)
                        
                    !computed concentration would only occur after retention time (volume / flow)
                    if (RetentionTime < Me%ExtVar%DT) then
                        Property%Concentration(CurrReservoir%Position) = Concentration
                            
                    else                            
                        Property%Concentration(CurrReservoir%Position) = Property%Concentration(CurrReservoir%Position)   &
                                                                          + (ConcDif * Me%ExtVar%DT / RetentionTime)
                    endif                               
                        
                endif  
                                
                Property => Property%Next

            enddo                
            !DNInflow Flows - always positive (computed from volume in node)
            CurrReservoir%VolumeNew = CurrReservoir%VolumeNew + Me%ReservoirInflows(CurrReservoir%Position) * Me%ExtVar%DT
            
            
            
            !TODO : Surface fluxes concentration change            
            !surface fluxes
            VolumeNew = CurrReservoir%VolumeNew + CurrReservoir%SurfaceFluxes * Me%ExtVar%DT
            !only remove what is available
            if (CurrReservoir%SurfaceFluxes .lt. 0.0 .and. VolumeNew .lt. 0.0) then
                !       m3/s            =       m3    /    s
                CurrReservoir%SurfaceFluxes = - CurrReservoir%VolumeNew / Me%ExtVar%DT
                VolumeNew = 0.0
            endif               
            CurrReservoir%VolumeNew = VolumeNew        

            
            CurrReservoir => CurrReservoir%Next
        enddo        
        
        
        !Actualize Volumes and concentration due to reservoir discharges (do not include imposed outflow)
        do iDis = 1, Me%nDischarges                   
            
            !avoid DN discharges tath can be present at the same time
            if (Me%DischargesActive(iDis)) then     
                
                ReservoirID = Me%ReservoirDischargeLink(iDis)

                call FindReservoir(ReservoirID, CurrReservoir, Found)
            
                ProcessDischargeHere = .true.
                !verify if outflow (save value and do not process here)
                if (CurrReservoir%Management%ImposedOutflow) then
            
                    !this discharge is one imposed outflow discharge?
                    call GetIsReservoirOutflow(Me%ObjDischarges, iDis, IsImposedOutflow, STAT = STAT_CALL)
                    if (STAT_CALL/=SUCCESS_) stop 'ModuleReservoir - UpdateReservoirStateVariables - ERR06'
                
                    !save value and do not process here but in outflow routine
                    if (IsImposedOutflow .and. Me%ReservoirDischargeFlow(iDis) .lt. 0.0) then
                        
                        !do not account new volume here. it will be accounted in routine compute outflow
                        ProcessDischargeHere = .false.
                        
                        !Outflow is positive
                        CurrReservoir%Outflow = CurrReservoir%Outflow - Me%ReservoirDischargeFlow(iDis)                        
                    
                        !!update outflows to DN
                        !Me%ReservoirOutflows(CurrReservoir%Position) = Me%ReservoirOutflows(CurrReservoir%Position)   &
                        !                                                - Me%ReservoirDischargeFlow(iDis)                    
                    endif
            
                endif                
                
                !outflows are not processed here but in routine of outflow c0mputation. they do not chane conc
                if (ProcessDischargeHere) then
                
                    !Update of volume can only be after DischargeProperty (because math is made with old volume + discharge volume)
                    VolumeNew = CurrReservoir%VolumeNew + Me%ReservoirDischargeFlow(iDis) * Me%ExtVar%DT
            
                    !only remove what is available
                    if (Me%ReservoirDischargeFlow(iDis).lt. 0.0 .and. VolumeNew .lt. 0.0) then
                        !       m3/s            =       m3    /    s
                        Me%ReservoirDischargeFlow(iDis) = - CurrReservoir%VolumeNew / Me%ExtVar%DT
                        VolumeNew = 0.0
                    endif                    
                        
                    !Just to acount in sum matrixes (positive and negative)
                    CurrReservoir%Discharges = CurrReservoir%Discharges + Me%ReservoirDischargeFlow(iDis)
            
            
                    nullify (Property)
                    Property => Me%FirstProperty
                    iProp = 0
                    do while (associated (Property))
                
                        if (Property%ComputeOptions%Discharges) then 
                    
                            iProp = iProp + 1
                    
                            Concentration = DischargeProperty (Me%ReservoirDischargeFlow(iDis), & 
                                                               Me%ReservoirDischargeConc(iDis, iProp), &
                                                               CurrReservoir,  Property, Property%IScoefficient, .false.) 
                    
                            !instant mixing in case of positive discharge, if neative it does not matter
                            if (Me%ReservoirDischargeFlow(iDis) .le. 0.0     &
                                 .or. Me%ComputeOptions%PropertyComputeMethod == InstantMixing_) then
                    
                                Property%Concentration(CurrReservoir%Position) = Concentration
                    
                            !complete mixing would only occur at retention time
                            else if (Me%ComputeOptions%PropertyComputeMethod == RetentionTimeMixing_) then                       
                        
                                ConcDif = Concentration - Property%Concentration(CurrReservoir%Position)
                        
                                ![s] = [m3] / [m3/s]
                                RetentionTime = CurrReservoir%VolumeOld / Me%ReservoirDischargeFlow(iDis)
                        
                                !computed concentration would only occur after retention time (volume / flow)
                                if (RetentionTime < Me%ExtVar%DT) then
                                    Property%Concentration(CurrReservoir%Position) = Concentration
                            
                                else                            
                                    Property%Concentration(CurrReservoir%Position) = &
                                        Property%Concentration(CurrReservoir%Position) &
                                        + (ConcDif * Me%ExtVar%DT / RetentionTime)
                                endif                                                           
                        
                            endif
                    
                        end if
                                
                        Property => Property%Next

                    enddo
            
                    !only after all propety updates that used old volume
                    CurrReservoir%VolumeNew = VolumeNew 
                    
                endif
                
            endif
        enddo        
        
        
        !!save global matrixes values   
        !CurrReservoir => Me%FirstReservoir
        !do while (associated (CurrReservoir))                    
        !
        !    
        !    call UpdateMatrixesFromVolume (CurrReservoir)
        !
        !    
        !    CurrReservoir => CurrReservoir%Next
        !enddo        

    end subroutine UpdateReservoirStateVariables
    
    !---------------------------------------------------------------------------    
            
    subroutine UpdateVolumes()
    
        !Local------------------------------------------------------------------
        type (T_Reservoir), pointer                 :: CurrReservoir

        
        
        CurrReservoir => Me%FirstReservoir
        do while (associated (CurrReservoir))                    
                  
        
            !Update Volume
            CurrReservoir%VolumeNew = CurrReservoir%VolumeNew - CurrReservoir%Outflow * Me%ExtVar%DT
            
            !volumes close to zero. use abs to not masquerade possible erros with negative volumes
            if (Abs(CurrReservoir%VolumeNew) < AlmostZero) then
                CurrReservoir%VolumeNew = 0.0
            endif        
            
            !Update water level at the end
            if (associated(CurrReservoir%Management%AccVolumeCurve)) then            
                CurrReservoir%WaterLevel = ComputeReservoirLevel(CurrReservoir)
            endif
            
            !Save in global matrix
            Me%ReservoirVolumes(CurrReservoir%Position)  = CurrReservoir%VolumeNew           
            Me%ReservoirPercFull(CurrReservoir%Position) = CurrReservoir%VolumeNew / CurrReservoir%MaxVolume * 100
            
            CurrReservoir%PercFull = CurrReservoir%VolumeNew / CurrReservoir%MaxVolume * 100    
        
        
            CurrReservoir => CurrReservoir%Next
        enddo         
    
    end subroutine UpdateVolumes
    
    !----------------------------------------------------------------------------
    
    real function DischargeProperty (DischargeFlow, DischargeConc, Reservoir,         &
                                   Property, ISCoef, Accumulate)
        !Arguments--------------------------------------------------------------
        real                                        :: DischargeFlow, DischargeConc
        type (T_Property), pointer                  :: Property
        type (T_Reservoir), pointer                 :: Reservoir
        real                                        :: ISCoef !, Concentration
        logical                                     :: Accumulate

        !Local------------------------------------------------------------------
        real                                        :: DischargeVolume
        real                                        :: OldMass, NewMass
        real                                        :: ISDischargeConc, ISConcentration

        
        ISDischargeConc = DischargeConc * ISCoef
        ISConcentration = Property%Concentration(Reservoir%Position) * ISCoef

        if (abs(DischargeFlow) > AllmostZero) then
            
            ![m3] = [s] * [m3/s]
            DischargeVolume  = dble(Me%ExtVar%DT)*dble(DischargeFlow)
            
            ![g] = [g/m3] * [m3]
            OldMass          = dble(ISConcentration) * Reservoir%VolumeNew            
        
            if      (DischargeFlow > 0.0) then

                !Explicit discharges input 
                NewMass          = OldMass + DischargeVolume * dble(ISDischargeConc)                                       

                ISConcentration = NewMass / (Reservoir%VolumeNew + DischargeFlow * Me%ExtVar%DT)

            elseif (DischargeFlow < 0.0 .and. Reservoir%VolumeNew > 0.0) then
                    
                !If the discharge flow is negative (Output) then the concentration
                !to consider is the concentration of the reservoir where the discharge
                !is located

                !If the property acculumlates in the water column 
                !(e.g particulate properties during removal) then the concentration will increase

                if (Accumulate) then
                    NewMass          = OldMass
                else
                    NewMass          = OldMass * (1.0 + DischargeVolume / Reservoir%VolumeNew)
                endif
                
                !if water remains
                if (abs(DischargeVolume) < Reservoir%VolumeNew) then
                   
                    ISConcentration    = NewMass / (Reservoir%VolumeNew + DischargeVolume)
                
                else   !if all water exits reservoi than accumulated mass needs to be accounted in bottom!
                    
                    ISConcentration  = 0.0
                    
                    if (Accumulate) then
                        ![kg/m2] = [kg/m2] + [g] * 1e-3 [kg/g] / m2
                        Property%BottomConc(Reservoir%Position) = Property%BottomConc(Reservoir%Position) +            &
                                                      (NewMass * 1e-3 /                        &
                                                       (Reservoir%BottomArea))
                    endif
                endif
            endif

        else
            
            !Do Nothing            

        endif
                    
        
        DischargeProperty = ISConcentration / ISCoef


    end function DischargeProperty

    !---------------------------------------------------------------------------    
    
    subroutine ComputeReservoirOutflows()

        !Arguments-----------------------------------------------------------_

        !Local-----------------------------------------------------------------
        type(T_Reservoir),           pointer    :: CurrentReservoir
        real                                    :: Outflow
		logical	                                :: GetOutflowFromOperation_
        real, dimension(6)  , target            :: AuxTime
        integer                                 :: NextMonth
        integer                                 :: SecondsUntilNextMonth
        real                                    :: CurrReservoir_MaxVol, CurrReservoir_NextMaxVol
        !Begin----------------------------------------------------------------

		! To get the year, month, day, hours, minutes and seconds of the date - Ana Oliveira
        call ExtractDate   (Me%ExtVar%ActualTime , AuxTime(1), AuxTime(2),         &
                                                   AuxTime(3), AuxTime(4),         &
                                                   AuxTime(5), AuxTime(6))
		
		if (int(AuxTime(2)) == 12) then
		 NextMonth = 1
		else
		 NextMonth = int(AuxTime(2)) + 1
		end if
        
        SecondsUntilNextMonth = NumberOfDaysInMonth(int(AuxTime(1)), int(AuxTime(2))) * 86400 - int(AuxTime(3)) * &
						    86400 + int(AuxTime(4)) * 3600 + int(AuxTime(5)) * 60 + int(AuxTime(6)) + 86400
        
        CurrentReservoir => Me%FirstReservoir
        
        do while (associated(CurrentReservoir))
            CurrReservoir_MaxVol = CurrentReservoir%MaxVolume_Monthly(int(AuxTime(2)))
            CurrReservoir_NextMaxVol = CurrentReservoir%MaxVolume_Monthly(NextMonth)
            
            ! If user defines one max volume per month a correction between consecutive months needs to be performed - Ana Oliveira
			if (CurrentReservoir%MaxVolumeType == MonthlyValues_) then
                if (CurrentReservoir%VolumeNew > CurrReservoir_NextMaxVol .and. &
                    CurrReservoir_NextMaxVol < CurrReservoir_MaxVol) then

                        CurrentReservoir%ExtraOutflow = (CurrReservoir_MaxVol - CurrReservoir_NextMaxVol) / SecondsUntilNextMonth
							
						GetOutflowFromOperation_ = .true.
                        
                else if (CurrentReservoir%VolumeNew < CurrentReservoir%MaxVolume .and. &
                            CurrentReservoir%VolumeNew > CurrReservoir_NextMaxVol .and. &
                            CurrReservoir_NextMaxVol < CurrReservoir_MaxVol) then

                            CurrentReservoir%ExtraOutflow = (CurrentReservoir%VolumeNew - CurrReservoir_NextMaxVol) / SecondsUntilNextMonth + Me%ReservoirInflows(CurrentReservoir%Position)
                                
                            GetOutflowFromOperation_ = .false.
                        
                else
                    CurrentReservoir%ExtraOutflow = 0.0
                endif
            else
                CurrentReservoir%ExtraOutflow = 0.0
            endif            
            
            !since reservoirs never go below (associated to discharger location) limit. 
            !in undefined MinVolume is zero
            if (CurrentReservoir%VolumeNew <= CurrentReservoir%MinVolume) then
                Outflow = 0.0 
         
            !if level is imposed
            else if (CurrentReservoir%Management%ON .and. CurrentReservoir%Management%ImposedLevel) then
                
                !outflow is the remainder balance betwewn all the fluxes already computed (VolumeNew) and 
                !the VolumeTarget read from file
                !negative outlflow (underpredicted inflows)
                if (CurrentReservoir%VolumeTarget >= CurrentReservoir%VolumeNew) then
                    Outflow = 0.0
                else                    
                    Outflow = (CurrentReservoir%VolumeNew - CurrentReservoir%VolumeTarget) / Me%ExtVar%DT
                endif            
            
            !if imposed by discharge verify computation
            else if (CurrentReservoir%Management%ON .and. CurrentReservoir%Management%ImposedOutflow) then  
                
                !computed already in discharges
                Outflow = CurrentReservoir%Outflow
                                  
            !operation
            else if (CurrentReservoir%Management%ON .and. CurrentReservoir%Management%ImposedOperation) then            
                                    
                if (GetOutflowFromOperation_) then
					!Compute outflow. uses environmental flow if curves defined above curr volume
					Outflow = GetOutflowFromOperation(CurrentReservoir)
				end if                                                     
                                
            !only minimum flow defined
            else if (CurrentReservoir%Management%ON) then
                
                Outflow = CurrentReservoir%Management%MinOutflow
                 
            !First verify if not managed 
            else if (.not. CurrentReservoir%Management%ON) then
                
                !type weir
                if (CurrentReservoir%IsWeir) then                   
                    Outflow = ComputeWeirFlow(CurrentReservoir)
                else                    
                    !type ditch, flow zero until maximum capacity, then all exits
                    Outflow = 0.0
                endif                        
                          
            endif
            
            
            !avoid under min volume and above max volume after outflow integration
            if (CurrentReservoir%VolumeNew - (Outflow * Me%ExtVar%DT) < CurrentReservoir%MinVolume) then
                Outflow = (CurrentReservoir%VolumeNew - CurrentReservoir%MinVolume) / Me%ExtVar%DT   
            else if (CurrentReservoir%VolumeNew - (Outflow * Me%ExtVar%DT) > CurrentReservoir%MaxVolume) then
                ![m3/s] = [m3] / [s]
                Outflow = (CurrentReservoir%VolumeNew - CurrentReservoir%MaxVolume) / Me%ExtVar%DT            
            endif

            !limit to available water
            Outflow = min(Outflow, CurrentReservoir%VolumeNew / Me%ExtVar%DT)
            
            !limit to positive and project flow
            Outflow = max(Outflow, 0.0)            
            if (Outflow > CurrentReservoir%Management%MaxOutflow) then
                Outflow = CurrentReservoir%Management%MaxOutflow
            endif
            
            !If the year of construction is bigger than the current year the outflow will be calculate 
            !by the operational curve, otherwise the outflow will bu equal the inflow and there is not 
            !accumulated volume - Ana Oliveira
            if (AuxTime(1) >= CurrentReservoir%ConstructionYear) then
                CurrentReservoir%Outflow = Outflow
            else
                CurrentReservoir%Outflow = Me%ReservoirInflows(CurrentReservoir%Position)
            endif

             !update outflows to DN - includes the ones from imposed discharge
            CurrentReservoir%Outflow = CurrentReservoir%Outflow + CurrentReservoir%ExtraOutflow
            Me%ReservoirOutflows(CurrentReservoir%Position) = Me%ReservoirOutflows(CurrentReservoir%Position)   &
                                                                + CurrentReservoir%Outflow + CurrentReservoir%ExtraOutflow
                                                                   
            CurrentReservoir => CurrentReservoir%Next
        enddo
                    

    end subroutine ComputeReservoirOutflows

    !--------------------------------------------------------------------------        

    real function ComputeWeirFlow(CurrentReservoir)
        !Arguments------------------------------------------------------------
        type(T_Reservoir),           pointer     :: CurrentReservoir

        !Local-----------------------------------------------------------------
        real                                     :: HeightAboveWeir
        real                                     :: Flow
        !---------------------------------------------------------------------
        
        !Q = cv * b * sqrt(2*g) * H^(1.5)            
        !Upstream of the Weir
        
        if (CurrentReservoir%WaterLevel /= null_real) then
            HeightAboveWeir = CurrentReservoir%WaterLevel  - CurrentReservoir%FlowOver%CrestLevel
            
            if (HeightAboveWeir > 0) then
                Flow = -sqrt(19.6) * CurrentReservoir%FlowOver%DischargeCoeficient *       &
                                        CurrentReservoir%FlowOver%WeirLength  * HeightAboveWeir ** 1.5
            else
                Flow = 0.
            endif
        else
            write(*,*) 'Can not define flow on reservoir ID : ',  CurrentReservoir%ID
            write(*,*) 'Because level is not defined'
            stop 'ComputeWeirFlow - ModuleResevoirs - ERR01'
        endif
        
        ComputeWeirFlow = Flow
    
    end function ComputeWeirFlow
    
    !-------------------------------------------------------------------------
    
    real function GetOutflowFromOperation(CurrentReservoir)
        !Arguments------------------------------------------------------------
        type(T_Reservoir),           pointer     :: CurrentReservoir

        !Local-----------------------------------------------------------------
        integer                                 :: i
        real                                    :: Outflow, ReservoirPercentageVolume
        real                                    :: PreviousCurvePercentageVolume, PreviousCurveOutflow, PreviousCurvePercMaxOutflow
        real                                    :: NextCurvePercentageVolume, NextCurveOutflow, NextCurvePercMaxOutflow
        real                                    :: ReservoirLevel, PercentInflow
        real                                    :: PreviousCurvePercInflow, NextCurvePercInflow
        real                                    :: PreviousCurveLevel, NextCurveLevel
        !Begin----------------------------------------------------------------
        
        Outflow = 0.0
        
        if (CurrentReservoir%Management%OperationType == Operation_PercVol_Outflow_) then
            
            ReservoirPercentageVolume     = CurrentReservoir%VolumeNew / CurrentReservoir%MaxVolume
            PreviousCurvePercentageVolume = CurrentReservoir%Management%OperationCurve(1, 1)
            PreviousCurveOutflow          = CurrentReservoir%Management%OperationCurve(1, 2)
            
            !reservoir empty. avoid all kinds of interpolations
            if (CurrentReservoir%VolumeNew < AllmostZero) then
                
                Outflow = 0.0
                
            !reservoir lower than first point
            else if (ReservoirPercentageVolume < PreviousCurvePercentageVolume) then
                
                !set environmental flow (if not defined is zero)
                Outflow = CurrentReservoir%Management%MinOutflow 
                    
            !reservoir higherthan last point
            else if (ReservoirPercentageVolume >     &
               CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 1)) then
                
                !set last ouflow
                Outflow = CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 2)
                
            else
                
                !go trough all points to find where belongs
do1:             do i = 2, CurrentReservoir%Management%OperationCurvePoints
                
                    NextCurvePercentageVolume = CurrentReservoir%Management%OperationCurve(i, 1)
                    NextCurveOutflow          = CurrentReservoir%Management%OperationCurve(i, 2)
                
                    if (ReservoirPercentageVolume >= PreviousCurvePercentageVolume &
                        .and. ReservoirPercentageVolume <= NextCurvePercentageVolume) then
                    
                        Outflow = LinearInterpolation(PreviousCurvePercentageVolume, PreviousCurveOutflow,       &
                                         NextCurvePercentageVolume, NextCurveOutflow, ReservoirPercentageVolume)   
                        exit do1
                        
                    endif
                
                    PreviousCurvePercentageVolume = NextCurvePercentageVolume
                    PreviousCurveOutflow          = NextCurveOutflow
                enddo do1
                
            endif
            
        else if (CurrentReservoir%Management%OperationType == Operation_PercVol_PercMaxOutflow_) then
            
            ReservoirPercentageVolume     = CurrentReservoir%VolumeNew / CurrentReservoir%MaxVolume
            PreviousCurvePercentageVolume = CurrentReservoir%Management%OperationCurve(1, 1)
            PreviousCurvePercMaxOutflow   = CurrentReservoir%Management%OperationCurve(1, 2)
            
            !reservoir empty. avoid all kinds of interpolations
            if (CurrentReservoir%VolumeNew < AllmostZero) then
                
                Outflow = 0.0
                
            !reservoir lower than first point
            else if (ReservoirPercentageVolume < PreviousCurvePercentageVolume) then
                
                !set environmental flow (if not defined is zero)
                Outflow = CurrentReservoir%Management%MinOutflow 
                    
            !reservoir higherthan last point
            else if (ReservoirPercentageVolume >     &
               CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 1)) then
                
                !set last ouflow
                Outflow = CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 2)
                
            else
                
                !go trough all points to find where belongs
do12:           do i = 2, CurrentReservoir%Management%OperationCurvePoints
                
                    NextCurvePercentageVolume = CurrentReservoir%Management%OperationCurve(i, 1)
                    NextCurvePercMaxOutflow   = CurrentReservoir%Management%OperationCurve(i, 2)
                
                    if (ReservoirPercentageVolume >= PreviousCurvePercentageVolume &
                        .and. ReservoirPercentageVolume <= NextCurvePercentageVolume) then
                    
                        Outflow = LinearInterpolation(PreviousCurvePercentageVolume, PreviousCurvePercMaxOutflow,       &
                                         NextCurvePercentageVolume, NextCurvePercMaxOutflow, ReservoirPercentageVolume)  &
                                    * CurrentReservoir%Management%MaxOutflow
                        exit do12
                        
                    endif
                
                    PreviousCurvePercentageVolume = NextCurvePercentageVolume
                    PreviousCurvePercMaxOutflow   = NextCurvePercMaxOutflow
                enddo do12
                
            endif            
                    
        elseif (CurrentReservoir%Management%OperationType == Operation_PercVol_PercInflow_) then
            
            ReservoirPercentageVolume     = CurrentReservoir%VolumeNew / CurrentReservoir%MaxVolume
            PreviousCurvePercentageVolume = CurrentReservoir%Management%OperationCurve(1, 1)
            PreviousCurvePercInflow       = 0.0
            
           !reservoir empty
            if (CurrentReservoir%VolumeNew < AllmostZero) then
                
                Outflow = 0.0
                
            !reservoir lower than first point
            else if (ReservoirPercentageVolume < PreviousCurvePercentageVolume) then
                
                !set environmental flow (if not defined is zero)
                Outflow = CurrentReservoir%Management%MinOutflow 
                    
            !reservoir higherthan last point
            else if (ReservoirPercentageVolume >      &
                  CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 1)) then
                
                !set last ouflow
                Outflow = CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 2)  &
                           * Me%ReservoirInflows(CurrentReservoir%Position)
                
            else            
            
                !go trough all points to find where belongs
do2:             do i = 2, CurrentReservoir%Management%OperationCurvePoints
                
                    NextCurvePercentageVolume = CurrentReservoir%Management%OperationCurve(i, 1)
                    NextCurvePercInflow       = CurrentReservoir%Management%OperationCurve(i, 2)
                
                    if (ReservoirPercentageVolume >= PreviousCurvePercentageVolume &
                        .and. ReservoirPercentageVolume <= NextCurvePercentageVolume) then
                    
                        PercentInflow = LinearInterpolation(PreviousCurvePercentageVolume, PreviousCurvePercInflow,       &
                                         NextCurvePercentageVolume, NextCurvePercInflow, ReservoirPercentageVolume)       
                              
                        Outflow = PercentInflow * Me%ReservoirInflows(CurrentReservoir%Position)
                        
                        exit do2
                        
                    endif
                
                    PreviousCurvePercentageVolume = NextCurvePercentageVolume
                    PreviousCurvePercInflow       = NextCurvePercInflow
                    
                enddo do2
                
            endif
            
        elseif (CurrentReservoir%Management%OperationType == Operation_Level_Outflow_) then
            
            ReservoirLevel                = ComputeReservoirLevel(CurrentReservoir)
            PreviousCurveLevel            = CurrentReservoir%Management%OperationCurve(1, 1)
            PreviousCurveOutflow          = CurrentReservoir%Management%OperationCurve(i, 2)
            
            !reservoir empty
            if (CurrentReservoir%VolumeNew < AllmostZero) then
                
                Outflow = 0.0
                
            !reservoir lower than first point
            else if (ReservoirLevel < PreviousCurveLevel) then
                
                !set environmental flow (if not defined is zero)
                Outflow = CurrentReservoir%Management%MinOutflow 
                    
            !reservoir higherthan last point
            else if (ReservoirLevel >      &
                CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 1)) then
                
                !set last ouflow
                Outflow = CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 2)
                
            else                        
            
                !go trough all points to find where belongs
do3:            do i = 2, CurrentReservoir%Management%OperationCurvePoints
                
                    NextCurveLevel            = CurrentReservoir%Management%OperationCurve(i, 1)
                    NextCurveOutflow          = CurrentReservoir%Management%OperationCurve(i, 2)
                
                    if (ReservoirLevel >= PreviousCurveLevel &
                        .and. ReservoirLevel <= NextCurveLevel) then
                    
                        Outflow = LinearInterpolation(PreviousCurveLevel, PreviousCurveOutflow,       &
                                         NextCurveLevel, NextCurveOutflow, ReservoirLevel)                             
                        
                        exit do3
                        
                    endif
                
                    PreviousCurveLevel    = NextCurveLevel
                    PreviousCurveOutflow  = NextCurveOutflow
                enddo do3
                
            endif
                    
        elseif (CurrentReservoir%Management%OperationType == Operation_Level_PercInflow_) then
            
            ReservoirLevel                = ComputeReservoirLevel(CurrentReservoir)
            PreviousCurveLevel            = CurrentReservoir%Management%OperationCurve(1, 1)
            PreviousCurvePercInflow       = CurrentReservoir%Management%OperationCurve(i, 2)
            
           !reservoir empty
            if (CurrentReservoir%VolumeNew < AllmostZero) then
                
                Outflow = 0.0
                
            !reservoir lower than first point
            else if (ReservoirLevel < PreviousCurveLevel) then
                
                !set environmental flow (if not defined is zero)
                Outflow = CurrentReservoir%Management%MinOutflow 
                    
            !reservoir higherthan last point
            else if (ReservoirLevel >       &
                CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 1)) then
                
                !set last ouflow
                Outflow = CurrentReservoir%Management%OperationCurve(CurrentReservoir%Management%OperationCurvePoints, 2)  &
                           * Me%ReservoirInflows(CurrentReservoir%Position)
                
            else            
                
            
                !go trough all points to find where belongs
do4:            do i = 2, CurrentReservoir%Management%OperationCurvePoints
                
                    NextCurveLevel            = CurrentReservoir%Management%OperationCurve(i, 1)
                    NextCurvePercInflow       = CurrentReservoir%Management%OperationCurve(i, 2)
                
                    if (ReservoirLevel >= PreviousCurveLevel &
                        .and. ReservoirLevel <= NextCurveLevel) then
                    
                        PercentInflow = LinearInterpolation(PreviousCurveLevel, PreviousCurvePercInflow,       &
                                         NextCurveLevel, NextCurvePercInflow, ReservoirLevel)       
                              
                        Outflow = PercentInflow * Me%ReservoirInflows(CurrentReservoir%Position)
                        
                        exit do4
                        
                    endif
                
                    PreviousCurveLevel            = NextCurveLevel
                    PreviousCurvePercInflow       = NextCurvePercInflow
                enddo do4
                   
            endif
            
        endif
        
        GetOutflowFromOperation = Outflow
    
    end function GetOutflowFromOperation
    
    !---------------------------------------------------------------------------
        
        
    subroutine ModifyOutput ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        logical                                     :: IsFinalFile
  
        if (Me%Output%TimeSerie) then            
            call WriteTimeSeries            
        endif
        
        if (Me%Output%HDF .and. Me%Output%HDFActive) then
            call WriteOutputHDF
        endif
        
        !Restart Output
        if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%ActualTime == Me%ExtVar%EndTime)) then
            if(Me%ExtVar%ActualTime >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                IsFinalFile = .false.
                call WriteFinalHDF(IsFinalFile)
                Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
            endif
        endif        
        
    end subroutine ModifyOutput
    
   !---------------------------------------------------------------------------
   
    subroutine WriteTimeSeries ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Reservoir    ), pointer                     :: CurrReservoir
        type (T_Property     ), pointer                     :: Property
        integer                                             :: i, j        
        integer                                             :: STAT_CALL
        

        i = 0
        j = 7 !base outputs
        CurrReservoir => Me%FirstReservoir
do1:    do while(associated(CurrReservoir))            
            
            i = i + 1
            
if1:        if (CurrReservoir%TimeSerie) then        
                
                Me%TimeSerie%DataLine (1) = CurrReservoir%VolumeNew
                Me%TimeSerie%DataLine (2) = CurrReservoir%PercFull
                Me%TimeSerie%DataLine (3) = CurrReservoir%WaterLevel
                Me%TimeSerie%DataLine (4) = CurrReservoir%DNInflow
                Me%TimeSerie%DataLine (5) = CurrReservoir%Discharges
                Me%TimeSerie%DataLine (6) = CurrReservoir%SurfaceFluxes
                Me%TimeSerie%DataLine (7) = CurrReservoir%Outflow


if2:            if (Me%HasProperties) then

                    Property => Me%FirstProperty
                    j = j + 1
                                        
                    do while (associated (Property))                

                        if (Property%ComputeOptions%TimeSerie) then
                            Me%TimeSerie%DataLine (j) = Property%Concentration(CurrReservoir%Position)   
                            j = j + 1

                        end if
                
                        Property => Property%Next

                    end do

                end if if2
                                
                call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),           &
                                        DataLine  = Me%TimeSerie%DataLine,      &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByNodes - ERR01'
            
            end if if1 

            CurrReservoir => CurrReservoir%Next
        
        end do do1
    
    end subroutine WriteTimeSeries
    
    !-------------------------------------------------------------------------
    
    subroutine WriteOutputHDF
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        real, dimension(6)  , target                :: AuxTime
        real, dimension(:)  , pointer               :: TimePointer       
        real, dimension(:,:), pointer               :: ReservoirVolume, ReservoirPercFull
        real, dimension(:,:), pointer               :: ReservoirDNInflow, ReservoirDischarges
        real, dimension(:,:), pointer               :: ReservoirSurfaceFluxes, ReservoirOutflow
        real, dimension(:,:), pointer               :: ReservoirConc, ReservoirBottomConc
        type(T_Reservoir), pointer                  :: CurrReservoir
        type(T_Property), pointer                   :: CurrProperty
        
        !Bounds
        ILB = Me%ExtVar%WorkSize2D%ILB
        IUB = Me%ExtVar%WorkSize2D%IUB

        JLB = Me%ExtVar%WorkSize2D%JLB
        JUB = Me%ExtVar%WorkSize2D%JUB


        if (Me%ExtVar%ActualTime >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

            !Writes current time
            call ExtractDate   (Me%ExtVar%ActualTime , AuxTime(1), AuxTime(2),         &
                                                AuxTime(3), AuxTime(4),                &
                                                AuxTime(5), AuxTime(6))
            TimePointer => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR01'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                 "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D      = TimePointer,                    &
                                 OutputNumber = Me%OutPut%NextOutPut,           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR02'

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR03'
            
            allocate(ReservoirVolume         (ILB:IUB, JLB:JUB))
            allocate(ReservoirPercFull       (ILB:IUB, JLB:JUB))
            allocate(ReservoirDNInflow       (ILB:IUB, JLB:JUB))
            allocate(ReservoirDischarges     (ILB:IUB, JLB:JUB))
            allocate(ReservoirSurfaceFluxes  (ILB:IUB, JLB:JUB))
            allocate(ReservoirOutflow        (ILB:IUB, JLB:JUB))
            allocate(ReservoirConc           (ILB:IUB, JLB:JUB))
            allocate(ReservoirBottomConc     (ILB:IUB, JLB:JUB))
            call SetMatrixValue (ReservoirVolume,        Me%ExtVar%WorkSize2D, 0.0)
            call SetMatrixValue (ReservoirPercFull,      Me%ExtVar%WorkSize2D, 0.0)
            call SetMatrixValue (ReservoirDNInflow,      Me%ExtVar%WorkSize2D, 0.0)
            call SetMatrixValue (ReservoirDischarges,    Me%ExtVar%WorkSize2D, 0.0)
            call SetMatrixValue (ReservoirSurfaceFluxes, Me%ExtVar%WorkSize2D, 0.0)
            call SetMatrixValue (ReservoirOutflow,       Me%ExtVar%WorkSize2D, 0.0)
            call SetMatrixValue (ReservoirConc,          Me%ExtVar%WorkSize2D, 0.0)
            call SetMatrixValue (ReservoirBottomConc,    Me%ExtVar%WorkSize2D, 0.0)
            
            !Create grid results
            CurrReservoir => Me%FirstReservoir
            do while (associated(CurrReservoir))
                
                ReservoirVolume        (CurrReservoir%GridI, CurrReservoir%GridJ) = CurrReservoir%VolumeNew
                ReservoirPercFull      (CurrReservoir%GridI, CurrReservoir%GridJ) = CurrReservoir%PercFull
                ReservoirDNInflow      (CurrReservoir%GridI, CurrReservoir%GridJ) = CurrReservoir%DNInflow
                ReservoirDischarges    (CurrReservoir%GridI, CurrReservoir%GridJ) = CurrReservoir%Discharges
                ReservoirSurfaceFluxes (CurrReservoir%GridI, CurrReservoir%GridJ) = CurrReservoir%SurfaceFluxes
                ReservoirOutflow       (CurrReservoir%GridI, CurrReservoir%GridJ) = CurrReservoir%Outflow
                
                CurrReservoir => CurrReservoir%Next
            enddo
            

            call HDF5WriteData   (Me%ObjHDF5, "//Results/volume",              &
                                  "volume", "m3",                              &
                                  Array2D      = ReservoirVolume,              &
                                  OutputNumber = Me%OutPut%NextOutPut,         &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR010'
       
            call HDF5WriteData   (Me%ObjHDF5, "//Results/percentage full",     &
                                  "percentage full", "%",                      &
                                  Array2D      = ReservoirPercFull,            &
                                  OutputNumber = Me%OutPut%NextOutPut,         &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR020'

            call HDF5WriteData   (Me%ObjHDF5, "//Results/inflow",              &
                                  "inflow", "m3/s",                            &
                                  Array2D      = ReservoirDNInflow,            &
                                  OutputNumber = Me%OutPut%NextOutPut,         &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR030'
            
            call HDF5WriteData   (Me%ObjHDF5, "//Results/discharges",          &
                                  "discharges", "m3/s",                        &
                                  Array2D      = ReservoirDischarges,          &
                                  OutputNumber = Me%OutPut%NextOutPut,         &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR040'
            
            call HDF5WriteData   (Me%ObjHDF5, "//Results/surface fluxes",      &
                                  "surface fluxes", "m3/s",                    &
                                  Array2D      = ReservoirSurfaceFluxes,       &
                                  OutputNumber = Me%OutPut%NextOutPut,         &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR050'
            
            call HDF5WriteData   (Me%ObjHDF5, "//Results/outflow",             &
                                  "outflow", "m3/s",                           &
                                  Array2D      = ReservoirOutflow,             &
                                  OutputNumber = Me%OutPut%NextOutPut,         &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR060'
           
            
            !Properties
            CurrProperty => Me%FirstProperty
            do while (associated(CurrProperty))
                
                if (CurrProperty%ComputeOptions%HDF) then
                    
                    CurrReservoir => Me%FirstReservoir
                    do while (associated(CurrReservoir))
                
                        ReservoirConc(CurrReservoir%GridI, CurrReservoir%GridJ) = &
                            CurrProperty%Concentration(CurrReservoir%Position) 
                                        
                        if(CurrProperty%ComputeOptions%BottomFluxes) then 
                        
                            ReservoirBottomConc(CurrReservoir%GridI,CurrReservoir%GridJ) = &
                                CurrProperty%BottomConc(CurrReservoir%Position)
                        
                        endif
                        CurrReservoir => CurrReservoir%Next
                    enddo
                    
                    call HDF5WriteData   (Me%ObjHDF5, "//Results/"//trim(CurrProperty%ID%Name),         &
                                            trim(CurrProperty%ID%Name), "mg/L",                           &
                                            Array2D      = ReservoirConc,                                 &
                                            OutputNumber = Me%OutPut%NextOutPut,                          &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR070'                    
                                        
                    if(CurrProperty%ComputeOptions%BottomFluxes) then                         
                        
                        call HDF5WriteData   (Me%ObjHDF5, "//Results/"//trim(CurrProperty%ID%Name)//"_Bottom",  &
                                                trim(CurrProperty%ID%Name)//"_Bottom", "kg/m2",                 &
                                                Array2D      = ReservoirBottomConc,                           &
                                                OutputNumber = Me%OutPut%NextOutPut,                          &
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR080'
                    endif                    
                    
                    
                endif
                
                CurrProperty => CurrProperty%Next
            enddo
            

            deallocate(ReservoirVolume         )
            deallocate(ReservoirPercFull       )
            deallocate(ReservoirDNInflow       )
            deallocate(ReservoirDischarges     )
            deallocate(ReservoirSurfaceFluxes  )
            deallocate(ReservoirOutflow        )
            deallocate(ReservoirConc           )
            deallocate(ReservoirBottomConc     )
            
            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR090'

            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

        endif    
    
    end subroutine WriteOutputHDF
    
    
    !---------------------------------------------------------------------------           
    
    subroutine WriteFinalHDF(IsFinalFile)
    
        !Arguments--------------------------------------------------------------
        logical                                     :: IsFinalFile
        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        character(LEN = PathLength)                 :: FileName
        integer                                     :: ObjHDF5           = 0
        integer                                     :: ObjHorizontalGrid = 0
        integer                                     :: ObjGridData       = 0
        integer                                     :: ObjBasinGeometry  = 0        
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        type (T_Time)                               :: Actual             
        real                                        :: Total_Mass_Created
        character (Len = StringLength)              :: str_mass_created, string_to_be_written     
        real, dimension(:,:), pointer               :: ReservoirVolume
        real, dimension(:,:), pointer               :: ReservoirConc, ReservoirBottomConc
        real, dimension(:,:), pointer               :: ReservoirMassCreated
        integer, dimension(:,:), pointer            :: BasinPoints       => null()
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB       
        real, pointer, dimension(:, :)              :: GridData
        type (T_Reservoir), pointer                 :: CurrReservoir        
        !Begin------------------------------------------------------------------
        
            
        
       !Build all Topography stuff
        
        !Constructs Horizontal Grid
        call ConstructHorizontalGrid(ObjHorizontalGrid, Me%ExtVar%TopographicFile, &
                                        STAT = STAT_CALL)           
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR010'

        !Constructs GridData
        call ConstructGridData      (ObjGridData, ObjHorizontalGrid,           &
                                        FileName = Me%ExtVar%TopographicFile,  &
                                        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR020'

        !Constructs BasinGeometry
        call ConstructBasinGeometry (ObjBasinGeometry, ObjGridData,            &
                                        ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHdf5Output - ModuleReservoirs - ERR030'

        !Gets BasinPoints
        call GetBasinPoints         (ObjBasinGeometry, BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR040'        
        
        !Gets the size of the grid
        call GetHorizontalGridSize (ObjHorizontalGrid,                               &
                                    Size     = Me%ExtVar%Size2D,                     &
                                    WorkSize = Me%ExtVar%WorkSize2D,                 &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR050'                           
            
        !Bounds
        WorkILB = Me%ExtVar%WorkSize2D%ILB 
        WorkIUB = Me%ExtVar%WorkSize2D%IUB 

        WorkJLB = Me%ExtVar%WorkSize2D%JLB 
        WorkJUB = Me%ExtVar%WorkSize2D%JUB              
        
        
        
        
        
        !Checks if it's at the end of the run 
        !or !if it's supposed to overwrite the final HDF file
        !if ((Me%ExtVar%Now == Me%ExtVar%EndTime) .or. Me%Output%RestartOverwrite) then
        if (IsFinalFile .or. Me%Output%RestartOverwrite) then

            Filename = trim(Me%Files%FinalFile)

        else

            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%ExtVar%ActualTime))//".fin")

        endif
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5,                                                     &
                            trim(filename),                                              &
                            HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalFile - ModuleReservoirs - ERR10'

        
        
        
        !Time
        Actual   = Me%ExtVar%ActualTime
         
        call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
        !Writes Time
        TimePtr => AuxTime
        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleReservoirs - ERR11'

        call HDF5WriteData  (ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                             Array1D = TimePtr, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleReservoirs - ERR12'
                                                
        
        

        !Write the Horizontal Grid
        call WriteHorizontalGrid(ObjHorizontalGrid, ObjHDF5,       &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR070'

        !Gets a pointer to GridData
        call GetGridData      (ObjGridData, GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR080'

        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR090'

        !Writes the GridData
        call HDF5WriteData   (ObjHDF5, "/Grid", "Bathymetry", "m",            &
                              Array2D = GridData,                                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR100'

        !Writes the WaterPoints
        call HDF5WriteData   (ObjHDF5, "/Grid", "BasinPoints", "-",         &
                              Array2D = BasinPoints,            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR110'

        
        
        !Output volume and conc
        
        allocate(ReservoirVolume         (WorkILB:WorkIUB, WorkJLB:WorkJUB))
        allocate(ReservoirConc           (WorkILB:WorkIUB, WorkJLB:WorkJUB))
        allocate(ReservoirBottomConc     (WorkILB:WorkIUB, WorkJLB:WorkJUB))
        allocate(ReservoirMassCreated    (WorkILB:WorkIUB, WorkJLB:WorkJUB))
        call SetMatrixValue (ReservoirVolume,        Me%ExtVar%WorkSize2D, 0.0)
        call SetMatrixValue (ReservoirConc,          Me%ExtVar%WorkSize2D, 0.0)
        call SetMatrixValue (ReservoirBottomConc,    Me%ExtVar%WorkSize2D, 0.0)
        call SetMatrixValue (ReservoirMassCreated,   Me%ExtVar%WorkSize2D, 0.0)
            
        !Create grid results
        CurrReservoir => Me%FirstReservoir
        do while (associated(CurrReservoir))
                
            ReservoirVolume        (CurrReservoir%GridI, CurrReservoir%GridJ) = CurrReservoir%VolumeNew
                
            CurrReservoir => CurrReservoir%Next
        enddo
            

        call HDF5WriteData   (ObjHDF5, "//Results/volume",              &
                                "volume", "m3",                              &
                                Array2D      = ReservoirVolume,              &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR010'        
            
      
        
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))

            !Sets limits for next write operations
            call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                                  WorkJUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR090'

            CurrReservoir => Me%FirstReservoir
            do while (associated(CurrReservoir))
                
                ReservoirConc(CurrReservoir%GridI, CurrReservoir%GridJ) = PropertyX%Concentration(CurrReservoir%Position)
                                                                           
                if(PropertyX%ComputeOptions%BottomFluxes) then 
                        
                    ReservoirBottomConc(CurrReservoir%GridI,CurrReservoir%GridJ)= PropertyX%BottomConc(CurrReservoir%Position)

                endif
                        
                if (PropertyX%ComputeOptions%MinConcentration .and. Me%ExtVar%ActualTime == Me%ExtVar%EndTime) then

                    ReservoirMassCreated(CurrReservoir%GridI,CurrReservoir%GridJ)= PropertyX%MassCreated(CurrReservoir%Position)
                    
                endif                        
                        
                        
                CurrReservoir => CurrReservoir%Next
            enddo           
            
            
            call HDF5WriteData   (ObjHDF5, "//Results/"//trim(PropertyX%ID%Name),         &
                                    trim(PropertyX%ID%Name), "mg/L",                           &
                                    Array2D      = ReservoirConc,                                 &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR070'                    
                                        
            if(PropertyX%ComputeOptions%BottomFluxes) then                         
                        
                call HDF5WriteData   (ObjHDF5, "//Results/"//trim(PropertyX%ID%Name),         &
                                        trim(PropertyX%ID%Name), "kg/m2",                          &
                                        Array2D      = ReservoirConc,                                 &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutputHDF - ModuleReservoirs - ERR080'
            endif
                        
            if (PropertyX%ComputeOptions%MinConcentration .and. Me%ExtVar%ActualTime == Me%ExtVar%EndTime) then
                            
                call HDF5WriteData   (ObjHDF5,                                        &
                                        "/Results/"//trim(PropertyX%ID%Name)//" Mass Created",& 
                                        "Property Mass Created",                        &
                                        "g",                                            &
                                        Array2D      = ReservoirMassCreated,            &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleReservoirs - ERR10.5'

                !g/1000 = kg, to avoid big numbers
                Total_Mass_Created = SUM(PropertyX%MassCreated)/1000

                write(str_mass_created, '(f20.8)') Total_Mass_Created
      
                string_to_be_written = 'Due to MinConcentration Reservoirs Total mass (kg) created on property ' //&
                                        trim(adjustl(adjustr(PropertyX%ID%Name)))//' = ' //&
                                        trim(adjustl(adjustr(str_mass_created))) 
            
                !Writes total mass created to "Error_and_Messages.log" file
                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)
                    
            endif
            

            PropertyX => PropertyX%Next

        enddo                        
        
        
        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR120'
       
        
        deallocate(ReservoirVolume         )
        deallocate(ReservoirConc           )      
        deallocate(ReservoirBottomConc     )   
        deallocate(ReservoirMassCreated    )   
        
        !Unget all stuff
        !Ungets the GridData
        call UngetGridData (ObjGridData, GridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR130'        
                    
        !Ungets BasinPoints
        call UngetBasin             (ObjBasinGeometry, BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR150'                     
        
        
        call KillBasinGeometry      (ObjBasinGeometry,   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR160'

        call KillGridData           (ObjGridData,        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR170'

        call KillHorizontalGrid     (ObjHorizontalGrid,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR180'            
        
        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalHDF - ModuleReservoirs - ERR0190'         
    
    end subroutine WriteFinalHDF
    
    !---------------------------------------------------------------------------
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillReservoirs(ObjReservoirsID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjReservoirsID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, i, STAT_CALL              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           
        logical                             :: IsFinalFile
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjReservoirsID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mReservoirs_,  Me%InstanceID)
            
            IsFinalFile = .true.
            call WriteFinalHDF(IsFinalFile)            

            if (nUsers == 0) then

                !Kills the TimeSerie
                do i = 1, Me%TimeSerie%nNodes
                    call KillTimeSerie(Me%TimeSerie%ObjTimeSerie(i), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillReservoirs - ModuleReservoirs - ERR01'
                end do    

                
                if (Me%OutPut%HDF) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillReservoirs - ModuleReservoirs  - ERR02'
                endif
                
                !DN May also construct
                if (Me%ComputeOptions%Discharges) then
                    call Kill_Discharges(Me%ObjDischarges, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillReservoirs - ModuleReservoirs - ERR02a'
                endif                  
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillReservoirs - ModuleReservoirs - ERR04'
                                
                !Deallocates Instance
                call DeallocateInstance ()

                ObjReservoirsID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillReservoirs
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Reservoirs), pointer          :: AuxObjReservoirs
        type (T_Reservoirs), pointer          :: PreviousObjReservoirs

        !Updates pointers
        if (Me%InstanceID == FirstObjReservoirs%InstanceID) then
            FirstObjReservoirs => FirstObjReservoirs%Next
        else
            PreviousObjReservoirs => FirstObjReservoirs
            AuxObjReservoirs      => FirstObjReservoirs%Next
            do while (AuxObjReservoirs%InstanceID /= Me%InstanceID)
                PreviousObjReservoirs => AuxObjReservoirs
                AuxObjReservoirs      => AuxObjReservoirs%Next
            enddo

            !Now update linked list
            PreviousObjReservoirs%Next => AuxObjReservoirs%Next

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

    subroutine Ready (ObjReservoirs_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjReservoirs_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjReservoirs_ID > 0) then
            call LocateObjReservoirs (ObjReservoirs_ID)
            ready_ = VerifyReadLock (mRESERVOIRS_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjReservoirs (ObjReservoirs_ID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjReservoirs_ID

        !Local-----------------------------------------------------------------

        Me => FirstObjReservoirs
        do while (associated (Me))
            if (Me%InstanceID == ObjReservoirs_ID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleReservoirs - LocateObjReservoirs - ERR01'

    end subroutine LocateObjReservoirs

    !--------------------------------------------------------------------------

end module ModuleReservoirs