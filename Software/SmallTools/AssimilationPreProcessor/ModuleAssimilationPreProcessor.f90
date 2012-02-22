!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : AssimilationPreProcessor
! MODULE        : AssimilationPreProcessor
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2007
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to create Covariance HDF5 file from HDF5 files
!
!------------------------------------------------------------------------------

!DataFile
!
!   <BeginHDF5File>                                                                 
!   NAME                    : char                  [-]         !Name of HDF5 file with data for covariance 
!   <EndHDF5File>                                               !(specify one block for each HDF5 file to use)
!                                                               !(same block is used for any type of file)
!
!   START_TIME              : YYYY MM DD HH MM SS   [-]         !Start date of time window for covariance
!   END_TIME                : YYYY MM DD HH MM SS   [-]         !End date of time window for covariance
!
!   HDF5_MAP_ITEM           : char                  [-]         !HDF5 mapping item name
!
!   3D_HDF5                 : 0/1                   [0]         !Input HDF5 files are 3D
!                                                               !0 - 2D ; 1 - 3D
!
!   METHOD                  : 1                     [1]         !Assimilation method which will use covariance
!                                                               !1 - SEEK
!
!   STATECOV_RANK           : integer               [-]         !Rank of covariance matrix (no. EOFs)
!
!   NORMALIZATION           : 0/1                   [0]         !Normalization of states prior covariance calculation
!
!   DECIMATION_FACTOR       : integer               [0]         !Number in a row of states discarded for calculation
!                                                               !(e.g. if 3 then only 1 state in 4 is considered)
!
!   OUTPUTFILENAME          : char                  [-]         !Name of produced HDF5 file with analysis results
!
!   CYCLIC_BOUNDARY         : 0/1                   [0]         !Cyclic boundary of domain
!
!   CYCLIC_DIRECTION        : 1/2                   [-]         !Direction of cyclic boundary (1 - XX, 2 - YY)
!
!   STATE_RECONSTRUCTION    : 0/1                   [0]         !States reconstruction with the specified covariance
!                                                               !rank
!
!   STATE_OUTPUTFILENAME    : char                  [-]         !Name of produced HDF5 file with reconstructed states
!
!   <beginproperty>                                            
!   NAME                    : char                  [-]         !Property name (should be equal to the one specified 
!                                                               !in ModuleGlobalData)
!
!   UNITS                   : char                  [-]         !Property units
!
!   DIMENSION               : 2D/3D                 [-]         !Property rank (should be consistent with MOHID
!                                                               !convention)
!
!   HDF_GROUP               : char                  [-]         !Path of the HDF5 group in HDF5 file for property
!
!   STATE_WINDOW            : 4/6*integer   [4/6*FillValueInt]  !ilb, jlb, iub, jub (2D), klb, kub (3D)
!
!   TYPE_ZUV                : U/V/Z                 [Z]         !Type of grid where property is defined in input data: 
!                                                               !Z = cell center, 
!                                                               !U = cell faces U, V = cell faces V
!   <endproperty>                                               !(specify one block for each property)
!

Module ModuleAssimilationPreProcessor

    use ModuleGlobalData
    use ModuleTime,                     only : T_Time, operator(+),                 &
                                               StartComputeTime,                    &
                                               GetComputeCurrentTime,               &
                                               operator(.LT.), operator(.GT.),      &
                                               operator(.EQ.), operator(.GE.),      &
                                               operator(-), operator(.NE.),         &
                                               SetDate, operator(.LE.),             &
                                               ActualizeCurrentTime, ExtractDate
    use ModuleTimeSerie,                only : StartTimeSerieInput,                 &
                                               GetTimeSerieDataMatrix, KillTimeSerie
    use ModuleEnterData,                only : ConstructEnterData, KillEnterData,   &
                                               GetData, ExtractBlockFromBuffer,     &
                                               Block_Unlock, WriteDataLine
    use ModuleHDF5
    use ModuleFunctions,                only : Check_Hydrodynamic_Property,         &
                                               Check_Water_Property,                &
                                               SetMatrixValue, ConstructPropertyID

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartAssimPreProcessor
    private ::      ConstructAssimPreProcessor
    private ::          ReadKeywords
    private ::              ReadGlobalData
    private ::              ReadProperties
    private ::                  ConstructProperty
    private ::                  AddProperty
    private ::              ReadHDF5FileName
    private ::                  AddHDF5File
    private ::                  ConstructHDF5File
    private ::          OpenAndClassHDF5Files
    private ::              HDF5TimeInstant
    private ::              HDF5Evaluator
    private ::                  SamePeriodHDF5
    private ::              AddCovHDF5File
    private ::              KillIndividualHDF5File
    private ::              ObtainInstantsTimes
    private ::              AddCovInstantsTimes
    private ::          OpenOutputFile
    private ::              ConstructHDF5Grid
    private ::              ConstructStateDimension
    private ::              KillGridFields
    private ::              StartExpCoefTimeSeries
    private ::                  AllocateTimeSerieBuffer
    private ::          AllocateVariables

    !Selector
    
    !Modifier
    public  :: ModifyAssimPreProcessor
    private ::      StateStatisticsCalculation
    private ::          OpenAndReadHDF5File
    private ::              ReadPropertyFields
    private ::                  ConvertFieldToFaces
    private ::          StatisticsCalculation2D
    private ::          StatisticsCalculation3D
    private ::          KillIndividualPropertyFields
    private ::          KillIndividualPropertyStats
    private ::              CenterField
    private ::          WriteNormFactorToOutput
    private ::      CovarianceCalculation
    private ::          StatesCovarianceCalculation
    private ::              OpenAndReadHDF5OneInstant
    private ::      EOFAnalysis
    private ::          NormalRand
    private ::              UniformRand01
    private ::          IEigenPower
    private ::              ProductMatVec
    private ::          ExpansionCoefCalculation
    private ::              WriteExpCoefTimeSerie
    private ::                  WriteBufferToFile
    private ::              KillExpCoefTimeSeries
    private ::          TranslateEOF
    private ::          WriteEOFToOutput
    private ::      StateReconstruction
    private ::          WriteStateToOutput

    !Destructor
    public  :: KillAssimPreProcessor
    private ::      KillIndividualFullHDF5File
    private ::      KillIndividualProperty                                                     

    !Management
    
    !Interfaces----------------------------------------------------------------

    !Parameter-----------------------------------------------------------------

    !Data assimilation methods
    integer, parameter                              :: SEEK                 = 1

    !Block
    character(LEN = StringLength), parameter        :: block_prop_begin     = '<beginproperty>'
    character(LEN = StringLength), parameter        :: block_prop_end       = '<endproperty>'
    character(LEN = StringLength), parameter        :: block_hdf5_begin     = '<BeginHDF5File>'
    character(LEN = StringLength), parameter        :: block_hdf5_end       = '<EndHDF5File>'

    !Property dimensions 
    integer, parameter                              :: Dim_2D               = 2
    integer, parameter                              :: Dim_3D               = 3

    !Direction
    integer, parameter                              :: DirectionX_          = 1
    integer, parameter                              :: DirectionY_          = 2

    !Types---------------------------------------------------------------------

    type       T_Field
        character(len=StringLength)                 :: Name
        character(len=StringLength)                 :: Units
        integer                                     :: IDNumber             = 0
        real, dimension(:,:  ),     pointer         :: Values2D
        real, dimension(:,:,: ),    pointer         :: Values3D 
        type(T_Field),              pointer         :: Next                 => null()
    end type  T_Field

    type       T_Grid
        character(len=StringLength)                 :: Name
        character(len=StringLength)                 :: AditionalName
        character(len=StringLength)                 :: Units
        real, dimension(:,:  )      , pointer       :: RealValues2D
        real, dimension(:,:,: )     , pointer       :: RealValues3D
        integer,    dimension(:,:,:), pointer       :: IntegerValues3D
        integer,    dimension(:,:  ), pointer       :: IntegerValues2D
        integer                                     :: Position        
    end type  T_Grid

    type T_SimpleStatistics       
        real, dimension(:, :, :), pointer           :: Average
        real, dimension(:, :, :), pointer           :: SquareAverage
        real, dimension(:, :, :), pointer           :: StandardDeviation

        real, dimension(:, :   ), pointer           :: Average2D
        real, dimension(:, :   ), pointer           :: SquareAverage2D
        real, dimension(:, :   ), pointer           :: StandardDeviation2D

        real                                        :: RunPeriod
        type (T_Time)                               :: LastCalculation
    end type T_SimpleStatistics

    type       T_Property
        type (T_PropertyID)                         :: ID                   
        integer                                     :: Dim                  = null_int
        integer                                     :: TypeZUV              = null_int
        integer                                     :: Rank
        character(len=PathLength)                   :: Group
        
        logical                                     :: ConvertToFaces       = .false.
        
        type (T_Size2D)                             :: Window2D
        type (T_Size3D)                             :: Window
        
        integer                                     :: ModuleType           = null_int
        integer                                     :: FirstStatePosition
        integer                                     :: StateVarNumber

        type (T_SimpleStatistics)                   :: Statistics
        real                                        :: AverStandardDev

        type(T_Field), pointer                      :: FirstField
        type(T_Field), pointer                      :: CurrentField

        character(StringLength)                     :: LastName             = null_str
        integer                                     :: LastRank

        logical                                     :: CyclicBoundary       = .false.

        type (T_Property), pointer                  :: Next, Prev           => null()
    end type T_Property

    type       T_CovarianceTime
        type(T_Time)                                :: Time
        type(T_CovarianceTime), pointer             :: Next                 => null()
    end type  T_CovarianceTime

    type       T_ExpansCoefTS
        character(len=PathLength)                   :: Name                 = null_str
        integer                                     :: ObjTimeSerie         = 0
        integer                                     :: BufferCount          = null_int
        real(8), dimension(:), pointer              :: TimeSerieData
        type (T_Time), dimension(:), pointer        :: TimeBuffer
        real, dimension(:,:), pointer               :: TSDataMatrix
    end type  T_ExpansCoefTS

    type       T_CyclicBoundary
        logical                                     :: ON
        integer                                     :: Direction
    end type T_CyclicBoundary

    type       T_HDF5File
        integer                                     :: HDFID                = 0
        character(len=StringLength)                 :: Name
        type(T_Time)                                :: StartTime
        type(T_Time)                                :: EndTime
        type(T_Time)                                :: StartFieldTime
        type(T_Time)                                :: EndFieldTime
        integer                                     :: Rank
        integer                                     :: NumberOfInstants
        type(T_Time), dimension(:), pointer         :: InstantsArray 
        integer                                     :: StartInstant, EndInstant
        type(T_CovarianceTime), pointer             :: FirstInstantTime
        logical                                     :: WaterProp = .false.
        logical                                     :: HydroProp = .false.
        type(T_HDF5File), pointer                   :: Next, Prev           => null()
    end type  T_HDF5File

    private :: T_AssimPreProcessor
    type       T_AssimPreProcessor
        type (T_Size3D)                             :: Size, WorkSize

        character(PathLength)                       :: DataFile
        character(len=PathLength)                   :: OutputFileName
        integer                                     :: ObjCovHDF5           = 0
        type(T_Time)                                :: StartTime
        type(T_Time)                                :: EndTime
        type(T_Time)                                :: CurrentTime
        type(T_Time)                                :: LastStatCalculation
        logical                                     :: HDF5_3D
        integer                                     :: Method
        logical                                     :: Normalization
        integer                                     :: StateCovRank
        integer                                     :: DecimationFactor
        type(T_CyclicBoundary)                      :: CyclicBoundary
        logical                                     :: StateReconstruction
        character(len=PathLength)                   :: StateRecFileName
        integer                                     :: ObjStateRecHDF5

        integer                                     :: StateVarNumber
        integer                                     :: PropertiesNumber
        integer                                     :: StatesNumber
        type(T_CovarianceTime), pointer             :: FirstCovarianceTime
        integer                                     :: CurrentStateVar
        integer                                     :: CurrentState

        type (T_Property), pointer                  :: FirstProperty
        type (T_Property), pointer                  :: LastProperty
        logical                                     :: Properties2D         = .false.
        logical                                     :: PropertiesInFaces    = .false.

        type (T_HDF5File), pointer                  :: FirstHDF5File
        type (T_HDF5File), pointer                  :: FirstHydroCovHDF5File
        type (T_HDF5File), pointer                  :: LastHydroCovHDF5File
        type (T_HDF5File), pointer                  :: FirstWaterCovHDF5File
        type (T_HDF5File), pointer                  :: LastWaterCovHDF5File
        integer                                     :: NumberHydroCovHDF5   = 0
        integer                                     :: NumberWaterCovHDF5   = 0
        real                                        :: RegularDT
        logical                                     :: VariableDT

        real(8), dimension(:, :),  pointer          :: StateDeviation
        real, dimension(:),  pointer                :: AverageState
        real, dimension(:),  pointer                :: StandardDeviation
        real, dimension(:, :),  pointer             :: Covariance
        real(8), dimension(:, :),  pointer          :: LMatrix
        real, dimension(:),  pointer                :: UVector
        real                                        :: CovTrace
        real(8), dimension(:, :),  pointer          :: StatesLMatrix
        type(T_ExpansCoefTS), dimension(:), pointer :: ExpansionCoef
        real, dimension(:), pointer                 :: StateVector

        integer                                     :: MaxBufferSize        = null_int
        integer                                     :: BufferSize           = null_int

        logical                                     :: FullCovariance = .false. !just test
        logical                                     :: FullEOFAnalysis = .false. !just test

        type(T_Grid)                                :: Bathymetry
        type(T_Grid)                                :: ConnectionX
        type(T_Grid)                                :: ConnectionY
        type(T_Grid)                                :: Mapping
        type(T_Grid)                                :: MappingFacesU
        type(T_Grid)                                :: MappingFacesV
        type(T_Grid)                                :: Latitude
        type(T_Grid)                                :: Longitude
        logical                                     :: AditionalMap

        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjTime              = 0

        type(T_AssimPreProcessor), pointer          :: Next                 => null()
    end type  T_AssimPreProcessor

    !Global Module Variables
    type (T_AssimPreProcessor), pointer             :: Me                   => null()

    integer                                         :: mAssimPreProcessor_  = 0 !just to compile

    !--------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartAssimPreProcessor(ObjAssimPreProcessorID, DataFile)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjAssimPreProcessorID 
        character(PathLength), intent(IN)               :: DataFile    

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        !Assures nullification of the global variable
        allocate(Me)

        !Returns ID
        ObjAssimPreProcessorID          = 1

        !Atribute the name of data file            
        Me%DataFile = DataFile

        call ConstructAssimPreProcessor

        !----------------------------------------------------------------------

    end subroutine StartAssimPreProcessor
 
    !--------------------------------------------------------------------------

    subroutine ConstructAssimPreProcessor

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real                                            :: DT
        type(T_Time)                                    :: RunCovBeginTime
        type (T_HDF5File),       pointer                :: ObjHDF5File
        integer                                         :: STAT_CALL

        !------------------------------------------------------------------------

        !Read input file
        call ReadKeywords

        !Construct HDF5 lists
        call OpenAndClassHDF5Files

        !Construct output
        call OpenOutputFile

        !Construct time reducing DT from start stat time
        DT = - Me%RegularDT 
        !(because there isn't a minus DT subroutine)

        if (associated(Me%FirstHydroCovHDF5File)) then

            RunCovBeginTime = Me%FirstHydroCovHDF5File%InstantsArray(1) +           &
                              (Me%DecimationFactor + 1)*DT

            ObjHDF5File => Me%LastHydroCovHDF5File
        else

            RunCovBeginTime = Me%FirstWaterCovHDF5File%InstantsArray(1) +           &
                              (Me%DecimationFactor + 1)*DT

            ObjHDF5File => Me%LastWaterCovHDF5File
        endif

        call StartComputeTime(Me%ObjTime, RunCovBeginTime,                          &
                              ObjHDF5File%InstantsArray                             &
                              (ObjHDF5File%NumberOfInstants),                       &
                              Me%RegularDT, Me%VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
        stop 'ConstructAssimPreProcessor - ModuleAssimilationPreProcessor - ERR01'

        !Actualized the time
        call GetComputeCurrentTime(Me%ObjTime,                                      &
                                   Me%CurrentTime, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructAssimPreProcessor - ModuleAssimilationPreProcessor - ERR02'

        !Allocate global variables
        call AllocateVariables

        !----------------------------------------------------------------------

    end subroutine ConstructAssimPreProcessor
 
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL

        !------------------------------------------------------------------------

        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadKeywords - ModuleAssimilationPreProcessor - ERR01'

        call ReadGlobalData

        call ReadProperties

        call ReadHDF5FileName

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadKeywords - ModuleAssimilationPreProcessor - ERR02'

        !----------------------------------------------------------------------

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine ReadGlobalData

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag

        !------------------------------------------------------------------------

        ! Obtain the start and end times for the covariance' calculation
        ! Start Time
        call GetData(Me%StartTime, Me%ObjEnterData, iflag,                      &
                     keyword      = 'START_TIME',                               &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR01'   

        ! End Time 
        call GetData(Me%EndTime,   Me%ObjEnterData, iflag,                      &
                     keyword      = 'END_TIME',                                 &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR02'   

        ! Verifies consistency in time variables
        if (Me%EndTime .lt. Me%StartTime) then
            write (*,*) 'Covariance End Time is BEFORE Covariance Start Time'
            write (*,*) 'Module :','AssimilationPreProcessor'
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR03'
        endif

        if (Me%EndTime .eq. Me%StartTime) then
            write (*,*) 'Covariance End Time is EQUAL Covariance Start Time'
            write (*,*) 'Module :','AssimilationPreProcessor'
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR04'
        endif

        ! Obtain HDF5 file's Map variable name
        call GetData(Me%Mapping%Name, Me%ObjEnterData, iflag,                   &
                     keyword      = 'HDF5_MAP_ITEM',                            &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR05'            

        ! Obtain HDF5 file's dimension (3D or 2D)
        call GetData(Me%HDF5_3D, Me%ObjEnterData, iflag,                        &
                     keyword      = '3D_HDF5',                                  &
                     default      = .false.,                                    &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR06'            

        ! Obtain cyclic boundary condition switch
        call GetData(Me%CyclicBoundary%ON, Me%ObjEnterData, iflag,              &
                     keyword      = 'CYCLIC_BOUNDARY',                          &
                     default      = .false.,                                    &
                     SearchType   = FromFile,                                   &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR07a'            

        if (Me%CyclicBoundary%ON) then
            call GetData(Me%CyclicBoundary%Direction, Me%ObjEnterData, iflag,   &
                         keyword    = 'CYCLIC_DIRECTION',                       & 
                         SearchType = FromFile,                                 &
                         ClientModule ='ModuleAssimilationPreProcessor',        &
                         STAT       = STAT_CALL)            

            if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                          &
                stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR07b'            
        endif

        ! Get data assimilation method the covariance is intended for 
        call GetData(Me%Method,                                                 &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'METHOD',                                     &
                     Default    = SEEK,                                         &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR08'

        ! Get normalization button state
        call GetData(Me%Normalization,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'NORMALIZATION',                              &
                     Default    = .false.,                                      &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR09'

        ! Get decimation factor
        call GetData(Me%DecimationFactor,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'DECIMATION_FACTOR',                          &
                     Default    = 0,                                            &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR10'

        if (Me%Method == SEEK) then

            ! Get covariance reduced rank
            call GetData(Me%StateCovRank,                                       &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromFile,                                 &
                         keyword    = 'STATECOV_RANK',                          &
                         ClientModule = 'ModuleAssimilationPreProcessor',       &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                          &
                stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR11'
        endif

        ! Get state reconstruction button state
        call GetData(Me%StateReconstruction,                                    &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'STATE_RECONSTRUCTION',                       &
                     Default    = .false.,                                      &
                     ClientModule = 'ModuleAssimilationPreProcessor',           &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR12'

        !The maximum BufferSize is set here to 0.1Mb (for each expansion coef.)
        !This lets perform 25000 outputs to the buffer (considering each output of 4 bytes)
        call GetData(Me%MaxBufferSize,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      ='MAX_BUFFER_SIZE',                           &
                     Default      = 100000,                                     &
                     ClientModule ='ModuleAssimilationPreProcessor',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ReadGlobalData - ModuleAssimilationPreProcessor - ERR13'   

        !----------------------------------------------------------------------

    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------

    subroutine ReadProperties

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: ClientNumber
        logical                                         :: BlockFound
        type(T_Property),   pointer                     :: NewProperty
        logical                                         :: AtLeastOneBlock = .false.

        !------------------------------------------------------------------------

        nullify (Me%FirstProperty)
        nullify (Me%LastProperty)

        ! Obtain state properties for covariance calculation
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                        block_prop_begin, block_prop_end,       &
                                        BlockFound, STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  

                    AtLeastOneBlock = .true.

                    ! Construct a New Property 
                    call ConstructProperty(NewProperty)

                    ! Add new Property to the Assimilation List 
                    call AddProperty(NewProperty)
                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)                                 &
                        stop 'ReadProperties - ModuleAssimilationPreProcessor - ERR01'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)                                     &
                    stop 'ReadProperties - ModuleAssimilationPreProcessor - ERR02'
            end if cd1
        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No property block is indicated in input file: '
            write(*,*) trim(block_prop_begin)
            write(*,*) trim(block_prop_end)
            stop 'ReadProperties - ModuleAssimilationPreProcessor - ERR03'
        end if       

        !----------------------------------------------------------------------

    end subroutine ReadProperties

    !--------------------------------------------------------------------------

    subroutine ConstructProperty(NewProperty)

        !Arguments---------------------------------------------------------------
        type(T_Property),   pointer                     :: NewProperty

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: iflag
        integer                                         :: STAT_CALL
        character(len=StringLength)                     :: Char_TypeZUV
        integer, dimension(:), pointer                  :: aux
        !integer                                         :: WorkILB, WorkIUB, WorkJLB
        !integer                                         :: WorkJUB, WorkKLB, WorkKUB

        !------------------------------------------------------------------------

        !Allocates new property
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR01'
        nullify  (NewProperty%Next)

        !Get keywords from property block
        !ID properties
        call ConstructPropertyID (NewProperty%ID, Me%ObjEnterData, FromBlock)

        !The other
        call GetData(NewProperty%Dim, Me%ObjEnterData, iflag,                   &
                     keyword        = 'DIMENSION',                              &  
                     default        = Dim_3D,                                   &
                     SearchType     = FromBlock,                                &
                     ClientModule   = 'ModuleSequentialAssimilation',           &
                     STAT           = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR02'

        if (NewProperty%Dim /= Dim_3D) then
            if (NewProperty%Dim /= Dim_2D) then
                write(*,*)  
                write(*,*) 'Property must be either 3D or 2D: ',                &
                            trim(NewProperty%ID%Name)
                stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR03'

            elseif (NewProperty%ID%IDNumber /= WaterLevel_) then
                write(*,*)  
                write(*,*) 'Only property Water Level is registered as 2D in MOHID Water'
                stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR04'
            endif
        endif       

        ! Obtain parameter group
        call GetData(NewProperty%Group,                                         &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromBlock,                                  &
                     keyword      = 'HDF_GROUP',                                &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR05'

        call GetData(Char_TypeZUV, Me%ObjEnterData, iflag,                      &
                     keyword        = 'TYPE_ZUV',                               &  
                     SearchType     = FromBlock,                                &
                     ClientModule   = 'ModuleSequentialAssimilation',           &
                     default        = "Z",                                      &
                     STAT           = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR06'

        NewProperty%TypeZUV  = TranslateTypeZUV(Char_TypeZUV)

        if (NewProperty%TypeZUV == TypeZ_ .and.                                 &
           (NewProperty%ID%IDNumber == VelocityU_ .or.                          &
            NewProperty%ID%IDNumber == VelocityV_)) then
            !check if to convert to faces
            call GetData(NewProperty%ConvertToFaces,                            &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromBlock,                                &
                         keyword    = 'CONVERT_TO_FACES',                       &
                         Default    = .false.,                                  &
                         ClientModule = 'ModuleAssimilationPreProcessor',       &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR06b'
        endif

        if (((NewProperty%TypeZUV == TypeU_) .or. (NewProperty%TypeZUV == TypeV_) &
            .or. (NewProperty%TypeZUV == TypeZ_ .and. NewProperty%ConvertToFaces)) &
            .and. .not. Me%PropertiesInFaces) Me%PropertiesInFaces = .true.

        if (NewProperty%Dim == Dim_2D) then

            allocate (aux(4))

            call GetData(aux,                                                   &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromBlock,                                &
                         keyword    = 'STATE_WINDOW',                           &
                         Default    = FillValueInt,                             &
                         ClientModule ='ModuleSequentialAssimilation',          &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR07'

            if (iflag == 4) then

                NewProperty%Window2D%ILB = aux(1) 
                NewProperty%Window2D%IUB = aux(2)
                NewProperty%Window2D%JLB = aux(3) 
                NewProperty%Window2D%JUB = aux(4)
            else

                write(*,*)  
                write(*,*) 'Spatial window not specified for property: ',       &
                            trim(NewProperty%ID%Name)
                stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR08'
            endif

        else !Dim_3D
            
            allocate (aux(6))

            call GetData(aux,                                                   &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromBlock,                                &
                         keyword    = 'STATE_WINDOW',                           &
                         Default    = FillValueInt,                             &
                         ClientModule ='ModuleSequentialAssimilation',          &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR09'

            if (iflag == 6) then

                NewProperty%Window%ILB = aux(1) 
                NewProperty%Window%IUB = aux(2)
                NewProperty%Window%JLB = aux(3) 
                NewProperty%Window%JUB = aux(4)
                NewProperty%Window%KLB = aux(5) 
                NewProperty%Window%KUB = aux(6) 
            else
                write(*,*)  
                write(*,*) 'Spatial window not specified for property: ',       &
                            trim(NewProperty%ID%Name)
                stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR10'
            endif
        endif

        deallocate (aux)

        !Find if property is from Hydrodynamic or WaterProperties modules
        !(only hydrodynamic and water properties are assumed possible state components)
        if (Check_Hydrodynamic_Property(NewProperty%ID%IDNumber)) then

            NewProperty%ModuleType = 1
            
        else if (Check_Water_Property(NewProperty%ID%IDNumber)) then

            NewProperty%ModuleType = 2

        else
            write (*,*)  
            write (*,*)'State property isnt from Hydrodynamic or WaterProperties modules:'
            write (*,*) trim(NewProperty%ID%Name)
            stop 'ConstructProperty - ModuleAssimilationPreProcessor - ERR11'
        end if
        !(this classification of properties is required so that two lists of hdf5 are
        !constructed)

    end subroutine ConstructProperty

    !--------------------------------------------------------------------------

    subroutine AddProperty(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),   pointer         :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the list a new property
        if (.not.associated(Me%FirstProperty)) then
            Me%PropertiesNumber = 1
            Me%FirstProperty                => NewProperty
            Me%LastProperty                 => NewProperty
        else
            NewProperty%Prev                => Me%LastProperty
            Me%LastProperty%Next            => NewProperty
            Me%LastProperty                 => NewProperty
            Me%PropertiesNumber             = Me%PropertiesNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine AddProperty 

    !--------------------------------------------------------------------------

    subroutine ReadHDF5FileName

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag
        type (T_HDF5File),       pointer                :: NewHDF5File
        integer                                         :: ClientNumber
        logical                                         :: BlockFound
        logical                                         :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

        nullify(Me%FirstHDF5File)

        !Read input hdf5s (these go all to the same list for now)
do1 :   do
          
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = block_hdf5_begin,             &
                                        block_end       = block_hdf5_end,               &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)

cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                
                    AtLeastOneBlock = .true.

                    call AddHDF5File                     (NewHDF5File, Me%FirstHDF5File)

                    call ConstructHDF5File               (NewHDF5File)

                    nullify(NewHDF5File)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                                  & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadHDF5FileName - ModuleAssimilationPreProcessor - ERR02'

                    exit do1

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadHDF5FileName - ModuleAssimilationPreProcessor - ERR03'
            else cd1
                stop 'ReadHDF5FileName - ModuleAssimilationPreProcessor - ERR04'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No HDF5 file block is indicated in input file: '
            write(*,*) trim(block_hdf5_begin)
            write(*,*) trim(block_hdf5_end)
            stop 'ReadHDF5FileName - ModuleAssimilationPreProcessor - ERR05'
        end if       

        !Read output HDF5 file name
        call GetData(Me%OutputFileName, Me%ObjEnterData, iflag,                         &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModuleSequentialAssimilation',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                      &
        stop 'ReadHDF5FileName - ModuleAssimilationPreProcessor - ERR06'   

        if (Me%StateReconstruction) then
            !Read output states reconstruction HDF5 file name
            call GetData(Me%StateRecFileName, Me%ObjEnterData, iflag,                   &
                         keyword      = 'STATE_OUTPUTFILENAME',                         &
                         SearchType   = FromFile,                                       &
                         ClientModule = 'ModuleSequentialAssimilation',                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                                  &
            stop 'ReadHDF5FileName - ModuleAssimilationPreProcessor - ERR07'   
        endif

    end subroutine ReadHDF5FileName

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5File(NewHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),      pointer           :: NewHDF5File

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------


        ! Obtain HDF5 file name
        call GetData(NewHDF5File%Name,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'NAME',                                             &
                     ClientModule = 'ModuleSequentialAssimilation',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5File - ModuleAssimilationPreProcessor - ERR01'

    end subroutine ConstructHDF5File

    !--------------------------------------------------------------------------

    subroutine AddHDF5File(ObjHDF5File, FirstHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),     pointer           :: ObjHDF5File
        type (T_HDF5File),     pointer           :: FirstHDF5File

        !Local-----------------------------------------------------------------
        type (T_HDF5File),     pointer           :: PreviousHDF5File => null()
        type (T_HDF5File),     pointer           :: NewHDF5File

        !Begin-----------------------------------------------------------------

        !Allocates new HDF5 file
        allocate (NewHDF5File)
        nullify  (NewHDF5File%Next)

        !Insert new HDF5 file into list and makes current point to it
        if (.not. associated(FirstHDF5File)) then
            FirstHDF5File            => NewHDF5File
            ObjHDF5File              => NewHDF5File
        else
            PreviousHDF5File         => FirstHDF5File
            ObjHDF5File              => FirstHDF5File%Next
            do while (associated(ObjHDF5File))
                PreviousHDF5File     => ObjHDF5File
                ObjHDF5File          => ObjHDF5File%Next
            enddo
            ObjHDF5File              => NewHDF5File
            PreviousHDF5File%Next    => NewHDF5File
        end if

    end subroutine AddHDF5File

    !--------------------------------------------------------------------------

    subroutine OpenAndClassHDF5Files

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_HDF5File),     pointer              :: HDF5FileX        => null()
        type (T_HDF5File),     pointer              :: HydroHDF5FileX   => null()
        type (T_HDF5File),     pointer              :: WaterHDF5FileX   => null()
        type (T_HDF5File),     pointer              :: HDF5ToKill       => null()
        logical                                     :: exist
        logical                                     :: FirstTimeHydro, FirstTimeWater
        integer                                     :: HDF5_READ
        integer                                     :: STAT_CALL
        type (T_Property), pointer                  :: ObjProperty      => null()
        logical                                     :: Relevant, GroupExist
        character(len=StringLength)                 :: PropertyName
        logical                                     :: AtLeastOneProp = .false.
        type(T_Time), dimension(:), pointer         :: AuxInstantsArray 
        integer                                     :: CurrentInstant, AuxNumberInstants
        real                                        :: DT
        real                                        :: LastDT, AuxDT
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        integer                                     :: Count
        real, parameter                             :: AuxTypeReal = 8

        !Begin-----------------------------------------------------------------

        !Go through input HDF5 files and check their relevance:
        ! - time window
        ! - state properties
        !Add relevant ones to appropriate list in appropriate order

        FirstTimeHydro = .true.
        FirstTimeWater = .true.

        HDF5FileX => Me%FirstHDF5File

        !In a DO cycle open all HDF5 files provided by the user
        do while (associated(HDF5FileX))

            !Verifies if file exists
            inquire(FILE = HDF5FileX%Name, EXIST = exist)
            if (.not. exist) then
                write(*,*)'HDF5 file does not exist:'//trim(HDF5FileX%Name)
                stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR01'
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name),                  &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR02'

            !Obtain start and end times of HDF5 file
            !(obtain number of instants) 
            call GetHDF5GroupNumberOfItems(HDF5FileX%HDFID, "/Time",                    &
                                           HDF5FileX%NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                & 
            stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR03'

            !(obtain HDF5 start time)
            HDF5FileX%StartTime = HDF5TimeInstant(1, HDF5FileX)

            !(obtain HDF5 end time)
            HDF5FileX%EndTime = HDF5TimeInstant(HDF5FileX%NumberOfInstants, HDF5FileX)

            !Get info about the rank and variables present
            !(only data needed for checking are obtained from each file)
            AtLeastOneProp = .false.

            ObjProperty => Me%FirstProperty

            do while(associated(ObjProperty))

                call GetHDF5GroupExist (HDF5FileX%HDFID, ObjProperty%Group, GroupExist)

                !check if file contains parameter required
                if (GroupExist) then  

                    AtLeastOneProp = .true.

                    if (ObjProperty%ModuleType == 1) then
                        
                        HDF5FileX%HydroProp = .true.
                    else

                        HDF5FileX%WaterProp = .true.
                    endif

                    !get field ID, Rank
                    call GetHDF5GroupID(HDF5FileX%HDFID, ObjProperty%Group,             &
                                        1, PropertyName, ObjProperty%ID%Units,          &
                                        ObjProperty%Rank,                               &
                                        STAT = STAT_CALL)                                
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR04'

                    if (Me%HDF5_3D .and. (.not. Me%Properties2D)) then
                        if (ObjProperty%Rank == 2) Me%Properties2D = .true.
                    endif

                    if ((HDF5FileX%HydroProp .and. (.not. FirstTimeHydro)) .or.         &
                        (HDF5FileX%WaterProp .and. (.not. FirstTimeWater))) then
                        !check if name of property is consistent with the previous (same group)
                        if (PropertyName .NE. ObjProperty%LastName) then
                            write(*,*)'HDF5 file does not contain property required:'   &
                                       //trim(HDF5FileX%Name)
                            write(*,*)'Property required:'//trim(ObjProperty%ID%Name)
                        stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR05'
                        end if

                        !check if rank of property is consistent with the previous
                        if (ObjProperty%Rank .NE. ObjProperty%LastRank) then
                            write(*,*)'File rank not consistent to previous rank:'      &
                                       //trim(HDF5FileX%Name)
                            write(*,*)'Property:'//trim(ObjProperty%ID%Name)
                        stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR06'
                        end if
                    end if 

                    ObjProperty%LastRank = ObjProperty%Rank
                    ObjProperty%LastName = PropertyName   
                end if         

                ObjProperty => ObjProperty%Next
            end do

            if (.not. AtLeastOneProp) then

                write(*,*)'HDF5 file does not contain any property required:'           &
                           //trim(HDF5FileX%Name)
                write(*,*)'Remove file block from input file'
                stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR07'
            else

                !Check if HDF5 file is relevant for covariance calculation in terms of time
                call HDF5Evaluator(HDF5FileX, Relevant)

                !If HDF5 file is relevant then obtain key properties and instants
                if (Relevant) then

                    !Get useful time information from file:
                    !Set instant array for this file 
                    allocate(AuxInstantsArray(1:HDF5FileX%NumberOfInstants))

                    !Fill array with instants
                    do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                        AuxInstantsArray(CurrentInstant) =                              &
                            HDF5TimeInstant(CurrentInstant, HDF5FileX)
                    end do

                    !Get start and end instants for this file
                    !select time window begin
                    do CurrentInstant = 1, HDF5FileX%NumberOfInstants

                        if (AuxInstantsArray(CurrentInstant)                            &
                            .ge. HDF5FileX%StartFieldTime) then

                            HDF5FileX%StartInstant = CurrentInstant
                            HDF5FileX%StartFieldTime = HDF5TimeInstant(CurrentInstant,  & 
                                                                       HDF5FileX)              
                            exit
                        
                        end if
                    end do
    
                    !select time window end
                    do CurrentInstant = HDF5FileX%StartInstant,                         &
                            HDF5FileX%NumberOfInstants, (Me%DecimationFactor + 1)

                        if (AuxInstantsArray(CurrentInstant)                            &
                            .eq. HDF5FileX%EndFieldTime) then

                            HDF5FileX%EndInstant = (CurrentInstant)
                            HDF5FileX%EndFieldTime = HDF5TimeInstant(CurrentInstant,    &
                                                                     HDF5FileX)
                            exit

                        elseif (AuxInstantsArray(CurrentInstant)                        &
                                .gt. HDF5FileX%EndFieldTime) then
                            !%EndInstant takes into account the decimation factor

                            HDF5FileX%EndInstant = (CurrentInstant -                    &
                                                    (Me%DecimationFactor + 1))
                            HDF5FileX%EndFieldTime = HDF5TimeInstant(CurrentInstant -   &
                                                    (Me%DecimationFactor + 1), HDF5FileX)
                            exit 

                        end if
                    end do

                    if (HDF5FileX%EndInstant == 0) then
                        HDF5FileX%EndInstant = (CurrentInstant -                        &
                                               (Me%DecimationFactor + 1))
                        HDF5FileX%EndFieldTime = HDF5TimeInstant(CurrentInstant -       &
                                                 (Me%DecimationFactor + 1), HDF5FileX)
                    endif
        
                    AuxNumberInstants = HDF5FileX%NumberOfInstants

                    HDF5FileX%NumberOfInstants = 1 + (HDF5FileX%EndInstant -            &
                                                 HDF5FileX%StartInstant)/               &
                                                 (Me%DecimationFactor + 1)

                    allocate(HDF5FileX%InstantsArray(1:HDF5FileX%NumberOfInstants))

                    Count = 0
                    do CurrentInstant = HDF5FileX%StartInstant, HDF5FileX%EndInstant,   &
                                       (Me%DecimationFactor + 1)
                        Count = Count + 1
                        HDF5FileX%InstantsArray(Count) = AuxInstantsArray(CurrentInstant)
                    enddo

                    !Calculate DT
                    if (HDF5FileX%NumberOfInstants .ge. 2) then
                        DT = HDF5FileX%InstantsArray(2) - HDF5FileX%InstantsArray(1)
                    else 
                        write(*,*) 'HDF5 file with less than 2 time instants:'                
                        write(*,*) HDF5FileX%Name
                    end if

                    !Calculate DT for consistency checking
                    if (AuxNumberInstants .ge. 3) then
                        AuxDT = AuxInstantsArray(3) - AuxInstantsArray(2)
                        Me%RegularDT = AuxDT
                    else
                        !Check if this DT is equal to the DT of last file: they must equal!
                        if ((AuxDT .NE. LastDT) .and. ((.not. FirstTimeHydro) .or.      &
                            (.not. FirstTimeWater)) .and. (AuxNumberInstants .ge. 3)) then
                            write(*,*) 'HDF5 file has not the same time interval as first'                
                            write(*,*) 'HDF5 file: '//trim(HDF5FileX%Name)
                            write(*,*) 'Required time interval: ',                      &
                                        LastDT/(Me%DecimationFactor + 1)
                            write(*,*) 'All HDF5 files must have same time interval.' 
                        stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR08'  
                        end if 
                        LastDT = AuxDT
                        Me%RegularDT = DT
                    end if

                    !Add file to list of relevant files
                    if (HDF5FileX%HydroProp) then 
                        call AddCovHDF5File(HDF5FileX, Me%FirstHydroCovHDF5File,        &
                                            Me%LastHydroCovHDF5File, Me%NumberHydroCovHDF5)
                        !next run of cycle (next property) is not the first one                
                        FirstTimeHydro = .false.
                    endif

                    if (HDF5FileX%WaterProp .and. .not. HDF5FileX%HydroProp) then 
                        call AddCovHDF5File(HDF5FileX, Me%FirstWaterCovHDF5File,        &
                                            Me%LastWaterCovHDF5File, Me%NumberWaterCovHDF5)
                        !next run of cycle (next property) is not the first one                
                        FirstTimeWater = .false.
                    endif
         
                    deallocate(AuxInstantsArray)
                    nullify(AuxInstantsArray)

                end if
            endif

            !Close HDF5
            call killhdf5(HDF5FileX%HDFID)

            HDF5ToKill  => HDF5FileX

            HDF5FileX   => HDF5FileX%Next

            Me%FirstHDF5File => Me%FirstHDF5File%Next

            !Deallocate HDF5 (if relevant a copy is in list of relevant files)
            call KillIndividualHDF5File(HDF5ToKill)

        end do

        if (FirstTimeHydro .and. FirstTimeWater) then

            !No HDF5 file is provided with suitable data
            call ExtractDate(Me%StartTime, Year, Month, Day, Hour,                      &  
                             Minute, Second)        
            write(*,*) 'Data lacking from'      
100         format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'to'      
            call ExtractDate(Me%EndTime, Year, Month, Day, Hour, Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'There are not requested data in the HDF5 files provided.'      
            stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR09'                

        else

            !Check if number of Hydro and Water files is equal if both types exist
            !(this is first test to assure same time instants are present in both types)
            !(saves more expensive checking operations)
            if ((Me%NumberHydroCovHDF5 /= 0) .and. (Me%NumberWaterCovHDF5 /= 0)) then
                if (Me%NumberHydroCovHDF5 /= Me%NumberWaterCovHDF5) then
                    write(*,*) 'Number of hydrodynamic properties HDF5 files different'      
                    write(*,*) 'from number of water properties HDF5 files.'      
                    stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR10'                
                endif
            endif

            !assume VariableDT
            Me%VariableDT = .True.

            !For each HDF5 file needed for the covariance calculation
            HydroHDF5FileX => Me%FirstHydroCovHDF5File
            WaterHDF5FileX => Me%FirstWaterCovHDF5File
            nullify(HDF5FileX)

            Me%StatesNumber = 0

            do while(associated(HydroHDF5FileX) .or. associated(WaterHDF5FileX))
                !the .OR. is for the case only one file type is present

                if (associated(HydroHDF5FileX) .and. associated(WaterHDF5FileX)) then
                    !Check if relevant instants are the same in both files
                    !(just check start and end times)
                    if ((HydroHDF5FileX%StartFieldTime /=                               &
                         WaterHDF5FileX%StartFieldTime) .or.                            &
                        (HydroHDF5FileX%EndFieldTime /=                                 &
                         WaterHDF5FileX%EndFieldTime)) then

                        write(*,*)'Hydrodynamic and Water properties files not having'
                        write(*,*)'same relevant time instants:'                        &
                                   //trim(HydroHDF5FileX%Name)//' '                     &
                                   //trim(WaterHDF5FileX%Name)                          
                        write(*,*)'Property:'//trim(ObjProperty%ID%Name)
                        stop 'OpenAndClassHDF5Files - ModuleAssimilationPreProcessor - ERR11'
                    endif

                    HDF5FileX => HydroHDF5FileX !only one type is required

                    HydroHDF5FileX => HydroHDF5FileX%Next            
                    WaterHDF5FileX => WaterHDF5FileX%Next

                else if (associated(HydroHDF5FileX)) then

                    HDF5FileX => HydroHDF5FileX
                    HydroHDF5FileX => HydroHDF5FileX%Next            
                else

                    HDF5FileX => WaterHDF5FileX
                    WaterHDF5FileX => WaterHDF5FileX%Next
                endif

                !Get instants' times and put them in list
                call ObtainInstantsTimes(HDF5FileX)               
                 
                call AddCovInstantsTimes(HDF5FileX)
                !(instant times are the same if both Hydro and Water files exist)

            end do
        endif

        if (Me%FirstCovarianceTime%Time > Me%StartTime) then
            write(*,*)
            write(*,*) 'Data lacking from'      
            call ExtractDate(Me%StartTime, Year, Month, Day, Hour, Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
            write(*,*) 'to'      
            call ExtractDate(Me%FirstCovarianceTime%Time, Year, Month, Day, Hour,       &  
                             Minute, Second)
            write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
        endif   

        !Calculates the size of the buffer for time series of expansion coef.
        if (Me%StatesNumber * AuxTypeReal > Me%MaxBufferSize) then
            Me%BufferSize  = int(Me%MaxBufferSize / (AuxTypeReal))
        else
            Me%BufferSize  = Me%StatesNumber
        endif

    end subroutine OpenAndClassHDF5Files

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant, ObjHDF5File)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        type(T_HDF5File), pointer               :: ObjHDF5File
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjHDF5File%HDFID, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = ObjHDF5File%HDFID,                        &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
        stop 'HDF5TimeInstant - ModuleAssimilationPreProcessor - ERR01'

        call SetDate(HDF5TimeInstant, Year  = TimeVector(1),                            &
                     Month  = TimeVector(2), Day      = TimeVector(3),                  &
                     Hour   = TimeVector(4), Minute   = TimeVector(5),                  &
                     Second = TimeVector(6))

        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------

    subroutine HDF5Evaluator(HDF5FileX, Relevant)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        logical                                     :: Relevant

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_HDF5File), pointer                   :: OtherFile      

        !Begin-----------------------------------------------------------------

        Relevant = .FALSE.

        !Check if there is a file with same period
        call SamePeriodHDF5(HDF5FileX, OtherFile, STAT_CALL)

        if (STAT_CALL .EQ. SUCCESS_) then 
            write(*,*)'Same period is contained in two HDF5 files:'
            write(*,*) trim(HDF5FileX%Name)
            write(*,*) trim(OtherFile%Name)
            stop 'HDF5Evaluator - ModuleAssimilationPreProcessor - ERR01'
        end if

        !See if the file is to be considered for statistics 
        !Check start and end time
        if ((HDF5FileX%EndTime >= Me%StartTime) .and.                                   &
            (HDF5FileX%EndTime <= Me%EndTime)) then

            !End statistics time is between start and end times of file
            if (HDF5FileX%StartTime < Me%StartTime) then

                !Start statistics time is after start time of file
                HDF5FileX%StartFieldTime = Me%StartTime
                HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                Relevant = .TRUE.
            else 

                !Start statistics time is before start time of file
                HDF5FileX%StartFieldTime = HDF5FileX%StartTime
                HDF5FileX%EndFieldTime = HDF5FileX%EndTime

                Relevant = .TRUE.
            end if
        else if ((HDF5FileX%StartTime >= Me%StartTime) .and.                            &
                 (HDF5FileX%StartTime <= Me%EndTime)) then

            !End statistics time is before end time of file
            HDF5FileX%StartFieldTime = HDF5FileX%StartTime
            HDF5FileX%EndFieldTime = Me%EndTime

            Relevant = .TRUE.

        else if ((HDF5FileX%StartTime < Me%StartTime) .and.                             &
                 (HDF5FileX%EndTime > Me%EndTime)) then

            !Statistics period is contained in file
            HDF5FileX%StartFieldTime = Me%StartTime
            HDF5FileX%EndFieldTime = Me%EndTime

            Relevant = .TRUE.
        end if

    end subroutine HDF5Evaluator

    !--------------------------------------------------------------------------

    subroutine ObtainInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                 :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        type (T_CovarianceTime), pointer          :: NewTime
        type (T_CovarianceTime), pointer          :: ObjTime        => null()
        type (T_CovarianceTime), pointer          :: PreviousTime   => null()
        integer                                   :: CurrentInstant
        integer                                   :: HDF5_READ
        integer                                   :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (ObjHDF5File%HDFID,                                  &
                            trim(ObjHDF5File%Name),                             &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ObtainInstantsTimes - ModuleAssimilationPreProcessor - ERR01'
        
        do CurrentInstant = ObjHDF5File%StartInstant, ObjHDF5File%EndInstant,   &
                            (Me%DecimationFactor + 1)

            !Allocates new instance
            allocate (NewTime)
            nullify  (NewTime%Next)

            NewTime%Time = HDF5TimeInstant(CurrentInstant, ObjHDF5File)

            !Insert New Instance into list and makes Current point to it
            if (.not. associated(ObjHDF5File%FirstInstantTime)) then
            !FirstField should be the same for all HDF5 files 
                ObjHDF5File%FirstInstantTime    => NewTime
                ObjTime                         => NewTime
            else
                PreviousTime                    => ObjHDF5File%FirstInstantTime
                ObjTime                         => ObjHDF5File%FirstInstantTime%Next
                do while (associated(ObjTime))
                    PreviousTime                => ObjTime
                    ObjTime                     => ObjTime%Next
                enddo
                ObjTime                         => NewTime
                PreviousTime%Next               => NewTime
            endif

        end do 

        call killhdf5(ObjHDF5File%HDFID) 

    end subroutine ObtainInstantsTimes

    !--------------------------------------------------------------------------

    subroutine AddCovHDF5File(HDF5FileX, FirstCovHDF5File, LastCovHDF5File, Number)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX, FirstCovHDF5File
        type(T_HDF5File), pointer                   :: LastCovHDF5File
        integer                                     :: Number

        !Local-----------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileAux      => null()
        type(T_HDF5File), pointer                   :: PreviousHDF5File => null()
        type(T_HDF5File), pointer                   :: LastHDF5File     => null()

        !Begin-----------------------------------------------------------------

        if (.not. associated(FirstCovHDF5File)) then

            call CreateCovHDF5File(FirstCovHDF5File, HDF5FileX)
            LastCovHDF5File => FirstCovHDF5File

            nullify(FirstCovHDF5File%Next)
            nullify(FirstCovHDF5File%Prev)

            !current file was added to list
            Number = 1
        else

            if (HDF5FileX%StartFieldTime < FirstCovHDF5File%StartFieldTime) then
                !current file should be the first file in list

                !save the list starting
                LastHDF5File => FirstCovHDF5File

                !make the first element in the list of relevant files equal to current file
                call CreateCovHDF5File(HDF5FileAux, HDF5FileX)
                FirstCovHDF5File => HDF5FileAux

                FirstCovHDF5File%Next => LastHDF5File
                LastHDF5File%Prev => FirstCovHDF5File
                
                nullify(FirstCovHDF5File%Prev)

                !current file was added to list
                Number = Number + 1

                !Adjust end field time if equal to next file start field time
                if ((LastHDF5File%StartFieldTime - FirstCovHDF5File%EndFieldTime)   &
                    .lt. Me%RegularDT) then

                    call AdjustHDF5EndInstant(FirstCovHDF5File)

                endif

            else
                !check next files in list

                !locate previous file in the first file
                PreviousHDF5File => FirstCovHDF5File                  

                do while(associated(PreviousHDF5File))

                    if (.not. associated(PreviousHDF5File%Next)) then
                        !current file is the last file in the list of relevant files
                        call CreateCovHDF5File(HDF5FileAux, HDF5FileX)
                        LastCovHDF5File => HDF5FileAux

                        PreviousHDF5File%Next => LastCovHDF5File
                        LastCovHDF5File%Prev => PreviousHDF5File

                        nullify(LastCovHDF5File%Next)

                        !current file was added to list
                        Number = Number + 1

                        !Adjust end field time if equal to next file start field time
                        if ((LastCovHDF5File%StartFieldTime -                       &
                            PreviousHDF5File%EndFieldTime) .lt. Me%RegularDT) then

                            call AdjustHDF5EndInstant(PreviousHDF5File)

                        endif

                        exit

                    else

                        !check if current file should be located before the next file
                        if (HDF5FileX%StartFieldTime < PreviousHDF5File%Next%StartFieldTime) then
                            !current file should be located before next file

                            !save the previous list begining in PreviousHDF5File%Next
                            LastHDF5File => PreviousHDF5File%Next

                            call CreateCovHDF5File(HDF5FileAux, HDF5FileX)

                            LastHDF5File%Prev => HDF5FileAux
                            HDF5FileAux%Next => LastHDF5File

                            PreviousHDF5File%Next => HDF5FileAux
                            HDF5FileAux%Prev => PreviousHDF5File

                            !current file was added to list
                            Number = Number + 1

                            !Adjust end field time if equal to next file start field time
                            if ((LastHDF5File%StartFieldTime -                      &
                                HDF5FileAux%EndFieldTime) .lt. Me%RegularDT) then

                                call AdjustHDF5EndInstant(HDF5FileAux)

                            endif

                            if ((HDF5FileAux%StartFieldTime -                       &
                                PreviousHDF5File%EndFieldTime) .lt. Me%RegularDT) then

                                call AdjustHDF5EndInstant(PreviousHDF5File)

                            endif

                            exit

                        end if
                    end if

                    !check next file in list
                    PreviousHDF5File => PreviousHDF5File%Next     
                end do
            end if
        end if

    end subroutine AddCovHDF5File

    !--------------------------------------------------------------------------

    subroutine SamePeriodHDF5(HDF5FileX, HDF5FileAux, STAT)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        type(T_HDF5File), pointer                   :: HDF5FileAux
        integer, intent(OUT)                        :: STAT

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        STAT = UNKNOWN_

        HDF5FileAux => Me%FirstHDF5File

        do while (associated(HDF5FileAux))

            if ((HDF5FileX%Name .NE. HDF5FileAux%Name) .and. ((HDF5FileX%HydroProp      &
                .and. HDF5FileAux%HydroProp) .or. (HDF5FileX%WaterProp                  &
                .and. HDF5FileAux%WaterProp))) then 
                !(not the same HDF5 file or different group type)

                !Check if the same period is in more than one file
                if (((HDF5FileX%StartTime >= HDF5FileAux%StartTime)                     &
                    .and. (HDF5FileX%StartTime <= HDF5FileAux%EndTime))                 &
                    .or. ((HDF5FileX%EndTime >= HDF5FileAux%StartTime)                  &
                    .and. (HDF5FileX%EndTime <= HDF5FileAux%EndTime))                   &
                    .or. ((HDF5FileX%StartTime <= HDF5FileAux%StartTime)                &
                    .and. (HDF5FileX%EndTime >= HDF5FileAux%EndTime))) then

                    if ((HDF5FileX%StartTime /= HDF5FileAux%EndTime) .and.              &
                        (HDF5FileAux%StartTime /= HDF5FileX%EndTime)) then

                        STAT = SUCCESS_
                    
                        exit

                    endif
                end if
            end if

            HDF5FileAux => HDF5FileAux%Next

        end do

    end subroutine SamePeriodHDF5

    !--------------------------------------------------------------------------

    subroutine CreateCovHDF5File(HDF5FileNew, HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX
        type(T_HDF5File), pointer                   :: HDF5FileNew

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        !This subroutine atributes the values of the fields of a HDF5File to another HDF5File
        allocate(HDF5FileNew)

        HDF5FileNew%Name             =  HDF5FileX%Name
        HDF5FileNew%StartTime        =  HDF5FileX%StartTime
        HDF5FileNew%EndTime          =  HDF5FileX%EndTime
        HDF5FileNew%StartFieldTime   =  HDF5FileX%StartFieldTime
        HDF5FileNew%EndFieldTime     =  HDF5FileX%EndFieldTime
        HDF5FileNew%Rank             =  HDF5FileX%Rank
        HDF5FileNew%NumberOfInstants =  HDF5FileX%NumberOfInstants
        HDF5FileNew%StartInstant     =  HDF5FileX%StartInstant
        HDF5FileNew%EndInstant       =  HDF5FileX%EndInstant
        HDF5FileNew%InstantsArray    => HDF5FileX%InstantsArray
        HDF5FileNew%FirstInstantTime => HDF5FileX%FirstInstantTime

    end subroutine CreateCovHDF5File

    !--------------------------------------------------------------------------

    subroutine AdjustHDF5EndInstant(HDF5FileX)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX

        !Local-----------------------------------------------------------------
        integer                                     :: HDF5_READ
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        HDF5FileX%EndInstant = HDF5FileX%EndInstant - 1
        
        if (HDF5FileX%EndInstant .lt. HDF5FileX%StartInstant) then

            write(*,*)'HDF5 file with no relevant data:'
            write(*,*) trim(HDF5FileX%Name)
            write(*,*)'Remove correspondent HDF5 block from input file.'
            stop 'AdjustHDF5EndInstant - ModuleAssimilationPreProcessor - ERR01'
        else

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name), HDF5_READ,       &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'AdjustHDF5EndInstant - ModuleAssimilationPreProcessor - ERR02'

            HDF5FileX%EndFieldTime = HDF5TimeInstant(HDF5FileX%EndInstant, HDF5FileX)
            
            call killhdf5(HDF5FileX%HDFID) 

            HDF5FileX%NumberOfInstants = HDF5FileX%NumberOfInstants - 1
        endif

    end subroutine AdjustHDF5EndInstant

    !--------------------------------------------------------------------------

    subroutine AddCovInstantsTimes(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File), pointer                 :: ObjHDF5File

        !Local-----------------------------------------------------------------
        type (T_CovarianceTime), pointer           :: NewTime
        type (T_CovarianceTime), pointer           :: ObjTime           => null()
        type (T_CovarianceTime), pointer           :: PreviousTime      => null()
        type (T_CovarianceTime), pointer           :: ObjInstantTime    => null()

        !Begin-----------------------------------------------------------------

        ObjInstantTime => ObjHDF5File%FirstInstantTime

        do while (associated(ObjInstantTime))

            !Allocates new instance
            allocate (NewTime)
            nullify  (NewTime%Next)

            NewTime%Time = ObjInstantTime%Time

            !Insert New Instance into list and makes Current point to it
            if (.not. associated(Me%FirstCovarianceTime)) then
            !First instance should be the same for all HDF5 files 
                Me%FirstCovarianceTime  => NewTime
                ObjTime                 => NewTime
            else
                PreviousTime            => Me%FirstCovarianceTime
                ObjTime                 => Me%FirstCovarianceTime%Next
                do while (associated(ObjTime))
                    PreviousTime        => ObjTime
                    ObjTime             => ObjTime%Next
                enddo
                ObjTime                 => NewTime
                PreviousTime%Next       => NewTime
            endif

            Me%StatesNumber = Me%StatesNumber + 1

            ObjInstantTime => ObjInstantTime%Next          

        end do 

    end subroutine AddCovInstantsTimes

    !--------------------------------------------------------------------------

    subroutine KillIndividualHDF5File(HDF5ToDispose)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5ToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        deallocate(HDF5ToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'KillIndividualHDF5File - ModuleAssimilationPreProcessor - ERR01'
        nullify(HDF5ToDispose)

    end subroutine KillIndividualHDF5File

    !--------------------------------------------------------------------------   

    subroutine OpenOutputFile

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        real,    dimension(:), pointer              :: TimePtr
        real,    dimension(6), target               :: AuxTime
        type(T_HDF5File), pointer                   :: HDF5FileX
        integer                                     :: neof
        character(len=PathLength)                   :: AuxTSFileName        = null_str    
        character(StringLength)                     :: AuxNum
        integer                                     :: state
        type(T_CovarianceTime), pointer             :: ObjCovarianceTime

        !Begin-----------------------------------------------------------------

        !Get grid values for HDF5 (from first relevant file)
        call ConstructHDF5Grid

        !Calculate number of state variables using HDF5 WaterPoints
        call ConstructStateDimension

        !Gets File Access Code
        call GetHDF5FileAccess (HDF5_CREATE = HDF5_CREATE)

        !Opens output HDF File
        call ConstructHDF5 (Me%ObjCovHDF5, trim(Me%OutputFileName), HDF5_CREATE,        &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR01'

        !Create output HDF5 file
        call ConstructOutputHDF5(Me%ObjCovHDF5)

        !Write Last Time to know the end of covariance period
        if (associated(Me%LastHydroCovHDF5File)) then

            HDF5FileX => Me%LastHydroCovHDF5File
        else

            HDF5FileX => Me%LastWaterCovHDF5File
        endif
 
        call ExtractDate   (HDF5FileX%EndFieldTime,                                     &
                            AuxTime(1), AuxTime(2), AuxTime(3),                         &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%ObjCovHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR02'

        call HDF5WriteData  (Me%ObjCovHDF5, "/Time",                                    &
                             "Time", "YYYY/MM/DD HH:MM:SS",                             &
                             Array1D = TimePtr,                                         &
                             OutputNumber = 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR03'

        call HDF5FlushMemory(Me%ObjCovHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR04'      

        if (Me%StateReconstruction) then
            !Opens state reconstructed output HDF file
            call ConstructHDF5(Me%ObjStateRecHDF5, trim(Me%StateRecFileName),           &
                               HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR05'

            !Create output state HDF5 file
            call ConstructOutputHDF5(Me%ObjStateRecHDF5)

            !Write state times
            ObjCovarianceTime => Me%FirstCovarianceTime
            state = 1

            do while(associated(ObjCovarianceTime))
                call ExtractDate(ObjCovarianceTime%Time,                                &
                                 AuxTime(1), AuxTime(2), AuxTime(3),                    &
                                 AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime

                call HDF5SetLimits  (Me%ObjStateRecHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR06'

                call HDF5WriteData  (Me%ObjStateRecHDF5, "/Time",                       &
                                     "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                     Array1D = TimePtr,                                 &
                                     OutputNumber = state, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR07'

                ObjCovarianceTime => ObjCovarianceTime%Next
                state = state + 1
            enddo

            call HDF5FlushMemory(Me%ObjStateRecHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenOutputFile - ModuleAssimilationPreProcessor - ERR08'      
        endif

        !Kill all grid fields except mapping
        call KillGridFields

        !Open expansion coef. time series files
        allocate(Me%ExpansionCoef(1:Me%StateCovRank))
        
        do neof = 1, Me%StateCovRank
            !Construct time serie name
            write(AuxNum, fmt=*) neof
            AuxTSFileName = "ExpansionCoef"//trim(adjustl(AuxNum))
            Me%ExpansionCoef(neof)%Name = trim(AuxTSFileName)//".ects"

            call StartExpCoefTimeSeries(Me%ExpansionCoef(neof))
        enddo

    end subroutine OpenOutputFile

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Grid

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_HDF5File), pointer                  :: HDF5FileX
        integer                                     :: Rank
        integer, dimension(7)                       :: Dimensions
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        integer                                     :: GridVarNumber, n
        character(len=StringLength)                 :: GridVariableName
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: auxi, auxj, k, i, j

        !Begin-----------------------------------------------------------------

        !Extract from first relevant file: grid, bathymetry, mapping
        !(it is assumed that these are equal in all files)
        if (associated(Me%FirstHydroCovHDF5File)) then
            !(hydrodynamic state is assumed more frequent than strict water properties)
        
            HDF5FileX => Me%FirstHydroCovHDF5File

        else

            HDF5FileX => Me%FirstWaterCovHDF5File

        endif

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (HDF5FileX%HDFID, trim(HDF5FileX%Name),                      &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR01'

        Me%ConnectionX%Name = trim("ConnectionX")
        Me%ConnectionY%Name = trim("ConnectionY")
        Me%Latitude%Name = trim("Latitude")
        Me%Longitude%Name = trim("Longitude")
        Me%Bathymetry%Name = trim("Bathymetry")

        !Get position of grid variables to acess later
        call GetHDF5GroupNumberOfItems(HDF5FileX%HDFID, "/Grid", GridVarNumber,         & 
                                       STAT = STAT_CALL)
        do n = 1, GridVarNumber
            call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                         &
                            n, GridVariableName,                                        &
                            STAT = STAT_CALL)                                
            if (GridVariableName == trim(Me%Bathymetry%Name)) then
                Me%Bathymetry%Position = n
            elseif (GridVariableName == trim(Me%ConnectionX%Name)) then
                Me%ConnectionX%Position = n
            elseif (GridVariableName == trim(Me%ConnectionY%Name)) then
                Me%ConnectionY%Position = n
            elseif (GridVariableName == trim(Me%Latitude%Name)) then
                Me%Latitude%Position = n
            elseif (GridVariableName == trim(Me%Longitude%Name)) then
                Me%Longitude%Position = n
            elseif (GridVariableName == trim(Me%Mapping%Name)) then
                Me%Mapping%Position = n
            endif
        end do

        !Get grid:
        !Mapping dimensions and units
        call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                             &
                            Me%Mapping%Position, trim(Me%Mapping%Name),                 &
                            Me%Mapping%Units, Rank,                                     &
                            Dimensions,                                                 &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_)                                                    &  
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR02'

        !Get file dimensions
        Me%WorkSize%ILB = 1
        Me%WorkSize%IUB = Dimensions(1)
        Me%WorkSize%JLB = 1
        Me%WorkSize%JUB = Dimensions(2)
        if (Me%HDF5_3D) then           
            Me%WorkSize%KLB = 1
        else 
            Me%WorkSize%KLB = 0
        endif
        Me%WorkSize%KUB = Dimensions(3)
        
        Me%Size%ILB = Me%WorkSize%ILB
        Me%Size%JLB = Me%WorkSize%JLB
        Me%Size%IUB = Me%WorkSize%IUB + 1 
        Me%Size%JUB = Me%WorkSize%JUB + 1
        Me%Size%KLB = Me%WorkSize%KLB
        Me%Size%KUB = Me%WorkSize%KUB

        !Allocate variable data matrixes
        nullify(Me%ConnectionX%RealValues2D, Me%ConnectionY%RealValues2D)
        nullify(Me%Longitude%RealValues2D, Me%Latitude%RealValues2D,                    & 
                Me%Bathymetry%RealValues2D)
        
        allocate(Me%ConnectionX%RealValues2D(Me%Size%ILB:Me%Size%IUB,                   &
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ConnectionY%RealValues2D(Me%Size%ILB:Me%Size%IUB,                   & 
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%Longitude%RealValues2D(Me%Size%ILB:Me%Size%IUB,                     & 
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%Latitude%RealValues2D(Me%Size%ILB:Me%Size%IUB,                      & 
                 Me%Size%JLB:Me%Size%JUB))
        allocate(Me%Bathymetry%RealValues2D(Me%Size%ILB:Me%Size%IUB,                    &   
                 Me%Size%JLB:Me%Size%JUB))
                 !(allocate always with size, but read and write with adequate!)
        
        if (Me%HDF5_3D) then 
            nullify(Me%Mapping%IntegerValues3D) 
            allocate(Me%Mapping%IntegerValues3D(Me%Size%ILB:Me%Size%IUB,                &
                     Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
        else 
            nullify(Me%Mapping%IntegerValues2D) !, MapPoints2D)
            allocate(Me%Mapping%IntegerValues2D(Me%Size%ILB:Me%Size%IUB,                &
                     Me%Size%JLB:Me%Size%JUB))
                     !(allocate always with size!)
        endif

        !Get other grid variable units
        call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                             &
                            Me%Bathymetry%Position, trim(Me%Bathymetry%Name),           &
                            Me%Bathymetry%Units,                                        &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_)                                                    &  
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR03'

        call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                             &
                            Me%ConnectionX%Position, trim(Me%ConnectionX%Name),         &
                            Me%ConnectionX%Units,                                       &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_)                                                    &  
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR04'

        call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                             &
                            Me%ConnectionY%Position, trim(Me%ConnectionY%Name),         &
                            Me%ConnectionY%Units,                                       &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_)                                                    &  
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR05'

        call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                             &
                            Me%Longitude%Position, trim(Me%Longitude%Name),             &
                            Me%Longitude%Units,                                         &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_)                                                    &  
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR06'

        call GetHDF5GroupID(HDF5FileX%HDFID, trim("/Grid"),                             &
                            Me%Latitude%Position, trim(Me%Latitude%Name),               &
                            Me%Latitude%Units,                                          &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_)                                                    &  
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR07'

        !Read grid values
        !connections and coordinates
        call HDF5SetLimits (HDF5FileX%HDFID, Me%WorkSize%ILB,                           &
                            Me%WorkSize%IUB+1, Me%WorkSize%JLB,                         &
                            Me%WorkSize%JUB+1, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    & 
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR08'

        call HDF5ReadData(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%ConnectionX%Name),                                    &
                          Array2D      = Me%ConnectionX%RealValues2D,                   &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR09'

        call HDF5ReadData(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%ConnectionY%Name),                                    &
                          Array2D      = Me%ConnectionY%RealValues2D,                   &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR10'

        call HDF5ReadData(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%Latitude%Name),                                       &
                          Array2D      = Me%Latitude%RealValues2D,                      &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR11'

        call HDF5ReadData(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%Longitude%Name),                                      &
                          Array2D      = Me%Longitude%RealValues2D,                     &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR12'

        !bathymetry
        call HDF5SetLimits (HDF5FileX%HDFID, Me%WorkSize%ILB,                           &
                            Me%WorkSize%IUB, Me%WorkSize%JLB,Me%WorkSize%JUB,           &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    & 
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR13'

        call HDF5ReadData(HDF5FileX%HDFID, "/Grid",                                     &
                          trim(Me%Bathymetry%Name),                                     &
                          Array2D      = Me%Bathymetry%RealValues2D,                    &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
        stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR14'

        !mapping
        if (Me%HDF5_3D) then 
            call HDF5SetLimits (HDF5FileX%HDFID, Me%WorkSize%ILB,                       &
                                Me%WorkSize%IUB, Me%WorkSize%JLB,Me%WorkSize%JUB,       &
                                Me%WorkSize%KLB,Me%WorkSize%KUB,                        &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                & 
            stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR15'
            
            call HDF5ReadData(HDF5FileX%HDFID, "/Grid",                                 &
                              trim(Me%Mapping%Name),                                    &
                              Array3D      = Me%Mapping%IntegerValues3D,                &
                              STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR16'

            if (Me%Properties2D) then
                !create a 2D mapping to allow 2D property visualization

                Me%Mapping%AditionalName = trim("MappingPoints2D")

                nullify(Me%Mapping%IntegerValues2D)
                allocate(Me%Mapping%IntegerValues2D(Me%Size%ILB:Me%Size%IUB,            &
                         Me%Size%JLB:Me%Size%JUB))

                Me%Mapping%IntegerValues2D =                                            & 
                            Me%Mapping%IntegerValues3D(Me%WorkSize%ILB:                 &
                            Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB,            &
                            Me%Size%KUB)
                !(assume that relevant mapping is for the upper layer)

                Me%AditionalMap = .true.

                write(*,*) 'Aditional grid variable MappingPoints2D created'
                write(*,*) 'based on upper 3D layer mapping.'

            endif
        else 

            call HDF5ReadData(HDF5FileX%HDFID, "/Grid",                                 &
                              trim(Me%Mapping%Name),                                    &
                              Array2D      = Me%Mapping%IntegerValues2D,                &
                              STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ConstructHDF5Grid - ModuleAssimilationPreProcessor - ERR17'
            
        endif

        call killhdf5(HDF5FileX%HDFID)           

        if (Me%PropertiesInFaces) then
            !Create Me%MappingFacesU 
            nullify(Me%MappingFacesU%IntegerValues3D) 
            allocate(Me%MappingFacesU%IntegerValues3D(Me%Size%ILB:Me%Size%IUB,          &
                     Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
            Me%MappingFacesU%IntegerValues3D(:,:,:) = 0

            !Create Me%MappingFacesV
            nullify(Me%MappingFacesV%IntegerValues3D) 
            allocate(Me%MappingFacesV%IntegerValues3D(Me%Size%ILB:Me%Size%IUB,          &
                     Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
            Me%MappingFacesV%IntegerValues3D(:,:,:) = 0

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB
            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB

            auxj = JLB
            auxi = ILB

            do k = KLB, KUB

                !Me%MappingFacesU
                !Boundary faces (JLB, JUB+1)
                auxj = JLB
                do i = ILB, IUB
                    if (Me%Mapping%IntegerValues3D(i,auxj,k) == 1) then

                        Me%MappingFacesU%IntegerValues3D(i,auxj,k) = 1
                    endif
                enddo
                auxj = JUB
                do i = ILB, IUB
                    if (Me%Mapping%IntegerValues3D(i,auxj,k) == 1) then

                        Me%MappingFacesU%IntegerValues3D(i,auxj + 1,k) = 1
                    endif
                enddo
                !Interior faces
                do j = JLB + 1, JUB
                    do i = ILB, IUB

                        if ((Me%Mapping%IntegerValues3D(i,j-1,k) == 1) .and.            &
                            (Me%Mapping%IntegerValues3D(i,j,k) == 1)) then

                            Me%MappingFacesU%IntegerValues3D(i,j,k) = 1
                        endif
                    enddo
                enddo

                !Me%MappingFacesV
                do j = JLB, JUB
                    !Boundary faces
                    auxi = ILB
                    if (Me%Mapping%IntegerValues3D(auxi,j,k) == 1) then

                        Me%MappingFacesV%IntegerValues3D(auxi,j,k) = 1
                    endif
                    auxi = IUB
                    if (Me%Mapping%IntegerValues3D(auxi,j,k) == 1) then

                        Me%MappingFacesV%IntegerValues3D(auxi + 1,j,k) = 1
                    endif
                    !Interior faces
                    do i = ILB + 1, IUB

                        if ((Me%Mapping%IntegerValues3D(i-1,j,k) == 1) .and.            &
                            (Me%Mapping%IntegerValues3D(i,j,k) == 1)) then

                            Me%MappingFacesV%IntegerValues3D(i,j,k) = 1
                        endif
                    enddo
                enddo

            enddo

        endif

    end subroutine ConstructHDF5Grid

    !--------------------------------------------------------------------------

    subroutine ConstructStateDimension

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: ObjProperty
        integer                                     :: i, j, k
        !integer                                     :: auxi, auxj
        integer                                     :: KLB, KUB, JLB, JUB, ILB, IUB
        integer, dimension(:,:,:), pointer          :: PropertyMap  => null()

        !Begin-----------------------------------------------------------------

        !Obtain state variables according with user-defined window and WaterPoints
        !(only WaterPoints are state variables)
        !(only faces adjacent to two WaterPoints and boundary faces are state variables)

        !This checks if properties have at least one state variable, else an error ocours

        Me%StateVarNumber = 0

        ObjProperty => Me%FirstProperty

        do while(associated(ObjProperty))

            ObjProperty%FirstStatePosition = Me%StateVarNumber + 1

            ObjProperty%StateVarNumber = 0

            if (ObjProperty%Rank == 2) then !2D property
                
                JLB = ObjProperty%Window2D%JLB
                JUB = ObjProperty%Window2D%JUB
                ILB = ObjProperty%Window2D%ILB
                IUB = ObjProperty%Window2D%IUB

                do j = JLB, JUB

                    do i = ILB, IUB

                        if (Me%Mapping%IntegerValues2D(i,j) == WaterPoint) then

                            ObjProperty%StateVarNumber = ObjProperty%StateVarNumber + 1

                            Me%StateVarNumber = Me%StateVarNumber + 1
                        endif
                    enddo
                enddo

            else !3D property

                KLB = ObjProperty%Window%KLB
                KUB = ObjProperty%Window%KUB
                JLB = ObjProperty%Window%JLB
                JUB = ObjProperty%Window%JUB
                ILB = ObjProperty%Window%ILB
                IUB = ObjProperty%Window%IUB

                select case(ObjProperty%TypeZUV)
                    case(TypeZ_)

                        if (.not. ObjProperty%ConvertToFaces) then

                            PropertyMap => Me%Mapping%IntegerValues3D

                        elseif (ObjProperty%ID%IDNumber == VelocityU_) then

                            PropertyMap => Me%MappingFacesU%IntegerValues3D

                            !check if cyclic boundary effective this direction this property
                            if (Me%CyclicBoundary%ON .and.                              &
                                Me%CyclicBoundary%Direction ==                          &
                                DirectionX_ .and. JLB == Me%WorkSize%JLB .and.          &
                                JUB == Me%WorkSize%JUB) then

                                ObjProperty%CyclicBoundary = .true.
                                !(this logical will be required afterwards)
                            endif

                            JUB = JUB + 1

                        else !(assumed velocity V)

                            PropertyMap => Me%MappingFacesV%IntegerValues3D

                            !check if cyclic boundary effective this direction this property
                            if (Me%CyclicBoundary%ON .and.                              &
                                Me%CyclicBoundary%Direction ==                          &
                                DirectionY_ .and. ILB == Me%WorkSize%ILB .and.          &
                                IUB == Me%WorkSize%IUB) then

                                ObjProperty%CyclicBoundary = .true.
                                !(this logical will be required afterwards)
                            endif

                            IUB = IUB + 1
                        
                        endif

                    case(TypeU_)

                        PropertyMap => Me%MappingFacesU%IntegerValues3D

                        !check if cyclic boundary effective this direction this property
                        if (Me%CyclicBoundary%ON .and.                                  &
                            Me%CyclicBoundary%Direction ==                              &
                            DirectionX_ .and. JLB == Me%WorkSize%JLB .and.              &
                            JUB == Me%WorkSize%JUB) then

                            ObjProperty%CyclicBoundary = .true.
                            !(this logical will be required afterwards)
                        endif

                        JUB = JUB + 1

                    case(TypeV_)

                        PropertyMap => Me%MappingFacesV%IntegerValues3D

                        !check if cyclic boundary effective this direction this property
                        if (Me%CyclicBoundary%ON .and.                                  &
                            Me%CyclicBoundary%Direction ==                              &
                            DirectionY_ .and. ILB == Me%WorkSize%ILB .and.              &
                            IUB == Me%WorkSize%IUB) then

                            ObjProperty%CyclicBoundary = .true.
                            !(this logical will be required afterwards)
                        endif

                        IUB = IUB + 1

                end select

                !if (ObjProperty%TypeZUV == TypeZ_ .and.                                 &
                !    .not. ObjProperty%ConvertToFaces) then !Z property

                !    do k = KLB, KUB
                !        do j = JLB, JUB
                !            do i = ILB, IUB

                !                if (Me%Mapping%IntegerValues3D(i,j,k) == 1) then

                !                    ObjProperty%StateVarNumber =                        &
                !                        ObjProperty%StateVarNumber + 1
                !                    Me%StateVarNumber = Me%StateVarNumber + 1
                !                endif
                !            enddo
                !        enddo   
                !    enddo

                !elseif (ObjProperty%TypeZUV == TypeU_ .or.                              &
                !        (ObjProperty%ConvertToFaces .and.                               &
                !        (ObjProperty%ID%IDNumber == VelocityU_))) then !U property

                !    !check if cyclic boundary is effective this direction this property
                !    if (Me%CyclicBoundary%ON .and. Me%CyclicBoundary%Direction ==       &
                !        DirectionX_ .and. JLB == Me%WorkSize%JLB .and.                  &
                !        JUB == Me%WorkSize%JUB) then

                !        ObjProperty%CyclicBoundary = .true.
                !        !(this logical will be required afterwards)
                !    endif

                !    do k = KLB, KUB
                !        !Boundary faces
                !        auxj = JLB
                !        do i = ILB, IUB
                !            if (Me%Mapping%IntegerValues3D(i,auxj,k) == 1) then

                !                ObjProperty%StateVarNumber =                            &
                !                    ObjProperty%StateVarNumber + 1
                !                Me%StateVarNumber = Me%StateVarNumber + 1

                !                !Me%MappingFacesU%IntegerValues3D(i,auxj,k) = 1
                !            endif
                !        enddo
                !        auxj = JUB
                !        do i = ILB, IUB
                !            if (Me%Mapping%IntegerValues3D(i,auxj,k) == 1) then

                !                ObjProperty%StateVarNumber =                            &
                !                    ObjProperty%StateVarNumber + 1
                !                Me%StateVarNumber = Me%StateVarNumber + 1

                !                !Me%MappingFacesU%IntegerValues3D(i,auxj + 1,k) = 1
                !            endif
                !        enddo

                !        !Interior faces
                !        do j = JLB + 1, JUB
                !            do i = ILB, IUB

                !                if ((Me%Mapping%IntegerValues3D(i,j-1,k) == 1) .and.    &
                !                    (Me%Mapping%IntegerValues3D(i,j,k) == 1)) then

                !                    ObjProperty%StateVarNumber =                        &
                !                        ObjProperty%StateVarNumber + 1
                !                    Me%StateVarNumber = Me%StateVarNumber + 1

                !                    !Me%MappingFacesU%IntegerValues3D(i,j,k) = 1
                !                endif
                !            enddo
                !        enddo
                !    enddo

                !else !V property or Z property and convert to faces (velocity V)

                !    !check if cyclic boundary effective in this direction in this property
                !    if (Me%CyclicBoundary%ON .and. Me%CyclicBoundary%Direction ==       &
                !        DirectionY_ .and. ILB == Me%WorkSize%ILB .and.                  &
                !        IUB == Me%WorkSize%IUB) then

                !        ObjProperty%CyclicBoundary = .true.
                !        !(this logical will be required afterwards)
                !    endif

                !    do k = KLB, KUB
                !        do j = JLB, JUB
                !            !Boundary faces
                !            auxi = ILB
                !            if (Me%Mapping%IntegerValues3D(auxi,j,k) == 1) then
                           
                !                ObjProperty%StateVarNumber =                            &
                !                    ObjProperty%StateVarNumber + 1
                !                Me%StateVarNumber = Me%StateVarNumber + 1
                            
                !                !Me%MappingFacesV%IntegerValues3D(auxi,j,k) = 1
                !            endif
                !            auxi = IUB
                !            if (Me%Mapping%IntegerValues3D(auxi,j,k) == 1) then
                           
                !                ObjProperty%StateVarNumber =                            &
                !                    ObjProperty%StateVarNumber + 1
                !                Me%StateVarNumber = Me%StateVarNumber + 1
                            
                !                !Me%MappingFacesV%IntegerValues3D(auxi + 1,j,k) = 1
                !            endif

                !            !Interior faces
                !            do i = ILB + 1, IUB

                !                if ((Me%Mapping%IntegerValues3D(i-1,j,k) == 1) .and.    &
                !                    (Me%Mapping%IntegerValues3D(i,j,k) == 1)) then

                !                    ObjProperty%StateVarNumber =                        &
                !                        ObjProperty%StateVarNumber + 1
                !                    Me%StateVarNumber = Me%StateVarNumber + 1

                !                    !Me%MappingFacesV%IntegerValues3D(i,j,k) = 1
                !                endif
                !            enddo
                !        enddo
                !    enddo
                !endif

                do k = KLB, KUB

                    do j = JLB, JUB

                        do i = ILB, IUB

                            if (PropertyMap(i,j,k) == WaterPoint) then

                                ObjProperty%StateVarNumber = ObjProperty%StateVarNumber + 1

                                Me%StateVarNumber = Me%StateVarNumber + 1
                            endif
                        enddo
                    enddo
                enddo

            endif

            if (ObjProperty%StateVarNumber == 0) then
                !Property without water points!
                write(*,*)
                write(*,*)'Property without water points in user-defined state window!'
                write(*,*)'Property:'//trim(ObjProperty%ID%Name)
                stop 'ConstructStateDimension - ModuleAssimilationPreProcessor - ERR01'
            endif

            ObjProperty => ObjProperty%Next
        enddo

    end subroutine ConstructStateDimension

    !------------------------------------------------------------------------

    subroutine ConstructOutputHDF5(ObjHDF5)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjHDF5

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Set grid
        call HDF5SetLimits(ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB+1,                 &
                           Me%WorkSize%JLB, Me%WorkSize%JUB+1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR02'

        call HDF5WriteData(ObjHDF5, "/Grid", Me%ConnectionX%Name, Me%ConnectionX%Units, &
                           Array2D = Me%ConnectionX%RealValues2D,                       &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR03'

        call HDF5WriteData(ObjHDF5, "/Grid", Me%ConnectionY%Name, Me%ConnectionY%Units, &
                           Array2D = Me%ConnectionY%RealValues2D,                       &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR04'

        call HDF5WriteData(ObjHDF5, "/Grid", Me%Latitude%Name, Me%Latitude%Units,       &
                           Array2D = Me%Latitude%RealValues2D,                          &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR05'

        call HDF5WriteData(ObjHDF5, "/Grid", Me%Longitude%Name, Me%Longitude%Units,     &
                           Array2D = Me%Longitude%RealValues2D,                         &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR06'

        !Set bathymetry          
        call HDF5SetLimits(ObjHDF5, Me%WorkSize%ILB,Me%WorkSize%IUB, Me%WorkSize%JLB,   &
                           Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR07'

        call HDF5WriteData(ObjHDF5, "/Grid", Me%Bathymetry%Name, Me%Bathymetry%Units,   &
                           Array2D      = Me%Bathymetry%RealValues2D,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR08'

        !Set map
        if (Me%HDF5_3D) then 
            call HDF5SetLimits(ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,               &
                               Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                               Me%WorkSize%KLB, Me%WorkSize%KUB,                        &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &      
                stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR09'

            call HDF5WriteData(ObjHDF5, "/Grid", trim(Me%Mapping%Name),                 &
                               trim(Me%Mapping%Units),                                  & 
                               Array3D      = Me%Mapping%IntegerValues3D,               &
                               STAT         = STAT_CALL)      
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR10'

            if (Me%AditionalMap) then

                call HDF5SetLimits(ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,           &
                                   Me%WorkSize%JLB, Me%WorkSize%JUB,                    &
                                   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &      
                stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR11'

                call HDF5WriteData(ObjHDF5, "/Grid", trim(Me%Mapping%AditionalName),    &
                                   trim(Me%Mapping%Units),                              & 
                                   Array2D      = Me%Mapping%IntegerValues2D,           &
                                   STAT         = STAT_CALL)      
                if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR12'

            endif 

        else
            call HDF5WriteData(ObjHDF5, "/Grid", trim(Me%Mapping%Name),                 &
                               trim(Me%Mapping%Units),                                  & 
                               Array2D      = Me%Mapping%IntegerValues2D,               &
                               STAT         = STAT_CALL)    
            if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructOutputHDF5 - ModuleAssimilationPreProcessor - ERR13'
        endif

    end subroutine ConstructOutputHDF5

    !------------------------------------------------------------------------

    subroutine KillGridFields

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Connection X
        deallocate(Me%ConnectionX%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'KillGridFields - ModuleAssimilationPreProcessor - ERR01'  
        nullify(Me%ConnectionX%RealValues2D)

        !Connection Y
        deallocate(Me%ConnectionY%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'KillGridFields - ModuleAssimilationPreProcessor - ERR02'  
        nullify(Me%ConnectionY%RealValues2D)

        !Latitude 
        deallocate(Me%Latitude%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'KillGridFields - ModuleAssimilationPreProcessor - ERR03'  
        nullify(Me%Latitude%RealValues2D)

        !Longitude 
        deallocate(Me%Longitude%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'KillGridFields - ModuleAssimilationPreProcessor - ERR04'  
        nullify(Me%Longitude%RealValues2D)

        !Bathymetry 
        deallocate(Me%Bathymetry%RealValues2D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'KillGridFields - ModuleAssimilationPreProcessor - ERR05'  
        nullify(Me%Bathymetry%RealValues2D)

    end subroutine KillGridFields

    !------------------------------------------------------------------------

    subroutine StartExpCoefTimeSeries(ObjExpansionCoef)

        !Arguments-------------------------------------------------------------
        type (T_ExpansCoefTS)                       :: ObjExpansionCoef

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call AllocateTimeSerieBuffer(ObjExpansionCoef)

        ! Write header for output time serie files
        call UnitsManager (ObjExpansionCoef%ObjTimeSerie, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'StartExpCoefTimeSeries - ModuleAssimilationPreProcessor - ERR01'

        open(ObjExpansionCoef%ObjTimeSerie,                                             &
            file = ObjExpansionCoef%Name, Status = 'unknown')
        
        !write information header
        call WriteDataLine(ObjExpansionCoef%ObjTimeSerie,                               & 
            "Time Serie created by AssimilationPreProcessor.exe")      

        write(ObjExpansionCoef%ObjTimeSerie, *)

        call WriteDataLine(ObjExpansionCoef%ObjTimeSerie, "SERIE_INITIAL_DATA",         &
                           Me%FirstCovarianceTime%Time)        
        call WriteDataLine(ObjExpansionCoef%ObjTimeSerie, "TIME_UNITS", "SECONDS")

        write(ObjExpansionCoef%ObjTimeSerie, *)
        write(ObjExpansionCoef%ObjTimeSerie, *)

        write(ObjExpansionCoef%ObjTimeSerie, fmt=1000)

1000    format(1x,"     Seconds   YY  MM  DD  HH  MM       SS", 3x, "expansion_coef.")

        write(ObjExpansionCoef%ObjTimeSerie, *) '<BeginTimeSerie>'

    end subroutine StartExpCoefTimeSeries

    !--------------------------------------------------------------------------

    subroutine AllocateTimeSerieBuffer(ObjExpansionCoef)

        !Arguments-------------------------------------------------------------
        type (T_ExpansCoefTS)                       :: ObjExpansionCoef

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        ObjExpansionCoef%BufferCount = 0

        !Allocates the TimeSerie Data Buffer
        allocate(ObjExpansionCoef%TimeSerieData(Me%BufferSize), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'AllocateTimeSerieBuffer - ModuleAssimilationPreProcessor - ERR01'
        ObjExpansionCoef%TimeSerieData = null_real

        !Allocates the TimeSerie Time Buffer
        allocate(ObjExpansionCoef%TimeBuffer(Me%BufferSize), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'AllocateTimeSerieBuffer - ModuleAssimilationPreProcessor - ERR02'

    end subroutine AllocateTimeSerieBuffer

    !------------------------------------------------------------------------

    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: ObjProperty
        integer                                     :: ILB, IUB
        integer                                     :: JLB, JUB
        integer                                     :: KLB, KUB

        !Begin-----------------------------------------------------------------

        !Property fields
        ObjProperty => Me%FirstProperty

        do while(associated(ObjProperty))

            ObjProperty%Statistics%RunPeriod         = 0.
            ObjProperty%Statistics%LastCalculation   = Me%CurrentTime

            select case (ObjProperty%Rank)

                case(2)

                    ILB = ObjProperty%Window2D%ILB
                    IUB = ObjProperty%Window2D%IUB
                    JLB = ObjProperty%Window2D%JLB
                    JUB = ObjProperty%Window2D%JUB

                    !Fields allocated for the space window (regardless of WaterPoints)
                    allocate (ObjProperty%Statistics%Average2D (ILB:IUB, JLB:JUB))

                    ObjProperty%Statistics%Average2D = 0.

                    if (Me%Normalization) then

                        allocate (ObjProperty%Statistics%SquareAverage2D                &
                                  (ILB:IUB, JLB:JUB))
                        allocate (ObjProperty%Statistics%StandardDeviation2D            &
                                  (ILB:IUB, JLB:JUB))

                        ObjProperty%Statistics%SquareAverage2D     = 0.
                        ObjProperty%Statistics%StandardDeviation2D = 0.
                    endif

                case(3)
            
                    if (ObjProperty%TypeZUV == TypeZ_ .and.                             &
                        .not. ObjProperty%ConvertToFaces) then !Z property

                        IUB = ObjProperty%Window%IUB
                        JUB = ObjProperty%Window%JUB
                    
                    elseif (ObjProperty%TypeZUV == TypeU_ .or.                          &
                        (ObjProperty%TypeZUV == TypeZ_ .and.                            &
                        ObjProperty%ConvertToFaces .and.                                &
                        ObjProperty%ID%IDNumber == VelocityU_)) then

                        IUB = ObjProperty%Window%IUB
                        JUB = ObjProperty%Window%JUB + 1
                    
                    else !(TypeV_ .or. (TypeZ_ .and. ConvertToFaces .and. VelocityV_)
                        
                        IUB = ObjProperty%Window%IUB + 1
                        JUB = ObjProperty%Window%JUB
                    endif
                    
                    ILB = ObjProperty%Window%ILB
                    JLB = ObjProperty%Window%JLB
                    KLB = ObjProperty%Window%KLB
                    KUB = ObjProperty%Window%KUB

                    !Fields allocated for the space window (regardless of WaterPoints)
                    allocate (ObjProperty%Statistics%Average (ILB:IUB, JLB:JUB, KLB:KUB))

                    ObjProperty%Statistics%Average = 0.

                    if (Me%Normalization) then

                        allocate (ObjProperty%Statistics%SquareAverage                  &
                                  (ILB:IUB, JLB:JUB, KLB:KUB))
                        allocate (ObjProperty%Statistics%StandardDeviation              &
                                  (ILB:IUB, JLB:JUB, KLB:KUB))

                        ObjProperty%Statistics%SquareAverage     = 0.
                        ObjProperty%Statistics%StandardDeviation = 0.
                    endif

            case default 

                stop 'AllocateVariables - ModuleAssimilationPreProcessor - ERR01'

            end select

            ObjProperty => ObjProperty%Next

        enddo

        !Covariance variables
        if (Me%FullCovariance) then
            allocate (Me%StateDeviation (1:Me%StateVarNumber, 1:Me%StatesNumber))
            Me%StateDeviation (:,:) = 0.
        endif

        allocate (Me%AverageState (1:Me%StateVarNumber))
        Me%AverageState (:) = 0.

        allocate (Me%StandardDeviation (1:Me%StateVarNumber))
        Me%StandardDeviation (:) = 0.

        if (Me%FullCovariance) then
            allocate (Me%Covariance (1:Me%StateVarNumber, 1:Me%StateVarNumber))
        else
            allocate (Me%Covariance (1:Me%StatesNumber*(Me%StatesNumber+1)/2,1:1))
        endif
        Me%Covariance (:,:) = 0.

        !State covariance
        if (Me%Method == SEEK) then
            allocate (Me%UVector (1:Me%StateCovRank))
            Me%UVector          (:)   = FillValueReal

            if (.not. Me%FullEOFAnalysis) then
                allocate (Me%LMatrix (1:Me%StateVarNumber, 1:1))
                allocate (Me%StatesLMatrix (1:Me%StatesNumber, 1:Me%StateCovRank))
                Me%LMatrix      (:,:)   = 0.
                Me%StatesLMatrix(:,:)   = 0.
            else
                allocate (Me%LMatrix (1:Me%StateVarNumber, 1:Me%StateCovRank))
                Me%LMatrix      (:,:)   = 0.
            endif
        endif

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine ModifyAssimPreProcessor

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        call StateStatisticsCalculation

        call CovarianceCalculation

        call EOFAnalysis

        if (Me%StateReconstruction) then

            call StateReconstruction

        endif

    end subroutine ModifyAssimPreProcessor

    !--------------------------------------------------------------------------

    subroutine StateStatisticsCalculation

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_CovarianceTime), pointer            :: ObjCovarianceTime
        logical                                     :: FirstFile
        type(T_HDF5File)       , pointer            :: ObjHDF5HydroFile, ObjHDF5WaterFile
        type(T_HDF5File)       , pointer            :: ObjHDF5File
        logical                                     :: Running
        real                                        :: DT
        type (T_Property), pointer                  :: ObjProperty, PropertyToKill
        integer                                     :: Count, i, j
        logical                                     :: LastTime = .false.

        !----------------------------------------------------------------------

        ObjCovarianceTime => Me%FirstCovarianceTime

        !Cycle each relevant HDF5 file
        FirstFile =.true.

        !For each HDF5 file needed for covariance
        ObjHDF5HydroFile => Me%FirstHydroCovHDF5File
        ObjHDF5WaterFile => Me%FirstWaterCovHDF5File

        Me%CurrentStateVar = 0
        Me%CurrentState = 0

        write(*,*)
        write(*,*)'Calculating statistics from historical states...'

        do while(associated(ObjHDF5HydroFile) .or. associated(ObjHDF5WaterFile))
            !if both lists are associated the number of instants is equal 

            !Open and read relevant data
            call OpenAndReadHDF5File(FirstFile, ObjHDF5HydroFile, ObjHDF5WaterFile)

            Me%CurrentTime  = ObjCovarianceTime%Time

            if (FirstFile) then
                FirstFile = .false. !next file is not first

                call ActualizeCurrentTime(Me%ObjTime, Me%RegularDT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                             &
                stop 'StateStatisticsCalculation - ModuleAssimilationPreProcessor - ERR01'
            endif

            !Cycle each Property
            Running      = .true.

            do while (Running)

                Me%CurrentState = Me%CurrentState + 1

                if (Me%CurrentState == (Me%StatesNumber - 1)) then
                    LastTime = .true.
                endif

                !For each Property calculate statistic 
                ObjProperty => Me%FirstProperty

                Me%CurrentStateVar = 0

                Count = Count + 1 

                if (Count .EQ. 1) then 

                    ObjProperty%CurrentField => ObjProperty%FirstField 

                end if

                do while(associated(ObjProperty))

                    select case (ObjProperty%Rank)

                        case(2)

                            call StatisticsCalculation2D(ObjProperty)
 
                        case(3)

                            call StatisticsCalculation3D(ObjProperty)

                    case default 

                stop 'StateStatisticsCalculation - ModuleAssimilationPreProcessor - ERR02'
            
                    end select

                    ObjProperty%CurrentField => ObjProperty%CurrentField%Next               

                    ObjProperty => ObjProperty%Next               

                end do

                if (associated(ObjCovarianceTime%Next)) then

                    DT = ObjCovarianceTime%Next%Time - Me%CurrentTime

                    ObjCovarianceTime => ObjCovarianceTime%Next

                end if

                !Actualization of time            
                Me%CurrentTime = Me%CurrentTime + DT            

                call ActualizeCurrentTime(Me%ObjTime, DT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                             &
                stop 'StateStatisticsCalculation - ModuleAssimilationPreProcessor - ERR03'

                !Running dependent of the last time of file
                if (associated(ObjHDF5HydroFile)) then
                    ObjHDF5File => ObjHDF5HydroFile
                else
                    ObjHDF5File => ObjHDF5WaterFile
                endif

                if (Me%CurrentTime <= ObjHDF5File%EndFieldTime) then
                    Running = .true.
                else
                    Running = .false.
                end if

            end do

            Count = 0

            !Kill field space for each Property 
            ObjProperty => Me%FirstProperty

            do while(associated(ObjProperty))  

                PropertyToKill => ObjProperty
                call KillIndividualPropertyFields(ObjProperty)

                if (LastTime) then
                    call KillIndividualPropertyStats(ObjProperty)
                endif

                ObjProperty    => ObjProperty%Next
            end do

            if (associated(ObjHDF5HydroFile)) then
                ObjHDF5HydroFile => ObjHDF5HydroFile%Next
            endif

            if (associated(ObjHDF5WaterFile)) then
                ObjHDF5WaterFile => ObjHDF5WaterFile%Next
            endif

        end do

        if (Me%FullCovariance) then
            !Calculate deviation from average
            do j = 1, Me%StatesNumber
            do i = 1, Me%StateVarNumber

                Me%StateDeviation(i,j) = Me%StateDeviation(i,j) - Me%AverageState(i)

            end do
            end do
        endif

        if (Me%Normalization) then

            ObjProperty => Me%FirstProperty

            do while(associated(ObjProperty))  
                !Calculate property average standard deviation
                
                ObjProperty%AverStandardDev = 0.

                Count = 0

                do i = ObjProperty%FirstStatePosition,                                  &
                       (ObjProperty%FirstStatePosition - 1) + ObjProperty%StateVarNumber

                    if (Me%StandardDeviation(i) /= 0.) then
                    
                        Count = Count + 1

                        ObjProperty%AverStandardDev = (ObjProperty%AverStandardDev      &
                                                      + Me%StandardDeviation(i))
                    endif

                end do

                ObjProperty%AverStandardDev = ObjProperty%AverStandardDev/real(Count)              

                call WriteNormFactorToOutput(ObjProperty%AverStandardDev,               &
                                             ObjProperty%ID%Name)

                if (Me%FullCovariance) then
                    !Recalculate state deviation
                    if (ObjProperty%AverStandardDev /= 0.) then
                        do j = 1, Me%StatesNumber
                        do i = 1, ObjProperty%FirstStatePosition,                       &
                                  (ObjProperty%FirstStatePosition - 1) +                &
                                  ObjProperty%StateVarNumber

                            Me%StateDeviation(i,j) = Me%StateDeviation(i,j) /           &
                                                     ObjProperty%AverStandardDev

                        end do
                        end do
                    endif
                endif

                ObjProperty    => ObjProperty%Next
            end do

            deallocate(Me%StandardDeviation, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'StateStatisticsCalculation - ModuleSequentialAssimilation - ERR04'
            nullify (Me%StandardDeviation)
        endif

    end subroutine StateStatisticsCalculation

    !--------------------------------------------------------------------------

    subroutine  OpenAndReadHDF5File(FirstFile, ObjHDF5HydroFile, ObjHDF5WaterFile)

        !Arguments-------------------------------------------------------------
        logical                                     :: FirstFile
        type(T_HDF5File), pointer                   :: ObjHDF5HydroFile
        type(T_HDF5File), pointer                   :: ObjHDF5WaterFile

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty
        type(T_HDF5File), pointer                   :: ObjHDF5File
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second 

        !Begin-----------------------------------------------------------------

        !Objective: read relevant data for statistics from file

        !getting parameter fields data
        ObjProperty => Me%FirstProperty

        do while(associated(ObjProperty))

            if (FirstFile) then
            !it is the first file accessed

                nullify(ObjProperty%FirstField)
                !nullify have to be here because each new parameter has to be 
                !nullified for the first file
            
            end if

            !Allocates first field for this parameter
            allocate (ObjProperty%FirstField)
            nullify (ObjProperty%FirstField)

            if (ObjProperty%ModuleType == 1) then

                ObjHDF5File => ObjHDF5HydroFile

            else

                ObjHDF5File => ObjHDF5WaterFile
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (ObjHDF5File%HDFID, trim(ObjHDF5File%Name),              & 
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                write(*,*) 'HDF5 file cannot be opened'//ObjHDF5File%Name                
                stop 'OpenAndReadHDF5File - ModuleAssimilationPreProcessor - ERR01'
            end if

            !Read fields data
            write(*,*)'Reading data from HDF5 file:'//trim(ObjHDF5File%Name) 

            call ReadPropertyFields(ObjHDF5File, ObjProperty)

            call killhdf5(ObjHDF5File%HDFID)

            ObjProperty%CurrentField => ObjProperty%FirstField

            ObjProperty => ObjProperty%Next

        end do

        !check if there are data lacking after the last data from last file 
        if (.not. associated(ObjHDF5File%Next)) then
            !(last file)
            if (ObjHDF5File%EndFieldTime < (Me%EndTime)) then
                !no DT is considered here
100             format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
                call ExtractDate(ObjHDF5File%EndFieldTime, Year, Month, Day, Hour,      &  
                                 Minute, Second)
                write(*,*) 'Data lacking from'      
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write(*,*) 'to'      
                call ExtractDate(Me%EndTime, Year, Month, Day, Hour, Minute, Second)
                write(*,fmt=100) Year, Month, Day, Hour, Minute, Second                    
            end if
        end if

    end subroutine OpenAndReadHDF5File

    !--------------------------------------------------------------------------

    subroutine ReadPropertyFields(ObjHDF5File, ObjProperty)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: ObjHDF5File
        type(T_Property), pointer                   :: ObjProperty
          
        !Local-----------------------------------------------------------------
        integer                                     :: CurrentInstant
        type (T_Field), pointer                     :: NewField
        integer                                     :: Count = 1
        integer                                     :: STAT_CALL
        type (T_Field), pointer                     :: AuxField
        
        !Begin-----------------------------------------------------------------

        if (ObjProperty%ConvertToFaces) then
            !Allocates AuxField
            allocate (AuxField)
            nullify  (AuxField%Next)

            nullify  (AuxField%Values3D)
            allocate (AuxField%Values3D(Me%WorkSize%ILB:Me%WorkSize%IUB,                &
                                        Me%WorkSize%JLB:Me%WorkSize%JUB,                &
                                        Me%WorkSize%KLB:Me%WorkSize%KUB))
        endif

        !read/copy property fields
        write(*,*)'Reading '//trim(ObjProperty%ID%Name)//' fields'

        do CurrentInstant = ObjHDF5File%StartInstant, ObjHDF5File%EndInstant,           &
                            (Me%DecimationFactor + 1)
               
            !construct new fields
            call AddField(ObjProperty%FirstField, NewField)
            NewField%IDNumber = Count
            Count = Count + 1

            !get field ID, Rank and Dimensions
            !(this needs to be done for each instant)
            call GetHDF5GroupID(ObjHDF5File%HDFID, ObjProperty%Group,                   &
                                CurrentInstant, NewField%Name,                          &
                                NewField%Units, ObjProperty%Rank,                       &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
            stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR01'
               
            NewField%Units = trim(NewField%Units)

            !get field values
            select case (ObjProperty%Rank)

                case(2)
                ! The HDF5 file contains 2D data

                    !allocate field
                    nullify (NewField%Values2D)
                    allocate(NewField%Values2D(Me%Size%ILB:Me%Size%IUB,                 &
                                               Me%Size%JLB:Me%Size%JUB))
                       
                    call HDF5SetLimits (ObjHDF5File%HDFID, Me%WorkSize%ILB,             &
                                        Me%WorkSize%IUB, Me%WorkSize%JLB,               &
                                        Me%WorkSize%JUB,                                &
                                        STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR02'
        
                    !read field
                    call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group,             &
                                      trim(NewField%Name),                              &
                                      Array2D      = NewField%Values2D,                 &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR03'

                case(3)
                ! The HDF5 file contains 3D data

                    !allocate field
                    nullify (NewField%Values3D)
                    allocate(NewField%Values3D(Me%Size%ILB:Me%Size%IUB,                 &
                                               Me%Size%JLB:Me%Size%JUB,                 &
                                               Me%Size%KLB:Me%Size%KUB))

                    select case (ObjProperty%TypeZUV)

                        case (TypeZ_)
                            !eventual convert to faces is done afterwards  
                            call HDF5SetLimits  (ObjHDF5File%HDFID, Me%WorkSize%ILB,    &
                                                 Me%WorkSize%IUB, Me%WorkSize%JLB,      &
                                                 Me%WorkSize%JUB,                       &
                                                 Me%WorkSize%KLB,Me%WorkSize%KUB,       &
                                                 STAT = STAT_CALL)                
                            if (STAT_CALL .NE. SUCCESS_)                                & 
                        stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR04'

                            if (.not. ObjProperty%ConvertToFaces) then
                                !read field
                                call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group, &
                                                  trim(NewField%Name),                  &
                                                  Array3D      = NewField%Values3D,     &
                                                  STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                        stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR05'
                            
                            else !ConvertToFaces
                                !read field
                                call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group, &
                                                  trim(NewField%Name),                  &
                                                  Array3D      = AuxField%Values3D,     &
                                                  STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                        stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR06'

                                call ConvertFieldToFaces(AuxField, NewField,            &
                                                         ObjProperty%ID%IDNumber,       &
                                                         ObjProperty%CyclicBoundary)
                            endif

                        case (TypeU_)
                            call HDF5SetLimits  (ObjHDF5File%HDFID, Me%WorkSize%ILB,    &
                                                 Me%WorkSize%IUB, Me%WorkSize%JLB,      &
                                                 Me%WorkSize%JUB + 1,                   &
                                                 Me%WorkSize%KLB,Me%WorkSize%KUB,       &
                                                 STAT = STAT_CALL)                
                            if (STAT_CALL .NE. SUCCESS_)                                & 
                        stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR07'

                            !read field
                            call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group,     &
                                              trim(NewField%Name),                      &
                                              Array3D      = NewField%Values3D,         &
                                              STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR08'

                        case (TypeV_)
                            call HDF5SetLimits  (ObjHDF5File%HDFID, Me%WorkSize%ILB,    &
                                                 Me%WorkSize%IUB + 1, Me%WorkSize%JLB,  &
                                                 Me%WorkSize%JUB,                       &
                                                 Me%WorkSize%KLB,Me%WorkSize%KUB,       &
                                                 STAT = STAT_CALL)                
                            if (STAT_CALL .NE. SUCCESS_)                                & 
                        stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR09'

                            !read field
                            call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group,     &
                                              trim(NewField%Name),                      &
                                              Array3D      = NewField%Values3D,         &
                                              STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR10'

                    end select
 
            case default 
                    
                stop 'ReadPropertyFields - ModuleAssimilationPreProcessor - ERR11'
                    
            end select
        end do

        if (ObjProperty%ConvertToFaces) then
            !Deallocates AuxField
            deallocate(AuxField%Values3D)
            deallocate(AuxField)
        endif

    end subroutine ReadPropertyFields

    !--------------------------------------------------------------------------

    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                     :: FirstField
        type (T_Field), pointer                     :: ObjField
        
        !Local-----------------------------------------------------------------
        type (T_Field), pointer                     :: NewField
        type (T_Field), pointer                     :: PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstField)) then
        !FirstField should be the same for all HDF5 files 
            FirstField              => NewField
            ObjField                => NewField
        else
            PreviousField           => FirstField
            ObjField                => FirstField%Next
            do while (associated(ObjField))
                PreviousField       => ObjField
                ObjField            => ObjField%Next
            enddo
            ObjField                => NewField
            PreviousField%Next      => NewField
        endif

    end subroutine AddField

    !--------------------------------------------------------------------------

    subroutine ConvertFieldToFaces (CenteredField, FacesField, IDNumber, &
                                    CyclicBoundary)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                     :: CenteredField
        type (T_Field), pointer                     :: FacesField
        integer                                     :: IDNumber
        logical                                     :: CyclicBoundary
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: KLB, KUB, JLB, JUB, ILB, IUB

        !Begin-----------------------------------------------------------------
        
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        FacesField%Values3D(:,:,:) = 0.

        select case (IDNumber)

            case(VelocityU_)

                if (CyclicBoundary) then

                    do k = KLB,KUB
                    do i = ILB,IUB

                        if (Me%Mapping%IntegerValues3D(i,JUB-1,k) == WaterPoint .and.   &
                            Me%Mapping%IntegerValues3D(i,JUB  ,k) == WaterPoint) then
                            FacesField%Values3D(i,JLB,k) =                              &
                                (CenteredField%Values3D(i,JUB-1,k) +                    &
                                 CenteredField%Values3D(i,JUB,k))/2.
                        endif
                        if (Me%Mapping%IntegerValues3D(i,JLB,k) == WaterPoint .and.     &
                            Me%Mapping%IntegerValues3D(i,JLB+1,k) == WaterPoint) then
                            FacesField%Values3D(i,JUB+1,k) =                            &
                                (CenteredField%Values3D(i,JLB,k) +                      &
                                 CenteredField%Values3D(i,JLB+1,k))/2.
                        endif
                    
                    enddo
                    enddo

                endif

                do k = KLB,KUB
                do j = JLB + 1,JUB
                do i = ILB,IUB

                    if (Me%Mapping%IntegerValues3D(i,j-1,k) == WaterPoint .and.         &
                        Me%Mapping%IntegerValues3D(i,j  ,k) == WaterPoint) then
                        FacesField%Values3D(i,j,k) =                                    &
                            (CenteredField%Values3D(i,j-1,k) +                          &
                             CenteredField%Values3D(i,j,k))/2.
                    endif
                
                enddo
                enddo
                enddo

            case(VelocityV_)

                if (CyclicBoundary) then

                    do k = KLB,KUB
                    do j = JLB,JUB

                        if (Me%Mapping%IntegerValues3D(IUB-1,j,k) == WaterPoint .and.   &
                            Me%Mapping%IntegerValues3D(IUB,j,k) == WaterPoint) then
                            FacesField%Values3D(ILB,j,k) =                              &
                                (CenteredField%Values3D(IUB-1,j,k) +                    &
                                 CenteredField%Values3D(IUB,j,k))/2.
                        endif
                        if (Me%Mapping%IntegerValues3D(ILB,j,k) == WaterPoint .and.     &
                            Me%Mapping%IntegerValues3D(ILB+1,j,k) == WaterPoint) then
                            FacesField%Values3D(IUB+1,j,k) =                            &
                                (CenteredField%Values3D(ILB,j,k) +                      &
                                 CenteredField%Values3D(ILB+1,j,k))/2.
                        endif

                    enddo
                    enddo

                endif

                do k = KLB,KUB
                do j = JLB,JUB
                do i = ILB + 1,IUB

                    if (Me%Mapping%IntegerValues3D(i,j-1,k) == WaterPoint .and.         &
                        Me%Mapping%IntegerValues3D(i,j  ,k) == WaterPoint) then
                        FacesField%Values3D(i,j,k) =                                    &
                            (CenteredField%Values3D(i,j-1,k) +                          &
                             CenteredField%Values3D(i,j,k))/2.
                    endif

                enddo
                enddo
                enddo

        endselect

    end subroutine ConvertFieldToFaces

    !--------------------------------------------------------------------------

    subroutine StatisticsCalculation2D(ObjProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        real                                        :: DT, DX
        integer                                     :: AuxCurrentStateVar

        !----------------------------------------------------------------------

        !Get limits of property state window
        ILB = ObjProperty%Window2D%ILB
        IUB = ObjProperty%Window2D%IUB
        JLB = ObjProperty%Window2D%JLB
        JUB = ObjProperty%Window2D%JUB
       
        !Time since last calculation 
        DT  = Me%CurrentTime - ObjProperty%Statistics%LastCalculation

cd1:    if (DT>0) then

            AuxCurrentStateVar = Me%CurrentStateVar

            !Loops
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then

                    Me%CurrentStateVar = Me%CurrentStateVar + 1

                    if (ObjProperty%CurrentField%Values2D(i, j) >= FillValueReal / 2.)  &
                        then

                        !Average
                        ObjProperty%Statistics%Average2D (i, j) =                       &
                            (ObjProperty%Statistics%Average2D (i, j) *                  &
                             ObjProperty%Statistics%RunPeriod +                         &
                             ObjProperty%CurrentField%Values2D (i, j) * DT) /           &
                                (ObjProperty%Statistics%RunPeriod + DT)

                        Me%AverageState(Me%CurrentStateVar) =                           &
                            ObjProperty%Statistics%Average2D (i, j)

                        if (Me%FullCovariance) then
                            !Copy state to StateDeviation matrix
                            Me%StateDeviation(Me%CurrentStateVar, Me%CurrentState) =    &
                                ObjProperty%CurrentField%Values2D (i, j)
                        endif
                    endif
                endif
            enddo
            enddo

            if (Me%Normalization) then

                !Loops
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then

                        AuxCurrentStateVar = AuxCurrentStateVar + 1

                        if (ObjProperty%CurrentField%Values2D(i, j) >=                  &
                            FillValueReal / 2.) then

                            !Square Average
                            ObjProperty%Statistics%SquareAverage2D (i, j) =             &
                                (ObjProperty%Statistics%SquareAverage2D (i, j) *        &
                                 ObjProperty%Statistics%RunPeriod +                     &
                                 ObjProperty%CurrentField%Values2D (i, j)**2 * DT) /    &
                                    (ObjProperty%Statistics%RunPeriod + DT)

                            !Standard deviation
                            DX = ObjProperty%Statistics%SquareAverage2D (i, j) -        &
                                 ObjProperty%Statistics%Average2D       (i, j) ** 2

                            DX = abs(DX) 

                            ObjProperty%Statistics%StandardDeviation2D(i, j) = sqrt(DX)

                            Me%StandardDeviation(AuxCurrentStateVar) =                  &
                                ObjProperty%Statistics%StandardDeviation2D(i, j)

                        endif
                    endif
                enddo
                enddo
            endif

            !Updates Time
            ObjProperty%Statistics%RunPeriod       = ObjProperty%Statistics%RunPeriod + DT
            ObjProperty%Statistics%LastCalculation = Me%CurrentTime
      
        endif cd1

    end subroutine StatisticsCalculation2D

    !--------------------------------------------------------------------------

    subroutine StatisticsCalculation3D(ObjProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        real                                        :: DT, DX
        integer                                     :: AuxCurrentStateVar

        !----------------------------------------------------------------------

        !Get limits of property state window (defined as worksize)
        ILB = ObjProperty%Window%ILB
        IUB = ObjProperty%Window%IUB
        JLB = ObjProperty%Window%JLB
        JUB = ObjProperty%Window%JUB
        KLB = ObjProperty%Window%KLB
        KUB = ObjProperty%Window%KUB

        !Time since last calculation 
        DT  = Me%CurrentTime - ObjProperty%Statistics%LastCalculation

cd1:    if (DT>0) then

            AuxCurrentStateVar = Me%CurrentStateVar

            if (ObjProperty%TypeZUV == TypeZ_ .and.                                     &
                .not. ObjProperty%ConvertToFaces) then

                !Loops
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then

                        Me%CurrentStateVar = Me%CurrentStateVar + 1

                        if (ObjProperty%CurrentField%Values3D(i, j, k) >=               &
                            FillValueReal / 2.) then

                            !Average
                            ObjProperty%Statistics%Average (i, j, k) =                  &
                                (ObjProperty%Statistics%Average (i, j, k) *             &
                                 ObjProperty%Statistics%RunPeriod      +                &
                                 ObjProperty%CurrentField%Values3D (i, j, k) * DT) /    &
                                    (ObjProperty%Statistics%RunPeriod + DT)

                            Me%AverageState(Me%CurrentStateVar) =                       &
                                ObjProperty%Statistics%Average (i, j, k)

                            if (Me%FullCovariance) then
                                !Copy state to StateDeviation matrix
                                Me%StateDeviation(Me%CurrentStateVar,Me%CurrentState) = &
                                    ObjProperty%CurrentField%Values3D (i, j, k)
                            endif
                        endif
                    endif
                enddo
                enddo
                enddo

                if (Me%Normalization) then

                    !Loops
                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB
                        if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then

                            AuxCurrentStateVar = AuxCurrentStateVar + 1

                            if (ObjProperty%CurrentField%Values3D(i, j, k) >=           &
                                FillValueReal / 2.) then

                                !Square Average
                                ObjProperty%Statistics%SquareAverage (i, j, k) =        &
                                    (ObjProperty%Statistics%SquareAverage(i, j, k) *    &
                                     ObjProperty%Statistics%RunPeriod      +            &
                                     (ObjProperty%CurrentField%Values3D(i, j, k)**2)*DT)/ &
                                     (ObjProperty%Statistics%RunPeriod + DT)

                                !Standard deviation
                                DX = ObjProperty%Statistics%SquareAverage(i, j, k) -    &
                                     ObjProperty%Statistics%Average(i, j, k) ** 2

                                DX = abs(DX) 

                                ObjProperty%Statistics%StandardDeviation(i, j, k) =     &
                                    sqrt(DX)

                                Me%StandardDeviation(AuxCurrentStateVar) =              &
                                    ObjProperty%Statistics%StandardDeviation(i, j, k)
                            endif
                        endif
                    enddo
                    enddo
                    enddo

                endif

            elseif (ObjProperty%TypeZUV == TypeU_ .or.                                  &
                    (ObjProperty%ConvertToFaces .and.                                   &
                    ObjProperty%ID%IDNumber == VelocityU_)) then

                !Loops
                do k = KLB, KUB
                do j = JLB, JUB + 1
                do i = ILB, IUB
                    if (Me%MappingFacesU%IntegerValues3D(i, j, k) == WaterPoint) then

                        Me%CurrentStateVar = Me%CurrentStateVar + 1

                        if (ObjProperty%CurrentField%Values3D(i, j, k) >=               &
                            FillValueReal / 2.) then

                            !Average
                            ObjProperty%Statistics%Average (i, j, k) =                  &
                                (ObjProperty%Statistics%Average (i, j, k) *             &
                                 ObjProperty%Statistics%RunPeriod      +                &
                                 ObjProperty%CurrentField%Values3D (i, j, k) * DT) /    &
                                    (ObjProperty%Statistics%RunPeriod + DT)

                            Me%AverageState(Me%CurrentStateVar) =                       &
                                ObjProperty%Statistics%Average (i, j, k)

                            if (Me%FullCovariance) then
                                !Copy state to StateDeviation matrix
                                Me%StateDeviation(Me%CurrentStateVar,Me%CurrentState) = &
                                    ObjProperty%CurrentField%Values3D (i, j, k)
                            endif
                        endif
                    endif
                enddo
                enddo
                enddo

                if (Me%Normalization) then

                    !Loops
                    do k = KLB, KUB
                    do j = JLB, JUB + 1
                    do i = ILB, IUB
                        if (Me%MappingFacesU%IntegerValues3D(i,j,k) == WaterPoint) then

                            AuxCurrentStateVar = AuxCurrentStateVar + 1

                            if (ObjProperty%CurrentField%Values3D(i, j, k) >=           &
                                FillValueReal / 2.) then

                                !Square Average
                                ObjProperty%Statistics%SquareAverage (i, j, k) =        &
                                    (ObjProperty%Statistics%SquareAverage(i, j, k) *    &
                                     ObjProperty%Statistics%RunPeriod      +            &
                                     (ObjProperty%CurrentField%Values3D(i,j,k)**2)*DT)  &
                                     / (ObjProperty%Statistics%RunPeriod + DT)

                                !Standard deviation
                                DX = ObjProperty%Statistics%SquareAverage(i, j, k) -    &
                                     ObjProperty%Statistics%Average(i, j, k) ** 2

                                DX = abs(DX) 

                                ObjProperty%Statistics%StandardDeviation(i, j, k) =     &
                                    sqrt(DX)

                                Me%StandardDeviation(AuxCurrentStateVar) =              &
                                    ObjProperty%Statistics%StandardDeviation(i, j, k)
                            endif
                        endif
                    enddo
                    enddo
                    enddo
                endif

            else !Type V

                !Loops
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB + 1
                    if (Me%MappingFacesV%IntegerValues3D(i, j, k) == WaterPoint) then

                        Me%CurrentStateVar = Me%CurrentStateVar + 1

                        if (ObjProperty%CurrentField%Values3D(i, j, k) >=               &
                            FillValueReal / 2.) then

                            !Average
                            ObjProperty%Statistics%Average (i, j, k) =                  &
                                (ObjProperty%Statistics%Average (i, j, k) *             &
                                 ObjProperty%Statistics%RunPeriod      +                &
                                 ObjProperty%CurrentField%Values3D (i, j, k) * DT) /    &
                                    (ObjProperty%Statistics%RunPeriod + DT)

                            Me%AverageState(Me%CurrentStateVar) =                       &
                                ObjProperty%Statistics%Average (i, j, k)

                            if (Me%FullCovariance) then
                                !Copy state to StateDeviation matrix
                                Me%StateDeviation(Me%CurrentStateVar,Me%CurrentState) = &
                                    ObjProperty%CurrentField%Values3D (i, j, k)
                            endif
                        endif
                    endif
                enddo
                enddo
                enddo

                if (Me%Normalization) then

                    !Loops
                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB + 1
                        if (Me%MappingFacesV%IntegerValues3D(i,j,k) == WaterPoint) then

                            AuxCurrentStateVar = AuxCurrentStateVar + 1

                            if (ObjProperty%CurrentField%Values3D(i, j, k) >=           &
                                FillValueReal / 2.) then

                                !Square Average
                                ObjProperty%Statistics%SquareAverage (i, j, k) =        &
                                    (ObjProperty%Statistics%SquareAverage(i, j, k) *    &
                                     ObjProperty%Statistics%RunPeriod      +            &
                                     (ObjProperty%CurrentField%Values3D(i,j,k)**2)*DT)  &
                                     / (ObjProperty%Statistics%RunPeriod + DT)

                                !Standard deviation
                                DX = ObjProperty%Statistics%SquareAverage(i, j, k) -    &
                                     ObjProperty%Statistics%Average(i, j, k) ** 2

                                DX = abs(DX) 

                                ObjProperty%Statistics%StandardDeviation(i, j, k) =     &
                                    sqrt(DX)

                                Me%StandardDeviation(AuxCurrentStateVar) =              &
                                    ObjProperty%Statistics%StandardDeviation(i, j, k)
                            endif
                        endif
                    enddo
                    enddo
                    enddo
                endif

            endif
            !Updates Time
            ObjProperty%Statistics%RunPeriod = ObjProperty%Statistics%RunPeriod + DT
            ObjProperty%Statistics%LastCalculation = Me%CurrentTime
        endif cd1

    end subroutine StatisticsCalculation3D

    !--------------------------------------------------------------------------

    subroutine KillIndividualPropertyFields(PropertyToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Field), pointer                      :: FieldToKill, CurrentField

        !Begin-----------------------------------------------------------------

        CurrentField => PropertyToDispose%FirstField

        do while(associated(CurrentField))

            FieldToKill => CurrentField
            CurrentField => CurrentField%Next

            if (associated(FieldToKill%Values2D)) then
                deallocate(FieldToKill%Values2D, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillIndividualPropertyFields - ModuleAssimilationPreProcessor - ERR01'  
                nullify(FieldToKill%Values2D)
            end if

            if (associated(FieldToKill%Values3D)) then
                deallocate(FieldToKill%Values3D, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'KillIndividualPropertyFields - ModuleAssimilationPreProcessor - ERR02'  
                nullify(FieldToKill%Values3D)
            end if

        end do 

    end subroutine KillIndividualPropertyFields

    !--------------------------------------------------------------------------

    subroutine KillIndividualPropertyStats(ObjProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty

        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer               :: AuxField2D
        real, dimension(:,:,:), pointer             :: AuxField
        real, dimension(:,:,:), pointer             :: AuxFieldFaces
        real, dimension(:,:,:), pointer             :: PropFieldCenter
        type (T_Size3D)                             :: PropFacesSize
        type (T_Size3D)                             :: FacesSize
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        select case (ObjProperty%Rank)

            case(2)

                allocate(AuxField2D(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:   &
                                    Me%WorkSize%JUB))
                AuxField2D(:,:) = FillValueReal
                call SetMatrixValue(AuxField2D, ObjProperty%Window2D,                   &
                                    ObjProperty%Statistics%Average2D)

                !Write average field to output file
                call HDF5SetLimits(Me%ObjCovHDF5, Me%WorkSize%ILB,                      &
                                    Me%WorkSize%IUB, Me%WorkSize%JLB,                   &
                                    Me%WorkSize%JUB, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            & 
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR01'

                call HDF5WriteData(Me%ObjCovHDF5,                                       &
                                   "/Statistics/"//trim(ObjProperty%ID%Name)//"/Average", & 
                                   "Average", ObjProperty%ID%Units,                     &
                                   Array2D = AuxField2D,                                &
                                   OutputNumber = 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR02'            

                !deallocate
                deallocate(ObjProperty%Statistics%Average2D)
                nullify(ObjProperty%Statistics%Average2D)

                if (Me%Normalization) then

                    !deallocate
                    deallocate(ObjProperty%Statistics%SquareAverage2D)
                    nullify(ObjProperty%Statistics%SquareAverage2D)

                    AuxField2D(:,:) = FillValueReal !allocated before for average
                    call SetMatrixValue(AuxField2D, ObjProperty%Window2D,               &
                                        ObjProperty%Statistics%StandardDeviation2D)

                    !Write standard deviation field to output file
                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                                "/Statistics/"//trim(ObjProperty%ID%Name)//"/StandDev", &
                                "StandDev", "-", Array2D = AuxField2D,                  &
                                OutputNumber = 1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR03'            

                    !deallocate
                    deallocate(ObjProperty%Statistics%StandardDeviation2D)
                    nullify(ObjProperty%Statistics%StandardDeviation2D)
                endif

                deallocate(AuxField2D)
                nullify(AuxField2D)

            case(3)
                
                allocate(AuxField(Me%WorkSize%ILB:Me%WorkSize%IUB,                      &
                                  Me%WorkSize%JLB:Me%WorkSize%JUB,                      &
                                  Me%WorkSize%KLB:Me%WorkSize%KUB))
                AuxField(:,:,:) = FillValueReal                     

                if (ObjProperty%TypeZUV == TypeZ_ .and.                                 &
                    .not. ObjProperty%ConvertToFaces) then            

                    call SetMatrixValue(AuxField, ObjProperty%Window,                   &
                                        ObjProperty%Statistics%Average)

                elseif (ObjProperty%TypeZUV == TypeU_ .or.                              &
                        (ObjProperty%ConvertToFaces .and.                               &
                        ObjProperty%ID%IDNumber == VelocityU_)) then

                    allocate(PropFieldCenter(                                           &
                             ObjProperty%Window%ILB:ObjProperty%Window%IUB,             &
                             ObjProperty%Window%JLB:ObjProperty%Window%JUB,             &
                             ObjProperty%Window%KLB:ObjProperty%Window%KUB))
                    PropFieldCenter(:, :, :) = FillValueReal                     

                    call CenterField(PropFieldCenter,                                   &
                                     ObjProperty%Statistics%Average,                    &
                                     ObjProperty%ID%IDNumber, ObjProperty%Window)

                    call SetMatrixValue(AuxField, ObjProperty%Window, PropFieldCenter)

                    !(write with Size)
                    FacesSize%ILB = Me%WorkSize%ILB
                    FacesSize%IUB = Me%WorkSize%IUB + 1
                    FacesSize%JLB = Me%WorkSize%JLB
                    FacesSize%JUB = Me%WorkSize%JUB + 1
                    FacesSize%KLB = Me%WorkSize%KLB
                    FacesSize%KUB = Me%WorkSize%KUB

                    PropFacesSize%ILB = ObjProperty%Window%ILB
                    PropFacesSize%IUB = ObjProperty%Window%IUB
                    PropFacesSize%JLB = ObjProperty%Window%JLB
                    PropFacesSize%JUB = ObjProperty%Window%JUB + 1
                    PropFacesSize%KLB = ObjProperty%Window%KLB
                    PropFacesSize%KUB = ObjProperty%Window%KUB

                    allocate(AuxFieldFaces(FacesSize%ILB:FacesSize%IUB,                 &
                                           FacesSize%JLB:FacesSize%JUB,                 &
                                           FacesSize%KLB:FacesSize%KUB))
                    AuxFieldFaces(:,:,:) = FillValueReal

                    call SetMatrixValue(AuxFieldFaces, PropFacesSize,                   &
                                        ObjProperty%Statistics%Average)

                    !write faces average
                    call HDF5SetLimits(Me%ObjCovHDF5, FacesSize%ILB,                    &
                                       FacesSize%IUB, FacesSize%JLB,                    &
                                       FacesSize%JUB, FacesSize%KLB,                    &
                                       FacesSize%KUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR04'

                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                         "/Statistics/"//trim(ObjProperty%ID%Name)//"/AverageFaces",    & 
                                       "AverageFaces",                                  &
                                       ObjProperty%ID%Units,                            &
                                       Array3D = AuxFieldFaces,                         &
                                       OutputNumber = 1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR05' 

                    !deallocate
                    deallocate(AuxFieldFaces)
                    nullify(AuxFieldFaces)

                else !TypeV_ or ObjProperty%ConvertToFaces and VelocityV_

                    allocate(PropFieldCenter(                                           &
                             ObjProperty%Window%ILB:ObjProperty%Window%IUB,             &
                             ObjProperty%Window%JLB:ObjProperty%Window%JUB,             &
                             ObjProperty%Window%KLB:ObjProperty%Window%KUB))
                    PropFieldCenter(:, :, :) = FillValueReal                      

                    call CenterField(PropFieldCenter,                                   &
                                     ObjProperty%Statistics%Average,                    &
                                     ObjProperty%ID%IDNumber, ObjProperty%Window)

                    call SetMatrixValue(AuxField, ObjProperty%Window, PropFieldCenter)

                    !(write with Size)
                    FacesSize%ILB = Me%WorkSize%ILB
                    FacesSize%IUB = Me%WorkSize%IUB + 1
                    FacesSize%JLB = Me%WorkSize%JLB
                    FacesSize%JUB = Me%WorkSize%JUB + 1
                    FacesSize%KLB = Me%WorkSize%KLB
                    FacesSize%KUB = Me%WorkSize%KUB

                    PropFacesSize%ILB = ObjProperty%Window%ILB
                    PropFacesSize%IUB = ObjProperty%Window%IUB + 1
                    PropFacesSize%JLB = ObjProperty%Window%JLB
                    PropFacesSize%JUB = ObjProperty%Window%JUB
                    PropFacesSize%KLB = ObjProperty%Window%KLB
                    PropFacesSize%KUB = ObjProperty%Window%KUB

                    allocate(AuxFieldFaces(FacesSize%ILB:FacesSize%IUB,                 &
                                           FacesSize%JLB:FacesSize%JUB,                 &
                                           FacesSize%KLB:FacesSize%KUB))
                    AuxFieldFaces(:,:,:) = FillValueReal

                    call SetMatrixValue(AuxFieldFaces, PropFacesSize,                   &
                                        ObjProperty%Statistics%Average)

                    !write faces average
                    call HDF5SetLimits(Me%ObjCovHDF5, FacesSize%ILB,                    &
                                       FacesSize%IUB, FacesSize%JLB,                    &
                                       FacesSize%JUB, FacesSize%KLB,                    &
                                       FacesSize%KUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR06'

                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                         "/Statistics/"//trim(ObjProperty%ID%Name)//"/AverageFaces",    & 
                                       "AverageFaces",                                  &
                                       ObjProperty%ID%Units,                            &
                                       Array3D = AuxFieldFaces,                         &
                                       OutputNumber = 1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR07' 

                    !deallocate
                    deallocate(AuxFieldFaces)
                    nullify(AuxFieldFaces)
                endif

                !deallocate
                deallocate(ObjProperty%Statistics%Average)
                nullify(ObjProperty%Statistics%Average)

                !Write average field to output file
                call HDF5SetLimits(Me%ObjCovHDF5, Me%WorkSize%ILB,                      &
                                   Me%WorkSize%IUB, Me%WorkSize%JLB,                    &
                                   Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB,   &
                                   STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            & 
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR08'

                call HDF5WriteData(Me%ObjCovHDF5,                                       &
                                   "/Statistics/"//trim(ObjProperty%ID%Name)//"/Average", & 
                                   "Average", ObjProperty%ID%Units,                     &
                                   Array3D = AuxField,                                  &
                                   OutputNumber = 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR09' 

                if (Me%Normalization) then

                    !deallocate
                    deallocate(ObjProperty%Statistics%SquareAverage)
                    nullify(ObjProperty%Statistics%SquareAverage)

                    AuxField(:,:,:) = FillValueReal !allocated before for average

                    if (ObjProperty%TypeZUV == TypeZ_ .and.                             &
                        .not. ObjProperty%ConvertToFaces) then            

                        call SetMatrixValue(AuxField, ObjProperty%Window,               &
                                            ObjProperty%Statistics%StandardDeviation)

                    elseif (ObjProperty%TypeZUV == TypeU_ .or.                          &
                            (ObjProperty%ConvertToFaces .and.                           &
                            ObjProperty%ID%IDNumber == VelocityU_)) then

                        PropFieldCenter(:,:,:) = FillValueReal !allocated before

                        call CenterField(PropFieldCenter,                               &
                                         ObjProperty%Statistics%StandardDeviation,      &
                                         ObjProperty%ID%IDNumber, ObjProperty%Window)

                        call SetMatrixValue(AuxField, ObjProperty%Window,               &
                                            PropFieldCenter)

                    else !TypeV_ or ObjProperty%ConvertToFaces and VelocityV_

                        PropFieldCenter(:,:,:) = FillValueReal !allocated before

                        call CenterField(PropFieldCenter,                               &
                                         ObjProperty%Statistics%StandardDeviation,      &
                                         ObjProperty%ID%IDNumber, ObjProperty%Window)

                        call SetMatrixValue(AuxField, ObjProperty%Window,               &
                                            PropFieldCenter)
                    endif

                    !deallocate
                    deallocate(ObjProperty%Statistics%StandardDeviation)
                    nullify(ObjProperty%Statistics%StandardDeviation)

                    !Write standard deviation field to output file
                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                                "/Statistics/"//trim(ObjProperty%ID%Name)//"/StandDev", &
                                "StandDev","-",Array3D = AuxField,                      &
                                OutputNumber = 1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR10'            

                endif

                deallocate(AuxField)
                nullify(AuxField)

                if (associated(PropFieldCenter)) then
                    deallocate(PropFieldCenter)
                    nullify(PropFieldCenter)
                endif
        end select

        call HDF5FlushMemory(Me%ObjCovHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'KillIndividualPropertyStats - ModuleAssimilationPreProcessor - ERR11'

    end subroutine KillIndividualPropertyStats

    !--------------------------------------------------------------------------

    subroutine WriteNormFactorToOutput(PropStandDevValue, PropName)

        !Arguments-------------------------------------------------------------
        real                                            :: PropStandDevValue
        character(len=StringLength)                     :: PropName

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer                     :: AuxValue
        integer                                         :: STAT_CALL

        !----------------------------------------------------------------------

        allocate(AuxValue(1:1))

        !Write normalization factor
        call HDF5SetLimits(Me%ObjCovHDF5, 1, 1, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                        & 
        stop 'WriteNormFactorToOutput - ModuleAssimilationPreProcessor - ERR01'

        AuxValue(1) = 1./PropStandDevValue

        call HDF5WriteData(Me%ObjCovHDF5,                                   &
                           "/AssimilationData/NormFactor/"                  &
                           //trim(PropName),                                & 
                           PropName, "-", Array1D = AuxValue,               &
                           OutputNumber = 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
        stop 'WriteNormFactorToOutput - ModuleAssimilationPreProcessor - ERR02'  

        deallocate(AuxValue)

    end subroutine WriteNormFactorToOutput

    !--------------------------------------------------------------------------

    subroutine CovarianceCalculation

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i

        !----------------------------------------------------------------------

        if (Me%FullCovariance) then

            write(*,*)
            write(*,*) 'Calculating covariance matrix...'

            !Calculate covariance
            if (Me%FullEOFAnalysis) then
                Me%Covariance = MATMUL(Me%StateDeviation, TRANSPOSE(Me%StateDeviation)) &
                                /(real(Me%StatesNumber))
            else
                Me%Covariance = MATMUL(TRANSPOSE(Me%StateDeviation), Me%StateDeviation) &
                                /(real(Me%StatesNumber))
            endif

            deallocate(Me%StateDeviation)
            nullify(Me%StateDeviation)

            Me%CovTrace = 0.

            !Calculate matrix trace
            do i = 1, Me%StateVarNumber

                Me%CovTrace = Me%CovTrace + Me%Covariance(i,i)
            end do
        else
            call StatesCovarianceCalculation
        endif

    end subroutine CovarianceCalculation

    !--------------------------------------------------------------------------

    subroutine StatesCovarianceCalculation

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, ii
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: statei, statej
        real(8)                                     :: s, tmp
        type(T_HDF5File)        , pointer           :: ObjHDF5HydroFile, ObjHDF5WaterFile
        type(T_HDF5File)        , pointer           :: ObjHDF5File, AuxObjHDF5File
        type(T_HDF5File)        , pointer           :: AuxObjHDF5HydroFile
        type(T_HDF5File)        , pointer           :: AuxObjHDF5WaterFile
        logical                                     :: LastFile
        real, dimension(:)      , pointer           :: AuxVector
        integer                                     :: CurrentInstant, AuxCurrentInstant
        type(T_Property)        , pointer           :: ObjProperty
        real                                        :: AuxStandardDev

        !----------------------------------------------------------------------

        !StatesCovarianceCalculation calculates Me%Covariance as a 1D column vector
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output,
        !
        ! Local parameters:
        !
        !

        allocate(AuxVector(1:Me%StateVarNumber))

        !Cycle each relevant HDF5 file
        LastFile =.true.

        !For each HDF5 file needed for covariance
        ObjHDF5HydroFile => Me%LastHydroCovHDF5File
        ObjHDF5WaterFile => Me%LastWaterCovHDF5File
        
        if (Me%NumberHydroCovHDF5 > 0) then
            ObjHDF5File => ObjHDF5HydroFile
        else
            ObjHDF5File => ObjHDF5WaterFile
        endif
        
        CurrentInstant = ObjHDF5File%EndInstant

        Me%CurrentState = 0

        ii = Me%StatesNumber*(Me%StatesNumber + 1)/2
        Me%CovTrace = 0.

        write(*,*)
        write(*,*)'Calculating covariance from historical states...'

do1:    do statei = Me%StatesNumber, 1, -1

            write(*,*) 'Processing state: ', statei

            s = 0.
            Me%CurrentStateVar = 0

            if (LastFile .and. (statei == Me%StatesNumber)) then
                !Read property fields from file
                call OpenAndReadHDF5OneInstant(CurrentInstant, ObjHDF5HydroFile,        &
                                               ObjHDF5WaterFile)
                LastFile = .false.
            endif

            !For each Property
            ObjProperty => Me%FirstProperty

            do while(associated(ObjProperty))

                if (Me%Normalization) then
                    AuxStandardDev = ObjProperty%AverStandardDev
                else
                    AuxStandardDev = 1
                endif

                !Cycle the property's variables
                select case (ObjProperty%Rank)

                    case(2)
                        !Get limits of property state window
                        ILB = ObjProperty%Window2D%ILB
                        IUB = ObjProperty%Window2D%IUB
                        JLB = ObjProperty%Window2D%JLB
                        JUB = ObjProperty%Window2D%JUB
                        
                        do j = JLB, JUB
                        do i = ILB, IUB
                            if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then

                                Me%CurrentStateVar = Me%CurrentStateVar + 1

                                if (ObjProperty%FirstField%Values2D(i, j)               &
                                    >= FillValueReal / 2.) then

                                    tmp = (ObjProperty%FirstField%Values2D (i, j) -     &
                                           Me%AverageState(Me%CurrentStateVar))/        &
                                           AuxStandardDev
                                    AuxVector(Me%CurrentStateVar) = tmp
                                    s        = s + tmp**2
                                endif
                            endif
                        enddo
                        enddo

                        deallocate(ObjProperty%FirstField%Values2D)
                        nullify(ObjProperty%FirstField%Values2D)

                    case(3)
                        !Get limits of property state window
                        ILB = ObjProperty%Window%ILB
                        IUB = ObjProperty%Window%IUB
                        JLB = ObjProperty%Window%JLB
                        JUB = ObjProperty%Window%JUB
                        KLB = ObjProperty%Window%KLB
                        KUB = ObjProperty%Window%KUB

                        if (ObjProperty%TypeZUV == TypeZ_ .and.                         &
                            .not. ObjProperty%ConvertToFaces) then

                            do k = KLB, KUB
                            do j = JLB, JUB
                            do i = ILB, IUB
                                if (Me%Mapping%IntegerValues3D(i,j,k) == WaterPoint)    &
                                    then
                                    Me%CurrentStateVar = Me%CurrentStateVar + 1

                                    if (ObjProperty%FirstField%Values3D(i, j, k)        &
                                        >= FillValueReal / 2.) then

                                        tmp = (ObjProperty%FirstField%                  &
                                               Values3D(i,j,k) -                        &
                                            Me%AverageState(Me%CurrentStateVar))/       &
                                               AuxStandardDev
                                        AuxVector(Me%CurrentStateVar) = tmp
                                        s        = s + tmp**2
                                    endif
                                endif
                            enddo
                            enddo
                            enddo

                        elseif (ObjProperty%TypeZUV == TypeU_ .or.                      &
                                (ObjProperty%ConvertToFaces .and.                       &
                                ObjProperty%ID%IDNumber == VelocityU_)) then

                            do k = KLB, KUB
                            do j = JLB, JUB + 1
                            do i = ILB, IUB
                                if (Me%MappingFacesU%IntegerValues3D(i,j,k) ==          &
                                    WaterPoint) then
                                    Me%CurrentStateVar = Me%CurrentStateVar + 1

                                    if (ObjProperty%FirstField%Values3D(i, j, k)        &
                                        >= FillValueReal / 2.) then

                                        tmp = (ObjProperty%FirstField%                  &
                                               Values3D(i,j,k) -                        &
                                            Me%AverageState(Me%CurrentStateVar))/       &
                                               AuxStandardDev
                                        AuxVector(Me%CurrentStateVar) = tmp
                                        s        = s + tmp**2
                                    endif
                                endif
                            enddo
                            enddo
                            enddo

                        else !Type V .or. ConvertToFaces .and. VelocityV_

                            do k = KLB, KUB
                            do j = JLB, JUB
                            do i = ILB, IUB + 1
                                if (Me%MappingFacesV%IntegerValues3D(i,j,k) ==          &
                                    WaterPoint) then
                                    Me%CurrentStateVar = Me%CurrentStateVar + 1

                                    if (ObjProperty%FirstField%Values3D(i, j, k)        &
                                        >= FillValueReal / 2.) then

                                        tmp = (ObjProperty%FirstField%                  &
                                               Values3D(i,j,k) -                        &
                                            Me%AverageState(Me%CurrentStateVar))/       &
                                               AuxStandardDev
                                        AuxVector(Me%CurrentStateVar) = tmp
                                        s        = s + tmp**2
                                    endif
                                endif
                            enddo
                            enddo
                            enddo
                        endif

                        deallocate(ObjProperty%FirstField%Values3D)
                        nullify(ObjProperty%FirstField%Values3D)

                end select

                ObjProperty => ObjProperty%Next               
            end do

            tmp                     = s/real(Me%StatesNumber)
            Me%Covariance(ii,1)     = tmp
            Me%CovTrace             = Me%CovTrace + tmp

            ii = ii - statei

            AuxObjHDF5HydroFile => Me%FirstHydroCovHDF5File
            AuxObjHDF5WaterFile => Me%FirstWaterCovHDF5File
        
            if (Me%NumberHydroCovHDF5 > 0) then
                AuxObjHDF5File => AuxObjHDF5HydroFile
            else
                AuxObjHDF5File => AuxObjHDF5WaterFile
            endif

            AuxCurrentInstant = AuxObjHDF5File%StartInstant

do2:        do statej = 1, statei-1

                call OpenAndReadHDF5OneInstant(AuxCurrentInstant, AuxObjHDF5HydroFile,  &
                                               AuxObjHDF5WaterFile)

                s = 0.
                Me%CurrentStateVar = 0

                !For each Property
                ObjProperty => Me%FirstProperty

                do while(associated(ObjProperty))

                    if (Me%Normalization) then
                        AuxStandardDev = ObjProperty%AverStandardDev
                    else
                        AuxStandardDev = 1
                    endif

                    !Cycle the property's variables
                    select case (ObjProperty%Rank)

                        case(2)
                            !Get limits of property state window
                            ILB = ObjProperty%Window2D%ILB
                            IUB = ObjProperty%Window2D%IUB
                            JLB = ObjProperty%Window2D%JLB
                            JUB = ObjProperty%Window2D%JUB

                            do j = JLB, JUB
                            do i = ILB, IUB
                                if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then

                                    Me%CurrentStateVar = Me%CurrentStateVar + 1

                                    if (ObjProperty%FirstField%Values2D(i, j)           &
                                        >= FillValueReal / 2.) then

                                        s = s + ((ObjProperty%FirstField%Values2D (i,j) &
                                            - Me%AverageState(Me%CurrentStateVar))*     &
                                            AuxVector(Me%CurrentStateVar))/             &
                                            AuxStandardDev
                                    endif
                                endif
                            enddo
                            enddo

                            if (statej /= statei-1) then
                                deallocate(ObjProperty%FirstField%Values2D)
                                nullify(ObjProperty%FirstField%Values2D)
                            endif
                        case(3)
                            !Get limits of property state window
                            ILB = ObjProperty%Window%ILB
                            IUB = ObjProperty%Window%IUB
                            JLB = ObjProperty%Window%JLB
                            JUB = ObjProperty%Window%JUB
                            KLB = ObjProperty%Window%KLB
                            KUB = ObjProperty%Window%KUB

                            if (ObjProperty%TypeZUV == TypeZ_ .and.                     &
                                .not. ObjProperty%ConvertToFaces) then

                                do k = KLB, KUB
                                do j = JLB, JUB
                                do i = ILB, IUB
                                    if (Me%Mapping%IntegerValues3D(i, j, k) ==          &
                                        WaterPoint) then

                                        Me%CurrentStateVar = Me%CurrentStateVar + 1

                                        if (ObjProperty%FirstField%Values3D(i, j, k)    &
                                            >= FillValueReal / 2.) then

                                            s = s + ((ObjProperty%FirstField%           &
                                                Values3D(i,j,k) -                       &
                                            Me%AverageState(Me%CurrentStateVar))*       &
                                                AuxVector(Me%CurrentStateVar))/         &
                                                AuxStandardDev
                                        endif
                                    endif
                                enddo
                                enddo
                                enddo

                            elseif (ObjProperty%TypeZUV == TypeU_ .or.                  &
                                    (ObjProperty%ConvertToFaces .and.                   &
                                    ObjProperty%ID%IDNumber == VelocityU_)) then

                                do k = KLB, KUB
                                do j = JLB, JUB + 1
                                do i = ILB, IUB
                                    if (Me%MappingFacesU%IntegerValues3D(i,j,k) ==      &
                                        WaterPoint) then
                                        Me%CurrentStateVar = Me%CurrentStateVar + 1

                                        if (ObjProperty%FirstField%Values3D(i, j, k)    &
                                            >= FillValueReal / 2.) then

                                            s = s + ((ObjProperty%FirstField%           &
                                                Values3D(i,j,k) -                       &
                                            Me%AverageState(Me%CurrentStateVar))*       &
                                                AuxVector(Me%CurrentStateVar))/         &
                                                AuxStandardDev
                                        endif
                                    endif
                                enddo
                                enddo
                                enddo

                            else !Type V .or. ConvertToFaces .and. VelocityV_

                                do k = KLB, KUB
                                do j = JLB, JUB
                                do i = ILB, IUB + 1
                                    if (Me%MappingFacesV%IntegerValues3D(i,j,k) ==      &
                                        WaterPoint) then
                                        Me%CurrentStateVar = Me%CurrentStateVar + 1

                                        if (ObjProperty%FirstField%Values3D(i, j, k)    &
                                            >= FillValueReal / 2.) then

                                            s = s + ((ObjProperty%FirstField%           &
                                                Values3D(i,j,k) -                       &
                                            Me%AverageState(Me%CurrentStateVar))*       &
                                                AuxVector(Me%CurrentStateVar))/         &
                                                AuxStandardDev
                                        endif
                                    endif
                                enddo
                                enddo
                                enddo
                            endif
                            if (statej /= statei-1) then
                                deallocate(ObjProperty%FirstField%Values3D)
                                nullify(ObjProperty%FirstField%Values3D)
                            endif
                    end select

                    ObjProperty => ObjProperty%Next               
                end do

                Me%Covariance(ii+statej,1)   = s/real(Me%StatesNumber)

                AuxCurrentInstant = AuxCurrentInstant + (Me%DecimationFactor + 1)
                if (AuxCurrentInstant > AuxObjHDF5File%EndInstant) then

                    if (Me%NumberHydroCovHDF5 > 0) then
                        AuxObjHDF5HydroFile => AuxObjHDF5HydroFile%Next
                        AuxObjHDF5File      => AuxObjHDF5HydroFile
                    endif

                    if (Me%NumberWaterCovHDF5 > 0) then
                        AuxObjHDF5WaterFile => AuxObjHDF5WaterFile%Next
                        if (Me%NumberHydroCovHDF5 == 0) then
                            AuxObjHDF5File  => AuxObjHDF5WaterFile
                        endif
                    endif

                    if (.not. (associated(AuxObjHDF5File))) exit do2

                    AuxCurrentInstant = AuxObjHDF5File%StartInstant
                endif
            enddo do2

            CurrentInstant = CurrentInstant - (Me%DecimationFactor + 1)
            if (CurrentInstant < ObjHDF5File%StartInstant) then

                if (Me%NumberHydroCovHDF5 > 0) then
                    ObjHDF5HydroFile => ObjHDF5HydroFile%Prev
                    ObjHDF5File      => ObjHDF5HydroFile
                endif

                if (Me%NumberWaterCovHDF5 > 0) then
                    ObjHDF5WaterFile => ObjHDF5WaterFile%Prev
                    if (Me%NumberHydroCovHDF5 == 0) then
                        ObjHDF5File  => ObjHDF5WaterFile
                    endif
                endif

                if (.not. (associated(ObjHDF5File))) exit do1

                CurrentInstant = ObjHDF5File%EndInstant
            endif
        enddo do1

        deallocate(AuxVector)

    end subroutine StatesCovarianceCalculation

    !--------------------------------------------------------------------------

    subroutine OpenAndReadHDF5OneInstant(Instant, ObjHDF5HydroFile, ObjHDF5WaterFile)

        !Arguments-------------------------------------------------------------
        integer                                     :: Instant
        type(T_HDF5File), pointer                   :: ObjHDF5HydroFile
        type(T_HDF5File), pointer                   :: ObjHDF5WaterFile

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty
        type(T_HDF5File), pointer                   :: ObjHDF5File
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        type (T_Field), pointer                     :: AuxField

        !Begin-----------------------------------------------------------------

        !getting parameter fields data
        ObjProperty => Me%FirstProperty

        do while(associated(ObjProperty))

            if ((ObjProperty%ModuleType == 1) .or. (Me%NumberWaterCovHDF5 == 0)) then

                ObjHDF5File => ObjHDF5HydroFile
            else

                ObjHDF5File => ObjHDF5WaterFile
            endif

            if (ObjProperty%ConvertToFaces) then
                !Allocates AuxField
                allocate (AuxField)
                nullify  (AuxField%Next)

                nullify  (AuxField%Values3D)
                allocate (AuxField%Values3D(Me%WorkSize%ILB:Me%WorkSize%IUB,            &
                                            Me%WorkSize%JLB:Me%WorkSize%JUB,            &
                                            Me%WorkSize%KLB:Me%WorkSize%KUB))
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            !Open HDF5 file
            call ConstructHDF5 (ObjHDF5File%HDFID, trim(ObjHDF5File%Name),              & 
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                write(*,*) 'HDF5 file cannot be opened'//ObjHDF5File%Name                
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR01'
            end if

            call GetHDF5GroupID(ObjHDF5File%HDFID, ObjProperty%Group, Instant,          &
                                ObjProperty%FirstField%Name, STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                                &  
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR02'

            !get field values
            select case (ObjProperty%Rank)

                case(2)
                    ! The HDF5 file contains 2D data
                    allocate(ObjProperty%FirstField%Values2D(Me%Size%ILB:Me%Size%IUB,   &
                         Me%Size%JLB:Me%Size%JUB))
                     
                    call HDF5SetLimits (ObjHDF5File%HDFID, Me%WorkSize%ILB,             &
                                        Me%WorkSize%IUB, Me%WorkSize%JLB,               &
                                        Me%WorkSize%JUB,                                &
                                        STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR03'
        
                    !read field
                    call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group,             &
                                      trim(ObjProperty%FirstField%Name),                &
                                      Array2D = ObjProperty%FirstField%Values2D,        &
                                      STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR04'

                case(3)
                    ! The HDF5 file contains 3D data
                    allocate(ObjProperty%FirstField%Values3D(Me%Size%ILB:Me%Size%IUB,   &
                             Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB ))

                    select case (ObjProperty%TypeZUV)

                        case (TypeZ_)
                            !eventual convert to faces is done afterwards
                            call HDF5SetLimits  (ObjHDF5File%HDFID, Me%WorkSize%ILB,    &
                                                 Me%WorkSize%IUB, Me%WorkSize%JLB,      &
                                                 Me%WorkSize%JUB,                       &
                                                 Me%WorkSize%KLB,Me%WorkSize%KUB,       &
                                                 STAT = STAT_CALL) 
                            if (STAT_CALL .NE. SUCCESS_)                                & 
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR05'

                            if (.not. ObjProperty%ConvertToFaces) then
                                !read field
                                call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group, &
                                            trim(ObjProperty%FirstField%Name),          &
                                            Array3D = ObjProperty%FirstField%Values3D,  &
                                            STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR06'

                            else !ConvertToFaces
                                !read field
                                call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group, &
                                                  trim(ObjProperty%FirstField%Name),    &
                                                  Array3D      = AuxField%Values3D,     &
                                                  STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR07'

                                call ConvertFieldToFaces(AuxField,ObjProperty%FirstField, &
                                                         ObjProperty%ID%IDNumber,       &
                                                         ObjProperty%CyclicBoundary)
                            endif

                        case (TypeU_)
                            call HDF5SetLimits  (ObjHDF5File%HDFID, Me%WorkSize%ILB,    &
                                                 Me%WorkSize%IUB, Me%WorkSize%JLB,      &
                                                 Me%WorkSize%JUB + 1,                   &
                                                 Me%WorkSize%KLB,Me%WorkSize%KUB,       &
                                                 STAT = STAT_CALL) 
                            if (STAT_CALL .NE. SUCCESS_)                                & 
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR08'

                            !read field
                            call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group,     &
                                              trim(ObjProperty%FirstField%Name),        &
                                              Array3D = ObjProperty%FirstField%Values3D, &
                                              STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR09'

                        case (TypeV_)
                            call HDF5SetLimits  (ObjHDF5File%HDFID, Me%WorkSize%ILB,    &
                                                 Me%WorkSize%IUB + 1, Me%WorkSize%JLB,  &
                                                 Me%WorkSize%JUB,                       &
                                                 Me%WorkSize%KLB,Me%WorkSize%KUB,       &
                                                 STAT = STAT_CALL) 
                            if (STAT_CALL .NE. SUCCESS_)                                & 
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR10'

                            !read field
                            call HDF5ReadData(ObjHDF5File%HDFID, ObjProperty%Group,     &
                                              trim(ObjProperty%FirstField%Name),        &
                                              Array3D = ObjProperty%FirstField%Values3D, &
                                              STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                stop 'OpenAndReadHDF5OneInstant - ModuleAssimilationPreProcessor - ERR11'

                    end select
 
            end select

            call killhdf5(ObjHDF5File%HDFID)

            if (ObjProperty%ConvertToFaces) then
                !Deallocates AuxField
                deallocate(AuxField%Values3D)
                deallocate(AuxField)
            endif

            ObjProperty => ObjProperty%Next

        end do

    end subroutine OpenAndReadHDF5OneInstant

    !--------------------------------------------------------------------------

    subroutine EOFAnalysis

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: idum, ieig, k, iter
        real(8)                                         :: tmp, angle
        real(8), dimension(:),      pointer             :: AuxValue
        integer                                         :: STAT_CALL

        integer, parameter                              :: miter = 10000
        real, parameter                                 :: prec = 1e-8
        
        !----------------------------------------------------------------------

        !EOFAnalysis performs the EOF analysis of Me%Covariance
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output,
        !
        ! Local parameters:
        !
        !

        idum = 1
        tmp  = 0.

        write(*,*)
        write(*,*) 'Performing EOF analysis...' 

        do ieig = 1, Me%StateCovRank
            if (Me%FullEOFAnalysis) then
                do k = 1, Me%StateVarNumber !nstatevars
                    Me%LMatrix(k,ieig) = NormalRand(idum)
                end do
            else
                do k = 1, Me%StatesNumber !nstates
                    Me%StatesLMatrix(k,ieig) = NormalRand(idum)
                end do
            endif

            call IEigenPower(ieig, miter - 1, ieig*prec, angle, iter)

            tmp = tmp + Me%UVector(ieig)
       
            write(*,*)
            write(*,*) 'EOF ', ieig, ' :'
            write(*,*) 'Number of iterations for estimation: ', iter
            write(*,*) 'EigenValue: ', Me%UVector(ieig)
            write(*,*) 'Fraction of variance explained (inertia): ',                &
                        Me%UVector(ieig)/Me%CovTrace
            write(*,*) 'Cumulative inertia: ', tmp/Me%CovTrace
            write(*,*) 'Precision: ', angle

            !Write EOF inertia to output HDF5
            allocate(AuxValue(1:1))
            AuxValue(1) = Me%UVector(ieig)/Me%CovTrace

            call HDF5SetLimits(Me%ObjCovHDF5, 1, 1, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            & 
            stop 'EOFAnalysis - ModuleAssimilationPreProcessor - ERR01'        

            call HDF5WriteData(Me%ObjCovHDF5, "/AssimilationData/Inertia",          &
                               "inertia", "-", Array1D = AuxValue,                  &
                               OutputNumber = ieig, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
            stop 'EOFAnalysis - ModuleAssimilationPreProcessor - ERR02'  

            deallocate(AuxValue)

            if (.not. Me%FullEOFAnalysis) then

                !calculate expansion coefficient for EOF
                call ExpansionCoefCalculation(ieig)

                !translate EOF in states space to state variables space and output
                call TranslateEOF(ieig)

                call WriteEOFToOutput(ieig)
            endif

        end do

        if (Me%FullEOFAnalysis) then

            call WriteOutput !This needs to be completed with eigenvalues
        endif

        if (.not. Me%StateReconstruction) call killhdf5(Me%ObjCovHDF5)

    end subroutine EOFAnalysis

    !--------------------------------------------------------------------------

    real function UniformRand01(idum)

        !Arguments-------------------------------------------------------------
        integer                                         :: idum

        !Local-----------------------------------------------------------------
        integer                                         :: k
        integer, parameter                              :: IA = 16807
        integer, parameter                              :: IM = 2147483647
        integer, parameter                              :: IQ=127773
        integer, parameter                              :: IR=2836
        real, parameter                                 :: AM=1./IM

        !----------------------------------------------------------------------

        !UniformRand01 returns a random number uniformely distributed in interval [0, 1]
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output, real UniformRand01, random real in the range 0 to 1.
        !
        ! Local parameters:
        !
        !

        k = idum/IQ
        idum = IA*(idum - k*IQ) - IR*k
        if (idum .lt. 0) then 
            idum = idum + IM
        endif

        UniformRand01 = AM*idum

    end function UniformRand01

    !--------------------------------------------------------------------------

    real function NormalRand(idum)

        !Arguments-------------------------------------------------------------
        integer                                         :: idum

        !Local-----------------------------------------------------------------
        integer , save                                  :: iset
        real    , save                                  :: gset
        real                                            :: v1, v2, fac, rsq

        !----------------------------------------------------------------------

        !NormalRand returns a random number following normal distribution N(0,1)
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output, real NormalRand, random real with distribution N(0,1)
        !
        ! Local parameters:
        !
        !

        DATA     iset /0/ !initializes iset


        if (iset .eq. 0) then
!10          v1 = 2.0 * UniformRand01(idum) - 1.0
            !following the same as(Fortran 77): if ((rsq .ge. 1.0).or.(rsq .eq. 0)) goto 10
            do
                v1 = 2.0 * UniformRand01(idum) - 1.0
                v2 = 2.0 * UniformRand01(idum) - 1.0
            
                rsq = v1*v1 + v2*v2
                if ((rsq .lt. 1.0) .and. (rsq /= 0)) exit
            end do
            
            fac = sqrt(- 2.0*alog(rsq)/rsq) !alog = neperian logarithm
            iset = 1
            gset = v2*fac
            NormalRand = v1*fac
        else
            iset = 0
            NormalRand = gset
        endif

    end function NormalRand

    !--------------------------------------------------------------------------

    subroutine IEigenPower(i, itermax, eps, angle, iter)

        !Arguments-------------------------------------------------------------
        integer,    intent(IN)                          :: i, itermax
        real,       intent(IN)                          :: eps
        real(8)                                         :: angle
        integer,    intent(OUT)                         :: iter

        !Local-----------------------------------------------------------------
        real(8), dimension(:), pointer                  :: vec2, vecr
        real(8), dimension(:,:), pointer                :: AuxLMatrix
        real(8), dimension(:), pointer                  :: AuxVector
        integer                                         :: j, k
        real(8)                                         :: tmp, a1, a2, b2, b, t, t2
        integer                                         :: n !nstates or nstatevars

        !----------------------------------------------------------------------

        !IEigenPower returns i-th eigenvector and eigenvalue of a covariance matrix
        !(calculated using the Power Method and iterations)
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input: i = eigen order, itermax = maximum number of
        !           iterations, eps = stop criterium, angle
        !
        !
        !    Input/output,
        !
        !    Output, eigvector = i-th eigen vector (EOF), eigvalue = i-th eigen value
        !
        ! Local parameters:
        !
        !

        !Allocate aux vectors:
        if (Me%FullEOFAnalysis) then
            n = Me%StateVarNumber
            AuxLMatrix => Me%LMatrix
        else
            n = Me%StatesNumber
            AuxLMatrix => Me%StatesLMatrix
        endif

        allocate(vecr(1:n))
        allocate(vec2(1:n))

        !1. Make new i-th eigen vector ortogonal with the others (which are normalized)
        do j = 1,i - 1
            tmp = 0
            do k = 1, n
                tmp = tmp + AuxLMatrix(k,i)*AuxLMatrix(k,j)
            end do
            do k = 1, n
                AuxLMatrix(k,i) = AuxLMatrix(k,i) - tmp*AuxLMatrix(k,j)
            end do
        end do

        b2 = 0.

        do k = 1, n
            b2 = b2 + AuxLMatrix(k,i)**2
        end do

        b = dsqrt(b2)
  
        do k = 1, n
            AuxLMatrix(k,i) = AuxLMatrix(k,i)/b
        end do

        if (Me%FullCovariance) then
            vecr(:) = MATMUL(Me%Covariance, AuxLMatrix(:,i))

        else
            allocate(AuxVector(1:n))
            AuxVector(:) = AuxLMatrix(:,i)
            
            call ProductMatVec(Me%Covariance, AuxVector, vecr, n)
            deallocate(AuxVector)
            nullify(AuxVector)
        endif

        a1 = 0.
        do k = 1, n
            a1 = a1 + vecr(k)*AuxLMatrix(k,i)
        end do

        do k = 1, n
            vecr(k) = vecr(k) - a1*AuxLMatrix(k,i)
        end do

        iter = 1

        !2. Make vecr ortogonal with the other eigen vectors (which are normalized)
 do1 :  do 
            do j = 1, i - 1 !500
                tmp = 0
                do k = 1, n
                    tmp = tmp + vecr(k)*AuxLMatrix(k,j)

                end do
                do k = 1, n
                    vecr(k) = vecr(k) - tmp*AuxLMatrix(k,j)

                end do
            end do

            b2 = 0.
        
            do k = 1, n
                b2 = b2 + vecr(k)**2
            end do

            b = dsqrt(b2)

            !stoping test
            angle = b/a1
        
            if (angle .lt. eps) exit do1 !goto 600

            iter = iter + 1

            if (iter .gt. itermax) exit do1 !goto 600

            do k = 1, n
                vec2(k) = vecr(k)/b
            end do

            if (Me%FullCovariance) then
                vecr = MATMUL(Me%Covariance, vec2)
            else
                call ProductMatVec(Me%Covariance, vec2, vecr, n)
            endif

            a2 = 0
            do k = 1, n
                a2 = a2 + vec2(k)*vecr(k)
            end do

            t = 2*b/(dabs(a1 - a2) + dsqrt((a1 - a2)**2 + 4*b2))

            if ((a1.le.a2) .and. ((a1.ne. a2).or.(a1+a2.lt.0))) then
                t = -t
            endif

            t2 = 1 + t**2

            if (a1**2 .gt. a2**2) then
                a1 = (a1 + (t**2)*a2 + 2*t*b)/t2 
                t2 = dsqrt(t2)
                tmp = t/t2
            
                do k = 1, n
                   vecr(k) = (vecr(k) - a2*vec2(k) - b*AuxLMatrix(k,i))*tmp
                   AuxLMatrix(k,i) = (AuxLMatrix(k,i) + t*vec2(k))/t2
                end do
            else
                a1 = (a2 + (t**2)*a1 - 2*t*b)/t2
                t2 = dsqrt(t2)
            
                do k = 1, n
                   vecr(k) = (vecr(k) - a2*vec2(k) - b*AuxLMatrix(k,i))/t2
                   AuxLMatrix(k,i) = (vec2(k) - t*AuxLMatrix(k,i))/t2
                end do
            endif
        !goto 500
        end do do1

        tmp = dsqrt(a1*a1 + b2) !600
      
        do k = 1, n
            AuxLMatrix(k,i) = (a1*AuxLMatrix(k,i) + vecr(k))/tmp
        end do

        Me%UVector(i) = a1

        !Deallocate aux vectors      
        deallocate(vecr)
        deallocate(vec2)
        nullify(vecr)
        nullify(vec2)

    end subroutine IEigenPower

    !--------------------------------------------------------------------------

    subroutine ProductMatVec(matrix, vector, product, n)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer                   :: matrix
        real(8), dimension(:), pointer                  :: vector, product
        integer                                         :: n

        !Local-----------------------------------------------------------------
        integer                                         :: ii, i, ij, j
        real(8)                                         :: s

        !----------------------------------------------------------------------

        !ProductMatVec returns the product of a symmetric matrix of rank n, arranged
        !in a 1D vector of size n*(n+1)/2.
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input, matrix = vector of size n*(n+1)/2, representing a symmetric matrix of
        !           rank n, vector = vector of size n.
        !
        !    Input/output,
        !
        !    Output, product = product of matrix and vector.
        !
        ! Local parameters:
        !
        !

        ii = 1

        do i = 1,n
            s = 0
            ij = ii
            do j = 1,i-1
                s = s + matrix(ij,1)*vector(j)
                ij = ij + 1
            enddo
            do j = i, n
                s = s + matrix(ij,1)*vector(j)
                ij = ij + j
            enddo
            product(i) = s
            ii = ii + i
        enddo

    end subroutine ProductMatVec

    !--------------------------------------------------------------------------

    subroutine TranslateEOF(ieig)

        !Arguments-------------------------------------------------------------
        integer                                     :: ieig

        !Local-----------------------------------------------------------------
        integer                                     :: state
        integer                                     :: i, j, k
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        real(8)                                     :: s
        type(T_HDF5File)        , pointer           :: ObjHDF5HydroFile, ObjHDF5WaterFile
        type(T_HDF5File)        , pointer           :: ObjHDF5File
        integer                                     :: CurrentInstant
        type(T_Property)        , pointer           :: ObjProperty

        !----------------------------------------------------------------------

        !TranslateEOFs translates the EOFs from states space to states variables space
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output,
        !
        ! Local parameters:
        !
        !

        Me%LMatrix(:,1) = 0.

        ObjHDF5HydroFile => Me%FirstHydroCovHDF5File
        ObjHDF5WaterFile => Me%FirstWaterCovHDF5File

        if (Me%NumberHydroCovHDF5 > 0) then
            ObjHDF5File => ObjHDF5HydroFile
        else
            ObjHDF5File => ObjHDF5WaterFile
        endif

        CurrentInstant = ObjHDF5File%StartInstant

do2:    do state = 1, Me%StatesNumber

            call OpenAndReadHDF5OneInstant(CurrentInstant, ObjHDF5HydroFile,            &
                                           ObjHDF5WaterFile)

            s = (1./sqrt(real(Me%StatesNumber)))*Me%StatesLMatrix(state,ieig) !*          &
                !(1./sqrt(Me%UVector(ieig)))

            Me%CurrentStateVar = 0

            !For each Property
            ObjProperty => Me%FirstProperty

            do while(associated(ObjProperty))

                !Cycle the property's variables
                select case (ObjProperty%Rank)

                    case(2)
                        !Get limits of property state window
                        ILB = ObjProperty%Window2D%ILB
                        IUB = ObjProperty%Window2D%IUB
                        JLB = ObjProperty%Window2D%JLB
                        JUB = ObjProperty%Window2D%JUB

                        do j = JLB, JUB
                        do i = ILB, IUB
                            if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then

                                Me%CurrentStateVar = Me%CurrentStateVar + 1
                                if (ObjProperty%FirstField%Values2D(i, j)               &
                                    >= FillValueReal / 2.) then

                                    Me%LMatrix(Me%CurrentStateVar,1) =                  &
                                        + Me%LMatrix(Me%CurrentStateVar,1)              &
                                        + s*((ObjProperty%FirstField%Values2D(i, j)     &
                                        - Me%AverageState(Me%CurrentStateVar))) !         &
                                        !/ObjProperty%AverStandardDev

                                endif
                            endif
                        enddo
                        enddo

                        deallocate(ObjProperty%FirstField%Values2D)
                    case(3)
                        !Get limits of property state window
                        ILB = ObjProperty%Window%ILB
                        IUB = ObjProperty%Window%IUB
                        JLB = ObjProperty%Window%JLB
                        JUB = ObjProperty%Window%JUB
                        KLB = ObjProperty%Window%KLB
                        KUB = ObjProperty%Window%KUB

                        if (ObjProperty%TypeZUV == TypeZ_ .and.                         &
                            .not. ObjProperty%ConvertToFaces) then

                            do k = KLB, KUB
                            do j = JLB, JUB
                            do i = ILB, IUB
                                if (Me%Mapping%IntegerValues3D(i, j, k) ==              &
                                    WaterPoint) then

                                    Me%CurrentStateVar = Me%CurrentStateVar + 1
                                    
                                    if (ObjProperty%FirstField%Values3D(i, j, k)        &
                                        >= FillValueReal / 2.) then

                                        Me%LMatrix(Me%CurrentStateVar,1) =              &
                                            Me%LMatrix(Me%CurrentStateVar,1)            &
                                            + s*                                        &
                                            (ObjProperty%FirstField%Values3D(i,j,k)     &
                                            - Me%AverageState(Me%CurrentStateVar)) !      &
                                            !/ObjProperty%AverStandardDev

                                    endif
                                endif
                            enddo
                            enddo
                            enddo

                        elseif (ObjProperty%TypeZUV == TypeU_ .or.                      &
                                (ObjProperty%ConvertToFaces .and.                       &
                                ObjProperty%ID%IDNumber == VelocityU_)) then

                            do k = KLB, KUB
                            do j = JLB, JUB + 1
                            do i = ILB, IUB
                                if (Me%MappingFacesU%IntegerValues3D(i, j, k) ==        &
                                    WaterPoint) then

                                    Me%CurrentStateVar = Me%CurrentStateVar + 1
                                    
                                    if (ObjProperty%FirstField%Values3D(i, j, k)        &
                                        >= FillValueReal / 2.) then

                                        Me%LMatrix(Me%CurrentStateVar,1) =              &
                                            Me%LMatrix(Me%CurrentStateVar,1)            &
                                            + s*                                        &
                                            (ObjProperty%FirstField%Values3D(i,j,k)     &
                                            - Me%AverageState(Me%CurrentStateVar)) !      &
                                            !/ObjProperty%AverStandardDev

                                    endif
                                endif
                            enddo
                            enddo
                            enddo

                        else !TypeV_ .or. (ConvertToFaces .and. VelocityV_)

                            do k = KLB, KUB
                            do j = JLB, JUB
                            do i = ILB, IUB + 1
                                if (Me%MappingFacesV%IntegerValues3D(i, j, k) ==        &
                                    WaterPoint) then

                                    Me%CurrentStateVar = Me%CurrentStateVar + 1
                                    
                                    if (ObjProperty%FirstField%Values3D(i, j, k)        &
                                        >= FillValueReal / 2.) then

                                        Me%LMatrix(Me%CurrentStateVar,1) =              &
                                            Me%LMatrix(Me%CurrentStateVar,1)            &
                                            + s*                                        &
                                            (ObjProperty%FirstField%Values3D(i,j,k)     &
                                            - Me%AverageState(Me%CurrentStateVar)) !      &
                                            !/ObjProperty%AverStandardDev

                                    endif
                                endif
                            enddo
                            enddo
                            enddo

                        endif
                        
                        deallocate(ObjProperty%FirstField%Values3D)
                end select

                ObjProperty => ObjProperty%Next               
            end do

            CurrentInstant = CurrentInstant + (Me%DecimationFactor + 1)
            if (CurrentInstant > ObjHDF5File%EndInstant) then

                if (Me%NumberHydroCovHDF5 > 0) then
                    ObjHDF5HydroFile => ObjHDF5HydroFile%Next
                    ObjHDF5File      => ObjHDF5HydroFile
                endif

                if (Me%NumberWaterCovHDF5 > 0) then
                    ObjHDF5WaterFile => ObjHDF5WaterFile%Next
                    if (Me%NumberHydroCovHDF5 == 0) then
                        ObjHDF5File  => ObjHDF5WaterFile
                    endif
                endif

                if (.not. (associated(ObjHDF5File))) exit do2

                CurrentInstant = ObjHDF5File%StartInstant
            endif
        enddo do2

    end subroutine TranslateEOF

    !--------------------------------------------------------------------------

    subroutine WriteEOFToOutput(neof)

        !Arguments-------------------------------------------------------------
        integer                                         :: neof

        !Local-----------------------------------------------------------------
        real, dimension(:,:),       pointer             :: AuxField2D, AbsAuxField2D
        real, dimension(:,:,:),     pointer             :: AuxField, AuxFieldFaces
        real, dimension(:,:,:),     pointer             :: AbsAuxField, AbsAuxFieldFaces
        real(8), dimension(:),      pointer             :: AuxValue
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: KLB, KUB
        integer                                         :: i, j, k
        integer                                         :: CurrentStateVar
        type(T_Property),           pointer             :: ObjProperty
        integer                                         :: STAT_CALL
        real(8)                                         :: EOFnorm

        !----------------------------------------------------------------------

        !(If property defined in faces it involves obtaining EOF values 
        !for cell center)
        !(EOF is defined for all WaterPoints)

        !Write EOF as property fields
        CurrentStateVar = 0

        EOFnorm = 0.

        !Cycle each property
        ObjProperty => Me%FirstProperty

        do while(associated(ObjProperty))

            select case (ObjProperty%Rank)

                case(2)

                    !field for relative EOF
                    allocate(AuxField2D(Me%WorkSize%ILB:Me%WorkSize%IUB,                &
                                        Me%WorkSize%JLB:Me%WorkSize%JUB))
                    AuxField2D(:,:) = FillValueReal

                    !field for absolute EOF
                    allocate(AbsAuxField2D(Me%WorkSize%ILB:Me%WorkSize%IUB,             &
                                        Me%WorkSize%JLB:Me%WorkSize%JUB))
                    AbsAuxField2D(:,:) = FillValueReal

                    ILB = ObjProperty%Window2D%ILB
                    IUB = ObjProperty%Window2D%IUB
                    JLB = ObjProperty%Window2D%JLB
                    JUB = ObjProperty%Window2D%JUB

                    do j = JLB, JUB
                    do i = ILB, IUB

                        if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then
                            
                            CurrentStateVar = CurrentStateVar + 1

                            AuxField2D(i,j) = Me%LMatrix(CurrentStateVar, 1)
                            
                            EOFnorm = EOFnorm + Me%LMatrix(CurrentStateVar, 1)**2

                            AbsAuxField2D(i,j) = AuxField2D(i,j) !*                       &
                                                 !sqrt(Me%UVector(neof))*                &
                                                 !ObjProperty%AverStandardDev
                        endif

                    enddo
                    enddo

                    call HDF5SetLimits(Me%ObjCovHDF5, Me%WorkSize%ILB,                  &
                                       Me%WorkSize%IUB, Me%WorkSize%JLB,                &
                                       Me%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR01'

                    !relative EOF
                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                                       "/Statistics/"//trim(ObjProperty%ID%Name)//"/EOF", &
                                       "eof", "-", Array2D = AuxField2D,                &
                                       OutputNumber = neof, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR02'            

                    deallocate(AuxField2D)
                    nullify(AuxField2D)

                    !absolute EOF
                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                            "/Statistics/"//trim(ObjProperty%ID%Name)//"/AbsoluteEOF",  &
                                       "absolute eof", ObjProperty%ID%Units,            &
                                       Array2D = AbsAuxField2D, OutputNumber = neof,    &
                                       STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR03'            

                    deallocate(AbsAuxField2D)
                    nullify(AbsAuxField2D)

                case(3)

                    ILB = ObjProperty%Window%ILB
                    IUB = ObjProperty%Window%IUB
                    JLB = ObjProperty%Window%JLB
                    JUB = ObjProperty%Window%JUB
                    KLB = ObjProperty%Window%KLB
                    KUB = ObjProperty%Window%KUB

                    !allocate a generic field (cell centered)
                    !relative EOF
                    allocate(AuxField(Me%WorkSize%ILB:Me%WorkSize%IUB,                  &
                                      Me%WorkSize%JLB:Me%WorkSize%JUB,                  &
                                      Me%WorkSize%KLB:Me%WorkSize%KUB))
                    AuxField(:,:,:) = FillValueReal

                    !absolute EOF 
                    allocate(AbsAuxField(Me%WorkSize%ILB:Me%WorkSize%IUB,               &
                                         Me%WorkSize%JLB:Me%WorkSize%JUB,               &
                                         Me%WorkSize%KLB:Me%WorkSize%KUB))
                    AbsAuxField(:,:,:) = FillValueReal

                    if (ObjProperty%TypeZUV == TypeZ_ .and.                             &
                        .not. ObjProperty%ConvertToFaces) then

                        do k = KLB, KUB
                        do j = JLB, JUB
                        do i = ILB, IUB
                            if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then

                                CurrentStateVar = CurrentStateVar + 1

                                AuxField(i,j,k) = Me%LMatrix(CurrentStateVar, 1)

                                EOFnorm = EOFnorm + Me%LMatrix(CurrentStateVar, 1)**2

                                AbsAuxField(i,j,k) = AuxField(i,j,k) !*                   &
                                                    !sqrt(Me%UVector(neof))*             &
                                                    !ObjProperty%AverStandardDev
                            endif
                        enddo
                        enddo
                        enddo

                    elseif (ObjProperty%TypeZUV == TypeU_ .or.                          &
                            (ObjProperty%ConvertToFaces .and.                           &
                            ObjProperty%ID%IDNumber == VelocityU_)) then

                        !relative EOF
                        allocate(AuxFieldFaces(Me%Size%ILB:Me%Size%IUB,                 &
                                               Me%Size%JLB:Me%Size%JUB,                 &
                                               Me%Size%KLB:Me%Size%KUB))
                        AuxFieldFaces(:,:,:) = FillValueReal

                        !absolute EOF 
                        allocate(AbsAuxFieldFaces(Me%Size%ILB:Me%Size%IUB,              &
                                                  Me%Size%JLB:Me%Size%JUB,              &
                                                  Me%Size%KLB:Me%Size%KUB))
                        AbsAuxFieldFaces(:,:,:) = FillValueReal

                        do k = KLB, KUB
                        do j = JLB, JUB + 1
                        do i = ILB, IUB
                            if (Me%MappingFacesU%IntegerValues3D(i, j, k) == 1) then

                                CurrentStateVar = CurrentStateVar + 1

                                AuxFieldFaces(i,j,k) = Me%LMatrix(CurrentStateVar, 1)

                                EOFnorm = EOFnorm + Me%LMatrix(CurrentStateVar, 1)**2

                                AbsAuxFieldFaces(i,j,k) = AuxFieldFaces(i,j,k) !*         &
                                                          !sqrt(Me%UVector(neof))*       &
                                                          !ObjProperty%AverStandardDev
                                !all state variables have an eof value (non fillvalue)
                            endif
                        enddo
                        enddo
                        enddo

                        !relative EOF
                        call CenterField(AuxField, AuxFieldFaces,                       &
                                         ObjProperty%ID%IDNumber, Me%WorkSize)

                        deallocate(AuxFieldFaces)
                        nullify(AuxFieldFaces)

                        !absolute EOF
                        call CenterField(AbsAuxField, AbsAuxFieldFaces,                 &
                                         ObjProperty%ID%IDNumber, Me%WorkSize)

                        deallocate(AbsAuxFieldFaces)
                        nullify(AbsAuxFieldFaces)

                    else !TypeV_ .or. (ConvertToFaces .and. VelocityV_)

                        !relative EOF
                        allocate(AuxFieldFaces(Me%Size%ILB:Me%Size%IUB,                 &
                                               Me%Size%JLB:Me%Size%JUB,                 &
                                               Me%Size%KLB:Me%Size%KUB))
                        AuxFieldFaces(:,:,:) = FillValueReal

                        !absolute EOF
                        allocate(AbsAuxFieldFaces(Me%Size%ILB:Me%Size%IUB,              &
                                                  Me%Size%JLB:Me%Size%JUB,              &
                                                  Me%Size%KLB:Me%Size%KUB))
                        AbsAuxFieldFaces(:,:,:) = FillValueReal

                        do k = KLB, KUB
                        do j = JLB, JUB
                        do i = ILB, IUB + 1
                            if (Me%MappingFacesV%IntegerValues3D(i, j, k) == 1) then

                                CurrentStateVar = CurrentStateVar + 1

                                AuxFieldFaces(i,j,k) = Me%LMatrix(CurrentStateVar, 1)

                                EOFnorm = EOFnorm + Me%LMatrix(CurrentStateVar, 1)**2

                                AbsAuxFieldFaces(i,j,k) = AuxFieldFaces(i,j,k) !*         &
                                                          !sqrt(Me%UVector(neof))*       &
                                                          !ObjProperty%AverStandardDev
                                !all state variables have an eof value (non fillvalue)
                            endif
                        enddo
                        enddo
                        enddo

                        !relative EOF
                        call CenterField(AuxField, AuxFieldFaces,                       &
                                         ObjProperty%ID%IDNumber, Me%WorkSize)

                        deallocate(AuxFieldFaces)
                        nullify(AuxFieldFaces)

                        !absolute EOF
                        call CenterField(AbsAuxField, AbsAuxFieldFaces,                 &
                                         ObjProperty%ID%IDNumber, Me%WorkSize)

                        deallocate(AbsAuxFieldFaces)
                        nullify(AbsAuxFieldFaces)

                    endif

                    call HDF5SetLimits(Me%ObjCovHDF5, Me%WorkSize%ILB,                  &
                                       Me%WorkSize%IUB, Me%WorkSize%JLB,                &
                                       Me%WorkSize%JUB, Me%WorkSize%KLB,                &
                                       Me%WorkSize%KUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR04'

                    !relative EOF
                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                           "/Statistics/"//trim(ObjProperty%ID%Name)//"/EOF","eof","-", &
                                       Array3D = AuxField, OutputNumber = neof,         &
                                       STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR05'

                    deallocate(AuxField)
                    nullify(AuxField)

                    !absolute EOF
                    call HDF5WriteData(Me%ObjCovHDF5,                                   &
                            "/Statistics/"//trim(ObjProperty%ID%Name)//"/AbsoluteEOF",  &
                                       "absolute eof",ObjProperty%ID%Units,             &
                                       Array3D = AbsAuxField, OutputNumber = neof,      &
                                       STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR06'

                    deallocate(AbsAuxField)
                    nullify(AbsAuxField)

            end select

            ObjProperty => ObjProperty%Next
        end do

        !Write EOF as state fields
        call HDF5SetLimits(Me%ObjCovHDF5, 1, Me%StateVarNumber, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    & 
        stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR07'

        allocate(AuxValue(1:Me%StateVarNumber))
        AuxValue(:) = Me%LMatrix(:,1)

        call HDF5WriteData(Me%ObjCovHDF5,                                               &
                           "/AssimilationData/EOF", "eof", "-",                         &
                           Array1D = AuxValue, OutputNumber = neof, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
        stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR08'  

        deallocate(AuxValue)

        !Write corresponding EigenValue
        call HDF5SetLimits(Me%ObjCovHDF5, 1, 1, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    & 
        stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR09'

        allocate(AuxValue(1:1))
        AuxValue(1) = Me%UVector(neof)

        call HDF5WriteData(Me%ObjCovHDF5,                                               &
                           "/AssimilationData/EigenValue",                              & 
                           "eigenvalue",                                                & 
                           "-",                                                         &
                           Array1D = AuxValue, OutputNumber = neof, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
        stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR10'  

        deallocate(AuxValue)

        call HDF5FlushMemory(Me%ObjCovHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteEOFToOutput - ModuleAssimilationPreProcessor - ERR11'

        EOFnorm = sqrt(EOFnorm)
        write(*,*) 'EOF norm:', EOFnorm

    end subroutine WriteEOFToOutput

    !--------------------------------------------------------------------------

    subroutine CenterField(PropCenter, PropFaces, IDNumber, Size)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :),   pointer             :: PropCenter, PropFaces
        integer                                         :: IDNumber
        type (T_Size3D)                                 :: Size

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: KLB, KUB
        integer                                         :: i, j, k

        !----------------------------------------------------------------------

        !Shorten
        ILB = Size%ILB 
        IUB = Size%IUB 
        JLB = Size%JLB 
        JUB = Size%JUB 
        KLB = Size%KLB 
        KUB = Size%KUB 

        if (IDNumber == VelocityU_) then

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then

                    !Check for non state faces (assume MISSING_NULL)
                    if (PropFaces(i,j,k) < FillValueReal / 2.) PropFaces(i,j,k) = 0.
                    if (PropFaces(i,j+1,k) < FillValueReal / 2.) PropFaces(i,j+1,k) = 0.

                    PropCenter(i, j, k) = (PropFaces(i, j, k) +                         &
                                           PropFaces(i, j+1, k)) / 2.
                else

                    PropCenter(i, j, k) = 0.
                endif
            enddo
            enddo
            enddo

        elseif (IDNumber == VelocityV_) then

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then

                    !Check for non state faces (assume MISSING_NULL)
                    if (PropFaces(i,j,k) < FillValueReal / 2.) PropFaces(i,j,k) = 0.
                    if (PropFaces(i+1,j,k) < FillValueReal / 2.) PropFaces(i+1,j,k) = 0.

                    PropCenter(i, j, k) = (PropFaces(i, j, k) +                         &
                                           PropFaces(i+1, j, k)) / 2.
                else

                    PropCenter(i, j, k) = 0.
                endif
            enddo
            enddo
            enddo
        else

            stop 'CenterField - ModuleAssimilationPreProcessor - ERR01'  

        endif

    end subroutine CenterField

    !--------------------------------------------------------------------------

    subroutine WriteEigenValuesToOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer                     :: AuxValue
        integer                                         :: STAT_CALL

        !----------------------------------------------------------------------

        write(*,*)
        write(*,*) 'Writing eigenvalues to output HDF5...'

        !Write eigenvalue
        call HDF5SetLimits(Me%ObjCovHDF5, 1, 1, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                        & 
        stop 'WriteEigenValuesToOutput - ModuleAssimilationPreProcessor - ERR01'

        allocate(AuxValue(1:1))

        AuxValue(1) = Me%UVector(1)

        call HDF5WriteData(Me%ObjCovHDF5,                                   &
                           "/AssimilationData/EigenValue",                  & 
                           "eigenvalue", "-", Array1D = AuxValue,           &
                           OutputNumber = 1,                                &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'WriteEigenValuesToOutput - ModuleAssimilationPreProcessor - ERR02'  

        call HDF5FlushMemory(Me%ObjCovHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'WriteEigenValuesToOutput - ModuleAssimilationPreProcessor - ERR03'

    end subroutine WriteEigenValuesToOutput

    !--------------------------------------------------------------------------

    subroutine ExpansionCoefCalculation(ieig)

        !Arguments-------------------------------------------------------------
        integer                                     :: ieig

        !Local-----------------------------------------------------------------
        real(8), dimension(:), pointer              :: InstantExpansionCoef
        real(8), dimension(:), pointer              :: LVector
        type(T_CovarianceTime), pointer             :: ObjCovarianceTime
        integer                                     :: state

        !----------------------------------------------------------------------

        write(*,*)
        write(*,*) 'Calculating expansion coefficient...'

        allocate(InstantExpansionCoef(1:Me%StatesNumber))
        InstantExpansionCoef(:) = 0.

        allocate(LVector(1:Me%StatesNumber))
        LVector(:) = Me%StatesLMatrix(:,ieig)

        !calculate expansion coefficient
        if (.not. Me%FullCovariance) then
            call ProductMatVec(Me%Covariance, LVector, InstantExpansionCoef,            &
                               Me%StatesNumber)
            !InstantExpansionCoef = InstantExpansionCoef*(1./sqrt(Me%UVector(ieig)))
            InstantExpansionCoef = InstantExpansionCoef*(1./(Me%UVector(ieig)))
        else
            !InstantExpansionCoef = MATMUL(Me%Covariance, LVector)*                    &
            !                       (1./sqrt(Me%UVector(ieig)))
            InstantExpansionCoef = MATMUL(Me%Covariance, LVector)*                      &
                                   (1./(Me%UVector(ieig)))
        endif

        deallocate(LVector)

        !write expansion coefficient to time series file
        !(one time series file per each eof, so that memory is minimized)
        state = 0
        ObjCovarianceTime => Me%FirstCovarianceTime

        do while(associated(ObjCovarianceTime))
            Me%CurrentTime = ObjCovarianceTime%Time
            state = state + 1

            call WriteExpCoefTimeSerie(InstantExpansionCoef(state), Me%ExpansionCoef(ieig))

            ObjCovarianceTime => ObjCovarianceTime%Next
        enddo

        deallocate(InstantExpansionCoef)

        call KillExpCoefTimeSeries(Me%ExpansionCoef(ieig))

    end subroutine ExpansionCoefCalculation

    !--------------------------------------------------------------------------

    subroutine WriteExpCoefTimeSerie(InstantValue, ObjExpansionCoef)

        !Arguments-------------------------------------------------------------
        real(8)                                     :: InstantValue
        type (T_ExpansCoefTS)                       :: ObjExpansionCoef

        !Local-----------------------------------------------------------------
        integer                                     :: IBC
        
        !Begin-----------------------------------------------------------------

        !(this is made similarly to time serie subroutines)
        !(the buffer is TimeSerieData in each location)

        !Write time and data value to buffer
        !Increments the internal buffer count
        ObjExpansionCoef%BufferCount = ObjExpansionCoef%BufferCount + 1 

        !Shorten variable
        IBC = ObjExpansionCoef%BufferCount

        !Stores the current time
        ObjExpansionCoef%TimeBuffer(IBC) = Me%CurrentTime

        !Stores the data value
        ObjExpansionCoef%TimeSerieData(IBC) = InstantValue

        !If buffer is full then write to file
        if (ObjExpansionCoef%BufferCount  == Me%BufferSize) then
           call WriteBufferToFile(ObjExpansionCoef)
           ObjExpansionCoef%BufferCount = 0
        endif

    end subroutine WriteExpCoefTimeSerie

    !--------------------------------------------------------------------------

    subroutine WriteBufferToFile(ObjExpansionCoef)

        !Arguments-------------------------------------------------------------
        type (T_ExpansCoefTS)                       :: ObjExpansionCoef

        !Local-----------------------------------------------------------------
        integer                                     :: iB, unit
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        
        !Begin-----------------------------------------------------------------

        unit = ObjExpansionCoef%ObjTimeSerie

        do iB = 1, ObjExpansionCoef%BufferCount

            !Writes date in the form of 1999 08 08 23 59 23
            call ExtractDate(ObjExpansionCoef%TimeBuffer(iB), Year = Year,      &
                             Month = Month, Day = Day, Hour = Hour,             &
                             Minute = Minute, Second = Second)

            !Writes time in seconds since the beginning of the output
            !Writes value in the buffer in exp format
            write(unit, fmt=1001) ObjExpansionCoef%TimeBuffer(iB) -             &
                                  Me%FirstCovarianceTime%Time, int(Year),       &
                                  int(Month), int(Day), int(Hour), int(Minute), &
                                  Second, ObjExpansionCoef%TimeSerieData(iB)
        enddo

    1001 format(f13.2, 1x, i4, 2x, i2, 2x, i2, 2x, i2, 2x, i2, 2x, f7.4,        & 
                2x, e20.12e3)

    end subroutine WriteBufferToFile

    !--------------------------------------------------------------------------

    subroutine KillExpCoefTimeSeries(ObjExpansionCoef)

        !Arguments-------------------------------------------------------------
        type (T_ExpansCoefTS)                       :: ObjExpansionCoef

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Disposes buffer
        call WriteBufferToFile(ObjExpansionCoef)

        !Closes file
        write(ObjExpansionCoef%ObjTimeSerie, *) '<EndTimeSerie>'

        call UnitsManager(ObjExpansionCoef%ObjTimeSerie, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)
            write(*,*)'Error closing Time Serie data file', trim(ObjExpansionCoef%Name)
            write(*,*)'KillExpCoefTimeSeries - ModuleAssimilationPreProcessor - WRN01'
        endif

        deallocate(ObjExpansionCoef%TimeSerieData, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'KillExpCoefTimeSeries - ModuleAssimilationPreProcessor - ERR01'
        nullify(ObjExpansionCoef%TimeSerieData)

        deallocate(ObjExpansionCoef%TimeBuffer, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'KillExpCoefTimeSeries - ModuleAssimilationPreProcessor - ERR02'
        nullify(ObjExpansionCoef%TimeBuffer)

    end subroutine KillExpCoefTimeSeries

    !--------------------------------------------------------------------------

    subroutine StateReconstruction

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: state, ieig, statevar
        type(T_CovarianceTime), pointer             :: ObjCovarianceTime
        integer                                     :: STAT_CALL
        real(8), dimension(:), pointer              :: LVector
        integer, pointer, dimension(:)              :: ColumnsToRead
        type (T_Property), pointer                  :: ObjProperty

        !----------------------------------------------------------------------

        write(*,*)
        write(*,*)'Reconstructing states...'

        allocate(LVector(1:Me%StateVarNumber))
        allocate(Me%StateVector(1:Me%StateVarNumber))
        allocate(ColumnsToRead(1:1))
        ColumnsToRead(1) = 8

        state = 0
        ObjCovarianceTime => Me%FirstCovarianceTime

        do while(associated(ObjCovarianceTime))
            
            Me%CurrentTime = ObjCovarianceTime%Time
            state = state + 1
            Me%StateVector(:) = 0.

            do ieig = 1, Me%StateCovRank

                if (state == 1) then
                    Me%ExpansionCoef(ieig)%ObjTimeSerie = 0
   
                    !read expansion coef.
                    call StartTimeSerieInput(Me%ExpansionCoef(ieig)%ObjTimeSerie,       &
                                             Me%ExpansionCoef(ieig)%Name,               &
                                             ColumnsToRead = ColumnsToRead,             &
                                             STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'StateReconstruction - ModuleAssimilationPreProcessor - ERR01'

                    call GetTimeSerieDataMatrix(Me%ExpansionCoef(ieig)%ObjTimeSerie,    &
                                                Me%ExpansionCoef(ieig)%TSDataMatrix,    &
                                                STAT = STAT_CALL)
                endif

                !read eof
                LVector(:) = FillValueReal

                call HDF5SetLimits (Me%ObjCovHDF5, 1, Me%StateVarNumber, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'StateReconstruction - ModuleSequentialAssimilation - ERR02'
                    
                call HDF5ReadData(Me%ObjCovHDF5, "/AssimilationData/EOF",               &
                                  "eof",                                                &
                                  Array1D = LVector,                                    &
                                  OutputNumber = ieig,                                  &
                                  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'StateReconstruction - ModuleSequentialAssimilation - ERR03'

                ObjProperty => Me%FirstProperty

                !cycle variables
                do statevar = 1, Me%StateVarNumber
                    
                    if ((statevar - ObjProperty%FirstStatePosition + 1)                     &
                        > ObjProperty%StateVarNumber) then
                        ObjProperty => ObjProperty%Next
                    endif

                    !partial sum
                    Me%StateVector(statevar) = Me%StateVector(statevar) +               &
                        Me%ExpansionCoef(ieig)%TSDataMatrix(state,2)*LVector(statevar) !  &
                        !*ObjProperty%AverStandardDev
                enddo
            enddo

            !sum average
            Me%StateVector = Me%StateVector*(sqrt(real(Me%StatesNumber)))               &
                             + Me%AverageState
            !Me%StateVector = Me%StateVector               &
            !                 + Me%AverageState

            !write state to output file
            call WriteStateToOutput(state)

            ObjCovarianceTime => ObjCovarianceTime%Next
        enddo

        deallocate(LVector)
        deallocate(Me%StateVector)
        deallocate(ColumnsToRead)

        call killhdf5(Me%ObjCovHDF5)
        call killhdf5(Me%ObjStateRecHDF5)

        !kill expansion coef. DataMatrix and time series
        do ieig = 1, Me%StateCovRank
            nullify(Me%ExpansionCoef(ieig)%TSDataMatrix)
            call KillTimeSerie(Me%ExpansionCoef(ieig)%ObjTimeSerie)
        enddo

    end subroutine StateReconstruction

    !--------------------------------------------------------------------------

    subroutine WriteStateToOutput(state)

        !Arguments-------------------------------------------------------------
        integer                                         :: state

        !Local-----------------------------------------------------------------
        real, dimension(:,:),       pointer             :: AuxField2D
        real, dimension(:,:,:),     pointer             :: AuxField, AuxFieldFaces
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: KLB, KUB
        integer                                         :: i, j, k
        integer                                         :: CurrentStateVar
        type(T_Property),           pointer             :: ObjProperty
        integer                                         :: STAT_CALL

        !----------------------------------------------------------------------

        !(If property defined in faces it involves obtaining state values 
        !for cell center)
        !(State is defined for all WaterPoints)

        !Write state as property fields
        CurrentStateVar = 0

        !Cycle each property
        ObjProperty => Me%FirstProperty

        do while(associated(ObjProperty))

            select case (ObjProperty%Rank)

                case(2)

                    allocate(AuxField2D(Me%WorkSize%ILB:Me%WorkSize%IUB,                &
                                        Me%WorkSize%JLB:Me%WorkSize%JUB))
                    AuxField2D(:,:) = FillValueReal

                    ILB = ObjProperty%Window2D%ILB
                    IUB = ObjProperty%Window2D%IUB
                    JLB = ObjProperty%Window2D%JLB
                    JUB = ObjProperty%Window2D%JUB

                    do j = JLB, JUB
                    do i = ILB, IUB

                        if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then
                            
                            CurrentStateVar = CurrentStateVar + 1

                            AuxField2D(i,j) = Me%StateVector(CurrentStateVar)

                        endif

                    enddo
                    enddo

                    call HDF5SetLimits(Me%ObjStateRecHDF5, Me%WorkSize%ILB,             &
                                       Me%WorkSize%IUB, Me%WorkSize%JLB,                &
                                       Me%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'WriteStateToOutput - ModuleAssimilationPreProcessor - ERR01'

                    call HDF5WriteData(Me%ObjStateRecHDF5,                              &
                                       "/Results/"//trim(ObjProperty%ID%Name),          &
                                       trim(ObjProperty%ID%Name), "-",                  &
                                       Array2D = AuxField2D,                            &
                                       OutputNumber = state, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'WriteStateToOutput - ModuleAssimilationPreProcessor - ERR02'            

                    deallocate(AuxField2D)
                    nullify(AuxField2D)

                case(3)

                    ILB = ObjProperty%Window%ILB
                    IUB = ObjProperty%Window%IUB
                    JLB = ObjProperty%Window%JLB
                    JUB = ObjProperty%Window%JUB
                    KLB = ObjProperty%Window%KLB
                    KUB = ObjProperty%Window%KUB

                    !allocate a generic field (cell centered)
                    allocate(AuxField(Me%WorkSize%ILB:Me%WorkSize%IUB,                  &
                                      Me%WorkSize%JLB:Me%WorkSize%JUB,                  &
                                      Me%WorkSize%KLB:Me%WorkSize%KUB))
                    AuxField(:,:,:) = FillValueReal

                    if (ObjProperty%TypeZUV == TypeZ_ .and.                             &
                        .not. ObjProperty%ConvertToFaces) then

                        do k = KLB, KUB
                        do j = JLB, JUB
                        do i = ILB, IUB
                            if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then

                                CurrentStateVar = CurrentStateVar + 1

                                AuxField(i,j,k) = Me%StateVector(CurrentStateVar)
                            endif
                        enddo
                        enddo
                        enddo

                    elseif (ObjProperty%TypeZUV == TypeU_ .or.                          &
                            (ObjProperty%ConvertToFaces .and.                           &
                            ObjProperty%ID%IDNumber == VelocityU_)) then

                        allocate(AuxFieldFaces(Me%WorkSize%ILB:Me%WorkSize%IUB,         &
                                               Me%WorkSize%JLB:Me%WorkSize%JUB + 1,     &
                                               Me%WorkSize%KLB:Me%WorkSize%KUB))
                        AuxFieldFaces(:,:,:) = FillValueReal

                        do k = KLB, KUB
                        do j = JLB, JUB + 1
                        do i = ILB, IUB
                            if (Me%MappingFacesU%IntegerValues3D(i, j, k) == 1) then

                                CurrentStateVar = CurrentStateVar + 1

                                AuxFieldFaces(i,j,k) =                                  &
                                    Me%StateVector(CurrentStateVar)
                                !all state variables have an non fillvalue
                            endif
                        enddo
                        enddo
                        enddo

                        call CenterField(AuxField, AuxFieldFaces,                       &
                                         ObjProperty%ID%IDNumber, Me%WorkSize)

                        deallocate(AuxFieldFaces)
                        nullify(AuxFieldFaces)

                    else !TypeV_ or (ConvertToFaces and VelocityV_)

                        allocate(AuxFieldFaces(Me%WorkSize%ILB:Me%WorkSize%IUB + 1,     &
                                               Me%WorkSize%JLB:Me%WorkSize%JUB,         &
                                               Me%WorkSize%KLB:Me%WorkSize%KUB))
                        AuxFieldFaces(:,:,:) = FillValueReal

                        do k = KLB, KUB
                        do j = JLB, JUB
                        do i = ILB, IUB + 1
                            if (Me%MappingFacesV%IntegerValues3D(i, j, k) == 1) then

                                CurrentStateVar = CurrentStateVar + 1

                                AuxFieldFaces(i,j,k) =                                  &
                                    Me%StateVector(CurrentStateVar)
                                !all state variables have an non fillvalue
                            endif
                        enddo
                        enddo
                        enddo

                        call CenterField(AuxField, AuxFieldFaces,                       &
                                         ObjProperty%ID%IDNumber, Me%WorkSize)

                        deallocate(AuxFieldFaces)
                        nullify(AuxFieldFaces)                         
                    endif

                    call HDF5SetLimits(Me%ObjStateRecHDF5, Me%WorkSize%ILB,             &
                                       Me%WorkSize%IUB, Me%WorkSize%JLB,                &
                                       Me%WorkSize%JUB, Me%WorkSize%KLB,                &
                                       Me%WorkSize%KUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        & 
                    stop 'WriteStateToOutput - ModuleAssimilationPreProcessor - ERR03'

                    call HDF5WriteData(Me%ObjStateRecHDF5,                              &
                                       "/Results/"//trim(ObjProperty%ID%Name),          &
                                       trim(ObjProperty%ID%Name),"-",                   &
                                       Array3D = AuxField, OutputNumber = state,        &
                                       STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'WriteStateToOutput - ModuleAssimilationPreProcessor - ERR04'

                    deallocate(AuxField)
                    nullify(AuxField)
            end select

            ObjProperty => ObjProperty%Next
        end do

        call HDF5FlushMemory(Me%ObjStateRecHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteStateToOutput - ModuleAssimilationPreProcessor - ERR09'

    end subroutine WriteStateToOutput

    !--------------------------------------------------------------------------

    subroutine WriteOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:,:),   pointer                 :: AuxField2D
        real, dimension(:,:,:), pointer                 :: AuxField
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: KLB, KUB
        integer                                         :: i, j, k, neof
        integer                                         :: CurrentStateVar
        type(T_Property), pointer                       :: ObjProperty
        integer                                         :: STAT_CALL

        !----------------------------------------------------------------------

        write(*,*)
        write(*,*) 'Writing EOF property fields to output HDF5...'

        !Cycle each EOF
        do neof = 1, Me%StateCovRank

            CurrentStateVar = 0

            !Cycle each property
            ObjProperty => Me%FirstProperty

            do while(associated(ObjProperty))

                select case (ObjProperty%Rank)

                    case(2)

                        allocate(AuxField2D(Me%WorkSize%ILB:Me%WorkSize%IUB,    &
                                 Me%WorkSize%JLB:Me%WorkSize%JUB))

                        ILB = ObjProperty%Window2D%ILB
                        IUB = ObjProperty%Window2D%IUB
                        JLB = ObjProperty%Window2D%JLB
                        JUB = ObjProperty%Window2D%JUB

                        do j = JLB, JUB
                        do i = ILB, IUB

                            if (Me%Mapping%IntegerValues2D(i, j) == WaterPoint) then
                                
                                CurrentStateVar = CurrentStateVar + 1

                                AuxField2D(i,j) = Me%LMatrix(CurrentStateVar, neof)

                            endif
                        enddo
                        enddo

                        call HDF5SetLimits(Me%ObjCovHDF5, Me%WorkSize%ILB,      &
                                            Me%WorkSize%IUB, Me%WorkSize%JLB,   &
                                            Me%WorkSize%JUB, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                            & 
                        stop 'WriteOutput - ModuleAssimilationPreProcessor - ERR01'

                        call HDF5WriteData(Me%ObjCovHDF5,                       &
                                           "/EOF/"//trim(ObjProperty%ID%Name),  & 
                                           ObjProperty%ID%Name,                 & 
                                           ObjProperty%ID%Units,                &
                                           Array2D = AuxField2D,                &
                                           OutputNumber = neof,                 &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'WriteOutput - ModuleAssimilationPreProcessor - ERR02'            

                        deallocate(AuxField2D)
                        nullify(AuxField2D)

                    case(3)

                        allocate(AuxField(Me%WorkSize%ILB:Me%WorkSize%IUB,      &
                                 Me%WorkSize%JLB:Me%WorkSize%JUB,               &
                                 Me%WorkSize%KLB:Me%WorkSize%KUB))

                        ILB = ObjProperty%Window%ILB
                        IUB = ObjProperty%Window%IUB
                        JLB = ObjProperty%Window%JLB
                        JUB = ObjProperty%Window%JUB
                        KLB = ObjProperty%Window%KLB
                        KUB = ObjProperty%Window%KUB

                        do k = KLB, KUB
                        do j = JLB, JUB
                        do i = ILB, IUB

                            if (Me%Mapping%IntegerValues3D(i, j, k) == WaterPoint) then

                                CurrentStateVar = CurrentStateVar + 1

                                AuxField(i,j,k) = Me%LMatrix(CurrentStateVar, neof)

                            endif
                        enddo
                        enddo
                        enddo

                        call HDF5SetLimits(Me%ObjCovHDF5, Me%WorkSize%ILB,      &
                                           Me%WorkSize%IUB, Me%WorkSize%JLB,    &
                                           Me%WorkSize%JUB, Me%WorkSize%KLB,    &
                                           Me%WorkSize%KUB, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                            & 
                        stop 'WriteOutput - ModuleAssimilationPreProcessor - ERR03'

                        call HDF5WriteData(Me%ObjCovHDF5,                       &
                                           "/EOF/"//trim(ObjProperty%ID%Name),  & 
                                           ObjProperty%ID%Name,                 & 
                                           ObjProperty%ID%Units,                &
                                           Array3D = AuxField,                  &
                                           OutputNumber = neof,                 &
                                           STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                              &
                            stop 'WriteOutput - ModuleAssimilationPreProcessor - ERR04'

                        deallocate(AuxField)
                        nullify(AuxField)

                end select

                ObjProperty => ObjProperty%Next
            end do

            call HDF5FlushMemory(Me%ObjCovHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'WriteOutput - ModuleAssimilationPreProcessor - ERR05'
        end do

    end subroutine WriteOutput

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine KillAssimPreProcessor

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5FileX, HDF5FileToKill
        type(T_Property), pointer                   :: PropertyX, PropertyToKill
        integer                                     :: STAT_CALL

        !------------------------------------------------------------------------

        !Kill the mapping field
        if (associated(Me%Mapping%IntegerValues3D)) then
            deallocate(Me%Mapping%IntegerValues3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR01'  
            nullify(Me%Mapping%IntegerValues3D)
        endif

        if (associated(Me%Mapping%IntegerValues2D)) then
            deallocate(Me%Mapping%IntegerValues2D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR02'  
            nullify(Me%Mapping%IntegerValues2D)
        endif

        if (associated(Me%MappingFacesU%IntegerValues3D)) then
            deallocate(Me%MappingFacesU%IntegerValues3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR03'  
            nullify(Me%MappingFacesU%IntegerValues3D)
        endif

        if (associated(Me%MappingFacesV%IntegerValues3D)) then
            deallocate(Me%MappingFacesV%IntegerValues3D, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR04'  
            nullify(Me%MappingFacesV%IntegerValues3D)
        endif

        !Kill assimilation variables
        deallocate(Me%Covariance, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR05'
        nullify(Me%Covariance)

        deallocate(Me%AverageState, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR06'
        nullify(Me%AverageState)

        deallocate(Me%UVector, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR07'
        nullify(Me%UVector)

        deallocate(Me%LMatrix, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR08'
        nullify(Me%LMatrix)

        if (associated(Me%StatesLMatrix)) then
            deallocate(Me%StatesLMatrix, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'KillAssimPreProcessor - ModuleAssimilationPreProcessor - ERR09'
            nullify(Me%StatesLMatrix)
        endif

        !Kill all HDF5Files from the list of the ones relevant
        HDF5FileX=> Me%FirstHydroCovHDF5File
        
        do while(associated(HDF5FileX))  
            HDF5FileToKill  => HDF5FileX
            HDF5FileX       => HDF5FileX%Next
            call KillIndividualFullHDF5File(HDF5FileToKill)
        end do
        nullify(Me%FirstHydroCovHDF5File)
        nullify(Me%LastHydroCovHDF5File)

        HDF5FileX=> Me%FirstWaterCovHDF5File
        
        do while(associated(HDF5FileX))  
            HDF5FileToKill  => HDF5FileX
            HDF5FileX       => HDF5FileX%Next
            call KillIndividualFullHDF5File(HDF5FileToKill)
        end do
        nullify(Me%FirstWaterCovHDF5File)
        nullify(Me%LastWaterCovHDF5File)

        !Kill all properties from the list
        PropertyX=> Me%FirstProperty

        do while(associated(PropertyX))  

            PropertyToKill  => PropertyX
            PropertyX       => PropertyX%Next
            call KillIndividualProperty(PropertyToKill)
        end do
        nullify(Me%FirstProperty)

        deallocate(Me)
        nullify(Me)

        !------------------------------------------------------------------------

    end subroutine KillAssimPreProcessor

    !--------------------------------------------------------------------------

    subroutine KillIndividualFullHDF5File(HDF5ToDispose)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: HDF5ToDispose

        !Local-----------------------------------------------------------------
        type(T_CovarianceTime), pointer             :: CovTimeX, CovTimeToKill
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Deallocate the InstantsArray
        if (associated(HDF5ToDispose%InstantsArray)) then 
            deallocate(HDF5ToDispose%InstantsArray, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillIndividualFullHDF5File - ModuleAssimilationPreProcessor - ERR01' 
            nullify(HDF5ToDispose%InstantsArray)
        end if

        !Deallocate the covariance time list
        CovTimeX=> HDF5ToDispose%FirstInstantTime
        
        do while(associated(CovTimeX))  
            CovTimeToKill   => CovTimeX
            CovTimeX        => CovTimeX%Next
            deallocate(CovTimeToKill, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillIndividualFullHDF5File - ModuleAssimilationPreProcessor - ERR02' 
            nullify(CovTimeToKill)
        end do

        deallocate(HDF5ToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'KillIndividualFullHDF5File - ModuleAssimilationPreProcessor - ERR03'
        nullify(HDF5ToDispose)

    end subroutine KillIndividualFullHDF5File

    !--------------------------------------------------------------------------    

    subroutine KillIndividualProperty(PropertyToDispose)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyToDispose

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !The fields of the parameters have already been killed
        deallocate(PropertyToDispose, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                        &
        stop 'KillIndividualProperty - ModuleAssimilationPreProcessor - ERR01'
        nullify(PropertyToDispose)

    end subroutine KillIndividualProperty

    !------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

end module ModuleAssimilationPreProcessor

