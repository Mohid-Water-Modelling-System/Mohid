!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : SequentialAssimilation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2007
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module which manages the sequential data assimilation operations
!
!------------------------------------------------------------------------------
!
!DataFile
!   METHOD                      : int               [1]         !Sequential data assimilation
!                                                               !method: 1= SEEK
!
!   DT                          : real              [-]         !Time step for assimilation
!                                                               !(seconds) starting at initial
!                                                               !time
!
!   INITIAL_STATECOV_FILE       : char              [-]         !Path to the initial covariance
!                                                               !HDF5 file
!
!   STATECOV_EVOLUTION          : 0/1               [0]         !Evolve (1) or not (0) the
!                                                               !covariance structure over time
!
!   STATECOV_EVOLUTION_SPLIT    : 0/1               [1]         !Calculate new EOF every model
!                                                               !split time step 
!
!   STATECOV_RANK               : int               [-]         !Rank of the initial covariance
!                                                               !subspace
!
!   LINEARIZATION_FACTOR        : real              [0.5]       !Model linearization factor used
!                                                               !to evolve the EOF in SEEK
!
!   FORGETTING_FACTOR           : real              [1.0]       !Forgetting factor model present
!                                                               !and past states: 1= perfect model
!
!   OBJECTIVE_ANALYSIS          : 0/1               [0]         !Perform (1) or not (0) objective
!                                                               !analysis at initial time
!
!   OUTPUT_TIME                 : sec. sec. sec.    [-]         !Output time for analysis state in HDF5: 
!                                                               !first, second and interval
!
!   EOF_OUTPUT_TIME             : sec. sec. sec.    [-]         !Output time for EOFs in HDF5:
!                                                               !first, second and interval
!
!   ERROR_TIME_SERIE            : 0/1               [0]         !Output time serie with forecast and
!                                                               !analysis error
!
!   TIME_SERIE                  : 0/1               [0]         !Output time series
!
!   TIME_SERIE_LOCATION         : char              [-]         !Path to time serie location file
!
!   RENORMALIZATION             : 0/1               [0]         !Perform renormalization of correction
!                                                               !basis after each analysis
!
!   REORTHONORMALIZATION        : 0/1               [0]         !Perform periodical reorthonormalization 
!                                                               !of correction basis (after evolution)
!
!   REORTHONORMALIZATION_DT     : real              [-]         !Time step for reorthonormalization of
!                                                               !correction basis (after evolution step)
!
!   <beginproperty>                                             !Block of assimilation state property
!   NAME                        : char              [-]         !Property name
!   UNITS                       : char              [-]         !Property units
!   DIMENSION                   : 2D/3D             [-]         !Property rank (consistent with MOHID
!                                                               !convention)
!   TYPE_ZUV                    : Z/U/V             [Z]         !Type of grid where property is defined: 
!                                                               !Z = cell center, 
!                                                               !U = cell faces U, V = cell faces V
!   STATE_WINDOW                : 4/6*int   [4/6*FillValueInt]  !Spatial window for state:
!                                                               !ilb, jlb, iub, jub (2D), klb, kub (3D)
!   MEASURE                     : 0/1               [0]         !Consider measurements (1) or not for
!                                                               !the state property
!
!   <<begin_measure>>                                           !Block of measurement
!   FILE_IN_TIME                : char          [timeserie]     !Type of file containing measured values
!   DT                          : real              [-]         !Time step of measurements
!   VARIANCE                    : real              [-]         !Variance of measurement error (constant)
!   TIME_TOLERANCE              : real              [0.0]       !Tolerance assumed between assimilation
!                                                               !and measurement times (past and future)
!   FILENAME                    : char              [-]         !Path to measurement file
!   DATA_COLUMN                 : int               [-]         !Measurement column in time serie 
!                                                               !measurement file
!   LOCALIZATION_I              : int               [-]         !Time serie localization i
!   LOCALIZATION_J              : int               [-]         !Time serie localization j
!   LOCALIZATION_K              : int               [-]         !Time serie localization k
!   <<end_measure>>
!   <endproperty>


Module ModuleSequentialAssimilation

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData,         only : GetGridData, UngetGridData
    use ModuleHorizontalGrid,   only : WriteHorizontalGrid, GetXYCellZ
    use ModuleTimeSerie,        only : StartTimeSerieInput, GetTimeSerieValue,          &
                                       KillTimeSerie, StartTimeSerie,                   &
                                       GetNumberOfTimeSeries, GetTimeSerieLocation,     &
                                       CorrectsCellsTimeSerie, WriteTimeSerie,          &
                                       WriteTimeSerieLine
    use ModuleHDF5
    use ModuleHydrodynamic,     only : GetWaterLevel, GetHorizontalVelocity,            &
                                       GetOldHorizontalVelocity,                        &
                                       GetWaterFluxes, GetVerticalVelocity,             &
                                       GetSubModelFluxes, GetChezyVelUV,                &
                                       GetHydroNeedsFather, CopyChezyVelUV,             &
                                       CopyWaterLevel, CopyHorizontalVelocity,          &
                                       CopyVerticalVelocity, CopyWaterFluxes,           &
                                       CopySubModelFluxes, UnGetHydrodynamic,           &
                                       SetWaterLevel, SetHorizontalVelocity,            &
                                       SetVerticalVelocity, SetWaterFluxes,             &
                                       SetSubModelFluxes, SetChezyVelUV,                &
                                       ReSetHydrodynamicProperties,                     &
                                       GetCyclicBoundary, GetHydroSeqAssimilation,      &
                                       PointToHydroState, NullifyHydroStatePointer
    use ModuleWaterProperties,  only : GetConcentration, CopyConcentration,             &
                                       SetConcentration, SetDensity, SetSigma,          &
                                       ReSetConcentration, ReSetDensity, CopyDensity,   &
                                       GetWaterPropertiesNumber,                        &
                                       GetWaterPropertiesIDArray, GetDensity, GetSigma, &
                                       UnGetWaterProperties, GetWaterSeqAssimilation,   &
                                       ConstructPropertiesIDArray, PointToDensity,      &
                                       PointToConcentration, NullifyDensityPointer,     &
                                       NullifyConcentrationPointer
    use ModuleFunctions,        only : SetMatrixValue, ConstructPropertyID,             &
                                       InvSingularDiagMatrix2D, CholLinSystemSolver,    &
                                       CholeskyFactorization, EOFAnalysis,              &
                                       Check_Hydrodynamic_Property, Check_Water_Property
    use ModuleGeometry,         only : GetGeometrySize, GetGeometryDistances,           &
                                       GetGeometryAreas, GetGeometryVolumes,            &
                                       GetGeometryWaterColumn, CopyGeometryDistances,   &
                                       CopyGeometryAreas, CopyGeometryVolumes,          &
                                       CopyGeometryWaterColumn, SetGeometryDistances,   &
                                       SetGeometryAreas, SetGeometryWaterColumn,        &
                                       SetGeometryVolumes, UnGetGeometry,               &
                                       ReSetGeometry, PointToGeometryState,             &
                                       NullifyGeometryStatePointer, GetLayer4Level
    use ModuleMap,              only : GetWaterPoints3D, GetOpenPoints3D, GetWetFaces,  &
                                       GetComputeFaces3D, GetImposedNormalFaces,        &
                                       GetImposedTangentialFaces, UnGetMap 

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartSequentialAssimilation
    private ::      AllocateInstance
    private ::      ConstructSequentialAssimilation
    private ::          Read_SeqAssimilation_Files_Name
    private ::          ConstructOptions
    private ::              ConstructProperty
    private ::                  ConstructMeasurement
    private ::                  ConstructMeasureStateLoc
    private ::                  AddMeasurement
    private ::              AddProperty
    private ::          ConstructLog
    private ::          ConstructSeqAssimVariables
    private ::              AllocateVariables
    private ::                  AddWaterProperty
    private ::              ConstructInitialStateCov
    private ::          ObjectiveAnalysis
    private ::              InitialUCalculationSEEK
    private ::              ConstructModelState
    private ::                  ActualizeOpenFaces
    private ::                  AddPropertyToState
    private ::          ConstructOutPutTime
    private ::          OpenHDF5OutPutFile

    !Selector
    public  :: GetSeqAssimilationOptions
    public  :: GetSeqAssimilationTime

    !Modifier
    public  :: SetModelInitialState
    private ::      CopyModelResultToFS
    private ::      CopyFullModelState
    private ::      DisturbFullModelState
    private ::      SetDisturbedStateInModel
    private ::      CopyUndisturbedStateToModel
    public  :: GetModelResult
    public  :: CovarianceCalculationSEEK
    public  :: ModifySequentialAssimilation
    private ::      ReadMeasures
    private ::      UCalculationSEEK
    private ::      SEEKAnalysis
    private ::      CopyModelAnalysedState
    private ::      SeqAssimilationOutPut

    !Destructor
    public  :: KillSequentialAssimilation
    private ::      Deassociate_External_Modules
    private ::      DeallocateVariables
    private ::          KillFullState
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjSequentialAssimilation 

    !Interfaces----------------------------------------------------------------

    !Parameter-----------------------------------------------------------------

    !Data assimilation methods
    integer, parameter                              :: SEEK                 = 1

    !Block
    character(LEN = StringLength), parameter        :: block_begin          = '<beginproperty>'
    character(LEN = StringLength), parameter        :: block_end            = '<endproperty>'
    character(LEN = StringLength), parameter        :: begin_measure        = '<<begin_measure>>'
    character(LEN = StringLength), parameter        :: end_measure          = '<<end_measure>>'

    !Property dimensions 
    integer, parameter                              :: Dim_2D               = 2
    integer, parameter                              :: Dim_3D               = 3

    !Direction
    integer, parameter                              :: DirectionX_          = 1
    integer, parameter                              :: DirectionY_          = 2

    !Measurement file type
    integer, parameter                              :: TimeSerie            = 1
    integer, parameter                              :: ProfileTimeSerie     = 2
    integer, parameter                              :: HDF                  = 3

    !Types---------------------------------------------------------------------

    type T_TimeSerie
        integer                                     :: ID                   = 0
        character(PathLength)                       :: FileName             = null_str
        integer                                     :: DataColumn           = null_int
    end type T_TimeSerie

    type T_Localization
        integer                                     :: I
        integer                                     :: J
        integer                                     :: K
    end type T_Localization

    private :: T_Property
    type       T_Property
        type (T_PropertyID)                         :: ID                   = null_int
        integer                                     :: Dim                  = null_int
        integer                                     :: TypeZUV              = null_int
        
        type (T_Size2D)                             :: Window2D
        type (T_Size3D)                             :: Window
        
        integer                                     :: ModuleType           = null_int
        logical                                     :: Measure              = .false.
        integer                                     :: MeasuresNumber       = 0
        integer                                     :: FirstStatePosition   = null_int !initialization: Carina

        real, dimension (:, :), pointer             :: Field2D              => null()
        real, dimension (:, :, :), pointer          :: Field                => null()
        real(8), dimension (:, :, :), pointer       :: FieldR8              => null()
        real, dimension (:, :, :), pointer          :: CenteredField        => null()
        real(8), dimension (:, :, :), pointer       :: CenteredFieldR8      => null()

        !logical                                     :: CyclicBoundary       = .false.

        type (T_Property), pointer                  :: Next           => null()
		type (T_Property), pointer                  :: Prev           => null()
    end type T_Property

    private :: T_Measure
    type       T_Measure
        type (T_PropertyID)                         :: ID                   = null_int  !initialization: Carina
        integer                                     :: FileType             = null_int  !initialization: Carina
        real                                        :: DT                   = null_real !initialization: Carina
        real                                        :: Variance             = null_real !initialization: Carina
        integer                                     :: StateLocation        = FillValueInt
        
        type (T_TimeSerie)                          :: TimeSerie

        type (T_Localization)                       :: Localization

        character(PathLength)                       :: FileName            = null_str  !initialization: Carina

        logical                                     :: ValueExists         = .false.   !initialization: Carina
        real                                        :: Tolerance           = null_real !initialization: Carina

        type (T_Measure), pointer                   :: Next           => null()
		type (T_Measure), pointer                   :: Prev           => null()
    end type T_Measure

    private :: T_StateProp
    type       T_StateProp
        integer                                     :: IDNumber           = null_int  !initialization: Carina
        real, dimension (:, :, :), pointer          :: Field              => null()  !initialization: Carina

        type (T_StateProp), pointer                 :: Next         => null()
        type (T_StateProp), pointer                 :: Prev         => null()
 end type T_StateProp

    private :: T_FullState
    type       T_FullState
        !ModuleHydrodynamic properties
        real,    dimension (:, :),    pointer       :: WaterLevel         => null()  !initialization: Carina
        real,    dimension (:, :, :), pointer       :: VelocityU          => null()  !initialization: Carina
        real,    dimension (:, :, :), pointer       :: VelocityV          => null()  !initialization: Carina
        real,    dimension (:, :, :), pointer       :: VelocityUOld       => null()  !initialization: Carina
        real,    dimension (:, :, :), pointer       :: VelocityVOld       => null()  !initialization: Carina
        real,    dimension (:, :, :), pointer       :: VelocityAcross     => null()  !initialization: Carina
        real(8), dimension (:, :, :), pointer       :: WaterFluxX         => null()  !initialization: Carina
        real(8), dimension (:, :, :), pointer       :: WaterFluxY         => null()  !initialization: Carina
        real(8), dimension (:, :, :), pointer       :: WaterFluxZ         => null()  !initialization: Carina
        real(8), dimension(:,:,:)   , pointer       :: SubModelqX         => null()  !initialization: Carina
        real(8), dimension(:,:,:)   , pointer       :: SubModelqY         => null()  !initialization: Carina
        real,    dimension(:,:)     , pointer       :: ChezyVelUV         => null()  !initialization: Carina

        !ModuleWaterProperties
        real, pointer, dimension(:,:,:)             :: Density            => null()  !initialization: Carina
        real, pointer, dimension(:,:,:)             :: SigmaDensity       => null()  !initialization: Carina
        type (T_StateProp), pointer                 :: FirstWaterProperty => null()  !initialization: Carina
        type (T_StateProp), pointer                 :: LastWaterProperty  => null()  !initialization: Carina

        !ModuleGeometry properties
        real,    dimension(:, :, :), pointer        :: SZZ                => null()  !initialization: Carina
        real,    dimension(:, :, :), pointer        :: DWZ                => null()  !initialization: Carina
		real,    dimension(:, :, :), pointer        :: DUZ                => null()  !initialization: Carina
		real,    dimension(:, :, :), pointer        :: DVZ                => null()  !initialization: Carina
        real,    dimension(:, :, :), pointer        :: DZZ                => null()  !initialization: Carina
        real,    dimension(:, :, :), pointer        :: AreaU              => null()  !initialization: Carina
		real,    dimension(:, :, :), pointer        :: AreaV              => null()  !initialization: Carina
        real,    dimension(:, :, :), pointer        :: ZCellCenter        => null()  !initialization: Carina
        real,    dimension(:, :), pointer           :: WaterColumnU      => null()  !initialization: Carina
		real,    dimension(:, :), pointer           :: WaterColumnV       => null()  !initialization: Carina
        real,    dimension(:, :), pointer           :: WaterColumnZ       => null()  !initialization: Carina
        real(8), dimension(:, :, :), pointer        :: VolumeZ           => null()  !initialization: Carina
		real(8), dimension(:, :, :), pointer        :: VolumeU           => null()  !initialization: Carina
        real(8), dimension(:, :, :), pointer        :: VolumeV           => null()  !initialization: Carina
		real(8), dimension(:, :, :), pointer        :: VolumeV           => null()  !initialization: Carina
    end type T_FullState

    private :: T_Files
    type       T_Files
         character(len=StringLength)                :: ConstructData = null_str !initialization: Carina
         character(len=StringLength)                :: OutPutFields  = null_str !initialization: Carina
    end type T_Files

    private :: T_OutPut
    type       T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime => null()  !initialization: Carina
         integer                                    :: NextOutPut = null_int  !initialization: Carina
		 integer                                    :: Number = null_int  !initialization: Carina
         logical                                    :: HDF5ON = .false.
         logical                                    :: TimeSerieON = .false.
         logical                                    :: ErrorTimeSerieON = .false.
    end type T_OutPut

    type       T_ErrorTimeSerie
        integer                                     :: TimeSerie = 0
        type (T_ErrorTimeSerie), pointer            :: Next
    end type   T_ErrorTimeSerie

    !type       T_CyclicBoundary
    !    logical                                     :: ON
    !    integer                                     :: Direction
    !end type T_CyclicBoundary

    private :: T_External
    type T_External
        integer, dimension(:,:,:), pointer          :: WaterPoints          => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: OpenPoints           => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: ComputeFacesU        => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: ComputeFacesV        => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: ImposedNormFacesU    => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: ImposedNormFacesV    => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: ImposedTangFacesU    => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: ImposedTangFacesV    => null()  !initialization: Carina
    end type T_External

    private :: T_SequentialAssimilation
    type       T_SequentialAssimilation
        integer                                     :: InstanceID = null_int  !initialization: Carina
        type (T_Size3D)                             :: Size, WorkSize

        type(T_Files)                               :: Files
        type(T_Time)                                :: CurrentTime
        type(T_Time)                                :: EndTime
        type(T_OutPut)                              :: OutPut
        type(T_OutPut)                              :: EOF_OutPut
        type(T_ErrorTimeSerie), pointer             :: FirstErrorTimeSerie
        !type(T_CyclicBoundary)                      :: CyclicBoundary

        logical                                     :: HydroSeqAssim = .false.
        logical                                     :: WaterSeqAssim = .false.

        !Data assimilation options
        type (T_Time)                               :: AssimilationTime
        integer                                     :: Method = null_int  !initialization: Carina
        real                                        :: DT     = null_real  !initialization: Carina
        character(len=StringLength)                 :: InitialStateCovFile = null_str  !initialization: Carina
        integer                                     :: StateCovRank = null_int  !initialization: Carina
        logical                                     :: ObjectiveAnalysis = .false. !initialization: Carina
        logical                                     :: MeasureErrorVar = .true.
        real                                        :: ForgettingFactor   = null_real  !initialization: Carina
        logical                                     :: StateCovEvolution = .false. !initialization: Carina
        real                                        :: StateCov_LinFactor  = null_real  !initialization: Carina
        logical                                     :: StateCov_Evolution_Split = .false. !initialization: Carina
        logical                                     :: ReNormalization = .false. !initialization: Carina
        logical                                     :: ReOrthonormalization = .false. !initialization: Carina
        real                                        :: ReOrthoNorm_DT = null_real  !initialization: Carina
        type (T_Time)                               :: ReOrthonormTime

        !Data assimilation variables
        real(8), dimension (:, :), pointer          :: State            => null()  !initialization: Carina
        integer, dimension (:, :), pointer          :: OpenPointState   => null()  !initialization: Carina
        real(8), dimension (:, :), pointer          :: MeasuredState    => null()  !initialization: Carina
        real,    dimension (:, :), pointer          :: MeasuresInvCov   => null()  !initialization: Carina
        real,    dimension (:, :), pointer          :: ObservOperator   => null()  !initialization: Carina
        !real, dimension (:, :), pointer              :: KalmanGain
        real(8), dimension (:, :), pointer          :: LMatrix          => null()  !initialization: Carina
        real, dimension (:, :), pointer             :: InvUMatrix       => null()  !initialization: Carina

        !Full state variables (for model covariance evolution)
        type (T_FullState)                          :: FS
        type (T_FullState), dimension (:), pointer  :: DFS => null()  !initialization: Carina
        integer                                     :: FullWaterPropNumber
        integer, dimension(:), pointer              :: PropertiesIDArray   => null()  !initialization: Carina
        logical                                     :: HydroSubModel       = .false.
        logical                                     :: PerturbeState       = .true.
        !logical                                     :: CalculateEOF         = .false.

        !Property list
        type(T_Property), pointer                   :: FirstProperty  => null()  !initialization: Carina
        type(T_Property), pointer                   :: LastProperty  => null()  !initialization: Carina
        integer                                     :: PropertiesNumber     = FillValueInt
        integer                                     :: StateVarNumber       = 0 a
        logical                                     :: PropertiesInFaces    = .false.

        !Measurement list
        type(T_Measure), pointer                    :: FirstMeasure  => null()  !initialization: Carina
        type(T_Measure), pointer                    :: LastMeasure   => null()  !initialization: Carina
        integer                                     :: MeasuresNumber       = FillValueInt
        integer                                     :: CurrentMeasuresNumber = null_int  !initialization: Carina
        integer                                     :: PreviousMeasuresNumber = null_int  !initialization: Carina

        !Mapping variables
        type(T_External)                            :: External_Var
        integer, dimension(:,:,:), pointer          :: WaterFacesU  => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: WaterFacesV  => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: OpenFacesU   => null()  !initialization: Carina
        integer, dimension(:,:,:), pointer          :: OpenFacesV   => null()  !initialization: Carina

        !Instance of other modules
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjCovHDF5           = 0
        integer                                     :: ObjHydrodynamic      = 0
        integer                                     :: ObjWaterProperties   = 0
        integer                                     :: ObjGeometry          = 0
        integer                                     :: ObjMap               = 0
        integer                                     :: ObjHDF5              = 0
        integer                                     :: ObjGridData          = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjTimeSerie         = 0

        type(T_SequentialAssimilation), pointer     :: Next                     => null()
    end type  T_SequentialAssimilation

    !Global Module Variables
    type (T_SequentialAssimilation), pointer        :: FirstObjSeqAssimilation  => null()
    type (T_SequentialAssimilation), pointer        :: Me                       => null()

    !--------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartSequentialAssimilation(SequentialAssimilationID,            &
                                           GridDataID, HorizontalGridID,        &
                                           TimeID,                              &
                                           HydrodynamicID, WaterPropertiesID,   &
                                           GeometryID, MapID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: SequentialAssimilationID
        integer                                         :: GridDataID
        integer                                         :: HorizontalGridID
        integer                                         :: TimeID
        integer                                         :: HydrodynamicID
        integer                                         :: WaterPropertiesID
        integer                                         :: GeometryID
        integer                                         :: MapID
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSequentialAssimilation_)) then
            nullify (FirstObjSeqAssimilation)
            call RegisterModule (mSequentialAssimilation_) 
        endif

        call Ready(SequentialAssimilationID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Associates External Instances
            Me%ObjTime            = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjGridData        = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalGrid  = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID) 
            Me%ObjHydrodynamic    = AssociateInstance (mHYDRODYNAMIC_,   HydrodynamicID  )
            Me%ObjWaterProperties = AssociateInstance (mWATERPROPERTIES_,               &
                                                       WaterPropertiesID )
            Me%ObjGeometry        = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap             = AssociateInstance (mMAP_,            MapID           )

            call ConstructSequentialAssimilation

            !Returns ID
            SequentialAssimilationID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0

            stop 'ModuleSequentialAssimilation - StartSequentialAssimilation - ERR01' 

        end if cd0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartSequentialAssimilation
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_SequentialAssimilation), pointer                         :: NewObjSeqAssimilation
        type (T_SequentialAssimilation), pointer                         :: PreviousObjSeqAssimilation

        !------------------------------------------------------------------------

        !Allocates new instance
        allocate (NewObjSeqAssimilation)
        nullify  (NewObjSeqAssimilation%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSeqAssimilation)) then
            FirstObjSeqAssimilation         => NewObjSeqAssimilation
            Me                    => NewObjSeqAssimilation
        else
            PreviousObjSeqAssimilation      => FirstObjSeqAssimilation
            Me                    => FirstObjSeqAssimilation%Next
            do while (associated(Me))
                PreviousObjSeqAssimilation  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjSeqAssimilation
            PreviousObjSeqAssimilation%Next => NewObjSeqAssimilation
        endif

        Me%InstanceID = RegisterNewInstance (mSequentialAssimilation_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ConstructSequentialAssimilation

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL

        !------------------------------------------------------------------------

        !Get size
        call GetGeometrySize(Me%ObjGeometry, Size = Me%Size,                    &
                             WorkSize = Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR01'

        !Get times
        call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR02'

        call GetComputeTimeLimits(Me%ObjTime, EndTime = Me%EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR03'

        !Get WaterPoints
        call GetWaterPoints3D(Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR04'

        !!Get cyclic boundary
        !call GetCyclicBoundary(Me%ObjHydrodynamic, Me%CyclicBoundary%ON,        &
        !                       Me%CyclicBoundary%Direction, STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)                                              &
        !    stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR04b'

        call Read_SeqAssimilation_Files_Name

        !Construct enter data 
        call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData,        &
                            STAT = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR05'

        call ConstructOptions

        if (Me%MeasuresNumber > 0) then

            call ConstructLog

            !Construct data assimilation variables
            call ConstructSeqAssimVariables

            call ConstructOutPutTime

            !Opens the seq. assimilation results HDF5 file
            if (Me%OutPut%HDF5ON .or. Me%EOF_OutPut%HDF5ON) call OpenHDF5OutPutFile

            !Construct the results Time Serie
            if (Me%OutPut%TimeSerieON) call ConstructTimeSerie

            !Construct the error Time Serie
            if (Me%OutPut%ErrorTimeSerieON) call ConstructErrorTimeSerie

            call UnGetMap (Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR06'

            !Perform objective analysis
            if (Me%ObjectiveAnalysis) then
                call ObjectiveAnalysis
            endif

            !Initialize Me%AssimilationTime
            Me%AssimilationTime = Me%CurrentTime + Me%DT

            if (Me%ReOrthonormalization) then

                Me%ReOrthonormTime = Me%CurrentTime + Me%ReOrthonorm_DT

            endif

        else

            write(*,*)  
            write(*,*) 'No measurements are specified for sequential assimilation.'
            stop 'ConstructSequentialAssimilation - ModuleSequentialAssimilation - ERR07'

        endif

    end subroutine ConstructSequentialAssimilation

    !--------------------------------------------------------------------------

    subroutine ConstructWaterFaces

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: auxi, auxj
        integer                                         :: i, j, k
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB

        !------------------------------------------------------------------------

        !Calculates the WaterFaces based on the WaterPoints
        !(cell faces confining with two WaterPoints or boundary faces)

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        allocate(Me%WaterFacesU(ILB:IUB,JLB:JUB + 1,KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructWaterFaces - ModuleSequentialAssimilation - ERR01'

        allocate(Me%WaterFacesV(ILB:IUB + 1,JLB:JUB,KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructWaterFaces - ModuleSequentialAssimilation - ERR02'

        allocate(Me%OpenFacesU(ILB:IUB,JLB:JUB + 1,KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructWaterFaces - ModuleSequentialAssimilation - ERR03'

        allocate(Me%OpenFacesV(ILB:IUB + 1,JLB:JUB,KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructWaterFaces - ModuleSequentialAssimilation - ERR04'

        !By default all values are zero
        Me%WaterFacesU      = 0
        Me%WaterFacesV      = 0
        Me%OpenFacesU       = 0
        Me%OpenFacesV       = 0

        auxj = JLB
        auxi = ILB

        do k = KLB, KUB

            !WaterFacesU

            !Boundary faces (JLB, JUB+1)
            auxj = JLB
            do i = ILB, IUB
                if (Me%External_Var%WaterPoints(i,auxj,k) == 1) then

                    Me%WaterFacesU(i,auxj,k) = 1
                endif
            enddo
            auxj = JUB
            do i = ILB, IUB
                if (Me%External_Var%WaterPoints(i,auxj,k) == 1) then

                    Me%WaterFacesU(i,auxj + 1,k) = 1
                endif
            enddo
            !Interior faces
            do j = JLB + 1, JUB
                do i = ILB, IUB

                    if ((Me%External_Var%WaterPoints(i,j-1,k) == 1) .and.           &
                        (Me%External_Var%WaterPoints(i,j,k) == 1)) then

                        Me%WaterFacesU(i,j,k) = 1
                    endif
                enddo
            enddo

            !do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            !    if (Me%External_Var%WaterPoints(i,auxj,k) == 1) then

            !        Me%WaterFacesU(i,auxj,k) = 1
            !    endif
            !enddo

            !do j = Me%WorkSize%JLB + 1, Me%WorkSize%JUB + 1

            !    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            !        if ((Me%External_Var%WaterPoints(i,j-1,k) == 1) .or.            &
            !            (Me%External_Var%WaterPoints(i,j,k) == 1)) then

            !            Me%WaterFacesU(i,j,k) = 1
            !        endif
            !    enddo
            !enddo

            !WaterFacesV

            do j = JLB, JUB
                !Boundary faces
                auxi = ILB
                if (Me%External_Var%WaterPoints(auxi,j,k) == 1) then

                    Me%WaterFacesV(auxi,j,k) = 1
                endif
                auxi = IUB
                if (Me%External_Var%WaterPoints(auxi,j,k) == 1) then

                    Me%WaterFacesV(auxi + 1,j,k) = 1
                endif
                !Interior faces
                do i = ILB + 1, IUB

                    if ((Me%External_Var%WaterPoints(i-1,j,k) == 1) .and.           &
                        (Me%External_Var%WaterPoints(i,j,k) == 1)) then

                        Me%WaterFacesV(i,j,k) = 1
                    endif
                enddo
            enddo

            !do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            !    if (Me%External_Var%WaterPoints(auxi,j,k) == 1) then
                   
            !        Me%WaterFacesV(auxi,j,k) = 1
            !    endif
            !    do i = Me%WorkSize%ILB + 1, Me%WorkSize%IUB + 1

            !        if ((Me%External_Var%WaterPoints(i-1,j,k) == 1) .or.            &
            !            (Me%External_Var%WaterPoints(i,j,k) == 1)) then

            !            Me%WaterFacesV(i,j,k) = 1
            !        endif
            !    enddo
            !enddo
        enddo

    end subroutine ConstructWaterFaces

    !--------------------------------------------------------------------------

    subroutine Read_SeqAssimilation_Files_Name

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        character(len = StringLength)                   :: Message
        logical                                         :: exist

        !------------------------------------------------------------------------

        !Opens the Sequential Assimilation data file 
        ! ---> ASCII file used to construct a new sequential assimilation
        Message   ='ASCII file used to construct a new sequential assimilation.'
        Message   = trim(Message)

        call ReadFileName('SEQASSIM_DAT', Me%Files%ConstructData,               &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'Read_SeqAssimilation_Files_Name - ModuleSequentialAssimilation - ERR01' 

        inquire(FILE   = Me%Files%ConstructData,                                &
                EXIST  = exist,                                                 &
                IOSTAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Read_SeqAssimilation_Files_Name - ModuleSequentialAssimilation - ERR02' 
        if (.NOT. exist)                                                        &
            stop 'Read_SeqAssimilation_Files_Name - ModuleSequentialAssimilation - ERR03' 

        ! ---> File in HDF format where is written instant fields of assimilation analysis
        Message   ='Instant fields of state properties analysis in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SEQASSIM_HDF', Me%Files%OutPutFields,                &
                           Message = Message, TIME_END = Me%EndTime,            &
                           Extension = 'sqt', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'Read_SeqAssimilation_Files_Name - ModuleSequentialAssimilation - ERR04' 

    end subroutine Read_SeqAssimilation_Files_Name

    !--------------------------------------------------------------------------

    subroutine ConstructOptions

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------
        logical                                         :: exist

        !Local-------------------------------------------------------------------
        integer                                         :: FromFile
        integer                                         :: iflag
        integer                                         :: ClientNumber
        integer                                         :: STAT_CALL
        logical                                         :: BlockFound
        type(T_Property),   pointer                     :: NewProperty

        !------------------------------------------------------------------------

        !Check if sequential data assimilation is commanded in modules
        !Hydrodynamic
        call GetHydroSeqAssimilation(Me%ObjHydrodynamic, Me%HydroSeqAssim, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR01'
        !WaterProperties
        call GetWaterSeqAssimilation(Me%ObjWaterProperties, Me%WaterSeqAssim, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR02'

        !Get number of water properties specified in ModuleWaterProperties
        call GetWaterPropertiesNumber(Me%ObjWaterProperties,                    &
                                      Me%FullWaterPropNumber, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR03'

        if (Me%WaterSeqAssim) then
            if (Me%FullWaterPropNumber > 0) then
                !Construct water properties IDs array
                call ConstructPropertiesIDArray(Me%ObjWaterProperties, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'ConstructOptions - ModuleSequentialAssimilation - ERR04'

                !Obtain IDs of water properties in ModuleWaterProperties
                call GetWaterPropertiesIDArray(Me%ObjWaterProperties,           &
                                               Me%PropertiesIDArray, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'ConstructOptions - ModuleSequentialAssimilation - ERR05'
            else
                write(*,*)  
                write(*,*) 'No water properties specified in ModuleWaterProperties but'
                write(*,*) 'Sequential Assimilation is commanded in ModuleWaterProperties.'
                stop 'ConstructOptions - ModuleSequentialAssimilation - ERR06'
            endif
        endif

        call GetExtractType(FromFile = FromFile)

        call GetData(Me%Method,                                                 &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'METHOD',                                     &
                     Default    = SEEK,                                         &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR07'

        if (Me%Method /= SEEK) then
            write(*,*)  
            write(*,*) 'Data assimilation method not implemented.'
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR08'
        endif

        call GetData(Me%DT,                                                     &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'DT',                                         &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR09'

        call GetData(Me%InitialStateCovFile,                                    &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'INITIAL_STATECOV_FILE',                      &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR10'

        inquire (file=trim(Me%InitialStateCovFile), exist = exist)
        if (.not. exist) then
            write(*,*)
            write(*,*)'Could not find covariance file '//trim(Me%InitialStateCovFile)
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR11'
        endif

        call GetData(Me%StateCovRank,                                           &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'STATECOV_RANK',                              &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR12'

        call GetData(Me%ObjectiveAnalysis,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'OBJECTIVE_ANALYSIS',                         &
                     Default    = .false.,                                      &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR13'

        if (Me%Method == SEEK) then
            call GetData(Me%ForgettingFactor,                                   &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromFile,                                 &
                         keyword    = 'FORGETTING_FACTOR',                      &
                         default    = 1.0,                                      &
                         ClientModule = 'ModuleSequentialAssimilation',         &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructOptions - ModuleSequentialAssimilation - ERR14'
        endif

        call GetData(Me%StateCovEvolution,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'STATECOV_EVOLUTION',                         &
                     Default    = .false.,                                      &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR15'

        call GetData(Me%Output%TimeSerieON,                                     &
                     Me%ObjEnterData, iflag,                                    &
                     keyword    = 'TIME_SERIE',                                 &
                     Default    = .false.,                                      &
                     SearchType = FromFile,                                     &
                     ClientModule ='ModuleSequentialAssimilation',              &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR15b'

        call GetData(Me%Output%ErrorTimeSerieON,                                &
                     Me%ObjEnterData, iflag,                                    &
                     keyword      ='ERROR_TIME_SERIE',                          &
                     Default    = .false.,                                      &
                     SearchType   = FromFile,                                   &
                     ClientModule ='ModuleSequentialAssimilation',              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR15c'

        if (Me%StateCovEvolution .and. (Me%Method == SEEK)) then

            call GetData(Me%StateCov_LinFactor,                                 &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromFile,                                 &
                         keyword    = 'LINEARIZATION_FACTOR',                   &
                         Default    = 0.5,                                      &
                         ClientModule = 'ModuleSequentialAssimilation',         &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructOptions - ModuleSequentialAssimilation - ERR16'

            call GetData(Me%StateCov_Evolution_Split,                           &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromFile,                                 &
                         keyword    = 'STATECOV_EVOLUTION_SPLIT',               &
                         Default    = .false.,                                   &
                         ClientModule = 'ModuleSequentialAssimilation',         &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructOptions - ModuleSequentialAssimilation - ERR16a'

            call GetData(Me%ReNormalization,                                    &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromFile,                                 &
                         keyword    = 'RENORMALIZATION',                        &
                         Default    = .false.,                                  &
                         ClientModule = 'ModuleSequentialAssimilation',         &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructOptions - ModuleSequentialAssimilation - ERR16b'

            call GetData(Me%ReOrthonormalization,                               &
                         Me%ObjEnterData, iflag,                                &
                         SearchType = FromFile,                                 &
                         keyword    = 'REORTHONORMALIZATION',                   &
                         Default    = .false.,                                  &
                         ClientModule = 'ModuleSequentialAssimilation',         &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructOptions - ModuleSequentialAssimilation - ERR16c'

            if (Me%ReOrthonormalization) then
                call GetData(Me%ReOrthoNorm_DT,                                 &
                             Me%ObjEnterData, iflag,                            &
                             SearchType = FromFile,                             &
                             keyword    = 'REORTHONORMALIZATION_DT',            &
                             ClientModule = 'ModuleSequentialAssimilation',     &
                             STAT       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'ConstructOptions - ModuleSequentialAssimilation - ERR16d'
            endif

        endif

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                        block_begin, block_end, BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                    
                    ! Construct a New Property 
                    call ConstructProperty(NewProperty, ClientNumber)

                    ! Add new Property to the Assimilation List 
                    call AddProperty     (NewProperty)
                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)                                 &
                        stop 'ConstructOptions - ModuleSequentialAssimilation - ERR17'
                        
                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)                                     &
                    stop 'ConstructOptions - ModuleSequentialAssimilation - ERR18'
            end if cd1
        end do do1

        if (Me%WaterSeqAssim .and. (Me%FullWaterPropNumber > 0)) then
            call UngetWaterProperties(Me%ObjWaterProperties,                    &
                                      Me%PropertiesIDArray, STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
            stop 'ConstructOptions - ModuleSequentialAssimilation - ERR19'
        endif

    end subroutine ConstructOptions

    !--------------------------------------------------------------------------

    subroutine ConstructProperty(NewProperty, ClientNumber)

        !Arguments---------------------------------------------------------------
        type(T_Property),   pointer                     :: NewProperty
        integer                                         :: ClientNumber

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: iflag
        integer                                         :: STAT_CALL
        character(len=StringLength)                     :: Char_TypeZUV
        integer, dimension(:), pointer                  :: aux
        integer                                         :: WorkILB, WorkIUB, WorkJLB
        integer                                         :: WorkJUB, WorkKLB, WorkKUB
        logical                                         :: BlockFound
        type(T_Measure),   pointer                      :: NewMeasure
        integer                                         :: i, j, k
        integer, dimension(:,:,:), pointer              :: PropertyMap  => null()

        !------------------------------------------------------------------------

        !Allocates new property
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR01'

        !Get keywords from property block
        !ID properties
        call ConstructPropertyID (NewProperty%ID, Me%ObjEnterData, FromBlock)

        !The other
        call GetData(NewProperty%Dim, Me%ObjEnterData, iflag,                       &
                     keyword        = 'DIMENSION',                                  &  
                     default        = Dim_3D,                                       &
                     SearchType     = FromBlock,                                    &
                     ClientModule   = 'ModuleSequentialAssimilation',               &
                     STAT           = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR02'

        if (NewProperty%Dim /= Dim_3D) then
            if (NewProperty%Dim /= Dim_2D) then
                write(*,*)  
                write(*,*) 'Property must be either 3D or 2D: ',                    &
                            trim(NewProperty%ID%Name)
                stop 'ConstructProperty - ModuleSequentialAssimilation - ERR03'

            elseif (NewProperty%ID%IDNumber /= WaterLevel_) then
                write(*,*)  
                write(*,*) 'Only property Water Level is registered as 2D in MOHID Water'
                stop 'ConstructProperty - ModuleSequentialAssimilation - ERR04'
            endif
        endif         

        call GetData(Char_TypeZUV, Me%ObjEnterData, iflag,                          &
                     keyword        = 'TYPE_ZUV',                                   &  
                     SearchType     = FromBlock,                                    &
                     ClientModule   = 'ModuleSequentialAssimilation',               &
                     default        = "Z",                                          &
                     STAT           = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR05'

        NewProperty%TypeZUV  = TranslateTypeZUV(Char_TypeZUV)

        if ((NewProperty%Dim /= Dim_3D) .and. (NewProperty%TypeZUV /= TypeZ_))      &
            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR06'

        if (((NewProperty%TypeZUV == TypeU_) .or. (NewProperty%TypeZUV == TypeZ_))  &
           .and. .not. Me%PropertiesInFaces) then
            Me%PropertiesInFaces = .true.
            call ConstructWaterFaces
        endif

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        if (NewProperty%Dim == Dim_2D) then

            allocate (aux(4))

            call GetData(aux,                                                       &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType = FromBlock,                                    &
                         keyword    = 'STATE_WINDOW',                               &
                         Default    = FillValueInt,                                 &
                         ClientModule ='ModuleSequentialAssimilation',              &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ConstructProperty - ModuleSequentialAssimilation - ERR07'

            if (iflag == 4) then

                NewProperty%Window2D%ILB = aux(1)
                NewProperty%Window2D%IUB = aux(2)
                NewProperty%Window2D%JLB = aux(3) 
                NewProperty%Window2D%JUB = aux(4)
            else

                write(*,*)  
                write(*,*) 'Spatial window not specified for property: ',           &
                            trim(NewProperty%ID%Name)
                write(*,*) 'Assumed full domain...'

                NewProperty%Window2D%ILB = WorkILB
                NewProperty%Window2D%IUB = WorkIUB
                NewProperty%Window2D%JLB = WorkJLB 
                NewProperty%Window2D%JUB = WorkJUB
            endif

        else !Dim_3D

            allocate (aux(6))

            call GetData(aux,                                                       &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType = FromBlock,                                    &
                         keyword    = 'STATE_WINDOW',                               &
                         Default    = FillValueInt,                                 &
                         ClientModule ='ModuleSequentialAssimilation',              &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ConstructProperty - ModuleSequentialAssimilation - ERR08'

            if (iflag == 6) then

                NewProperty%Window%ILB = aux(1)
                if (NewProperty%TypeZUV == TypeV_) then
                    NewProperty%Window%IUB = aux(2) + 1
                else
                    NewProperty%Window%IUB = aux(2)
                endif
                NewProperty%Window%JLB = aux(3) 
                if (NewProperty%TypeZUV == TypeU_) then
                    NewProperty%Window%JUB = aux(4) + 1
                else
                    NewProperty%Window%JUB = aux(4)
                endif
                NewProperty%Window%KLB = aux(5) 
                NewProperty%Window%KUB = aux(6)
            else
                write(*,*)  
                write(*,*) 'Spatial window not specified for property: ',           &
                            trim(NewProperty%ID%Name)
                write(*,*) 'Assumed full domain...'

                NewProperty%Window%ILB = WorkILB 
                if (NewProperty%TypeZUV == TypeV_) then
                    NewProperty%Window%IUB = WorkIUB + 1
                else
                    NewProperty%Window%IUB = WorkIUB
                endif
                NewProperty%Window%JLB = WorkJLB 
                if (NewProperty%TypeZUV == TypeU_) then
                    NewProperty%Window%JUB = WorkJUB + 1
                else
                    NewProperty%Window%JUB = WorkJUB
                endif
                NewProperty%Window%KLB = WorkJLB 
                NewProperty%Window%KUB = WorkJUB
            endif

            !if (NewProperty%TypeZUV == TypeU_) then

                !if (Me%CyclicBoundary%ON .and. Me%CyclicBoundary%Direction ==       &
                !    DirectionX_ .and. NewProperty%Window%JLB == WorkJLB .and.       &
                !    NewProperty%Window%JUB == WorkJUB + 1) then

                !    NewProperty%CyclicBoundary = .true.
                !endif

            !elseif (NewProperty%TypeZUV == TypeV_) then

                !if (Me%CyclicBoundary%ON .and. Me%CyclicBoundary%Direction ==       &
                !    DirectionY_ .and. NewProperty%Window%ILB == WorkILB .and.       &
                !    NewProperty%Window%IUB == WorkIUB + 1) then

                !    NewProperty%CyclicBoundary = .true.
                !endif
            !endif

        endif

        deallocate (aux)

        !Allocate property field
        if (NewProperty%Dim == Dim_2D) then

            allocate(NewProperty%Field2D(NewProperty%Window2D%ILB:                  &
                                         NewProperty%Window2D%IUB,                  &
                                         NewProperty%Window2D%JLB:                  &
                                         NewProperty%Window2D%JUB))

        else !Dim_3D

            if ((NewProperty%ID%IDNumber == WaterFluxX_) .or.                       &
                (NewProperty%ID%IDNumber == WaterFluxY_)) then

                allocate(NewProperty%FieldR8(NewProperty%Window%ILB:                &
                                             NewProperty%Window%IUB,                &
                                             NewProperty%Window%JLB:                &
                                             NewProperty%Window%JUB,                &
                                             NewProperty%Window%KLB:                &
                                             NewProperty%Window%KUB))

                if (NewProperty%TypeZUV /= TypeZ_ .and. (Me%Output%TimeSerieON .or. &
                    Me%Output%HDF5ON)) then

                    allocate(NewProperty%CenteredFieldR8(NewProperty%Window%ILB:    &
                                                         NewProperty%Window%IUB,    &
                                                         NewProperty%Window%JLB:    &
                                                         NewProperty%Window%JUB,    &
                                                         NewProperty%Window%KLB:    &
                                                         NewProperty%Window%KUB))
                endif
            else

                allocate(NewProperty%Field(NewProperty%Window%ILB:                  &
                                           NewProperty%Window%IUB,                  &
                                           NewProperty%Window%JLB:                  &
                                           NewProperty%Window%JUB,                  &
                                           NewProperty%Window%KLB:                  &
                                           NewProperty%Window%KUB))

                if (NewProperty%TypeZUV /= TypeZ_ .and. (Me%Output%TimeSerieON .or. &
                    Me%Output%HDF5ON)) then

                    allocate(NewProperty%CenteredField(NewProperty%Window%ILB:      &
                                                       NewProperty%Window%IUB,      &
                                                       NewProperty%Window%JLB:      &
                                                       NewProperty%Window%JUB,      &
                                                       NewProperty%Window%KLB:      &
                                                       NewProperty%Window%KUB))
                endif
            endif      
        endif

        !Find if there are measurements for this property
        call GetData(NewProperty%Measure,                                           &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType = FromBlock,                                        &
                     keyword    = 'MEASURE',                                        &
                     Default    = .false.,                                          &
                     ClientModule = 'ModuleSequentialAssimilation',                 &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR09'

        !Find if property is from Hydrodynamic or WaterProperties modules
        if (Check_Hydrodynamic_Property(NewProperty%ID%IDNumber)) then
            
            if (Me%HydroSeqAssim) then
                NewProperty%ModuleType = 1
            else
                write (*,*)  
                write (*,*)'State hydrodynamic property but sequential assimilation not'
                write (*,*)'commanded in the ModuleHydrodynamic:'
                write (*,*) trim(NewProperty%ID%Name)
                stop 'ConstructProperty - ModuleSequentialAssimilation - ERR10'
            endif

        else if (Check_Water_Property(NewProperty%ID%IDNumber)) then

            if (Me%WaterSeqAssim) then

                !check if property is specified in WaterProperties input file (in ID list)
                if (Me%FullWaterPropNumber > 0 .and.                                &
                    PropertyInWaterPropertiesList(NewProperty%ID%IDNumber)) then
                    !property is specified in ModuleWaterProperties
                    NewProperty%ModuleType = 2
                else
                    write (*,*)  
                    write (*,*)'State water property absent from ModuleWaterProperties:'
                    write (*,*) trim(NewProperty%ID%Name)
                    stop 'ConstructProperty - ModuleSequentialAssimilation - ERR11'
                endif                                        
            else
                write (*,*)  
                write (*,*)'State water property but sequential assimilation not'
                write (*,*)'commanded in ModuleWaterProperties:'
                write (*,*) trim(NewProperty%ID%Name)
                stop 'ConstructProperty - ModuleSequentialAssimilation - ERR11'
            endif
        else
            write (*,*)  
            write (*,*)'State property is not from Hydrodynamic or WaterProperties modules:'
            write (*,*) trim(NewProperty%ID%Name)
            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR12'
        end if

        if (NewProperty%Measure) then

            !Find measurement blocks
do1 :       do
                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,           &
                                           begin_measure, end_measure, BlockFound,  &
                                           STAT = STAT_CALL)

cd1 :           if (STAT_CALL .EQ. SUCCESS_ ) then    
cd2 :               if (BlockFound) then                                                  

                        ! Construct a new measurement
                        call ConstructMeasurement(NewMeasure, NewProperty%Dim)

                        NewMeasure%ID%Name = NewProperty%ID%Name

                        if (NewProperty%Dim == Dim_2D) then

                            ! Construct measurement state location
                            call ConstructMeasureStateLoc(NewMeasure%Localization,  &
                                                          NewProperty%Dim,          &
                                                          NewProperty%TypeZUV,      &
                                                          PropertyWindow2D =        &
                                                          NewProperty%Window2D,     &
                                                          MeasureStateLocation =    &
                                                          NewMeasure%StateLocation) 

                            if (NewMeasure%StateLocation == FillValueInt) then
                                write(*,*)
                                write(*,*) 'Property: '//trim(NewProperty%ID%Name)  
                                write(*,*) 'Measure localization outside model domain:'
                                write(*,*) 'Localization I: ', NewMeasure%Localization%I
                                write(*,*) 'Localization J: ', NewMeasure%Localization%J
                            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR13'
                            endif

                        else !Dim_3D
  
                            ! Construct measurement state location
                            call ConstructMeasureStateLoc(NewMeasure%Localization,  &
                                                          NewProperty%Dim,          &
                                                          NewProperty%TypeZUV,      &
                                                          PropertyWindow =          &
                                                          NewProperty%Window,       &
                                                          MeasureStateLocation =    &
                                                          NewMeasure%StateLocation) 

                            if (NewMeasure%StateLocation == FillValueInt) then
                                write(*,*)
                                write(*,*) 'Property: '//trim(NewProperty%ID%Name)  
                                write(*,*) 'Measure localization outside model domain:'
                                write(*,*) 'Localization I: ', NewMeasure%Localization%I
                                write(*,*) 'Localization J: ', NewMeasure%Localization%J
                                write(*,*) 'Localization K: ', NewMeasure%Localization%K
                            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR14'
                            endif

                        endif

                        ! Add new measurement to list 
                        call AddMeasurement(NewMeasure)

                        NewProperty%MeasuresNumber = NewProperty%MeasuresNumber + 1
                        
                    else
                        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)                                 &
                            stop 'ConstructProperty - ModuleSequentialAssimilation - ERR15'
                        exit do1    !No more blocks
                    end if cd2

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                    write(*,*)  
                    write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    if(STAT_CALL .ne. SUCCESS_)                                     &
                        stop 'ConstructProperty - ModuleSequentialAssimilation - ERR16'
                end if cd1
            end do do1
        endif

        NewProperty%FirstStatePosition = Me%StateVarNumber + 1

        !Actualize the total model state variables number
        if (NewProperty%Dim == Dim_2D) then
            do j = NewProperty%Window2D%JLB, NewProperty%Window2D%JUB
                do i = NewProperty%Window2D%ILB, NewProperty%Window2D%IUB
                    if (Me%External_Var%WaterPoints(i,j,Me%WorkSize%KUB) == 1) then
                        !found state variable!
                        
                        Me%StateVarNumber = Me%StateVarNumber + 1
                    endif
                enddo
            enddo

        else !Dim_3D

            select case(NewProperty%TypeZUV)
                case(TypeZ_)
                    PropertyMap => Me%External_Var%WaterPoints

                case(TypeU_)
                    PropertyMap => Me%WaterFacesU

                case(TypeV_)
                    PropertyMap => Me%WaterFacesV
            end select

            do k = NewProperty%Window%KLB, NewProperty%Window%KUB
                do j = NewProperty%Window%JLB, NewProperty%Window%JUB
                    do i = NewProperty%Window%ILB, NewProperty%Window%IUB
                        if (PropertyMap(i,j,k) == 1) then
                            !found state variable!
                            
                            Me%StateVarNumber = Me%StateVarNumber + 1
                        endif
                    enddo
                enddo
            enddo
        endif

    end subroutine ConstructProperty

    !--------------------------------------------------------------------------

    logical function PropertyInWaterPropertiesList(Property)

        !Arguments-------------------------------------------------------------
        integer, intent (IN)                            :: Property

        !Local-----------------------------------------------------------------
        integer                                         :: nwp

        !----------------------------------------------------------------------

        !Find if property is in ModuleWaterProperties
        PropertyInWaterPropertiesList = .false.

        do nwp =1, Me%FullWaterPropNumber
            if (Property == Me%PropertiesIDArray(nwp)) then
                PropertyInWaterPropertiesList = .true.
                exit
            endif
        enddo       

    end function PropertyInWaterPropertiesList

    !--------------------------------------------------------------------------

    subroutine ConstructMeasurement(NewMeasure, PropertyDimension)

        !Arguments---------------------------------------------------------------
        type(T_Measure),   pointer                      :: NewMeasure
        integer                                         :: PropertyDimension

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: iflag
        integer                                         :: STAT_CALL
        character(len=StringLength)                     :: AuxString
        integer                                         :: AuxValueInt
        real                                            :: AuxValueReal

        !------------------------------------------------------------------------

        !Allocates new measurement
        allocate (NewMeasure, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR01'

        !Get keywords from measurement block
        !Get file type
        call GetData(AuxString, Me%ObjEnterData,  iflag,                        &
                     SearchType     = FromBlockInBlock,                         &
                     keyword        = 'FILE_IN_TIME',                           &
                     default        = "Timeserie",                              &
                     ClientModule   = 'ModuleSequentialAssimilation',           &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR02'

        select case (trim(adjustl(AuxString)))
            case ("Hdf",        "HDF",          "hdf")
                NewMeasure%FileType    = HDF
                write(*,*)
                write(*,*)'Not implemented option for keyword FILE_IN_TIME'
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR03'

            case ("Timeserie", "TIMESERIE", "timeserie", "TimeSerie")
                NewMeasure%FileType    = TimeSerie

            case ("Profile_Timeserie", "PROFILE_TIMESERIE",                     &
                  "profile_timeserie", "Profile_TimeSerie")
                NewMeasure%FileType    = ProfileTimeSerie
                write(*,*)
                write(*,*)'Not implemented option for keyword FILE_IN_TIME'
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR04'

            case default
                write(*,*)
                write(*,*)'Invalid option for keyword FILE_IN_TIME'
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR05'
        end select

        !Get measurement time interval
        call GetData(NewMeasure%DT, Me%ObjEnterData, iflag,                     &
                     SearchType = FromBlockInBlock,                             &
                     keyword    = 'DT',                                         &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR06'
        
        if (NewMeasure%DT > Me%DT) then
            write(*,*)
            write(*,*)'Measurement DT larger that specified assimilation DT.'
            stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR07'
        else
            !Check if assimilation interval is multiple of measurement DT
            AuxValueInt = Me%DT/NewMeasure%DT
            AuxValueReal = Me%DT/NewMeasure%DT
            if ((AuxValueReal/AuxValueInt) /= 1.) then
                write(*,*)
                write(*,*)'Specified assimilation DT is not multiple of measurement DT.'
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR08'
            endif
        endif

        if (Me%MeasureErrorVar) then
            !Get measurement error variance 
            call GetData(NewMeasure%Variance, Me%ObjEnterData, iflag,           &
                         SearchType = FromBlockInBlock,                         &
                         keyword    = 'VARIANCE',                               &
                         ClientModule = 'ModuleSequentialAssimilation',         &
                         STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR09'

            if (iflag == 0) then
        
                if (Me%ObjectiveAnalysis .and. Me%MeasureErrorVar) then
                    Me%MeasureErrorVar = .false.
                    write(*,*)
                    write(*,*)                                                  &
                        'At least one measurement without error variance specification.'
                    write(*,*)                                                  &
                        'Objective Analysis is ON: error variance is not considered.'
                else
                    write(*,*)
                    write(*,*)                                                  &
                        'Assimilation measurement without error variance specification.'
                    write(*,*)'Objective Analysis is OFF: specify error variance.'
                    stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR10'
                endif
            endif
        endif

        !Get tolerance in measurement time 
        call GetData(NewMeasure%Tolerance, Me%ObjEnterData, iflag,              &
                     SearchType = FromBlockInBlock,                             &
                     keyword    = 'TIME_TOLERANCE',                             &
                     default    = 0.0,                                          &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR11'

        !Get measurement file name
        call GetData(NewMeasure%FileName,                                       &
                     Me%ObjEnterData , iflag,                                   &
                     SearchType   = FromBlockInBlock,                           &
                     keyword      = 'FILENAME',                                 &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR12'

        if (NewMeasure%FileType == TimeSerie) then

            !Get time series column
            call GetData(NewMeasure%TimeSerie%DataColumn,                       &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromBlockInBlock,                       &
                         keyword      = 'DATA_COLUMN',                          &
                         ClientModule = 'ModuleSequentialAssimilation',         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR13'

            !Get time series location
            !Searches for the Localization I
            call GetData(NewMeasure%Localization%I,                             &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromBlockInBlock,                       &
                         keyword      ='LOCALIZATION_I',                        &
                         ClientModule ='ModuleSequentialAssimilation',          &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .or. iflag /= 1)                        &
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR14'

            !Searches for the Localization J
            call GetData(NewMeasure%Localization%J,                             &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromBlockInBlock,                       &
                         keyword      ='LOCALIZATION_J',                        &
                         ClientModule ='ModuleSequentialAssimilation',          &
                         STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .or. iflag /= 1)                        &
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR15'

            if (PropertyDimension == Dim_3D) then
                !Searches for the Localization K
                call GetData(NewMeasure%Localization%K,                         &
                             Me%ObjEnterData, iflag,                            &
                             SearchType   = FromBlockInBlock,                   &
                             keyword      ='LOCALIZATION_K',                    &
                             ClientModule ='ModuleSequentialAssimilation',      &
                             STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_ .or. iflag /= 1)                    &
                    stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR16'
            endif

            call StartTimeSerieInput(NewMeasure%TimeSerie%ID,                   &
                                     NewMeasure%FileName,                       &
                                     Me%ObjTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructMeasurement - ModuleSequentialAssimilation - ERR17'

        elseif (NewMeasure%FileType == HDF) then

            !...

        endif

    end subroutine ConstructMeasurement

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

    subroutine AddMeasurement(NewMeasure)

        !Arguments-------------------------------------------------------------
        type(T_Measure),   pointer          :: NewMeasure

        !----------------------------------------------------------------------

        ! Add to the list a new property
        if (.not.associated(Me%FirstMeasure)) then
            Me%MeasuresNumber = 1
            Me%FirstMeasure                 => NewMeasure
            Me%LastMeasure                  => NewMeasure
        else
            NewMeasure%Prev                 => Me%LastMeasure
            Me%LastMeasure%Next             => NewMeasure
            Me%LastMeasure                  => NewMeasure
            Me%MeasuresNumber               = Me%MeasuresNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine AddMeasurement 

    !--------------------------------------------------------------------------

    subroutine ConstructLog

        !Local-----------------------------------------------------------------
        type(T_Property),   pointer         :: CurrentProperty

        !----------------------------------------------------------------------

        write(*, *)"------------------ SEQUENTIAL ASSIMILATION -----------------"
        write(*, *)
        write(*, *)"Num of Properties      : ", Me%PropertiesNumber
        write(*, *)
        write(*, *)"Num of Measurements    : ", Me%MeasuresNumber
        write(*, *)

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))

            write(*, *)"Property               : ", trim(CurrentProperty%ID%Name)
            write(*, *)"---Num of Measurements : ", CurrentProperty%MeasuresNumber
            write(*, *)

            CurrentProperty=>CurrentProperty%Next
        enddo

    end subroutine ConstructLog

    !----------------------------------------------------------------------

    subroutine ConstructSeqAssimVariables

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        type(T_Measure),   pointer                      :: ObjMeasure
        integer                                         :: MeasureCount

        !------------------------------------------------------------------------

        Me%PreviousMeasuresNumber = Me%MeasuresNumber

        !Allocate data assimilation variables
        call AllocateVariables 

        !Obtain initial state covariance structure 
        call ConstructInitialStateCov

        ObjMeasure => Me%FirstMeasure
        MeasureCount = 0

        do while (associated (ObjMeasure))
            MeasureCount = MeasureCount + 1

            if (Me%MeasureErrorVar) then
                !Construction of measurement error covariance matrix
                Me%MeasuresInvCov(MeasureCount, MeasureCount) = 1/ObjMeasure%Variance
            endif

            !Construction of observation operator
            Me%ObservOperator(MeasureCount, ObjMeasure%StateLocation) = 1.0 

            ObjMeasure => ObjMeasure%Next
        enddo

    end subroutine ConstructSeqAssimVariables

    !----------------------------------------------------------------------

    Subroutine AllocateVariables

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: neof, nwp
        integer                                         :: STAT_CALL
        type (T_StateProp), pointer                     :: NewWaterProperty
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        !integer, dimension(:), pointer                  :: PropertiesIDArray

        !------------------------------------------------------------------------

        !State
        allocate (Me%State (1:Me%StateVarNumber,1:1))
        !State OpenPoints
        allocate (Me%OpenPointState (1:Me%StateVarNumber,1:1))
        !Measured state
        allocate (Me%MeasuredState (1:Me%MeasuresNumber,1:1))
        !Observation operator
        allocate (Me%ObservOperator (1:Me%MeasuresNumber, 1:Me%StateVarNumber))
        !Kalman gain
        !allocate (Me%KalmanGain (1:Me%StateVarNumber, 1:Me%MeasuresNumber))

        Me%State            (:,:)   = FillValueReal
        Me%OpenPointState   (:,:)   = 0.
        Me%MeasuredState    (:,:)   = FillValueReal
        Me%ObservOperator   (:,:)   = 0.
        !Me%KalmanGain     (:,:)   = 0.

        if (Me%MeasureErrorVar) then
            !Measurement error covariance matrix
            allocate (Me%MeasuresInvCov (1:Me%MeasuresNumber, 1:Me%MeasuresNumber))
            Me%MeasuresInvCov   (:,:)   = 0.
        endif

        if (Me%Method == SEEK) then
            !allocate error reduced space variables    
            allocate (Me%LMatrix (1:Me%StateVarNumber, 1:Me%StateCovRank))
            allocate (Me%InvUMatrix (1:Me%StateCovRank, 1:Me%StateCovRank))

            Me%LMatrix      (:,:)   = FillValueReal
            Me%InvUMatrix   (:,:)   = 0.

            if (Me%StateCovEvolution) then
                !get logical for submodel in ModuleHydrodynamic
                call GetHydroNeedsFather(Me%ObjHydrodynamic, Me%HydroSubModel, STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'AllocateVariables - ModuleSequentialAssimilation - ERR01'

                !allocate full state matrixes
                allocate (Me%DFS (1:Me%StateCovRank)) !One full state per each EOF

                ILB = Me%Size%ILB 
                IUB = Me%Size%IUB 

                JLB = Me%Size%JLB 
                JUB = Me%Size%JUB 

                KLB = Me%Size%KLB 
                KUB = Me%Size%KUB 

                !ModuleHydrodynamic
                allocate (Me%FS%WaterLevel      (ILB:IUB, JLB:JUB))
                allocate (Me%FS%VelocityU       (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%VelocityV       (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%VelocityUOld    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%VelocityVOld    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%VelocityAcross  (ILB:IUB, JLB:JUB, KLB:KUB)) !or cartesian?
                allocate (Me%FS%WaterFluxX      (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%WaterFluxY      (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%WaterFluxZ      (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%ChezyVelUV      (ILB:IUB, JLB:JUB)) 

                Me%FS%WaterLevel        (:,:)   = FillValueReal
                Me%FS%VelocityU         (:,:,:) = FillValueReal
                Me%FS%VelocityV         (:,:,:) = FillValueReal
                Me%FS%VelocityUOld      (:,:,:) = FillValueReal
                Me%FS%VelocityVOld      (:,:,:) = FillValueReal
                Me%FS%VelocityAcross    (:,:,:) = FillValueReal !or cartesian?
                Me%FS%WaterFluxX        (:,:,:) = dble(FillValueReal)
                Me%FS%WaterFluxY        (:,:,:) = dble(FillValueReal)
                Me%FS%WaterFluxZ        (:,:,:) = dble(FillValueReal)
                Me%FS%ChezyVelUV        (:,:)   = 0.

                !Do next only if submodel
                if (Me%HydroSubModel) then
                    allocate (Me%FS%SubModelqX      (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%FS%SubModelqY      (ILB:IUB, JLB:JUB, KLB:KUB))

                    Me%FS%SubModelqX    (:,:,:) = FillValueReal
                    Me%FS%SubModelqY    (:,:,:) = FillValueReal
                endif

                !Point to hydrodynamic state module memory space
                call PointToHydroState(Me%ObjHydrodynamic, STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'AllocateVariables - ModuleSequentialAssimilation - ERR02'

                !ModuleWaterProperties
                allocate (Me%FS%Density         (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%SigmaDensity    (ILB:IUB, JLB:JUB, KLB:KUB))

                Me%FS%Density           (:,:,:) = FillValueReal
                Me%FS%SigmaDensity      (:,:,:) = FillValueReal

                !Point to density state module memory space
                call PointToDensity(Me%ObjWaterProperties, STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'AllocateVariables - ModuleSequentialAssimilation - ERR03'

                nullify(Me%FS%FirstWaterProperty)
                nullify(Me%FS%LastWaterProperty)

                !ModuleGeometry
                allocate (Me%FS%SZZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%DWZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%DUZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%DVZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%DZZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%AreaU        (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%AreaV        (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%ZCellCenter  (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%WaterColumnU (ILB:IUB, JLB:JUB))
                allocate (Me%FS%WaterColumnV (ILB:IUB, JLB:JUB))
                allocate (Me%FS%WaterColumnZ (ILB:IUB, JLB:JUB))
                allocate (Me%FS%VolumeZ      (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%VolumeU      (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%VolumeV      (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Me%FS%VolumeZOld   (ILB:IUB, JLB:JUB, KLB:KUB))

                Me%FS%SZZ               (:,:,:) = FillValueReal
                Me%FS%DWZ               (:,:,:) = FillValueReal
                Me%FS%DUZ               (:,:,:) = FillValueReal
                Me%FS%DVZ               (:,:,:) = FillValueReal
                Me%FS%DZZ               (:,:,:) = FillValueReal
                Me%FS%AreaU             (:,:,:) = FillValueReal
                Me%FS%AreaV             (:,:,:) = FillValueReal
                Me%FS%ZCellCenter       (:,:,:) = FillValueReal
                Me%FS%WaterColumnU      (:,:)   = FillValueReal
                Me%FS%WaterColumnV      (:,:)   = FillValueReal
                Me%FS%WaterColumnZ      (:,:)   = FillValueReal
                Me%FS%VolumeZ           (:,:,:) = FillValueReal
                Me%FS%VolumeU           (:,:,:) = FillValueReal
                Me%FS%VolumeV           (:,:,:) = FillValueReal
                Me%FS%VolumeZOld        (:,:,:) = FillValueReal

                !Point to geometry state module memory space
                call PointToGeometryState(Me%ObjGeometry, STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'AllocateVariables - ModuleSequentialAssimilation - ERR04'

                do neof = 1, Me%StateCovRank

                    !ModuleHydrodynamic
                    allocate (Me%DFS(neof)%WaterLevel       (ILB:IUB, JLB:JUB))
                    allocate (Me%DFS(neof)%VelocityU        (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%VelocityV        (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%VelocityUOld     (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%VelocityVOld     (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%VelocityAcross   (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%WaterFluxX       (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%WaterFluxY       (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%WaterFluxZ       (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%ChezyVelUV       (ILB:IUB, JLB:JUB)) 

                    Me%DFS(neof)%WaterLevel     (:,:)   = FillValueReal
                    Me%DFS(neof)%VelocityU      (:,:,:) = FillValueReal
                    Me%DFS(neof)%VelocityV      (:,:,:) = FillValueReal
                    Me%DFS(neof)%VelocityUOld   (:,:,:) = FillValueReal
                    Me%DFS(neof)%VelocityVOld   (:,:,:) = FillValueReal
                    Me%DFS(neof)%VelocityAcross (:,:,:) = FillValueReal
                    Me%DFS(neof)%WaterFluxX     (:,:,:) = dble(FillValueReal)
                    Me%DFS(neof)%WaterFluxY     (:,:,:) = dble(FillValueReal)
                    Me%DFS(neof)%WaterFluxZ     (:,:,:) = dble(FillValueReal)
                    Me%DFS(neof)%ChezyVelUV     (:,:)   = 0.

                    !Do next only if submodel
                    if (Me%HydroSubModel) then
                        allocate (Me%DFS(neof)%SubModelqX   (ILB:IUB, JLB:JUB, KLB:KUB))
                        allocate (Me%DFS(neof)%SubModelqY   (ILB:IUB, JLB:JUB, KLB:KUB))

                        Me%DFS(neof)%SubModelqX (:,:,:) = FillValueReal
                        Me%DFS(neof)%SubModelqY (:,:,:) = FillValueReal
                    endif

                    !ModuleWaterProperties
                    allocate (Me%DFS(neof)%Density          (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%SigmaDensity     (ILB:IUB, JLB:JUB, KLB:KUB))

                    Me%DFS(neof)%Density        (:,:,:) = FillValueReal
                    Me%DFS(neof)%SigmaDensity   (:,:,:) = FillValueReal

                    nullify(Me%DFS(neof)%FirstWaterProperty)
                    nullify(Me%DFS(neof)%LastWaterProperty)

                    !ModuleGeometry
                    allocate (Me%DFS(neof)%SZZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%DWZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%DUZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%DVZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%DZZ          (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%AreaU        (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%AreaV        (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%ZCellCenter  (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%WaterColumnU (ILB:IUB, JLB:JUB))
                    allocate (Me%DFS(neof)%WaterColumnV (ILB:IUB, JLB:JUB))
                    allocate (Me%DFS(neof)%WaterColumnZ (ILB:IUB, JLB:JUB))
                    allocate (Me%DFS(neof)%VolumeZ      (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%VolumeU      (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%VolumeV      (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate (Me%DFS(neof)%VolumeZOld   (ILB:IUB, JLB:JUB, KLB:KUB))

                    Me%DFS(neof)%SZZ            (:,:,:) = FillValueReal
                    Me%DFS(neof)%DWZ            (:,:,:) = FillValueReal
                    Me%DFS(neof)%DUZ            (:,:,:) = FillValueReal
                    Me%DFS(neof)%DVZ            (:,:,:) = FillValueReal
                    Me%DFS(neof)%DZZ            (:,:,:) = FillValueReal
                    Me%DFS(neof)%AreaU          (:,:,:) = FillValueReal
                    Me%DFS(neof)%AreaV          (:,:,:) = FillValueReal
                    Me%DFS(neof)%ZCellCenter    (:,:,:) = FillValueReal
                    Me%DFS(neof)%WaterColumnU   (:,:)   = FillValueReal
                    Me%DFS(neof)%WaterColumnV   (:,:)   = FillValueReal
                    Me%DFS(neof)%WaterColumnZ   (:,:)   = FillValueReal
                    Me%DFS(neof)%VolumeZ        (:,:,:) = FillValueReal
                    Me%DFS(neof)%VolumeU        (:,:,:) = FillValueReal
                    Me%DFS(neof)%VolumeV        (:,:,:) = FillValueReal
                    Me%DFS(neof)%VolumeZOld     (:,:,:) = FillValueReal

                enddo

                if (Me%FullWaterPropNumber > 0) then

                    !nullify(Me%FS%FirstWaterProperty)
                    !nullify(Me%FS%LastWaterProperty)

                    if (.not. Me%WaterSeqAssim) then 
                        !Construct water properties IDs array
                        call ConstructPropertiesIDArray(Me%ObjWaterProperties, STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'AllocateVariables - ModuleSequentialAssimilation - ERR05'
                    endif

                    !Obtain IDs of water properties in ModuleWaterProperties
                    call GetWaterPropertiesIDArray(Me%ObjWaterProperties,               &
                                                   Me%PropertiesIDArray, STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'AllocateVariables - ModuleSequentialAssimilation - ERR06'

                    !Construct water properties list and allocate water properties variables
                    do nwp = 1, Me%FullWaterPropNumber

                        ! Construct a New Property ID from ModuleWaterProperties
                        allocate (NewWaterProperty, STAT = STAT_CALL)            
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'AllocateVariables - ModuleSequentialAssimilation - ERR07'
            
                        NewWaterProperty%IDNumber = Me%PropertiesIDArray(nwp)

                        ! Add new Property to the State List 
                        call AddWaterProperty (NewWaterProperty,                        &
                                               Me%FS%FirstWaterProperty,                &
                                               Me%FS%LastWaterProperty)         

                        allocate (NewWaterProperty%Field (ILB:IUB, JLB:JUB, KLB:KUB))

                        NewWaterProperty%Field (:,:,:)  = FillValueReal

                        !Point to water property state module memory space
                        call PointToConcentration(Me%ObjWaterProperties,                &
                                                  NewWaterProperty%IDNumber,            &
                                                  STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'AllocateVariables - ModuleSequentialAssimilation - ERR08'
                    enddo

                    do neof = 1, Me%StateCovRank

                        !nullify(Me%DFS(neof)%FirstWaterProperty)
                        !nullify(Me%DFS(neof)%LastWaterProperty)

                        do nwp = 1, Me%FullWaterPropNumber

                            ! Construct a New Property ID from ModuleWaterProperties
                            allocate (NewWaterProperty, STAT = STAT_CALL)            
                            if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'AllocateVariables - ModuleSequentialAssimilation - ERR09'
            
                            NewWaterProperty%IDNumber = Me%PropertiesIDArray(nwp)

                            ! Add new Property to the State List 
                            call AddWaterProperty (NewWaterProperty,                    &
                                                   Me%DFS(neof)%FirstWaterProperty,     &
                                                   Me%DFS(neof)%LastWaterProperty)                
                
                            allocate (NewWaterProperty%Field (ILB:IUB, JLB:JUB, KLB:KUB))
            
                            NewWaterProperty%Field (:,:,:)  = FillValueReal
                        enddo
                    enddo

                    call UngetWaterProperties(Me%ObjWaterProperties,                    &
                                              Me%PropertiesIDArray, STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'AllocateVariables - ModuleSequentialAssimilation - ERR10'
                endif
            endif
        endif

      !----------------------------------------------------------------------

    end subroutine AllocateVariables   

    !----------------------------------------------------------------------

    Subroutine AddWaterProperty (NewWaterProperty, FirstWaterProperty, LastWaterProperty)

        !Arguments---------------------------------------------------------------
        type (T_StateProp), pointer                 :: NewWaterProperty
        type (T_StateProp), pointer                 :: FirstWaterProperty
        type (T_StateProp), pointer                 :: LastWaterProperty

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        ! Add to the list a new property
        if (.not.associated(FirstWaterProperty)) then
            FirstWaterProperty              => NewWaterProperty
            LastWaterProperty               => NewWaterProperty
        else
            NewWaterProperty%Prev           => LastWaterProperty
            LastWaterProperty%Next          => NewWaterProperty
            LastWaterProperty               => NewWaterProperty
        end if 

      !----------------------------------------------------------------------

    end subroutine AddWaterProperty

    !----------------------------------------------------------------------

    subroutine ConstructInitialStateCov

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------
        integer                                         :: HDF5_READ

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        character(len=StringLength)                     :: AuxEOFName 
        !character(len=StringLength)                     ::  AuxUName
        integer, dimension(7)                           :: Dimensions
        integer                                         :: Rank, NumberEOF, neof
        !real, dimension(:), pointer                     :: EigenValue
        real(8), dimension (:), pointer                 :: EOFField
        !real, dimension(:,:), pointer                   :: UMatrix

        !------------------------------------------------------------------------

        !Open initial state covariance HDF5
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        call ConstructHDF5 (Me%ObjCovHDF5, trim(Me%InitialStateCovFile),        &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR01'

        !(according to method)
        if (Me%Method == SEEK) then
            !Read general information
            
            !Check number of EOFs
            call GetHDF5GroupNumberOfItems (Me%ObjCovHDF5,                      &
                                            "/AssimilationData/EOF", NumberEOF, &
                                            STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR02'

            if (NumberEOF < Me%StateCovRank) then
                write(*,*)  
                write(*,*) 'Initial state covariance file does not contain specified EOFs.'
                stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR03'
            else
                
                !(if EOF number is correct it is assumed that eigenvalules number is also
                !correct, because of PreProcessor algorithm)

                allocate(EOFField(1:Me%StateVarNumber))
                !allocate(EigenValue(1:1))
                !allocate(UMatrix(1:Me%StateCovRank,1:Me%StateCovRank))
                !UMatrix(:,:) = 0.

                do neof = 1, Me%StateCovRank 
                    
                    !Get EOFs
                    EOFField(:) = FillValueReal

                    !Get field ID
                    call GetHDF5GroupID(Me%ObjCovHDF5,                          & 
                                        "/AssimilationData/EOF",                &
                                        neof, AuxEOFName, Rank = Rank,          &
                                        Dimensions = Dimensions,                &
                                        STAT = STAT_CALL)                                
                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR04'

                    !Check consistency with previously inputed data
                    if ((Rank /= 1) .or. (Dimensions(1) /= Me%StateVarNumber)) then
                        write(*,*)  
                        write(*,*) 'Initial EOFs do not have correct dimension.'
                    stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR05'
                    endif

                    !Read EOFs
                    call HDF5SetLimits (Me%ObjCovHDF5, 1, Me%StateVarNumber,    &
                                        STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                &
                stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR06'
                    
                    call HDF5ReadData(Me%ObjCovHDF5, "/AssimilationData/EOF",   &
                                      AuxEOFName,                               &
                                      Array1D = EOFField,                       &
                                      STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR07'

                    Me%LMatrix(:,neof) = EOFField(:)

                    !(eigenvalues matrix is assumed identity)
                    !!Get eigenvalues
                    !!Get field ID
                    !call GetHDF5GroupID(Me%ObjCovHDF5,                          & 
                    !                    "/AssimilationData/EigenValue",         &
                    !                    neof, AuxUName, Rank = Rank,            &
                    !                    Dimensions = Dimensions,                &
                    !                STAT = STAT_CALL)                                
                    !if (STAT_CALL .NE. SUCCESS_)                                &
                !stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR08'

                    !!Check consistency with previously input data
                    !if ((Rank /= 1) .and. (Dimensions(1) /= 1)) then
                    !    write(*,*)  
                    !    write(*,*) 'Eigen value field does not have correct dimension.'
                    !stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR09'
                    !endif

                    !call HDF5SetLimits (Me%ObjCovHDF5, 1, 1, STAT = STAT_CALL)
                    !if (STAT_CALL .NE. SUCCESS_)                                &
                    !stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR10'

                    !call HDF5ReadData(Me%ObjCovHDF5, "/AssimilationData/EigenValue", &
                    !                  AuxUName, Array1D = EigenValue,           &
                    !                  STAT = STAT_CALL)
                    !if (STAT_CALL /= SUCCESS_)                                  &
                    !stop 'ConstructInitialStateCov - ModuleSequentialAssimilation - ERR11'

                    !UMatrix(neof, neof) = EigenValue(1)
                    Me%InvUMatrix(neof, neof) = 1.
                enddo
            endif

            !!Invert U matrix
            !call InvSingularDiagMatrix2D(UMatrix, Me%StateCovRank, Me%InvUMatrix) 

            !deallocate(UMatrix)
            deallocate(EOFField)
            !deallocate(EigenValue)

        else
            !Other methods not yet implemented

        endif

    end subroutine ConstructInitialStateCov

    !----------------------------------------------------------------------

    subroutine ConstructMeasureStateLoc(MeasureLocalization, PropertyDim,           &
                                        PropertyType, PropertyWindow2D,             &
                                        PropertyWindow, MeasureStateLocation) 

        !Arguments---------------------------------------------------------------
        type (T_Localization)                           :: MeasureLocalization
        integer                                         :: PropertyDim
        integer                                         :: PropertyType
        type (T_Size2D), optional                       :: PropertyWindow2D
        type (T_Size3D), optional                       :: PropertyWindow
        integer                                         :: MeasureStateLocation

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: k, i, j
        integer                                         :: AuxStateLocation
        integer, dimension(:,:,:), pointer              :: PropertyMap  => null()

        !------------------------------------------------------------------------

        MeasureStateLocation = FillValueInt

        AuxStateLocation = Me%StateVarNumber !this comes from the previous property

        if (PropertyDim == Dim_2D) then

            if (((MeasureLocalization%J >= PropertyWindow2D%JLB) .and.          &
                (MeasureLocalization%J <= PropertyWindow2D%JUB)) .and.          &
                ((MeasureLocalization%I >= PropertyWindow2D%ILB) .and.          &
                (MeasureLocalization%I <= PropertyWindow2D%IUB))) then
                !Measurement within domain

                !Cycle property spatial window and count variables
do1:            do j = PropertyWindow2D%JLB, PropertyWindow2D%JUB
                    do i = PropertyWindow2D%ILB, PropertyWindow2D%IUB
                        if (Me%External_Var%WaterPoints(i,j,Me%WorkSize%KUB) == 1) then
                            
                            !every WaterPoint/Face is a state variable
                            AuxStateLocation = AuxStateLocation + 1

                            if ((MeasureLocalization%J == j) .and.              &
                                (MeasureLocalization%I == i)) then

                                !Measurement state location found
                                MeasureStateLocation = AuxStateLocation
                                exit do1

                            endif
                        endif
                    enddo
                enddo do1
            endif

       elseif (PropertyDim == Dim_3D) then

            if (((MeasureLocalization%K >= PropertyWindow%KLB) .and.            &
                (MeasureLocalization%K <= PropertyWindow%KUB)) .and.            &
                ((MeasureLocalization%J >= PropertyWindow%JLB) .and.            &
                (MeasureLocalization%J <= PropertyWindow%JUB)) .and.            &
                ((MeasureLocalization%I >= PropertyWindow%ILB) .and.            &
                (MeasureLocalization%I <= PropertyWindow%IUB))) then
                !Measurement within domain

                select case(PropertyType)
                    case(TypeZ_)
                        PropertyMap => Me%External_Var%WaterPoints

                    case(TypeU_)
                        PropertyMap => Me%WaterFacesU

                    case(TypeV_)
                        PropertyMap => Me%WaterFacesV
                end select

                !Cycle property spatial window and count variables
do2:            do k = PropertyWindow%KLB, PropertyWindow%KUB
                    do j = PropertyWindow%JLB, PropertyWindow%JUB
                        do i = PropertyWindow%ILB, PropertyWindow%IUB
                            if (PropertyMap(i,j,k) == 1) then

                                !every WaterPoint/Face is a state variable                   
                                AuxStateLocation = AuxStateLocation + 1

                                if ((MeasureLocalization%K == k) .and.          &
                                    (MeasureLocalization%J == j) .and.          &
                                    (MeasureLocalization%I == i)) then

                                    !Measurement state location found
                                    MeasureStateLocation = AuxStateLocation
                                    exit do2

                                endif
                            endif
                        enddo
                    enddo
                enddo do2
            endif
        
        else

            stop 'ConstructMeasureStateLoc - ModuleSequentialAssimilation - ERR01'

        endif

        !(this subroutine explains the state filling method: klb->kub, 
        ! jlb->jub, ilb->iub)

    end subroutine ConstructMeasureStateLoc

    !----------------------------------------------------------------------

    subroutine ObjectiveAnalysis 

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real                                            :: Year, Month, Day  
        real                                            :: Hour, Minute, Second

        !------------------------------------------------------------------------

        call ReadMeasures

        !Perform data assimilation only if
        if (Me%CurrentMeasuresNumber >= Me%StateCovRank) then

            if (Me%Method == SEEK) then

                call InitialUCalculationSEEK

            else
                !...
            endif

            call ConstructModelState

            if (Me%Method == SEEK) then

                call SEEKAnalysis

            else
                !...
            endif

            call CopyModelAnalysedState

            if (Me%OutPut%HDF5ON .or. Me%OutPut%TimeSerieON) call SeqAssimilationOutPut

        else

            write(*,*)
            write(*,*) 'Number of measurements, ',                              &
                       Me%CurrentMeasuresNumber, ' is insuficcient.'  
            write(*,*) 'Data assimilation is not performed for time :'

100         format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)

            call ExtractDate(Me%CurrentTime, Year, Month, Day, Hour,            &
                             Minute, Second)
            write(*,fmt=100)Year, Month, Day, Hour, Minute, Second
        endif               

    end subroutine ObjectiveAnalysis

    !--------------------------------------------------------------------------

    subroutine ConstructOutPutTime

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        type(T_Time)                                    :: BeginTime, EndTime
        integer                                         :: AuxInt
        real                                            :: ModelDT

        !------------------------------------------------------------------------

        if (Me%ObjectiveAnalysis) then
            BeginTime = Me%CurrentTime
        else
            BeginTime = Me%CurrentTime + Me%DT
        endif

        AuxInt = (Me%EndTime - Me%CurrentTime)/Me%DT
        if (Me%EndTime > (Me%CurrentTime + AuxInt*Me%DT)) then
            EndTime = Me%CurrentTime + AuxInt*Me%DT
        else
            EndTime = Me%EndTime
        endif

        call GetOutPutTime(Me%ObjEnterData,                                     &
                           CurrentTime   = BeginTime,                           &
                           EndTime       = EndTime,                             &
                           keyword       = 'OUTPUT_TIME',                       &
                           SearchType    = FromFile,                            &
                           OutPutsTime   = Me%OutPut%OutTime,                   &
                           OutPutsOn     = Me%OutPut%HDF5ON,                    &
                           STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructOutPutTime - ModuleSequentialAssimilation - ERR01'

        if (Me%StateCovEvolution .and. (Me%Method == SEEK)) then

            call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructOutPutTime - ModuleSequentialAssimilation - ERR02' 

            BeginTime = Me%CurrentTime

            AuxInt = (Me%EndTime - Me%CurrentTime)/ModelDT
            if (Me%EndTime > (Me%CurrentTime + AuxInt*ModelDT)) then
                EndTime = Me%CurrentTime + AuxInt*ModelDT
            else
                EndTime = Me%EndTime
            endif

            call GetOutPutTime(Me%ObjEnterData,                                 &
                               CurrentTime   = BeginTime,                       &
                               EndTime       = EndTime,                         &
                               keyword       = 'EOF_OUTPUT_TIME',               &
                               SearchType    = FromFile,                        &
                               OutPutsTime   = Me%EOF_OutPut%OutTime,           &
                               OutPutsOn     = Me%EOF_OutPut%HDF5ON,            &
                               STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructOutPutTime - ModuleSequentialAssimilation - ERR03'

        endif

        if (Me%OutPut%HDF5ON) then 

            Me%OutPut%NextOutPut = 1

        endif

        if (Me%EOF_OutPut%HDF5ON) then 

            Me%EOF_OutPut%NextOutPut = 1

        endif

    end subroutine ConstructOutPutTime 

    !--------------------------------------------------------------------------

    subroutine OpenHDF5OutPutFile

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: HDF5_CREATE
        real, pointer, dimension(:, :)                  :: Bathymetry
        integer                                         :: WorkILB, WorkIUB
        integer                                         :: WorkJLB, WorkJUB
        integer                                         :: WorkKLB, WorkKUB

        !------------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        !Gets a pointer to Bathymetry
        call GetGridData      (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR01'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5,                                    &
                                 trim(Me%Files%OutPutFields)//"5",              &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,            &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR03'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Bathymetry,                             &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR04'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",        &
                              Array3D = Me%External_Var%WaterPoints,                                   &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR05'

        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB+1, WorkJLB,          &
                              WorkJUB+1, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR06'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR07'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR08'

        !call output of EOF
        if (Me%EOF_OutPut%HDF5ON) call EOFOutPut

        !Ungets the Bathymetry
        call UngetGridData (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OpenHDF5OutPutFile - ModuleSequentialAssimilation - ERR09'

    end subroutine OpenHDF5OutPutFile 

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON
        integer                                             :: STAT_CALL
        integer                                             :: iflag, dn, Id, Jd 
        integer                                             :: TimeSerieNumber, i, j
        character(len=PathLength)                           :: TimeSerieLocationFile
        type(T_Property), pointer                           :: ObjProperty
        integer                                             :: Count
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Begin-----------------------------------------------------------------
       
        !Allocates PropertyList
        allocate(PropertyList(Me%PropertiesNumber), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR01' 

        !Fills up PropertyList
        ObjProperty => Me%FirstProperty
        Count = 0

        do while (associated (ObjProperty))
            Count = Count + 1
            PropertyList(Count) = trim(GetPropertyName (ObjProperty%ID%IDNumber))

            ObjProperty => ObjProperty%Next
        enddo

        do i=1,Me%PropertiesNumber
            do j=1,len_trim(PropertyList(i))
                if (PropertyList(i)(j:j)==' ') PropertyList(i)(j:j)='_'
            enddo
        enddo

        call GetData(TimeSerieLocationFile,                                     &
                     Me%ObjEnterData,iflag,                                     &
                     SearchType   = FromFile,                                   &
                     keyword      = 'TIME_SERIE_LOCATION',                      &
                     ClientModule = 'ModuleSequentialAssimilation',             &
                     Default      = Me%Files%ConstructData,                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR02' 

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                        &
                            TimeSerieLocationFile,                              &
                            PropertyList, "srsa",                               &
                            WaterPoints3D = Me%External_Var%WaterPoints,        &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR03' 

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR04' 

        !Corrects if necessary the cell based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR05'  

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                      &  
                                      CoordX   = CoordX,                        &
                                      CoordY   = CoordY,                        & 
                                      CoordON  = CoordON,                       &
                                      STAT     = STAT_CALL)
            if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd,   &
                                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR06'

                if (Id < 0 .or. Jd < 0)                                         &
                    stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR07'

                if (Me%External_Var%WaterPoints(Id,JD,Me%WorkSize%KUB) == WaterPoint) &
                    then

                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd,    &
                                                STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR08'
                else

                    stop 'ConstructTimeSerie - ModuleSequentialAssimilation - ERR09'

                endif
            endif

        enddo
        
    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------

    subroutine ConstructErrorTimeSerie

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: nProperties, imeas
        integer                                             :: ObjTimeSerie
        character(len=StringLength)                         :: MeasureName
        type (T_ErrorTimeSerie), pointer                    :: PrevPointer
        type (T_ErrorTimeSerie), pointer                    :: Aux_Pointer
        integer                                             :: STAT_CALL
        character(StringLength), pointer, dimension(:)      :: PropertyList

        !Begin-----------------------------------------------------------------
       
        !First each time series has two properties: forecast error, analysis error
        nProperties = 2

        !Allocates Time Serie PropertyList
        allocate(PropertyList(nProperties), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructErrorTimeSerie - ModuleSequentialAssimilation - ERR01'

        !Fills up PropertyList
        PropertyList(1) = "forecast_error"
        PropertyList(2) = "analysis_error"

        !Nullifies first 
        nullify(Me%FirstErrorTimeSerie)

        !For each measure allocate a time serie
        do imeas = 1, Me%MeasuresNumber

            write(MeasureName,'(i3)') imeas
          
            !Starts a new Time Serie
            ObjTimeSerie = 0
            call StartTimeSerie(ObjTimeSerie,                                   &
                                Me%ObjTime,                                     &
                                Me%Files%ConstructData,                         &
                                PropertyList, "srsae",                          &
                                ResultFileName = trim(MeasureName),             &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructErrorTimeSerie - ModuleSequentialAssimilation - ERR02'

            !Insert New Time Serie to List
            if (.not. associated(Me%FirstErrorTimeSerie)) then
                allocate(Me%FirstErrorTimeSerie)
                Me%FirstErrorTimeSerie%TimeSerie = ObjTimeSerie
                nullify(Me%FirstErrorTimeSerie%Next)
            else
                PrevPointer             => Me%FirstErrorTimeSerie
                Aux_Pointer             => Me%FirstErrorTimeSerie%Next
                do while (associated(Aux_Pointer))
                    PrevPointer         => Aux_Pointer
                    Aux_Pointer         => Aux_Pointer%Next
                enddo
                allocate(Aux_Pointer)
                Aux_Pointer%TimeSerie   =  ObjTimeSerie
                PrevPointer%Next        => Aux_Pointer
                nullify(Aux_Pointer%Next)
            endif
        enddo

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'ConstructErrorTimeSerie - ModuleSequentialAssimilation - ERR03'
        
    end subroutine ConstructErrorTimeSerie

    !----------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    !--------------------------------------------------------------------------

    subroutine GetSeqAssimilationTime(SequentialAssimilationID,                 &
                                      AssimilationTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: SequentialAssimilationID
        type(T_Time),      intent(OUT)              :: AssimilationTime
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SequentialAssimilationID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            AssimilationTime = Me%AssimilationTime

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                      &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetSeqAssimilationTime

    !--------------------------------------------------------------------------
    
    subroutine GetSeqAssimilationOptions(SequentialAssimilationID,              &
                                         StateCovEvolution, StateCovRank, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: SequentialAssimilationID
        logical,           intent(OUT)  :: StateCovEvolution
        integer,           intent(OUT)  :: StateCovRank
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SequentialAssimilationID, ready_)  
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                     &
            (ready_ == READ_LOCK_ERR_)) then

            StateCovEvolution = Me%StateCovEvolution

            StateCovRank = Me%StateCovRank

            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetSeqAssimilationOptions

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine SetModelInitialState(CycleNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: CycleNumber

        !Local-----------------------------------------------------------------
        !real, pointer, dimension(:, :)              :: ObjArray2D => null()
        !real, pointer, dimension(:, :,:)            :: ObjArray3D => null()

        !----------------------------------------------------------------------

        if (CycleNumber == 1) then

            !Get full model state from ModuleHydrodynamic, ModuleWaterProperties and 
            !ModuleGeometry to FS variables

            call CopyModelResultToFS(Me%FS)

        endif

        if (CycleNumber /= Me%StateCovRank + 1) then

            !Copy from FS variables to respective DFS variables (CycleNumber signals EOF)
            if (Me%PerturbeState .or. Me%StateCov_Evolution_Split) then 

                call CopyFullModelState(CycleNumber)

                !Disturb respective DFS variables (only assimilation state variables)
                call DisturbFullModelState(CycleNumber)

            endif

            !SET ModuleHydrodynamic, ModuleWaterProperties and ModuleGeometry 
            !respective DFS variables
            call SetDisturbedStateInModel(CycleNumber)
            !(SET => Model works on memory in current module)

        else

            !COPY FS variables to ModuleHydrodynamic and ModuleWaterProperties variables
            call CopyUndisturbedStateToModel
            !(COPY => Model works on memory in Hydrodynamic and WaterProperties modules)

            !State should only be perturbed again if new covariance is calculated
            Me%PerturbeState = .false.

        endif

    end subroutine SetModelInitialState

    !--------------------------------------------------------------------------

    subroutine CopyFullModelState(CycleNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: CycleNumber

        !Local-----------------------------------------------------------------
        type (T_StateProp), pointer                 :: ObjWaterPropertyFS
        type (T_StateProp), pointer                 :: ObjWaterPropertyDFS

        !----------------------------------------------------------------------

        !ModuleHydrodynamic variables
        !>>> Water Level
        Me%DFS(CycleNumber)%WaterLevel(:,:)         = Me%FS%WaterLevel(:,:)
        !>>> Velocity U
        Me%DFS(CycleNumber)%VelocityU(:,:,:)        = Me%FS%VelocityU(:,:,:)
        !>>> Velocity V
        Me%DFS(CycleNumber)%VelocityV(:,:,:)        = Me%FS%VelocityV(:,:,:)
        !>>> Velocity U Old
        Me%DFS(CycleNumber)%VelocityUOld(:,:,:)     = Me%FS%VelocityUOld(:,:,:)
        !>>> Velocity V Old
        Me%DFS(CycleNumber)%VelocityVOld(:,:,:)     = Me%FS%VelocityVOld(:,:,:)
        !>>> Velocity W
        Me%DFS(CycleNumber)%VelocityAcross(:,:,:)   = Me%FS%VelocityAcross(:,:,:)
        !>>> Water Flux X
        Me%DFS(CycleNumber)%WaterFluxX(:,:,:)       = Me%FS%WaterFluxX(:,:,:)
        !>>> Water Flux Y
        Me%DFS(CycleNumber)%WaterFluxY(:,:,:)       = Me%FS%WaterFluxY(:,:,:)
        !>>> Water Flux Z
        Me%DFS(CycleNumber)%WaterFluxZ(:,:,:)       = Me%FS%WaterFluxZ(:,:,:)
        !>>> ChezyVelUV
        Me%DFS(CycleNumber)%ChezyVelUV(:,:)         = Me%FS%ChezyVelUV(:,:)

        if (Me%HydroSubModel) then
            !>>> SubModel fluxes
            Me%DFS(CycleNumber)%SubModelqX(:,:,:)   = Me%FS%SubModelqX(:,:,:)
            Me%DFS(CycleNumber)%SubModelqY(:,:,:)   = Me%FS%SubModelqY(:,:,:)
        endif

        !ModuleWaterProperties variables
        !>>> Density
        Me%DFS(CycleNumber)%Density(:,:,:)          = Me%FS%Density(:,:,:)
        !>>> Sigma Density
        Me%DFS(CycleNumber)%SigmaDensity(:,:,:)     = Me%FS%SigmaDensity(:,:,:)

        ObjWaterPropertyFS => Me%FS%FirstWaterProperty
        ObjWaterPropertyDFS => Me%DFS(CycleNumber)%FirstWaterProperty

        do while (associated(ObjWaterPropertyFS)) !Same association number in both lists

            ObjWaterPropertyDFS%Field(:,:,:)        = ObjWaterPropertyFS%Field(:,:,:)

            ObjWaterPropertyFS  => ObjWaterPropertyFS%Next
            ObjWaterPropertyDFS => ObjWaterPropertyDFS%Next

        enddo

        !ModuleGeometry variables
        !>>> Distances
        Me%DFS(CycleNumber)%SZZ(:,:,:)              = Me%FS%SZZ(:,:,:)
        Me%DFS(CycleNumber)%DWZ(:,:,:)              = Me%FS%DWZ(:,:,:)
        Me%DFS(CycleNumber)%DUZ(:,:,:)              = Me%FS%DUZ(:,:,:)
        Me%DFS(CycleNumber)%DVZ(:,:,:)              = Me%FS%DVZ(:,:,:)
        Me%DFS(CycleNumber)%DZZ(:,:,:)              = Me%FS%DZZ(:,:,:)
        Me%DFS(CycleNumber)%ZCellCenter(:,:,:)      = Me%FS%ZCellCenter(:,:,:)
        !>>> Areas
        Me%DFS(CycleNumber)%AreaU(:,:,:)            = Me%FS%AreaU(:,:,:)
        Me%DFS(CycleNumber)%AreaV(:,:,:)            = Me%FS%AreaV(:,:,:)
        !>>> WaterColumn
        Me%DFS(CycleNumber)%WaterColumnZ(:,:)       = Me%FS%WaterColumnZ(:,:)
        Me%DFS(CycleNumber)%WaterColumnU(:,:)       = Me%FS%WaterColumnU(:,:)
        Me%DFS(CycleNumber)%WaterColumnV(:,:)       = Me%FS%WaterColumnV(:,:)
        !>>> Volume
        Me%DFS(CycleNumber)%VolumeZ(:,:,:)          = Me%FS%VolumeZ(:,:,:)
        Me%DFS(CycleNumber)%VolumeZOld(:,:,:)       = Me%FS%VolumeZOld(:,:,:)
        Me%DFS(CycleNumber)%VolumeU(:,:,:)          = Me%FS%VolumeU(:,:,:)
        Me%DFS(CycleNumber)%VolumeV(:,:,:)          = Me%FS%VolumeV(:,:,:)

    end subroutine CopyFullModelState

    !--------------------------------------------------------------------------

    subroutine DisturbFullModelState(CycleNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: CycleNumber

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real, dimension (:, :), pointer             :: AuxMatrixDFS2D
        real, dimension (:, :, :), pointer          :: AuxMatrixDFS
        real(8), dimension (:, :, :), pointer       :: AuxMatrixDFSR8
        integer                                     :: i,j,k, Count
        type(T_Property), pointer                   :: ObjProperty
        type (T_StateProp) , pointer                :: ObjDFSProperty
        integer, dimension(:,:,:), pointer          :: PropertyMap      => null()
        !integer, dimension(:,:,:), pointer          :: PropertyOpenMap  => null()

        !----------------------------------------------------------------------

        !(only ModuleHydrodynamic or ModuleWaterProperties variables are potentially
        !disturbed, because are the only possible state variables)

        !Get WaterPoints
        call GetWaterPoints3D(Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR01'

        !Get OpenPoints
        call GetOpenPoints3D(Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR02'

        ObjProperty => Me%FirstProperty

        Count = 0 !indicator of property position in EOF vector

        do while (associated (ObjProperty))

            !Select state property matrixes
            if (ObjProperty%ModuleType == 1) then !ModuleHydrodynamic property

case1 :         select case(ObjProperty%ID%IDNumber)
                    case(WaterLevel_    )

                        AuxMatrixDFS2D => Me%DFS(CycleNumber)%WaterLevel 

                    case(VelocityU_     )

                        AuxMatrixDFS => Me%DFS(CycleNumber)%VelocityU 

                    case(VelocityV_     ) 

                        AuxMatrixDFS => Me%DFS(CycleNumber)%VelocityV 

                    case(VelocityW_     )

                        AuxMatrixDFS => Me%DFS(CycleNumber)%VelocityAcross

                    case(WaterFluxX_    )

                        AuxMatrixDFSR8 => Me%DFS(CycleNumber)%WaterFluxX

                    case(WaterFluxY_    )

                        AuxMatrixDFSR8 => Me%DFS(CycleNumber)%WaterFluxY

                    case default

                        stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR03' 

                end select case1

            else !ModuleWaterProperties property

                !All water properties are assumed 3D!

                !Cycle water properties list and find the correct matrix to point to
                ObjDFSProperty => Me%DFS(CycleNumber)%FirstWaterProperty

                do while (associated (ObjDFSProperty))

                    if (ObjProperty%ID%IDNumber == ObjDFSProperty%IDNumber) then
                        !Property found in full state!

                        AuxMatrixDFS => ObjDFSProperty%Field
                    
                        exit !cycle proceeds only till properties found 

                    endif

                    ObjDFSProperty => ObjDFSProperty%Next
                enddo

                !Density variables are not suitable for data assimilation state

            endif

            if (.not. associated(AuxMatrixDFS) .and.                            &
                .not. associated(AuxMatrixDFS2D) .and.                          &
                .not. associated(AuxMatrixDFSR8))                               &
                stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR04' 

            !Cycle all property's variables
            if (ObjProperty%Dim == Dim_2D) then !Typically this is Water Level

                do j = ObjProperty%Window2D%JLB, ObjProperty%Window2D%JUB
                    do i = ObjProperty%Window2D%ILB, ObjProperty%Window2D%IUB
                        
                        if (Me%External_Var%WaterPoints(i,j,Me%WorkSize%KUB) == 1) then
                            !found state variable!

                            Count = Count + 1 !Count advances per state variable

                            !if (Me%External_Var%OpenPoints(i,j,Me%WorkSize%KUB) == 1) then

                                !Calculate new disturbed property value
                                AuxMatrixDFS2D(i,j) = AuxMatrixDFS2D(i,j) +     &
                                                      Me%StateCov_LinFactor*    &
                                                      Me%LMatrix(Count, CycleNumber)
                            !endif
                        endif
                    enddo
                enddo   

                nullify(AuxMatrixDFS2D)

            else !(Dim_3D)

                select case(ObjProperty%TypeZUV)
                    case(TypeZ_)
                        PropertyMap     => Me%External_Var%WaterPoints
                        !PropertyOpenMap => Me%External_Var%OpenPoints

                    case(TypeU_)
                        !if (.not. associated(Me%External_Var%WetFacesU)) then
                        !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU, &
                        !                     Me%External_Var%WetFacesV, STAT = STAT_CALL)
                        !    if (STAT_CALL /= SUCCESS_)                          &
                        !stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR05'
                        !endif
                        PropertyMap     => Me%WaterFacesU
                        !PropertyOpenMap => Me%External_Var%WetFacesU

                    case(TypeV_)
                        !if (.not. associated(Me%External_Var%WetFacesV)) then
                        !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU, &
                        !                     Me%External_Var%WetFacesV, STAT = STAT_CALL)
                        !    if (STAT_CALL /= SUCCESS_)                          &
                        !stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR06'
                        !endif
                        PropertyMap     => Me%WaterFacesV
                        !PropertyOpenMap => Me%External_Var%WetFacesV
                end select

                if ((ObjProperty%ID%IDNumber == WaterFluxX_) .or.               &
                    (ObjProperty%ID%IDNumber == WaterFluxY_)) then

                    do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                        do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                            do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB

                                if (PropertyMap(i,j,k) == 1) then
                                    !found state variable!

                                    Count = Count + 1 !Count advances per state variable
                            
                                    !if (PropertyOpenMap(i,j,k) == 1) then

                                        !Calculate new disturbed property value
                                        AuxMatrixDFSR8(i,j,k) =                 &
                                                        AuxMatrixDFSR8(i,j,k) + &
                                                        Me%StateCov_LinFactor*  &
                                                        Me%LMatrix(Count, CycleNumber)
                                    !endif
                                endif
                            enddo
                        enddo
                    enddo

                    nullify(AuxMatrixDFSR8)

                else

                    do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                        do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                            do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB

                                if (PropertyMap(i,j,k) == 1) then
                                    !found state variable!

                                    Count = Count + 1 !Count advances per state variable
                            
                                    !if (PropertyOpenMap(i,j,k) == 1) then

                                        !Calculate new disturbed property value
                                        AuxMatrixDFS(i,j,k) =                   &
                                                        AuxMatrixDFS(i,j,k) +   &
                                                        Me%StateCov_LinFactor*  &
                                                        Me%LMatrix(Count, CycleNumber)
                                    !endif
                                endif
                            enddo
                        enddo
                    enddo

                    nullify(AuxMatrixDFS)

                endif
            endif

            ObjProperty => ObjProperty%Next

        enddo

        call UnGetMap (Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR05'

        call UnGetMap (Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR06'

        !if (associated(Me%External_Var%WetFacesU)) then
        !    call UnGetMap (Me%ObjMap, Me%External_Var%WetFacesU, STAT = STAT_CALL)
        !    if (STAT_CALL /= SUCCESS_)                                          &
        !        stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR09'
        !endif

        !if (associated(Me%External_Var%WetFacesV)) then
        !    call UnGetMap (Me%ObjMap, Me%External_Var%WetFacesV, STAT = STAT_CALL)
        !    if (STAT_CALL /= SUCCESS_)                                          &
        !        stop 'DisturbFullModelState - ModuleSequentialAssimilation - ERR10'
        !endif

    end subroutine DisturbFullModelState

    !--------------------------------------------------------------------------

    subroutine SetDisturbedStateInModel(CycleNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: CycleNumber

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_StateProp), pointer                 :: ObjWaterProperty

        !----------------------------------------------------------------------

        !This is SetModelInitialState subroutine backwards (instead of Get is used Set)
        !(All full state variables must be set because state is possibly fully altered by
        !previous CycleNumber item)

        !ModuleHydrodynamic variables
        !>>> Water Level
        call SetWaterLevel(Me%ObjHydrodynamic, Me%DFS(CycleNumber)%WaterLevel,  &
                           STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR01'
        !>>> Velocity U
        call SetHorizontalVelocity(Me%ObjHydrodynamic,                          &
                                   Velocity_U = Me%DFS(CycleNumber)%VelocityU,  &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR02'
        !>>> Velocity V
        call SetHorizontalVelocity(Me%ObjHydrodynamic,                          &
                                   Velocity_V = Me%DFS(CycleNumber)%VelocityV,  &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR03'
        !>>> Velocity U Old
        call SetHorizontalVelocity(Me%ObjHydrodynamic,                          &
                                   OldVelocity_U = Me%DFS(CycleNumber)%VelocityUOld, &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR04'
        !>>> Velocity V Old
        call SetHorizontalVelocity(Me%ObjHydrodynamic,                          &
                                   OldVelocity_V = Me%DFS(CycleNumber)%VelocityVOld, &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR05'
        !>>> Velocity W
        call SetVerticalVelocity(Me%ObjHydrodynamic,                            &
                          Velocity_Across = Me%DFS(CycleNumber)%VelocityAcross, &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR06'
        !>>> Water Flux X
        call SetWaterFluxes(Me%ObjHydrodynamic,                                 &
                            WaterFluxX = Me%DFS(CycleNumber)%WaterFluxX,        &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR07'
        !>>> Water Flux Y
        call SetWaterFluxes(Me%ObjHydrodynamic,                                 &
                            WaterFluxY = Me%DFS(CycleNumber)%WaterFluxY,        &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR08'
        !>>> Water Flux Z
        call SetWaterFluxes(Me%ObjHydrodynamic,                                 &
                            WaterFluxZ = Me%DFS(CycleNumber)%WaterFluxZ,        &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR09'
        !>>> ChezyVelUV
        call SetChezyVelUV(Me%ObjHydrodynamic,                                  &
                           ChezyVelUV = Me%DFS(CycleNumber)%ChezyVelUV,         &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR10'

        if (Me%HydroSubModel) then
            !>>> SubModel fluxes
            call SetSubModelFluxes(Me%ObjHydrodynamic,                          &
                                   SubModelqX = Me%DFS(CycleNumber)%SubModelqX, &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR11'
            call SetSubModelFluxes(Me%ObjHydrodynamic,                          &
                                   SubModelqY = Me%DFS(CycleNumber)%SubModelqY, &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR12'
        endif

        !ModuleWaterProperties variables
        !>>> Density
        call SetDensity(Me%ObjWaterProperties, Me%DFS(CycleNumber)%Density, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR13'
        !>>> Sigma Density
        call SetSigma(Me%ObjWaterProperties, Me%DFS(CycleNumber)%SigmaDensity, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR14'
        !>>> Water properties
        ObjWaterProperty => Me%DFS(CycleNumber)%FirstWaterProperty

        do while (associated(ObjWaterProperty))

            call SetConcentration(Me%ObjWaterProperties, ObjWaterProperty%Field, &
                                  ObjWaterProperty%IDNumber, STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR15'

            ObjWaterProperty => ObjWaterProperty%Next

        enddo

        !ModuleGeometry variables
        !>>> Distances
        call SetGeometryDistances(Me%ObjGeometry, Me%DFS(CycleNumber)%SZZ,      &
                                  Me%DFS(CycleNumber)%DWZ,                      &
                                  Me%DFS(CycleNumber)%DUZ, Me%DFS(CycleNumber)%DVZ, &
                                  Me%DFS(CycleNumber)%DZZ,                      &
                                  Me%DFS(CycleNumber)%ZCellCenter, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR16'

        !>>> Areas
        call SetGeometryAreas(Me%ObjGeometry, Me%DFS(CycleNumber)%AreaU,        &
                              Me%DFS(CycleNumber)%AreaV, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR17'

        !>>> WaterColumn
        call SetGeometryWaterColumn(Me%ObjGeometry,                             &
                                    Me%DFS(CycleNumber)%WaterColumnZ,           &
                                    Me%DFS(CycleNumber)%WaterColumnU,           &
                                    Me%DFS(CycleNumber)%WaterColumnV, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR18'

        !>>> Volumes
        call SetGeometryVolumes(Me%ObjGeometry, Me%DFS(CycleNumber)%VolumeZ,    &
                                Me%DFS(CycleNumber)%VolumeZOld,                 &
                                Me%DFS(CycleNumber)%VolumeU,                    &
                                Me%DFS(CycleNumber)%VolumeV, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'SetDisturbedStateInModel - ModuleSequentialAssimilation - ERR19'

    end subroutine SetDisturbedStateInModel

    !--------------------------------------------------------------------------

    subroutine CopyUndisturbedStateToModel

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Size2D)                             :: AuxSize2D
        integer                                     :: STAT_CALL
        type (T_StateProp), pointer                 :: ObjWaterProperty

        !----------------------------------------------------------------------

        !1. Point Hydrodynamic module variables to original space in this module
        call ReSetHydrodynamicProperties(Me%ObjHydrodynamic, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR01'

        AuxSize2D%ILB = Me%Size%ILB 
        AuxSize2D%IUB = Me%Size%IUB 
        AuxSize2D%JLB = Me%Size%JLB 
        AuxSize2D%JUB = Me%Size%JUB 

        !2. Copy Me%FS values to Hydrodynamic module variables
        !(reestablish the original variable values)
        !Copy subroutines are individualized because they are used in case for analysis too
        call CopyWaterLevel(Me%ObjHydrodynamic, AuxSize2D, Me%FS%WaterLevel, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR02'
        
        call CopyHorizontalVelocity(Me%ObjHydrodynamic, Me%Size,                &
                                    Velocity_U = Me%FS%VelocityU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR03'

        call CopyHorizontalVelocity(Me%ObjHydrodynamic, Me%Size,                &
                                    Velocity_V = Me%FS%VelocityV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR04'

        call CopyHorizontalVelocity(Me%ObjHydrodynamic, Me%Size,                &
                                    OldVelocity_U = Me%FS%VelocityUOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR05'

        call CopyHorizontalVelocity(Me%ObjHydrodynamic, Me%Size,                &
                                    OldVelocity_V = Me%FS%VelocityVOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR06'

        call CopyVerticalVelocity(Me%ObjHydrodynamic, Me%Size,                  &
                                  Velocity_Across = Me%FS%VelocityAcross,       &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR07'

        call CopyWaterFluxes(Me%ObjHydrodynamic, Me%Size,                       &
                             WaterFluxX = Me%FS%WaterFluxX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR08'

        call CopyWaterFluxes(Me%ObjHydrodynamic, Me%Size,                       &
                             WaterFluxY = Me%FS%WaterFluxY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR09'

        call CopyWaterFluxes(Me%ObjHydrodynamic, Me%Size,                       &
                             WaterFluxZ = Me%FS%WaterFluxZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR10'

        call CopyChezyVelUV(Me%ObjHydrodynamic, AuxSize2D,                      &
                            ChezyVelUV = Me%FS%ChezyVelUV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR11'

        if (Me%HydroSubModel) then
            call CopySubModelFluxes(Me%ObjHydrodynamic, Me%Size,                &
                                    SubModelqX = Me%FS%SubModelqX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR12'

            call CopySubModelFluxes(Me%ObjHydrodynamic, Me%Size,                &
                                    SubModelqY = Me%FS%SubModelqY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR13'
        endif

        !The same for WaterProperties module
        call ReSetDensity(Me%ObjWaterProperties, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR14'

        call CopyDensity(Me%ObjWaterProperties, Me%Size, Me%FS%Density,         &
                         Me%FS%SigmaDensity, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR15'

        !(only IDs are required, instead of Me%FS could be used Me%DFS)
        ObjWaterProperty => Me%FS%FirstWaterProperty

        do while (associated(ObjWaterProperty))

            call ReSetConcentration(Me%ObjWaterProperties,                      &
                                    ObjWaterProperty%IDNumber, STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR16'

            call CopyConcentration(Me%ObjWaterProperties, Me%Size,              &
                                   ObjWaterProperty%Field,                      &
                                   ObjWaterProperty%IDNumber, STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR17'

            ObjWaterProperty => ObjWaterProperty%Next

        enddo

        !The same for Geometry module
        call ReSetGeometry(Me%ObjGeometry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR18'

        call CopyGeometryDistances(Me%ObjGeometry, Me%FS%SZZ, Me%FS%DWZ,        &
                                   Me%FS%DUZ, Me%FS%DVZ, Me%FS%DZZ,             &
                                   Me%FS%ZCellCenter, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR19'

        call CopyGeometryAreas(Me%ObjGeometry, Me%FS%AreaU, Me%FS%AreaV, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR20'

        call CopyGeometryVolumes(Me%ObjGeometry, Me%FS%VolumeZ,                 &
                                 Me%FS%VolumeZOld, Me%FS%VolumeU, Me%FS%VolumeV, &
                                 STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR21'

        call CopyGeometryWaterColumn(Me%ObjGeometry, Me%FS%WaterColumnZ,        &
                               Me%FS%WaterColumnU, Me%FS%WaterColumnV, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'CopyUndisturbedStateToModel - ModuleSequentialAssimilation - ERR22'

    end subroutine CopyUndisturbedStateToModel

    !--------------------------------------------------------------------------

    subroutine GetModelResult(CycleNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: CycleNumber

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        if (CycleNumber /= Me%StateCovRank + 1) then

            !Get ModuleHydrodynamic and ModuleWaterProperties results for respective DFS
            !variables
            call CopyModelResultToFS(Me%DFS(CycleNumber))

        else

            !Get ModuleHydrodynamic and ModuleWaterProperties results for FS variables
            call CopyModelResultToFS(Me%FS)

        endif

    end subroutine GetModelResult

    !--------------------------------------------------------------------------

    subroutine CopyModelResultToFS(ObjFullState)

        !Arguments-------------------------------------------------------------
        type (T_FullState)                          :: ObjFullState

        !Local-----------------------------------------------------------------
        real, dimension (:, :), pointer             :: AuxField2D
        real, dimension (:, :, :), pointer          :: AuxField
        real(8), dimension (:, :, :), pointer       :: AuxFieldR8
        integer                                     :: STAT_CALL
        type (T_StateProp), pointer                 :: ObjWaterProperty

        !----------------------------------------------------------------------

        !Get full model state from ModuleHydrodynamic and ModuleWaterProperties to FS
        !variables

        !ModuleHydrodynamic variables
        !>>> Water Level
        call GetWaterLevel(Me%ObjHydrodynamic, AuxField2D, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR01'

        ObjFullState%WaterLevel(:,:) = AuxField2D(:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR02'

        !>>> Velocity U
        call GetHorizontalVelocity(Me%ObjHydrodynamic,                      &
                                   Velocity_U = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR03'

        ObjFullState%VelocityU(:,:,:) = AuxField(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR04'
        
        !>>> Velocity V
        call GetHorizontalVelocity(Me%ObjHydrodynamic,                      &
                                   Velocity_V = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR05'

        ObjFullState%VelocityV(:,:,:) = AuxField(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR06'

        !>>> Velocity U Old
        call GetOldHorizontalVelocity(Me%ObjHydrodynamic,                   &
                                      OldVelocity_U = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR07'

        ObjFullState%VelocityUOld(:,:,:) = AuxField(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR08'
        
        !>>> Velocity V Old
        call GetOldHorizontalVelocity(Me%ObjHydrodynamic,                   &
                                   OldVelocity_V = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR09'

        ObjFullState%VelocityVOld(:,:,:) = AuxField(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR10'

        !>>> Velocity Across
        call GetVerticalVelocity(Me%ObjHydrodynamic,                        &
                                 Velocity_Across = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR11'

        ObjFullState%VelocityAcross(:,:,:) = AuxField(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR12'

        !>>> Water Flux X
        call GetWaterFluxes(Me%ObjHydrodynamic,                             &
                            WaterFluxX = AuxFieldR8, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR13'

        ObjFullState%WaterFluxX(:,:,:) = AuxFieldR8(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxFieldR8, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR14'

        !>>> Water Flux Y
        call GetWaterFluxes(Me%ObjHydrodynamic,                             &
                            WaterFluxY = AuxFieldR8, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR15'

        ObjFullState%WaterFluxY(:,:,:) = AuxFieldR8(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxFieldR8, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR16'

        !>>> Water Flux Z
        call GetWaterFluxes(Me%ObjHydrodynamic,                             &
                            WaterFluxZ = AuxFieldR8, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR17'

        ObjFullState%WaterFluxZ(:,:,:) = AuxFieldR8(:,:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxFieldR8, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR18'

        !>>> ChezyVelUV
        call GetChezyVelUV(Me%ObjHydrodynamic,                              &
                           ChezyVelUV = AuxField2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR19'

        ObjFullState%ChezyVelUV(:,:) = AuxField2D(:,:)

        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR20'

        if (Me%HydroSubModel) then

            !>>> SubModel fluxes
            call GetSubModelFluxes(Me%ObjHydrodynamic,                      &
                                   SubModelqX = AuxFieldR8, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR21'

            ObjFullState%SubModelqX(:,:,:) = AuxFieldR8(:,:,:)

            call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxFieldR8, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR22'

            call GetSubModelFluxes(Me%ObjHydrodynamic,                      &
                                   SubModelqY = AuxFieldR8, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR23'

            ObjFullState%SubModelqY(:,:,:) = AuxFieldR8(:,:,:)

            call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxFieldR8, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR24'

        endif

        !ModuleWaterProperties variables
        !>>> Density
        call GetDensity(Me%ObjWaterProperties, AuxField, Me%CurrentTime, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR25'

        ObjFullState%Density(:,:,:) = AuxField(:,:,:)

        call UngetWaterProperties(Me%ObjWaterProperties, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR26'

        !>>> SigmaDensity
        call GetSigma(Me%ObjWaterProperties, AuxField, Me%CurrentTime, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR27'

        ObjFullState%SigmaDensity(:,:,:) = AuxField(:,:,:)

        call UngetWaterProperties(Me%ObjWaterProperties, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR28'

        !>>> Properties
        ObjWaterProperty => ObjFullState%FirstWaterProperty

        do while (associated(ObjWaterProperty))

            call GetConcentration(Me%ObjWaterProperties, AuxField,          &
                                  ObjWaterProperty%IDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR29'

            ObjWaterProperty%Field(:,:,:) = AuxField(:,:,:)

            call UngetWaterProperties(Me%ObjWaterProperties, AuxField, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR30'

            ObjWaterProperty => ObjWaterProperty%Next

        enddo

        !ModuleGeometry variables
        !>>> Distances
        !SZZ
        call GetGeometryDistances(Me%ObjGeometry, SZZ = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR31'

        ObjFullState%SZZ(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR32'

        !DWZ
        call GetGeometryDistances(Me%ObjGeometry, DWZ = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR33'

        ObjFullState%DWZ(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR34'

        !DUZ
        call GetGeometryDistances(Me%ObjGeometry, DUZ = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR35'

        ObjFullState%DUZ(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR36'

        !DVZ
        call GetGeometryDistances(Me%ObjGeometry, DVZ = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR37'

        ObjFullState%DVZ(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR38'

        !DZZ
        call GetGeometryDistances(Me%ObjGeometry, DZZ = AuxField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR39'

        ObjFullState%DZZ(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR40'

        !ZCellCenter
        call GetGeometryDistances(Me%ObjGeometry, ZCellCenter = AuxField,   &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR41'

        ObjFullState%ZCellCenter(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR42'

        !AreaU
        call GetGeometryAreas(Me%ObjGeometry, AreaU = AuxField,             &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR43'

        ObjFullState%AreaU(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR44'

        !AreaV
        call GetGeometryAreas(Me%ObjGeometry, AreaV = AuxField,             &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR45'

        ObjFullState%AreaV(:,:,:) = AuxField(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR46'

        !VolumeZ
        call GetGeometryVolumes(Me%ObjGeometry, VolumeZ = AuxFieldR8,       &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR47'

        ObjFullState%VolumeZ(:,:,:) = AuxFieldR8(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxFieldR8, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR48'

        !VolumeZOld
        call GetGeometryVolumes(Me%ObjGeometry, VolumeZOld = AuxFieldR8,    &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR49'

        ObjFullState%VolumeZOld(:,:,:) = AuxFieldR8(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxFieldR8, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR50'

        !VolumeU
        call GetGeometryVolumes(Me%ObjGeometry, VolumeU = AuxFieldR8,       &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR51'

        ObjFullState%VolumeU(:,:,:) = AuxFieldR8(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxFieldR8, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR52'

        !VolumeV
        call GetGeometryVolumes(Me%ObjGeometry, VolumeV = AuxFieldR8,       &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR53'

        ObjFullState%VolumeV(:,:,:) = AuxFieldR8(:,:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxFieldR8, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR54'

        !WaterColumnU
        call GetGeometryWaterColumn(Me%ObjGeometry, WaterColumnU = AuxField2D, &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR55'

        ObjFullState%WaterColumnU(:,:) = AuxField2D(:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField2D, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR56'

        !WaterColumnV
        call GetGeometryWaterColumn(Me%ObjGeometry, WaterColumnV = AuxField2D, &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR57'

        ObjFullState%WaterColumnV(:,:) = AuxField2D(:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField2D, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR58'

        !WaterColumnZ
        call GetGeometryWaterColumn(Me%ObjGeometry, WaterColumn = AuxField2D, &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR59'

        ObjFullState%WaterColumnZ(:,:) = AuxField2D(:,:)

        call  UnGetGeometry(Me%ObjGeometry, AuxField2D, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'CopyModelResultToFS - ModuleSequentialAssimilation - ERR60'

    end subroutine CopyModelResultToFS

    !--------------------------------------------------------------------------

    subroutine CovarianceCalculationSEEK

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: neof, i, j, k
        integer                                     :: Count
        real, dimension (:, :), pointer             :: AuxMatrixDFS2D, AuxMatrixFS2D
        real, dimension (:, :, :), pointer          :: AuxMatrixDFS, AuxMatrixFS
        real(8), dimension (:, :, :), pointer       :: AuxMatrixDFSR8, AuxMatrixFSR8
        type (T_StateProp) , pointer                :: ObjDFSProperty, ObjFSProperty
        type(T_Property), pointer                   :: ObjProperty
        integer, dimension(:,:,:), pointer          :: PropertyMap      => null()
        !integer, dimension(:,:,:), pointer          :: PropertyOpenMap  => null()
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        !STAT_ = UNKNOWN_

        !call Ready(SequentialAssimilationID, ready_)

        !if (ready_ .EQ. IDLE_ERR_) then

            if (.not. Me%PerturbeState .or. Me%StateCov_Evolution_Split) then

                !Get WaterPoints
                call GetWaterPoints3D(Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                &
                    stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR01'

                !Get OpenPoints
                call GetOpenPoints3D(Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                &
                    stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR02'

                call CopyModelResultToFS(Me%FS)

                do neof = 1, Me%StateCovRank !Cycle the EOFs

                    ObjProperty => Me%FirstProperty

                    Count = 0 !indicator of property position in EOF vector

                    do while (associated (ObjProperty))

                        !Select state property matrixes
                        if (ObjProperty%ModuleType == 1) then !ModuleHydrodynamic property

        case1 :             select case(ObjProperty%ID%IDNumber)
                                case(WaterLevel_    )

                                    AuxMatrixDFS2D => Me%DFS(neof)%WaterLevel 
                                    AuxMatrixFS2D => Me%FS%WaterLevel
 
                                case(VelocityU_     )

                                    AuxMatrixDFS => Me%DFS(neof)%VelocityU 
                                    AuxMatrixFS => Me%FS%VelocityU

                                case(VelocityV_     ) 

                                    AuxMatrixDFS => Me%DFS(neof)%VelocityV 
                                    AuxMatrixFS => Me%FS%VelocityV

                                case(VelocityW_     )

                                    AuxMatrixDFS => Me%DFS(neof)%VelocityAcross
                                    AuxMatrixFS => Me%FS%VelocityAcross

                                case(WaterFluxX_    )

                                    AuxMatrixDFSR8 => Me%DFS(neof)%WaterFluxX
                                    AuxMatrixFSR8 => Me%FS%WaterFluxX

                                case(WaterFluxY_    )

                                    AuxMatrixDFSR8 => Me%DFS(neof)%WaterFluxY
                                    AuxMatrixFSR8 => Me%FS%WaterFluxY

                                case default

                            stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR03' 

                            end select case1

                        else !ModuleWaterProperties property

                            !All water properties are assumed 3D!

                            !Cycle water properties list and find the correct matrix to point to
                            ObjDFSProperty => Me%DFS(neof)%FirstWaterProperty
                            ObjFSProperty => Me%FS%FirstWaterProperty

                            do while (associated (ObjDFSProperty) .and. associated (ObjFSProperty))

                                if ((ObjProperty%ID%IDNumber == ObjDFSProperty%IDNumber)    &
                                    .and.                                                   &
                                    (ObjProperty%ID%IDNumber == ObjFSProperty%IDNumber)) then

                                    AuxMatrixDFS => ObjDFSProperty%Field
                                    AuxMatrixFS => ObjFSProperty%Field
                        
                                    exit !cycle proceeds only till properties found 

                                endif

                                ObjDFSProperty => ObjDFSProperty%Next
                                ObjFSProperty => ObjFSProperty%Next
                            enddo

                            if (.not. associated(AuxMatrixFS))                              &
                            stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR04' 

                        endif

                        !Cycle all property's variables
                        if (ObjProperty%Dim == Dim_2D) then

                            do j = ObjProperty%Window2D%JLB, ObjProperty%Window2D%JUB
                                do i = ObjProperty%Window2D%ILB, ObjProperty%Window2D%IUB

                                    if (Me%External_Var%WaterPoints(i,j,Me%WorkSize%KUB)    &
                                        == 1) then
                                        !found state variable!

                                        Count = Count + 1 !Count advances per state variable

                                        !if (Me%External_Var%OpenPoints(i,j,Me%WorkSize%KUB) &
                                        !    == 1) then

                                            !Calculate new EOF property value
                                            Me%LMatrix(Count,neof) = (AuxMatrixDFS2D(i,j) - &
                                                                      AuxMatrixFS2D(i,j))/  &
                                                                      Me%StateCov_LinFactor
                                        !else
                                        !    Me%LMatrix(Count,neof) = 0.
                                        !endif
                                    endif
                                enddo
                            enddo   

                            AuxMatrixDFS2D(:,:) = FillValueReal
                            if (neof == Me%StateCovRank) AuxMatrixFS2D(:,:) = FillValueReal

                            nullify(AuxMatrixDFS2D)
                            nullify(AuxMatrixFS2D)

                        else !(Dim_3D)

                            select case(ObjProperty%TypeZUV)
                                case(TypeZ_)
                                    PropertyMap     => Me%External_Var%WaterPoints
                                    !PropertyOpenMap => Me%External_Var%OpenPoints

                                case(TypeU_)
                                    !if (.not. associated(Me%External_Var%WetFacesU)) then
                                    !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU, &
                                    !                     Me%External_Var%WetFacesV, STAT = STAT_CALL)
                                    !    if (STAT_CALL /= SUCCESS_)                          &
                            !stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR05'
                                    !endif
                                    PropertyMap     => Me%WaterFacesU
                                    !PropertyOpenMap => Me%External_Var%WetFacesU

                                case(TypeV_)
                                    !if (.not. associated(Me%External_Var%WetFacesV)) then
                                    !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU, &
                                    !                     Me%External_Var%WetFacesV, STAT = STAT_CALL)
                                    !    if (STAT_CALL /= SUCCESS_)                          &
                            !stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR06'
                                    !endif
                                    PropertyMap     => Me%WaterFacesV
                                    !PropertyOpenMap => Me%External_Var%WetFacesV
                            end select

                            if ((ObjProperty%ID%IDNumber == WaterFluxX_) .or.               &
                                (ObjProperty%ID%IDNumber == WaterFluxY_)) then

                                do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                                    do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                                        do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB

                                            if (PropertyMap(i,j,k) == 1) then
                                                !found state variable!

                                                Count = Count + 1 !Count advances per state variable
                            
                                                !if (PropertyOpenMap(i,j,k) == 1) then

                                                    !Calculate new EOF property value
                                                    Me%LMatrix(Count, neof) =               &
                                                                   (AuxMatrixDFSR8(i,j,k) - &
                                                                   AuxMatrixFSR8(i,j,k))/   &
                                                                   Me%StateCov_LinFactor
                                                !else
                                                !    Me%LMatrix(Count,neof) = 0.
                                                !endif
                                            endif
                                        enddo
                                    enddo
                                enddo

                                AuxMatrixDFSR8(:,:,:) = FillValueReal
                                if (neof == Me%StateCovRank) AuxMatrixFSR8(:,:,:) = FillValueReal

                                nullify(AuxMatrixDFSR8)
                                nullify(AuxMatrixFSR8)

                            else

                                do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                                    do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                                        do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB

                                            if (PropertyMap(i,j,k) == 1) then
                                                !found state variable!

                                                Count = Count + 1 !Count advances per state variable
                            
                                                !if (PropertyOpenMap(i,j,k) == 1) then

                                                    !Calculate new EOF property value
                                                    Me%LMatrix(Count, neof) =               &
                                                                    (AuxMatrixDFS(i,j,k) -  &
                                                                    AuxMatrixFS(i,j,k))/    &
                                                                    Me%StateCov_LinFactor
                                                !else
                                                !    Me%LMatrix(Count,neof) = 0.
                                                !endif
                                            endif
                                        enddo
                                    enddo
                                enddo

                                AuxMatrixDFS(:,:,:) = FillValueReal
                                if (neof == Me%StateCovRank) AuxMatrixFS(:,:,:) = FillValueReal

                                nullify(AuxMatrixDFS)
                                nullify(AuxMatrixFS)

                            endif
                        endif

                        ObjProperty => ObjProperty%Next

                    enddo
                enddo

                !if (Me%ReOrthonormalization .and. Me%CurrentTime >= Me%ReOrthonormTime) then

                !    call ReOrthonormalization

                !    Me%ReOrthonormTime = Me%CurrentTime + Me%ReOrthonorm_DT

                !endif

                !call output of EOF
                if (Me%EOF_OutPut%HDF5ON) call EOFOutPut

                call UnGetMap (Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR07'

                call UnGetMap (Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR08'

                !if (associated(Me%External_Var%WetFacesU)) then
                !    call UnGetMap (Me%ObjMap, Me%External_Var%WetFacesU, STAT = STAT_CALL)
                !    if (STAT_CALL /= SUCCESS_)                                              &
                !        stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR09'
                !endif

                !if (associated(Me%External_Var%WetFacesV)) then
                !    call UnGetMap (Me%ObjMap, Me%External_Var%WetFacesV, STAT = STAT_CALL)
                !    if (STAT_CALL /= SUCCESS_)                                              &
                !        stop 'CovarianceCalculationSEEK - ModuleSequentialAssimilation - ERR10'
                !endif

                Me%PerturbeState = .true.

            !else

            !    Me%PerturbeState = .false.
            endif

        !    STAT_ = SUCCESS_
        !else               
        !    STAT_ = ready_
        !end if

        !if (present(STAT)) STAT = STAT_

    end subroutine CovarianceCalculationSEEK

    !--------------------------------------------------------------------------

    subroutine ReOrthonormalization

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:, :),  pointer             :: AuxCov
        real, dimension (:,:), pointer              :: AuxChol
        real, dimension (:,:), pointer              :: AuxU
        real, dimension (:,:), pointer              :: AuxI
        real(8), dimension (:,:), pointer           :: Aux1, Aux2
        real(8), dimension (:,:), pointer           :: Aux3, Aux4, Aux5
        real(8), dimension (:,:), pointer           :: AuxNeofLMatrix, AuxLMatrix
        integer                                     :: neof

        !----------------------------------------------------------------------

        allocate (AuxCov(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (AuxChol(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (AuxU(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (AuxI(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (Aux1(1:Me%StateCovRank,1:1))
        allocate (Aux2(1:Me%StateCovRank,1:1))
        allocate (Aux3(1:Me%StateCovRank,1:Me%StateVarNumber))
        allocate (Aux4(1:Me%StateVarNumber,1:Me%StateCovRank))
        allocate (Aux5(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (AuxNeofLMatrix(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (AuxLMatrix(1:Me%StateVarNumber,1:Me%StateCovRank))

        AuxCov (:,:)    = FillValueReal
        AuxChol (:,:)   = 0.
        AuxU (:,:)      = FillValueReal
        AuxI (:,:)      = 0.
        Aux1 (:,:)      = FillValueReal
        Aux2 (:,:)      = FillValueReal
        Aux3 (:,:)      = FillValueReal
        Aux4 (:,:)      = FillValueReal
        Aux5(:,:)       = FillValueReal
        AuxNeofLMatrix (:,:) = 0.

        !1. Calculate Cholesky decomposition of U
        !Calculate identity matrix
        do neof = 1, Me%StateCovRank 
            AuxI(neof,neof) = 1.
        enddo

        !Calculate U from invU*U = I, U = inv(invU), I = identity matrix
        !invU is assumed a symmetric positive definite matrix (theoretically is)
        do neof = 1, Me%StateCovRank
            Aux1(:,1) = AuxI(:,neof)
            Aux2(:,1) = FillValueReal
            call CholLinSystemSolver(Me%InvUMatrix, Aux1, Me%StateCovRank,      &
                                     Aux2)
            AuxU(:,neof) = Aux2(:,1)
        enddo

        !calculate the Cholesky factorization of InvU
        call CholeskyFactorization(AuxU, Me%StateCovRank, AuxChol)

        !2. Calculate the covariance matrix (in neof space)
        Aux3 = MATMUL(TRANSPOSE(AuxChol),TRANSPOSE(Me%LMatrix))
        Aux4 = MATMUL(Me%LMatrix,AuxChol)

        AuxCov = MATMUL(Aux3, Aux4)

        !3. Perform EOF analysis (in neof space)
        AuxU(:,:) = 0.
        call EOFAnalysis (AuxCov, AuxNeofLMatrix, AuxU, Me%StateCovRank,        &
                          Me%StateCovRank)

        !4. Translate EOFs for state variables space        
        AuxLMatrix(:,:) = Me%LMatrix(:,:)
        Aux5 = MATMUL(AuxChol, AuxNeofLMatrix)
        Me%LMatrix = MATMUL(AuxLMatrix, Aux5) 
        !(already multiplied by square root of new U)
        
        !set InvU as identity matrix
        Me%InvUMatrix(:,:) = 0.
        do neof = 1, Me%StateCovRank 

            Me%InvUMatrix(neof,neof) = 1.

        enddo

        deallocate (AuxCov)
        deallocate (AuxChol)
        deallocate (AuxU)
        deallocate (AuxI)
        deallocate (Aux1)
        deallocate (Aux2)
        deallocate (Aux3)
        deallocate (Aux4)
        deallocate (Aux5)
        deallocate (AuxNeofLMatrix)
        deallocate (AuxLMatrix)

    end subroutine ReOrthonormalization

    !--------------------------------------------------------------------------

    subroutine ModifySequentialAssimilation(SequentialAssimilationID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: SequentialAssimilationID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        real                                        :: Year, Month, Day, Hour 
        real                                        :: Minute, Second
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SequentialAssimilationID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ModifySequentialAssimilation - ModuleSequentialAssimilation - ERR01'

            !Verify existences and check measured state
            call ReadMeasures

            if (Me%CurrentMeasuresNumber /= Me%PreviousMeasuresNumber) then

                !Construct new measure covariance matrix and observations operator

            endif

            !Perform data assimilation only if
            if (Me%CurrentMeasuresNumber >= Me%StateCovRank) then

                if (Me%Method == SEEK) then

                    call UCalculationSEEK

                else
                    !...
                endif

                call ConstructModelState

                if (Me%Method == SEEK) then

                    call SEEKAnalysis

                    if (Me%StateCovEvolution) then
                      
                        if (Me%ReNormalization) then

                            call ReNormalization

                        endif

                        if (Me%ReOrthonormalization .and.                       &
                            Me%CurrentTime >= Me%ReOrthonormTime) then

                            call ReOrthonormalization

                            Me%ReOrthonormTime = Me%CurrentTime + Me%ReOrthonorm_DT

                        endif

                    endif

                else
                    !...
                endif

                call CopyModelAnalysedState

                !call output of analysis
                if (Me%OutPut%HDF5ON .or. Me%OutPut%TimeSerieON)                &
                    call SeqAssimilationOutPut

            else

                write(*,*)
                write(*,*) 'Number of measurements, ',                          &
                           Me%CurrentMeasuresNumber, ' is insuficcient.'  
                write(*,*) 'Data assimilation is not performed for time :'

100             format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)

                call ExtractDate(Me%CurrentTime, Year, Month, Day, Hour,        &
                                 Minute, Second)
                write(*,fmt=100)Year, Month, Day, Hour, Minute, Second

            endif               

            !Actualize Me%AssimilationTime
            Me%AssimilationTime = Me%CurrentTime + Me%DT

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifySequentialAssimilation

    !--------------------------------------------------------------------------

    subroutine ReadMeasures

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Measure),   pointer                  :: ObjMeasure
        type(T_Time)                                :: Time1
        real                                        :: Value1
        type(T_Time)                                :: Time2
        real                                        :: Value2
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        !Cycle measurement list and access each measurement values
        ObjMeasure => Me%FirstMeasure

        Me%CurrentMeasuresNumber = 0

        do while (associated (ObjMeasure))

            !resets the logical variable
            ObjMeasure%ValueExists = .false.

            if (ObjMeasure%FileType == TimeSerie) then

                !Get the two nearest time series values
                call GetTimeSerieValue(ObjMeasure%TimeSerie%ID, Me%CurrentTime, &
                                       ObjMeasure%TimeSerie%DataColumn, Time1,  &
                                       Value1, Time2, Value2, TimeCycle, STAT_CALL)

                !Find which is closer to current time
                if ((Me%CurrentTime - Time1) >= (Time2 - Me%CurrentTime)) then

                    !Check tolerance for time 2:
                    if ((Time2 - Me%CurrentTime) <= ObjMeasure%Tolerance) then

                        !Signals the existence of measurement (in a logical)
                        ObjMeasure%ValueExists = .true.

                        !Actualize Me%CurrentMeasureNumber
                        Me%CurrentMeasuresNumber = Me%CurrentMeasuresNumber + 1

                        !Read measurement value to Me%MeasuredState
                        Me%MeasuredState(Me%CurrentMeasuresNumber,1) = Value2
                        
                    endif
                else

                    if ((Me%CurrentTime - Time1) <= ObjMeasure%Tolerance) then

                        !Signals the existence of measurement (in a logical)
                        ObjMeasure%ValueExists = .true.

                        !Actualize Me%CurrentMeasureNumber
                        Me%CurrentMeasuresNumber = Me%CurrentMeasuresNumber + 1

                        !Read measurement value to Me%MeasuredState
                        Me%MeasuredState(Me%CurrentMeasuresNumber,1) = Value1

                    endif
                endif
            else if (ObjMeasure%FileType == HDF) then

                !...

            else

                !Not yet implemented measurement file options

            endif

            ObjMeasure => ObjMeasure%Next
        enddo

        !# Me%MeasuredState is only used till the Me%CurrentMeasureNumber line #

    end subroutine ReadMeasures

    !--------------------------------------------------------------------------

    subroutine InitialUCalculationSEEK

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real(8), dimension (:, :), pointer          :: Aux1, Aux2, Aux3

        !----------------------------------------------------------------------

        !allocate auxiliary variables
        allocate(Aux1(1:Me%MeasuresNumber, 1:Me%StateCovRank))
        allocate(Aux2(1:Me%MeasuresNumber, 1:Me%StateCovRank))
        allocate(Aux3(1:Me%StateCovRank, 1:Me%MeasuresNumber))

        Aux1 (:,:)    = FillValueReal
        Aux2 (:,:)    = FillValueReal
        Aux3 (:,:)    = FillValueReal

        !invU = L'*ObsOper'*inv(R)*(ObsOper*L)
        !Me%InvUMatrix = TRANSPOSE(Me%LMatrix)*TRANSPOSE(Me%ObservOperator)*     &
        !                Me%MeasuresInvCov*(Me%ObservOperator*Me%LMatrix)
        
        Aux1 = MATMUL(Me%ObservOperator, Me%LMatrix)
        Aux2 = MATMUL(Me%MeasuresInvCov, Aux1)
        Aux3 = MATMUL(TRANSPOSE(Me%LMatrix), TRANSPOSE(Me%ObservOperator))
        Me%InvUMatrix = MATMUL(Aux3, Aux2)

        deallocate(Aux1)
        deallocate(Aux2)
        deallocate(Aux3)

    end subroutine InitialUCalculationSEEK

    !--------------------------------------------------------------------------

    subroutine UCalculationSEEK

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real(8), dimension (:, :), pointer          :: Aux1, Aux2, Aux3, Aux4

        !----------------------------------------------------------------------

        !allocate auxiliary variables
        allocate(Aux1(1:Me%MeasuresNumber, 1:Me%StateCovRank))
        allocate(Aux2(1:Me%MeasuresNumber, 1:Me%StateCovRank))
        allocate(Aux3(1:Me%StateCovRank, 1:Me%MeasuresNumber))
        allocate(Aux4(1:Me%StateCovRank, 1:Me%StateCovRank))

        Aux1 (:,:)    = FillValueReal
        Aux2 (:,:)    = FillValueReal
        Aux3 (:,:)    = FillValueReal
        Aux4 (:,:)    = FillValueReal

        !invU = ffactor*invU + L'*ObsOper'*inv(R)*(ObsOper*L)
        !Me%InvUMatrix = Me%ForgettingFactor*Me%InvUMatrix +                     &
        !                TRANSPOSE(Me%LMatrix)*TRANSPOSE(Me%ObservOperator)*     &
        !                Me%MeasuresInvCov*(Me%ObservOperator*Me%LMatrix)
        
        Aux1 = MATMUL(Me%ObservOperator, Me%LMatrix)
        Aux2 = MATMUL(Me%MeasuresInvCov, Aux1)
        Aux3 = MATMUL(TRANSPOSE(Me%LMatrix), TRANSPOSE(Me%ObservOperator))
        Aux4 = MATMUL(Aux3, Aux2)

        Me%InvUMatrix =  Me%ForgettingFactor*Me%InvUMatrix + Aux4

        !# If Me%ForgettingFactor = 1. then model error is assumed null (perfect model)

        deallocate(Aux1)
        deallocate(Aux2)
        deallocate(Aux3)
        deallocate(Aux4)

    end subroutine UCalculationSEEK

    !--------------------------------------------------------------------------

    subroutine ConstructModelState

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty
        integer                                     :: STAT_CALL !, STAT_CALL2
        real, dimension (:, :), pointer             :: AuxField2D
        real, dimension (:, :, :), pointer          :: AuxField
        real(8), dimension (:, :, :), pointer       :: AuxFieldR8

        !----------------------------------------------------------------------

        !Get WaterPoints
        call GetWaterPoints3D(Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ConstructModelState - ModuleSequentialAssimilation - ERR01'

        !Get OpenPoints (these are get by default)
        call GetOpenPoints3D(Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ConstructModelState - ModuleSequentialAssimilation - ERR02'

        !Actualize OpenFaces
        if (associated(Me%WaterFacesU)) then
            call ActualizeOpenFaces
        endif

        !Assess assimilation state property list
        ObjProperty => Me%FirstProperty

        do while (associated (ObjProperty))

            if (ObjProperty%ModuleType == 1) then !ModuleHydrodynamic property

                !All hydrodynamic properties are assumed 3D except water level that is 2D!

case1 :         select case(ObjProperty%ID%IDNumber)
                    case(WaterLevel_    )

                        call GetWaterLevel(Me%ObjHydrodynamic, AuxField2D, STAT_CALL)
 
                    case(VelocityU_     )

                        call GetHorizontalVelocity(Me%ObjHydrodynamic,              &
                                                   Velocity_U = AuxField,           &
                                                   STAT = STAT_CALL)
                        
                        !if (.not. associated(Me%External_Var%WetFacesU)) then
                        !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU,  &
                        !                     Me%External_Var%WetFacesV, STAT = STAT_CALL2)
                        !endif

                    case(VelocityV_     ) 

                        call GetHorizontalVelocity(Me%ObjHydrodynamic,              &
                                                   Velocity_V = AuxField,           &
                                                   STAT = STAT_CALL)

                        !if (.not. associated(Me%External_Var%WetFacesV)) then
                        !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU,  &
                        !                     Me%External_Var%WetFacesV, STAT = STAT_CALL2)
                        !endif

                    case(VelocityW_     )

                        call GetVerticalVelocity(Me%ObjHydrodynamic,                &
                                                 Velocity_Across = AuxField,        &
                                                 STAT = STAT_CALL)

                    case(WaterFluxX_    )

                        call GetWaterFluxes(Me%ObjHydrodynamic,                     &
                                            WaterFluxX = AuxFieldR8, STAT = STAT_CALL)

                        !if (.not. associated(Me%External_Var%WetFacesU)) then
                        !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU,  &
                        !                     Me%External_Var%WetFacesV, STAT = STAT_CALL2)
                        !endif

                    case(WaterFluxY_    )

                        call GetWaterFluxes(Me%ObjHydrodynamic,                     &
                                            WaterFluxY = AuxFieldR8, STAT = STAT_CALL)

                        !if (.not. associated(Me%External_Var%WetFacesV)) then
                        !    call GetWetFaces(Me%ObjMap, Me%External_Var%WetFacesU,  &
                        !                     Me%External_Var%WetFacesV, STAT = STAT_CALL2)
                        !endif

                    case default

                        stop 'ConstructModelState - ModuleSequentialAssimilation - ERR03' 

                end select case1

                if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'ConstructModelState - ModuleSequentialAssimilation - ERR04'

                !if (STAT_CALL2 /= SUCCESS_)                                         &
                !    stop 'ConstructModelState - ModuleSequentialAssimilation - ERR05'

                if (ObjProperty%Dim == Dim_2D) then

                    call SetMatrixValue(ObjProperty%Field2D, ObjProperty%Window2D,  &
                                        AuxField2D)

                    call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField2D,          &
                                           STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'ConstructModelState - ModuleSequentialAssimilation - ERR06'

                else

                    if ((ObjProperty%ID%IDNumber == WaterFluxX_) .or.               &
                        (ObjProperty%ID%IDNumber == WaterFluxY_)) then

                        call SetMatrixValue(ObjProperty%FieldR8,                    &
                                            ObjProperty%Window, AuxFieldR8)

                        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxFieldR8,      &
                                               STAT = STAT_CALL) 
                        if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'ConstructModelState - ModuleSequentialAssimilation - ERR07'

                    else
                    
                        call SetMatrixValue(ObjProperty%Field, ObjProperty%Window,  &
                                            AuxField)

                        call UnGetHydrodynamic(Me%ObjHydrodynamic, AuxField,        &
                                               STAT = STAT_CALL) 
                        if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'ConstructModelState - ModuleSequentialAssimilation - ERR08'

                    endif
                endif
            else !ModuleWaterProperties property

                !All water properties are assumed 3D!
                call GetConcentration(Me%ObjWaterProperties, AuxField,              &
                                      ObjProperty%ID%IDNumber, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'ConstructModelState - ModuleSequentialAssimilation - ERR09'

                call SetMatrixValue(ObjProperty%Field, ObjProperty%Window, AuxField)

                call UnGetWaterProperties(Me%ObjWaterProperties, AuxField, STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'ConstructModelState - ModuleSequentialAssimilation - ERR10'
            endif

            call AddPropertyToState(ObjProperty)

            !!Nullify parameters' field
            !if (ObjProperty%Dim == Dim_2D) then

            !    call SetMatrixValue(ObjProperty%Field2D, ObjProperty%Window2D,      &
            !                        FillValueReal)
            !else

            !    if ((ObjProperty%ID%IDNumber == WaterFluxX_) .or.                   &
            !        (ObjProperty%ID%IDNumber == WaterFluxY_)) then

            !        call SetMatrixValue(ObjProperty%FieldR8, ObjProperty%Window,    &
            !                            dble(FillValueReal))

            !    else

            !        call SetMatrixValue(ObjProperty%Field, ObjProperty%Window,      &
            !                            FillValueReal)

            !    endif

            !endif

            ObjProperty => ObjProperty%Next

        enddo

        call UnGetMap (Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructModelState - ModuleSequentialAssimilation - ERR11'

        call UnGetMap (Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ConstructModelState - ModuleSequentialAssimilation - ERR12'

        !if (associated(Me%External_Var%WetFacesU)) then
        !    call UnGetMap (Me%ObjMap, Me%External_Var%WetFacesU, STAT = STAT_CALL)
        !    if (STAT_CALL /= SUCCESS_)                                              &
        !        stop 'ConstructModelState - ModuleSequentialAssimilation - ERR13'
        !endif

        !if (associated(Me%External_Var%WetFacesV)) then
        !    call UnGetMap (Me%ObjMap, Me%External_Var%WetFacesV, STAT = STAT_CALL)
        !    if (STAT_CALL /= SUCCESS_)                                              &
        !        stop 'ConstructModelState - ModuleSequentialAssimilation - ERR14'
        !endif

    end subroutine ConstructModelState

    !--------------------------------------------------------------------------

    subroutine SEEKAnalysis

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real(8), dimension (:,:), pointer           :: AuxY, AuxX
        real(8), dimension (:,:), pointer           :: Aux1, Aux2, Aux3, Aux4
        real(8), dimension (:,:), pointer           :: Aux5, Aux6
        integer                                     :: STAT_CALL
        integer                                     :: imeas
        real, dimension(:), pointer                 :: DataBuffer
        type (T_ErrorTimeSerie), pointer            :: ErrorTimeSerie

        !----------------------------------------------------------------------

        !allocate aux variables
        allocate (Aux1 (1:Me%MeasuresNumber,1:1))
        allocate (Aux2 (1:Me%MeasuresNumber,1:1))
        allocate (Aux3 (1:Me%MeasuresNumber,1:1))
        allocate (Aux4 (1:Me%StateVarNumber,1:1))

        allocate (AuxY (1:Me%StateCovRank,1:1))
        allocate (AuxX (1:Me%StateCovRank,1:1))

        Aux1 (:,:)    = FillValueReal
        Aux2 (:,:)    = FillValueReal
        Aux3 (:,:)    = FillValueReal
        Aux4 (:,:)    = FillValueReal

        AuxY (:,:)    = FillValueReal
        AuxX (:,:)    = FillValueReal

        !Calculate Y = L'*ObservOper'*inv(R)*(Measures - ObservOper*State)
        !AuxY = TRANSPOSE(Me%LMatrix)*TRANSPOSE(Me%ObservOperator)*Me%MeasuresInvCov*    &
        !       (Me%MeasuredState - Me%ObservOperator*Me%State)

        Aux1 = MATMUL(Me%ObservOperator, Me%State)
        Aux2 = Me%MeasuredState - Aux1
        Aux3 = MATMUL(Me%MeasuresInvCov, Aux2)
        Aux4 = MATMUL(TRANSPOSE(Me%ObservOperator), Aux3)

        AuxY = MATMUL(TRANSPOSE(Me%LMatrix), Aux4)

        !Calculate X from U*X = Y, U = invU
        !invU is assumed a symmetric positive definite matrix (theoretically is)
        call CholLinSystemSolver(Me%InvUMatrix, AuxY, Me%StateCovRank, AuxX)

        !Calculate State = State + L*X (correction only made in OpenPoints)
        Me%State = Me%State + Me%OpenPointState*MATMUL(Me%LMatrix, AuxX)

        if (Me%OutPut%ErrorTimeSerieON) then

            !Allocates DataBuffer
            allocate(DataBuffer(2), STAT = STAT_CALL)
            if (STAT_CALL /= 0)                                                     &
                stop 'SEEKAnalysis - ModuleSequentialAssimilation - ERR01'
        
            ErrorTimeSerie => Me%FirstErrorTimeSerie

            allocate (Aux5 (1:Me%MeasuresNumber,1:1))
            allocate (Aux6 (1:Me%MeasuresNumber,1:1))

            Aux5 (:,:)    = FillValueReal
            Aux6 (:,:)    = FillValueReal

            Aux5 = MATMUL(Me%ObservOperator, Me%State)
            Aux6 = Me%MeasuredState - Aux5

            do imeas = 1, Me%MeasuresNumber

                DataBuffer(1) = Aux2(imeas,1)
                DataBuffer(2) = Aux6(imeas,1)

                call WriteTimeSerieLine(ErrorTimeSerie%TimeSerie, DataBuffer,       &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SEEKAnalysis - ModuleSequentialAssimilation - ERR02'
                ErrorTimeSerie => ErrorTimeSerie%Next
            enddo

            deallocate(DataBuffer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'SEEKAnalysis - ModuleSequentialAssimilation - ERR03'

            deallocate (Aux5)
            deallocate (Aux6)

        endif

        !deallocate aux variables
        deallocate (AuxY)
        deallocate (AuxX)

        deallocate (Aux1)
        deallocate (Aux2)
        deallocate (Aux3)
        deallocate (Aux4)

        !Reference: Hoteit, Ibrahim, 2001, Filtres de Kalman Réduits Efficaces pour
        !           l'Assimilation de Données en Oceanographie, Thèse de Docteur en 
        !           Mathématiques Appliquées de l'Université de Joseph Fourrier

    end subroutine SEEKAnalysis

    !--------------------------------------------------------------------------

    subroutine AddPropertyToState(ObjProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, Position
        integer, dimension(:,:,:), pointer          :: PropertyMap      => null()
        integer, dimension(:,:,:), pointer          :: PropertyOpenMap  => null()

        !----------------------------------------------------------------------

        !Add property field to State vector
        Position = ObjProperty%FirstStatePosition

        if (ObjProperty%Dim == Dim_2D) then

            do j = ObjProperty%Window2D%JLB, ObjProperty%Window2D%JUB
                do i = ObjProperty%Window2D%ILB, ObjProperty%Window2D%IUB
                    if (Me%External_Var%WaterPoints(i,j,Me%WorkSize%KUB) == 1) then
                        !found state variable!

                        Me%State(Position,1) = ObjProperty%Field2D(i,j)

                        Me%OpenPointState(Position,1) =                             &
                            Me%External_Var%WaterPoints(i,j,Me%WorkSize%KUB)
                        !(if 2D assumed to be water level)

                        Position = Position + 1
                    endif
                enddo
            enddo

        else !(Dim_3D)

            select case(ObjProperty%TypeZUV)
                case(TypeZ_)
                    PropertyMap     => Me%External_Var%WaterPoints
                    PropertyOpenMap => Me%External_Var%OpenPoints

                case(TypeU_)
                    PropertyMap     => Me%WaterFacesU
                    PropertyOpenMap => Me%OpenFacesU

                case(TypeV_)
                    PropertyMap     => Me%WaterFacesV
                    PropertyOpenMap => Me%OpenFacesV
            end select

            if ((ObjProperty%ID%IDNumber == WaterFluxX_) .or.                       &
                (ObjProperty%ID%IDNumber == WaterFluxY_)) then

                do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                    do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                        do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB
                            if (PropertyMap(i,j,k) == 1) then
                                !found state variable!

                                Me%State(Position,1) = ObjProperty%FieldR8(i,j,k)

                                Me%OpenPointState(Position,1) = PropertyOpenMap(i,j,k)

                                Position = Position + 1
                            endif
                        enddo
                    enddo
                enddo

            else

                do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                    do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                        do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB
                            if (PropertyMap(i,j,k) == 1) then
                                !found state variable!

                                Me%State(Position,1) = ObjProperty%Field(i,j,k)

                                Me%OpenPointState(Position,1) = PropertyOpenMap(i,j,k)

                                Position = Position + 1
                            endif
                        enddo
                    enddo
                enddo

            endif

        endif

        !(this subroutine explains the state filling method: klb->kub, 
        ! jlb->jub, ilb->iub)

    end subroutine AddPropertyToState

    !--------------------------------------------------------------------------

    subroutine ActualizeOpenFaces

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        !integer                                         :: auxi, auxj
        integer                                         :: i, j, k
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: STAT_CALL

        !----------------------------------------------------------------------

        !Calculates the OpenFaces based on ComputeFaces, ImposedNormalFaces and
        !ImposedTangentialFaces

        call GetComputeFaces3D(Me%ObjMap, Me%External_Var%ComputeFacesU,            &
                               Me%External_Var%ComputeFacesV, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR01'

        call GetImposedNormalFaces(Me%ObjMap, Me%External_Var%ImposedNormFacesU,    &
                                   Me%External_Var%ImposedNormFacesV, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR02'

        call GetImposedTangentialFaces(Me%ObjMap, Me%External_Var%ImposedTangFacesU, &
                                   Me%External_Var%ImposedTangFacesV, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR03'

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !By default all values are zero
        Me%OpenFacesU       = 0
        Me%OpenFacesV       = 0

        !auxj = JLB
        !auxi = ILB

        do k = KLB, KUB

            !OpenFacesU
            do j = JLB, JUB + 1
                do i = ILB, IUB
                    if (Me%External_Var%ComputeFacesU(i,j,k) == 1 .or.              &
                        Me%External_Var%ImposedNormFacesU(i,j,k) == 1 .or.          &
                        Me%External_Var%ImposedTangFacesU(i,j,k) == 1) then
                    
                        Me%OpenFacesU(i,j,k) = 1

                    endif
                enddo
            enddo

            !OpenFacesV
            do j = JLB, JUB
                do i = ILB, IUB + 1
                    if (Me%External_Var%ComputeFacesV(i,j,k) == 1 .or.              &
                        Me%External_Var%ImposedNormFacesV(i,j,k) == 1 .or.          &
                        Me%External_Var%ImposedTangFacesV(i,j,k) == 1) then
                    
                        Me%OpenFacesV(i,j,k) = 1

                    endif
                enddo
            enddo

        enddo

        call UnGetMap (Me%ObjMap, Me%External_Var%ComputeFacesU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR04'

        call UnGetMap (Me%ObjMap, Me%External_Var%ComputeFacesV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR05'

        call UnGetMap (Me%ObjMap, Me%External_Var%ImposedNormFacesU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR06'

        call UnGetMap (Me%ObjMap, Me%External_Var%ImposedNormFacesV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR07'

        call UnGetMap (Me%ObjMap, Me%External_Var%ImposedTangFacesU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR08'

        call UnGetMap (Me%ObjMap, Me%External_Var%ImposedTangFacesV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ActualizeOpenFaces - ModuleSequentialAssimilation - ERR09'

    end subroutine ActualizeOpenFaces

    !--------------------------------------------------------------------------

    subroutine ReNormalization

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        !real(8), dimension (:,:), pointer           :: Aux1
        real, dimension (:,:), pointer              :: AuxC
        real, dimension (:,:), pointer              :: AuxU
        real, dimension (:,:), pointer              :: AuxI
        real(8), dimension (:,:), pointer           :: Aux1, Aux2
        real(8), dimension (:,:), pointer           :: AuxLMatrix
        integer                                     :: neof

        !----------------------------------------------------------------------

        allocate (AuxC(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (AuxU(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (AuxI(1:Me%StateCovRank,1:Me%StateCovRank))
        allocate (Aux1(1:Me%StateCovRank,1:1))
        allocate (Aux2(1:Me%StateCovRank,1:1))
        allocate (AuxLMatrix(1:Me%StateVarNumber,1:Me%StateCovRank))
        AuxC (:,:) = 0.
        AuxU (:,:) = FillValueReal
        AuxI (:,:) = 0.
        Aux1 (:,:) = FillValueReal
        Aux2 (:,:) = FillValueReal
        AuxLMatrix (:,:) = 0.

        !Calculate identity matrix
        do neof = 1, Me%StateCovRank 
            AuxI(neof,neof) = 1.
        enddo

        !Calculate U from invU*U = I, U = inv(invU), I = identity matrix
        !invU is assumed a symmetric positive definite matrix (theoretically is)
        do neof = 1, Me%StateCovRank
            Aux1(:,1) = AuxI(:,neof)
            Aux2(:,1) = FillValueReal
            call CholLinSystemSolver(Me%InvUMatrix, Aux1, Me%StateCovRank,  &
                                    Aux2)
            AuxU(:,neof) = Aux2(:,1)
        enddo

        !calculate the Cholesky factorization of InvU
        call CholeskyFactorization(AuxU, Me%StateCovRank, AuxC)

        AuxLMatrix(:,:) = Me%LMatrix(:,:)

        !calculate the new LMatrix
        Me%LMatrix = MATMUL(AuxLMatrix, AuxC)

        !set InvU as identity matrix
        Me%InvUMatrix(:,:) = 0.
        do neof = 1, Me%StateCovRank 

            Me%InvUMatrix(neof,neof) = 1.

        enddo

        deallocate (AuxC)
        deallocate (AuxU)
        deallocate (AuxI)
        deallocate (Aux1)
        deallocate (Aux2)
        deallocate (AuxLMatrix)

    end subroutine ReNormalization

    !--------------------------------------------------------------------------

    subroutine CopyModelAnalysedState

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: ObjProperty
        integer                                     :: i, j, k, Position
        integer                                     :: STAT_CALL
        integer, dimension(:,:,:), pointer          :: PropertyMap  => null()
        !integer, dimension(:, :), pointer           :: OpenMap2D

        !----------------------------------------------------------------------

        !Get WaterPoints
        call GetWaterPoints3D(Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR01'

        !!Get OpenPoints
        !call GetOpenPoints3D(Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
        !if (STAT_CALL .NE. SUCCESS_)                                                &
        !    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR02'

        !Assess assimilation state property list
        ObjProperty => Me%FirstProperty

        do while (associated (ObjProperty))

            Position = ObjProperty%FirstStatePosition

            if (ObjProperty%Dim == Dim_2D) then

                !allocate (OpenMap2D)
                !OpenMap2D = 0

                !1. Copy analysed state to property field (must be nullified before)
                do j = ObjProperty%Window2D%JLB, ObjProperty%Window2D%JUB
                    do i = ObjProperty%Window2D%ILB, ObjProperty%Window2D%IUB
                        if (Me%External_Var%WaterPoints(i,j,Me%WorkSize%KUB) == 1) then
                            !found state variable!

                            ObjProperty%Field2D(i,j) = Me%State(Position,1)

                            !OpenMap2D = Me%External_Var%OpenPoints(i,j,Me%WorkSize%KUB)

                            Position = Position + 1
                        endif
                    enddo
                enddo

                !2. Copy from property field to ModuleHydrodynamic and ModuleWaterProperties
                if (ObjProperty%ID%IDNumber == WaterLevel_) then

                    call CopyWaterLevel(Me%ObjHydrodynamic, ObjProperty%Window2D,   &
                                        ObjProperty%Field2D, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR03'
                endif

                !3. Nullify property field
                if (.not. Me%OutPut%HDF5ON) then
                    call SetMatrixValue(ObjProperty%Field2D, ObjProperty%Window2D,  &
                                        FillValueReal)
                endif

                !deallocate(OpenMap2D)
                !nullify(OpenMap2D)

            else

                select case(ObjProperty%TypeZUV)
                    case(TypeZ_)
                        PropertyMap => Me%External_Var%WaterPoints

                    case(TypeU_)
                        PropertyMap => Me%WaterFacesU

                    case(TypeV_)
                        PropertyMap => Me%WaterFacesV
                end select

                if ((ObjProperty%ID%IDNumber == WaterFluxX_) .or.                   &
                    (ObjProperty%ID%IDNumber == WaterFluxY_)) then

                    do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                        do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                            do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB
                                if (PropertyMap(i,j,k) == 1) then
                                    !found state variable!

                                    ObjProperty%FieldR8(i,j,k) = Me%State(Position,1)

                                    Position = Position + 1
                                endif
                            enddo
                        enddo
                    enddo

                else

                    do k = ObjProperty%Window%KLB, ObjProperty%Window%KUB
                        do j = ObjProperty%Window%JLB, ObjProperty%Window%JUB
                            do i = ObjProperty%Window%ILB, ObjProperty%Window%IUB
                                if (PropertyMap(i,j,k) == 1) then
                                    !found state variable!

                                    if ((ObjProperty%ModuleType /= 2) .or.          &
                                        Me%State(Position,1) >= 0.) then

                                        ObjProperty%Field(i,j,k) = Me%State(Position,1)
                                    else
                                        !Assume all water properties concentration positive
                                        ObjProperty%Field(i,j,k) = 0.
                                    endif

                                    Position = Position + 1

                                endif
                            enddo
                        enddo
                    enddo

                endif

                !2. Copy from property field to ModuleHydrodynamic and ModuleWaterProperties
                if (ObjProperty%ModuleType == 1) then !ModuleHydrodynamic property

case1 :             select case(ObjProperty%ID%IDNumber)

                        case(VelocityU_     )

                            call CopyHorizontalVelocity(Me%ObjHydrodynamic,         &
                                                ObjProperty%Window,                 &
                                                Velocity_U = ObjProperty%Field,     &
                                                !MapMatrix = Me%OpenFacesU,          &
                                                STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR04'

                        case(VelocityV_     ) 

                            call CopyHorizontalVelocity(Me%ObjHydrodynamic,         &
                                                ObjProperty%Window,                 &
                                                Velocity_V = ObjProperty%Field,     &
                                                !MapMatrix = Me%OpenFacesV,          &
                                                STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR05'

                        case(VelocityW_     )

                            call CopyVerticalVelocity(Me%ObjHydrodynamic,           &
                                                      ObjProperty%Window,           &
                                            Velocity_Across = ObjProperty%Field,    &
                                            !MapMatrix = Me%External_Var%OpenPoints, &
                                                      STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR06'

                        case(WaterFluxX_    )

                            call CopyWaterFluxes(Me%ObjHydrodynamic,                &
                                                 ObjProperty%Window,                &
                                                 WaterFluxX = ObjProperty%FieldR8,  &
                                                 !MapMatrix = Me%OpenFacesU,         &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR07'

                        case(WaterFluxY_    )

                            call CopyWaterFluxes(Me%ObjHydrodynamic,                &
                                                 ObjProperty%Window,                &
                                                 WaterFluxY = ObjProperty%FieldR8,  &
                                                 !MapMatrix = Me%OpenFacesV,         &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR08'
                    end select case1

                else !ModuleWaterProperties property

                    call CopyConcentration(Me%ObjWaterProperties,                   &
                                           ObjProperty%Window,                      &
                                           ObjProperty%Field,                       &
                                           ObjProperty%ID%IDNumber, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR09'

                endif

                !3. Nullify property field
                if (.not. Me%OutPut%HDF5ON) then
                    if ((ObjProperty%ID%IDNumber == WaterFluxX_) .or.               &
                        (ObjProperty%ID%IDNumber == WaterFluxY_)) then

                        call SetMatrixValue(ObjProperty%FieldR8, ObjProperty%Window, &
                                            dble(FillValueReal))
                    else

                        call SetMatrixValue(ObjProperty%Field, ObjProperty%Window,  &
                                            FillValueReal)
                    endif
                endif
            endif

            ObjProperty => ObjProperty%Next

        enddo

        call UnGetMap (Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR10'

        !call UnGetMap (Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)                                                  &
        !    stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR11'

    end subroutine CopyModelAnalysedState

    !--------------------------------------------------------------------------

    subroutine SeqAssimilationOutPut

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Property), pointer                   :: ObjProperty
        logical                                     :: FirstTime
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WorkILB, WorkIUB, WorkJLB
        integer                                     :: WorkJUB, WorkKLB, WorkKUB
        integer                                     :: OutPutNumber
        !real, dimension (:, :, :), pointer          :: CenteredField
        !real(8), dimension (:, :, :), pointer       :: CenteredFieldR8

        !----------------------------------------------------------------------

        !Get WaterPoints
        call GetWaterPoints3D(Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR01'

        if (Me%OutPut%HDF5ON) then 

            !Bounds
            ILB = Me%Size%ILB 
            IUB = Me%Size%IUB 
            JLB = Me%Size%JLB 
            JUB = Me%Size%JUB 
            KLB = Me%Size%KLB 
            KUB = Me%Size%KUB

            WorkILB = Me%WorkSize%ILB 
            WorkIUB = Me%WorkSize%IUB 
            WorkJLB = Me%WorkSize%JLB 
            WorkJUB = Me%WorkSize%JUB 
            WorkKLB = Me%WorkSize%KLB 
            WorkKUB = Me%WorkSize%KUB 

            ObjProperty => Me%FirstProperty

            FirstTime = .true.        

            OutPutNumber = Me%OutPut%NextOutPut

    TOut:   if (Me%CurrentTime >= Me%OutPut%OutTime(OutPutNumber)) then

    PropX:      do while (associated(ObjProperty))

    First:          if (FirstTime) then 

                        !Writes current time
                        call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2),     &
                                            AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))
                        TimePtr => AuxTime

                        !Get OpenPoints
                        call GetOpenPoints3D(Me%ObjMap, Me%External_Var%OpenPoints,     &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR02'

                        !Writes Time
                        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR03'

                        call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",               &
                                             "YYYY/MM/DD HH:MM:SS", Array1D = TimePtr,  &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR04'

                        !!Writes SZZ
                        !call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                        !                     WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
                        !if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR04'

                        !call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",    &
                        !                     "m", Array3D = Me%ExternalVar%SZZ,            &
                        !                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        !if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR05'

                        !Writes OpenPoints
                        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,              &
                                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,        &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR05'

                        call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints", &
                                             "-", Array3D = Me%External_Var%OpenPoints, &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR06'
               
                        call UnGetMap (Me%ObjMap, Me%External_Var%OpenPoints, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR07'

                        FirstTime = .false.
                    endif First

                    select case (ObjProperty%Dim)

                        case(Dim_2D)

                            call HDF5SetLimits(Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,   &
                                               WorkJUB, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  & 
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR08'

                            call HDF5WriteData(Me%ObjHDF5,                              &
                                               "/Results/"//trim(ObjProperty%ID%Name),  &
                                               trim(ObjProperty%ID%Name),               &
                                               ObjProperty%ID%Units,                    &
                                               Array2D = ObjProperty%Field2D,           &
                                               OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR09'            

                            !call SetMatrixValue(ObjProperty%Field2D, ObjProperty%Window2D, &
                            !                    FillValueReal)
                        case(Dim_3D)

                            call HDF5SetLimits(Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,   &
                                               WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  & 
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR10'

                            if (ObjProperty%TypeZUV == TypeZ_) then

                                call HDF5WriteData(Me%ObjHDF5,                          &
                                                   "/Results/"//trim(ObjProperty%ID%Name), &
                                                   trim(ObjProperty%ID%Name),           &
                                                   ObjProperty%ID%Units,                &
                                                   Array3D = ObjProperty%Field,         &
                                                   OutputNumber = OutPutNumber,         &
                                                   STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)                              &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR11' 

                            else

                                if ((ObjProperty%ID%IDNumber /= WaterFluxX_) .and.      &
                                    (ObjProperty%ID%IDNumber /= WaterFluxY_)) then

                                    !allocate(CenteredField(WorkILB:WorkIUB,             &
                                    !         WorkJLB:WorkJUB,WorkKLB:WorkKUB))
                                    !CenteredField(:, :, :) = FillValueReal 

                                    call CenterField(ObjProperty%CenteredField,         &
                                                     ObjProperty%Field,                 &
                                                     ObjProperty%TypeZUV,               &
                                                     ObjProperty%Window)
                                                     !ObjProperty%Window,                &
                                                     !ObjProperty%CyclicBoundary)

                                    call HDF5WriteData(Me%ObjHDF5,                      &
                                                       "/Results/"//trim(ObjProperty%ID%Name), &
                                                       trim(ObjProperty%ID%Name),       &
                                                       ObjProperty%ID%Units,            &
                                                       Array3D = ObjProperty%CenteredField, &
                                                       OutputNumber = OutPutNumber,     &
                                                       STAT = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_)                          &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR12' 

                                    !deallocate(CenteredField)
                                    !nullify(CenteredField)

                                    !call SetMatrixValue(ObjProperty%Field,              &
                                    !                    ObjProperty%Window,             &
                                    !                    FillValueReal)
                                else

                                    !allocate(CenteredFieldR8(WorkILB:WorkIUB,           &
                                    !         WorkJLB:WorkJUB,WorkKLB:WorkKUB))
                                    !CenteredFieldR8(:, :, :) = FillValueReal 

                                    call CenterFieldR8(ObjProperty%CenteredFieldR8,     &
                                                     ObjProperty%FieldR8,               &
                                                     ObjProperty%TypeZUV,               &
                                                     ObjProperty%Window)
                                                     !ObjProperty%Window,                &
                                                     !ObjProperty%CyclicBoundary)

                                    call HDF5WriteData(Me%ObjHDF5,                      &
                                                     "/Results/"//trim(ObjProperty%ID%Name), &
                                                     trim(ObjProperty%ID%Name),         &
                                                     ObjProperty%ID%Units,              &
                                                     Array3D = ObjProperty%CenteredFieldR8, &
                                                     OutputNumber = OutPutNumber,       &
                                                     STAT = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_)                          &
                            stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR13' 

                                    !deallocate(CenteredFieldR8)
                                    !nullify(CenteredFieldR8)

                                    !call SetMatrixValue(ObjProperty%FieldR8,            &
                                    !                    ObjProperty%Window,             &
                                    !                    dble(FillValueReal))
                                endif


                            endif
                    end select

                    !Writes everything to disk
                    call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'SeqAssimilationOutPut - ModuleSequentialAssimilation - ERR14'

                    ObjProperty => ObjProperty%Next

                enddo PropX

                Me%OutPut%NextOutPut = OutPutNumber + 1

            endif  TOut    

            nullify(ObjProperty)

        endif 

        if (Me%OutPut%TimeSerieON) then 

            call OutPutTimeSeries

        endif

        call UnGetMap (Me%ObjMap, Me%External_Var%WaterPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'CopyModelAnalysedState - ModuleSequentialAssimilation - ERR15'

    end subroutine SeqAssimilationOutPut

    !--------------------------------------------------------------------------

    subroutine CenterField(PropCenter, PropFaces, TypeZUV, Size) !, CyclicBoundary)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :),   pointer             :: PropCenter, PropFaces
        integer                                         :: TypeZUV
        type (T_Size3D)                                 :: Size
        !logical                                         :: CyclicBoundary

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

        select case (TypeZUV)

            case (TypeU_)

                do k = KLB, KUB
                do j = JLB, JUB - 1
                do i = ILB, IUB
                    if (Me%External_Var%WaterPoints(i, j, k) == WaterPoint) then

                        !Check for non state faces (assume MISSING_NULL)
                        if (PropFaces(i,j,k) < FillValueReal/2.) PropFaces(i,j,k) = 0.
                        if (PropFaces(i,j+1,k) < FillValueReal/2.) PropFaces(i,j+1,k) = 0.

                        PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
                                               PropFaces(i, j+1, k)) / 2.
                    else

                        !PropCenter(i, j, k) = 0.
                    endif
                enddo
                enddo
                enddo

            case (TypeV_)

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB - 1
                    if (Me%External_Var%WaterPoints(i, j, k) == WaterPoint) then

                        !Check for non state faces (assume MISSING_NULL)
                        if (PropFaces(i,j,k) < FillValueReal/2.) PropFaces(i,j,k) = 0.
                        if (PropFaces(i+1,j,k) < FillValueReal/2.) PropFaces(i+1,j,k) = 0.

                        PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
                                               PropFaces(i+1, j, k)) / 2.
                    else

                        !PropCenter(i, j, k) = 0.
                    endif
                enddo
                enddo
                enddo

        end select

        !if (TypeZUV == TypeU_) then

            !if (.not. CyclicBoundary) then

            !    do k = KLB, KUB
            !    do j = JLB, JUB - 1
            !    do i = ILB, IUB
            !        if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then
               
            !            PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
            !                                   PropFaces(i, j+1, k)) / 2.
            !        else

            !            PropCenter(i, j, k) = 0.
            !        endif
            !    enddo
            !    enddo
            !    enddo

            !else !(cyclic boundary)

            !    do k = KLB, KUB
            !    do j = JLB, JUB - 2
            !    do i = ILB, IUB
            !        if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then
               
            !            PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
            !                                   PropFaces(i, j+1, k)) / 2.
            !        else

            !            PropCenter(i, j, k) = 0.
            !        endif
            !    enddo
            !    enddo

            !    do i = ILB, IUB
            !        if (Me%External_Var%OpenPoints(i, JUB - 1, k) == OpenPoint) then

            !            if (Me%External_Var%OpenPoints(i, JLB, k) == OpenPoint) then

            !                PropCenter(i, JUB - 1, k) = (PropFaces(i, JUB - 1, k) +     &
            !                                        PropFaces(i, JLB, k)) / 2.
            !            else

            !                PropCenter(i, JUB - 1, k) = PropFaces(i, JUB - 1, k) / 2.
            !                !(assumed null value for property at JLB face)
            !            endif
            !        else

            !            PropCenter(i, JUB - 1, k) = 0.
            !        endif
            !    enddo
            !    enddo
            !endif

        !else !TypeV_

            !if (.not. CyclicBoundary) then

            !    do k = KLB, KUB
            !    do j = JLB, JUB
            !    do i = ILB, IUB - 1
            !        if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then

            !            PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
            !                                   PropFaces(i+1, j, k)) / 2.
            !        else

            !            PropCenter(i, j, k) = 0.
            !        endif
            !    enddo
            !    enddo
            !    enddo

            !else !(cyclic boundary)

            !    do k = KLB, KUB
            !    do j = JLB, JUB
            !    do i = ILB, IUB - 2
            !        if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then
               
            !            PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
            !                                   PropFaces(i+1,j, k)) / 2.
            !        else

            !            PropCenter(i, j, k) = 0.
            !        endif
            !    enddo
            !    enddo

            !    do j = JLB, JUB
            !        if (Me%External_Var%OpenPoints(IUB - 1, j, k) == OpenPoint) then
                    
            !            if (Me%External_Var%OpenPoints(ILB, j, k) == OpenPoint) then
               
            !                PropCenter(IUB - 1, j, k) = (PropFaces(IUB - 1, j, k) +     &
            !                                        PropFaces(ILB, j, k)) / 2.
            !            else

            !                PropCenter(IUB - 1, j, k) = PropFaces(IUB - 1, j, k) / 2.
            !                !(assumed null value for property at ILB face)
            !            endif
            !        else

            !            PropCenter(IUB - 1, j, k) = 0.
            !        endif
            !    enddo
            !    enddo
            !endif
        !endif

    end subroutine CenterField

    !--------------------------------------------------------------------------

    subroutine CenterFieldR8(PropCenter, PropFaces, TypeZUV, Size) !, CyclicBoundary)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :, :),   pointer          :: PropCenter, PropFaces
        integer                                         :: TypeZUV
        type (T_Size3D)                                 :: Size
        !logical                                         :: CyclicBoundary

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

        select case (TypeZUV)

            case (TypeU_)

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then

                        !Check for non state faces (assume MISSING_NULL)
                        if (PropFaces(i,j,k) < FillValueReal/2.) PropFaces(i,j,k) = 0.
                        if (PropFaces(i,j+1,k) < FillValueReal/2.) PropFaces(i,j+1,k) = 0.

                        PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
                                               PropFaces(i, j+1, k)) / 2.
                    else

                        PropCenter(i, j, k) = 0.
                    endif
                enddo
                enddo
                enddo

            case (TypeV_)

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then

                        !Check for non state faces (assume MISSING_NULL)
                        if (PropFaces(i,j,k) < FillValueReal/2.) PropFaces(i,j,k) = 0.
                        if (PropFaces(i+1,j,k) < FillValueReal/2.) PropFaces(i+1,j,k) = 0.

                        PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
                                               PropFaces(i+1, j, k)) / 2.
                    else

                        PropCenter(i, j, k) = 0.
                    endif
                enddo
                enddo
                enddo

        end select

        !if (TypeZUV == TypeU_) then

        !    if (.not. CyclicBoundary) then

        !        do k = KLB, KUB
        !        do j = JLB, JUB - 1
        !        do i = ILB, IUB
        !            if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then
               
        !                PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
        !                                       PropFaces(i, j+1, k)) / 2.
        !            else

        !                PropCenter(i, j, k) = 0.
        !            endif
        !        enddo
        !        enddo
        !        enddo

        !    else !(cyclic boundary)

        !        do k = KLB, KUB
        !        do j = JLB, JUB - 2
        !        do i = ILB, IUB
        !            if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then
               
        !                PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
        !                                       PropFaces(i, j+1, k)) / 2.
        !            else

        !                PropCenter(i, j, k) = 0.
        !            endif
        !        enddo
        !        enddo

        !        do i = ILB, IUB
        !            if (Me%External_Var%OpenPoints(i, JUB - 1, k) == OpenPoint) then
                    
        !                if (Me%External_Var%OpenPoints(i, JLB, k) == OpenPoint) then

        !                    PropCenter(i, JUB - 1, k) = (PropFaces(i, JUB - 1, k) + &
        !                                            PropFaces(i, JLB, k)) / 2.
        !                else

        !                    PropCenter(i, JUB - 1, k) = PropFaces(i, JUB - 1, k) / 2.
        !                    !(assumed null value for property at JLB face)
        !                endif
        !            else

        !                PropCenter(i, JUB - 1, k) = 0.
        !            endif
        !        enddo
        !        enddo
        !    endif

        !else !TypeV_

        !    if (.not. CyclicBoundary) then

        !        do k = KLB, KUB
        !        do j = JLB, JUB
        !        do i = ILB, IUB - 1
        !            if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then

        !                PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
        !                                       PropFaces(i+1, j, k)) / 2.
        !            else

        !                PropCenter(i, j, k) = 0.
        !            endif
        !        enddo
        !        enddo
        !        enddo

        !    else !(cyclic boundary)

        !        do k = KLB, KUB
        !        do j = JLB, JUB
        !        do i = ILB, IUB - 2
        !            if (Me%External_Var%OpenPoints(i, j, k) == OpenPoint) then
               
        !                PropCenter(i, j, k) = (PropFaces(i, j, k) +                     &
        !                                       PropFaces(i+1,j, k)) / 2.
        !            else

        !                PropCenter(i, j, k) = 0.
        !            endif
        !        enddo
        !        enddo

        !        do j = JLB, JUB
        !            if (Me%External_Var%OpenPoints(IUB - 1, j, k) == OpenPoint) then
                    
        !                if (Me%External_Var%OpenPoints(ILB, j, k) == OpenPoint) then
               
        !                    PropCenter(IUB - 1, j, k) = (PropFaces(IUB - 1, j, k) +     &
        !                                                 PropFaces(ILB, j, k)) / 2.
        !                else

        !                    PropCenter(IUB - 1, j, k) = PropFaces(IUB - 1, j, k) / 2.
        !                    !(assumed null value for property at ILB face)
        !                endif
        !            else

        !                PropCenter(IUB - 1, j, k) = 0.
        !            endif
        !        enddo
        !        enddo
        !    endif
        !endif

    end subroutine CenterFieldR8

    !--------------------------------------------------------------------------

    subroutine OutPutTimeSeries

        !Arguments-------------------------------------------------------------
        

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                            :: DepthLevel
        integer                                         :: STAT_CALL, TimeSerieNumber, dn
        integer                                         :: id, jd, kd
        logical                                         :: DepthON
        type(T_Property), pointer                       :: ObjProperty

        !----------------------------------------------------------------------

        !Corrects if necessary the cell of the time serie based in the time serie depth
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR01'  

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                          &  
                                      LocalizationI = id,                           &
                                      LocalizationJ = jd,                           &
                                      DepthLevel    = DepthLevel,                   &
                                      DepthON       = DepthON,                      & 
                                      STAT          = STAT_CALL)
            if (DepthON) then

                if (Id < 0 .or. Jd < 0)                                             &
                    stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR02'

                kd = GetLayer4Level(Me%ObjGeometry, id, jd, DepthLevel, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                          &
                    stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR03'

                if (Me%External_Var%WaterPoints(id, jd, kd) == WaterPoint) then

                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn,  k = kd,       &
                                                STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR04'
                else

                    stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR05'

                endif
            endif

        enddo

        ObjProperty => Me%FirstProperty

        do while (associated(ObjProperty))

            select case (ObjProperty%Dim)

                case(Dim_2D)

                    call WriteTimeSerie(Me%ObjTimeSerie,                            &
                                        Data2D = ObjProperty%Field2D,               &
                                        STAT   = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR06'

                case(Dim_3D)

                    if (ObjProperty%TypeZUV == TypeZ_) then

                        call WriteTimeSerie(Me%ObjTimeSerie,                        &
                                            Data3D = ObjProperty%Field,             &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                  &
                            stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR07'
                    else

                        if ((ObjProperty%ID%IDNumber /= WaterFluxX_) .and.          &
                            (ObjProperty%ID%IDNumber /= WaterFluxY_)) then

                            call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                                Data3D = ObjProperty%CenteredField, &
                                                STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                                stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR08'
                        else

                            call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                                Data3D_8 = ObjProperty%CenteredFieldR8, &
                                                STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                                stop 'OutPutTimeSeries - ModuleSequentialAssimilation - ERR09'
                        endif
                    endif

            end select

            ObjProperty => ObjProperty%Next

        enddo
  
    end subroutine OutPutTimeSeries

    !--------------------------------------------------------------------------

    subroutine EOFOutPut

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Property), pointer                   :: ObjProperty
        logical                                     :: FirstTime
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        real, dimension(:,:),       pointer         :: AuxField2D
        real, dimension(:,:,:),     pointer         :: AuxField, AuxFieldFaces
        integer                                     :: neof
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WorkILB, WorkIUB, WorkJLB
        integer                                     :: WorkJUB, WorkKLB, WorkKUB
        integer                                     :: SizeILB, SizeIUB, SizeJLB
        integer                                     :: SizeJUB, SizeKLB, SizeKUB
        integer                                     :: i, j, k
        integer                                     :: OutPutNumber
        integer, dimension(:,:,:), pointer          :: PropertyMap  => null()
        integer                                     :: CurrentStateVar
        character(len=StringLength)                 :: EOFName

        !----------------------------------------------------------------------

        !Bounds
        SizeILB = Me%Size%ILB 
        SizeIUB = Me%Size%IUB 
        SizeJLB = Me%Size%JLB 
        SizeJUB = Me%Size%JUB 
        SizeKLB = Me%Size%KLB 
        SizeKUB = Me%Size%KUB

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'EOFOutPut - ModuleSequentialAssimilation - ERR01'

        FirstTime = .true.        

        OutPutNumber = Me%EOF_OutPut%NextOutPut

TOut:   if (Me%CurrentTime >= Me%EOF_OutPut%OutTime(OutPutNumber)) then

            do neof = 1, Me%StateCovRank !Cycle the EOFs

                write(EOFName,'(i3)') neof

                ObjProperty => Me%FirstProperty

                !Write EOF as property fields
                CurrentStateVar = 0

    PropX:      do while (associated(ObjProperty))

    First:          if (FirstTime) then 

                        !Writes current time
                        call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2), &
                                            AuxTime(3), AuxTime(4), AuxTime(5),     &
                                            AuxTime(6))
                        TimePtr => AuxTime

                        !Writes Time
                        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                  &
                            stop 'EOFOutPut - ModuleSequentialAssimilation - ERR02'

                        call HDF5WriteData(Me%ObjHDF5, "/EOFs/Time", "Time",        &
                                           "YYYY/MM/DD HH:MM:SS", Array1D = TimePtr, &
                                           OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                  &
                            stop 'EOFOutPut - ModuleSequentialAssimilation - ERR03'

                        FirstTime = .false.
                    endif First

                    select case (ObjProperty%Dim)

                        case(Dim_2D)

                            !field for EOF
                            allocate(AuxField2D(WorkILB:WorkIUB,WorkJLB:WorkJUB))
                            AuxField2D(:,:) = FillValueReal

                            ILB = ObjProperty%Window2D%ILB
                            IUB = ObjProperty%Window2D%IUB
                            JLB = ObjProperty%Window2D%JLB
                            JUB = ObjProperty%Window2D%JUB

                            do j = JLB, JUB
                            do i = ILB, IUB

                                if (Me%External_Var%WaterPoints(i, j, WorkKUB)      &
                                    == WaterPoint) then
                                
                                    CurrentStateVar = CurrentStateVar + 1                            

                                    AuxField2D(i,j) = Me%LMatrix(CurrentStateVar, neof)
                            
                                endif

                            enddo
                            enddo

                            call HDF5SetLimits(Me%ObjHDF5, WorkILB,                 &
                                               WorkIUB, WorkJLB, WorkJUB,           &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                            & 
                            stop 'EOFOutPut - ModuleSequentialAssimilation - ERR04'

                            call HDF5WriteData(Me%ObjHDF5,                          &
                                "/EOFs/"//trim(ObjProperty%ID%Name)//"/"//trim(adjustl(EOFName)), &
                                               "eof", "-", Array2D = AuxField2D,    &
                                        OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                            stop 'EOFOutPut - ModuleSequentialAssimilation - ERR05'            

                            deallocate(AuxField2D)
                            nullify(AuxField2D)

                        case(Dim_3D)

                            ILB = ObjProperty%Window%ILB
                            IUB = ObjProperty%Window%IUB
                            JLB = ObjProperty%Window%JLB
                            JUB = ObjProperty%Window%JUB
                            KLB = ObjProperty%Window%KLB
                            KUB = ObjProperty%Window%KUB

                            !allocate a generic field (cell centered)
                            allocate(AuxField(WorkILB:WorkIUB, WorkJLB:WorkJUB,     &
                                              WorkKLB:WorkKUB))
                            AuxField(:,:,:) = FillValueReal

                            if (ObjProperty%TypeZUV == TypeZ_) then

                                do k = KLB, KUB
                                do j = JLB, JUB
                                do i = ILB, IUB
                                    if (Me%External_Var%WaterPoints(i, j, k) ==     &
                                        WaterPoint) then

                                        CurrentStateVar = CurrentStateVar + 1

                                        AuxField(i,j,k) =                           &
                                            Me%LMatrix(CurrentStateVar, neof)

                                    endif
                                enddo
                                enddo
                                enddo

                            else
                            
                                if (ObjProperty%TypeZUV == TypeU_) then

                                    PropertyMap => Me%WaterFacesU

                                else

                                    PropertyMap => Me%WaterFacesV

                                endif

                                allocate(AuxFieldFaces(SizeILB:SizeIUB,             &
                                                       SizeJLB:SizeJUB,             &
                                                       SizeKLB:SizeKUB))
                                AuxFieldFaces(:,:,:) = FillValueReal

                                do k = KLB, KUB
                                do j = JLB, JUB
                                do i = ILB, IUB
                                    if (PropertyMap(i, j, k) == 1) then

                                        CurrentStateVar = CurrentStateVar + 1

                                        AuxFieldFaces(i,j,k) =                      &
                                            Me%LMatrix(CurrentStateVar, neof)

                                    endif
                                enddo
                                enddo
                                enddo

                                call CenterField(AuxField, AuxFieldFaces,           &
                                                 ObjProperty%TypeZUV, ObjProperty%Window)
                                deallocate(AuxFieldFaces)
                                nullify(AuxFieldFaces)

                            endif

                            call HDF5SetLimits(Me%ObjHDF5, WorkILB,                 &
                                               WorkIUB, WorkJLB, WorkJUB,           &
                                               WorkKLB, WorkKUB, STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                            & 
                            stop 'EOFOutPut - ModuleSequentialAssimilation - ERR06'

                            call HDF5WriteData(Me%ObjHDF5,                          &
                                "/EOFs/"//trim(ObjProperty%ID%Name)//"/"//trim(adjustl(EOFName)), &
                                               "eof", "-", Array3D = AuxField,      &
                                        OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                              &
                            stop 'EOFOutPut - ModuleSequentialAssimilation - ERR07'            

                            deallocate(AuxField)
                            nullify(AuxField)

                    end select

                    !Writes everything to disk
                    call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'EOFOutPut - ModuleSequentialAssimilation - ERR08'

                    ObjProperty => ObjProperty%Next

                enddo PropX

            enddo

            Me%EOF_OutPut%NextOutPut = OutPutNumber + 1

        endif  TOut    

        nullify(ObjProperty)

    end subroutine EOFOutPut

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine KillSequentialAssimilation(SequentialAssimilationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: SequentialAssimilationID              
        integer, optional, intent(OUT)              :: STAT

        !External----------------------------------------------------------------
        integer                                     :: ready_              

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_, nUsers
        integer                                     :: STAT_CALL
        type (T_StateProp), pointer                 :: ObjWaterProperty
        type (T_ErrorTimeSerie), pointer            :: ErrorTimeSerie
        integer                                     :: imeas

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SequentialAssimilationID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSequentialAssimilation_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%OutPut%HDF5ON) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                stop 'KillSequentialAssimilation - ModuleSequentialAssimilation - ERR01'
                endif

                if (Me%OutPut%TimeSerieON) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'KillSequentialAssimilation - ModuleSequentialAssimilation - ERR02'
                endif

                if (Me%OutPut%ErrorTimeSerieON) then

                    ErrorTimeSerie => Me%FirstErrorTimeSerie
                    do imeas = 1, Me%MeasuresNumber

                        call KillTimeSerie(ErrorTimeSerie%TimeSerie, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'KillSequentialAssimilation - ModuleSequentialAssimilation - ERR03'
                        ErrorTimeSerie => ErrorTimeSerie%Next

                    enddo
                endif

                if (Me%StateCovEvolution) then
                    !Nullify aux state pointers
                    !>>> Hydrodynamic
                    call NullifyHydroStatePointer(Me%ObjHydrodynamic, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                stop 'KillSequentialAssimilation - ModuleSequentialAssimilation - ERR04'

                    !>>> WaterProperties
                    call NullifyDensityPointer(Me%ObjWaterProperties, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                stop 'KillSequentialAssimilation - ModuleSequentialAssimilation - ERR05'

                    !(only IDs are required, instead of Me%FS could be used Me%DFS)
                    ObjWaterProperty => Me%FS%FirstWaterProperty
                    do while (associated(ObjWaterProperty)) 

                        call NullifyConcentrationPointer(Me%ObjWaterProperties,         &
                                                         ObjWaterProperty%IDNumber,     &
                                                         STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                stop 'KillSequentialAssimilation - ModuleSequentialAssimilation - ERR06'

                        ObjWaterProperty => ObjWaterProperty%Next
                    enddo

                    call NullifyGeometryStatePointer(Me%ObjGeometry, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                stop 'KillSequentialAssimilation - ModuleSequentialAssimilation - ERR07'
                endif

                call Deassociate_External_Modules

                call DeallocateVariables

                !Deallocates Instance
                call DeallocateInstance ()

                SequentialAssimilationID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillSequentialAssimilation
        
    !------------------------------------------------------------------------

    Subroutine Deassociate_External_Modules

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: nUsers

        !------------------------------------------------------------------------

        nUsers = DeassociateInstance (mTIME_,           Me%ObjTime)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR01'

        nUsers = DeassociateInstance (mGRIDDATA_,       Me%ObjGridData)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR02'

        nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR03'
        
        nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR04'

        nUsers = DeassociateInstance (mMAP_,            Me%ObjMap)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR05'

        nUsers = DeassociateInstance (mHYDRODYNAMIC_,   Me%ObjHydrodynamic)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR06'

        nUsers = DeassociateInstance (mWATERPROPERTIES_, Me%ObjWaterProperties)
        if (nUsers == 0) stop 'KillHydrodynamic - ModuleHydrodynamic - ERR07'
        
    End Subroutine Deassociate_External_Modules

    !------------------------------------------------------------------------

    subroutine DeallocateVariables

        !Arguments---------------------------------------------------------------


        !External----------------------------------------------------------------


        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: fs
        type(T_Property), pointer                   :: ObjProperty
        type(T_Measure),   pointer                  :: ObjMeasure

        !------------------------------------------------------------------------

        !Deallocate map faces variables
        if (associated(Me%WaterFacesU)) then

            deallocate (Me%WaterFacesU, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR01'
            nullify (Me%WaterFacesU)

            deallocate (Me%WaterFacesV, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR02'
            nullify (Me%WaterFacesV)

            deallocate (Me%OpenFacesU, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR03'
            nullify (Me%OpenFacesU)

            deallocate (Me%OpenFacesV, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR04'
            nullify (Me%OpenFacesV)

        endif

        !Deallocate sequential assimilation filter variables
        deallocate (Me%State, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR05'
        nullify (Me%State) 

        deallocate (Me%OpenPointState, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR06'
        nullify (Me%State) 

        deallocate (Me%MeasuredState, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR07'
        nullify (Me%MeasuredState) 

        deallocate (Me%MeasuresInvCov, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR08'
        nullify (Me%MeasuresInvCov)

        deallocate (Me%ObservOperator, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR09'
        nullify (Me%ObservOperator)

        if (Me%Method == SEEK) then

            deallocate (Me%LMatrix, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR10'
            nullify (Me%LMatrix)

            deallocate (Me%InvUMatrix, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR11'
            nullify (Me%InvUMatrix)

            if (Me%StateCovEvolution) then

                !Deallocate full state variables
                call KillFullState(Me%FS)

                do fs = 1, Me%StateCovRank

                    call KillFullState(Me%DFS(fs))

                enddo

                deallocate (Me%DFS, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR12'
                nullify (Me%DFS)
            endif
        else
            !...
        endif

        !Deallocate state property fields
        ObjProperty => Me%FirstProperty

        do while (associated(ObjProperty))

            if (associated(ObjProperty%Field)) then

                deallocate (ObjProperty%Field, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR13'
                nullify (ObjProperty%Field)

            elseif (associated(ObjProperty%Field2D)) then

                deallocate (ObjProperty%Field2D, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR14'
                nullify (ObjProperty%Field2D)

            else

                deallocate (ObjProperty%FieldR8, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR15'
                nullify (ObjProperty%FieldR8)
            endif

            ObjProperty => ObjProperty%Next
        enddo

        nullify (Me%FirstProperty, Me%LastProperty)

        !Deallocate measurements
        ObjMeasure => Me%FirstMeasure

        do while (associated(ObjMeasure))

            call KillTimeSerie(ObjMeasure%TimeSerie%ID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'DeallocateVariables - ModuleSequentialAssimilation - ERR16'
            ObjMeasure => ObjMeasure%Next
        enddo

        nullify (Me%FirstMeasure, Me%LastMeasure)

    end subroutine DeallocateVariables

    !------------------------------------------------------------------------

    subroutine KillFullState(ObjFullState)

        !Arguments---------------------------------------------------------------
        type (T_FullState)                          :: ObjFullState

        !External----------------------------------------------------------------


        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_StateProp) , pointer                :: ObjFSProperty

        !------------------------------------------------------------------------

        !Deallocate full state variables
        deallocate (ObjFullState%WaterLevel, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR01'
        nullify (ObjFullState%WaterLevel)

        deallocate (ObjFullState%VelocityU, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR02'
        nullify (ObjFullState%VelocityU)

        deallocate (ObjFullState%VelocityV, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR03'
        nullify (ObjFullState%VelocityV)

        deallocate (ObjFullState%VelocityUOld, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR02'
        nullify (ObjFullState%VelocityUOld)

        deallocate (ObjFullState%VelocityVOld, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR03'
        nullify (ObjFullState%VelocityVOld)

        deallocate (ObjFullState%VelocityAcross, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR04'
        nullify (ObjFullState%VelocityAcross)

        deallocate (ObjFullState%WaterFluxX, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR05'
        nullify (ObjFullState%WaterFluxX)

        deallocate (ObjFullState%WaterFluxY, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR06'
        nullify (ObjFullState%WaterFluxY)

        deallocate (ObjFullState%WaterFluxZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR07'
        nullify (ObjFullState%WaterFluxZ)

        deallocate (ObjFullState%ChezyVelUV, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR08'
        nullify (ObjFullState%ChezyVelUV)

        if (Me%HydroSubModel) then
            deallocate (ObjFullState%SubModelqX, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'KillFullState - ModuleSequentialAssimilation - ERR09'
            nullify (ObjFullState%SubModelqX)

            deallocate (ObjFullState%SubModelqY, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'KillFullState - ModuleSequentialAssimilation - ERR10'
            nullify (ObjFullState%SubModelqY)
        endif

        deallocate (ObjFullState%Density, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR11'
        nullify (ObjFullState%Density)

        deallocate (ObjFullState%SigmaDensity, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR12'
        nullify (ObjFullState%SigmaDensity)

        if (Me%FullWaterPropNumber > 0) then
            !Deallocate water properties field
            ObjFSProperty => ObjFullState%FirstWaterProperty

            do while (associated(ObjFSProperty))

                deallocate(ObjFSProperty%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'KillFullState - ModuleSequentialAssimilation - ERR13'
                nullify   (ObjFSProperty%Field)

                ObjFSProperty => ObjFSProperty%Next
            enddo

            nullify (ObjFullState%FirstWaterProperty,ObjFullState%LastWaterProperty)
        endif

        deallocate (ObjFullState%SZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR14'
        nullify (ObjFullState%SZZ)

        deallocate (ObjFullState%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR15'
        nullify (ObjFullState%DWZ)

        deallocate (ObjFullState%DUZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR16'
        nullify (ObjFullState%DUZ)

        deallocate (ObjFullState%DVZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR17'
        nullify (ObjFullState%DVZ)

        deallocate (ObjFullState%DZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR18'
        nullify (ObjFullState%DZZ)

        deallocate (ObjFullState%ZCellCenter, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR19'
        nullify (ObjFullState%ZCellCenter)

        deallocate (ObjFullState%AreaU, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR20'
        nullify (ObjFullState%AreaU)

        deallocate (ObjFullState%AreaV, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR21'
        nullify (ObjFullState%AreaV)

        deallocate (ObjFullState%WaterColumnZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR22'
        nullify (ObjFullState%WaterColumnZ)

        deallocate (ObjFullState%WaterColumnU, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR23'
        nullify (ObjFullState%WaterColumnU)

        deallocate (ObjFullState%WaterColumnV, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR24'
        nullify (ObjFullState%WaterColumnV)

        deallocate (ObjFullState%VolumeZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR25'
        nullify (ObjFullState%VolumeZ)

        deallocate (ObjFullState%VolumeZOld, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR26'
        nullify (ObjFullState%VolumeZOld)

        deallocate (ObjFullState%VolumeU, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR27'
        nullify (ObjFullState%VolumeU)

        deallocate (ObjFullState%VolumeV, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'KillFullState - ModuleSequentialAssimilation - ERR28'
        nullify (ObjFullState%VolumeV)

    end subroutine KillFullState

    !------------------------------------------------------------------------

    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_SequentialAssimilation), pointer          :: AuxObjSeqAssimilation
        type (T_SequentialAssimilation), pointer          :: PreviousObjSeqAssimilation

        !------------------------------------------------------------------------

        !Updates pointers
        if (Me%InstanceID == FirstObjSeqAssimilation%InstanceID) then
            FirstObjSeqAssimilation => FirstObjSeqAssimilation%Next
        else
            PreviousObjSeqAssimilation => FirstObjSeqAssimilation
            AuxObjSeqAssimilation      => FirstObjSeqAssimilation%Next
            do while (AuxObjSeqAssimilation%InstanceID /= Me%InstanceID)
                PreviousObjSeqAssimilation => AuxObjSeqAssimilation
                AuxObjSeqAssimilation      => AuxObjSeqAssimilation%Next
            enddo

            !Now update linked list
            PreviousObjSeqAssimilation%Next => AuxObjSeqAssimilation%Next
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

    subroutine Ready (ObjSequentialAssimilation_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSequentialAssimilation_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjSequentialAssimilation_ID > 0) then
            call LocateObjSequentialAssimilation (ObjSequentialAssimilation_ID)
            ready_ = VerifyReadLock (mSequentialAssimilation_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjSequentialAssimilation (ObjSequentialAssimilationID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSequentialAssimilationID

        !Local-----------------------------------------------------------------

        Me => FirstObjSeqAssimilation
        do while (associated (Me))
            if (Me%InstanceID == ObjSequentialAssimilationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleSequentialAssimilation - LocateObjSequentialAssimilation - ERR01'

    end subroutine LocateObjSequentialAssimilation

    !--------------------------------------------------------------------------

end module ModuleSequentialAssimilation









