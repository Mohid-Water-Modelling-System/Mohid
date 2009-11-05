!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Lagrangian
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jun 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Three dimensional lagrangian tracers model
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

Module ModuleLagrangian

!BOP
!
! !MODULE: ModuleLagrangian

!    !DESCRIPTION: 
!     This model is responsible Create, Move and Destroy Lagrangian Tracers
 
! !REVISION HISTORY: 
!    Jul2001   Frank Braunschweig  New Implementation
!
! !FILES USED:   
!  'nomfich.dat' file is where the name files are given .
!  The keywords relevant for the hydrodynamic model are:
!    PARTIC_DATA : data file
!    PARTIC_INI  : restart file
!    PARTIC_FIN  : final output file
!    PARTIC_HDF  : transiente output file
!              
!
! !SEE ALSO:    
!  http://www.mohid.com
!


!
!DataFile
!
!   DT_PARTIC               : sec.                      [DT_Model]              !Particle Time Step    
!   OUTPUT_TIME             : sec. sec. sec.            []                      !Output Time
!   RESTART_FILE_OUTPUT_TIME: sec. sec. sec.            []                      !Output Time to write restart files
!   RESTART_FILE_OVERWRITE  : 0/1                       [1]                     !Overwrite intermediate restart files    
!   PARTIC_BOX              : char                      []                      !Particle Box definition file
!   MONITOR_BOX             : char                      []                      !Particle Monitoring box
!   MONITOR_BOX_PROP_MASS   : char                      []                      !Name of property to monitor Mass in a box
!   MONITOR_BOX_MASS_FRACTION: 1/2                      [1]                     !Mass OutPut Integration Type in box
!                                                                                1 - Arithmetic  2 - Geometric
!   MONITOR_BOX_MIN_CONC    : real                      [1]                     ! Minimum Tracer Conc. Value for Av.Concentration 
!                                                                               !and Contamination Probab. in Each Monitoring Box
!   EULERIAN_MONITOR_BOX    : char                      []                      !Eulerian results monitoring box
!   ASSOCIATE_BEACH_PROB    : 0/1                       [0]                     !Associates Beaching Probabilities
!   DEFAULT_BEACHING_PROB   : real                      [0.5]                   !Outbox Beaching Probability
!   BEACHING_LIMIT          : real                      [5.0]                   !Maximum distance between particles and coast
                                                                                !for particle beaching
!   BEACHING_BOX_FILENAME   : char                      []                      !Beaching Probability Box definition file
!   BOXES_BEACHING_PROB     : real(Boxes Number)        []                      !List of Inbox Beaching Probability
!   OUTPUT_CONC             : 1/2                       [1]                     !OutPut Integration Type
!   OUTPUT_MAX_TRACER       : 0/1                       [0]                     !Checks if the users wants to output the maximum 
                                                                                 ! tracer concentration in each cell
!   OVERLAY_VELOCITY        : 0/1                       [0]                     !If a adicional velocity field is to be added
!                                                                                1 - Maximum 2 - Mean
!
!
!<BeginOrigin>
!   ORIGIN_NAME             : char                      [Origin_xx]             !Name of the origin
!   OLD                     : 0/1                       [0]                     !Old Origin
!   GROUP_ID                : integer                   [1]                     !Group to which belong Origin
!   EMISSION_SPATIAL        : Point/Accident/Box        [-]                     !Spatial emission type
!   EMISSION_TEMPORAL       : Continuous/Instantaneous  [-]                     !Temporal emission type
!   DT_EMIT                 : sec                       [DT_PARTIC]             !Interval between continuous emissions
!   START_PARTIC_EMIT       : YYYY MM DD HH MM SS       [BeginModel]            !InitialData of the emission
!   STOP_PARTIC_EMIT        : YYYY MM DD HH MM SS       [EndModel]              !FinalData of the emission    
!   NBR_PARTIC              : int                       [1]                     !Number of Particles in each emission
!   FLOW                    : real                      [-]                     !Flow associated to point origin 
!   FLOW_VARIABLE           : 0/1                       [0]                     !Check if the user wants a variable water flow
!   DISCHARGE_FILE          : char                      [ ]                     !Name of the time serie input where is defined 
                                                                                !the variable flow
!   FLOW_COLUMN             : int                       [ ]                     !Column of the time serie input where is defined
                                                                                !a variable flow
!                                                                                (continous emission)
!   POINT_VOLUME            : real                      [-]                     !Volume of instantanous emission
!                                                                                (point or accident)
!   TVOL200                 : real                      [-]                     !Time to double volume
!   VOLUME_INCREASE         : Double/Velocity                                   !Type of function to calculate volume increase
!   VOLFAC                  : real                      [10]                    !Factor of how many times a particle can grow
!   SPLIT_PART              : 0/1                       [0]                     !Split big particles
!   KILL_LAND_PARTICLES     : 0/1                       [0]                     !Kills particles which are located in a Waterpoint
!                                                                               !which is not a OpenPoint                   
!   ACCIDENT_METHOD         : Fay/Thickness             [Fay]                   !Way to calculate initial area of accident
!   FLOAT                   : 0/1                       [0]                     !Floating particle
!       THICKNESS_METERS    : meters                    []                      !Thickness of floating particles
!                                                                               ! (just Boxes / Accident Initialization)
!   MOVEMENT                : SullivanAllen/NotRandom   []                      !Aleatoria Horizontal Movement
!       VARVELHX            :                           [0.2]                   !
!       VARVELH             :                           [0.0]                   !
!   TURB_V                  : Constant/Profile          []                      !Vertical turbulence parameterization
!       VARVELVX            :                           [0.0]
!       VARVELV             :                           [0.0]
!   ADVECTION               : 0/1                       [1]                     !Move Particle due horizontal velocity    
!   WINDCOEF                : real                      [0.03]                  !Wind transfer Coeficient
!   WINDXY                  : real real                 [0.0 0.0]               !If this keyword is defined than the the wind 
                                                                                !the wind velocity defined in the atmosphere 
                                                                                !module is override nad the wind use by the 
                                                                                !tracers is this one 
!   SEDIMENTATION           : Stokes/Imposed            []                      !Sedimentation type
!       D50                 : real  (mm)                [0.002]                 !Stokes
!       SED_VELOCITY        : real                      []                      !Sedimentation Velocity
!       MIN_SED_VELOCITY    : real                      [0.0]                   !Minimum Sedimention velocity
!   POSITION_METERS         : meters meters             []                      !X and Y Position of the origin in meters 
!   POSITION_CELLS          : Cell Cell                 []                      !X and Y Position of the origin in grid cells
!   POSITION_COORDINATES    : Xcoord Ycoord             []                      !X and Y Position of the origin in the grid        
                                                                                !coordinates
!   DEPTH_METERS            : meters (from surface)     []                      !Depth of emission relativ to surface
!   DEPTH_CELLS             : Cell                      []                      !Depth in Cells (from bottom)
!   MAINTAIN_RELATIVE_POSITION : logical                []                      !Check is the user wants to maintain 
                                                                                !the vertical relative position of the origin
!   INCRP                   : int                       [1]                     !Increment of grid cells to fill Boxes
!   BOX_NUMBER              : int                       []                      !Number of box to associate to origin
!   BOXVOLINIC              : real                      [VolumeOfBox]           !Initial Volume of a particle in the box
!   SW_EXTINCTION_COEF      : real                      [1/20]                  !Short Wave Extintion factor
!   SW_PERCENTAGE           : real                      [0.6]                   !Parcentage of shortwave radiation
!   THEORIC_AREA            : 0/1                       [0]                     !Uses Theoric Area for Oil Processes
!   MOVING_ORIGIN           : 0                         [0]                     !A moving origin (Just Point Emission)
!   MOVING_ORIGIN_FILE      : char                      []                      !FileName with Trajectory of origin
!   MOVING_ORIGIN_UNITS     : char                      [Cells]                 !Time Serie given in meters or in cells
!   MOVING_ORIGIN_COLUMN_X  : int                       []                      !Column in time serie file - X Coordinate
!   MOVING_ORIGIN_COLUMN_Y  : int                       []                      !Column in time serie file - Y Coordinate
!   COMPUTE_AGE             : logical                   [0]                     !This logical option allows to compute the 

!   COMPUTE_AVERAGE_POSITION: logical                   [0]                     !This logical option allows to compute the 
!                                                                                average position and the radius of influence of 
!                                                                                of the tracers of each origin
!   COMPUTE_PLUME           : logical                   [0]                     !Computes Particle Plume due density gradients
!   COMPUTE_BUOYANCY        : logical                   [0]                     !Computes Particle vertical velocity evolution
                                                                                !due to density gradients
!   DENSITY_METHOD          : int                       [3]                     !Way to calculate particle density 
                                                                                !(1-LeendertseState_, 
                                                                                ! 2-UNESCOState_)
!   COEF_INITIAL_MIXING     : real                      []                      !Coefficient use to control volume increase due to 
!   JET_DATA_FILE           : char                      [********.***]          !Data file where the jet properties are defined
!   JET_DT                  : real                      [600]                   !Time interval between actualizations of the jet 
!                                                                               !properties
!   OIL_BEACHING            : 0/1                       [0]                     !Beaching Process 
!                                                                               !velocities gradientes between the tracer and the 
!                                                                               !avewage flow
!   DEPOSITION              : logical                   []                      !Checks if the tracers can deposited
!       TAU_ERO             : real                      [Pa]                    !Critical shear stress of erosion
!       TAU_DEP             : real                      [Pa]                    !Critical shear stress of deposition
!       BOTTOM_DISTANCE     : real                      [m]                     !Distance from bottom below which the tracer can 
!                                                                               !sediment
!       TIME_DECAY          : real                      [s]                     !Decay time use to compute a relxation term that 
                                                                                !makes the critical shear stress of erosion tend to 
                                                                                !the average tracer erosion rate
                                                                                !of the cell where the tracer is deposited. 
                                                                                !This is use to compute 
                                                                                !the shadow effect of the large sediments over the 
                                                                                !smaller ones. 
!       BOTTOM_EMISSION     : logical                   []                      !Checks if the tracers are emited from the bottom
!       EROSION_RATE        : real                      [g/m2/s]                !Rate of tracers erosion

!   PARTITION_WATER            : logical                []                      !Checks if the tracers has two phases 
!                                                                               !(adsorbe and dissolved) in the water column.
!       PARTITION_COEF_WATER   : real                   []                      ! partition coefficent in the water column
!       PARTITION_RATE_WATER   : real                   [s-1]                   !Rate of transfer between the two phases.
!       PARTITION_COUPLE_WATER : real                   [M/L^3]                 !Concentration of the dissolved phase. 
!                                                                               !The dissolved phase is admitted with a constant 
!                                                                               !concentration

!   PARTITION_SED              : logical                []                      !Checks if the tracers has two phases 
!                                                                               !(adsorbe and dissolved) in the sediment
!       PARTITION_COEF_SED     : real                   []                      ! partition coefficent in the sediment
!       PARTITION_RATE_SED     : real                   [s-1]                   !Rate of transfer between the two phases.
!       PARTITION_COUPLE_SED   : real                   [M/L^3]                 !Concentration of the dissolved phase in the 
!                                                                               !intersticial water. The dissolved phase is 
!                                                                               !admitted with a constant concentration.
! WQM_DATA_FILE             : character                 []                      !Data File of the WQM module
!
!<<BeginProperty>>
!   NAME                    : char                      []                      !Name of the property
!   UNITS                   : char                      []                      !Units of the property
!   CONCENTRATION           : real                      []                      !Concentration of the property
!   CONC_VARIABLE           : 0/1                       [0]                     !Check if the user wants a variable concentration
!   CONC_COLUMN             : int                       [ ]                     !Column of the time serie input where is defined 
!                                                                               !a variable concentration
!   GEOMETRIC               : 0/1                       [0]                     !Check if user wants concentration 
!                                                                                based on geometric mean 
!   NOWQM                   : logical                   [0]                     ! To compute age without running moduleWQM
!   MIN_CONCENTRATION       : real                      []
!   AMBIENT_CONC            : real                      [0.0]                   !Ambient concentration
!   T90_VARIABLE            : logical                   [0]                     !Check if the user wants to compute T90 function 
!                                                                                of ambient properties: salinity,temperature,light
!   T90_VAR_METHOD          : integer                   [1]                     !Fecal decay according to Canteras et al. (1995)
!                                                       [2]                     !Fecal decay according to Chapra (1997)
!   T90                     : real                      [7200]                  !Coliform Decay rate (s)
!   STATISTICS              : logical                   [0]                     !Wheter to calculate or not the statistic
!   STATISTICS_FILE         : char                      []                      !File name with the statistics definition
!   STATISTICS_LAG          : logical                   [0]                     !Do a frequency analysis tracer by tracer.
!<<EndProperty>>
!!  parameter from the Module Oil
!
!<<EndOil>>
!<EndOrigin>

!
!
!
!
! to do 
!   - SPLIT_PART
!   - GridThickness -> One per origin or one per Group????
!

    use ModuleGlobalData
    use ModuleTriangulation,    only : ConstructTriangulation, GetNumberOfBoundaryNodes,&
                                       GetBoundaryNodes, KillTriangulation
    use ModuleTime                  
    use ModuleHDF5
    use ModuleFunctions,        only : SigmaUNESCO, SigmaLeendertse, SigmaUNESCOPressureCorrection,&
                                       InterpolateValueInTime, RodaXY,                      &
                                       ComputeT90_Chapra, ComputeT90_Canteras,              &
                                       GetDataOnlineString, SetMatrixValue, TimeToString, ChangeSuffix
    use ModuleEnterData,        only : ReadFileName, ConstructEnterData, GetData,           &
                                       ExtractBlockFromBuffer, ExtractBlockFromBlock,       &
                                       Block_Unlock, GetOutPutTime, RewindBuffer,           &
                                       GetKeywordFromLine, GetFullBufferLine,               &
                                       ReplaceFullBufferLine,                               &
                                       GetBlockSize, KillEnterData
    use ModuleWaterQuality,     only : StartWaterQuality, WaterQuality, GetDTWQM,           &
                                       GetWQPropIndex, KillWaterQuality
    use ModuleGridData,         only : GetGridData, GetMaximumValue, UngetGridData
    use ModuleTimeSerie,        only : StartTimeSerie, StartTimeSerieInput, WriteTimeSerie, &
                                       WriteTimeSerieLine, GetTimeSerieValue, KillTimeSerie
    use ModuleLightExtinction,  only : ConstructLightExtinction, ModifyLightExtinctionField,&
                                       GetLightExtinctionOptions, KillLightExtinction,      &
                                       GetShortWaveExtinctionField, UnGetLightExtinction,   &
                                       GetLongWaveExtinctionCoef, GetRadiationPercentages
    use ModuleHorizontalMap,    only : GetBoundaries, UnGetHorizontalMap
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, WriteHorizontalGrid,              &
                                       UnGetHorizontalGrid, GetGridOrigin, GetGridCoordType,&
                                       GetCoordTypeList, GetGridAngle, GetCheckDistortion,  &
                                       LocateCell, GetDefineCellsMap, GetGridLatitudeLongitude
    use ModuleAssimilation,     only : StartAssimilation, GetAssimilationField,             &
                                       UnGetAssimilation, KillAssimilation
    use ModuleGeometry,         only : GetGeometrySize, GetGeometryWaterColumn,             &
                                       GetGeometryDistances, GetGeometryVolumes,            &
                                       GetGeometryKFloor, UnGetGeometry
    use ModuleMap,              only : GetWaterPoints3D, GetLandPoints3D, GetOpenPoints3D,  &
                                       GetComputeFaces3D, UngetMap             
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes, UngetBoxDif,&
                                       BoxDif, KillBoxDif        
    use ModuleTurbulence,       only : GetMixingLengthVertical, GetMixingLengthHorizontal,  &
                                       UngetTurbulence       
    use ModuleHydrodynamic,     only : StartHydrodynamic,                                   &
                                       KillHydrodynamic, GetHorizontalVelocity,             &
                                       GetVerticalVelocity, UngetHydrodynamic
    use ModuleWaterProperties,  only : Construct_WaterProperties, WaterPropertyExists,      &
                                       GetConcentration, GetDensity, GetSigma,              &
                                       UngetWaterProperties, GetSPM, KillWaterProperties,   &
                                       GetWaterPropertiesSubModulesID, GetDensityOptions,   &
                                       GetFiltrationRate, SetLagrangianSinksSources
    use ModuleStatistic,        only : ConstructStatistic, ModifyStatistic, KillStatistic,  &
                                       GetStatisticClassesNumber, GetStatisticClasses,      &
                                       UnGetStatistic
    use ModuleOil
    use ModuleJet,              only : Construct_Jet, GetPlumeTemperature, GetOutPutMatrix, &
                                       GetPlumeLocation, GetPlumeVelocity, GetPlumeDilution,&
                                       GetPlumeDensity, GetPlumeMixingHorLength,            &
                                       GetPlumeThickness, GetPlumeSalinity, ModifyJet,      &
                                       UnGetJet, KillJet
#ifndef _WAVES_
    use ModuleWaves
#endif

!#ifdef _CGI_
!    use dflib
!#endif

    implicit none 

    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructLagrangian
    private ::      ConstructGlobalVariables
    private ::      ConstructParticleGrid
    private ::      ConstructOrigins
    private ::          ConstructEmissionType
    private ::          ConstructParticOil
    private ::      VerifyOriginProperties
    private ::          VerifyPropertyList
    private ::      VerifyBeachingProbabilities
    private ::      ReadFinalPartic
    private ::      MergeOldWithNewOrigins
    private ::      ConstructEmission
    private ::          ConstructEmissionTime
    private ::      ConstructTimeSeries
    private ::      ConstructHDF5Output
    private ::      ConstructLog
    private ::      ConstructParticStatistic

    private ::          AllocateNewOrigin
    private ::          InsertOriginToList
    private ::          DeleteOrigin

    private ::          AllocateNewProperty
    private ::          InsertPropertyToList
    private ::          DeleteProperty

    private ::          AllocateNewParticle
    private ::          InsertParticleToList
    private ::          ActualizeJetProperties
    private ::          GiveJetPropertiesToParticle
    private ::          DeleteParticle

    private ::      ConstructParticLightExtinction


    !Modifier
    public  :: ModifyLagrangian
    private ::      ParticleEmission
    private ::          EmissionBox
    private ::          EmissionPoint
    private ::          EmissionAccident
    private ::      FillGridThickness
    private ::      FillGridConcentration
    private ::      ParticleDensity
    private ::      VerifyParticleBeaching
    private ::      MovePartic
    private ::          MoveParticHorizontal
    private ::              MoveParticVertical
    private ::                  WD_
    private ::          VerifyRessuspension
    private ::          VerifyDeposition
    private ::      VerifyLandParticles
    private ::      PurgeParticles
    private ::      VolumeVariation
    private ::      Dilution
    private ::          GetAmbientConcentration
    private ::      PropertiesEvolution
    private ::      ColiformDecay
    private ::          ComputeT90
    private ::      PartitionDecay
    private ::      ComputeAreaVolume
    private ::      InternalParticOil
    private ::          OilGridConcentration
    private ::      NewParticleMass
    private ::      NewParticleAge
    private ::      MonitorParticle
    private ::      ModifyParticStatistic
    private ::          ComputeStatisticsLag
    private ::          ActualizesTauErosionGrid

    private ::      LightEvolution

    private ::      ParticleOutput
    private ::          WriteGridConcentration
    private ::          WriteOilGridThickness
    private ::          WriteOilGridConcentration

    private ::      OutputRestartFile


    private :: Convert_XY_CellIJ
    private :: Convert_CellIJ_XY
    private :: Convert_CellIJ_IJ
    private :: Convert_Z_CellK
    private :: Convert_CellK_Z
    private :: Convert_CellK_K

    private :: Search_Property

    !Selector
    public  :: SetLagrangianShear
    public  :: SetLagrangianAtmPressure
    public  :: SetLagrangianWind
    public  :: GetLagrangianAirOptions
    public  :: SetLagrangianSolarRadiation

    !Destructor
    public  :: KillLagrangian
    private ::      WriteFinalPartic
    private ::      DeallocateOriginList
    private ::      KillParticleStatistic
    private ::          WriteFrequencyLag
    private ::      KillLight

    !Managment
    private :: ReadLockExternalVar
    private :: ReadLockEulerianDensity
    private :: ReadUnLockExternalVar
    private :: ReadUnLockEulerianDensity
    private :: Ready


    !Parameter
    character(LEN = StringLength), parameter    :: block_begin          = '<BeginOrigin>'
    character(LEN = StringLength), parameter    :: block_end            = '<EndOrigin>'
    character(LEN = StringLength), parameter    :: property_begin       = '<<BeginProperty>>'
    character(LEN = StringLength), parameter    :: property_end         = '<<EndProperty>>'
    character(LEN = StringLength), parameter    :: oil_begin            = '<<BeginOil>>'
    character(LEN = StringLength), parameter    :: oil_end              = '<<EndOil>>'

    character(LEN = StringLength), parameter    :: block_begin_clone    = '<BeginOrigin_Clone>'
    character(LEN = StringLength), parameter    :: block_end_clone      = '<EndOrigin_Clone>'


    character(LEN = StringLength), parameter    :: Char_Continuous      = 'Continuous'
    character(LEN = StringLength), parameter    :: Char_Instantaneous   = 'Instantaneous'
    character(LEN = StringLength), parameter    :: Char_Point           = 'Point'
    character(LEN = StringLength), parameter    :: Char_Accident        = 'Accident'
    character(LEN = StringLength), parameter    :: Char_Box             = 'Box'

    character(LEN = StringLength), parameter    :: Char_SullivanAllen   = 'SullivanAllen'
    character(LEN = StringLength), parameter    :: Char_NotRandom       = 'NotRandom'
    character(LEN = StringLength), parameter    :: Char_Profile         = 'Profile'
    character(LEN = StringLength), parameter    :: Char_Constant        = 'Constant'

    character(LEN = StringLength), parameter    :: Char_Stokes          = 'Stokes'
    character(LEN = StringLength), parameter    :: Char_Imposed         = 'Imposed'

    character(LEN = StringLength), parameter    :: Char_Fay             = 'Fay'
    character(LEN = StringLength), parameter    :: Char_Thickness       = 'Thickness'

    character(LEN = StringLength), parameter    :: Char_Double          = 'Double'
    character(LEN = StringLength), parameter    :: Char_Velocity        = 'Velocity'

    character(LEN = StringLength), parameter    :: Char_Cells           = 'Cells'
    character(LEN = StringLength), parameter    :: Char_Meters          = 'Meters'

    integer, parameter                          :: Continuous_          = 1
    integer, parameter                          :: Instantaneous_       = 2
    integer, parameter                          :: Point_               = 1
    integer, parameter                          :: Accident_            = 2
    integer, parameter                          :: Box_                 = 3

    !Aleat movement
    integer, parameter                          :: SullivanAllen_       = 1
    integer, parameter                          :: NotRandom_           = 2

    !Standard deviation
    integer, parameter                          :: VerticalTurbConstant = 1
    integer, parameter                          :: VerticalTurb         = 2

    !Sedimentation
    integer, parameter                          :: Stokes_              = 1
    integer, parameter                          :: Imposed_             = 2

    !accident
    integer, parameter                          :: Fay_                 = 1
    integer, parameter                          :: Thickness_           = 2

    !GridType
    integer, parameter                          :: Grid1D               = 1
    integer, parameter                          :: Grid2D               = 2

    !TVolType
    integer, parameter                          :: Velocity_            = 1
    integer, parameter                          :: Double_              = 2

    !OutputIntegrationType
    integer, parameter                          :: Maximum              = 1
    integer, parameter                          :: Mean                 = 2

    !Relative position 
    integer, parameter                          :: Cells                = 1
    integer, parameter                          :: Meters               = 2

    !T90 Calc Method
    integer, parameter                          :: Canteras             = 1
    integer, parameter                          :: Chapra               = 2


    !Online Emission options 
    integer, parameter                          :: ParticleOne          = 1
    integer, parameter                          :: Particle100          = 2
    
    ! Monitorization
    integer, parameter                      :: Arithmetic           = 1
    integer, parameter                      :: Geometric            = 2


    !Internal states
    type T_State
        logical                                 :: ContCalc             = OFF
        logical                                 :: WQM                  = OFF
        logical                                 :: Larvae               = OFF
        logical                                 :: Wind                 = OFF
        logical                                 :: Sedimentation        = OFF
        logical                                 :: Deposition           = OFF
        logical                                 :: VariableGeom         = OFF
        logical                                 :: ShearVel             = OFF
        logical                                 :: Oil                  = OFF
        logical                                 :: FCF                  = OFF
        logical                                 :: Monitor              = OFF
        logical                                 :: MonitorPropMass      = OFF
        logical                                 :: EulerianMonitor      = OFF
        logical                                 :: Partition            = OFF
        logical                                 :: AssociateBeachProb   = OFF
        logical                                 :: HaveBeachingProbBox  = OFF
        logical                                 :: Statistics           = OFF
        logical                                 :: Age                  = OFF
        logical                                 :: ComputePlume         = OFF
        logical                                 :: PlumeShear           = OFF
        logical                                 :: FarFieldBuoyancy     = OFF
        logical                                 :: T90Variable          = OFF
        logical                                 :: Filtration           = OFF
    end type T_State

    !IO
    type T_Files
        character(PathLength)                   :: ConstructData
        character(PathLength)                   :: Initial
        character(PathLength)                   :: Final
        character(PathLength)                   :: TransientHDF
        character(PathLength)                   :: BoxDataFile
        character(PathLength)                   :: MonitorBox
        character(PathLength)                   :: EulerianMonitorBox
        character(PathLength)                   :: BeachingBoxFileName    
    end type T_Files


!#ifdef _CGI_
    type T_Online
        real,  dimension(:,:), pointer          :: StartDate
        real,  dimension(:),   pointer          :: X, Y 
        real,  dimension(:),   pointer          :: WindCoef
        character(Len=23)                       :: TimeStamp
    end type T_Online
!#endif        

    !Output
    type       T_OutPut
         type (T_Time), pointer, dimension(:)   :: OutTime, RestartOutTime
         integer                                :: NextOutPut, NextRestartOutPut
         integer                                :: OutPutConcType = Maximum
         logical                                :: ConcMaxTracer   
         logical                                :: OriginEnvelope
         logical                                :: Write_, WriteRestartFile = .false. 
         logical                                :: RestartOverwrite
    end type T_OutPut

    type T_OverLay
        logical                                 :: Overlay
        real, dimension(:, :, :), pointer       :: VelUFinal
        real, dimension(:, :, :), pointer       :: VelVFinal
    end type T_OverLay

    type T_Monitorization
        real(8), dimension(:),    pointer       :: SurfaceBoxVolume
        real(8), dimension(:),    pointer       :: InstBoxVolume
        real(8), dimension(:),    pointer       :: InstBoxMass
        real(8), dimension(:, :), pointer       :: InstMassByOrigin
        real(8), dimension(:, :), pointer       :: InstVolumeByOrigin
        real(8), dimension(:),    pointer       :: IntgBoxVolume
        real(8), dimension(:, :), pointer       :: IntgVolumeByOrigin
        integer, dimension(:),    pointer       :: NumberOfCellsPerBox
        integer, dimension(:, :), pointer       :: NumberOfCellsFromOrigin
        integer, dimension(:), pointer          :: ObjTimeSerie
        character(StringLength)                 :: MassProperty
        integer                                 :: MassPropertyID           = null_int

        real(8), dimension(:),    pointer       :: InstBoxLogMass
        real(8), dimension(:),    pointer       :: InstBoxConc
        integer, dimension(:), pointer          :: NumberOfTracers
        real(8), dimension(:, :), pointer       :: InstBoxMassFractionByOrigin
        real(8), dimension(:, :), pointer       :: InstLogMassByOrigin
        integer, dimension(:, :), pointer       :: NumberOfTracersFromOrigin
        integer                                 :: MassFractionType         = Arithmetic
        integer                                 :: EulerianMonitorBoxType   = Arithmetic
        real                                    :: MonitorBox_TracerMinConc = null_real
        real                                    :: ContaminationDepth       = null_real
        real(8), dimension(:), pointer          :: ContaminationProbability 
        real(8), dimension(:), pointer          :: AverageBoxContaminatedConc  
        integer, dimension(:), pointer          :: NbrBoxContaminatedTracers 
        real(8), dimension(:), pointer          :: VolBoxContaminatedTracers 

    end type T_Monitorization
    
    type T_EulerianMonitor
        real(8), dimension(:, :, :), pointer    :: Mass
    end type T_EulerianMonitor

    type T_Filtration
        real,    dimension(:, :, :), pointer    :: RelativeMassFilter, MassFiltered
        type (T_Time)                           :: Next
    end type T_Filtration

    !Defines a generic Position
    type T_Position
        integer                                 :: I                        = null_int
        integer                                 :: J                        = null_int
        integer                                 :: K                        = null_int
        real                                    :: X                        = null_real
        real                                    :: Y                        = null_real
        real                                    :: Z                        = null_real
        real                                    :: CellI                    = null_real
        real                                    :: CellJ                    = null_real
        real                                    :: CellK                    = null_real
        integer                                 :: DepthDefinition          = null_int
        logical                                 :: MaintainRelative         = .false.
        real                                    :: Depth                    = null_real
    end type T_Position

    !Defines the movement of a particle
    type T_Movement

        !Horizontal
        integer                                 :: MovType                  = NotRandom_
        real                                    :: VarVelHX                 = 0.2
        real                                    :: VarVelH                  = 0.0

        !Vertical
        integer                                 :: StandardDeviationType    = VerticalTurbConstant
        real                                    :: VarVelVX                 = 0.0
        real                                    :: VarVelV                  = 0.0

        real                                    :: WindTransferCoef         = 0.03
        logical                                 :: WindOriginON             = .false.
        real                                    :: WindX, WindY

        !Sediment stuff
        integer                                 :: SedimentationType        = null_int
        real                                    :: SedVel                   = null_real
        real                                    :: MinSedVel                = null_real
        real                                    :: D50                      = null_real
        logical                                 :: Sediment                 =.FALSE.

        real                                    :: ThicknessMeters          = null_real

        !movement options
        logical                                 :: Float                    = OFF !Floating particle, does not mix with water
        logical                                 :: Advection                = ON
        logical                                 :: KillLandParticles        = OFF
        
        !Volumes / Areas
        real                                    :: TVOL200                  = null_real  !Time to double volume.
        integer                                 :: TVolType                 = Double_
        logical                                 :: SPLIT_PART               = OFF
        real                                    :: VOLFAC                   = null_real  
                                                  !VOLFAC * (Inicial Particle Volume) -> Particle die.

        !Plume definition
        integer                                 :: DensityMethod            = null_int
        logical                                 :: CorrecPress              = OFF
        !real                                    :: ParticleSigmaDensity     = null_real
        real                                    :: CoefInitialMixing        = null_real
        real                                    :: InitialVelocityU         = null_real
        real                                    :: InitialVelocityV         = null_real
        real                                    :: InitialVelocityW         = null_real
        real                                    :: JetDT                    = null_real
        real                                    :: JetSalinity              = null_real
        real                                    :: JetTemperature           = null_real
        real                                    :: PlumeDilution            = null_real
        type (T_Time)                           :: NextJetActualization
        character(PathLength)                   :: JetDataFile
        character(PathLength)                   :: JetFileOut
        integer                                 :: JetUnit
        integer                                 :: ObjJet                   = 0

    end type T_Movement

    !Defines the particle Geometry
    type T_ParticleGeometry
        real                                    :: VolumeOil                = null_real
        real                                    :: Volume                   = null_real
        real                                    :: InitialVolume            = null_real
        real                                    :: VolVar                   = null_real
    end type T_ParticleGeometry

    !Defines the parameters necessary to compute the fluxes associated to a 
    !partition coefficient
    type T_Partition
        logical                                 :: ON                       = OFF
        real                                    :: Coefficient              = null_real 
        real                                    :: TransferRate             = null_real 
        !By default is assumed that the particle property is particulated and is 
        !transfering mass to a dissolvde property that in the environment was an
        !aproximated constant value
        real                                    :: CoupleProp               = null_real 
    end type T_Partition


    !One Property
    type  T_Property
        integer                                 :: ID                       = null_int
        character(StringLength)                 :: Name
        character(StringLength)                 :: units
        real                                    :: concentration            = null_real 
        real                                    :: min_concentration        = 0.0       !Used to cut values in GridConcentration
        real                                    :: AmbientConcentration     = null_real
        logical                                 :: HaveAmbientConcentration = OFF
        logical                                 :: T90Variable              = OFF
        integer                                 :: T90Var_Method            = null_int
        real                                    :: T90                      = 7200.     !Coliform Bacteria decay time
        type(T_Partition)                       :: SedimentPartition
        type(T_Partition)                       :: WaterPartition
        logical                                 :: Statistics               = .false.
        character(len=PathLength)               :: StatisticsFile
        logical                                 :: StatisticsLag
        real, dimension(:,:,:,:), pointer       :: FrequencyLag
        integer                                 :: nClassesLag
        real, dimension(:,:    ), pointer       :: ClassesLag
        integer                                 :: Statistic1_ID            = 0
        integer                                 :: Statistic2_ID            = 0
        logical                                 :: ConcVariable
        integer                                 :: ConcColumn
        logical                                 :: NoWQM                    =.false.
        type (T_Property), pointer              :: Next
        real                                    :: ExtinctionParameter
        logical                                 :: Filtration               = OFF
        integer                                 :: Geometric                = null_int
    end type T_Property


    !Particle List
    type T_Partic
        integer                                 :: ID
        type (T_Position)                       :: Position
        type (T_ParticleGeometry)               :: Geometry
        real, dimension(:), pointer             :: Concentration
        real, dimension(:), pointer             :: Mass
        real                                    :: TpercursoX               = -null_real
        real                                    :: TpercursoY               = -null_real
        real                                    :: TpercursoZ               = -null_real
        real                                    :: UD_old                   = null_real
        real                                    :: VD_old                   = null_real
        real                                    :: WD_old                   = null_real
        logical                                 :: KillPartic               = OFF
        logical                                 :: Freazed                  = OFF
        logical                                 :: Beached                  = OFF
        logical                                 :: Deposited                = OFF
        real                                    :: TauErosion               = null_real
        real                                    :: ErosionRateProbability   = null_real
        real                                    :: U, V, W                  = null_real
        real                                    :: RelU, RelV               = null_real
        real                                    :: SDU, SDV                 = null_real
        real                                    :: SigmaDensity             = null_real
        type (T_Partic), pointer                :: Next                     => null()
        type (T_Partic), pointer                :: Prev                     => null()
        real                                    :: Age                      = 0.
        real                                    :: Radiation                = 0.
        real                                    :: ShortWaveExt             = null_real
        real                                    :: T90                      = null_real
    end type T_Partic

    !Particle deposition
    type  T_Deposition
        real                                    :: BottomDistance           = null_real
        real                                    :: TauErosion               = null_real
        real                                    :: ErosionRate              = null_real
        real                                    :: TauDeposition            = null_real
        real                                    :: Tdecay                   = null_real
        logical                                 :: BottomEmission           = OFF
    end type T_Deposition


    !Origin list
    type T_Origin
        character(StringLength)                 :: Name
        integer                                 :: ID
        integer                                 :: GroupID                  = 1
        logical                                 :: Old                      = OFF
        type (T_State)                          :: State
        integer                                 :: EmissionSpatial    
        integer                                 :: EmissionTemporal
        real                                    :: DT_Emit
        real                                    :: Flow
        integer                                 :: NbrParticlesIteration    = 1
        integer                                 :: AccidentMethod
        type (T_Time)                           :: AccidentTime
        logical                                 :: AccidentFinished         = .false.
        real                                    :: PointVolume
        logical                                 :: FlowVariable
        integer                                 :: FlowColumn
        character(PathLength)                   :: DischargeFile
        integer                                 :: TimeSerieInput           = 0
        type (T_Time)                           :: StartEmission
        type (T_Time)                           :: StopEmission
        type (T_Time)                           :: NextEmission
        integer                                 :: INCRP                    = 1
        integer                                 :: BoxNumber                = 1
        real                                    :: ParticleBoxVolume
        type (T_Movement)                       :: Movement
        type (T_Deposition)                     :: Deposition
        type (T_Position)                       :: Position
        integer                                 :: nParticle
        integer                                 :: nProperties
        integer                                 :: NextParticID
        real                                    :: Photoinhibition
        real                                    :: ShortWavePercentage
        logical                                 :: UseTheoricArea
        real                                    :: AreaTotal
        real                                    :: VolumeTotal
        real                                    :: VolumeOilTotal
        real                                    :: VolTotOilBeached
        real                                    :: VolTotBeached
        logical                                 :: AveragePositionON        = .false. 
        real                                    :: CoefRadius
        real, dimension(:, :), pointer          :: GridThickness
        real, dimension(:, :), pointer          :: OilGridConcentration
        logical                                 :: MovingOrigin
        character(StringLength)                 :: MovingOriginUnits
        integer                                 :: MovingOriginColumnX
        integer                                 :: MovingOriginColumnY
        character(PathLength)                   :: MovingOriginFile
        character(PathLength)                   :: WQM_DataFile
        integer                                 :: ObjTimeSerie   = 0
        integer                                 :: WaterQualityID = 0
        type (T_Time)                           :: NextWQMCompute
        real                                    :: DTWQM
        integer                                 :: ObjOil = 0
        logical                                 :: ComputeAge
        type (T_Property), pointer              :: FirstProperty  => null()
        type (T_Partic), pointer                :: FirstPartic    => null()
        type (T_Origin), pointer                :: Next           => null()
        logical                                 :: Beaching       = OFF
        logical                                 :: Filtration     = OFF
    end type T_Origin

    type T_ParticleGrid
        integer                                 :: GridType     = null_int
        real, dimension(:),    pointer          :: ParticX
        real, dimension(:),    pointer          :: ParticY
        real, dimension(:, :), pointer          :: ParticXX
        real, dimension(:, :), pointer          :: ParticYY
        real, dimension(:, :), pointer          :: AreaCell
    end type T_ParticleGrid

    type T_Light
        integer, dimension(:      ), pointer    :: ObjLightExtinction
        real,    dimension(:,:,:,:), pointer    :: TopRadiationCells
        logical                                 :: Compute      = OFF
    end type T_Light

    !ExternalVar
    type T_ExternalVar

        !Time
        type(T_Time)                            :: Now
        type(T_Time     )                       :: BeginTime
        type(T_Time     )                       :: EndTime

        real                                    :: RunPeriod

        !ObjBathymetry
        real,    dimension(:, : ), pointer      :: Bathymetry
        real,    dimension(:    ), pointer      :: XX
        real,    dimension(:    ), pointer      :: YY

        real,    dimension(:, : ), pointer      :: XX_IE
        real,    dimension(:, : ), pointer      :: YY_IE

        !ObjHorizontalMObj
        integer, pointer, dimension(:,:  )      :: BoundaryPoints2D

  
        !ObjGeometry
        type (T_Size3D)                         :: Size
        type (T_Size3D)                         :: WorkSize
        real,    pointer, dimension(:,:  )      :: WaterColumn
        real,    pointer, dimension(:,:,:)      :: SZZ, DWZ, ZCellCenter, DWZ_Xgrad, DWZ_Ygrad
        real(8), pointer, dimension(:,:,:)      :: VolumeZ
        integer, pointer, dimension(:,:  )      :: kFloor

        !ObjMap
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D
        integer, pointer, dimension(:,:,:)      :: LandPoints3D
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D
        integer, pointer, dimension(:,:,:)      :: ComputeFaces3D_U, ComputeFaces3D_V

        !ObjTurbulence
        real,    pointer, dimension(:,:,:)      :: Lupward
        real,    pointer, dimension(:,:,:)      :: Ldownward
        real,    pointer, dimension(:,:,:)      :: MixingLengthX    
        real,    pointer, dimension(:,:,:)      :: MixingLengthY 

        !ObjHydrodynamic
        real,    pointer, dimension(:,:,:)      :: Velocity_U
        real,    pointer, dimension(:,:,:)      :: Velocity_V
        real,    pointer, dimension(:,:,:)      :: Velocity_W

        !ObjInterfaceWaterAir
        real,    pointer, dimension(:,:  )      :: WindX
        real,    pointer, dimension(:,:  )      :: WindY
        real,    pointer, dimension(:,:  )      :: SurfaceRadiation
        real,    pointer, dimension(:,:  )      :: AtmPressure
        real,    pointer, dimension(:,:  )      :: WaveHeight, WavePeriod
 

        !ObjInterfaceSedimentWater
        real,    pointer, dimension(:,:  )      :: ShearVelocity
        real,    pointer, dimension(:,:  )      :: BottomStress


        !ObjOil
        real,    pointer, dimension(:,:  )      :: SpreadingVelocityX
        real,    pointer, dimension(:,:  )      :: SpreadingVelocityY
        real                                    :: DiffVelocity
        real                                    :: VWaterContent
        real                                    :: AreaTotal
        real                                    :: OilDensity
        real                                    :: MDispersed
        integer                                 :: ThicknessGradient, Fay, SpreadingMethod

        !ObjAssimilation
        real,    pointer, dimension(:,:,:)      :: OverLayU
        real,    pointer, dimension(:,:,:)      :: OverLayV

        !ObjWaterProperties
        real,    pointer, dimension(:,:,:)      :: Density
        real,    pointer, dimension(:,:,:)      :: SigmaDensity
    end type T_ExternalVar


    type T_Lagrangian
        integer                                 :: InstanceID
        type (T_State)                          :: State
        type (T_Files)                          :: Files
        type (T_OutPut)                         :: Output
        integer                                 :: nOrigins
        integer                                 :: nOldOrigins
        integer                                 :: nGroups
        logical                                 :: RunOnline
        integer, dimension(:    ),  pointer     :: GroupIDs
        integer, dimension(:    ),  pointer     :: nOriginsGroup
        type (T_Light)                          :: Light
        real,    dimension(:,:,:),  pointer     :: TauErosionGrid
        real,    dimension(:,:,:),  pointer     :: MassSedGrid
        real                                    :: DT_Partic
        real                                    :: BeachingLimit
        real                                    :: DefaultBeachingProbability
        real   , dimension(:,:,:),    pointer   :: BeachingProbability
        type(T_OverLay        )                 :: Overlay
        type(T_Monitorization )                 :: Monitor
        type(T_EulerianMonitor)                 :: EulerianMonitor
        type(T_Filtration     )                 :: Filtration
        type(T_Time)                            :: NextCompute
        type(T_Time)                            :: Now
        type(T_ParticleGrid)                    :: Grid
        type(T_Origin     ), pointer            :: FirstOrigin          => null()
        type(T_Origin     ), pointer            :: FirstOldOrigin       => null()
        type(T_ExternalVar)                     :: ExternalVar

!#ifdef _CGI_
        type(T_Online)                             :: Online    
!#endif      


        logical                                 :: WritesTimeSerie
        logical                                 :: RunOnlyMov2D         = .false.
        integer                                 :: ObjTimeSerie         = 0
        integer                                 :: ObjWaterProperties   = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: ObjGridData          = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjHorizontalMap     = 0 
        integer                                 :: ObjGeometry          = 0 
        integer                                 :: ObjMap               = 0 
        integer                                 :: ObjHydrodynamic      = 0 
        integer                                 :: ObjTurbulence        = 0
        integer                                 :: ObjWaves             = 0
        integer                                 :: ObjBoxDif            = 0
        integer                                 :: ObjMonBox            = 0
        integer                                 :: ObjEulMonBox         = 0
        integer                                 :: ObjBeachingProbBox   = 0
        integer                                 :: ObjAssimilation      = 0

        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjEnterDataClone    = 0
        integer                                 :: ObjEnterDataOriginal = 0

        type(T_Lagrangian     ), pointer        :: Next                 => null()
    end type T_Lagrangian

    !Global Module Variables
    type (T_Lagrangian), pointer                :: FirstLagrangian      => null()
    type (T_Lagrangian), pointer                :: Me                   => null()


    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine ConstructLagrangian(LagrangianID,                                         &
                                   TimeID,                                               &
                                   GridDataID,                                           &
                                   HorizontalGridID,                                     &
                                   HorizontalMapID,                                      &
                                   GeometryID,                                           &
                                   MapID,                                                &
                                   AssimilationID,                                       &
                                   HydrodynamicID,                                       &
                                   TurbulenceID,                                         &
                                   WavesID,                                              &
                                   WaterPropertiesID,                                    &
                                   STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        integer                                     :: TimeID
        integer                                     :: GridDataID
        integer                                     :: HorizontalGridID
        integer                                     :: HorizontalMapID
        integer                                     :: GeometryID
        integer                                     :: MapID
        integer                                     :: AssimilationID
        integer                                     :: HydrodynamicID
        integer                                     :: TurbulenceID
        integer                                     :: WavesID
        integer                                     :: WaterPropertiesID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL


        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mLagrangian_)) then
            nullify (FirstLagrangian)
            call RegisterModule (mLagrangian_) 
        endif

        call Ready(LagrangianID, ready_)

        if (ready_ .EQ. OFF_ERR_) then

            !Allocates a new Instance
            call AllocateInstance
            nullify (Me%FirstOrigin)

            !External Modules
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID           )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID       )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID  )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID       )
            Me%ObjMap            = AssociateInstance (mMAP_,            MapID            )
            Me%ObjHydrodynamic   = AssociateInstance (mHYDRODYNAMIC_,   HydrodynamicID   )
            Me%ObjTurbulence     = AssociateInstance (mTURBULENCE_,     TurbulenceID     )
            Me%ObjWaterProperties= AssociateInstance (mWATERPROPERTIES_,WaterPropertiesID)

            if(WavesID /= 0)then
                Me%ObjWaves      = AssociateInstance (mWAVES_,          WavesID)
            end if


            !Gets Time
            call GetComputeTimeLimits(Me%ObjTime,                             &
                                      BeginTime = Me%ExternalVar%BeginTime,   &
                                      EndTime   = Me%ExternalVar%EndTime,     &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangian - ModuleLagrangian - ERR02'

            ! Actualized the time
            call GetComputeCurrentTime(Me%ObjTime,                            &
                                       Me%ExternalVar%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangian - ModuleLagrangian - ERR03'


            !Gets Pointer to External modules
            call ReadLockExternalVar

            call ReadLockEulerianDensity

            ! Construct the variable common to all module  
            call ConstructGlobalVariables

            ! Constructs the Particle Grid
            call ConstructParticleGrid

            !Construct enter data 
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangian - ModuleLagrangian - ERR13'

            !Constructs the origin list 
            call ConstructOrigins

            if (Me%OverLay%Overlay) then

                if (AssimilationID == 0) then

                    !Associated the Assimilation
                    call StartAssimilation (Me%ObjAssimilation,                       &
                                            Me%ObjTime,                               &
                                            Me%ObjGridData,                           &
                                            Me%ObjHorizontalGrid,                     &
                                            Me%ObjHorizontalMap,                      &
                                            Me%ObjMap,                                &
                                            Me%ObjGeometry,                           &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangian - ModuleLagrangian - ERR08a'
                    
                    AssimilationID = Me%ObjAssimilation

                else
                    
                    Me%ObjAssimilation = AssociateInstance (mASSIMILATION_, AssimilationID)

                endif
            endif

            !Verifies the Origin Lists of all Origins
            call VerifyOriginProperties   
            
            !Verifies the Oil-Beaching Probabilities
            if (Me%State%AssociateBeachProb) then
                call VerifyBeachingProbabilities
            end if
 

            !Reads previuos data
            if (Me%State%ContCalc) then
                call ReadFinalPartic 
            end if

            !Starts the Light Extinction
            if (Me%State%WQM .or. Me%State%T90Variable) then
                call ConstructParticLightExtinction
            endif

            !Constructs the Time Series
            if (Me%WritesTimeSerie) then
                call ConstructTimeSeries  
            endif



            !Kills EnterData
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangian - ModuleLagrangian - ERR14'

            !Merges Old with New Origins
            call MergeOldWithNewOrigins

            !Moves Moving origins
            call MoveOrigin               

            !Constructs Emission for every origin
            call ConstructEmission        

            !Constructs the overlay velocity for every origin
            call ConstructOverlay         

            !Constructs the Monitoring 
            call ConstructMonitoring      

            !Starts the HDF Output
            if (Me%OutPut%Write_) call ConstructHDF5Output      

            !Starts the Statistic
            if (Me%State%Statistics) then
                call NewParticleMass            
                call ConstructParticStatistic   
            endif

            !Message to the User
            call ConstructLog             

            !Frees Pointer to External modules
            call ReadUnLockExternalVar   

            call ReadUnLockEulerianDensity

            !Returns ID
            LagrangianID    = Me%InstanceID
            STAT_           = SUCCESS_

        else

            stop 'ConstructLagrangian - ModuleLagrangian - ERR99' 

        endif


        if (present(STAT)) STAT = STAT_


    end subroutine ConstructLagrangian

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Lagrangian), pointer                :: NewLagrangian
        type (T_Lagrangian), pointer                :: PreviousLagrangian


        !Allocates new instance
        allocate (NewLagrangian)
        nullify  (NewLagrangian%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstLagrangian)) then
            FirstLagrangian    => NewLagrangian
            Me                 => NewLagrangian
        else
            PreviousLagrangian => FirstLagrangian
            Me                 => FirstLagrangian%Next
            do while (associated(Me))
                PreviousLagrangian  => Me
                Me                  => Me%Next
            enddo
            Me                      => NewLagrangian
            PreviousLagrangian%Next => NewLagrangian
        endif

        Me%InstanceID = RegisterNewInstance (mLAGRANGIAN_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(len = StringLength)               :: Message

        !----------------------------------------------------------------------

        !Initialize the origin ID   
        Me%nOrigins    = 0
        Me%nOldOrigins = 0

        !Initialize the origin list   
        nullify (Me%FirstOrigin   )
        nullify (Me%FirstOldOrigin)
 
        !Gets Size
        call GetGeometrySize(Me%ObjGeometry,                                  &
                             Size     = Me%ExternalVar%Size,                  &
                             WorkSize = Me%ExternalVar%WorkSize,              &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleLagrangian - ERR04'

        !Input data file
        Message   ='ASCII file used to construct lagrangian module'
        call ReadFileName('PARTIC_DATA', Me%Files%ConstructData,              &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleLagrangian - ERR05'



        !Transient HDF File
        Message   ='Instant fields of lagrangian particles in HDF format.'
        call ReadFileName('PARTIC_HDF', Me%Files%TransientHDF,                &
                           Message = Message,                                            &
                           TIME_END = Me%ExternalVar%EndTime,                 &
                           Extension = 'hdf', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleLagrangian - ERR06'


        !Final file
        Message   ='Final Particle Conditions.'
        call ReadFileName('PARTIC_FIN', Me%Files%Final,                       &
                           Message = Message,                                            &
                           TIME_END = Me%ExternalVar%EndTime,                 &
                           Extension = 'ptf', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleLagrangian - ERR07'


        !Initial Conditions
        Message   ='Initial Particle Conditions.'
        call ReadFileName('PARTIC_INI', Me%Files%Initial,                     &
                           Message = Message, TIME_END = Me%ExternalVar%Now,  &
                           Extension = 'ptf',STAT = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_  ) then
            write(*,*)  
            write(*,*) 'Initial file not found.'
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleLagrangian - ERR08'

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then
            write(*,*)  
            write(*,*) 'Keyword for the initial file not found in nomfich.dat.'
            write(*,*) 'ConstructGlobalVariables - ModuleLagrangian - WRN01'
            write(*,*)  

        else if (STAT_CALL .EQ. SUCCESS_             ) then
            continue
        else
            stop 'ConstructGlobalVariables - ModuleLagrangian - ERR09'
        end if cd1  

        !----------------------------------------------------------------------

    end subroutine ConstructGlobalVariables

    !--------------------------------------------------------------------------

    subroutine ConstructParticleGrid

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: i, j, STAT_CALL
        real,    dimension (:, :), pointer          :: XX_IE, YY_IE
        integer                                     :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG
        integer                                     :: GRID_COORD, CoordType, NLRD
        real,    dimension (:),    pointer          :: XX, YY
        logical                                     :: DistortionOn
        integer, dimension (:, :), pointer          :: DefineCellsMap



        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB

        !Gets Coordinate Types List
        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(Me%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR10'

        call GetCheckDistortion (Me%ObjHorizontalGrid, DistortionOn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR20'


        if (DistortionOn) then
            Me%Grid%GridType = Grid2D
        else
            if (CoordType == GEOG        ) Me%Grid%GridType = Grid2D
            if (CoordType == UTM         ) Me%Grid%GridType = Grid1D
            if (CoordType == MIL_PORT    ) Me%Grid%GridType = Grid1D 
            if (CoordType == SIMPLE_GEOG ) Me%Grid%GridType = Grid2D
            if (CoordType == GRID_COORD  ) Me%Grid%GridType = Grid1D
            if (CoordType == NLRD        ) Me%Grid%GridType = Grid1D
        endif


        if (Me%Grid%GridType == Grid1D) then

            nullify (Me%Grid%ParticXX)
            nullify (Me%Grid%ParticYY)

            !Gets XX and YY
            nullify(XX, YY)
            call GetHorizontalGrid (Me%ObjHorizontalGrid, XX = XX, YY = YY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR30'
            
            allocate(Me%Grid%ParticX(JLB:JUB))
            allocate(Me%Grid%ParticY(ILB:IUB))

            Me%Grid%ParticX(WS_JLB) = 0.
            Me%Grid%ParticY(WS_ILB) = 0.


            do j = WS_JLB+1, WS_JUB + 1
                Me%Grid%ParticX(j) = Me%Grid%ParticX(j-1) +  (XX(j) - XX(j-1))
            enddo

            do i = WS_ILB+1, WS_IUB + 1
                Me%Grid%ParticY(i) = Me%Grid%ParticY(i-1) + (YY(i) - YY(i-1))
            enddo


            call UnGetHorizontalGrid (Me%ObjHorizontalGrid, XX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR40'

            call UnGetHorizontalGrid (Me%ObjHorizontalGrid, YY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR50'

        else

            nullify (Me%Grid%ParticX)
            nullify (Me%Grid%ParticY)

            !Gets Horizontal Grid
            call GetHorizontalGrid(Me%ObjHorizontalGrid, XX_IE = XX_IE, YY_IE = YY_IE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ER60'

            call GetDefineCellsMap(Me%ObjHorizontalGrid, DefineCellsMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR70'

            allocate(Me%Grid%ParticXX(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR80'

            allocate(Me%Grid%ParticYY(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR90'


            Me%Grid%ParticXX(:, :) = 0.
            Me%Grid%ParticYY(:, :) = 0.

            do j = WS_JLB, WS_JUB
            do i = WS_ILB,   WS_IUB

                if (DefineCellsMap(i, j) == 1) then

                    Me%Grid%ParticXX(i  , j  ) = (XX_IE(i,  j  ) + XX_IE(i+1, j  )) / 2.
                    Me%Grid%ParticXX(i  , j+1) = (XX_IE(i,  j+1) + XX_IE(i+1, j+1)) / 2.
                    Me%Grid%ParticYY(i  , j  ) = (YY_IE(i,  j  ) + YY_IE(i  , j+1)) / 2.
                    Me%Grid%ParticYY(i+1, j  ) = (YY_IE(i+1,j  ) + YY_IE(i+1, j+1)) / 2.
                endif
            
            enddo
            enddo

    
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, XX_IE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR100'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, YY_IE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR110'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DefineCellsMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR120'

        endif


        !Calculates the area of each Grid Cell
        allocate (Me%Grid%AreaCell(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangian - ERR130'

        Me%Grid%AreaCell = null_real
        do j = WS_JLB, WS_JUB
        do i = WS_ILB, WS_IUB

            select case (Me%Grid%GridType)
            
            case (Grid1D)
            
                Me%Grid%AreaCell(i, j) = (Me%Grid%ParticX(j+1) - Me%Grid%ParticX(j)) * &
                                                    (Me%Grid%ParticY(i+1) - Me%Grid%ParticY(i)) 

            case (Grid2D)

                Me%Grid%AreaCell(i, j) = (Me%Grid%ParticXX(i, j+1) - Me%Grid%ParticXX(i, j)) * &
                                                    (Me%Grid%ParticYY(i+1, j) - Me%Grid%ParticYY(i, j)) 

            end select

        enddo
        enddo


    end subroutine ConstructParticleGrid    

    !--------------------------------------------------------------------------

    subroutine ConstructOrigins

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: flag
        integer                                     :: STAT_CALL
        integer                                     :: ClientNumber, ClientNumberClone
        real                                        :: DT, DT_PARTIC, MinError
        logical                                     :: BlockFound, PropertyFound
        type (T_Origin), pointer                    :: NewOrigin, OriginalOrigin
        character(3)                                :: Aux
        character(StringLength)                     :: EmissionSpatial
        character(StringLength)                     :: EmissionTemporal
        character(PathLength)                       :: String, String2, RootPath
        real, dimension(:), allocatable             :: Aux2
        real, dimension(1:2)                        :: Position
        real                                        :: Depth, Xorig, Yorig, GridAngle
        logical                                     :: HaveOrigin, DistortionOn
        type (T_Property), pointer                  :: NewProperty
        integer                                     :: i, j, k
        logical                                     :: ret
        integer                                     :: PropertyID, DensityMethod
        logical                                     :: WP_HaveProperty, SedimentDefined = .false.
        logical                                     :: SalOK = .false., TempOK = .false., PressureCorrection
        logical                                     :: FoundCloneOrigin, ClonesExist
        integer                                     :: Nmax, no, NFirstExtraction
        logical                                     :: FirstExtraction

        !Begin-----------------------------------------------------------------


        call GetOutPutTime(Me%ObjEnterData,                               &
                           CurrentTime = Me%ExternalVar%Now,              &
                           EndTime     = Me%ExternalVar%EndTime,          &
                           keyword     = 'OUTPUT_TIME',                   &
                           SearchType  = FromFile,                        &
                           OutPutsTime = Me%OutPut%OutTime,               &
                           OutPutsOn   = Me%OutPut%Write_,                &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR10'

        if (Me%OutPut%Write_) then
            Me%OutPut%NextOutPut = 1
        endif

        call GetOutPutTime(Me%ObjEnterData,                                         &
                           CurrentTime = Me%ExternalVar%Now,                        &
                           EndTime     = Me%ExternalVar%EndTime,                    &
                           keyword     = 'RESTART_FILE_OUTPUT_TIME',                &
                           SearchType  = FromFile,                                  &
                           OutPutsTime = Me%OutPut%RestartOutTime,                  &
                           OutPutsOn   = Me%OutPut%WriteRestartFile,                &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR11'

        if(Me%OutPut%WriteRestartFile)then

            Me%OutPut%NextRestartOutput = 1

        end if 

        !Checks wether to overwrite the Restart File OR not
        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='RESTART_FILE_OVERWRITE',                            &
                     ClientModule ='ModuleLagrangian',                                  &
                     Default      = .true.,                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR12'

        !Output Concentration Type
        call GetData(Me%OutPut%OutPutConcType,                                          &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_CONC',                                       &
                     ClientModule ='ModuleLagrangian',                                  &
                     Default      = Maximum,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR20'

        !Checks if the users wants to output the maximum tracer concentration in each cell
        call GetData(Me%OutPut%ConcMaxTracer,                                           &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_MAX_TRACER',                                 &
                     ClientModule ='ModuleLagrangian',                                  &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR30'


        call GetData(Me%OutPut%OriginEnvelope,                                          & 
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_ORIGIN_ENVELOPE',                            &
                     ClientModule ='ModuleLagrangian',                                  &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR40'


        !Time series (will be valid to all the properties)
        call GetData(Me%WritesTimeSerie,                                                &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     Keyword      = 'TIME_SERIE',                                       &
                     ClientModule = 'ModuleLagrangian',                                 &
                     Default      = .false.,                                            &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR50'


        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR70'

        !Get time step for particle from file
        DT_PARTIC = null_real
        call GetData(DT_PARTIC,                                                         &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='DT_PARTIC',                                         &
                     ClientModule ='ModuleLagrangian',                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR80'
                          
        !Checks time steps
        if (flag .EQ. 0) then
            DT_PARTIC = DT
        else
            if (DT_PARTIC > DT) then
                if (amod(DT_PARTIC, DT) == 0.0) then
                    !OK
                else
                    write(*,*)'Particle Time step : ', DT_PARTIC
                    write(*,*)'Model    Time step : ', DT
                    write(*,*)'Particle Time step must multiply or submultiply of Model Time step'
                    stop      'ConstructOrigins - ModuleLagrangian - ERR90'
                endif
            elseif (DT > DT_PARTIC) then
                !if (amod(DT, DT_PARTIC) == 0.0) then
                    !OK
                !else
                !    write(*,*)'Particle Time step : ', DT_PARTIC
                !    write(*,*)'Model    Time step : ', DT
                !    write(*,*)'Particle Time step must multiply or submultiply of Model Time step'
                !    stop      'ConstructOrigins - ModuleLagrangian - ERR100'
                !endif
            
                !The run period must be a multiple of the model DT
                !The abs function is used, to avoid rounding erros
                !The old way was removed, to be able to run with Timesteps lower tehn 1 sec
                !Frank Dec - 2000
                MinError = min (abs(mod (DT, DT_PARTIC)),                           &
                                abs(DT_PARTIC - mod (DT, DT_PARTIC)))
                if (MinError >= 1.e-5) then
                    write(*,*)' Time step error DT_PARTIC - Run period must be a multiple of DT'
                    stop      'ConstructOrigins - ModuleLagrangian - ERR110.' 
                endif
            endif
        end if

        !Sets time step
        Me%DT_Partic   = DT_PARTIC
        Me%NextCompute = Me%ExternalVar%Now + Me%DT_Partic
        Me%Now         = Me%ExternalVar%Now
        
        !Boxes data file
        call GetData(Me%Files%BoxDataFile,                                              &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='PARTIC_BOX',                                        &
                     ClientModule ='ModuleLagrangian',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR120'
        if (flag == 1) then

            !Starts BoxDif
            call StartBoxDif(BoxDifID         = Me%ObjBoxDif,                           &
                             TimeID           = Me%ObjTime,                             &
                             HorizontalGridID = Me%ObjHorizontalGrid,                   &
                             BoxesFilePath    = Me%Files%BoxDataFile,                   &
                             WaterPoints3D    = Me%ExternalVar%Waterpoints3D,           &
                             Size3D           = Me%ExternalVar%Size,                    &
                             WorkSize3D       = Me%ExternalVar%WorkSize,                &
                             STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR130'
            
        endif

        !Boxes data file
        call GetData(Me%Files%MonitorBox,                                               &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='MONITOR_BOX',                                       &
                     ClientModule ='ModuleLagrangian',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR140'
        if (flag == 1) then

            !Starts BoxDif
            call StartBoxDif(BoxDifID         = Me%ObjMonBox,                           & 
                             TimeID           = Me%ObjTime,                             &
                             HorizontalGridID = Me%ObjHorizontalGrid,                   &
                             BoxesFilePath    = Me%Files%MonitorBox,                    &
                             WaterPoints3D    = Me%ExternalVar%Waterpoints3D,           &
                             Size3D           = Me%ExternalVar%Size,                    &
                             WorkSize3D       = Me%ExternalVar%WorkSize,                &
                             STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR150'

            Me%State%Monitor = .true.


            !Name of the Property
            call GetData(Me%Monitor%MassProperty,                                       &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='MONITOR_BOX_PROP_MASS',                         &
                         ClientModule ='ModuleLagrangian',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR160'

            if (flag == 1) then

                !Checks if the property is recognized by the model       
                if (.not. CheckPropertyName(trim(Me%Monitor%MassProperty), Me%Monitor%MassPropertyID)) then
                    write (*,*)'Unknown Property : ', trim(Me%Monitor%MassProperty)
                    stop 'ConstructOrigins - ModuleLagrangian - ERR170'
                end if

                Me%State%MonitorPropMass = .true.
            end if
            
            if (Me%State%MonitorPropMass) then
            
                !Output Concentration Type
                call GetData(Me%Monitor%MassFractionType,                                       &
                             Me%ObjEnterData,                                                   &
                             flag,                                                              &
                             SearchType   = FromFile,                                           &
                             keyword      ='MONITOR_BOX_MASS_FRACTION',                         &
                             ClientModule ='ModuleLagrangian',                                  &
                             Default      = Arithmetic,                                         &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR175'

                ! Minimum Tracer Concentration Value for Average Concentration 
                !and Contamination Probability in Each Monitoring Box

                call GetData(Me%Monitor%MonitorBox_TracerMinConc,                           &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType   = FromFile,                                       &
                             keyword      ='MONITOR_BOX_MIN_CONC',                          &
                             ClientModule ='ModuleLagrangian',                              &
                             Default      = 1.,                                             &
                             STAT         = STAT_CALL)
               if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR176'

               ! Depth considered to contamination probability
                call GetData(Me%Monitor%ContaminationDepth,                                 &
                             Me%ObjEnterData,                                               &
                             flag,                                                          &
                             SearchType   = FromFile,                                       &
                             keyword      ='MONITOR_BOX_CONT_DEPTH',               &
                             ClientModule ='ModuleLagrangian',                              &
                             Default      = null_real,                                      &
                             STAT         = STAT_CALL)
               if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR177'
            
            end if

        endif

        !Eulerian boxes data file
        call GetData(Me%Files%EulerianMonitorBox,                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='EULERIAN_MONITOR_BOX',                              &
                     ClientModule ='ModuleLagrangian',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR180'
        
        if (flag == 1) Me%State%EulerianMonitor = .true.            

        !Associates Beaching Probabilities 
        call GetData(Me%State%AssociateBeachProb,                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='ASSOCIATE_BEACH_PROB',                              &
                     Default      = .false.,                                            &
                     ClientModule ='ModuleLagrangian',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR190'

        if (Me%State%AssociateBeachProb) then

            !Maximum distance between particles and coast for particle beaching 
            call GetData(Me%BeachingLimit,                                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='BEACHING_LIMIT',                                &
                         Default      = 5.0,                                            &
                         ClientModule ='ModuleLagrangian',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR200'

            !Outbox Beaching Probability
            call GetData(Me%DefaultBeachingProbability,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='DEFAULT_BEACHING_PROB',                         &
                         Default      = 0.5,                                            &
                         ClientModule ='ModuleLagrangian',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR210'


            !Beaching Probabilities Boxes data file
            call GetData(Me%Files%BeachingBoxFileName,                                  &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='BEACHING_BOX_FILENAME',                         &
                         ClientModule ='ModuleLagrangian',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR220'
            if (flag == 1) then

                !Starts BoxDif
                call StartBoxDif(BoxDifID         = Me%ObjBeachingProbBox,              &
                                 TimeID           = Me%ObjTime,                         &
                                 HorizontalGridID = Me%ObjHorizontalGrid,               &
                                 BoxesFilePath    = Me%Files%BeachingBoxFileName,       &
                                 WaterPoints3D    = Me%ExternalVar%Waterpoints3D,       &
                                 Size3D           = Me%ExternalVar%Size,                &
                                 WorkSize3D       = Me%ExternalVar%WorkSize,            &
                                 STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR230'

                Me%State%HaveBeachingProbBox = .true.
        
            endif

        endif


        !OVERLAY_VELOCITY
        call GetData(Me%OverLay%Overlay,                                                &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OVERLAY_VELOCITY',                                  &
                     ClientModule ='ModuleLagrangian',                                  &  
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR240'


        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR250'

        FoundCloneOrigin = .false.
        ClonesExist      = .true.

        Nmax        = 1
        FirstExtraction = .true. 
        NFirstExtraction = 0


        !Run module lagrangian online
        call GetData(Me%RunOnline,                                                      &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='RUN_ONLINE',                                        &
                     ClientModule ='ModuleLagrangian',                                  &  
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR255'

!#ifdef _CGI_

        if (Me%RunOnline) then

            call ReadQueryString(Nmax)

        endif

!#endif

        !Number of times it read the lagrangian data looking for origins
        !Only when the _CGI_ option is on is able to read several times 
        !the origin blocks
SB:     do no = 1, Nmax

DW:     do  

            if (FoundCloneOrigin) then

                BlockFound      = .true. 

            else 

                call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                            ClientNumber,                                    &
                                            block_begin, block_end,                          &
                                            BlockFound,                                      &
                                            STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR260'

            endif

BF:         if (BlockFound) then 

                if (FirstExtraction) then 
                    NFirstExtraction = NFirstExtraction + 1
                endif
                    

                !Allocates a new origin
                call AllocateNewOrigin(NewOrigin)

                if (.not. FoundCloneOrigin) then
                    nullify(OriginalOrigin)
                    OriginalOrigin => NewOrigin
                endif

                !Gets its name    
                call GetData(NewOrigin%Name,                                            &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlock,                                  &
                             keyword      ='ORIGIN_NAME',                               &
                             ClientModule ='ModuleLagrangian',                          &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR270'

                if (flag == 0) then
                    Aux = ' '
                    write(Aux,'(i3)') Me%nOrigins + 1
                    NewOrigin%Name = trim(adjustl('Origin_'//trim(adjustl(Aux))))
                end if

                !Old Origin    
                call GetData(NewOrigin%Old,                                              &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='OLD',                                        &
                             ClientModule ='ModuleLagrangian',                           &
                             Default      = .false.,                                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR280'

                if (NewOrigin%Old) then
                    Me%State%ContCalc = .true.
                endif

                !Group ID
                call GetData(NewOrigin%GroupID,                                          &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='GROUP_ID',                                   &
                             ClientModule ='ModuleLagrangian',                           &
                             Default      = 1,                                           &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR290'


                !Gets the spatial emission type
                call GetData(EmissionSpatial,                                            &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='EMISSION_SPATIAL',                           &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR300'

                !Gets the temporal emission type
                call GetData(EmissionTemporal,                                           &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='EMISSION_TEMPORAL',                          &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR310'

                !Converts Strings to int
                call ConstructEmissionType(NewOrigin, EmissionSpatial, EmissionTemporal)

                !
ET:             if (NewOrigin%EmissionTemporal == Continuous_) then

                    !Gets the interval between emissions
                    call GetData(NewOrigin%DT_Emit,                                      &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='DT_EMIT',                                &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = Me%DT_Partic,                 &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR320'

                    !Gets the interval between emissions
                    call GetData(NewOrigin%StartEmission,                                &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='START_PARTIC_EMIT',                      &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = Me%ExternalVar%BeginTime,     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR330'

                    !Gets the interval between emissions
                    call GetData(NewOrigin%StopEmission,                                 &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='STOP_PARTIC_EMIT',                       &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = Me%ExternalVar%EndTime,       &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR340'


                    !Gets flow associated to a continuous emission
                    if (NewOrigin%EmissionSpatial == Point_) then


                        !Flow variable in time
                        call GetData(NewOrigin%FlowVariable,                             &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='FLOW_VARIABLE',                      &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR360'

                        if (NewOrigin%FlowVariable) then

                            !Discharge file name
                            call GetData(NewOrigin%DischargeFile,                            &
                                         Me%ObjEnterData,                                    &
                                         flag,                                               &
                                         SearchType   = FromBlock,                           &
                                         keyword      ='DISCHARGE_FILE',                     &
                                         ClientModule ='ModuleLagrangian',                   &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR370'


                            call StartTimeSerieInput(NewOrigin%TimeSerieInput,           &
                                                     NewOrigin%DischargeFile,            &
                                                     Me%ObjTime,                         &
                                                     STAT = STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR380'

                            !Discharge file name
                            call GetData(NewOrigin%FlowColumn,                               &
                                         Me%ObjEnterData,                                    &
                                         flag,                                               &
                                         SearchType   = FromBlock,                           &
                                         keyword      ='FLOW_COLUMN',                        &
                                         ClientModule ='ModuleLagrangian',                   &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR390'

                        else

                            !If not variable, get flow of the origin
                            call GetData(NewOrigin%Flow,                                     &
                                         Me%ObjEnterData,                                    &
                                         flag,                                               &
                                         SearchType   = FromBlock,                           &
                                         keyword      ='FLOW',                               &
                                         ClientModule ='ModuleLagrangian',                   &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR350'
                            if (flag /= 1) then
                                write(*,*)'Keyword FLOW not defined at origin :',trim(adjustl(NewOrigin%Name))
                                stop      'ConstructOrigins - ModuleLagrangian - ERR27'
                        
                        endif


                        endif



                        call GetData(NewOrigin%MovingOrigin,                             &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     default      = .false.,                             &
                                     keyword      ='MOVING_ORIGIN',                      &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR400'

                        if (NewOrigin%MovingOrigin) then

                            call GetData(NewOrigin%MovingOriginFile,                     &
                                         Me%ObjEnterData,                     &
                                         flag,                                           &
                                         SearchType   = FromBlock,                       &
                                         keyword      ='MOVING_ORIGIN_FILE',             &
                                         ClientModule ='ModuleLagrangian',               &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_ .or. flag /= 1) stop 'ConstructOrigins - ModuleLagrangian - ERR410'

                            call GetData(NewOrigin%MovingOriginUnits,                    &
                                         Me%ObjEnterData,                     &
                                         flag,                                           &
                                         default      = 'Cells',                         &
                                         SearchType   = FromBlock,                       &
                                         keyword      ='MOVING_ORIGIN_UNITS',            &
                                         ClientModule ='ModuleLagrangian',               &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR420'


                            call GetData(NewOrigin%MovingOriginColumnX,                  &
                                         Me%ObjEnterData,                     &
                                         flag,                                           &
                                         SearchType   = FromBlock,                       &
                                         keyword      ='MOVING_ORIGIN_COLUMN_X',         &
                                         ClientModule ='ModuleLagrangian',               &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_ .or. flag /= 1) stop 'ConstructOrigins - ModuleLagrangian - ERR430'

                            call GetData(NewOrigin%MovingOriginColumnY,                  &
                                         Me%ObjEnterData,                     &
                                         flag,                                           &
                                         SearchType   = FromBlock,                       &
                                         keyword      ='MOVING_ORIGIN_COLUMN_Y',         &
                                         ClientModule ='ModuleLagrangian',               &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_ .or. flag /= 1) stop 'ConstructOrigins - ModuleLagrangian - ERR440'


                            call StartTimeSerieInput(NewOrigin%ObjTimeSerie,             &
                                                     NewOrigin%MovingOriginFile,         &
                                                     Me%ObjTime,                         &
                                                     STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR450'

   
                        endif

                    endif

                else ET

                    if (NewOrigin%EmissionSpatial == Point_ .or.                         &
                        NewOrigin%EmissionSpatial == Accident_ ) then
                        call GetData(NewOrigin%PointVolume,                              &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='POINT_VOLUME',                       &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR460'
                        if (flag /= 1) then
                            write(*,*)'Keyword POINT_VOLUME not defined at origin :',trim(adjustl(NewOrigin%Name))
                            stop      'ConstructOrigins - ModuleLagrangian - ERR29'
                        endif
                    endif
               
                endif ET

                !Time to double volume
                call GetData(NewOrigin%Movement%TVOL200,                                 &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='TVOL200',                                    &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR470'

                if (flag == 1) then
                    NewOrigin%State%VariableGeom = ON
                endif

                
                !Split particle
                call GetData(NewOrigin%Movement%SPLIT_PART,                              &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      = 'SPLIT_PART',                                &
                             Default      = OFF,                                         &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR480'


                !Volume Factor when Particle dies
                call GetData(NewOrigin%Movement%VOLFAC,                                  &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      = 'VOLFAC',                                    &
                             Default      = 10.,                                         &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR490'

                !How to calculate volume increase
                call GetData(String2,                                                    &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      = 'VOLUME_INCREASE',                           &
                             Default      = Char_Double,                                 &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR500'

                select case (trim(adjustl(String2)))

                case (Char_Velocity)

                    NewOrigin%Movement%TVolType = Velocity_

                case (Char_Double)

                    NewOrigin%Movement%TVolType = Double_

                case default

                    write(*,*)'Invalid option for keyword VOLUME_INCREASE'
                    write(*,*)
                    stop 'ConstructOrigins - ModuleLagrangian - ERR510'

                end select


                !Floating Particle
                call GetData(NewOrigin%Movement%Float,                                   &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='FLOAT',                                      &
                             Default      = OFF,                                         &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR520'

            
                !Calculates Plume
                call GetData(NewOrigin%State%FarFieldBuoyancy,                           &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='COMPUTE_BUOYANCY',                           &
                             Default      = OFF,                                         &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR530'

                
                if (NewOrigin%State%FarFieldBuoyancy) then

                    Me%State%FarFieldBuoyancy = .true.

                    NewOrigin%Movement%InitialVelocityW = 0.

                endif

           
                !Calculates Plume
                call GetData(NewOrigin%State%ComputePlume,                               &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='COMPUTE_PLUME',                              &
                             Default      = OFF,                                         &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR540'

                !Plume Data
                if (NewOrigin%State%ComputePlume) then


                    Me%State%ComputePlume = .true.
                    
                    if (.not. NewOrigin%State%FarFieldBuoyancy) then
                        write(*,*) 'The MOHID JET is ON but'
                        write(*,*) 'no buoyancy in the far field is being compute'
                        write(*,*) 'ConstructOrigins - ModuleLagrangian - WRN10'
                    endif

                    !Momentum balance in horizontal direcction
                    !Increase of volume tracer due to shear effect
                    call GetData(NewOrigin%State%PlumeShear,                             &
                                 Me%ObjEnterData,                                        &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='PLUME_SHEAR',                            &
                                 Default      = .true.,                                  &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR550'


                                        
                    !Plume Friction
                    call GetData(NewOrigin%Movement%CoefInitialMixing,                   &
                                 Me%ObjEnterData,                                        &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='COEF_INITIAL_MIXING',                    &
                                 Default      = 1.0,                                     &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR560'


                    !Plume Jet Data file
                    call GetData(NewOrigin%Movement%JetDataFile,                         &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='JET_DATA_FILE',                          &
                                 Default      = "********.***",                          &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR570'
                    if (flag == 0)             stop 'ConstructOrigins - ModuleLagrangian - ERR580'

                    !Time interval for the actualization of Plume Jet properties
                    call GetData(NewOrigin%Movement%JetDT,                               &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='JET_DT',                                 &
                                 Default      = 600.0,                                   &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR590'


                    NewOrigin%Movement%NextJetActualization = Me%ExternalVar%BeginTime

                endif
                
                if (Me%State%FarFieldBuoyancy .or. Me%State%ComputePlume) then
                    !Density Evolution
                    call GetData(NewOrigin%Movement%DensityMethod,                       &
                                 Me%ObjEnterData,                                        &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='DENSITY_METHOD',                         &
                                 Default      = UNESCOState_,                            &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR600'
                    
                    call GetData(NewOrigin%Movement%CorrecPress,                        &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlock,                              &
                                 keyword      ='PRESSURE_CORRECTION',                   &
                                 Default      = .true.,                                 &
                                 ClientModule ='ModuleLagrangian',                      &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR610'

                    call GetDensityOptions(Me%ObjWaterProperties, DensityMethod,        &
                                           PressureCorrection, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR620'
                    
                    if (NewOrigin%Movement%DensityMethod/=DensityMethod) then
                        stop 'ConstructOrigins - ModuleLagrangian - ERR630'
                    endif
                    
                    if (.not. (NewOrigin%Movement%CorrecPress .EQV. PressureCorrection)) then
                        stop 'ConstructOrigins - ModuleLagrangian - ERR640'
                    endif                    
               
                endif


                !Type of Accident
AC:             if (NewOrigin%EmissionSpatial == Accident_) then

                    call GetData(NewOrigin%AccidentMethod,                               &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='ACCIDENT_METHOD',                        &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = Fay_,                                    &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR670'

                    call GetData(NewOrigin%AccidentTime,                                 &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='ACCIDENT_TIME',                          &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = Me%ExternalVar%BeginTime,     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR680'


                endif AC

                !Particle Thickness                
                if (NewOrigin%Movement%Float) then
                    Me%State%Wind = ON

                    if (NewOrigin%EmissionSpatial  == Box_ .or.                          &
                        (NewOrigin%EmissionSpatial == Accident_ .and.                    &
                         NewOrigin%AccidentMethod  == Thickness_)) then
                        call GetData(NewOrigin%Movement%ThicknessMeters,                 &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='THICKNESS_METERS',                   &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR690'
                        if (flag /= 1) then
                            write(*,*)'Keyword THICKNESS_METERS not defined at origin :',trim(adjustl(NewOrigin%Name))
                            stop      'ConstructOrigins - ModuleLagrangian - ERR38'
                        endif
                    endif

                endif


                !Movement type
                call GetData(String,                                                     &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='MOVEMENT',                                   &
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR700'

MO:             if (flag == 1) then

                    select case(trim(adjustl(String)))
                
                    case(Char_SullivanAllen)
                        NewOrigin%Movement%MovType = SullivanAllen_

                        ! ---> Definition of the Horizontal and Vertical variance in 
                        !      the form of a percentage of the average velocity      
                        ! UStandardDeviation = VarVelHX * Vel + VarVelH
                        call GetData(NewOrigin%Movement%VarVelHX,                        &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='VARVELHX',                           &
                                     ClientModule ='ModuleLagrangian',                   &  
                                     Default      = 0.2,                                 &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR710'


                        call GetData(NewOrigin%Movement%VarVelH,                         &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='VARVELH',                            &
                                     ClientModule ='ModuleLagrangian',                   &  
                                     Default      = 0.0,                                 &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR720'


                        call GetData(String2,                                            &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='TURB_V',                             &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR730'

TURB_V:                 if (flag == 1) then   

                            select case(trim(adjustl(String2)))
                            
                            case(Char_Constant)

                                NewOrigin%Movement%StandardDeviationType = VerticalTurbConstant

                                call GetData(NewOrigin%Movement%VarVelVX,                &
                                             Me%ObjEnterData,                 &
                                             flag,                                       &
                                             SearchType   = FromBlock,                   &
                                             keyword      ='VARVELVX',                   &
                                             ClientModule ='ModuleLagrangian',           &  
                                             Default      = 0.0,                         &
                                             STAT         = STAT_CALL)             
                                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR740'


                                call GetData(NewOrigin%Movement%VarVelV,                 &
                                             Me%ObjEnterData,                 &
                                             flag,                                       &
                                             SearchType   = FromBlock,                   &
                                             keyword      ='VARVELV',                    &
                                             ClientModule ='ModuleLagrangian',           &  
                                             Default      = 0.0,                         &
                                             STAT         = STAT_CALL)             
                                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR750'

                            case(Char_Profile )

                                NewOrigin%Movement%StandardDeviationType = VerticalTurb
                                
                                !Will need shear velocity    
                                Me%State%ShearVel = ON
                                
                                !Sets unused variables to dummy
                                NewOrigin%Movement%VarVelVX  = null_real
                                NewOrigin%Movement%VarVelV   = null_real

                            case default

                                write(*,*)'Invalid option for TURB_V'
                                stop 'ConstructOrigins - ModuleLagrangian - ERR760'

                            end select

                        else TURB_V
                        
                            write(*,*)'Keyword TURB_V not found'
                            stop 'ConstructOrigins - ModuleLagrangian - ERR780'

                        endif TURB_V

                    case(Char_NotRandom    )

                        NewOrigin%Movement%MovType = NotRandom_

                    case default

                        write(*,*)'Invalid horizontal movement keyword'
                        stop 'ConstructOrigins - ModuleLagrangian - ERR790'
                 
                    end select

                else MO

                    write(*,*)'Keyword MOVEMENT not defined at origin', trim(adjustl(NewOrigin%Name))
                    stop      'ConstructOrigins - ModuleLagrangian - ERR800'
                end if MO


                !Advection
                call GetData(NewOrigin%Movement%Advection,                               &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='ADVECTION',                                  &
                             ClientModule ='ModuleLagrangian',                           &  
                             Default      = ON,                                          &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR810'

                !Kill Land Particles
                call GetData(NewOrigin%Movement%KillLandParticles,                       &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='KILL_LAND_PARTICLES',                        &
                             ClientModule ='ModuleLagrangian',                           &  
                             Default      = OFF,                                         &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR820'



                !WINDCOEF
                call GetData(NewOrigin%Movement%WindTransferCoef,                        &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='WINDCOEF',                                   &
                             ClientModule ='ModuleLagrangian',                           &  
                             Default      = 0.03,                                        &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR830'


                allocate(Aux2(2))

                Aux2(1:2) = 0

                !WINDXY
                call GetData(Aux2,                                                       &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='WINDXY',                                     &
                             ClientModule ='ModuleLagrangian',                           &  
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR840'

                
                if (flag == 2) then
                    NewOrigin%Movement%WindOriginON = .true.
                else if (flag == 0) then
                    NewOrigin%Movement%WindOriginON = .false.
                endif

                NewOrigin%Movement%WindX = Aux2(1)
                NewOrigin%Movement%WindY = Aux2(2)

                deallocate(Aux2)

                !SEDIMENTATION
                call GetData(String,                                                     &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='SEDIMENTATION',                              &
                             ClientModule ='ModuleLagrangian',                           &  
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR850'

SE:             if (flag == 1) then
 
                    select case(trim(adjustl(String)))
                
                    case(Char_Stokes )
                    
                        NewOrigin%Movement%SedimentationType = Stokes_

                        call GetData(NewOrigin%Movement%D50,                             &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='D50',                                &
                                     default      = 0.002,                               & 
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR860'

                        NewOrigin%State%Sedimentation = ON

                    case(Char_Imposed)
                    
                        NewOrigin%Movement%SedimentationType = Imposed_

                        call GetData(NewOrigin%Movement%SedVel,                          &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='SED_VELOCITY',                       &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR870'
                        if (flag == 0) then
                            write(*,*)'Sedimentation velocity not defined, keyword SED_VELOCITY'
                            stop 'ConstructOrigins - ModuleLagrangian - ERR880'
                        endif

                        NewOrigin%State%Sedimentation     = ON

                    case default
                    
                        write(*,*)'Invalid Sedimentaion type, keyword SEDIMENTATION'
                        stop 'ConstructOrigins - ModuleLagrangian - ERR890'

                    end select

                end if SE

               !Particles with deposition 
                call GetData(NewOrigin%State%Deposition,                                 &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='DEPOSITION',                                 &
                             default      = .false.,                                     & 
                             ClientModule ='ModuleLagrangian',                           &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                               &
                    call SetError(FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR900')

                if (NewOrigin%State%Deposition .and. .not.NewOrigin%State%Sedimentation) &
                    call SetError(FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR910')


DE:             if (NewOrigin%State%Deposition) then

                    Me%State%Deposition = ON

                    call GetData(NewOrigin%Deposition%TauErosion,                        &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='TAU_ERO',                                &
                                 default      = 0.2,                                     &                                  
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                           &
                        call SetError(FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR920')

                    call GetData(NewOrigin%Deposition%TauDeposition,                     &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='TAU_DEP',                                &
                                 default      = 0.1,                                     &                                  
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                           &
                        call SetError(FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR930')

                   call GetData(NewOrigin%Deposition%BottomDistance,                    &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='BOTTOM_DISTANCE',                        &
                                 default      = 0.1,                                     &                                                                   
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                           &
                        call SetError(FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR940')

                   call GetData(NewOrigin%Deposition%Tdecay,                             &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='TIME_DECAY',                             &
                                 !defaul 2 days
                                 default      = 172800.,                                 &                                                                   
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                           &
                        call SetError(FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR950')

                   call GetData(NewOrigin%Deposition%BottomEmission,                     &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='BOTTOM_EMISSION',                        &
                                 !by default the particle are emitted in the water column
                                 default      = .false.,                                 &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                           &
                        call SetError(FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR960')

                   call GetData(NewOrigin%Deposition%ErosionRate,                        &
                                 Me%ObjEnterData,                                        &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='EROSION_RATE',                           &
                                 !by default the erosion rate in 5e-2 g/m2/s             
                                 !This value make sense if the concentration is in mg/l = g/m3
                                 default      = 5.e-2,                                   &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                           &
                        call SetError(FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR970')

                endif DE


                !Min Sedimentation velocity
                call GetData(NewOrigin%Movement%MinSedVel,                               &
                             Me%ObjEnterData,                                            &
                             flag,                                                       &
                             SearchType   =  FromBlock,                                  &
                             keyword      = 'MIN_SED_VELOCITY',                          &
                             default      =  0.0,                                        &
                             ClientModule = 'ModuleLagrangian',                          &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR980'


                !Reads parameter specific to cada Spatial emission type
PA:             if (NewOrigin%EmissionSpatial == Point_ .or.                             &
                    NewOrigin%EmissionSpatial == Accident_ ) then

                    !Gets the number of particles to emit
                    call GetData(NewOrigin%NbrParticlesIteration,                        &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='NBR_PARTIC',                             &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = 1,                                       &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR990'


                    !Horizontal position in meters
                    Position(:) = null_real
                    HaveOrigin  = .false.
                    call GetData(Position,                                               &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='POSITION_METERS',                        &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1000'

                    if (flag == 2) then

                        NewOrigin%Position%X = Position(1)
                        NewOrigin%Position%Y = Position(2)

                        call Convert_XY_CellIJ(NewOrigin%Position)
                        call Convert_CellIJ_IJ(NewOrigin%Position)

                        HaveOrigin                      = .true.

                    else
                    
                        !Horizontal position in grid cells
                        call GetData(Position,                                          &
                                     Me%ObjEnterData,                                   &
                                     flag,                                              &
                                     SearchType   = FromBlock,                          &
                                     keyword      ='POSITION_CELLS',                    &
                                     ClientModule ='ModuleLagrangian',                  &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1010'

                        if (flag == 2) then

                            NewOrigin%Position%CellI = Position(2)
                            NewOrigin%Position%CellJ = Position(1)

                            call Convert_CellIJ_XY(NewOrigin%Position)
                            call Convert_CellIJ_IJ(NewOrigin%Position)

                            HaveOrigin                      = .true.

                        else 
                        
                            call GetCheckDistortion (Me%ObjHorizontalGrid, DistortionOn)
                        
                            !Horizontal position in coordinates X, Y
                            call GetData(Position,                                  &
                                         Me%ObjEnterData,                           &
                                         flag,                                      &
                                         SearchType   = FromBlock,                  &
                                         keyword      ='POSITION_COORDINATES',      &
                                         ClientModule ='ModuleLagrangian',          &
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1020'

                            if (DistortionOn .and. flag > 0) stop 'ConstructOrigins - ModuleLagrangian - ERR1025'

                            !Gets Origin
                            call GetGridOrigin  (Me%ObjHorizontalGrid, Xorig, Yorig)

                            !Gets GridAngle
                            call GetGridAngle   (Me%ObjHorizontalGrid, GridAngle)

                            call RODAXY(-Xorig,  -Yorig, 0.0, Position(1), Position(2)) 

                            call RODAXY(0.0, 0.0, -GridAngle, Position(1), Position(2))  

                            NewOrigin%Position%X = Position(1)
                            NewOrigin%Position%Y = Position(2)

                            call Convert_Coord1D_CellIJ (NewOrigin%Position,            &
                                                         Me%ExternalVar%XX, Me%ExternalVar%YY)

                            call Convert_CellIJ_XY(NewOrigin%Position)
                            call Convert_CellIJ_IJ(NewOrigin%Position)

                            HaveOrigin                      = .true.

                        endif

                    endif

                    !Checks if a valid horizontal origin was found
                    if (.not. HaveOrigin) then

                        write(*,*)'No Valid Horizontal Location defined for ',trim(adjustl(NewOrigin%Name))
                        stop 'ConstructOrigins - ModuleLagrangian - ERR66'

                    endif

                    !Verifies the horizontal location of the origin
                    i = NewOrigin%Position%I
                    j = NewOrigin%Position%J
                    k = Me%ExternalVar%WorkSize%KUB
                    if (Me%ExternalVar%WaterPoints3D(i, j, k) /= WaterPoint) then
                        write(*,*)'Invalid Location defined for ',trim(adjustl(NewOrigin%Name))
                        write(*,*)'Point [i]:', NewOrigin%Position%I
                        write(*,*)'Point [j]:', NewOrigin%Position%J
                        write(*,*)'Is not a WaterPoint'
                        stop 'ConstructOrigins - ModuleLagrangian - ERR1030'
                    endif


                    !Vertical position
                    Depth      = null_real
                    HaveOrigin = .false.
                    call GetData(Depth,                                                  &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='DEPTH_METERS',                           &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)             
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1040'

                    if (flag == 1) then

                        NewOrigin%Position%Z = Depth

                        call Convert_Z_CellK (NewOrigin, NewOrigin%Position)
                        call Convert_CellK_K (NewOrigin%Position)

                        NewOrigin%Position%DepthDefinition = Meters

                        HaveOrigin                      = .true.

                    else

                        call GetData(Depth,                                              &
                                     Me%ObjEnterData,                         &
                                     flag,                                               &
                                     SearchType   = FromBlock,                           &
                                     keyword      ='DEPTH_CELLS',                        &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1050'

                        if (flag == 1) then
    
                            NewOrigin%Position%CellK = Depth

                            call Convert_CellK_Z (NewOrigin%Position)
                            call Convert_CellK_K (NewOrigin%Position)
                        
                            HaveOrigin                      = .true.

                            NewOrigin%Position%DepthDefinition = Cells

                        endif

                    endif

                    call GetData(NewOrigin%Position%MaintainRelative,                    &
                                 Me%ObjEnterData,                                        &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='MAINTAIN_RELATIVE_POSITION',             &
                                 Default      = .false.,                                 &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)             
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1060'

                    if (NewOrigin%Position%MaintainRelative) then
                    
                        NewOrigin%Position%Depth = Depth

                    endif

                    !Initial Position of floating Particle close to the surface
                    if (NewOrigin%Movement%Float) then

                        i = NewOrigin%Position%I
                        j = NewOrigin%Position%J
                        k = Me%ExternalVar%WorkSize%KUB
                        NewOrigin%Position%Z = Me%ExternalVar%SZZ(i, j, k)
                    
                        call Convert_Z_CellK  (NewOrigin, NewOrigin%Position)
                        call Convert_CellK_K  (NewOrigin%Position)

                        HaveOrigin = .true.

                    endif

                    !Checks if a valid vertical origin was found
                    if (.not. HaveOrigin) then

                        write(*,*)'No Valid Vertical Location defined for ',trim(adjustl(NewOrigin%Name))
                        stop 'ConstructOrigins - ModuleLagrangian - ERR1080'

                    endif

                    !Verifies the location of the origin
                    i = NewOrigin%Position%I
                    j = NewOrigin%Position%J
                    k = NewOrigin%Position%K
                    if (Me%ExternalVar%WaterPoints3D(i, j, k) /= WaterPoint) then
                        write(*,*)'Invalid Location defined for ',trim(adjustl(NewOrigin%Name))
                        write(*,*)'Point [i]:', NewOrigin%Position%I
                        write(*,*)'Point [j]:', NewOrigin%Position%J
                        write(*,*)'Point [k]:', NewOrigin%Position%K
                        write(*,*)'Is not a WaterPoint'
                        stop 'ConstructOrigins - ModuleLagrangian - ERR1090'
                    endif

                endif PA

BX:             if (NewOrigin%EmissionSpatial == Box_) then

                    call GetData(NewOrigin%INCRP,                                        &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='INCRP',                                  &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = 1,                                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1100'


                    !Volume of the particles inside the
                    !If not given OneCell _> One particle
                    call GetData(NewOrigin%ParticleBoxVolume,                            &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='BOXVOLINIC',                             &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = null_real,                               &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1110'

                    !Box Number
                    call GetData(NewOrigin%BoxNumber,                                    &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='BOX_NUMBER',                             &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1120'
                    if (flag /= 1) then
                        write(*,*)'Keyword BOX_NUMBER not defined at origin :', trim(adjustl(NewOrigin%Name))
                        stop 'ConstructOrigins - ModuleLagrangian - ERR1130'
                    endif

                endif BX

                
                !Searches for Properties
DOPROP:         do 
                    call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,            &
                                               property_begin, property_end,             &
                                               PropertyFound, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1140'
                    if (PropertyFound) then

                        !Allocates NewProperty
                        call AllocateNewProperty(NewProperty)

                        !Name of the Property
                        call GetData(NewProperty%Name,                                   &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='NAME',                               &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1150'
                        if (flag == 0) stop 'ConstructOrigins - ModuleLagrangian - ERR1160'


                        !Checks if the property is recognized by the model
                        if (.not. CheckPropertyName(trim(NewProperty%Name), NewProperty%ID)) then
                            write (*,*)'Unknown Property : ', trim(NewProperty%Name)
                            stop 'ConstructOrigins - ModuleLagrangian - ERR1170'
                        end if

                        !Units of the Property
                        call GetData(NewProperty%Units,                                  &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='UNITS',                              &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1180'

                        !Concentration of the Property
                        call GetData(NewProperty%Concentration,                          &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='CONCENTRATION',                      &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1190'
                        if (flag == 0) stop 'ConstructOrigins - ModuleLagrangian - ERR1200'

                        !Concentration variable in time
                        call GetData(NewProperty%ConcVariable,                           &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='CONC_VARIABLE',                      &
                                     default      = .false.,                             &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1210'

                        if (NewProperty%ConcVariable) then
   
                            !Concentration column
                            call GetData(NewProperty%ConcColumn,                             &
                                         Me%ObjEnterData,                                    &
                                         flag,                                               &
                                         SearchType   = FromBlockInBlock,                    &
                                         keyword      ='CONC_COLUMN',                        &
                                         ClientModule ='ModuleLagrangian',                   &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1220'
                            if (flag == 0) stop 'ConstructOrigins - ModuleLagrangian - ERR1230'

                        endif


                        !GEOMETRIC
                         call GetData(NewProperty%Geometric,                                 &
                                         Me%ObjEnterData,                                    &
                                         flag,                                               &
                                         SearchType   = FromBlockInBlock,                    &
                                         keyword      ='GEOMETRIC',                          &
                                         ClientModule ='ModuleLagrangian',                   &
                                         Default      = 0,                                   &    
                                        STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1230'

                        !No WQM
                         call GetData(NewProperty%NoWQM,                                     &
                                         Me%ObjEnterData,                                    &
                                         flag,                                               &
                                         SearchType   = FromBlockInBlock,                    &
                                         keyword      ='NOWQM',                              &
                                         ClientModule ='ModuleLagrangian',                   &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1240'


                        !MinimunConcentration 
                        call GetData(NewProperty%Min_concentration,                      &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='MIN_CONCENTRATION',                  &
                                     ClientModule ='ModuleLagrangian',                   &
                                     Default      = 0.0,                                 &    
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1250'

                        !AmbientConcentration 
                        call GetData(NewProperty%AmbientConcentration,                   &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='AMBIENT_CONC',                       &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1260'

                        if (flag == 1) then
                            NewProperty%HaveAmbientConcentration = ON
                        else
                            !Checks if the WP module have this property
                            WP_HaveProperty = WaterPropertyExists (                      &
                                                    Me%ObjWaterProperties,               &
                                                    NewProperty%ID, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1270'

                            if (.not. WP_HaveProperty) then
                                write(*,*)'No Ambient Concentration for the Property : ',trim(NewProperty%Name)
                                write(*,*)'Use the KeyWord                           : AMBIENT_CONC'
                                write(*,*)'OR define the property in the Eulerian Module'
                                stop      'ConstructOrigins - ModuleLagrangian - ERR1280'
                            endif
                        endif

    
                        !Extra Stuff
                        ret = CheckPropertyName (NewProperty%Name, PropertyID)
                        select case (PropertyID)

                        !T90
                        case (Fecal_Coliforms_)

                            call GetData(NewProperty%T90Variable,                        &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='T90_VARIABLE',                   &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = .false.,                         &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1290'


                            if (NewProperty%T90Variable) then
                                
                                if (NewProperty%T90Variable) Me%State%T90Variable = .true.
                                call GetData(NewProperty%T90Var_Method,                  &
                                             Me%ObjEnterData,                            &
                                             flag,                                       &
                                             SearchType   = FromBlockInBlock,            &
                                             keyword      ='T90_VAR_METHOD',             &
                                             ClientModule ='ModuleLagrangian',           &
                                            STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1300'

                                if (flag /= 1) then
                                    write(*,*)'Keyword T90_VAR_METHOD not defined at origin :',trim(adjustl(NewOrigin%Name))
                                    stop 'ConstructOrigins - ModuleLagrangian - ERR94b'
                                endif

                            endif


                            call GetData(NewProperty%T90,                                &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='T90',                            &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = 7200.,                           &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1310'


                        end select

                        call GetData(NewProperty%WaterPartition%ON,                      &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='PARTITION_WATER',                    &
                                     ClientModule ='ModuleLagrangian',                   &
                                     Default      =.false.,                              &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                       &
                            call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1320')

WP:                     if (NewProperty%WaterPartition%ON) then

                            Me%State%Partition = ON
                            NewOrigin%State%Partition     = ON

                            call GetData(NewProperty%WaterPartition%Coefficient,         &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='PARTITION_COEF_WATER',           &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = 0.90,                            &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                   &
                                call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1330')

                            call GetData(NewProperty%WaterPartition%TransferRate,        &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='PARTITION_RATE_WATER',           &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = 1e-3,                            &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                   &
                                call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1340')


                            call GetData(NewProperty%WaterPartition%CoupleProp,          &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='PARTITION_COUPLE_WATER',         &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = 0.,                              &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                   &
                                call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1350')

                        endif WP

                        call GetData(NewProperty%SedimentPartition%ON,                   &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='PARTITION_SED',                      & 
                                     ClientModule ='ModuleLagrangian',                   &
                                     Default      =.false.,                              &
                                     STAT         = STAT_CALL)                          
                        if (STAT_CALL /= SUCCESS_)                                       &
                            call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1360')

                        if (NewProperty%SedimentPartition%ON.and..not.NewOrigin%State%Deposition) &
                            call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1370')


SP:                     if (NewProperty%SedimentPartition%ON) then

                            Me%State%Partition = ON
                            NewOrigin%State%Partition     = ON

                            call GetData(NewProperty%SedimentPartition%Coefficient,      &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='PARTITION_COEF_SED',             &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = 0.98,                            &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                   &
                                call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1380')

                            call GetData(NewProperty%SedimentPartition%TransferRate,     &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='PARTITION_RATE_SED',             &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = 1e-4,                            &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                   &
                                call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1390')


                            call GetData(NewProperty%SedimentPartition%CoupleProp,       &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='PARTITION_COUPLE_SED',           &
                                         ClientModule ='ModuleLagrangian',               &
                                         Default      = 0.,                              &    
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                   &
                                call SetError (FATAL_, KEYWORD_, 'ConstructOrigins - ModuleLagrangian - ERR1400')

                        endif SP

                        !Statistics
                        call GetData(NewProperty%Statistics,                             &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='STATISTICS',                         &
                                     default      = .false.,                             &
                                     ClientModule ='ModuleLagrangian',                   &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1410'
                        
    
                        if (NewProperty%Statistics) then

                            !StatisticsFile
                            call GetData(NewProperty%StatisticsFile,                     &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='STATISTICS_FILE',                &
                                         ClientModule ='ModuleLagrangian',               &
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1420'
                            if (flag /= 1) then
                                write(*,*)'No Statistcis file definition for the Property : ',trim(NewProperty%Name)
                                write(*,*)'Use the KeyWord                                : STATISTICS_FILE'
                                write(*,*)'OR disable STATISTICS'
                                stop      'ConstructOrigins - ModuleLagrangian - ERR96'
                            endif

                            !Do frequency analysis considering the tracers concentration
                            call GetData(NewProperty%StatisticsLag,                      &
                                         Me%ObjEnterData,                                &
                                         flag,                                           &
                                         SearchType   = FromBlockInBlock,                &
                                         keyword      ='STATISTICS_LAG',                 &
                                         ClientModule ='ModuleLagrangian',               &
                                         default      = .false.,                         &
                                         STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1430'


                        endif

                        !This property has an extinction parameter. This parameter can be use 
                        !to compute the effect of this property in the light extinction
                        call GetData(NewProperty%ExtinctionParameter,                    &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='EXTINCTION_PARAMETER',               &
                                     ClientModule ='ModuleLagrangian',                   &
                                     default      = 1.,                                  &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1440'


                        !This property is being filter from the water column
                        call GetData(NewProperty%Filtration,                             &
                                     Me%ObjEnterData,                                    &
                                     flag,                                               &
                                     SearchType   = FromBlockInBlock,                    &
                                     keyword      ='FILTRATION',                         &
                                     ClientModule ='ModuleLagrangian',                   &
                                     default      = .false.,                             &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1450'

                        if (NewProperty%Filtration) NewOrigin%Filtration = .true. 


                        if (PropertyID == Sediment) SedimentDefined = .true.

                        if (NewOrigin%State%FarFieldBuoyancy .or. NewOrigin%State%ComputePlume) then

                            if (PropertyID == Salinity_   ) then
                                if (NewOrigin%State%ComputePlume)                       &
                                    NewOrigin%Movement%JetSalinity    = NewProperty%Concentration
                                SalOK = .true.
                            endif

                            if (PropertyID == Temperature_) then
                                if (NewOrigin%State%ComputePlume)                       &
                                    NewOrigin%Movement%JetTemperature = NewProperty%Concentration
                                TempOK = .true.
                            endif

                        endif

                        !Insert Property into list of properties
                        call InsertPropertyToList(NewOrigin, NewProperty, SetStates = .true.)

                    else

                        exit DOPROP

                    endif

                enddo DOPROP

                if (NewOrigin%Filtration) Me%State%Filtration = .true. 

                if ((NewOrigin%State%FarFieldBuoyancy  .or. Me%State%ComputePlume)      &
                     .and. .not.(TempOK .and. SalOK)) then
                    call SetError (FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR1460')
                endif
                if (NewOrigin%State%Deposition .and. .not. SedimentDefined)              &
                    call SetError (FATAL_, INTERNAL_, 'ConstructOrigins - ModuleLagrangian - ERR1470')



                !If neceassary Starts The WQM module for this origin
                if (NewOrigin%State%WQM) then
                    
                    !WQM Data File
                    call GetData(NewOrigin%WQM_DataFile,                                 &
                                 Me%ObjEnterData,                                        &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='WQM_DATA_FILE',                          &
                                 ClientModule ='ModuleLagrangian',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1480'

                    call StartWaterQuality(NewOrigin%WaterQualityID,                     &
                                           NewOrigin%WQM_DataFile,                       &
                                           STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1490'

                    NewOrigin%NextWQMCompute = Me%ExternalVar%Now
                    call GetDTWQM (NewOrigin%WaterQualityID, DTSecond = NewOrigin%DTWQM, &
                                   STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1500'

                endif

                !If neceassary Starts The Oil module for this origin
                if (NewOrigin%State%Oil) then

                    call GetData(NewOrigin%UseTheoricArea,                               &
                                 Me%ObjEnterData,                             &
                                 flag,                                                   &
                                 SearchType   = FromBlock,                               &
                                 keyword      ='THEORIC_AREA',                           &
                                 ClientModule ='ModuleLagrangian',                       &
                                 Default      = .false.,                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1510'

                    call ConstructParticOil (NewOrigin, ClientNumber)

                endif

                call GetData(NewOrigin%State%Age,                                        &
                             Me%ObjEnterData,                                 &
                             flag,                                                       &
                             SearchType   = FromBlock,                                   &
                             keyword      ='COMPUTE_AGE',                                &
                             ClientModule ='ModuleLagrangian',                           &
                             Default      = .false.,                                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1520'

                if (NewOrigin%State%Age) Me%State%Age = .true.


                if (NewOrigin%State%ComputePlume) then

                    call Construct_Jet(JetID       = NewOrigin%Movement%ObjJet,          &
                                       FileName    = NewOrigin%Movement%JetDataFile,     &
                                       PositionX   = NewOrigin%Position%X,               & 
                                       PositionY   = NewOrigin%Position%Y,               & 
                                       Flow        = NewOrigin%Flow,                     &
                                       Salinity    = NewOrigin%Movement%JetSalinity,     &
                                       Temperature = NewOrigin%Movement%JetTemperature,  &
                                       STAT        = STAT_CALL) 

                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1530'

                    call UnitsManager(NewOrigin%Movement%JetUnit, OPEN_FILE,             &
                                      STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ConstructOrigins - ModuleLagrangian - ERR1540'
                    endif

                    call ReadFileName("ROOT_SRT", RootPath, STAT = STAT_CALL)

                    NewOrigin%Movement%JetFileOut = trim(RootPath)//trim(NewOrigin%Name)//".jet"

                    open (file   = NewOrigin%Movement%JetFileOut,                           &
                          unit   = NewOrigin%Movement%JetUnit,                              & 
                          status = "unknown",                                               &
                          form   = "formatted", IOSTAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) then
                        write(*,*) 'Error opening time series file ',                       &
                                    trim(NewOrigin%Movement%JetFileOut)
                        stop 'ConstructOrigins - ModuleLagrangian - ERR1560'
                    endif

                    write(NewOrigin%Movement%JetUnit,'(A110)')"Year Month Day Hour Minutes "&
                    // "Seconds Dilution X Y Z Density Temperature "                        &
                    // "Salinity U V W MixingHorLength Thickness"
                    write(NewOrigin%Movement%JetUnit,*) '<BeginTimeSerie>'


                    if (NewOrigin%Old) then


                    endif

                endif


                call GetData(NewOrigin%Beaching,                                        &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlock,                                  &
                             keyword      ='BEACHING',                                  &
                             ClientModule ='ModuleLagrangian',                          &
                             Default      = OFF,                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1570'



                call GetData(NewOrigin%AveragePositionON,                               &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlock,                                  &
                             keyword      ='COMPUTE_AVERAGE_POSITION',                  &
                             ClientModule ='ModuleLagrangian',                          &
                             Default      = OFF,                                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1580'

                if (NewOrigin%AveragePositionON) then

                    call GetData(NewOrigin%CoefRadius,                                  &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlock,                              &
                                 keyword      ='COEF_RADIUS',                           &
                                 ClientModule ='ModuleLagrangian',                      &
                                 !90% of the particles inside the circle
                                 Default      = 1.645,                                  &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1590'

                endif
                


                !Insert New Origin to the list of origins    
                call InsertOriginToList (Me%FirstOrigin, NewOrigin, Me%nOrigins)

                if (ClonesExist) then
                    call CheckForOriginClones(OriginalOrigin%Name, ClientNumber,        &
                                              ClientNumberClone, FoundCloneOrigin, ClonesExist)
                endif
 
            else BF
            
                FirstExtraction = .false. 

                exit DW !No more blocks

            end if BF

        end do DW

        end do SB


        !Finished reading block -> unlocks block reading
        call Block_Unlock(Me%ObjEnterData,                                              &
                          ClientNumber,                                                 &
                          STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1600'

        if (Me%ObjEnterDataClone   /= 0)  then
            !Finished reading block -> unlocks block reading
            call Block_Unlock(Me%ObjEnterDataClone,                                     &
                              ClientNumberClone,                                        &
                              STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1610'

            call KillEnterData(Me%ObjEnterDataClone, STAT = STAT_CALL)  

            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1620'
        endif

        if (Me%ObjEnterDataOriginal /= 0) then
            call KillEnterData(Me%ObjEnterDataOriginal, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangian - ERR1630'
        endif

!#ifdef _CGI_

        if (Me%RunOnline) then

            call ChangeOriginOnline(DT, Nmax, NFirstExtraction)
        endif

!#endif

        nullify(OriginalOrigin, NewOrigin)

    end subroutine ConstructOrigins

    !--------------------------------------------------------------------------

    subroutine CheckForOriginClones(OriginName, ClientNumber, ClientNumberClone,        &
                                    FoundCloneOrigin, ClonesExist)


        !Arguments-------------------------------------------------------------
        character(len=*), intent(IN)                :: OriginName
        integer         , intent(IN)                :: ClientNumber
        integer         , intent(OUT)               :: ClientNumberClone
        logical         , intent(OUT)               :: FoundCloneOrigin, ClonesExist

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: OriginNameClone
        logical                                     :: BlockFound
        integer                                     :: STAT_CALL, flag
        integer                                     :: StartLine, EndLine 
        integer, save                               :: PreviousStart


        !Begin-----------------------------------------------------------------

        !first time check if exist "clone" blocks
i1:     if (Me%ObjEnterDataClone == 0) then

            !Construct enter data "clone"
            call ConstructEnterData(Me%ObjEnterDataClone, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR10'

            call ExtractBlockFromBuffer(Me%ObjEnterDataClone,                           &
                                        ClientNumberClone,                              &
                                        block_begin_Clone, block_end_Clone,             &
                                        BlockFound,                                     &
                                        STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR20'

i2:         if (BlockFound) then
    
                !Construct enter data original again
                call ConstructEnterData(Me%ObjEnterDataOriginal, Me%Files%ConstructData, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR30'

                PreviousStart  = FillValueInt

                call RewindBuffer(Me%ObjEnterDataClone,                                 &
                                  STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR40'


            else  i2

                ClonesExist = .false. 

                !Finished reading block -> unlocks block reading
                call Block_Unlock(Me%ObjEnterDataClone,                                 &
                                  ClientNumberClone,                                    &
                                  STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR50'

                call KillEnterData(Me%ObjEnterDataClone, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR60'

            endif i2

        endif i1

i3:     if (ClonesExist) then

d2:     do
 
            call ExtractBlockFromBuffer(Me%ObjEnterDataClone,                           &
                                        ClientNumberClone,                              &
                                        block_begin_Clone, block_end_Clone,             &
                                        BlockFound,                                     &
                                        STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR70'


i4:         if (BlockFound) then 

                call GetBlockSize(Me%ObjEnterDataClone, ClientNumberClone, StartLine, EndLine, STAT= STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR80'

i5:             if (StartLine > PreviousStart) then

                    call GetData(OriginNameClone,                                       &
                                 Me%ObjEnterDataClone,                                  &
                                 flag,                                                  &
                                 SearchType   = FromBlock,                              &
                                 keyword      ='CLONE',                                 &
                                 ClientModule ='ModuleLagrangian',                      &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR90'

                    if (flag == 0) stop 'CheckForOriginClones - ModuleLagrangian - ERR100'


i6:                if (trim(OriginNameClone) == trim(OriginName)) then
                        FoundCloneOrigin = .true.

                        call ReplaceKeywordsInClone(ClientNumber, StartLine, EndLine)

                        PreviousStart = StartLine

                        exit

                    endif i6

                endif i5

            else i4

                PreviousStart  = FillValueInt

                call RewindBuffer(Me%ObjEnterDataClone, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangian - ERR110'

                FoundCloneOrigin = .false.

                exit

            endif i4

        enddo d2

        endif i3

    end subroutine CheckForOriginClones

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine ReplaceKeywordsInClone(ClientNumber, StartClone, EndClone)


        !Arguments-------------------------------------------------------------
        integer                                     :: ClientNumber, StartClone, EndClone

        !Local-----------------------------------------------------------------
        character(len=PathLength)                   :: FullBufferLine
        character(len=StringLength)                 :: KeywordClone, KeywordOriginal
        logical                                     :: ExistKeywordClone, ExistKeywordOriginal
        integer                                     :: StartLine, EndLine
        integer                                     :: lineC, lineO
        integer                                     :: STAT_CALL


        !Begin-----------------------------------------------------------------

        call GetBlockSize(Me%ObjEnterData, ClientNumber, StartLine, EndLine, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangian - ERR10'

    

d1:     do lineC = StartClone, EndClone

            call GetKeywordFromLine(Me%ObjEnterDataClone, lineC, ExistKeywordClone, KeywordClone, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangian - ERR20'

i1:         if (ExistKeywordClone) then
                             
i2:             if (trim(KeywordClone) /= 'CLONE') then

d2:                 do lineO = StartLine, EndLine                

                        call GetKeywordFromLine(Me%ObjEnterDataOriginal, lineO, ExistKeywordOriginal, &
                                                KeywordOriginal, STAT = STAT_CALL)  
                        if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangian - ERR30'
                    
i3:                     if (ExistKeywordOriginal) then

i4:                         if (trim(KeywordClone) == trim(KeywordOriginal)) then

                                call GetFullBufferLine    (Me%ObjEnterDataClone, lineC, FullBufferLine, STAT = STAT_CALL)  
                                if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangian - ERR40'

                                call ReplaceFullBufferLine(Me%ObjEnterData,      lineO, FullBufferLine, STAT = STAT_CALL)  
                                if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangian - ERR50'

                                exit

                            endif i4
                        endif i3
                    enddo d2
                endif i2
            endif i1

        enddo d1


    end subroutine ReplaceKeywordsInClone

    !--------------------------------------------------------------------------

!#ifdef _CGI_

    !--------------------------------------------------------------------------
    subroutine ReadQueryString(Nmax)


        !Arguments-------------------------------------------------------------
        integer                                     :: Nmax

        !Local-----------------------------------------------------------------
        real,  dimension(:), pointer                :: AuxReal
        !character(Len=1024)                         :: Online_St, Online_St2
        integer                                     :: i !, j, k, no
        !integer                                     :: STAT_CALL
        !integer(2)                                  :: STAT2


        !Begin-----------------------------------------------------------------

        !if ('TIME_STAMP'/= Online_St(1:10))  stop 'Lagrangian - ReadQueryString - ERR20'

        !j = Scan(Online_St,'&')

        !1 - Year, 2 - Month, 3 - Day, 4 - Hours, 5 - Minutes, 6 - seconds, 7 - miliseconds
        !Me%Online%TimeStamp=Online_St(12:j-1)

        call GetDataOnlineString ('TIME_STAMP', CharData = Me%Online%TimeStamp)

        !Online_St=trim(Online_St(j+1:Len_Trim(Online_St)))

        !if ('LAG_NORIGINS' /= Online_St(1:12)) stop 'Lagrangian - ReadQueryString - ERR30'

        !j = Scan(Online_St,'&')
        !Online_St2=Online_St(14:j-1)
       
        !read(Online_St2,*,IOSTAT=STAT_CALL) Nmax
        !if (STAT_CALL /= SUCCESS_)  stop 'Lagrangian - ReadQueryString - ERR40'

        call GetDataOnlineString ('LAG_NORIGINS', IntData = Nmax)

        allocate(Me%Online%X(1:Nmax), Me%Online%Y(1:Nmax))
        allocate(Me%Online%StartDate(1:Nmax,6))

        allocate(AuxReal(1:2*Nmax))

        !Online_St=trim(Online_St(j+1:Len_Trim(Online_St)))

        !if ('LAG_XY'/= Online_St(1:6))  stop 'Lagrangian - ReadQueryString - ERR50'

        !Online_St2=' '
        !j = Scan(Online_St,'&')

        !Online_St2=Online_St(8:j-1)

        !do i=1,Len_Trim(Online_St2)
        !    if (Online_St2(i:i) =='_') Online_St2(i:i) = ' '
        !enddo

        !read(Online_St2,*,IOSTAT=STAT_CALL) (Me%Online%X(i),Me%Online%Y(i),i=1,Nmax)
        !if (STAT_CALL /= SUCCESS_)  stop 'Lagrangian - ReadQueryString - ERR60'

        call GetDataOnlineString ('LAG_XY', ArrayData = AuxReal)

        do i = 1, Nmax
            Me%Online%X(i) = AuxReal(2*i-1)
            Me%Online%Y(i) = AuxReal(2*i  )
        enddo

        deallocate(AuxReal)

        allocate(AuxReal(1:6*Nmax))

        call GetDataOnlineString ('LAG_START', ArrayData = AuxReal)
         
        !Online_St=trim(Online_St(j+1:Len_Trim(Online_St)))

        !if ('LAG_START'/= Online_St(1:9)) stop 'Lagrangian - ReadQueryString - ERR70'

        !j = Scan(Online_St,'&')

        !Online_St2=Online_St(11:j-1)

        
        do i=1,Nmax
            !1 - Year, 2 - Month, 3 - Day, 4 - Hours, 5 - Minutes, 6 - seconds, 7 - miliseconds
            !do k = 1, 6
                !no = scan(Online_St2,'_') 

                !if (no==0) no = len_trim(Online_St2)+1

                !read(Online_St2(1:no-1),*,IOSTAT=STAT_CALL) Me%Online%StartDate(i,k)

                !if (STAT_CALL /= SUCCESS_)  stop 'Lagrangian - ReadQueryString - ERR80'

                !if (no/=0) Online_St2 = trim(Online_St2(no+1:len_trim(Online_St2)))
            !enddo
            Me%Online%StartDate(i,1:6) = AuxReal(i*6-5:i*6)
        enddo

        deallocate(AuxReal)
        

        !Online_St=trim(Online_St(j+1:Len_Trim(Online_St)))

        !if ('WIND_COEF'/= Online_St(1:9))  stop 'Lagrangian - ReadQueryString - ERR90'

        !j = Scan(Online_St,'&')

        !Online_St2=Online_St(11:j-1)

        !do i=1,Len_Trim(Online_St2)
        !    if (Online_St2(i:i) =='_') Online_St2(i:i) = ' '
        !enddo

        allocate(Me%Online%WindCoef(1:Nmax))

        call GetDataOnlineString ('WIND_COEF', ArrayData = Me%Online%WindCoef)

        !read(Online_St2,*,IOSTAT=STAT_CALL) (Me%Online%WindCoef(i),i=1,Nmax)

        !if (STAT_CALL /= SUCCESS_)  stop 'Lagrangian - ReadQueryString - ERR100'

        !Online_St=trim(Online_St(j+1:Len_Trim(Online_St)))

        !Me%Online%ON            = .true. 

    end subroutine ReadQueryString

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ChangeOriginOnline(DT, Nmax, NFirstExtraction)


        !Arguments-------------------------------------------------------------
        real                                        :: DT
        integer                                     :: Nmax, NFirstExtraction

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: NewOrigin
        character(Len=PathLength)                   :: String
        real                                        :: Xorig, Yorig, GridAngle, Xaux, Yaux
        integer                                     :: i, j


        !Begin-----------------------------------------------------------------

!        if (Me%Online%ON) then

        Me%Files%TransientHDF =trim(Me%Files%TransientHDF)//'_'//trim(Me%Online%TimeStamp)//'.hdf'

        NewOrigin=> Me%FirstOrigin            
        
        do i=1, Nmax
        do j=1, NFirstExtraction

            !Gets Origin
            call GetGridOrigin  (Me%ObjHorizontalGrid, Xorig, Yorig)

            !Gets GridAngle
            call GetGridAngle   (Me%ObjHorizontalGrid, GridAngle)

            Xaux = Me%Online%X(i)
            Yaux = Me%Online%Y(i)

            call RODAXY(-Xorig, -Yorig,  0.0, Xaux, Yaux)  
            call RODAXY(0.0, 0.0, -GridAngle, Xaux, Yaux)  

            NewOrigin%Position%X = Xaux
            NewOrigin%Position%Y = Yaux

            call Convert_Coord1D_CellIJ (NewOrigin%Position,            &
                                         Me%ExternalVar%XX, Me%ExternalVar%YY)

            call Convert_CellIJ_XY(NewOrigin%Position)
            call Convert_CellIJ_IJ(NewOrigin%Position)  

            !Change Origin name
            write(String,*) i
            NewOrigin%Name = trim(NewOrigin%Name)//'_'//trim(adjustl(String))

            Me%Online%StartDate(i,6) = Me%Online%StartDate(i,6)

            call SetDate (NewOrigin%StartEmission, Me%Online%StartDate(i,1), Me%Online%StartDate(i,2), &
                                                   Me%Online%StartDate(i,3), Me%Online%StartDate(i,4), &
                                                   Me%Online%StartDate(i,5), Me%Online%StartDate(i,6))


            NewOrigin%StopEmission = NewOrigin%StartEmission + DT

            NewOrigin%Movement%WindTransferCoef = Me%Online%WindCoef(i)

            NewOrigin => NewOrigin%Next

        enddo
        enddo

!        endif

        deallocate(Me%Online%X, Me%Online%Y)
        deallocate(Me%Online%StartDate     )
        deallocate(Me%Online%WindCoef      )

    end subroutine ChangeOriginOnline

    !--------------------------------------------------------------------------
!#endif
    !--------------------------------------------------------------------------
    subroutine ConstructEmissionType(NewOrigin, EmissionSpatial, EmissionTemporal)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: NewOrigin
        character(len=*)                            :: EmissionSpatial
        character(len=*)                            :: EmissionTemporal

        !Local-----------------------------------------------------------------


        !Spatial Emission type
        select case (trim(adjustl(EmissionSpatial)))

        case (trim(adjustl(Char_Point)))

            NewOrigin%EmissionSpatial = Point_

        case (trim(adjustl(Char_Box)))

            NewOrigin%EmissionSpatial = Box_

        case (trim(adjustl(Char_Accident)))

            NewOrigin%EmissionSpatial = Accident_

        case default
            
            write(*,*)'Unknown Spatial Emission type ', trim(adjustl(EmissionSpatial))
            write(*,*)'Valid option [1] :',trim(adjustl(Char_Point))
            write(*,*)'Valid option [2] :',trim(adjustl(Char_Box))
            write(*,*)'Valid option [3] :',trim(adjustl(Char_Accident))
            write(*,*)'Origin           :',trim(adjustl(NewOrigin%Name))
            stop 'ConstructEmissionType - ModuleLagrangian - ERR01'

        end select


        !Temporal Emission type
        select case (trim(adjustl(EmissionTemporal)))

        case (trim(adjustl(Char_Continuous)))

            NewOrigin%EmissionTemporal = Continuous_

        case (trim(adjustl(Char_Instantaneous)))

            NewOrigin%EmissionTemporal = Instantaneous_

        case default

            write(*,*)'Unknown Temporal Emission type ', trim(adjustl(EmissionTemporal))
            write(*,*)'Valid option [1] :',trim(adjustl(Char_Continuous))
            write(*,*)'Valid option [2] :',trim(adjustl(Char_Instantaneous))
            write(*,*)'Origin           :',trim(adjustl(NewOrigin%Name))
            stop 'ConstructEmissionType - ModuleLagrangian - ERR01'

        end select

        
        !Verifies consistency
        if (NewOrigin%EmissionTemporal == Continuous_ .and. NewOrigin%EmissionSpatial == Accident_) then
            write(*,*)'Cant emit particle in continous as accident'
            write(*,*)'Origin :',trim(adjustl(NewOrigin%Name))
            stop 'ConstructEmissionType - ModuleLagrangian - ERR02'
        endif

        !If BoxDif isnt associated keyword Partic PARTIC_BOX wasnt found
        if (NewOrigin%EmissionSpatial == Box_ .and. Me%ObjBoxDif == 0) then
            write(*,*)'Cant emit particle in a box, once keyword PARTIC_BOX wasnt given'
            write(*,*)'Origin :',trim(adjustl(NewOrigin%Name))
            stop 'ConstructEmissionType - ModuleLagrangian - ERR03'
        endif


    end subroutine ConstructEmissionType

    !--------------------------------------------------------------------------

    subroutine ConstructEmission 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            !Constructs the Temporal emission 
            select case (CurrentOrigin%EmissionTemporal)

            case (Continuous_)

                call ConstructEmissionTime          (CurrentOrigin)

                if (Me%ExternalVar%Now >= CurrentOrigin%NextEmission) then

                    !Realizes the first emission
                    select case (CurrentOrigin%EmissionSpatial)

                    case (Box_)
            
                        call EmissionBox           (CurrentOrigin)

                    case (Point_)

                        call EmissionPoint         (CurrentOrigin)

                    ! The accident can not be continuous YET!!!!!
!                    case (Accident_)

!                        call EmissionAccident      (CurrentOrigin)

                    end select

                endif

            case (Instantaneous_)

                if (.not. CurrentOrigin%Old) then

                    !Realizes the first emission
                    select case (CurrentOrigin%EmissionSpatial)

                    case (Box_)
            
                        call EmissionBox           (CurrentOrigin)

                    case (Point_)

                        call EmissionPoint         (CurrentOrigin)

                    case (Accident_)

                        if (CurrentOrigin%AccidentTime ==                                &
                            Me%ExternalVar%BeginTime)                         &
                            call EmissionAccident  (CurrentOrigin)

                    end select

                endif

            end select
    
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine ConstructEmission

    !--------------------------------------------------------------------------

    subroutine ConstructOverlay 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: SizeILB, SizeIUB
        integer                                     :: SizeJLB, SizeJUB
        integer                                     :: SizeKLB, SizeKUB

        !Shorten Size
        SizeILB     = Me%ExternalVar%Size%ILB   
        SizeIUB     = Me%ExternalVar%Size%IUB   
        SizeJLB     = Me%ExternalVar%Size%JLB   
        SizeJUB     = Me%ExternalVar%Size%JUB   
        SizeKLB     = Me%ExternalVar%Size%KLB   
        SizeKUB     = Me%ExternalVar%Size%KUB   
        
        if (Me%OverLay%Overlay) then

            allocate (Me%OverLay%VelUFinal (SizeILB:SizeIUB, SizeJLB:SizeJUB, SizeKLB:SizeKUB))
            allocate (Me%OverLay%VelVFinal (SizeILB:SizeIUB, SizeJLB:SizeJUB, SizeKLB:SizeKUB))

        else
        
            nullify (Me%OverLay%VelUFinal)
            nullify (Me%OverLay%VelVFinal)

        endif


    end subroutine ConstructOverlay

    !--------------------------------------------------------------------------

    subroutine ConstructMonitoring 

        !Local-----------------------------------------------------------------
        integer                                         :: NumberOfBoxes, nB
        integer                                         :: NumberOfOrigins, iO
        integer                                         :: STAT_CALL
        character(StringLength), dimension(:), pointer  :: PropertyList
        character(StringLength)                         :: AuxChar
        type (T_Property), pointer                      :: CurrentProperty
        integer                                         :: iProp
        integer                                         :: ILB, IUB
        integer                                         :: JLB, JUB
        integer                                         :: KLB, KUB
 
        !Begin-----------------------------------------------------------------
        
        if (Me%State%Monitor) then

            call GetNumberOfBoxes(Me%ObjMonBox, NumberOfBoxes3D = NumberOfBoxes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMonitoring - ModuleLagrangian - ERR01'

            NumberOfOrigins = Me%nOrigins

            allocate (Me%Monitor%SurfaceBoxVolume           (NumberOfBoxes                  ))
            allocate (Me%Monitor%InstBoxVolume              (NumberOfBoxes                  ))
            allocate (Me%Monitor%InstVolumeByOrigin         (NumberOfBoxes, NumberOfOrigins ))
            allocate (Me%Monitor%IntgBoxVolume              (NumberOfBoxes                  ))
            allocate (Me%Monitor%IntgVolumeByOrigin         (NumberOfBoxes, NumberOfOrigins ))
            allocate (Me%Monitor%InstBoxMass                (NumberOfBoxes                  ))
            allocate (Me%Monitor%InstMassByOrigin           (NumberOfBoxes, NumberOfOrigins ))
            allocate (Me%Monitor%NumberOfCellsPerBox        (NumberOfBoxes                  ))
            allocate (Me%Monitor%NumberOfCellsFromOrigin    (NumberOfBoxes, NumberOfOrigins ))
            allocate (Me%Monitor%ObjTimeSerie               (NumberOfBoxes                  ))            

            allocate (Me%Monitor%InstBoxLogMass             (NumberOfBoxes                  ))
            allocate (Me%Monitor%InstBoxConc                (NumberOfBoxes                  ))
            allocate (Me%Monitor%NumberOfTracers            (NumberOfBoxes                  ))
            allocate (Me%Monitor%InstLogMassByOrigin        (NumberOfBoxes, NumberOfOrigins ))
            allocate (Me%Monitor%NumberOfTracersFromOrigin  (NumberOfBoxes, NumberOfOrigins ))
            allocate (Me%Monitor%InstBoxMassFractionByOrigin(NumberOfBoxes, NumberOfOrigins ))
            allocate (Me%Monitor%NbrBoxContaminatedTracers  (NumberOfBoxes                  ))
            allocate (Me%Monitor%VolBoxContaminatedTracers  (NumberOfBoxes                  ))
            allocate (Me%Monitor%AverageBoxContaminatedConc (NumberOfBoxes                  ))
            allocate (Me%Monitor%ContaminationProbability   (NumberOfBoxes                  ))

            Me%Monitor%IntgBoxVolume      = 0.
            Me%Monitor%IntgVolumeByOrigin = 0.
            Me%Monitor%ObjTimeSerie       = 0



            If (Me%State%MonitorPropMass) then
                allocate (PropertyList (2 + 4*(NumberOfOrigins + 1)))
            Else
                allocate (PropertyList (2*(NumberOfOrigins + 1)))
            End If

            do nB = 1, NumberOfBoxes

                PropertyList (1) = "Inst. Volume"
                do iO = 1, NumberOfOrigins
                    write(AuxChar, fmt=*)iO
                    PropertyList (iO+1) = "Inst. VolFrOrigin_"//trim(adjustl(AuxChar))
                enddo

                PropertyList ((NumberOfOrigins + 1) + 1) = "Intg. Volume"
                do iO = 1, NumberOfOrigins
                    write(AuxChar, fmt=*)iO
                    PropertyList ((NumberOfOrigins + 1) + 1 + iO) = "Intg. VolFrOrigin_"//trim(adjustl(AuxChar))
                enddo

                If (Me%State%MonitorPropMass) then

                    PropertyList (2 * (NumberOfOrigins + 1) + 1) = "Inst. Mass"
                    do iO = 1, NumberOfOrigins
                        write(AuxChar, fmt=*)iO
                        PropertyList (2 * (NumberOfOrigins + 1) + 1 + iO) = "Inst. MassFrOrigin_"//trim(adjustl(AuxChar))
                    enddo
                    
                    if (Me%Monitor%MassFractionType .EQ. Arithmetic) then
                        PropertyList (3 * (NumberOfOrigins + 1) + 1) = "Inst. Arit. Conc."
                        PropertyList (3 * (NumberOfOrigins + 1) + 2) = "Cont. Arit. Conc."
                    else
                        PropertyList (3 * (NumberOfOrigins + 1) + 1) = "Inst. Geom. Conc."
                        PropertyList (3 * (NumberOfOrigins + 1) + 2) = "Cont. Geom. Conc."
                    end if
                        PropertyList (3 * (NumberOfOrigins + 1) + 3) = "Contam. Prob."

                    do iO = 1, NumberOfOrigins
                        write(AuxChar, fmt=*)iO
                        PropertyList (3 * (NumberOfOrigins + 1) + 3 + iO) = "Inst. RelMassFrOrigin_"//trim(adjustl(AuxChar))
                    enddo


                End If

                write(AuxChar, fmt=*)nB
                call StartTimeSerie (Me%Monitor%ObjTimeSerie(nB),                           &
                                     Me%ObjTime,                                            &
                                     Me%Files%ConstructData,                                &
                                     PropertyList,                                          &
                                     Extension   = "LMB",                                   &
                                     ResultFileName = "MonitorBox_"//trim(adjustl(AuxChar)),&
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMonitoring - ModuleLagrangian - ERR02'

            enddo
            deallocate (PropertyList)

        else

            nullify  (Me%Monitor%SurfaceBoxVolume           )
            nullify  (Me%Monitor%InstBoxVolume              )
            nullify  (Me%Monitor%InstVolumeByOrigin         )
            nullify  (Me%Monitor%IntgBoxVolume              )
            nullify  (Me%Monitor%IntgVolumeByOrigin         )
            nullify  (Me%Monitor%InstBoxMass                )
            nullify  (Me%Monitor%InstMassByOrigin           )
            nullify  (Me%Monitor%NumberOfCellsPerBox        )       
            nullify  (Me%Monitor%NumberOfCellsFromOrigin    ) 
            nullify  (Me%Monitor%ObjTimeSerie               )
            nullify  (Me%Monitor%InstLogMassByOrigin        )
            nullify  (Me%Monitor%NumberOfTracersFromOrigin  )
            nullify  (Me%Monitor%InstBoxMassFractionByOrigin)
            nullify  (Me%Monitor%InstBoxLogMass             )         
            nullify  (Me%Monitor%InstBoxConc                )
            nullify  (Me%Monitor%NumberOfTracers            )
            nullify  (Me%Monitor%InstBoxMassFractionByOrigin)
            nullify  (Me%Monitor%InstLogMassByOrigin        )
            nullify  (Me%Monitor%NumberOfTracersFromOrigin  )
            nullify  (Me%Monitor%ContaminationProbability   )    
            nullify  (Me%Monitor%AverageBoxContaminatedConc )
            nullify  (Me%Monitor%NbrBoxContaminatedTracers  )
            nullify  (Me%Monitor%VolBoxContaminatedTracers  )
        endif

        if (Me%State%EulerianMonitor) then

            ILB     = Me%ExternalVar%Size%ILB   
            IUB     = Me%ExternalVar%Size%IUB   
            JLB     = Me%ExternalVar%Size%JLB   
            JUB     = Me%ExternalVar%Size%JUB   
            KLB     = Me%ExternalVar%Size%KLB   
            KUB     = Me%ExternalVar%Size%KUB   

            !This test is done for simply reason
            if (Me%nGroups > 1) then
                write(*,*)'Cannot write Eulerian Box Time Series for simulation with more then one Group'
                stop 'ConstructMonitoring - ModuleLagrangian - ERR40'
            endif

            !Counts Number of properties
            iProp = 0
            CurrentProperty => Me%FirstOrigin%FirstProperty
            do while (associated(CurrentProperty))
                iProp = iProp + 1    
                CurrentProperty => CurrentProperty%Next
            enddo

            allocate (PropertyList(1:iProp))

            !Counts Number of properties
            iProp = 0
            CurrentProperty => Me%FirstOrigin%FirstProperty
            do while (associated(CurrentProperty))
                iProp               = iProp + 1
                PropertyList(iProp) = "Lagrangian_"//trim(CurrentProperty%Name)
                CurrentProperty => CurrentProperty%Next
            enddo

            allocate(Me%EulerianMonitor%Mass(ILB:IUB, JLB:JUB, KLB:KUB))

            call SetMatrixValue(Me%EulerianMonitor%Mass, Me%ExternalVar%Size, &
                                dble(0.), Me%ExternalVar%Waterpoints3D)
            
            !Starts BoxDif
            call StartBoxDif(BoxDifID         = Me%ObjEulMonBox,                        &
                             TimeID           = Me%ObjTime,                             &
                             HorizontalGridID = Me%ObjHorizontalGrid,                   &
                             BoxesFilePath    = Me%Files%EulerianMonitorBox,            &
                             ScalarOutputList = PropertyList,                           &
                             WaterPoints3D    = Me%ExternalVar%Waterpoints3D,           &
                             Size3D           = Me%ExternalVar%Size,                    &
                             WorkSize3D       = Me%ExternalVar%WorkSize,                &
                             STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMonitoring - ModuleLagrangian - ERR050'
            
            deallocate (PropertyList)

        end if


    end subroutine ConstructMonitoring

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len = StringLength), dimension(:), pointer    :: PropertyList
        integer, pointer, dimension(:, :, :)                    :: WaterPoints3D
        type (T_Property), pointer                              :: CurrentProperty
        integer                                                 :: iProp
        integer                                                 :: STAT_CALL
        integer                                                 :: iflag
        character(len=PathLength)                               :: TimeSerieLocationFile


        !This test is done for simply reason
        if (Me%nGroups > 1) then
            write(*,*)'Cannot write Time Series for simulation with more then one Group'
            stop 'ConstructTimeSeries - ModuleLagrangian - ERR01'
        endif

        !Counts Number of properties
        iProp = 0
        CurrentProperty => Me%FirstOrigin%FirstProperty
        do while (associated(CurrentProperty))
            iProp = iProp + 1    

            CurrentProperty => CurrentProperty%Next
        enddo


        !Allocates Property List
        if (iProp > 0) then
            allocate(PropertyList(iProp))
        else
            Me%WritesTimeSerie = .false.
            write(*,*)'No Properties defined'
            write(*,*)'Particle Time Series Disabled'
            return
        endif
            

        !Fills Property List
        iProp = 0
        CurrentProperty => Me%FirstOrigin%FirstProperty
        do while (associated(CurrentProperty))
            iProp = iProp + 1    
            PropertyList(iProp) = trim(adjustl(CurrentProperty%Name))
            CurrentProperty => CurrentProperty%Next
        enddo


        !Gets the position of the water points in the Map Module
        call GetWaterPoints3D(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangian - ERR03'
        
        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData ,iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleLagrangian',                                 &
                     Default      = Me%Files%ConstructData,                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ConstructTimeSerie - ModuleLagrangian - ERR04' 

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                 &
                            TimeSerieLocationFile,                                       &
                            PropertyList, "srl",                                         &
                            WaterPoints3D = WaterPoints3D,                               &
                            STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleLagrangian - ERR05'

        !Ungets the WaterPoints 
        call UnGetMap (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangian - ERR06'

        deallocate (PropertyList)


    end subroutine ConstructTimeSeries

    !--------------------------------------------------------------------------

    subroutine VerifyOriginProperties 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Origin), pointer                    :: FirstInGroup, NextInGroup
        logical                                     :: Equal
        integer, dimension(:), allocatable          :: AuxGroups, AuxGroupsIDs
        integer                                     :: iO, i, ILB, IUB, JLB, JUB, KLB, KUB
        logical                                     :: HaveID
          
        !Checks how many different Groups exists
        Me%nGroups = 0
        allocate (AuxGroups    (Me%nOrigins))
        allocate (AuxGroupsIDs (Me%nOrigins))
        CurrentOrigin => Me%FirstOrigin
        iO            =  1
        do while (associated(CurrentOrigin))
        
            HaveID = .false.
            do i = 1, iO - 1
                if (AuxGroups(i) == CurrentOrigin%GroupID) then
                    HaveID = .true.
                endif
            enddo

            if (.not. HaveID) then
                Me%nGroups               = Me%nGroups + 1
                AuxGroupsIDs(Me%nGroups) = CurrentOrigin%GroupID
            endif

            AuxGroups(iO) =  CurrentOrigin%GroupID
            CurrentOrigin => CurrentOrigin%Next
            iO            =  iO + 1

        enddo

        ILB = Me%ExternalVar%Size%ILB
        IUB = Me%ExternalVar%Size%IUB

        JLB = Me%ExternalVar%Size%JLB
        JUB = Me%ExternalVar%Size%JUB

        KLB = Me%ExternalVar%Size%KLB
        KUB = Me%ExternalVar%Size%KUB

        if (Me%State%Deposition)  then
        
            allocate   (Me%TauErosionGrid (Me%nGroups, ILB:IUB, JLB:JUB))
            allocate   (Me%MassSedGrid    (Me%nGroups, ILB:IUB, JLB:JUB))

            Me%TauErosionGrid (:,:,:) = 0.
            Me%MassSedGrid    (:,:,:) = 0.

        endif

        if (Me%State%Filtration)  then

            allocate(Me%Filtration%RelativeMassFilter(ILB:IUB,JLB:JUB,KLB:KUB))

            allocate(Me%Filtration%MassFiltered(ILB:IUB,JLB:JUB,KLB:KUB))

            Me%Filtration%RelativeMassFilter(:,:,:) = 0. 
            Me%Filtration%MassFiltered      (:,:,:) = 0. 
            
            Me%Filtration%Next = Me%ExternalVar%Now

        endif


        !Origins per Group
        allocate   (Me%nOriginsGroup  (Me%nGroups))
        Me%nOriginsGroup = 0



        !Group IDs
        allocate   (Me%GroupIDs       (Me%nGroups))
        do i = 1, Me%nGroups
            Me%GroupIDs (i) = AuxGroupsIDs(i)
        enddo

        deallocate (AuxGroups   )
        deallocate (AuxGroupsIDs)


        !Checks if origins within the same group have the same properties
        do i = 1, Me%nGroups
        
            nullify (FirstInGroup, NextInGroup)
            CurrentOrigin => Me%FirstOrigin
            do while (associated (CurrentOrigin))

                if (CurrentOrigin%GroupID == Me%GroupIDs(i)) then

                    Me%nOriginsGroup (i) = Me%nOriginsGroup (i) + 1
                
                    if (.not. associated(FirstInGroup)) then
                        FirstInGroup => CurrentOrigin
                    else
                        NextInGroup  => CurrentOrigin
                        call VerifyPropertyList (FirstInGroup, NextInGroup, Equal)
                        if (.not. Equal) then
                            
                            write (*,*)'Origins in the same Group but with different properties'
                            write (*,*)'First   Origin :', trim(adjustl(FirstInGroup%Name))
                            write (*,*)'Current Origin :', trim(adjustl(NextInGroup%Name))
                            write (*,*)'Group ID       :', CurrentOrigin%GroupID
                            stop 'VerifyOriginProperties - ModuleLagrangian - ERR01'
                        endif
                    endif

                endif

                CurrentOrigin => CurrentOrigin%Next
            enddo

        enddo

    end subroutine VerifyOriginProperties

    !--------------------------------------------------------------------------

    subroutine MergeOldWithNewOrigins

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Origin), pointer                    :: CurrentOldOrigin
        logical                                     :: Equal

        CurrentOrigin => Me%FirstOrigin
        do while (associated(CurrentOrigin))

            if (CurrentOrigin%Old) then
            
                CurrentOldOrigin => Me%FirstOldOrigin
OldOrigin:      do while (associated(CurrentOldOrigin))

                    if (trim (adjustl(CurrentOldOrigin%Name)) ==                         &
                        trim (adjustl(CurrentOrigin%Name   ))) then

                        !Verifies Property List
                        call VerifyPropertyList (CurrentOldOrigin, CurrentOrigin, Equal)
                        if (.not. Equal) then
                            write (*,*) 'The Properties from the Old and Current Origin are different'
                            write (*,*) 'Old Origin Name: ',trim(CurrentOldOrigin%Name)
                            stop 'MergeOldWithNewOrigins - ModuleLagrangian - ERR02'
                        endif

                        CurrentOrigin%nParticle   =  CurrentOldOrigin%nParticle
                        CurrentOrigin%FirstPartic => CurrentOldOrigin%FirstPartic
                        CurrentOrigin%ObjOil      =  CurrentOldOrigin%ObjOil


                        nullify (CurrentOldOrigin%FirstPartic)

                        exit OldOrigin
                    endif

                    CurrentOldOrigin => CurrentOldOrigin%Next
                enddo OldOrigin

            endif
            
            CurrentOrigin => CurrentOrigin%Next
        enddo

        call DeallocateOriginList (Me%FirstOldOrigin, Me%nOldOrigins)

    end subroutine MergeOldWithNewOrigins

    !--------------------------------------------------------------------------

    subroutine VerifyPropertyList (OriginOne, OriginTwo, Equal)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: OriginOne
        type (T_Origin), pointer                    :: OriginTwo
        logical                                     :: Equal

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyOne
        type (T_Property), pointer                  :: PropertyTwo

        !Same amount of properties
        if (OriginOne%nProperties /= OriginTwo%nProperties) then
            Equal = .false.
            return
        endif

        PropertyOne => OriginOne%FirstProperty
        PropertyTwo => OriginTwo%FirstProperty
        do while (associated(PropertyOne) .and. associated(PropertyTwo))
            
            if (trim(PropertyOne%Name) /= trim(PropertyTwo%Name)) then
                Equal = .false.
                return
            endif

            PropertyOne => PropertyOne%Next
            PropertyTwo => PropertyTwo%Next

        enddo

        Equal = .true.

    end subroutine VerifyPropertyList

    !--------------------------------------------------------------------------

    subroutine EmissionBox (CurrentOrigin)
        
        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Local-----------------------------------------------------------------
        real                                        :: FACTOR
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k
        real                                        :: VOLCEL
        real                                        :: DPPX, DPPY, DPPZ
        integer                                     :: KPP, LP
        type (T_Partic), pointer                    :: NewParticle
        real                                        :: XI, YI, ZI
        real                                        :: auxX, auxY, auxZ, BoxThickness, DepthPartic
        integer, pointer, dimension(:,:,:)          :: Boxes
        integer                                     :: STAT_CALL, BoxCell

        
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB

        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB

        KLB = Me%ExternalVar%WorkSize%KLB
        KUB = Me%ExternalVar%WorkSize%KUB


        !Get the boxes
        call GetBoxes(Me%ObjBoxDif, Boxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'EmissionBox - ModuleLagrangian - ERR01'

                       
        !Distribui Floating particle
FLOAT:  if (CurrentOrigin%Movement%Float .or. CurrentOrigin%Deposition%BottomEmission) then

            FACTOR = 0.0

            if (CurrentOrigin%Movement%Float) then

                BoxThickness = CurrentOrigin%Movement%ThicknessMeters

            else if (CurrentOrigin%Deposition%BottomEmission) then

                BoxThickness = CurrentOrigin%Deposition%BottomDistance

            endif

            do j = JLB, JUB, CurrentOrigin%INCRP 
            do i = ILB, IUB, CurrentOrigin%INCRP 

OP1:            if (Me%ExternalVar%OpenPoints3D(i, j, KUB) == OpenPoint) then

                if (CurrentOrigin%Movement%Float) then
                
                    BoxCell     = Boxes(i, j, KUB)

                    DepthPartic = Me%ExternalVar%SZZ(i, j, KUB)

                else if (CurrentOrigin%Deposition%BottomEmission) then

                    BoxCell     = Boxes(i, j, Me%ExternalVar%KFloor(i, j))


                    DepthPartic = Me%ExternalVar%SZZ(i, j, Me%ExternalVar%KFloor(i, j)-1) +               &
                                               CurrentOrigin%Deposition%BottomDistance / 2.

                endif

BOX2D:          if (BoxCell == CurrentOrigin%BoxNumber) then

                if (CurrentOrigin%ParticleBoxVolume > 0.0) then

                    VOLCEL = Me%Grid%AreaCell(i, j) *                         &
                             BoxThickness + FACTOR

                    KPP    = INT(VOLCEL / CurrentOrigin%ParticleBoxVolume)

                    FACTOR = VOLCEL - CurrentOrigin%ParticleBoxVolume * KPP

                else
                                              
                    KPP = 1                         

                end if

                XI  = j - 1.0                                      
                YI  = i - 1.0                                       

                if (KPP >= 1) then
                    DPPX = 1.0 / KPP                                       
                    DPPY = 1.0 / KPP                                       
                
                    do LP = 1, KPP

                        call AllocateNewParticle (NewParticle, CurrentOrigin%nProperties, &
                                                  CurrentOrigin%NextParticID)

                        call RANDOM_NUMBER(auxX)
                        call RANDOM_NUMBER(auxY)

                        NewParticle%Position%CellJ = XI + DPPX * KPP * auxX
                        NewParticle%Position%CellI = YI + DPPY * KPP * auxY

                        call Convert_CellIJ_XY(NewParticle%Position)
                        call Convert_CellIJ_IJ(NewParticle%Position)

                        !Z position
                        NewParticle%Position%Z = DepthPartic

                        call Convert_Z_CellK  (CurrentOrigin, NewParticle%Position)
                        call Convert_CellK_K  (               NewParticle%Position)

                        NewParticle%Geometry%VolVar = 0.0

                        
                        if (CurrentOrigin%ParticleBoxVolume > 0.0) then
                            NewParticle%Geometry%Volume     = CurrentOrigin%ParticleBoxVolume
                        else
                            select case (Me%Grid%GridType)

                            case (Grid1D)
                                NewParticle%Geometry%Volume =                                           &
                                    (Me%Grid%ParticX(j+1) - Me%Grid%ParticX(j)) * &
                                    (Me%Grid%ParticY(i+1) - Me%Grid%ParticY(i)) * &
                                     BoxThickness
                            
                            case (Grid2D)
                                NewParticle%Geometry%Volume =                                                   &
                                    (Me%Grid%ParticXX(i, j+1) - Me%Grid%ParticXX(i, j)) * &
                                    (Me%Grid%ParticYY(i+1, j) - Me%Grid%ParticYY(i, j)) * &
                                     BoxThickness

                            end select

                        endif

                        !Stores initial Volume
                        NewParticle%Geometry%InitialVolume = NewParticle%Geometry%Volume
                        
                        if (CurrentOrigin%State%Oil)                                    &
                            NewParticle%Geometry%VolumeOil     = NewParticle%Geometry%Volume
                       
                        !Inserts Particle to list
                        call InsertParticleToList (CurrentOrigin, NewParticle, .true.)

                    enddo
                endif    

            endif BOX2D

            endif OP1

            enddo
            enddo

        else FLOAT
    
            FACTOR = 0.0
            do k = KLB, KUB
            do j = JLB, JUB, CurrentOrigin%INCRP 
            do i = ILB, IUB, CurrentOrigin%INCRP 

OP:         if ((Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) .and. &
                (Boxes(i, j, k) == CurrentOrigin%BoxNumber)) then
 
                if (CurrentOrigin%ParticleBoxVolume > 0.0) then                                 

                    VOLCEL = Me%ExternalVar%VolumeZ(I,J,K) + FACTOR    
                                                       
                    KPP    = INT(VOLCEL / CurrentOrigin%ParticleBoxVolume)       
                                                
                    FACTOR = VOLCEL - CurrentOrigin%ParticleBoxVolume * KPP  

                else
                                              
                    KPP = 1                         

                end if
                                                

                ! ---> Considera-se que a particula no ponto X=1,Y=1, Z=1, se encontra
                !      no centro da celula I=1, J=1 , K=1                             
                XI  = j - 1.0                                      
                YI  = i - 1.0                                       
                ZI  = k - 1.0                                       


                ! ---> Distribui aleatoriamente segundo uma distribuicao uniforme   
                !      os tracadores em cada celula X=XI+RANDOM(0:1) (ver NEWPART)  
               if (KPP .GE. 1) then
                    DPPX = 1.0 / KPP                                       
                    DPPY = 1.0 / KPP                                       
                    DPPZ = 1.0  
                
                    do LP = 1, KPP

                        call AllocateNewParticle (NewParticle, CurrentOrigin%nProperties, &
                                                  CurrentOrigin%NextParticID)

                        call RANDOM_NUMBER(auxX)
                        call RANDOM_NUMBER(auxY)

                        NewParticle%Position%CellJ = XI + DPPX * KPP * auxX
                        NewParticle%Position%CellI = YI + DPPY * KPP * auxY

                        call Convert_CellIJ_XY(NewParticle%Position)
                        call Convert_CellIJ_IJ(NewParticle%Position)

                        !Random Z position
                        call RANDOM_NUMBER(auxZ)
                        NewParticle%Position%CellK = ZI + DPPZ * auxZ 
                        call Convert_CellK_Z  (NewParticle%Position)
                        call Convert_CellK_K  (NewParticle%Position)

                        NewParticle%Geometry%VolVar     = 0.0
                        
                        if (CurrentOrigin%ParticleBoxVolume > 0.0) then
                            NewParticle%Geometry%Volume        = CurrentOrigin%ParticleBoxVolume
                        else
                            NewParticle%Geometry%Volume        =                             &
                                   Me%ExternalVar%VolumeZ(NewParticle%Position%I, &
                                                                     NewParticle%Position%J, &
                                                                     NewParticle%Position%K)
                        endif

                        !Stores initial Volume
                        NewParticle%Geometry%InitialVolume = NewParticle%Geometry%Volume

                        if (CurrentOrigin%State%Oil)                                    &
                            NewParticle%Geometry%VolumeOil     = NewParticle%Geometry%Volume

                        !Inserts Particle to list
                        call InsertParticleToList (CurrentOrigin, NewParticle, .true.)
                                                              
                    enddo
                endif
            endif OP
            enddo
            enddo
            enddo

        endif FLOAT


        !Unget The Boxes
        call UngetBoxDif(Me%ObjBoxDif, Boxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'EmissionBox - ModuleLagrangian - ERR01'


    end subroutine EmissionBox

    !--------------------------------------------------------------------------

    subroutine EmissionPoint (CurrentOrigin)
        
        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Local-----------------------------------------------------------------
        integer                                     :: nP
        type (T_Partic), pointer                    :: NewParticle

        !Actualizes Variable Floe
        if (CurrentOrigin%FlowVariable .and.                                             &
            CurrentOrigin%EmissionTemporal == Continuous_)  then
            call ActualizeOrigin (CurrentOrigin)
        endif

        do nP = 1, CurrentOrigin%NbrParticlesIteration

            !Allocates a new Particle
            call AllocateNewParticle (NewParticle, CurrentOrigin%nProperties,            &
                                      CurrentOrigin%NextParticID)
            
            !Sets Initial Position
            NewParticle%Position = CurrentOrigin%Position


            if (CurrentOrigin%Position%MaintainRelative .and. CurrentOrigin%Position%DepthDefinition == Cells) then

                NewParticle%Position%CellK = NewParticle%Position%Depth

                call Convert_CellK_Z (NewParticle%Position)
                call Convert_CellK_K (NewParticle%Position)


            endif
            

            NewParticle%Geometry%VolVar     = 0.0

            !Volume
            select case (CurrentOrigin%EmissionTemporal)

            case (Continuous_)

                !Volume of the Seed Particle
                NewParticle%Geometry%Volume = CurrentOrigin%Flow     *                   &
                                              CurrentOrigin%DT_EMIT  /                   &
                                              CurrentOrigin%NbrParticlesIteration

            case (Instantaneous_)

                !Volume of the Seed Particle
                NewParticle%Geometry%Volume = CurrentOrigin%PointVolume /                &
                                              CurrentOrigin%NbrParticlesIteration

            end select            

            !Stores initial Volume
            NewParticle%Geometry%InitialVolume = NewParticle%Geometry%Volume

            if (CurrentOrigin%State%Oil)                                                 &
               NewParticle%Geometry%VolumeOil     = NewParticle%Geometry%Volume

            !Inserts Particle to list
            call InsertParticleToList (CurrentOrigin, NewParticle, .true.)


            !If simulating a plume give initial velocity and density to the tracer
            if (CurrentOrigin%State%ComputePlume) then
                
                if (Me%Now >= CurrentOrigin%Movement%NextJetActualization) then

                    call ActualizeJetProperties(CurrentOrigin)

                    CurrentOrigin%Movement%NextJetActualization = CurrentOrigin%Movement%NextJetActualization + &
                                                                  CurrentOrigin%Movement%JetDT
                endif

                call GiveJetPropertiesToParticle (CurrentOrigin, NewParticle)

            endif

            if (CurrentOrigin%State%FarFieldBuoyancy) then

                NewParticle%W = CurrentOrigin%Movement%InitialVelocityW                

            endif
            
        enddo


    end subroutine EmissionPoint
    
    !--------------------------------------------------------------------------

    subroutine ActualizeOrigin (CurrentOrigin)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Gets Flow values that limited the current instant
        call GetTimeSerieValue(CurrentOrigin%TimeSerieInput,                     &
                               Me%ExternalVar%Now,                               &
                               CurrentOrigin%FlowColumn,                         &
                               Time1, Value1, Time2, Value2, TimeCycle,          &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeOrigin - ModuleLagrangian - ERR01.'

        if (TimeCycle) then
            CurrentOrigin%Flow = Value1
        else
            !Interpolates Value for current instant
            call InterpolateValueInTime(Me%ExternalVar%Now, Time1,               &
                                        Value1, Time2, Value2, CurrentOrigin%Flow)
        endif


        CurrentProperty => CurrentOrigin%FirstProperty
        do while (associated(CurrentProperty))

            if (CurrentProperty%ConcVariable) then
                call ActualizeConcentration (CurrentOrigin, CurrentProperty)
            endif

            CurrentProperty => CurrentProperty%Next
        enddo

        nullify(CurrentProperty)


    end subroutine ActualizeOrigin 

    !--------------------------------------------------------------------------
    
    subroutine ActualizeConcentration (CurrentOrigin, CurrentProperty)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Property), pointer                  :: CurrentProperty

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Gets Flow values that limited the current instant
        call GetTimeSerieValue(CurrentOrigin%TimeSerieInput,                     &
                               Me%ExternalVar%Now,                               &
                               CurrentProperty%ConcColumn,                       &
                               Time1, Value1, Time2, Value2, TimeCycle,          &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeOrigin - ModuleLagrangian - ERR01.'

        if (TimeCycle) then
            CurrentProperty%Concentration = Value1
        else
            !Interpolates Value for current instant
            call InterpolateValueInTime(Me%ExternalVar%Now, Time1,               &
                                        Value1, Time2, Value2,                   &
                                        CurrentProperty%Concentration)
        endif

    end subroutine ActualizeConcentration        

    !--------------------------------------------------------------------------

    subroutine ActualizeJetProperties(CurrentOrigin)
        
        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:), pointer           :: Temperature3D, Salinity3D
        real                                        :: PlumeDilution, PlumeX, PlumeY, PlumeZ,        &
                                                       PlumeDensity, PlumeTemperature, PlumeSalinity,&
                                                       PlumeU, PlumeV, PlumeW, PlumeMixingHorLength, &
                                                       PlumeThickness 
        real                                        :: Year, Month, Day, Hour, Minutes, Seconds
        integer                                     :: STAT_CALL, I, J, OutPutNumber
        logical                                     :: OutPutJet
        !Begin-----------------------------------------------------------------

        I = CurrentOrigin%Position%I
        J = CurrentOrigin%Position%J

        
        !Gets the temperature and the Salinity from the Eulerian model
        call GetConcentration(Me%ObjWaterProperties,                                    &
                              ConcentrationX    = Temperature3D,                        &
                              PropertyXIDNumber = Temperature_,                         &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR01'

        call GetConcentration(Me%ObjWaterProperties,                                    &
                              ConcentrationX    = Salinity3D,                           &
                              PropertyXIDNumber = Salinity_,                            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR02'

        OutPutNumber = Me%OutPut%NextOutPut
        if (Me%ExternalVar%Now >= Me%OutPut%OutTime(OutPutNumber)) then 
            OutPutJet = .true.
        else
            OutPutJet = .false.
        endif
            
        call ModifyJet(JetID          = CurrentOrigin%Movement%ObjJet,                  &
                       Salinity       = Salinity3D,                                     &
                       Temperature    = Temperature3D,                                  &
                       VelU           = Me%ExternalVar%Velocity_U,                      &
                       VelV           = Me%ExternalVar%Velocity_V,                      &
                       VelW           = Me%ExternalVar%Velocity_W,                      &
                       SZZ            = Me%ExternalVar%SZZ,                             &
                       I              = I,                                              &
                       J              = J,                                              &
                       BottomLayer    = Me%ExternalVar%kFloor(I, J),                    &
                       SurfaceLayer   = Me%ExternalVar%WorkSize%KUB,                    &
                       OutPutOK       = .true.,                                         &
                       STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR03'

        !UnGets the temperature and the Salinity from the Eulerian model
        call UnGetWaterProperties(Me%ObjWaterProperties,                                &
                                  Temperature3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR04'

        call UnGetWaterProperties(Me%ObjWaterProperties,                                &
                                  Salinity3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR05'


        call ExtractDate(Me%ExternalVar%Now, Year, Month, Day, Hour, Minutes, Seconds)


        call GetPlumeDilution(CurrentOrigin%Movement%ObjJet, PlumeDilution, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR07'

        call GetPlumeLocation(CurrentOrigin%Movement%ObjJet, PlumeX, PlumeY, PlumeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR08'

        call GetPlumeDensity(CurrentOrigin%Movement%ObjJet, PlumeDensity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR09'

        call GetPlumeTemperature(CurrentOrigin%Movement%ObjJet, PlumeTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR10'

        call GetPlumeSalinity(CurrentOrigin%Movement%ObjJet, PlumeSalinity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR11'

        call GetPlumeVelocity(CurrentOrigin%Movement%ObjJet, PlumeU, PlumeV, PlumeW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR12'

        call GetPlumeMixingHorLength(CurrentOrigin%Movement%ObjJet, PlumeMixingHorLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR13'

        call GetPlumeThickness(CurrentOrigin%Movement%ObjJet, PlumeThickness, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangian - ERR14'


        write(CurrentOrigin%Movement%JetUnit,'(7F8.0, 2F14.3, F12.3, F16.6, 2F8.4,3F12.5, F10.2, F8.3)') &
                                                Year, Month, Day, Hour, Minutes, Seconds, PlumeDilution, &
                                                PlumeX, PlumeY, PlumeZ, PlumeDensity, PlumeTemperature,  &
                                                PlumeSalinity, PlumeU, PlumeV, PlumeW,                   &
                                                PlumeMixingHorLength, PlumeThickness

    end subroutine ActualizeJetProperties
    !--------------------------------------------------------------------------


    subroutine GiveJetPropertiesToParticle(CurrentOrigin, NewParticle)
        
        !Arguments-------------------------------------------------------------
        type (T_Origin    ), pointer                :: CurrentOrigin
        type (T_Partic    ), pointer                :: NewParticle
        type (T_Property  ), pointer                :: CurrentProperty

        !Local-----------------------------------------------------------------
        real                                        :: AmbientConcentration
        integer                                     :: iP
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------


        call GetPlumeVelocity(JetID     = CurrentOrigin%Movement%ObjJet,                 &
                              PlumeVelU = CurrentOrigin%Movement%InitialVelocityU,       &
                              PlumeVelV = CurrentOrigin%Movement%InitialVelocityV,       &
                              PlumeVelW = CurrentOrigin%Movement%InitialVelocityW,       &
                              STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangian - ERR01'

        NewParticle%U = CurrentOrigin%Movement%InitialVelocityU
        NewParticle%V = CurrentOrigin%Movement%InitialVelocityV
        NewParticle%W = CurrentOrigin%Movement%InitialVelocityW

        call GetPlumeDilution(JetID         = CurrentOrigin%Movement%ObjJet,             &
                              PlumeDilution = CurrentOrigin%Movement%PlumeDilution,      &
                              STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangian - ERR02'

        !Initial volume Dilution 
        NewParticle%Geometry%Volume        = CurrentOrigin%Movement%PlumeDilution * NewParticle%Geometry%Volume  

        NewParticle%Geometry%VolVar = 0.0

        iP = 1

        CurrentProperty => CurrentOrigin%FirstProperty
        do while (associated(CurrentProperty))

            if      (CurrentProperty%ID == Salinity_) then

                call GetPlumeSalinity (JetID         = CurrentOrigin%Movement%ObjJet,    &
                                       PlumeSalinity = NewParticle%Concentration(iP),    &
                                       STAT          = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangian - ERR03'
               


            else if (CurrentProperty%ID == Temperature_) then
                
                call GetPlumeTemperature (JetID            = CurrentOrigin%Movement%ObjJet,&
                                          PlumeTemperature = NewParticle%Concentration(iP),&
                                          STAT             = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangian - ERR04'

            else

                !Properties dilution
                call GetAmbientConcentration (CurrentProperty,                           &
                                              NewParticle%Position,                      &
                                              AmbientConcentration)

            
                NewParticle%Concentration(iP) =(CurrentProperty%Concentration +                                     &
                                                AmbientConcentration * (CurrentOrigin%Movement%PlumeDilution - 1)) / &
                                                CurrentOrigin%Movement%PlumeDilution

            endif

            NewParticle%Mass         (iP) = NewParticle%Concentration(iP) *              &
                                            NewParticle%Geometry%Volume
  
            iP = iP + 1
            CurrentProperty => CurrentProperty%Next

        enddo

        !Location after initial dilution
        call GetPlumeLocation(JetID  = CurrentOrigin%Movement%ObjJet,                    &
                              PlumeX = NewParticle%Position%X,                           &
                              PlumeY = NewParticle%Position%Y,                           &
                              PlumeZ = NewParticle%Position%Z,                           &
                              STAT   = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangian - ERR05'

        call Convert_Z_CellK   (CurrentOrigin, NewParticle%Position)
        call Convert_CellK_K   (NewParticle%Position)

        call Convert_XY_CellIJ (NewParticle%Position)
        call Convert_CellIJ_IJ (NewParticle%Position)

        if (Me%ExternalVar%OpenPoints3D(NewParticle%Position%i,               &
                                                   NewParticle%Position%j,               &
                                                   NewParticle%Position%k) /= OpenPoint) then
            
            !Sets Horizontal Position equal to origin
            NewParticle%Position%X = CurrentOrigin%Position%X
            NewParticle%Position%Y = CurrentOrigin%Position%Y

            call Convert_XY_CellIJ (NewParticle%Position)
            call Convert_CellIJ_IJ (NewParticle%Position)

            write(*,*) 'The jet of Origin = ', trim(CurrentOrigin%Name), ' intersected land' 

        endif

    end subroutine GiveJetPropertiesToParticle

    subroutine EmissionAccident (CurrentOrigin)
        
        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Local-----------------------------------------------------------------
        type (T_Partic), pointer                    :: NewParticle
        real                                        :: AreaTotal, AreaParticle, AreaSum
        real                                        :: AccRad
        integer                                     :: II, JJ
        real                                        :: aux, ang
        real                                        :: Dist, XDif, YDif
        real                                        :: API
        integer                                     :: STAT_CALL
        integer                                     :: i, j, k, KUB
        logical                                     :: Emited
        type (T_Position)                           :: NewPosition


        !Cell of the Accident
        JJ = CurrentOrigin%position%J
        II = CurrentOrigin%position%I         

        
        select case (CurrentOrigin%AccidentMethod)

        case (Fay_)

            !Total Area given by the F_Fay Method
            if (CurrentOrigin%ObjOil == 0) then
                write(*,*)'Origin : ', trim(CurrentOrigin%Name)
                write(*,*)'The Fay Method can only be used with Oil Accident'
                write(*,*)'Use ACCIDENT_METHOD : Thickness'
                write(*,*)'or intialize the oil property'
            endif

            call GetOilAPI (CurrentOrigin%ObjOil, API = API, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'EmissionAccident - ModuleLagrangian - ERR01'
            
            AreaTotal = F_FayArea(VolInic = CurrentOrigin%PointVolume, API = API)

        case (Thickness_)

            AreaTotal = CurrentOrigin%PointVolume / CurrentOrigin%Movement%ThicknessMeters

        end select


        !Calcula a area ocupada por cada tracador
        AreaParticle = AreaTotal / float(CurrentOrigin%NbrParticlesIteration)

        !Allocates the first Particle
        call AllocateNewParticle (NewParticle, CurrentOrigin%nProperties,                &
                                  CurrentOrigin%NextParticID)

        !Stores first particle / from the center of the cell
        NewParticle%Position                = CurrentOrigin%Position
        NewParticle%Geometry%Volume         = CurrentOrigin%PointVolume /                &
                                              float(CurrentOrigin%NbrParticlesIteration)

        NewParticle%Geometry%VolVar         = 0.0

        !Stores initial Volume
        NewParticle%Geometry%InitialVolume = NewParticle%Geometry%Volume
     
        if (CurrentOrigin%State%Oil)                                    &
            NewParticle%Geometry%VolumeOil     = NewParticle%Geometry%Volume

        !Inserts Particle to list
        call InsertParticleToList (CurrentOrigin, NewParticle, .true.)


        !Radius of the Accident
        AccRad = sqrt(AreaTotal / PI)

        !Area Ocupied until know
        AreaSum = AreaParticle
        do while (abs(AreaTotal - (AreaParticle/2.)) > AreaSum)

            !Allocates new particle
            call AllocateNewParticle (NewParticle, CurrentOrigin%nProperties,            &
                                      CurrentOrigin%NextParticID)
                
            !While it doesnt hit a waterpoint, try emitting
            Emited =.FALSE.
            do while(.NOT. Emited)

                !Aleatory Angle [0, 2PI]
                call RANDOM_NUMBER(aux)
                Ang = aux * PI * 2.0

                !Aleatory Distance [0, AccRad]
                call RANDOM_NUMBER(aux)
                Dist = aux * AccRad

                !Distance from origin in meters
                Xdif = Dist * cos(Ang)
                Ydif = Dist * sin(Ang)

                !New Position
                NewPosition%X = CurrentOrigin%Position%X + Xdif
                NewPosition%Y = CurrentOrigin%Position%Y + Ydif
                NewPosition%Z = CurrentOrigin%Position%Z

                !Converts Horizontal Position
                call Convert_XY_CellIJ (NewPosition)
                call Convert_CellIJ_IJ (NewPosition)

                i   = NewPosition%I
                j   = NewPosition%J
                KUB = Me%ExternalVar%WorkSize%KUB
                
                !Stores the Particle
                if (Me%ExternalVar%OpenPoints3D(i, j, KUB) == OpenPoint) then


                    !If the vertical position of the Particle isnt a waterpoint, put it 
                    !close to the suface
                    call Convert_Z_CellK (CurrentOrigin, NewPosition)
                    call Convert_CellK_K (               NewPosition)
                    k = NewPosition%K
                    if (Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then
                        NewParticle%Position = NewPosition
                    else
                        NewPosition%Z = Me%ExternalVar%SZZ(i, j, KUB)
                        call Convert_Z_CellK (CurrentOrigin, NewPosition)
                        call Convert_CellK_K (               NewPosition)
                        NewParticle%Position = NewPosition
                    endif

                    NewParticle%Geometry%Volume = CurrentOrigin%PointVolume /            &
                                                  float(CurrentOrigin%NbrParticlesIteration)

                    NewParticle%Geometry%VolVar = 0.0

                    !Stores initial Volume
                    NewParticle%Geometry%InitialVolume = NewParticle%Geometry%Volume

                    if (CurrentOrigin%State%Oil)                                    &
                        NewParticle%Geometry%VolumeOil     = NewParticle%Geometry%Volume

                    !Inserts Particle to List
                    call InsertParticleToList (CurrentOrigin, NewParticle, .true.)

                    AreaSum = AreaSum + AreaParticle
                    Emited  =.true.

                else

                    Emited = .false.

                endif

            enddo

        enddo


        CurrentOrigin%AccidentFinished = .true.


    end subroutine EmissionAccident

    !--------------------------------------------------------------------------

    subroutine ConstructEmissionTime (CurrentOrigin)
        
        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Local-----------------------------------------------------------------

        !Next Emission time
        if (CurrentOrigin%Old) then
            CurrentOrigin%NextEmission = Me%ExternalVar%Now +                 &
                                         CurrentOrigin%DT_Emit
        else
            CurrentOrigin%NextEmission = Me%ExternalVar%Now
        endif

        if (CurrentOrigin%NextEmission < CurrentOrigin%StartEmission) then
            CurrentOrigin%NextEmission = CurrentOrigin%StartEmission
        endif

    end subroutine ConstructEmissionTime
    
    !--------------------------------------------------------------------------

    subroutine ConstructParticOil (NewOrigin, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Origin),     pointer                :: NewOrigin
        integer                                     :: ClientNumber

        !Local-----------------------------------------------------------------
        logical                                     :: OilSectionFound
        integer                                     :: STAT_CALL


        !Extracts the Oil Section from the Origin
        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,             &
                                   oil_begin, oil_end,                                   &
                                   OilSectionFound, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticOil - ModuleLagrangian - ERR00'


        if (OilSectionFound .and. .not. NewOrigin%Old) then


            !Starts Oil
            call StartOil(OilID             = NewOrigin%ObjOil,                      &     
                          TimeID            = Me%ObjTime,                            &     
                          EnterDataID       = Me%ObjEnterData,                       &
                          HorizontalGridID  = Me%ObjHorizontalGrid,                  &     
                          GeometryID        = Me%ObjGeometry,                        &     
                          MapID             = Me%ObjMap,                             &     
                          DT                = Me%DT_PARTIC,                          &     
                          ContCalc          = NewOrigin%Old,                         &
                          ExtractType       = FromBlockInBlock,                      &
                          STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticOil - ModuleLagrangian - ERR01'

        elseif (.not. OilSectionFound) then

            write (*,*)'Oil Section not defined for Origin : ', trim(NewOrigin%Name)
            stop       'ConstructParticOil - ModuleLagrangian - ERR03'

        endif
  
            !Allocates GridThickness
            allocate (NewOrigin%GridThickness(Me%ExternalVar%Size%ILB:            &
                                              Me%ExternalVar%Size%IUB,            &
                                              Me%ExternalVar%Size%JLB:            &
                                              Me%ExternalVar%Size%JUB),           &
                                              STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticOil - ModuleLagrangian - ERR04'


            !Allocates OilGridConcentration
            allocate (NewOrigin%OilGridConcentration(Me%ExternalVar%Size%ILB:            &
                                                     Me%ExternalVar%Size%IUB,            &
                                                     Me%ExternalVar%Size%JLB:            &
                                                     Me%ExternalVar%Size%JUB),           &
                                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticOil - ModuleLagrangian - ERR05'

    end subroutine ConstructParticOil

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine VerifyBeachingProbabilities

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,    pointer, dimension(:)              :: BoxesBeachingProbability
        integer, pointer, dimension(:,:,:)          :: BeachingProbBoxes
        integer                                     :: flag
        integer                                     :: NumberOfBoxes
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k, Box
        integer                                     :: STAT_CALL

        !-----------------------------------------------------------------------
        
        allocate (Me%BeachingProbability(Me%ExternalVar%Size%ILB:  &
                                                    Me%ExternalVar%Size%IUB,  &
                                                    Me%ExternalVar%Size%JLB:  &
                                                    Me%ExternalVar%Size%JUB,  &
                                                    Me%ExternalVar%Size%KLB:  &
                                                    Me%ExternalVar%Size%KUB), &
                                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangian - ERR01'
        
        if ( .NOT. Me%State%HaveBeachingProbBox) then

            Me%BeachingProbability = Me%DefaultBeachingProbability

        else if (Me%State%HaveBeachingProbBox) then



            !Shorten
            ILB    = Me%ExternalVar%WorkSize%ILB
            IUB    = Me%ExternalVar%WorkSize%IUB
            JLB    = Me%ExternalVar%WorkSize%JLB
            JUB    = Me%ExternalVar%WorkSize%JUB
            KLB    = Me%ExternalVar%WorkSize%KLB
            KUB    = Me%ExternalVar%WorkSize%KUB
        
            !Get the boxes
            call GetBoxes(Me%ObjBeachingProbBox, BeachingProbBoxes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangian - ERR02'

            call GetNumberOfBoxes(Me%ObjBeachingProbBox, NumberOfBoxes3D = NumberOfBoxes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangian - ERR03'

            allocate (BoxesBeachingProbability(NumberOfBoxes))

            call GetData(BoxesBeachingProbability,                                       &
                         Me%ObjEnterData,                                     &
                         flag,                                                           &
                         SearchType = FromFile,                                          &
                         keyword    = 'BOXES_BEACHING_PROB',                             &
                         ClientModule = 'ModuleLagrangian',                              &
                         STAT       = STAT_CALL)
            if       (STAT_CALL .EQ. SIZE_ERR_)  then
                write(*,*) 
                write(*,*) 'Error calling GetData.  '
                write(*,*) 'Number of box values is incorrect:'
                write(*,*) '    NumberOfBoxes =', NumberOfBoxes
                write(*,*) '    BoxesData   =', flag
                stop       'Subroutine VerifyBeachingProbabilities; Module ModuleLagrangian. ERR04.'

            else if ((STAT_CALL .NE. SIZE_ERR_) .AND.                      &
                     (STAT_CALL .NE. SUCCESS_)) then                                                                        
                stop 'Subroutine VerifyBeachingProbabilities; Module ModuleLagrangian. ERR05.'
            end if               

            if (flag==0) then
                write(*,*) 
                write(*,*) 'Error do not have the box beaching probability.'           
                stop       'Subroutine VerifyBeachingProbabilities; Module ModuleLagrangian. ERR06.'
            end if



            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                Box = BeachingProbBoxes(i, j, k)
            
                if (Box .GT. 0 ) then
            
                    Me%BeachingProbability (i,j,k) = BoxesBeachingProbability(Box)
                     
                else if (Box .LE. 0 ) then
            
                    Me%BeachingProbability (i,j,k) = Me%DefaultBeachingProbability
            
                endif
            enddo
            enddo
            enddo

            deallocate (BoxesBeachingProbability)

            !Unget The Boxes
            call UngetBoxDif(Me%ObjBeachingProbBox, BeachingProbBoxes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangian - ERR07'

        end if

    end subroutine VerifyBeachingProbabilities

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Output

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, pointer, dimension(:, :)              :: Bathymetry
        integer, pointer, dimension(:, :, :)        :: WaterPoints3D
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Bounds
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB

        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB

        KLB = Me%ExternalVar%WorkSize%KLB
        KUB = Me%ExternalVar%WorkSize%KUB

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5,                                  &
                                 trim(Me%Files%TransientHDF)//"5",            &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR05'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR06'

        !Gets a pointer to Bathymetry
        call GetGridData        (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR01'

        !Gets WaterPoints3D
        call GetWaterPoints3D   (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR03'


        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB,       &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR17'


        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",         &
                              Array2D = Bathymetry,                                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR18'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",      &
                              Array3D = WaterPoints3D,                                   &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR19'

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR23'

        !Ungets the Bathymetry
        call UngetGridData (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR24'

        !Ungets the WaterPoints
        call UnGetMap        (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangian - ERR25'

    end subroutine ConstructHDF5Output
  
    !--------------------------------------------------------------------------

    subroutine ConstructParticStatistic

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:, :, :, :), pointer        :: GridConc
        real, dimension(:, :, :, :), pointer        :: GridMaxTracer
        type (T_Property), pointer                  :: CurrentProperty
        real, dimension(:, :, :), pointer           :: GridConc3D 
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: iProp, STAT_CALL


        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        KLB    = Me%ExternalVar%Size%KLB
        KUB    = Me%ExternalVar%Size%KUB

        !This test is done for simply reason
        if (Me%nGroups > 1) then
            write(*,*)'Cannot calculate Statistics for simulation with more then one Group'
            stop 'ConstructParticStatistic - ModuleLagrangian - ERR01'
        endif


        !Allocates auxiliar variable
        allocate (GridConc3D (ILB:IUB, JLB:JUB, KLB:KUB         ))


        !Fills Grid concentration
        call FillGridConcentration  (1, GridConc, GridMaxTracer = GridMaxTracer)

        iProp = 1
        CurrentProperty => Me%FirstOrigin%FirstProperty
        do while (associated(CurrentProperty))

            if (CurrentProperty%Statistics) then

                GridConc3D = GridConc(:, :, :, iProp)

                call ConstructStatistic (CurrentProperty%Statistic1_ID,                     &
                                         ObjTime          = Me%ObjTime,                     &
                                         ObjHDF5          = Me%ObjHDF5,                     &
                                         Size             = Me%ExternalVar%Size,            &
                                         WorkSize         = Me%ExternalVar%WorkSize,        &
                                         DataFile         = CurrentProperty%StatisticsFile, &
                                         Name             = CurrentProperty%Name,           &
                                         STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'ConstructParticStatistic - ModuleLagrangian - ERR02'

                if (Me%OutPut%ConcMaxTracer) then

                    GridConc3D = GridMaxTracer(:, :, :, iProp)

                    call ConstructStatistic (CurrentProperty%Statistic2_ID,                 &
                                             ObjTime          = Me%ObjTime,                 &
                                             ObjHDF5          = Me%ObjHDF5,                 &
                                             Size             = Me%ExternalVar%Size,        &
                                             WorkSize         = Me%ExternalVar%WorkSize,    &
                                             DataFile         = CurrentProperty%StatisticsFile,&
                                             Name             = CurrentProperty%Name,       &
                                             GroupName        = "MaxTracer",                &
                                             STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                              &
                        stop 'ConstructParticStatistic - ModuleLagrangian - ERR03'

                endif

                if (CurrentProperty%StatisticsLag) then
                
                    call GetStatisticClassesNumber(CurrentProperty%Statistic1_ID,        &
                                                   CurrentProperty%nClassesLag, STAT= STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                           &
                        stop 'ConstructParticStatistic - ModuleLagrangian - ERR04'

                    allocate(CurrentProperty%FrequencyLag(ILB:IUB,JLB:JUB,KLB:KUB,1:CurrentProperty%nClassesLag))

                    CurrentProperty%FrequencyLag(:,:,:,:) = 0.

                endif


            endif
            CurrentProperty => CurrentProperty%Next
            iProp = iProp + 1
        enddo

        !Deallocates auxiliar variable
        deallocate (GridConc3D   )
        deallocate (GridConc     )
        if (Me%OutPut%ConcMaxTracer) deallocate (GridMaxTracer)

    end subroutine ConstructParticStatistic

    !--------------------------------------------------------------------------

    subroutine ConstructParticLightExtinction
        
        !Local-----------------------------------------------------------------
        integer                                     :: iGroup, STAT_CALL  
        integer                                     :: LightExtinctionID
          
        !Begin-----------------------------------------------------------------

        allocate (Me%Light%ObjLightExtinction(Me%nGroups))

        Me%Light%ObjLightExtinction(:) = 0

        allocate (Me%Light%TopRadiationCells(Me%ExternalVar%Size%ILB:Me%ExternalVar%Size%IUB, &
                                             Me%ExternalVar%Size%JLB:Me%ExternalVar%Size%JUB, &
                                             Me%ExternalVar%Size%KLB:Me%ExternalVar%Size%KUB, &
                                             1 : Me%nGroups))

        call GetWaterPropertiesSubModulesID(WaterPropertiesID      = Me%ObjWaterProperties, &
                                            LightExtinctionID      = LightExtinctionID,     &
                                            STAT                   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticLightExtinction - Lagrangian - ERR00'


        if(LightExtinctionID == 0)then

            Me%Light%Compute = ON
            
            do iGroup = 1, Me%nGroups

                call ConstructLightExtinction (Me%Light%ObjLightExtinction(iGroup),          &
                                               Me%ObjTime, Me%ObjEnterData,                  &
                                               Me%ExternalVar%Size, Me%ExternalVar%WorkSize, &
                                               STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructParticLightExtinction - Lagrangian - ERR01'
            enddo

        else

            Me%Light%Compute = OFF

            do iGroup = 1, Me%nGroups
                Me%Light%ObjLightExtinction(iGroup) = AssociateInstance(mLIGHTEXTINCTION_,LightExtinctionID)
            enddo

        end if




    end subroutine ConstructParticLightExtinction

    !--------------------------------------------------------------------------

    subroutine ConstructLog

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: iGroup


        write(*, *)"----------------------- LAGRANGIAN -----------------------"
        write(*, *)
        write(*, *)"Number of Origins : ", Me%nOrigins
        write(*, *)"Number of Groups  : ", Me%nGroups
        write(*, *)

        do iGroup = 1, Me%nGroups
            write(*, *)"GroupID           : ", Me%GroupIDs(iGroup)
            write(*, *)"---Number of Orig.: ", Me%nOriginsGroup(iGroup)
            write(*, *)
        enddo

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            write(*, *)"Origin            : ", trim(CurrentOrigin%Name)
            write(*, *)"---Group ID       : ", CurrentOrigin%GroupID
            write(*, *)"---Number of Part.: ", CurrentOrigin%nParticle
            write(*, *)"---Number of Prop.: ", CurrentOrigin%nProperties
            write(*, *)

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine ConstructLog

    !--------------------------------------------------------------------------

    subroutine AllocateNewOrigin (NewOrigin)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: NewOrigin

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        nullify  (NewOrigin)
        allocate (NewOrigin, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateNewOrigin - ModuleLagrangian - ERR01'
        nullify  (NewOrigin%Next)
        nullify  (NewOrigin%FirstPartic)
        nullify  (NewOrigin%FirstProperty)
        nullify  (NewOrigin%GridThickness)
        nullify  (NewOrigin%OilGridConcentration)

        NewOrigin%nProperties  = 0
        NewOrigin%nParticle    = 0
        NewOrigin%NextParticID = 1

    end subroutine AllocateNewOrigin

    !--------------------------------------------------------------------------

    subroutine InsertOriginToList(FirstOrigin, NewOrigin, nOrigins)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: FirstOrigin
        type (T_Origin), pointer                    :: NewOrigin
        integer                                     :: nOrigins

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin  => null()
        type (T_Origin), pointer                    :: PreviousOrigin => null()


        !Inserts a new origin to the list if origins
        if (.not. associated(FirstOrigin)) then
            FirstOrigin => NewOrigin
        else
            PreviousOrigin => FirstOrigin
            CurrentOrigin  => PreviousOrigin%Next
            do while (associated(CurrentOrigin))
                PreviousOrigin => CurrentOrigin
                CurrentOrigin  => PreviousOrigin%Next
            enddo
            PreviousOrigin%Next => NewOrigin
        endif

        !Increments origin number by one
        nOrigins           = nOrigins + 1 
        NewOrigin%ID       = nOrigins

    end subroutine InsertOriginToList

    !--------------------------------------------------------------------------

    subroutine DeleteOrigin (FirstOrigin, OriginToDelete, nOrigins)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: FirstOrigin
        type (T_Origin), pointer                    :: OriginToDelete
        integer                                     :: nOrigins

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin    => null()
        type (T_Origin), pointer                    :: NextOrigin       => null()
        type (T_Origin), pointer                    :: PreviousOrigin   => null()
        type (T_Partic), pointer                    :: CurrentParticle  => null()
        type (T_Property), pointer                  :: CurrentProperty  => null()
        integer                                     :: STAT_CALL

        nullify (PreviousOrigin)
        CurrentOrigin => FirstOrigin

        do 
            if (CurrentOrigin%ID == OriginToDelete%ID) then

                NextOrigin => CurrentOrigin%Next
 
                !Updates first origin
                if (OriginToDelete%ID == FirstOrigin%ID) then
                    FirstOrigin => NextOrigin
                endif

                !Kill oil
                if (OriginToDelete%ObjOil /= 0) then
                    call KillOil (OriginToDelete%ObjOil, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangian - ERR01'
                endif

                !Kill WaterQuality
                if (OriginToDelete%State%WQM) then
                    call KillWaterQuality (OriginToDelete%WaterQualityID, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangian - ERR02'
                endif

                !Kill TimeSerie
                if (OriginToDelete%MovingOrigin) then
                    call KillTimeSerie(OriginToDelete%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangian - ERR03'
                endif

                !Kills ObjJet
                if (OriginToDelete%Movement%ObjJet /= 0) then
                   

                    write(OriginToDelete%Movement%JetUnit,*) '<EndTimeSerie>'

                    call UnitsManager(OriginToDelete%Movement%JetUnit, CLOSE_FILE, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangian - ERR03a'

                    call KillJet (JetID = OriginToDelete%Movement%ObjJet,  STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangian - ERR04'
                endif

                
                if (OriginToDelete%FlowVariable) then
                    call KillTimeSerie(OriginToDelete%TimeSerieInput, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangian - ERR04a'
                endif



                !Kill PropertyList
                CurrentProperty => OriginToDelete%FirstProperty
                do while (associated(CurrentProperty))
                    call DeleteProperty (OriginToDelete, CurrentProperty)
                    CurrentProperty => OriginToDelete%FirstProperty
                enddo

                !Kill ParticleList
                CurrentParticle => OriginToDelete%FirstPartic
                do while (associated(CurrentParticle))
                    call DeleteParticle (OriginToDelete, CurrentParticle)
                    CurrentParticle => OriginToDelete%FirstPartic
                enddo
                    
                !Update Previous Next
                if (associated(PreviousOrigin)) then
                    PreviousOrigin%Next => NextOrigin
                endif

                !Deallocate Origin
                nullify       (OriginToDelete%Next)
                deallocate    (OriginToDelete)
                nullify       (OriginToDelete)

                !Decreases number of origins
                nOrigins = nOrigins - 1

                if (nOrigins == 0) then
                    nullify (FirstOrigin) 
                endif

                return
            endif

            PreviousOrigin => CurrentOrigin
            CurrentOrigin  => CurrentOrigin%Next

        enddo
        

    end subroutine DeleteOrigin

    !--------------------------------------------------------------------------

    subroutine AllocateNewProperty(NewProperty)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        nullify  (NewProperty)
        allocate (NewProperty, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateNewProperty - ModuleLagrangian - ERR01'
        nullify  (NewProperty%Next)

    end subroutine AllocateNewProperty

    !--------------------------------------------------------------------------

    subroutine InsertPropertyToList(Origin, NewProperty, SetStates)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: Origin
        type (T_Property), pointer                  :: NewProperty
        logical                                     :: SetStates

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty   => null()
        type (T_Property), pointer                  :: PreviousProperty  => null()

        !Inserts a new property to the list of properties
        if (.not. associated(Origin%FirstProperty)) then
            Origin%FirstProperty => NewProperty
        else
            PreviousProperty => Origin%FirstProperty
            CurrentProperty  => PreviousProperty%Next
            do while (associated(CurrentProperty))
                PreviousProperty => CurrentProperty
                CurrentProperty  => PreviousProperty%Next
            enddo
            PreviousProperty%Next => NewProperty
        endif

        Origin%nProperties      = Origin%nProperties + 1

        !Statistics    
        if (NewProperty%Statistics) then
            Origin%State%Statistics         = ON
            Me%State%Statistics  = ON
        endif

        if (SetStates) then

            select case (NewProperty%ID)

            case (Phytoplankton_                  )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON
                                            
            case (Zooplankton_                    )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON
                                            
            case (Larvae_                         )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON
                Origin%State%Larvae     = ON
                Me%State%Larvae         = ON

            case (Age_                            )
                if(.not.NewProperty%NoWQM) then
                 Origin%State%WQM        = ON
                 Me%State%WQM            = ON
                endif
                  
                                            

            case (PON_)
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (DONRefractory_                  )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON
           
            case (DONNon_Refractory_              )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON
           
            case (Ammonia_                        )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (Nitrate_                        )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (Nitrite_                        )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (BOD_                            )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (Oxygen_                         )
                if(.not.NewProperty%NoWQM) then
                Origin%State%WQM        = ON
                Me%State%WQM            = ON
                endif

            case (Ciliate_                        )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (Bacteria_                       )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (POP_                            )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (DOPRefractory_                  )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (DOPNon_Refractory_              )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (Inorganic_Phosphorus_           )
                Origin%State%WQM        = ON
                Me%State%WQM            = ON

            case (Oil_                            )
                Origin%State%Oil        = ON
                Me%State%Oil            = ON

            case (Fecal_Coliforms_                )
                Origin%State%FCF        = ON
                Me%State%FCF            = ON

            end select

        endif


    end subroutine InsertPropertyToList

    !--------------------------------------------------------------------------

    subroutine DeleteProperty (Origin, PropertyToDelete)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: Origin
        type (T_Property), pointer                  :: PropertyToDelete

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty    => null()
        type (T_Property), pointer                  :: PreviousProperty   => null()
        type (T_Property), pointer                  :: NextProperty       => null()

        nullify (PreviousProperty)
        CurrentProperty => Origin%FirstProperty

        do 
            if (CurrentProperty%ID == PropertyToDelete%ID) then

                NextProperty     => CurrentProperty%Next
 
                !Updates first origin
                if (PropertyToDelete%ID == Origin%FirstProperty%ID) then
                    Origin%FirstProperty => NextProperty
                endif

                !Updates Previous Next
                if (associated(PreviousProperty)) then
                    PreviousProperty%Next => NextProperty
                endif

                !Deallocate properties
                nullify       (PropertyToDelete%Next)
                deallocate    (PropertyToDelete)
                nullify       (PropertyToDelete)

                !Decreases number of properties
                Origin%nProperties = Origin%nProperties - 1

                if (Origin%nProperties == 0) then
                    nullify (Origin%FirstProperty) 
                endif

                return
            endif

            PreviousProperty  => CurrentProperty
            CurrentProperty   => CurrentProperty%Next

        enddo

    end subroutine DeleteProperty

    !--------------------------------------------------------------------------

    subroutine AllocateNewParticle (NewPartic, nProp, NextParticID)

        !Arguments-------------------------------------------------------------
        type (T_Partic), pointer                    :: NewPartic
        integer                                     :: nProp
        integer                                     :: NextParticID

        nullify  (NewPartic)
        allocate (NewPartic)
        nullify  (NewPartic%Next)
        nullify  (NewPartic%Prev)

        if (nProp > 0) then
            allocate(NewPartic%Concentration (nProp))
            allocate(NewPartic%Mass          (nProp))
        else
            nullify (NewPartic%Concentration        )
            nullify (NewPartic%Mass                 )
        endif

        NewPartic%KillPartic = OFF

        NewPartic%TpercursoX = abs(null_real)
        NewPartic%TpercursoY = abs(null_real)
        NewPartic%TpercursoZ = abs(null_real)

        NewPartic%ID         = NextParticID
        NextParticID         = NextParticID + 1


    end subroutine AllocateNewParticle

    !--------------------------------------------------------------------------

    subroutine InsertParticleToList (Origin, NewParticle, InitConcentration)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: Origin
        type (T_Partic), pointer                    :: NewParticle
        logical                                     :: InitConcentration

        !Local-----------------------------------------------------------------
        type (T_Partic), pointer                    :: CurrentParticle  => null()
        type (T_Partic), pointer                    :: PreviousParticle => null()
        type (T_Property), pointer                  :: CurrentProperty  => null()
        integer                                     :: iP

        !Initializes Mass and Concentration of the new particle
        if (InitConcentration) then

            iP = 1
            CurrentProperty => Origin%FirstProperty
            do while (associated(CurrentProperty))
                
                NewParticle%Concentration(iP) = CurrentProperty%Concentration
                NewParticle%Mass         (iP) = CurrentProperty%Concentration *          &
                                                NewParticle%Geometry%Volume
        
                iP = iP + 1
                CurrentProperty => CurrentProperty%Next
            enddo

        endif


        if (Origin%State%Deposition) then
        
            NewParticle%TauErosion = Origin%Deposition%TauErosion
            NewParticle%Deposited  = Origin%Deposition%BottomEmission

        endif


        !Inserts a new property to the list of properties
        if (.not. associated(Origin%FirstPartic)) then
            Origin%FirstPartic => NewParticle
        else
            PreviousParticle => Origin%FirstPartic
            CurrentParticle  => PreviousParticle%Next
            do while (associated(CurrentParticle))
                PreviousParticle => CurrentParticle
                CurrentParticle  => PreviousParticle%Next
            enddo
            PreviousParticle%Next => NewParticle
            NewParticle%Prev      => PreviousParticle
        endif

        Origin%nParticle = Origin%nParticle + 1


    end subroutine InsertParticleToList

    !--------------------------------------------------------------------------

    subroutine DeleteParticle (Origin, ParticleToDelete)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: Origin
        type (T_Partic), pointer                    :: ParticleToDelete

        !Local-----------------------------------------------------------------
        type (T_Partic), pointer                    :: CurrentPartic => null()
        type (T_Partic), pointer                    :: NextParticle  => null()
        type (T_Partic), pointer                    :: PrevParticle  => null()

        logical                                     :: ParticleDeleted
        !Begin-----------------------------------------------------------------

        ParticleDeleted = .false. 

        CurrentPartic => Origin%FirstPartic

d1:     do while (associated(CurrentPartic)) 
                     
i1:         if (CurrentPartic%ID == ParticleToDelete%ID) then


                PrevParticle => CurrentPartic%Prev
                NextParticle => CurrentPartic%Next

 
                !Updates foward pointer
                if (associated(CurrentPartic%Prev)) then
                    if (associated(NextParticle)) then
                        PrevParticle%Next => NextParticle
                    else
                        nullify(PrevParticle%Next)
                    endif
                endif

                !Updates backward pointer
                if (associated(CurrentPartic%Next)) then
                    if (associated(PrevParticle)) then
                        NextParticle%Prev => PrevParticle
                    else
                        nullify(NextParticle%Prev)
                    endif
                endif

                !Updates first origin
                if (ParticleToDelete%ID == Origin%FirstPartic%ID) then
                    Origin%FirstPartic => NextParticle
                endif


                !Delete Concentration and Mass Buffer
                if (associated(ParticleToDelete%Concentration)) then
                    deallocate(ParticleToDelete%Concentration)
                    nullify   (ParticleToDelete%Concentration)
                endif

                if (associated(ParticleToDelete%Mass         )) then
                    deallocate(ParticleToDelete%Mass         )
                    nullify   (ParticleToDelete%Mass         )
                endif

                !Deallocate Particle
!                nullify       (ParticleToDelete%Next)
!                nullify       (ParticleToDelete%Prev)
                deallocate    (ParticleToDelete)
                nullify       (ParticleToDelete, CurrentPartic)

                !Decreases number of particle
                Origin%nParticle = Origin%nParticle - 1

                if (Origin%nParticle == 0) then
                    nullify (Origin%FirstPartic) 
                endif

                ParticleDeleted = .true. 

                exit

            endif i1

            CurrentPartic => CurrentPartic%Next

        enddo d1

        if (.not. ParticleDeleted) then
            stop 'DeleteParticleDeleted - ModuleLagrangian - ERR10'
        endif

    end subroutine DeleteParticle
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !----------------------------------------------------------------------------
    
    subroutine SetLagrangianAtmPressure(LagrangianID, AtmPressure, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: LagrangianID
        real, pointer, dimension(:,:)               :: AtmPressure
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            Me%ExternalVar%AtmPressure => AtmPressure
                        
            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetLagrangianAtmPressure

    !----------------------------------------------------------------------------
    
    subroutine SetLagrangianSolarRadiation(LagrangianID, SolarRadiation, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: LagrangianID
        real, pointer, dimension(:,:)               :: SolarRadiation
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            Me%ExternalVar%SurfaceRadiation => SolarRadiation
                        
            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetLagrangianSolarRadiation

    !----------------------------------------------------------------------------
    
    subroutine SetLagrangianWind(LagrangianID, WindX, WindY, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: LagrangianID
        real, pointer, dimension(:,:)               :: WindX, WindY
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            if(Me%State%Wind)then
                Me%ExternalVar%WindX => WindX
                Me%ExternalVar%WindY => WindY
            end if


            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetLagrangianWind


    !----------------------------------------------------------------------

    subroutine SetLagrangianShear(LagrangianID, ShearStress, ShearVelocity, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: LagrangianID
        real, dimension(:,:), pointer   :: ShearStress, ShearVelocity
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_              
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            if(Me%State%Deposition) Me%ExternalVar%BottomStress  => ShearStress            
            if(Me%State%ShearVel  ) Me%ExternalVar%ShearVelocity => ShearVelocity 
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetLagrangianShear


    !----------------------------------------------------------------------

    subroutine GetLagrangianAirOptions(LagrangianID, Oil, Wind, WaterQuality, T90Variable, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        logical, optional, intent(OUT)              :: Oil, Wind, WaterQuality, T90Variable
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Oil             = Me%State%Oil
            Wind            = Me%State%Wind
            WaterQuality    = Me%State%WQM
            T90Variable     = Me%State%T90Variable

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_            


    end subroutine GetLagrangianAirOptions  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine ModifyLagrangian(LagrangianID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_    


        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call ReadLockExternalVar

            do while (Me%ExternalVar%Now >= Me%NextCompute)

                Me%Now = Me%NextCompute

                if (Me%OverLay%Overlay) then
                    call OverLayProcess
                endif

                call ReadLockEulerianDensity

                !Emits new particles in all origins
                call ParticleEmission       ()

                !Calculates the Particle Density
                call ParticleDensity        ()

                !Calculates the GridThickness
                if (Me%State%Oil) then                
                    call FillGridThickness      ()             
                endif

                !Moves the Origins
                call MoveOrigin             ()

                !Moves the Particles
                call MovePartic             ()

                !Removes Particles in intertidal areas
                call VerifyLandParticles    ()      

                !Eliminates Particles
                call PurgeParticles         ()

                !Increases Volume
                call VolumeVariation        ()

                !Eliminates Particles
                call PurgeParticles         ()

                !Verifies is a particle is beached
                if (Me%State%AssociateBeachProb) then
                    call VerifyParticleBeaching() 
                endif

                !Dilute Particle with Ambiente Concentration
                call Dilution                ()

                if (Me%State%WQM .or. Me%State%T90Variable) then
                    call LightEvolution      ()
                endif

                !Calculates Sinks and Sources due WQM
                if (Me%State%WQM) then
                    call PropertiesEvolution ()
                endif

                !Applies T90 to Coliform Particle
                if (Me%State%FCF) then
                    call ColiformDecay       ()
                endif

                if (Me%State%Partition) then
                    call PartitionDecay      ()
                endif

                !Internal Oil Processes
                if (Me%State%Oil) then
                    call ComputeAreaVolume  ()
                    call InternalParticOil   ()
                endif

                !Calculates the Mass of each Particle
                call NewParticleMass        ()

                !Compue Age evolution
                if (Me%State%Age)  then
                    call NewParticleAge ()
                endif

                !Monitorization
                if (Me%State%Monitor) then
                    call MonitorParticle    ()
                endif

                if (Me%State%EulerianMonitor)then
                    call MonitorEulerian    ()
                end if

                call ReadUnLockEulerianDensity


                if (Me%State%Filtration) then
                    call ActualizesMassFilterGrid()
                endif

                !Statistic
                if (Me%State%Statistics) then
                    call ModifyParticStatistic ()
                endif

                if (Me%State%Deposition) then
                    call ActualizesTauErosionGrid()
                endif

                !Writes Particle Output
                call ParticleOutput         ()

                !Writes Time Series !FM
                if (Me%WritesTimeSerie) then
                   call OutPut_TimeSeries   ()
                endif

                if(Me%Output%WriteRestartFile)then
                    call OutputRestartFile
                end if

                Me%NextCompute = Me%NextCompute + Me%DT_Partic

            enddo
            !endif

            call ReadUnLockExternalVar  ()

            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLagrangian

    !--------------------------------------------------------------------------
    subroutine OverLayProcess ()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: WorkSizeILB, WorkSizeIUB
        integer                                     :: WorkSizeJLB, WorkSizeJUB
        integer                                     :: WorkSizeKLB, WorkSizeKUB
        integer                                     :: i, j, k
        integer                                     :: STAT_CALL


        !Shorten
        WorkSizeILB = Me%ExternalVar%WorkSize%ILB
        WorkSizeIUB = Me%ExternalVar%WorkSize%IUB
        WorkSizeJLB = Me%ExternalVar%WorkSize%JLB
        WorkSizeJUB = Me%ExternalVar%WorkSize%JUB
        WorkSizeKLB = Me%ExternalVar%WorkSize%KLB
        WorkSizeKUB = Me%ExternalVar%WorkSize%KUB


        call GetAssimilationField(Me%ObjAssimilation,                 &
                                  ID          = VelocityU_,                      &
                                  Field3D     = Me%ExternalVar%OverLayU,  &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangian - ERR01' 

        call GetAssimilationField(Me%ObjAssimilation,                 &
                                  ID          = VelocityV_,                      &
                                  Field3D     = Me%ExternalVar%OverLayV,  &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangian - ERR02' 

        Me%OverLay%VelUFinal = 0.
        Me%OverLay%VelVFinal = 0.

        do k = WorkSizeKLB, WorkSizeKUB
        do j = WorkSizeJLB, WorkSizeJUB
        do i = WorkSizeILB, WorkSizeIUB

            Me%OverLay%VelUFinal(i, j, k) = (Me%ExternalVar%OverLayU   (i, j  , k)     +   &
                                                        Me%ExternalVar%OverLayU   (i, j+1, k))/2. +   &
                                                        Me%ExternalVar%Velocity_U (i, j  , k)

            Me%OverLay%VelVFinal(i, j, k) = (Me%ExternalVar%OverLayV   (i  , j, k)     +   &
                                                        Me%ExternalVar%OverLayV   (i+1, j, k))/2. +   &
                                                        Me%ExternalVar%Velocity_V (i  , j, k)
        enddo
        enddo
        enddo

        call UngetAssimilation (Me%ObjAssimilation, Me%ExternalVar%OverLayU, &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangian - ERR03' 

        call UngetAssimilation (Me%ObjAssimilation, Me%ExternalVar%OverLayV, &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangian - ERR04' 


    end subroutine OverLayProcess

    !--------------------------------------------------------------------------

    subroutine ParticleEmission ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin


        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%EmissionTemporal == Continuous_) then

                if (Me%Now >= CurrentOrigin%NextEmission .and.    &
                    Me%Now <= CurrentOrigin%StopEmission) then

                    select case (CurrentOrigin%EmissionSpatial)

                    case (Point_)

                        call EmissionPoint      (CurrentOrigin)

                    case (Box_)

                        call EmissionBox        (CurrentOrigin)

                    end select

                    CurrentOrigin%NextEmission  = CurrentOrigin%NextEmission +           &
                                                  CurrentOrigin%DT_EMIT

                endif

            endif

            if (CurrentOrigin%EmissionSpatial == Accident_) then
                
                if (Me%ExternalVar%Now >= CurrentOrigin%AccidentTime .and.    &
                    .not. CurrentOrigin%AccidentFinished) then
                
                    call EmissionAccident(CurrentOrigin)

                endif

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr


    end subroutine ParticleEmission

   !--------------------------------------------------------------------------

    subroutine FillGridThickness ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        real, dimension(:, :), pointer              :: SumParticVolumeCell
        integer                                     :: i, j, STAT_CALL

        
        allocate (SumParticVolumeCell(Me%ExternalVar%Size%ILB:                &
                                      Me%ExternalVar%Size%IUB,                &
                                      Me%ExternalVar%Size%JLB:                &
                                      Me%ExternalVar%Size%JUB),               &
                                      STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FillGridThickness - ModuleLagrangian- ERR01'

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))
        
            if (CurrentOrigin%State%Oil) then
                
                !Adds the volume of all particles in the every cell
                SumParticVolumeCell = 0.
                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J

                    SumParticVolumeCell(i, j) = SumParticVolumeCell(i, j) +              & 
                                                CurrentPartic%Geometry%Volume

                    CurrentPartic => CurrentPartic%Next
                enddo


                !Calculates the GridThickness
                CurrentOrigin%GridThickness = SumParticVolumeCell /                      &
                                              Me%Grid%AreaCell

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

        deallocate (SumParticVolumeCell, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FillGridThickness - ModuleLagrangian- ERR02'


    end subroutine FillGridThickness

    !--------------------------------------------------------------------------

    subroutine OilGridConcentration (CurrentOrigin, WaveHeight, WaterDensity)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        real                                        :: WaveHeight
        real                                        :: WaterDensity

        !Local-----------------------------------------------------------------
 !       real                                        :: Specific_MSurface
 !       real                                        :: Specific_MDispersed
 !       real                                        :: MixingDepth
        
        if (CurrentOrigin%State%Oil) then
        
            !Calculates the Concentration
            ! CurrentOrigin%OilGridConcentration = (Specific_MSurface + Specific_MDispersed) / (MixingDepth * WaterDensity)
            CurrentOrigin%OilGridConcentration   = 1e6 * ((CurrentOrigin%GridThickness * Me%ExternalVar%OilDensity) + & 
                              (Me%ExternalVar%MDispersed / (max(AllmostZero,Me%ExternalVar%AreaTotal)))) / & 
                              (1.5 * WaveHeight * WaterDensity)                      

        endif


    end subroutine OilGridConcentration

    !--------------------------------------------------------------------------

     subroutine ParticleDensity () 
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin),   pointer                  :: CurrentOrigin
        type (T_Partic),   pointer                  :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: TemperatureID, SalinityID
        real                                        :: T, S, Depth
        integer                                     :: iProp
        
        !Begin-----------------------------------------------------------------

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))
        
            if (CurrentOrigin%State%FarFieldBuoyancy) then

                iProp = 0

                !Finds Temperature and Salinity
                CurrentProperty => CurrentOrigin%FirstProperty
                do while (associated(CurrentProperty))

                    iProp = iProp + 1

                    select case (CurrentProperty%ID)
                        case (Temperature_)
                            TemperatureID = iProp
                                                            
                        case (Salinity_   )
                            SalinityID    = iProp
                    end select

                    CurrentProperty => CurrentProperty%Next

                enddo
      

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))
                    T = CurrentPartic%Concentration (TemperatureID)
                    S = CurrentPartic%Concentration (SalinityID   )
                    
                    if (CurrentOrigin%Movement%DensityMethod == UNESCOState_) then
                        CurrentPartic%SigmaDensity = SigmaUNESCO     (T, S)

                    elseif (CurrentOrigin%Movement%DensityMethod == LeendertseState_) then
                        CurrentPartic%SigmaDensity = SigmaLeendertse (T, S)

                    else 
                        write(*,*)'Invalid Density Method'
                        stop      'ParticleDensity - ModuleLagrangian - ERR10'
                    endif
                    CurrentPartic => CurrentPartic%Next
                enddo
                
                    CurrentPartic => CurrentOrigin%FirstPartic
                    do while (associated(CurrentPartic))
                        T = CurrentPartic%Concentration (TemperatureID)
                        S = CurrentPartic%Concentration (SalinityID   )
                        
                        Depth = - 1 * Me%ExternalVar%ZCellCenter(CurrentPartic%Position%I, &
                                                                 CurrentPartic%Position%J, &
                                                                 CurrentPartic%Position%K)

                        if (CurrentOrigin%Movement%DensityMethod == UNESCOState_) then
                            CurrentPartic%SigmaDensity = SigmaUNESCOPressureCorrection (T, S, Depth,CurrentPartic%SigmaDensity)
               
                        elseif (CurrentOrigin%Movement%DensityMethod == LeendertseState_) then
                            write(*,*)'Invalid Density Method'
                            stop      'ParticleDensity - ModuleLagrangian - ERR20'

                        else 
                            write(*,*)'Invalid Density Method'
                            stop      'ParticleDensity - ModuleLagrangian - ERR30'
                        endif
                        CurrentPartic => CurrentPartic%Next
                    enddo                

            endif
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr     
        
     
    end subroutine ParticleDensity
       
    !--------------------------------------------------------------------------

    subroutine ActualizesTauErosionGrid ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin),   pointer                  :: CurrentOrigin
        type (T_Partic),   pointer                  :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: i, j, ILB, IUB, JLB, JUB, &
                                                       iGroup, Sediment_ID
        logical                                     :: FoundSediment

        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB

        Me%TauErosionGrid(:, :, :) = 0.

        Me%MassSedGrid   (:, :, :) = 0.

Group:  do iGroup = 1, Me%nGroups

            CurrentOrigin => Me%FirstOrigin
CurrOr:     do while (associated(CurrentOrigin))
        

                !Sediment Property 
                Sediment_ID      = 0
                FoundSediment = .false.
                CurrentProperty => CurrentOrigin%FirstProperty
                do while (associated(CurrentProperty) .and. .not. FoundSediment)
                    Sediment_ID = Sediment_ID + 1
                    if (CurrentProperty%ID == Sediment) FoundSediment = .true.
                    CurrentProperty => CurrentProperty%Next
                enddo

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J

                    if (CurrentPartic%Deposited) then

                        Me%TauErosionGrid  (iGroup, i, j) = Me%TauErosionGrid  (iGroup, i, j) + &
                                                                       CurrentPartic%Mass(Sediment_ID)              * &
                                                                       CurrentOrigin%Deposition%TauErosion

                        Me%MassSedGrid(iGroup, i, j) =                        &
                            Me%MassSedGrid(iGroup, i, j) + CurrentPartic%Mass(Sediment_ID) 

                    endif


                    CurrentPartic => CurrentPartic%Next
                enddo

                CurrentOrigin => CurrentOrigin%Next

            enddo CurrOr

            !Compute the probability of each tracer of being eroded
            !function of the Origin Erosion rate 
            CurrentOrigin => Me%FirstOrigin
CurrOr1:    do while (associated(CurrentOrigin))
        

                !Sediment Property 
                Sediment_ID      = 0
                FoundSediment = .false.
                CurrentProperty => CurrentOrigin%FirstProperty
                do while (associated(CurrentProperty) .and. .not. FoundSediment)
                    Sediment_ID = Sediment_ID + 1
                    if (CurrentProperty%ID == Sediment) FoundSediment = .true.
                    CurrentProperty => CurrentProperty%Next
                enddo

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J

                    if (CurrentPartic%Deposited) then
                        
                        CurrentPartic%ErosionRateProbability =                           &
                            min (1., CurrentOrigin%Deposition%ErosionRate *              &
                                     Me%DT_Partic              /              &
                                     Me%MassSedGrid(iGroup, i, j))
                    endif
                    CurrentPartic => CurrentPartic%Next
                enddo

                CurrentOrigin => CurrentOrigin%Next

            enddo CurrOr1


        enddo Group

        do iGroup =   1, Me%nGroups
        do i      = ILB, IUB
        do j      = JLB, JUB

            if (Me%MassSedGrid(iGroup, i, j) > 0.) then

                Me%TauErosionGrid  (iGroup, i, j) =                           &
                    Me%TauErosionGrid  (iGroup, i, j)        /                & 
                    Me%MassSedGrid(iGroup, i, j)

            else

                Me%TauErosionGrid  (iGroup, i, j) = 0.

            endif

        enddo
        enddo
        enddo


    end subroutine ActualizesTauErosionGrid

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine ActualizesMassFilterGrid ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin),   pointer                  :: CurrentOrigin
        type (T_Partic),   pointer                  :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        real,   dimension(:,:,:), pointer           :: FiltrationRateX
        real                                        :: DTEulerian, AuxFilter
        integer                                     :: i, j, k, ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: FilterProp_ID, STAT_CALL

        !Begin-----------------------------------------------------------------

ifilt:  if (Me%ExternalVar%Now >= Me%Filtration%Next) then

            ILB = Me%ExternalVar%WorkSize%ILB
            IUB = Me%ExternalVar%WorkSize%IUB
            JLB = Me%ExternalVar%WorkSize%JLB
            JUB = Me%ExternalVar%WorkSize%JUB
            KLB = Me%ExternalVar%WorkSize%KLB
            KUB = Me%ExternalVar%WorkSize%KUB

            Me%Filtration%RelativeMassFilter(:,:,:) = 0. 
            Me%Filtration%MassFiltered      (:,:,:) = 0. 

            call GetFiltrationRate (Me%ObjWaterProperties, FiltrationRateX,             &
                                    DTEulerian, Fecal_Coliforms_, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop "Lagrangian - ActualizesMassFilterGrid - ERR10"

            if (DTEulerian < Me%DT_Partic) then
                stop "Lagrangian - ActualizesMassFilterGrid - ERR20"
            endif


            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint .and.           &
                    FiltrationRateX(i, j, k) > 0.) then
                    ![ ] = [1/T] * [T] / [L3]
                    Me%Filtration%RelativeMassFilter(i,j,k) = FiltrationRateX(i, j, k) * &
                                                              DTEulerian               

                    if (Me%Filtration%RelativeMassFilter(i,j,k) > 1.)                   &
                        stop "Lagrangian - ActualizesMassFilterGrid - ERR30"
                endif

            enddo
            enddo
            enddo

            CurrentOrigin => Me%FirstOrigin
CurrOr:     do while (associated(CurrentOrigin))
    

                !Filter Property 
                FilterProp_ID      = 0
                CurrentProperty => CurrentOrigin%FirstProperty
                do while (associated(CurrentProperty))
                    FilterProp_ID = FilterProp_ID + 1
                    if (CurrentProperty%Filtration) then
                        if (CurrentProperty%ID /= Fecal_Coliforms_) then
                            write(*,*) "The filtration process in the lagrangian model is only ready for fecal coliforms"
                            stop "Lagrangian - ActualizesMassFilterGrid - ERR40"
                        endif
                        exit
                    endif

                    CurrentProperty => CurrentProperty%Next
                enddo

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J
                    k = CurrentPartic%Position%K

                    if (Me%Filtration%RelativeMassFilter(i,j,k) > 0.) then

                        AuxFilter = Me%Filtration%RelativeMassFilter(i,j,k) * CurrentPartic%Mass(FilterProp_ID)

                        Me%Filtration%MassFiltered(i, j, k) = Me%Filtration%MassFiltered(i, j, k) + AuxFilter

                        CurrentPartic%Mass(FilterProp_ID) = CurrentPartic%Mass(FilterProp_ID) - AuxFilter
                    endif

                    CurrentPartic => CurrentPartic%Next
                enddo

                CurrentOrigin => CurrentOrigin%Next

            enddo CurrOr

            call UnGetWaterProperties (Me%ObjWaterProperties, FiltrationRateX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "Lagrangian - ActualizesMassFilterGrid - ERR60"

            call SetLagrangianSinksSources(Me%ObjWaterProperties, Fecal_Coliforms_,     &
                                           Me%Filtration%MassFiltered, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "Lagrangian - ActualizesMassFilterGrid - ERR50"

            Me%Filtration%Next = Me%Filtration%Next + DTEulerian

        endif ifilt


    end subroutine ActualizesMassFilterGrid

    !--------------------------------------------------------------------------

    subroutine VerifyParticleBeaching()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: i, j, k
        real                                        :: CellI, CellJ, CellK
        real                                        :: BalX, BalY
        real                                        :: Rand1

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))


IfBeaching: if (CurrentOrigin%Beaching) then
                    
                CurrentPartic => CurrentOrigin%FirstPartic
CurrPart:       do while (associated(CurrentPartic))

                    !Grid Cell of the particle
                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J
                    k = CurrentPartic%Position%K

                    !Cell Position
                    CellI = CurrentPartic%Position%CellI
                    CellJ = CurrentPartic%Position%CellJ
                    CellK = CurrentPartic%Position%CellK

                    !Fraction of the cell
                    BALX  = CellJ - int(CellJ)
                    BALY  = CellI - int(CellI)

IfParticNotBeached: if (.NOT. CurrentPartic%Beached) then

                        call RANDOM_NUMBER(Rand1)

                        select case (Me%Grid%GridType)
                     
                            case (Grid1D)

                                if ((Me%ExternalVar%OpenPoints3D(i, j + 1, k).NE. OpenPoint)  .AND.                       &
                                      ((1-BALX) *  (Me%Grid%ParticX(j+1) - Me%Grid%ParticX(j))                 &
                                       .LT. Me%BeachingLimit)) then                                   
                       
                                    if ((Me%BeachingProbability(i,j+1,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i,j+1,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON

              
                                else if ((Me%ExternalVar%OpenPoints3D(i, j - 1, k).NE. OpenPoint)  .AND.                  &
                                           (BALX *  (Me%Grid%ParticX(j+1) - Me%Grid%ParticX(j))                &
                                            .LT. Me%BeachingLimit)) then                                  
                                                      
                                    if ((Me%BeachingProbability(i,j-1,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i,j-1,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON
                                            
                                else if ((Me%ExternalVar%OpenPoints3D(i + 1, j, k).NE. OpenPoint)  .AND.                  &
                                           ((1-BALY) *  (Me%Grid%ParticY(i+1) - Me%Grid%ParticY(i))            &
                                            .LT. Me%BeachingLimit)) then                                   

                                    if ((Me%BeachingProbability(i+1,j,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i+1,j,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON

                                else if ((Me%ExternalVar%OpenPoints3D(i - 1, j, k).NE. OpenPoint)  .AND.                  &
                                           (BALY *  (Me%Grid%ParticY(i+1) - Me%Grid%ParticY(i))                &
                                            .LT. Me%BeachingLimit)) then                                  

                                    if ((Me%BeachingProbability(i-1,j,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i-1,j,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON
                    
                                end if
                                                                        
                            case (Grid2D)

                                if ((Me%ExternalVar%OpenPoints3D(i, j+1, k).NE. OpenPoint)  .AND.                         &
                                      ((1-BALX) *  (Me%Grid%ParticXX(i, j+1) - Me%Grid%ParticXX(i, j))         &
                                       .LT. Me%BeachingLimit)) then                                   

                                    if ((Me%BeachingProbability(i,j+1,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i,j+1,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON
               
                                else if ((Me%ExternalVar%OpenPoints3D(i, j-1, k).NE. OpenPoint)  .AND.                    &
                                           (BALX *  (Me%Grid%ParticXX(i, j+1) - Me%Grid%ParticXX(i, j))        &
                                            .LT. Me%BeachingLimit)) then                                  


                                    if ((Me%BeachingProbability(i,j-1,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i,j-1,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON
             
                                else if ((Me%ExternalVar%OpenPoints3D(i+1, j, k).NE. OpenPoint)  .AND.                    &
                                           ((1-BALY) *  (Me%Grid%ParticYY(i+1, j) - Me%Grid%ParticYY(i, j))    &
                                            .LT. Me%BeachingLimit)) then                                   

                                    if ((Me%BeachingProbability(i+1,j,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i+1,j,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON
               
                                else if ((Me%ExternalVar%OpenPoints3D(i-1, j, k).NE. OpenPoint)  .AND.                    &
                                           (BALY *  (Me%Grid%ParticYY(i+1, j) - Me%Grid%ParticYY(i, j))        &
                                            .LT. Me%BeachingLimit)) then                                  

                                    if ((Me%BeachingProbability(i-1,j,k) .GT. Rand1)                                      &
                                        .OR. (Me%BeachingProbability(i-1,j,k) .EQ. 1))                                    &  
                                       CurrentPartic%Beached = ON
                        
                                end if


                        end select


                    end if IfParticNotBeached
                
                    CurrentPartic => CurrentPartic%Next

                enddo CurrPart

            end if IfBeaching
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine VerifyParticleBeaching

    !--------------------------------------------------------------------------

    subroutine MoveOrigin ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2
        real                                        :: NewValueX, NewValueY
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL


        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%MovingOrigin) then

                !Gets Value arround current instant     
                call GetTimeSerieValue(CurrentOrigin%ObjTimeSerie,                       &
                                       Me%ExternalVar%Now,                    &
                                       CurrentOrigin%MovingOriginColumnX,                &
                                       Time1, Value1, Time2, Value2, TimeCycle,          &
                                       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MoveOrigin - ModuleLagrangian - ERR01'

                if (TimeCycle) then
                    NewValueX = Value1
                else
                    !Interpolates Value for current instant
                    call InterpolateValueInTime(Me%ExternalVar%Now, Time1,    &
                                                Value1, Time2, Value2, NewValueX)
                endif
                                

                !Gets Value arround current instant    
                call GetTimeSerieValue(CurrentOrigin%ObjTimeSerie,                       &
                                       Me%ExternalVar%Now,                    &
                                       CurrentOrigin%MovingOriginColumnY,                &
                                       Time1, Value1, Time2, Value2, TimeCycle,          &
                                       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MoveOrigin - ModuleLagrangian - ERR02'

                if (TimeCycle) then
                    NewValueY = Value1
                else
                    !Interpolates Value for current instant
                    call InterpolateValueInTime(Me%ExternalVar%Now, Time1,    &
                                                Value1, Time2, Value2, NewValueY)
                endif
                 
                if (trim(CurrentOrigin%MovingOriginUnits) == trim(Char_Cells)) then

                    CurrentOrigin%Position%CellI = NewValueY
                    CurrentOrigin%Position%CellJ = NewValueX

                    !Convert Coordinates
                    call Convert_CellIJ_XY(CurrentOrigin%Position)
                    call Convert_CellIJ_IJ(CurrentOrigin%Position)


                else
                
                    CurrentOrigin%Position%X = NewValueX
                    CurrentOrigin%Position%Y = NewValueY

                    !Convert Coordinates
                    call Convert_XY_CellIJ(CurrentOrigin%Position)
                    call Convert_CellIJ_IJ(CurrentOrigin%Position)


                endif

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine MoveOrigin

    !--------------------------------------------------------------------------

    subroutine MovePartic ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        real, dimension(:, :, :), pointer           :: Temperature3D
        real                                        :: TemperatureX, DensityX
        real                                        :: VolInic
        integer                                     :: STAT_CALL
        integer                                     :: i, j, k
        integer                                     :: ThicknessGradient, Fay
        integer                                     :: SpreadingMethod


        if (Me%State%Oil) then

            !Gets the temperature and the Density from the Eulerian model
            call GetConcentration(Me%ObjWaterProperties,                   &
                                  ConcentrationX    = Temperature3D,                 &
                                  PropertyXIDNumber = Temperature_,                  &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write (*,*)'The Lagrangian Module needs the Water Temperature from the eulerian module'
                stop 'MovePartic - ModuleLagrangian - ERR01'
            endif

        endif



        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))


            !Calls Oil Active Processes to calculate Spreading Velocity
            if (CurrentOrigin%State%Oil .and. CurrentOrigin%nParticle > 0) then

                i            = CurrentOrigin%Position%I
                j            = CurrentOrigin%Position%J
                k            = CurrentOrigin%Position%K
                TemperatureX = Temperature3D (i, j, k)
                DensityX     = Me%ExternalVar%Density (i, j, k)

                !Calculates area and total volume    
                call ComputeAreaVolume ()

                !If it is a ccident the volume to be considered is the 
                !CurrentOrigin%PointVolume, on the other hand its
                !
                select case (CurrentOrigin%EmissionTemporal)

                case (Continuous_)

                    VolInic = CurrentOrigin%Flow * (CurrentOrigin%StopEmission - CurrentOrigin%StartEmission)

                case (Instantaneous_)

                    VolInic = CurrentOrigin%PointVolume

                end select            

                !Calculates the Oil Active Processes    
                call OilActiveProcesses(CurrentOrigin%ObjOil,                            &
                                        CurrentOrigin%GridThickness,                     &
                                        WaterTemperature   = TemperatureX,               &
                                        WaterDensity       = DensityX,                   &
                                        VolInic            = VolInic,                    &
                                        DT                 = Me%DT_PARTIC,    &
                                        AreaTotal          = CurrentOrigin%AreaTotal,    &
                                        STAT               = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangian - ERR03'

                !Gets Oil Spreading Parameters
                call GetOilSpreadingList(ThicknessGradient = ThicknessGradient,          &
                                         Fay               = Fay)
                call GetOilSpreading(CurrentOrigin%ObjOil,                               &    
                                     SpreadingMethod       = SpreadingMethod,            &
                                     STAT                  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangian - ERR03a'

                
                if     (SpreadingMethod == ThicknessGradient) then

                    call GetOilSpreadingVelocity(CurrentOrigin%ObjOil,                      &
                         SpreadingVelocityX = Me%ExternalVar%SpreadingVelocityX, &
                         SpreadingVelocityY = Me%ExternalVar%SpreadingVelocityY, &
                         STAT               = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangian - ERR03b'

                elseif (SpreadingMethod == Fay              ) then

                    call GetOilSpreadingVelocity(CurrentOrigin%ObjOil,                   &
                         DiffVelocity   = Me%ExternalVar%DiffVelocity,        &
                         STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangian - ERR03c'
                endif

            endif

            !Moves the particles horizontaly
            call MoveParticHorizontal(CurrentOrigin, ThicknessGradient, Fay, SpreadingMethod)

            !Ungets Spreading Velocity
            if (CurrentOrigin%State%Oil .and. CurrentOrigin%nParticle > 0) then

                if     (SpreadingMethod == ThicknessGradient) then

                    call UnGetOil (CurrentOrigin%ObjOil,                                 &
                                   Me%ExternalVar%SpreadingVelocityX,         &
                                   STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangian - ERR03d'

                    call UnGetOil (CurrentOrigin%ObjOil,                                 &
                                   Me%ExternalVar%SpreadingVelocityY,         &
                                   STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangian - ERR03e'

                endif

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr


        !Ungets Concentration from the eulerian module
        if (Me%State%Oil) then
            call UngetWaterProperties (Me%ObjWaterProperties, Temperature3D,  &
                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangian - ERR04'
        endif


    end subroutine MovePartic

    !--------------------------------------------------------------------------

    subroutine MoveParticHorizontal (CurrentOrigin, ThicknessGradient, Fay, SpreadingMethod)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: ThicknessGradient, Fay
        integer                                     :: SpreadingMethod

        !Local-----------------------------------------------------------------
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: i, j, k
        real                                        :: CellI, CellJ, CellK
        real                                        :: BalX, BalY, BalZ
        real                                        :: U, V
        real                                        :: UINT, VINT
        real                                        :: UWind, VWind
        real                                        :: U1, V1
        real                                        :: UD, VD
        real                                        :: DX, DY
        real                                        :: UOil
        real                                        :: VOil
        real                                        :: UPlume
        real                                        :: VPlume
        integer                                     :: Layers
        real                                        :: EspSup, Esp, EspInf
        real, dimension(:, :, :), pointer           :: Velocity_U, Velocity_V
        real                                        :: VelModH
        real                                        :: UStandardDeviation, VStandardDeviation
        real                                        :: TlagrangeX, TlagrangeY
        real                                        :: RAND, R1, R2
        type (T_Position)                           :: NewPosition
        integer                                     :: NewI, NewJ
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        logical                                     :: PositionCorrected
        real                                        :: TauAux
        integer                                     :: ID_Group
        real                                        :: Radius, Area
        real                                        :: VolOld, VolNew, dVol
        real                                        :: Modulus, WindX, WindY
        real                                        :: GradDWx, GradDWy, Aux
        logical                                     :: NoIntU, NoIntV

        !Shorten Var
        if (Me%OverLay%OverLay) then
            Velocity_U => Me%OverLay%VelUFinal
            Velocity_V => Me%OverLay%VelVFinal
        else
            Velocity_U => Me%ExternalVar%Velocity_U
            Velocity_V => Me%ExternalVar%Velocity_V
        endif

        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB



        CurrentPartic => CurrentOrigin%FirstPartic
CP:     do while (associated (CurrentPartic))
            
BD:         if (CurrentPartic%Beached .or. CurrentPartic%Deposited) then

                !Grid Cell of the particle
                i        = CurrentPartic%Position%I
                j        = CurrentPartic%Position%J
                ID_Group = CurrentOrigin%GroupID

                CurrentPartic%Position  = CurrentPartic%Position

                if (CurrentPartic%Deposited) then

                    TauAux = CurrentPartic%TauErosion + (Me%TauErosionGrid(ID_Group, i, j) - &
                             CurrentPartic%TauErosion) * Me%DT_Partic / CurrentOrigin%Deposition%Tdecay

                    CurrentPartic%TauErosion = max (CurrentOrigin%Deposition%TauErosion, TauAux)


                    if (VerifyRessuspension(CurrentPartic)) then
                        CurrentPartic%Deposited = .false. 
                        !The effect of other sediment classes ends when the particle is ressuspended
                        CurrentPartic%TauErosion = CurrentOrigin%Deposition%TauErosion
                    endif

                endif
            
            else BD
                

                !Grid Cell of the particle
                i     = CurrentPartic%Position%I
                j     = CurrentPartic%Position%J
                k     = CurrentPartic%Position%K

                !Cell Postition
                CellI = CurrentPartic%Position%CellI
                CellJ = CurrentPartic%Position%CellJ
                CellK = CurrentPartic%Position%CellK

                !Fraction of the cell
                BALX  = CellJ - int(CellJ)
                BALY  = CellI - int(CellI)
                BALZ  = CellK - int(CellK)


                !Linear Interpolation to obtain the velocity of the Water
                U = LinearInterpolation(Velocity_U(i,  j  ,k), Velocity_U(i,  j+1,k), Balx)
                V = LinearInterpolation(Velocity_V(i,  j  ,k), Velocity_V(i+1,j  ,k), Baly)

                !Linear Interpolation to obtain the thickness gradient
                GradDWx = LinearInterpolation(Me%ExternalVar%DWZ_Xgrad(i,  j  ,k), Me%ExternalVar%DWZ_Xgrad(i,  j+1,k), Balx)
                GradDWy = LinearInterpolation(Me%ExternalVar%DWZ_Ygrad(i,  j  ,k), Me%ExternalVar%DWZ_Ygrad(i+1,j  ,k), Baly)

                !Floating particle 
MF:             if (CurrentOrigin%Movement%Float) then 

                    !Velocity due Water
                    UINT = U                        
                    VINT = V           

                    !Velocity due wind                      
                    if (CurrentOrigin%Movement%WindOriginON) then

                        WindX = CurrentOrigin%Movement%WindX
                        WindY = CurrentOrigin%Movement%WindY

                    else

                        WindX = Me%ExternalVar%WindX(i, j)
                        WindY = Me%ExternalVar%WindY(i, j)

                    endif


                    UWind = CurrentOrigin%Movement%WindTransferCoef * WindX
                            
                    VWind = CurrentOrigin%Movement%WindTransferCoef * WindY

                    !Plume velocity
                    UPlume = 0.
                    VPlume = 0.

                    if (CurrentOrigin%State%Oil) then

                        if      (SpreadingMethod == ThicknessGradient) then
                    
                            UOil = LinearInterpolation(                                      &
                                    Me%ExternalVar%SpreadingVelocityX(i, j  ),    &
                                    Me%ExternalVar%SpreadingVelocityX(i, j+1),    &
                                    Balx)

                            VOil = LinearInterpolation(                                      &
                                    Me%ExternalVar%SpreadingVelocityY(i, j  ),    &
                                    Me%ExternalVar%SpreadingVelocityY(i+1, j),    &
                                    Baly)
                                                
                        elseif  (SpreadingMethod == Fay              ) then

                            call RANDOM_NUMBER(R1)
                            call RANDOM_NUMBER(R2)

                            UOil = R1 * cos(2 * PI * R2) * Me%ExternalVar%DiffVelocity
                            VOil = R1 * sin(2 * PI * R2) * Me%ExternalVar%DiffVelocity
                            
                        else

                            UOil  = 0.0
                            VOil  = 0.0

                        endif

                     else

                        UOil  = 0.0
                        VOil  = 0.0
                
                    end if

                else MF

                    UWind = 0.0
                    VWind = 0.0

                    UOil  = 0.0
                    VOil  = 0.0


                    !Interpolate horizontal velocities in the vertical
                    Layers = Me%ExternalVar%WorkSize%KUB -                        &
                             Me%ExternalVar%WorkSize%KLB + 1 

                    if (Layers > 1) then !More than 1 layer

                        EspSup = Me%ExternalVar%DWZ(i, j, k+1)
                        Esp    = Me%ExternalVar%DWZ(i, j, k  )
                        EspInf = Me%ExternalVar%DWZ(i, j, k-1)

                        NoIntU = .false.
                        NoIntV = .false. 

                        if (Me%ExternalVar%ComputeFaces3D_U(i, j,   k) /= Covered .or.  &
                            Me%ExternalVar%ComputeFaces3D_U(i, j+1, k) /= Covered ) then
                           !No vertical interpolation
                            NoIntU = .true.
                        endif

                        if (Me%ExternalVar%ComputeFaces3D_V(i, j, k  ) /= Covered .or.  &
                            Me%ExternalVar%ComputeFaces3D_V(i+1, j, k) /= Covered ) then
                           !No vertical interpolation
                            NoIntV = .true.
                        endif


                        !Not Close to the bottom
                        if (k /= Me%ExternalVar%kFloor(i, j) .and. BALZ < 0.5) then
                        
                            if (NoIntU) then 
                                UINT = U
                            else

                                U1 = LinearInterpolation(Velocity_U(I,  J  ,K-1),       &
                                                         Velocity_U(I,  J+1,K-1), Balx)

                                UINT = 2.0 * (U1 * (0.5-BALZ) * Esp + U * BALZ * Esp +  &
                                       U * 0.5 * EspInf) / (Esp + EspInf)
                            endif

                            if (NoIntV) then 
                                VINT = V
                            else

                                V1 = LinearInterpolation(Velocity_V(I,  J  ,K-1),       &
                                                         Velocity_V(I+1,J  ,K-1), Baly)

                                VINT = 2.0 * (V1 * (0.5-BALZ) * Esp + V * BALZ * Esp +  &
                                       V * 0.5 * EspInf) / (Esp + EspInf)
                            endif

                        !Not Close to the surface
                        else if ((K /= Me%ExternalVar%WorkSize%KUB) .and. (BALZ > 0.5 )) then

                            if (NoIntU) then 
                                UINT = U
                            else
                                U1 = LinearInterpolation(Velocity_U(I,  J  ,K+1),       &
                                                         Velocity_U(I,  J+1,K+1), Balx)

                                UINT = 2.0 * (U * 0.5 * EspSup + U * (1.0-BALZ) * Esp + &
                                       U1 * (BALZ-0.5) * Esp) / (Esp + EspSup)
                            endif

                            if (NoIntV) then 
                                VINT = V
                            else
                                V1 = LinearInterpolation(Velocity_V(I,  J  ,K+1),       &
                                                         Velocity_V(I+1,J  ,K+1), Baly)

                                VINT = 2.0 * (V * 0.5 * EspSup + V * (1.0-BALZ) * Esp + &
                                       V1 * (BALZ-0.5) * Esp) / (Esp + EspSup)
                            endif

                        !Close to surface / Bottom
                        else

                            UINT = U      
                            VINT = V      

                        end if

                    !1 Layer
                    else 

                        UINT = U          
                        VINT = V               

                    end if

                end if MF


MT:             if (CurrentOrigin%Movement%MovType == SullivanAllen_) then 
            
                    VelModH  = abs(cmplx(U, V))

                    UStandardDeviation = CurrentOrigin%Movement%VarVelHX * VelModH +         &
                                         CurrentOrigin%Movement%VarVelH
                    VStandardDeviation = CurrentOrigin%Movement%VarVelHX * VelModH +         &
                                         CurrentOrigin%Movement%VarVelH

                    if (UStandardDeviation > 0.0) then
                        TlagrangeX = Me%ExternalVar%MixingLengthX(i, j, k) /      &
                                     UStandardDeviation
                    else
                        TlagrangeX = 0.0     
                    endif

                    if (CurrentPartic%TpercursoX >= TlagrangeX) then  
                        call random_number(RAND)
                        !SQRT(3.0)=1.732050808 
                        !UD                       = (1.0 - 2.0 * RAND) * 1.732050808 *        &
                        !                           UStandardDeviation 

                        UD                       = 1.732050808 * UStandardDeviation * RAND

                        call random_number(RAND)
                        Aux                      = (1.0 - 2.0 * RAND)
                        
                        if (Aux >= GradDWx)   UD = - UD 

                        CurrentPartic%TpercursoX = Me%DT_Partic
                        CurrentPartic%UD_old     = UD                                  
                                      
                    else
                        UD                       = CurrentPartic%UD_old                               
                        CurrentPartic%TpercursoX = CurrentPartic%TpercursoX + Me%DT_Partic
                    end if


                    if (VStandardDeviation > 0.0) then
                        TlagrangeY = Me%ExternalVar%MixingLengthY(i, j, k) /      &
                                     VStandardDeviation
                    else
                        TlagrangeY = 0.0
                    endif

                    if (CurrentPartic%TpercursoY >= TlagrangeY) then
                        call random_number(RAND)
                        !VD                       = (1.0 - 2.0 * RAND) * 1.732050808 *        &
                        !                           VStandardDeviation  

                        VD                       = 1.732050808 * VStandardDeviation  * RAND

                        call random_number(RAND)
                        Aux                      = (1.0 - 2.0 * RAND)
                        
                        if (Aux >= GradDWy)   VD = - VD 

                        CurrentPartic%TpercursoY = Me%DT_Partic
                        CurrentPartic%VD_old     = VD                                    
                    else
                        VD                       = CurrentPartic%VD_old                            
                        CurrentPartic%TpercursoY = CurrentPartic%TpercursoY + Me%DT_Partic
                    endif
                       
                else if (CurrentOrigin%Movement%MovType .EQ. NotRandom_    ) then MT

                    UD = 0.0   
                    VD = 0.0    

                end if MT

                !Velocity due plume
                if (CurrentOrigin%State%PlumeShear) then    
                    Radius = (0.75 * CurrentPartic%Geometry%Volume/Pi) ** 0.33333
                    Area   = Pi * Radius ** 2.
                    VolOld = CurrentPartic%Geometry%Volume - CurrentPartic%Geometry%VolVar
                    VolNew = CurrentPartic%Geometry%Volume
                    dVol   = CurrentPartic%Geometry%VolVar

                    !Acceleration due drag force
                    Modulus = sqrt((UINT - CurrentPartic%U)**2. + (VINT - CurrentPartic%V)**2.)

                    
                                       !Momentum diffusion by small scale turbulence
                    CurrentPartic%U =  CurrentPartic%U * VolOld / VolNew + dVol * UINT / VolNew 

                                       !Momentum diffusion by small scale turbulence
                    CurrentPartic%V =  CurrentPartic%V * VolOld / VolNew + dVol * VINT / VolNew 


                    CurrentPartic%RelU = CurrentPartic%U - UINT
                    CurrentPartic%RelV = CurrentPartic%V - VINT

                    CurrentPartic%SDU  = UStandardDeviation
                    CurrentPartic%SDV  = VStandardDeviation

                endif


                if (CurrentOrigin%Movement%Advection) then

                    if (CurrentOrigin%State%PlumeShear) then    
                        DX = CurrentPartic%U *  Me%DT_Partic
                        DY = CurrentPartic%V *  Me%DT_Partic
                        !large scale turbulence
                        DX = DX + UD *  Me%DT_Partic
                        DY = DY + VD *  Me%DT_Partic
                    else
                        DX = (UINT + UD + UWIND + UOIL) * Me%DT_Partic
                        DY = (VINT + VD + VWIND + VOIL) * Me%DT_Partic
                    endif
                else
                    DX =         UD                 * Me%DT_Partic
                    DY =         VD                 * Me%DT_Partic
                endif

            
                !New Position
                NewPosition%X = CurrentPartic%Position%X + DX
                NewPosition%Y = CurrentPartic%Position%Y + DY

                !Convert Coordinates
                call Convert_XY_CellIJ(NewPosition)
                call Convert_CellIJ_IJ(NewPosition)

                !Verifies new position
                NewI = NewPosition%i
                NewJ = NewPosition%j

                !If NewI or NewJ are outside the limits of domain, kill the particle
                if (NewI > WS_IUB .or. NewI < WS_ILB .or. NewJ > WS_JUB .or. NewJ < WS_JLB) then
            
                    CurrentPartic%KillPartic = ON 

                !If NewI or NewJ are outside the limits of domain, kill the particle
                elseif (Me%ExternalVar%BoundaryPoints2D(NewI, NewJ) == Boundary) then

                    CurrentPartic%KillPartic = ON

                !If it isnt a OpenPoint, donesnt move, reset TPercurso
                elseif  (Me%ExternalVar%OpenPoints3D(NewI, NewJ, k) /= OpenPoint)    then

                    !If a particle doesnt move the freazed state is ON
                    CurrentPartic%Freazed    = ON

                    CurrentPartic%TpercursoX = abs(null_real)
                    CurrentPartic%TpercursoY = abs(null_real)

                !Moves it 
                else

                    call MoveParticVertical  (CurrentOrigin, CurrentPartic, NewPosition, VelModH)

                    call Convert_Z_CellK (CurrentOrigin, NewPosition, PositionCorrected)
                    call Convert_CellK_K (               NewPosition)

                    !If PositionCorrected initialize the vertical random dispersion in the next step
                    if (PositionCorrected) then
                        CurrentPartic%TpercursoZ = abs(null_real)
                        CurrentPartic%W = 0.
                    endif

                    !If a particle moved the freazed state is OFF
                    CurrentPartic%Freazed   = OFF

                    !Stores New Position
                    CurrentPartic%Position  = NewPosition

                    if (CurrentOrigin%State%Deposition) then
                    
                        CurrentPartic%Deposited = VerifyDeposition(CurrentOrigin, CurrentPartic)

                    endif

                endif

            end if BD

            CurrentPartic => CurrentPartic%Next

        enddo CP


    end subroutine MoveParticHorizontal

    !--------------------------------------------------------------------------

    logical function VerifyRessuspension(CurrentPartic)

        !Arguments-------------------------------------------------------------
   
        type (T_Partic), pointer                    :: CurrentPartic
        !Local-----------------------------------------------------------------
        real, save                                  :: Ranval
        logical                                     :: VerifyShearStress, VerifyErosionRate
        integer                                     :: i, j

        i = CurrentPartic%Position%I
        j = CurrentPartic%Position%J

        VerifyRessuspension = OFF

        VerifyShearStress   = OFF
        VerifyErosionRate   = OFF

        !Partheniades, E., 1965. Erosion and deposition of cohesive soils.
        !J. Hydr. Div., ASCE 91 HY1 , 105–139.
     
        !The ressuspension of a tracer is function of the bottom shear stress
        !and of the erosion rate. The erosion rate (Er) is quantifiied in the form of
        !a probability that is equal to  min (1, Er * dt / MassSedTotal(i,j))
cd1:    if (Me%ExternalVar%BottomStress(CurrentPartic%Position%I,             &
                                                   CurrentPartic%Position%J) >=          &
                                                   CurrentPartic%TauErosion) then

            VerifyShearStress = ON

         endif cd1

        if (VerifyShearStress) then

            call random_number(Ranval)

            if (CurrentPartic%ErosionRateProbability > RanVal) VerifyErosionRate = ON


        endif

        if (VerifyShearStress .and. VerifyErosionRate) VerifyRessuspension = ON
       
    end function VerifyRessuspension

    !--------------------------------------------------------------------------

    logical function VerifyDeposition(CurrentOrigin, CurrentPartic)

        !Arguments-------------------------------------------------------------
   
        type (T_Origin),     pointer                :: CurrentOrigin
        type (T_Partic),     pointer                :: CurrentPartic
        !Local-----------------------------------------------------------------
        real                                        :: Aux
        real, save                                  :: Ranval
        integer                                     :: i, j, kbottom

        i = CurrentPartic%Position%I
        j = CurrentPartic%Position%J
        !k = CurrentPartic%Position%K

        kbottom = Me%ExternalVar%KFloor(i, j)

        VerifyDeposition = OFF

cd1:    if (CurrentOrigin%Deposition%BottomDistance    >=                               &
            (Me%ExternalVar%SZZ(i, j, kbottom -1)- CurrentPartic%Position%Z)) then

!Odd, N.V.M., Owen, M.W., 1972. A two-layer model of mud transport in the Thames estuary. 
!In: Proceedings. Institution of Civil Engineers, London, pp. 195-202.

cd2:        if (Me%ExternalVar%BottomStress(i,j) <                                      &
                CurrentOrigin%Deposition%TauDeposition) then
                            
                call random_number(Ranval)

                Aux = (CurrentOrigin%Deposition%TauDeposition -                         &
                       Me%ExternalVar%BottomStress(i,j)) /                              &
                       CurrentOrigin%Deposition%TauDeposition

                !if BottomStress = 0 => Aux = 1 and Aux is always > Randval
                !In this case the probability of the particle be deposited is 100%
                if (Aux > Ranval) VerifyDeposition = ON

            endif cd2

        endif cd1

    end function VerifyDeposition


    real function LinearInterpolation(Val1, Val2, Coef)
        
        !Arguments-------------------------------------------------------------
        real                                        :: Val1, Val2, Coef
            
        !Begin-----------------------------------------------------------------

        LinearInterpolation = Val1 * (1.0-Coef) +  Val2 * Coef

    end function LinearInterpolation    

    !--------------------------------------------------------------------------

    subroutine MoveParticVertical(CurrentOrigin, CurrentPartic,           &
                                  NewPosition, VelModH)

        !Arguments-------------------------------------------------------------
   
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        type (T_Position)                           :: NewPosition
        real, intent (IN )                          :: VelModH

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, KUB
        real                                        :: CellK, BALZ
        real                                        :: CompZ_Up, CompZ_Down
        real                                        :: CompZ1_Up, CompZ1_Down
        real                                        :: AuxCompMisturaZ_Up, AuxCompMisturaZ_Down
        real                                        :: EspSup, Esp, EspInf
        real                                        :: VELQZ, D50M, VQ1
        real                                        :: WStandardDeviation
        real                                        :: W, WD
        real                                        :: Radius, Area, VolOld, VolNew, dVol
        real                                        :: ai, dw, Cd, DeltaD, AuxW, dwt        

        !------------------------------------------------------------------------

        i       = CurrentPartic%Position%I
        j       = CurrentPartic%Position%J
        k       = CurrentPartic%Position%K

        CellK   = CurrentPartic%Position%CellK
        KUB     = Me%ExternalVar%WorkSize%KUB
        BALZ    = CellK - int(CellK)


MF:     if (CurrentOrigin%Movement%Float) then

            NewPosition%Z = Me%ExternalVar%SZZ(i, j, KUB)

        else MF

SA:         if (CurrentOrigin%Movement%MovType == SullivanAllen_) then

                CompZ_Up   = Me%ExternalVar%Lupward  (i, j, k)
                CompZ_Down = Me%ExternalVar%Ldownward(i, j, k)

                if ((CompZ_Up < 0.0) .OR. (CompZ_Down < 0.0)) then
                    CompZ_Up   = 0.0
                    CompZ_Down = 0.0
                end if
    
                if ((Me%ExternalVar%WorkSize%KUB -                                     &
                     Me%ExternalVar%WorkSize%KLB + 1) > 1) then 

                    EspSup = Me%ExternalVar%DWZ(i, j, k+1)
                    Esp    = Me%ExternalVar%DWZ(i, j, k  )
                    EspInf = Me%ExternalVar%DWZ(i, j, k-1)

                    if ((K /= Me%ExternalVar%WorkSize%KLB) .AND.              &
                        (BALZ < 0.5)) then       !Close to the bottom

                        CompZ1_Up   = Me%ExternalVar%Lupward  (i, j, k-1)
                        CompZ1_Down = Me%ExternalVar%Ldownward(i, j, k-1)

                        AuxCompMisturaZ_Up   = 2.0 * (CompZ1_Up * (0.5 - BALZ) * Esp +   &
                                               CompZ_Up * BALZ * Esp + CompZ_Up * 0.5 *  &
                                               EspInf) / (Esp + EspInf)

                        AuxCompMisturaZ_Down = 2.0 * (CompZ1_Down * (0.5 - BALZ) * Esp + &
                                               CompZ_Down * BALZ * Esp + CompZ_Down *    &
                                               0.5 * EspInf) / (Esp + EspInf)
                 
                    else if ((K /= Me%ExternalVar%WorkSize%KUB) .AND.                  &
                             (BALZ .GT. 0.5 )) then !Close to the surface

                        CompZ1_Up   = Me%ExternalVar%Lupward  (i, j, k+1)
                        CompZ1_Down = Me%ExternalVar%Ldownward(i, j, k+1)

                        AuxCompMisturaZ_Up   = 2.0 * (CompZ1_Up * (0.5 - BALZ) * Esp +   &
                                               CompZ_Up * BALZ * Esp + CompZ_Up * 0.5 *  &
                                               EspInf) / (Esp + EspInf)

                        AuxCompMisturaZ_Down = 2.0 * (CompZ1_Down * (0.5 - BALZ) * Esp + &
                                               CompZ_Down * BALZ * Esp + CompZ_Down *    &
                                               0.5 * EspInf) / (Esp + EspInf)
                    else

                        AuxCompMisturaZ_Up   = CompZ_Up
                        AuxCompMisturaZ_Down = CompZ_Down           

                    endif

                else !1 Layer

                    AuxCompMisturaZ_Up   = CompZ_Up
                    AuxCompMisturaZ_Down = CompZ_Down           

                endif

            endif SA



DP:         if (CurrentOrigin%State%Sedimentation) then

                if (CurrentOrigin%Movement%SedimentationType .EQ. Stokes_ ) then

                    D50M  =  CurrentOrigin%Movement%D50 / 1000.0 ! passagem para metros 
                    VQ1   =  0.4949 * LOG10(D50M)**2.0 + 2.1795 * LOG10(D50M) + 3.7394  
                    VELQZ =  -1.0 / (10.0**VQ1)                                          

                    if (abs(CurrentOrigin%Movement%MinSedVel) > ABS(VELQZ))              &
                        VELQZ =  -1. * CurrentOrigin%Movement%MinSedVel

                else if (CurrentOrigin%Movement%SedimentationType .EQ. Imposed_) then

                    VELQZ =  -1. * CurrentOrigin%Movement%SedVel

                else

                    VELQZ = 0.0

                end if

            else

                VELQZ = 0.0

            end if DP


MT:         if (CurrentOrigin%Movement%MovType .EQ. SullivanAllen_) then

                select case (CurrentOrigin%Movement%StandardDeviationType)
                
                case (VerticalTurbConstant) 

                    WStandardDeviation = CurrentOrigin%Movement%VarVelVX * VelModH +     &
                                         CurrentOrigin%Movement%VarVelV

                    WD = WD_(CurrentPartic, WStandardDeviation, AuxCompMisturaZ_Up,      &
                             AuxCompMisturaZ_Down, Me%DT_Partic)

                case (VerticalTurb)

                    WStandardDeviation = 1.0975 * Me%ExternalVar%ShearVelocity(i, j) 

                    WD = WD_(CurrentPartic, WStandardDeviation, AuxCompMisturaZ_Up,      &
                             AuxCompMisturaZ_Down, Me%DT_Partic)

                end select

            else if (CurrentOrigin%Movement%MovType .EQ. NotRandom_    ) then MT

                WD = 0.0

            end if MT


            if (CurrentOrigin%Movement%Advection) then

                if (Me%ExternalVar%WorkSize%KUB -                             &
                    Me%ExternalVar%WorkSize%KLB + 1 == 1) then

                    W = 0.

                else

                    W = Me%ExternalVar%Velocity_W(i, j, k  ) * (1.0 - BALZ) + &
                        Me%ExternalVar%Velocity_W(i, j, k+1) * BALZ 

                endif

            else
                
                W = 0.0

            end if

            !Velocity due plume
PL:         if (CurrentOrigin%State%FarFieldBuoyancy) then    
                Radius = (0.75 * CurrentPartic%Geometry%Volume/Pi) ** 0.33333
                Area   = 4 * Pi * Radius ** 2.
                VolOld = CurrentPartic%Geometry%Volume - CurrentPartic%Geometry%VolVar
                VolNew = CurrentPartic%Geometry%Volume
                dVol   = CurrentPartic%Geometry%VolVar

                !Acceleration due density gradient
                DeltaD = Me%ExternalVar%SigmaDensity(i, j, k) - CurrentPartic%SigmaDensity

                !Buoyancy is consider null for low density gradients 
                if (abs(DeltaD) > 0.5) then
                    ai     = Gravity * DeltaD / (Me%ExternalVar%SigmaDensity(i, j, k) + &
                                                 SigmaDensityReference)
                else                
                    ai     = 0.
                endif

                dw     = CurrentPartic%W - W

                ! Taken from CORMIX
                Cd = 0.055 

                AuxW = Cd * abs(dw) * Area / VolNew
                
                                   !Momentum diffusion by small scale turbulence

                if (ai > 0) then
                    CurrentPartic%W =  (CurrentPartic%W * VolOld / VolNew + dVol * W / VolNew + &
                                        !buoyancy + shear friction 
                                        Me%DT_Partic * (ai + W * AuxW)) / (1 + AuxW)
                else
                    CurrentPartic%W =   CurrentPartic%W * VolOld / VolNew + dVol * W / VolNew + &
                                        !buoyancy + shear friction 
                                        Me%DT_Partic * (ai + dw * AuxW)
                endif

                dwt = Me%ExternalVar%DWZ(i, j, k) / Me%DT_Partic

                if (abs(CurrentPartic%W) > dwt  .and. abs(CurrentPartic%W) > 10. * WD) then
                
                    CurrentPartic%W = CurrentPartic%W * 10. * abs(WD) / abs(CurrentPartic%W) 

                endif

                                                          !tracer vel. + large scale turbulence 
                NewPosition%Z = CurrentPartic%Position%Z - (CurrentPartic%W + WD) *  Me%DT_Partic

            else PL

                CurrentPartic%W    = 0.

                NewPosition%Z = CurrentPartic%Position%Z - (W + WD + VELQZ) *  Me%DT_Partic

            endif PL

        endif MF

        !------------------------------------------------------------------------

    end subroutine MoveParticVertical

    !--------------------------------------------------------------------------

    real function WD_(CurrentPartic, WStandardDeviation, AuxCompMisturaZ_Up,             &
                      AuxCompMisturaZ_Down, DT_Partic)

        !Arguments-------------------------------------------------------------
        type (T_Partic), pointer                    :: CurrentPartic
        real,    intent(IN)                         :: WStandardDeviation
        real,    intent(IN)                         :: AuxCompMisturaZ_Up
        real,    intent(IN)                         :: AuxCompMisturaZ_Down
        real,    intent(IN)                         :: DT_Partic

        !Local-----------------------------------------------------------------
        real                                        :: TlagrangeZ
        real                                        :: RAND


        !Calculate TlagrangeZ
        if (WStandardDeviation > 0.0) then
            if (CurrentPartic%WD_old > 0.0) then 

                TlagrangeZ = AuxCompMisturaZ_Up   / WStandardDeviation

            else

                TlagrangeZ = AuxCompMisturaZ_Down / WStandardDeviation

            end if 

        else 

            TlagrangeZ = 0.0

        end if 


        if (CurrentPartic%TpercursoZ >= TlagrangeZ) then
            call random_number(RAND)
                                                                   !SQRT(3.0)=1.732050808 
            WD_                      = (1.0 - 2.0 * RAND) * 1.732050808 * WStandardDeviation
            CurrentPartic%TpercursoZ = DT_Partic
            CurrentPartic%WD_old     = WD_       
                                                          
        else

            WD_                      = CurrentPartic%WD_old                             
            CurrentPartic%TpercursoZ = CurrentPartic%TpercursoZ + DT_Partic
        end if

        !----------------------------------------------------------------------

    end function WD_

    !--------------------------------------------------------------------------

    subroutine VolumeVariation ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: i, j, k
        real                                        :: CellI, CellJ
        real                                        :: Balx, Baly
        real                                        :: U, V, VelQ
        real                                        :: FKK
        real                                        :: VolOld
        real                                        :: PlumeCoef, RelativeVel, VarianceTurb
        real, dimension(:, :, :), pointer           :: Velocity_U, Velocity_V

        !Shorten Var
        if (Me%OverLay%OverLay) then
            Velocity_U => Me%OverLay%VelUFinal
            Velocity_V => Me%OverLay%VelVFinal
        else
            Velocity_U => Me%ExternalVar%Velocity_U
            Velocity_V => Me%ExternalVar%Velocity_V
        endif
        

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%State%VariableGeom) then

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

DB:                 if (.not. CurrentPartic%Deposited .and.                             &
                        .not. CurrentPartic%Beached   .and.                             &
                        .not. CurrentPartic%Freazed        ) then
    
                        !Grid Cell of the particle
                        i     = CurrentPartic%Position%I
                        j     = CurrentPartic%Position%J
                        k     = CurrentPartic%Position%K

                        !Cell Postition
                        CellI = CurrentPartic%Position%CellI
                        CellJ = CurrentPartic%Position%CellJ

                        !Fraction of the cell
                        BALX  = CellJ - int(CellJ)
                        BALY  = CellI - int(CellI)

                        !Linear Interpolation to obtain the velocity of the Water
                        U = LinearInterpolation(Velocity_U(i,  j  ,k), Velocity_U(i,  j+1,k), Balx)
                        V = LinearInterpolation(Velocity_V(i,  j  ,k), Velocity_V(i+1,j  ,k), Baly)

                        !Velocity ^ 2
                        VelQ    = U**2. + V**2.
                    
                        !ZKMIN    =  0.0
                        !HREF     = 10.0
                        !VSQRREF  = 1*1   ! quadrado da velocidade
                        !ZKVOL    = LOG(2.)/TVOL200
                        !FK = ZKMIN + ZKVOL * (H/HREF) * (VSQR/VSQRREF)
                        select case (CurrentOrigin%Movement%TVolType)

                        case (Double_)

                            FKK     = (log(2.)/ CurrentOrigin%Movement%TVOL200)

                        case (Velocity_)

                            FKK     = (log(2.)/ CurrentOrigin%Movement%TVOL200) *                &  !ZKVOL
                                      (Me%ExternalVar%WaterColumn(i, j)/10.0) *       &  !(H/HREF)
                                      (VelQ)                                                        !(VSQR/VSQRREF)

                        end select

                        if (CurrentOrigin%State%PlumeShear) then
                        
                            !The vertical velocity is not considerer because in cases
                            !where buoancy is very high the Vertical velocity 
                            !is artificially high 
                            !In the far field is assumed that the turbulence is mainly horizontal
                            RelativeVel  = CurrentPartic%RelU**2 + CurrentPartic%RelV**2
                            RelativeVel  = RelativeVel * CurrentOrigin%Movement%CoefInitialMixing
                            VarianceTurb = CurrentPartic%SDU**2  + CurrentPartic%SDV**2 

                            if (VarianceTurb> abs(null_real)) then
                                write(*,*)'CurrentPartic%SDU', CurrentPartic%SDU
                                write(*,*)'CurrentPartic%SDV', CurrentPartic%SDV
                                stop 'Lagrangian - VolumeVariation - ERR01'
                            endif 
  
                           ! 10% of the horizontal velocity
                            if (VarianceTurb == 0) then
                               ! 10% of the horizontal velocity
                                VarianceTurb = max(1e-3,0.01 * VelQ)
                            endif

                            PlumeCoef = (RelativeVel + VarianceTurb) / VarianceTurb

                            FKK = FKK * PlumeCoef

                        endif 

                        !Old Volume
                        VolOld  = CurrentPartic%Geometry%Volume
                    
                        !New Volume
                        CurrentPartic%Geometry%Volume = CurrentPartic%Geometry%Volume *      &
                                                        (1. + Me%DT_Partic * FKK) ! Explicito

                        !Volume Variation
                        CurrentPartic%Geometry%VolVar = CurrentPartic%Geometry%Volume - VolOld

                        !Verifies dead volume
                        if (CurrentPartic%Geometry%Volume >                                  &
                            CurrentPartic%Geometry%InitialVolume * CurrentOrigin%Movement%VOLFAC) then
                        
                            CurrentPartic%KillPartic = ON
                        endif

                    endif DB

                    CurrentPartic => CurrentPartic%Next

                enddo

            endif
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

        

    end subroutine VolumeVariation

    !--------------------------------------------------------------------------

    subroutine Dilution ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Property), pointer                  :: CurrentProperty
        type (T_Partic), pointer                    :: CurrentPartic
        real                                        :: Concentration
        integer                                     :: nProp
        real                                        :: OldVolume

        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%State%VariableGeom) then

                nProp = 1
                CurrentProperty => CurrentOrigin%FirstProperty
                do while (associated(CurrentProperty))

                    CurrentPartic => CurrentOrigin%FirstPartic
                    do while (associated(CurrentPartic))
                    
DB:                     if (.not. CurrentPartic%Deposited .and.                         &
                            .not. CurrentPartic%Beached   .and.                         &
                            .not. CurrentPartic%Freazed        ) then
                    

                            call GetAmbientConcentration (CurrentProperty,              &
                                                          CurrentPartic%Position,       &
                                                          Concentration)

                            if (CurrentPartic%Geometry%VolVar < 0) then

                                stop 'Dilution - ModuleLagrangian - ERR01'

                            endif

                            OldVolume = CurrentPartic%Geometry%Volume -                 &
                                        CurrentPartic%Geometry%VolVar

                            if (nprop <1 .or. nprop > 4 .or. CurrentPartic%Geometry%Volume <0) then
                                write(*,*) CurrentOrigin%name
                                write(*,*) CurrentProperty%ID
                                write(*,*) CurrentPartic%ID
                                write(*,*) nprop
                                write(*,*) CurrentPartic%Geometry%Volume 
                            endif

                            !Calculates New Concentration
                            CurrentPartic%Concentration(nProp) =                        &
                                (CurrentPartic%Concentration(nProp) * OldVolume +       &
                                 Concentration * CurrentPartic%Geometry%VolVar) /       &
                                 CurrentPartic%Geometry%Volume

                        endif DB

                        CurrentPartic => CurrentPartic%Next
                        
                    enddo
        


                    nProp = nProp + 1
                    CurrentProperty => CurrentProperty%Next
                enddo

            endif
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr


    end subroutine Dilution

    !--------------------------------------------------------------------------

    subroutine GetAmbientConcentration (Property, Position, Concentration)

        !Arguments-------------------------------------------------------------
   
        type (T_Property), pointer                  :: Property
        type (T_Position)                           :: Position
        real                                        :: Concentration

        !Local-----------------------------------------------------------------
        real, dimension(:, :, :), pointer           :: ConcentrationX
        integer                                     :: STAT_CALL

        !If the property have a given ambient concentration, use it, else
        !try to get it from the Waterproperties module
        if (Property%HaveAmbientConcentration) then
            Concentration = Property%AmbientConcentration
        else

            call GetConcentration(Me%ObjWaterProperties,       &
                                  ConcentrationX    = ConcentrationX,    &
                                  PropertyXIDNumber = Property%ID,       &
                                  STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                
                write(*,*)'AmbientConcentration not given for : ', trim(Property%Name)
                stop 'GetAmbientConcentration - ModuleLagrangian - ERR10'

            endif

            Concentration = ConcentrationX(Position%I, Position%J, Position%K)

            call UngetWaterProperties(Me%ObjWaterProperties,                            &
                                      ConcentrationX,                                   &
                                      STAT              = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                
                stop 'GetAmbientConcentration - ModuleLagrangian - ERR20'

            endif

        endif


    end subroutine GetAmbientConcentration

    !--------------------------------------------------------------------------

    subroutine VerifyLandParticles ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: i, j, k


        !Set KILL state of all particle which are in intertidal areas to ON

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))
                
            if (CurrentOrigin%Movement%KillLandParticles) then

                CurrentPartic  => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J
                    k = CurrentPartic%Position%K

                    if (Me%ExternalVar%OpenPoints3D(i, j, k)  /= OpenPoint .and. &
                        Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then

                        CurrentPartic%KillPartic = ON

                    endif

                    CurrentPartic => CurrentPartic%Next

                enddo

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr


    end subroutine VerifyLandParticles 

    !--------------------------------------------------------------------------

    subroutine PurgeParticles ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin    => null()
        type (T_Partic), pointer                    :: CurrentPartic    => null()
        type (T_Partic), pointer                    :: ParticToDelete   => null()
        type (T_Partic), pointer                    :: ParticToContinue => null()

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))
                
            CurrentPartic  => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                !Deletes all Particle with KillPartic State on
                if (CurrentPartic%KillPartic) then
                
                    ParticToContinue => CurrentPartic%Next
                    ParticToDelete   => CurrentPartic
 
                    if(associated(CurrentPartic%Prev))then
                        CurrentPartic%Prev%Next   => CurrentPartic%Next
                    else
                        CurrentOrigin%FirstPartic => CurrentPartic%Next
                    end if

                    if(associated(CurrentPartic%Next))then
                        CurrentPartic%Next%Prev   => CurrentPartic%Prev
                    end if

                    !Delete Concentration and Mass Buffer
                    if (associated(CurrentPartic%Concentration)) then
                        deallocate(CurrentPartic%Concentration)
                        nullify   (CurrentPartic%Concentration)
                    endif

                    if (associated(CurrentPartic%Mass         )) then
                        deallocate(CurrentPartic%Mass         )
                        nullify   (CurrentPartic%Mass         )
                    endif

                    deallocate    (CurrentPartic)
                    nullify       (CurrentPartic)

                    !Decreases number of particle
                    CurrentOrigin%nParticle = CurrentOrigin%nParticle - 1

                    if (CurrentOrigin%nParticle == 0) then
                        nullify (CurrentOrigin%FirstPartic) 
                    endif


                    !call DeleteParticle (CurrentOrigin, ParticToDelete)

                    !Restarts search
                    CurrentPartic  => ParticToContinue

                else

                    CurrentPartic => CurrentPartic%Next

                endif

            enddo

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine PurgeParticles

       
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine LightEvolution ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentParticle
        type (T_Property), pointer                  :: CurrentProperty
        real,    dimension(:, :, :, :), pointer     :: GridConc
        real,    dimension(:, :, :   ), pointer     :: Concentration
        real,    dimension(:, :, :   ), pointer     :: ShortWaveExtinctionField
        real                                        :: UnitsCoef, dh1, CenterRadiation, SWPercentage
        integer                                     :: iGroup, STAT_CALL, iProp, i, j, k
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB, Kbottom
        logical                                     :: NeedsParameters, NeedsConcentrations, &
                                                       NeedsPhyto, NeedsSPM, FoundPhyto, FoundSed
        !Begin-----------------------------------------------------------------

        !Shorten
        ILB    = Me%ExternalVar%WorkSize%ILB
        IUB    = Me%ExternalVar%WorkSize%IUB
        JLB    = Me%ExternalVar%WorkSize%JLB
        JUB    = Me%ExternalVar%WorkSize%JUB
        KLB    = Me%ExternalVar%WorkSize%KLB
        KUB    = Me%ExternalVar%WorkSize%KUB


        !Updates light
d1:     do iGroup = 1, Me%nGroups


i0:         if(Me%Light%Compute)then

                call GetLightExtinctionOptions(LightExtinctionID    = Me%Light%ObjLightExtinction(iGroup),& 
                                               NeedsParameters      = NeedsParameters,       &
                                               NeedsConcentrations  = NeedsConcentrations,   &
                                               NeedsPhyto           = NeedsPhyto,            &
                                               NeedsSPM             = NeedsSPM,              &
                                               STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'LightEvolution - ModuleLagrangian - ERR01'

                call GetRadiationPercentages (Me%Light%ObjLightExtinction(iGroup),           &
                                              SWPercentage = SWPercentage, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)stop 'LightEvolution - ModuleLagrangian - ERR01'


i1:             if(NeedsConcentrations)then

                    call FillGridConcentration(iGroup, GridConc)

                    CurrentOrigin => Me%FirstOrigin 
Catch:              do while (associated(CurrentOrigin))
                        if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                            CurrentProperty  => CurrentOrigin%FirstProperty
                            exit Catch
                        endif
                        CurrentOrigin => CurrentOrigin%Next
                    enddo Catch

                    iProp = 1

                    FoundPhyto = .false.
                    FoundSed   = .false.

dw1:                do while (associated(CurrentProperty))

                        if (CurrentProperty%ID == Phytoplankton_    ) FoundPhyto = .true.
                        if (CurrentProperty%ID == Cohesive_Sediment_) FoundSed   = .true.

i2:                     if ((NeedsPhyto.and.FoundPhyto).or.(NeedsSPM.and.FoundSed)) then

                            UnitsCoef = 1e-3
                         
                            Concentration => GridConc(ILB:IUB,JLB:JUB,KLB:KUB,iProp)

i3:                         if(NeedsParameters)then

                            
                                call ModifyLightExtinctionField(LightExtinctionID   = Me%Light%ObjLightExtinction(iGroup),&
                                                                WaterPoints3D       = Me%ExternalVar%WaterPoints3D,   &
                                                                CurrentTime         = Me%ExternalVar%Now,             &
                                                                PropertyID          = CurrentProperty%ID,             &
                                                                Concentration       = Concentration,                  &
                                                                UnitsCoef           = UnitsCoef,                      &
                                                                ExtinctionParameter = CurrentProperty%ExtinctionParameter,&
                                                                STAT                = STAT_CALL)
                                if (STAT_CALL/= SUCCESS_) stop 'LightEvolution - ModuleLagrangian - ERR02'

                            else i3

                                call ModifyLightExtinctionField(LightExtinctionID   = Me%Light%ObjLightExtinction(iGroup),&
                                                                WaterPoints3D       = Me%ExternalVar%WaterPoints3D,   &
                                                                CurrentTime         = Me%ExternalVar%Now,             &
                                                                PropertyID          = CurrentProperty%ID,             &
                                                                Concentration       = Concentration,                  &
                                                                UnitsCoef           = UnitsCoef,                      &
                                                                STAT                = STAT_CALL)
                                if (STAT_CALL/= SUCCESS_) stop 'LightEvolution - ModuleLagrangian - ERR03'

                            end if i3


                        endif i2

                        CurrentProperty => CurrentProperty%Next
                        iProp = iProp + 1

                    enddo dw1

                    if (NeedsPhyto.and..not.FoundPhyto) then

                        stop 'LightEvolution - ModuleLagrangian - ERR04'

                    endif
                
                    if (NeedsSPM.and..not. FoundSed) then

                        stop 'LightEvolution - ModuleLagrangian - ERR05'

                    endif

                    deallocate(GridConc)

                else i1

                    call ModifyLightExtinctionField(LightExtinctionID   = Me%Light%ObjLightExtinction(iGroup),&
                                                    WaterPoints3D       = Me%ExternalVar%WaterPoints3D,   &
                                                    CurrentTime         = Me%ExternalVar%Now,             &
                                                    STAT                = STAT_CALL)
                    if (STAT_CALL/= SUCCESS_) stop 'LightEvolution - ModuleLagrangian - ERR06'

                end if i1

            end if i0

            call GetShortWaveExtinctionField(Me%Light%ObjLightExtinction(iGroup), ShortWaveExtinctionField, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'LightEvolution - ModuleLagrangian - ERR07'



d2:         do j = JLB, JUB
d3:         do i = ILB, IUB

                
i4:             if (Me%ExternalVar%WaterPoints3D(i, j, KUB) == WaterPoint) then 
                    
                    kbottom = Me%ExternalVar%KFloor(i, j)

                    Me%Light%TopRadiationCells(i, j, KUB, iGroup) = Me%ExternalVar%SurfaceRadiation(i, j) * SWPercentage
                   
d4:                 do k= KUB-1, kbottom - 1, -1

                        Me%Light%TopRadiationCells(i, j, k, iGroup) = Me%Light%TopRadiationCells   (i,j,k+1, iGroup) * &
                                                                      exp(-ShortWaveExtinctionField(i,j,k+1)         * &
                                                                      Me%ExternalVar%DWZ(i, j, k + 1))

                        Me%Light%TopRadiationCells(i, j, k, iGroup)= max(Me%Light%TopRadiationCells(i, j, k, iGroup), 0.001)

                    end do d4

                end if i4

            end do d3
            end do d2
        
            CurrentOrigin => Me%FirstOrigin 
dw2:        do while (associated(CurrentOrigin))
                
i5:             if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then

                    !It is admited that the tracers thickness is always smaller than the cells thickness
                    CurrentParticle => CurrentOrigin%FirstPartic
dw3:                do while (associated(CurrentParticle))

                        i = CurrentParticle%Position%I
                        j = CurrentParticle%Position%J
                        k = CurrentParticle%Position%K

                        dh1 = CurrentParticle%Position%Z - Me%ExternalVar%SZZ(i, j, k)

                        !Radiation in the center of the particle (I2)
                        !I2 = I1*exp(-K * Thickness/2) <=> I1 = I2 / exp(-K * Thickness/2)
                        !I1 = top particle face
                        !I3 = bottom particle face
                        !I3 = I2*exp(-K * Thickness/2)
                        CenterRadiation =  Me%Light%TopRadiationCells(i, j, k, iGroup) *    &
                                           exp( - ShortWaveExtinctionField(i,j,k) * dh1)

                        !dh2 = CurrentParticle%Geometry%Thickness

                        !Iaverage = (I1 - I3) / (K * Thickness) = I2 / (K * Thickness) * 
                        !           (1/exp(-K * Thickness/2) - exp(-K * Thickness/2))
                        !Iaverage = I2 / (K * Thickness) * (exp(K * Thickness/2) - exp(-K * Thickness/2)) 
                        !AverageRadiation = CenterRadiation / ShortWaveExtinctionField(i,j,k) / dh2 * &
                        !                   (exp(  ShortWaveExtinctionField(i,j,k) * dh2/2)         - &
                        !                    exp(- ShortWaveExtinctionField(i,j,k) + dh2/2)) 

                        CurrentParticle%Radiation    = CenterRadiation  
                        CurrentParticle%ShortWaveExt = ShortWaveExtinctionField(i, j, k)

                        CurrentParticle => CurrentParticle%Next
                    enddo dw3


                endif i5
                CurrentOrigin => CurrentOrigin%Next
            enddo dw2

            call UnGetLightExtinction(Me%Light%ObjLightExtinction(iGroup), ShortWaveExtinctionField, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'LightEvolution - ModuleLagrangian - ERR08'
 

        enddo d1



    end subroutine LightEvolution

!--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine PropertiesEvolution ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        real, dimension(:), pointer                 :: TemperatureX
        real, dimension(:), pointer                 :: SalinityX
        real, dimension(:), pointer                 :: SWRadiationX
        real, dimension(:), pointer                 :: LightExtCoefX
        real, dimension(:), pointer                 :: Thickness
        real, dimension(:), pointer                 :: FishFoodX
        real, dimension(:, :), pointer              :: WQMMass
        real, dimension(:, :, :), pointer           :: Temperature3D, Salinity3D
        real, dimension(:, :, :), pointer           :: FishFood3D
        integer                                     :: nParticle
        integer                                     :: i, j, k, nProp
        integer                                     :: ZooWQM
        integer                                     :: LarvaeWQM
        integer                                     :: AgeWQM
        integer                                     :: PhytoWQM
        integer                                     :: AmmoniaWQM
        integer                                     :: NitrateWQM
        integer                                     :: NitriteWQM
        integer                                     :: DONRefractoryWQM
        integer                                     :: DONNonRefractoryWQM
        integer                                     :: PartOrganicNitrogenWQM
        integer                                     :: OxygenWQM
        integer                                     :: BODWQM
        integer                                     :: BacteriaWQM
        integer                                     :: CiliateWQM
        integer                                     :: DOPRefractoryWQM
        integer                                     :: DOPNonRefractoryWQM
        integer                                     :: PartOrganicPhosphorusWQM
        integer                                     :: InorganicPhosphorusWQM
        integer                                     :: STAT_CALL
        real                                        :: ScaleRelationHV, Diameter
        type (T_TIME)                               :: Actual


        !Begin-----------------------------------------------------------------


        !Actual Time
        Actual = Me%ExternalVar%Now

        !Gets the temperature and the Salinity from the Eulerian model
        call GetConcentration(Me%ObjWaterProperties,                           &
                              ConcentrationX    = Temperature3D,                         &
                              PropertyXIDNumber = Temperature_,                          &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR02'

        call GetConcentration(Me%ObjWaterProperties,                           &
                              ConcentrationX    = Salinity3D,                            &
                              PropertyXIDNumber = Salinity_,                             &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR03'


        !Gets the FishFood from the Eulerian model
        if (Me%State%Larvae) then
            call GetConcentration(Me%ObjWaterProperties,                        &
                                  ConcentrationX    = FishFood3D,                         &
                                  PropertyXIDNumber = FishFood_,                          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR003'
        end if

        !In this version the lagrangian model calls the WQM for each origin once,
        !so origins with different inicial conditions can be studied at the same
        !time.
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))
                
            if (CurrentOrigin%State%WQM .and. Actual .ge. CurrentOrigin%NextWQMCompute) then


                !Gets the indexes from the WQM
                call GetWQPropIndex(CurrentOrigin%WaterQualityID,                                       &      
                                  Zoo                             = ZooWQM,                             &      
                                  Larvae                          = LarvaeWQM,                          &      
                                  Age                             = AgeWQM,                             &      
                                  Phyto                           = PhytoWQM,                           &      
                                  Ammonia                         = AmmoniaWQM,                         &      
                                  Nitrate                         = NitrateWQM,                         &      
                                  Nitrite                         = NitriteWQM,                         &      
                                  DissOrganicNitrogenRefractory   = DONRefractoryWQM,                   &      
                                  DONNonRefractory                = DONNonRefractoryWQM,                &      
                                  PartOrganicNitrogen             = PartOrganicNitrogenWQM,             &      
                                  Oxygen                          = OxygenWQM,                          &      
                                  BOD                             = BODWQM,                             &      
                                  Bacteria                        = BacteriaWQM,                        &      
                                  Ciliate                         = CiliateWQM,                         &      
                                  DissOrganicPhosphorusRefractory = DOPRefractoryWQM,                   &      
                                  DOPNonRefractory                = DOPNonRefractoryWQM,                &      
                                  PartOrganicPhosphorus           = PartOrganicPhosphorusWQM,           &      
                                  InorganicPhosphorus             = InorganicPhosphorusWQM,             &      
                                  STAT                            = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR01'


                !Allocates the temporary matrixes
                allocate (TemperatureX        (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR05'

                allocate (SalinityX           (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR06'

                allocate (SWRadiationX        (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR07'

                allocate (LightExtCoefX       (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR07a'

                allocate (Thickness           (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR07b'

                allocate (WQMMass(CurrentOrigin%nProperties, CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR08'

                if (CurrentOrigin%State%Larvae)  then            
                    allocate (FishFoodX       (CurrentOrigin%nParticle), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR08b'
                end if


                CurrentPartic  => CurrentOrigin%FirstPartic
                nParticle      = 1
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J
                    k = CurrentPartic%Position%K

                    !Stores Temperature
                    TemperatureX (nParticle) = Temperature3D (i, j, k)

                    !Stores Salinity
                    SalinityX    (nParticle) = Salinity3D    (i, j, k)

                   !Stores FishFood
                    if (CurrentOrigin%State%Larvae)                                        &            
                        FishFoodX (nParticle) = FishFood3D (i, j, k)
                    

                    !
                    !Pina, 2001
                    !
                    !Radiation at the point of the Particle
                    ScaleRelationHV = 100. 
                    !The tracer geometry is addimited equal to a cilinder
                    Diameter = 2 * (ScaleRelationHV * CurrentPartic%Geometry%Volume / Pi) 
                    !The relation between the tracer Diameter and is thickness is consider 
                    !proportional to the relation of turbulence scales H/V. 
                    Thickness     (nParticle)   = Diameter / ScaleRelationHV 
                    SWRadiationX  (nParticle)   = CurrentPartic%Radiation
                    LightExtCoefX (nParticle)   = CurrentPartic%ShortWaveExt

                    !Stores Mass 
                    nProp = 1
                    CurrentProperty => CurrentOrigin%FirstProperty
                    do while (associated(CurrentProperty))

                        select case (CurrentProperty%ID)

                        case (Phytoplankton_                  )
                            WQMMass(PhytoWQM,            nParticle) = CurrentPartic%Concentration (nProp)
                                                            
                        case (Zooplankton_                    )
                            WQMMass(ZooWQM,              nParticle) = CurrentPartic%Concentration (nProp)
                                                            
                        case (Larvae_                         )
                            WQMMass(LarvaeWQM,           nParticle) = CurrentPartic%Concentration (nProp)

                        case (Age_                         )
                            WQMMass(AgeWQM,              nParticle) = CurrentPartic%Concentration (nProp)
 
                        case (PON_                            )
                            WQMMass(PartOrganicNitrogenWQM, nParticle) = CurrentPartic%Concentration (nProp)

                        case (DONRefractory_                  )
                            WQMMass(DONRefractoryWQM,    nParticle) = CurrentPartic%Concentration (nProp)
                           
                        case (DONNon_Refractory_              )
                            WQMMass(DONNonRefractoryWQM, nParticle) = CurrentPartic%Concentration (nProp)
                           
                        case (Ammonia_                        )
                            WQMMass(AmmoniaWQM,          nParticle) = CurrentPartic%Concentration (nProp)

                        case (Nitrate_                        )
                            WQMMass(NitrateWQM,          nParticle) = CurrentPartic%Concentration (nProp)

                        case (Nitrite_                        )
                            WQMMass(NitriteWQM,          nParticle) = CurrentPartic%Concentration (nProp)
     
                        case (BOD_                            )
                            WQMMass(BODWQM,              nParticle) = CurrentPartic%Concentration (nProp)

                        case (Oxygen_                         )
                            WQMMass(OxygenWQM,           nParticle) = CurrentPartic%Concentration (nProp)

                        case (Ciliate_                        )
                            WQMMass(CiliateWQM,          nParticle) = CurrentPartic%Concentration (nProp)

                        case (Bacteria_                       )
                            WQMMass(BacteriaWQM,         nParticle) = CurrentPartic%Concentration (nProp)

                        case (POP_                            )
                            WQMMass(PartOrganicPhosphorusWQM,nParticle) = CurrentPartic%Concentration (nProp)

                        case (DOPRefractory_                  )
                            WQMMass(DOPRefractoryWQM,    nParticle) = CurrentPartic%Concentration (nProp)

                        case (DOPNon_Refractory_              )
                            WQMMass(DOPNonRefractoryWQM, nParticle) = CurrentPartic%Concentration (nProp)

                        case (Inorganic_Phosphorus_           )
                            WQMMass(InorganicPhosphorusWQM, nParticle) = CurrentPartic%Concentration (nProp)
    
                        end select

                        nProp           = nProp + 1
                        CurrentProperty => CurrentProperty%Next
                    enddo

                    CurrentPartic => CurrentPartic%Next
                    nParticle     = nParticle + 1
                enddo
                
                if (CurrentOrigin%State%Larvae)then

                    !Runs the WaterQualityModel    
                    call WaterQuality(CurrentOrigin%WaterQualityID,                            &
                                      Salinity                      = SalinityX,               &
                                      Temperature                   = TemperatureX,            &
                                      ShortWaveRadiation            = SWRadiationX,            &
                                      LightExtCoefField             = LightExtCoefX,           &
                                      Thickness                     = Thickness,               &
                                      Mass                          = WQMMass,                 &
                                      WQArrayLB                     = 1,                       &
                                      WQArrayUB                     = CurrentOrigin%nParticle, &
                                      FishFood                      = FishFoodX,               &
                                      STAT                          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR09a'

                else

                    !Runs the WaterQualityModel    
                    call WaterQuality(CurrentOrigin%WaterQualityID,                            &
                                      Salinity                      = SalinityX,               &
                                      Temperature                   = TemperatureX,            &
                                      ShortWaveRadiation            = SWRadiationX,            &
                                      LightExtCoefField             = LightExtCoefX,           &
                                      Thickness                     = Thickness,               &
                                      Mass                          = WQMMass,                 &
                                      WQArrayLB                     = 1,                       &
                                      WQArrayUB                     = CurrentOrigin%nParticle, &
                                      STAT                          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR09b'

                end if
                !Restores Mass
                CurrentPartic  => CurrentOrigin%FirstPartic
                nParticle      = 1
                do while (associated(CurrentPartic))

                    !Stores Mass 
                    nProp           = 1
                    CurrentProperty => CurrentOrigin%FirstProperty
                    do while (associated(CurrentProperty))

                        select case (CurrentProperty%ID)

                        case (Phytoplankton_                  )
                            CurrentPartic%Concentration (nProp) = WQMMass(PhytoWQM,  nParticle)
                                                            
                        case (Zooplankton_                    )
                            CurrentPartic%Concentration (nProp) = WQMMass(ZooWQM,    nParticle)
                                                            
                        case (Larvae_                         )
                            CurrentPartic%Concentration (nProp) = WQMMass(LarvaeWQM, nParticle)
                                                            
                        case (Age_                            )
                            CurrentPartic%Concentration (nProp) = WQMMass(AgeWQM,    nParticle)

                        case (PON_   )
                            CurrentPartic%Concentration (nProp) = WQMMass(PartOrganicNitrogenWQM, nParticle)

                        case (DONRefractory_                  )
                            CurrentPartic%Concentration (nProp) = WQMMass(DONRefractoryWQM,    nParticle)
                           
                        case (DONNon_Refractory_              )
                            CurrentPartic%Concentration (nProp) = WQMMass(DONNonRefractoryWQM, nParticle)
                           
                        case (Ammonia_                        )
                            CurrentPartic%Concentration (nProp) = WQMMass(AmmoniaWQM, nParticle)

                        case (Nitrate_                        )
                            CurrentPartic%Concentration (nProp) = WQMMass(NitrateWQM, nParticle)

                        case (Nitrite_                        )
                            CurrentPartic%Concentration (nProp) = WQMMass(NitriteWQM, nParticle)
     
                        case (BOD_                            )
                            CurrentPartic%Concentration (nProp) = WQMMass(BODWQM,     nParticle)

                        case (Oxygen_                         )
                            CurrentPartic%Concentration (nProp) = WQMMass(OxygenWQM,  nParticle)

                        case (Ciliate_                        )
                            CurrentPartic%Concentration (nProp) = WQMMass(CiliateWQM, nParticle)

                        case (Bacteria_                       )
                            CurrentPartic%Concentration (nProp) = WQMMass(BacteriaWQM,nParticle)

                        case (POP_ )
                            CurrentPartic%Concentration (nProp) = WQMMass(PartOrganicPhosphorusWQM, nParticle)

                        case (DOPRefractory_                  )
                            CurrentPartic%Concentration (nProp) = WQMMass(DOPRefractoryWQM,         nParticle)

                        case (DOPNon_Refractory_              )
                            CurrentPartic%Concentration (nProp) = WQMMass(DOPNonRefractoryWQM,      nParticle)

                        case (Inorganic_Phosphorus_           )
                            CurrentPartic%Concentration (nProp) = WQMMass(InorganicPhosphorusWQM,   nParticle) 
    
                        end select

                        nProp           = nProp + 1
                        CurrentProperty => CurrentProperty%Next
                    enddo

                    CurrentPartic => CurrentPartic%Next
                    nParticle     = nParticle + 1
                enddo


                !Deallocates the temporary matrixes
                deallocate (TemperatureX,           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR10'
                nullify (TemperatureX) 

                deallocate (SalinityX,              STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR11'
                nullify (SalinityX)

                deallocate (SWRadiationX,           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR12'
                nullify (SWRadiationX)

                deallocate (LightExtCoefX,          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR12'
                nullify (LightExtCoefX)

                deallocate (Thickness,               STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR12'
                nullify (Thickness)

                deallocate (WQMMass,                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR13'
                nullify (WQMMass)

                if (CurrentOrigin%State%Larvae)  then            
                    deallocate (FishFoodX,           STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR113'
                    nullify (FishFoodX) 
                end if
                

                CurrentOrigin%NextWQMCompute = CurrentOrigin%NextWQMCompute + CurrentOrigin%DTWQM

            end if

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

        !Ungets Concentration from the eulerian module
        call UngetWaterProperties (Me%ObjWaterProperties, Temperature3D,       &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR14'

        call UngetWaterProperties (Me%ObjWaterProperties,  Salinity3D,         &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR15'
       
        if (Me%State%Larvae)  then            
            call UngetWaterProperties (Me%ObjWaterProperties, FishFood3D,      &
                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangian - ERR115'
        end if

    end subroutine PropertiesEvolution

    !--------------------------------------------------------------------------

    subroutine ColiformDecay ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: nProp, iSal, iTemp
        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%State%FCF) then

P1:             if (Me%State%T90Variable) then

                    CurrentProperty => CurrentOrigin%FirstProperty
                    nProp =   1
                    iSal  = -99
                    iTemp = -99
                    do while (associated(CurrentProperty))

                        if (CurrentProperty%ID == Salinity_   ) iSal  = nProp
                        if (CurrentProperty%ID == Temperature_) iTemp = nProp

                        nProp           =  nProp + 1
                        CurrentProperty => CurrentProperty%Next

                    enddo

                    if (iSal == -99) then

                        write(*,*) 'The origin ='//trim(CurrentOrigin%Name)// &
                                    ' must have a salinity property to compute a variable T90' 
                        stop 'ColiformDecay - ModuleLagrangian - ERR01'

                    endif

                    if (iTemp == -99) then

                        write(*,*) 'The origin ='//trim(CurrentOrigin%Name)// &
                                    ' must have a temperature property to compute a variable T90'
                        stop 'ColiformDecay - ModuleLagrangian - ERR02'

                    endif


                endif P1


                nProp = 1
                CurrentProperty => CurrentOrigin%FirstProperty
CurrProp:       do while (associated(CurrentProperty))

                    if (CurrentProperty%ID == Fecal_Coliforms_) then

                        CurrentPartic => CurrentOrigin%FirstPartic
                        do while (associated(CurrentPartic))

                            if (CurrentProperty%T90Variable) then

                                CurrentPartic%T90 = ComputeT90(CurrentPartic, iSal, iTemp,CurrentProperty%T90Var_Method)

                            else

                                CurrentPartic%T90 = CurrentProperty%T90

                            endif
                    
                            CurrentPartic%Concentration(nProp) = CurrentPartic%Concentration(nProp) / &
                                                                 (1.0 + Me%DT_Partic *     &
                                                                 (log(10.) / CurrentPartic%T90))

                            CurrentPartic => CurrentPartic%Next
                        enddo

                        exit CurrProp

                    else
                        nProp           =  nProp + 1
                        CurrentProperty => CurrentProperty%Next
                    endif
                enddo CurrProp


            endif
    
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine ColiformDecay

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    real function ComputeT90 (CurrentPartic, iSal, iTemp, Method)

        !Arguments-------------------------------------------------------------
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: iSal, iTemp, Method
        !Local-----------------------------------------------------------------
        real                                        :: Convert, Light, Sal, Temp

        !Begin-----------------------------------------------------------------

 
        Sal        = CurrentPartic%Concentration(iSal)
        Temp       = CurrentPartic%Concentration(iTemp)


!______Mortality model selection

        If (Method .eq. Canteras) then

            ComputeT90 = ComputeT90_Canteras (Temp,Sal,CurrentPartic%Radiation)
            

        elseif (Method .eq. Chapra) then

            !Converts W in ly/hr
            Convert    = 0.086325
            Light      = Convert * CurrentPartic%Radiation

            ComputeT90 = ComputeT90_Chapra (Temp, Sal, Light)  

        else

            write (*,*) 'T90 calculation method unknown'
            stop 'ComputeT90 - ModuleLagrangian - ERR1'
        
        endif


    end function ComputeT90

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine PartitionDecay ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        real                                        :: kd, rp, rd, ConcP, ConcD, Kdp, dt
        integer                                     :: nProp
        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%State%Partition) then

                nProp = 1
                CurrentProperty => CurrentOrigin%FirstProperty
CurrProp:       do while (associated(CurrentProperty))

                    if (CurrentProperty%WaterPartition%ON .or. CurrentProperty%SedimentPartition%ON) then

                        CurrentPartic => CurrentOrigin%FirstPartic
                        do while (associated(CurrentPartic))

                            if (CurrentProperty%WaterPartition%ON .and. .not. CurrentPartic%Deposited) then
                            
                                kd = CurrentProperty%WaterPartition%Coefficient 

                                rp = kd / (1 + kd)
                                rd =  1 / (1 + kd)
                                ConcP = CurrentPartic%Concentration(nProp)
                                ConcD = CurrentProperty%WaterPartition%CoupleProp
                              
                                Kdp   = CurrentProperty%WaterPartition%TransferRate
                                dt    = Me%DT_Partic

                                ConcP = ConcP + dt * Kdp * (rp * ConcD - rd * ConcP)

                                CurrentPartic%Concentration(nProp) = ConcP

                            endif

                            if (CurrentProperty%SedimentPartition%ON .and. CurrentPartic%Deposited) then
                    
                                kd = CurrentProperty%SedimentPartition%Coefficient 

                                rp = kd / (1 + kd)
                                rd =  1 / (1 + kd)
                                ConcP = CurrentPartic%Concentration(nProp)
                                ConcD = CurrentProperty%SedimentPartition%CoupleProp
                              
                                Kdp   = CurrentProperty%SedimentPartition%TransferRate
                                dt    = Me%DT_Partic

                                ConcP = ConcP + dt * Kdp * (rp * ConcD - rd * ConcP)

                                CurrentPartic%Concentration(nProp) = ConcP

                            endif


                            CurrentPartic => CurrentPartic%Next
                        enddo

                        exit CurrProp

                    else
                        nProp           =  nProp + 1
                        CurrentProperty => CurrentProperty%Next
                    endif
                enddo CurrProp


            endif
    
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine PartitionDecay

    !--------------------------------------------------------------------------

    subroutine ComputeAreaVolume ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        integer, dimension(:, :), pointer           :: AreaFlag
        integer                                     :: ILB, IUB, JLB, JUB, i, j
        real                                        :: VolInic
        real                                        :: API

        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        
        allocate (AreaFlag(ILB:IUB, JLB:JUB))
        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            CurrentOrigin%VolTotOilBeached = 0.
            CurrentOrigin%VolTotBeached    = 0.
            CurrentOrigin%AreaTotal        = 0.
            CurrentOrigin%VolumeTotal      = 0.
            CurrentOrigin%VolumeOilTotal   = 0.
            AreaFlag                       = 0

            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                i = CurrentPartic%Position%I
                j = CurrentPartic%Position%J

                if (CurrentPartic%Beached)    then

                    CurrentOrigin%VolTotBeached     = CurrentOrigin%VolTotBeached +         &
                                                      CurrentPartic%Geometry%Volume

                    CurrentOrigin%VolTotOilBeached  = CurrentOrigin%VolTotOilBeached +      &
                                                      CurrentPartic%Geometry%VolumeOil


                else if (.NOT. CurrentPartic%Beached) then

                    if (.not. CurrentOrigin%UseTheoricArea) then
                        if (AreaFlag(i, j) == 0) then
                            CurrentOrigin%AreaTotal = CurrentOrigin%AreaTotal +              &
                                                      Me%Grid%AreaCell (i, j)
                            AreaFlag(i, j)          = 1
                        endif
                    endif

                    CurrentOrigin%VolumeTotal   = CurrentOrigin%VolumeTotal +                &
                                                  CurrentPartic%Geometry%Volume

                    CurrentOrigin%VolumeOilTotal= CurrentOrigin%VolumeOilTotal +             &
                                                  CurrentPartic%Geometry%VolumeOil


                end if

                CurrentPartic => CurrentPartic%Next
            enddo


            if (CurrentOrigin%UseTheoricArea) then

                call GetOilAPI (CurrentOrigin%ObjOil, API = API)

                select case (CurrentOrigin%EmissionTemporal)
                case (Continuous_)
                    VolInic = CurrentOrigin%VolumeTotal
                case (Instantaneous_)
                    VolInic = CurrentOrigin%PointVolume
                end select            

                CurrentOrigin%AreaTotal = F_FayArea(VolInic = VolInic, API = API)

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

        deallocate (AreaFlag)

    end subroutine ComputeAreaVolume

    !--------------------------------------------------------------------------

    subroutine InternalParticOil ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic    
        real, dimension(:, :, :), pointer           :: Temperature3D
        real, dimension(:, :, :), pointer           :: SPM3D
        integer                                     :: i, j, k
        real                                        :: WaterTemperature, WaterDensity, SPM
        real                                        :: UWIND, VWIND, Wind
        real                                        :: AtmPressure
        real                                        :: WaveHeight, WavePeriod
        real                                        :: Factor
        real                                        :: VWaterContent
        real                                        :: MDispersed
        real                                        :: OilDensity
        real                                        :: AreaTotalOUT
        real                                        :: VolumeTotalOUT, VolOld
        integer                                     :: STAT_CALL

 
         
        !Gets the temperature, the Density and the SPM from the Eulerian model
        call GetConcentration(Me%ObjWaterProperties,                            &
                              ConcentrationX    = Temperature3D,                &
                              PropertyXIDNumber = Temperature_,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR01'

        call GetSPM          (Me%ObjWaterProperties,                            &
                              SPM               = SPM3D,                        &
                              STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR02'

#ifndef _WAVES_
        call GetWaves (WavesID    = Me%ObjWaves,                                &
                       WavePeriod = Me%ExternalVar%WavePeriod,                  &
                       WaveHeight = Me%ExternalVar%WaveHeight,                  &
                       STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR03'
#endif

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

        
            if (CurrentOrigin%State%Oil .and. CurrentOrigin%nParticle > 0) then

                i       = CurrentOrigin%Position%I
                j       = CurrentOrigin%Position%J
                k       = CurrentOrigin%Position%K
        
                UWIND   = Me%ExternalVar%WindX (i, j)
                VWIND   = Me%ExternalVar%WindY (i, j)
                Wind    = abs(cmplx(UWIND, VWIND))

                AtmPressure         = Me%ExternalVar%AtmPressure(i, j)
                
                WaveHeight          = Me%ExternalVar%WaveHeight (i, j)

                WavePeriod          = Me%ExternalVar%WavePeriod (i, j)

                WaterTemperature    = Temperature3D             (i, j, k)
                WaterDensity        = Me%ExternalVar%Density    (i, j, k)
                SPM                 = SPM3D                     (i, j, k)

                if (CurrentOrigin%UseTheoricArea) then
                    CurrentOrigin%AreaTotal = -1.
                endif

                !Runs Oil Internal Processes
                call OilInternalProcesses(CurrentOrigin%ObjOil,                                   &
                                          Wind                  = Wind,                           &
                                          AtmosphericPressure   = AtmPressure,                    &
                                          WaterTemperature      = WaterTemperature,               &
                                          WaterDensity          = WaterDensity,                   &
                                          SPM                   = SPM,                            &
                                          VWaterContent         = VWaterContent,                  &
                                          MDispersed            = MDispersed,                     &
                                          OilDensity            = OilDensity,                     &
                                          VolTotOilBeached      = CurrentOrigin%VolTotOilBeached, &
                                          VolTotBeached         = CurrentOrigin%VolTotBeached,    &
                                          VolumeTotalIN         = CurrentOrigin%VolumeOilTotal,   &   
                                          VolumeTotalOUT        = VolumeTotalOUT,                 &     
                                          AreaTotal             = CurrentOrigin%AreaTotal,        & 
                                          AreaTotalOUT          = AreaTotalOUT,                   & 
                                          WaveHeight            = WaveHeight,                     &
                                          WavePeriod            = WavePeriod,                     &
                                          STAT                  = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR07'


                Me%ExternalVar%VWaterContent = VWaterContent
                Me%ExternalVar%MDispersed    = MDispersed
                Me%ExternalVar%OilDensity    = OilDensity
                Me%ExternalVar%AreaTotal     = AreaTotalOUT

                call OilGridConcentration  (CurrentOrigin, WaveHeight, WaterDensity)       

                !Modifies OilVolume
                Factor                                  =  VolumeTotalOUT / (max(AllmostZero,CurrentOrigin%VolumeOilTotal))
                CurrentPartic                           => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    if (.NOT. CurrentPartic%Beached ) then
                        
                        !Old Volume
                        VolOld  = CurrentPartic%Geometry%Volume
                    
                        !New Volume
                        CurrentPartic%Geometry%VolumeOil = CurrentPartic%Geometry%VolumeOil * Factor
                        CurrentPartic%Geometry%Volume    = CurrentPartic%Geometry%VolumeOil / &
                                                           (1 - Me%ExternalVar%VWaterContent)

                        !Volume Variation
                        CurrentPartic%Geometry%VolVar = CurrentPartic%Geometry%Volume - VolOld


                    else if (CurrentPartic%Beached) then

                        CurrentPartic%Geometry%VolVar = 0.0
         
                    end if 
         
                    CurrentPartic => CurrentPartic%Next

                enddo

            endif    

            CurrentOrigin => CurrentOrigin%Next
        enddo CurrOr

#ifndef _WAVES_
        call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WavePeriod, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR08'

        call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WaveHeight, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR09'
#endif
        !Ungets Concentration from the eulerian module
        call UngetWaterProperties (Me%ObjWaterProperties, Temperature3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR10'

        call UngetWaterProperties (Me%ObjWaterProperties,  SPM3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangian - ERR11'

    end subroutine InternalParticOil

    !--------------------------------------------------------------------------

    subroutine NewParticleMass ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%nProperties >= 1) then

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    CurrentPartic%Mass(:) = CurrentPartic%Concentration(:) *             &
                                            CurrentPartic%Geometry%Volume

                    CurrentPartic => CurrentPartic%Next
                enddo

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine NewParticleMass

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine NewParticleAge ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: nprop

        
        

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%State%Age) then

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    CurrentPartic%Age = CurrentPartic%Age + Me%DT_Partic
                
                    if (.not.CurrentOrigin%State%WQM) then

                        nprop = 1

                        CurrentProperty => CurrentOrigin%FirstProperty
                        
                        do while (associated(CurrentProperty))

                          if (CurrentProperty%ID==Age_) then

                               CurrentPartic%Concentration(nprop) = CurrentPartic%Concentration(nprop) + Me%DT_Partic/3600.

                          endif
                         
                          nprop=nprop +1
                          
                          CurrentProperty=>CurrentProperty%Next  
                     
                        enddo
                    
                    endif
                  
                    CurrentPartic => CurrentPartic%Next
                
                enddo

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine NewParticleAge

    !--------------------------------------------------------------------------

    subroutine MonitorParticle ()


        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        integer, dimension(:,:,:), pointer          :: MonitorBoxes
        integer                                     :: NumberOfBoxes, NumberOfOrigins
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k, Box, iO
        integer                                     :: STAT_CALL
        real, dimension(:), pointer                 :: AuxReal
        integer                                     :: nProp
        real(8)                                     :: AuxReal2
        real(8), dimension(:), pointer              :: SumContaminatedTracerConc
        real(8), dimension(:,:), pointer            :: SumNegativeLogMassByOrigin
        real(8), dimension(:), pointer              :: SumNegativeLogMass
        !Shorten
        ILB    = Me%ExternalVar%WorkSize%ILB
        IUB    = Me%ExternalVar%WorkSize%IUB
        JLB    = Me%ExternalVar%WorkSize%JLB
        JUB    = Me%ExternalVar%WorkSize%JUB
        KLB    = Me%ExternalVar%WorkSize%KLB
        KUB    = Me%ExternalVar%WorkSize%KUB
        
        !Get the boxes
        call GetBoxes(Me%ObjMonBox, MonitorBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MonitorParticle - ModuleLagrangian - ERR01'

        call GetNumberOfBoxes(Me%ObjMonBox, NumberOfBoxes3D = NumberOfBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MonitorParticle - ModuleLagrangian - ERR01'

        NumberOfOrigins = Me%nOrigins


        !Calculates the area of each monitoring boxes
        Me%Monitor%SurfaceBoxVolume = 0.
        do j = JLB, JUB
        do i = ILB, IUB
            Box = MonitorBoxes(i, j, KUB)
            if (Box > 0 .and. Me%ExternalVar%OpenPoints3D(i, j, KUB) == OpenPoint) then
                Me%Monitor%SurfaceBoxVolume (Box) = Me%Monitor%SurfaceBoxVolume (Box)                                          &
                                                    + (min(Me%ExternalVar%DWZ(i, j, KUB), Me%Monitor%ContaminationDepth)       &
                                                    * (Me%ExternalVar%VolumeZ   (i, j, KUB) /  Me%ExternalVar%DWZ(i, j, KUB))) 
            endif
        enddo
        enddo      
        
        !Calculates the volume of each monitoring boxes
        Me%Monitor%InstBoxVolume = 0.
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            Box = MonitorBoxes(i, j, k)
            if (Box > 0 .and. Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then
                Me%Monitor%InstBoxVolume (Box) = Me%Monitor%InstBoxVolume (Box) + &
                                                            Me%ExternalVar%VolumeZ   (i, j, k)
            endif
        enddo
        enddo
        enddo

        !Integrates the volume of each monitoring box
        do Box = 1, NumberOfBoxes
            Me%Monitor%IntgBoxVolume (Box) = Me%Monitor%IntgBoxVolume (Box) + &
                                                        Me%Monitor%InstBoxVolume (Box)
        enddo



        !Calculates the volume contributed from a given origin to the volume of a 
        !monitoring cell
        Me%Monitor%InstVolumeByOrigin = 0.
        CurrentOrigin => Me%FirstOrigin
        do while (associated(CurrentOrigin))

            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                i  = CurrentPartic%Position%I
                j  = CurrentPartic%Position%J
                k  = CurrentPartic%Position%K
               
                !Number of box in which Particle is located
                Box = MonitorBoxes(i, j, k)

                if (Box > 0 .and. Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then
                    Me%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) =     &
                        Me%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) + &
                        dble(CurrentPartic%Geometry%Volume)
                endif

                CurrentPartic => CurrentPartic%Next
            enddo

            CurrentOrigin => CurrentOrigin%Next

        enddo

MMass:  If (Me%State%MonitorPropMass) then

            ! Calculates the mass contributed from a given origin to the volume of a 
            ! monitoring group of cells

            Me%Monitor%InstMassByOrigin          = 0.
            Me%Monitor%InstLogMassByOrigin       = 0.
            Me%Monitor%NumberOfTracersFromOrigin = 0
            Me%Monitor%ContaminationProbability    (:)   = 0.
            Me%Monitor%AverageBoxContaminatedConc  (:)   = 0.
            Me%Monitor%NbrBoxContaminatedTracers   (:)   = 0
            Me%Monitor%VolBoxContaminatedTracers   (:)   = 0.

            allocate (SumContaminatedTracerConc (NumberOfBoxes))            
            allocate (SumNegativeLogMassByOrigin (NumberOfBoxes, NumberOfOrigins))            
            allocate (SumNegativeLogMass (NumberOfBoxes))
            SumContaminatedTracerConc            = 0.            
            SumNegativeLogMassByOrigin           = 0.
            SumNegativeLogMass                   = 0.
            
            CurrentOrigin => Me%FirstOrigin
            do while (associated(CurrentOrigin))


                nProp = 1
                CurrentProperty => CurrentOrigin%FirstProperty
                do while (associated(CurrentProperty))

                    if (CurrentProperty%Name == Me%Monitor%MassProperty) then
    
                        CurrentPartic => CurrentOrigin%FirstPartic
                        do while (associated(CurrentPartic))
                            i  = CurrentPartic%Position%I
                            j  = CurrentPartic%Position%J
                            k  = CurrentPartic%Position%K
           
                            !Number of box in which Particle is located
                            Box = MonitorBoxes(i, j, k)

                            if (Box > 0 .and. Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then
                                Me%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID) =     &
                                Me%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID) + &
                                (CurrentPartic%Concentration(nProp) * dble(CurrentPartic%Geometry%Volume))
                            

                                AuxReal2 = CurrentPartic%Concentration(nProp) * dble(CurrentPartic%Geometry%Volume)
                                
                                if (AuxReal2 .GT. 1.) then
                                    Me%Monitor%InstLogMassByOrigin (Box, CurrentOrigin%ID) =     &
                                    Me%Monitor%InstLogMassByOrigin (Box, CurrentOrigin%ID) + log10(AuxReal2)
                                elseif (AuxReal2 .GT. 0. .and. AuxReal2 .LT.1) then
                                    SumNegativeLogMass(CurrentOrigin%ID) =                       &
                                    SumNegativeLogMassByOrigin(Box, CurrentOrigin%ID) + dabs(AuxReal2) 
                                end if

                                    Me%Monitor%NumberOfTracersFromOrigin(Box, CurrentOrigin%ID) =           &
                                    Me%Monitor%NumberOfTracersFromOrigin(Box, CurrentOrigin%ID) + 1 
                                    
                            
                                If ( CurrentPartic%Concentration(nProp) .GE. Me%Monitor%MonitorBox_TracerMinConc) then
                                    If (Me%Monitor%MassFractionType .EQ. Arithmetic) then
                                            SumContaminatedTracerConc(Box)=                                 &
                                            SumContaminatedTracerConc(Box) + CurrentPartic%Concentration(nProp)
                                        else
                                            SumContaminatedTracerConc(Box) =                                &
                                            SumContaminatedTracerConc(Box) +                                &
                                            log10(CurrentPartic%Concentration(nProp))                                           
                                    end if
                                    
                                    Me%Monitor%VolBoxContaminatedTracers(Box) =                             &
                                    Me%Monitor%VolBoxContaminatedTracers(Box) +                             &
                                    dble(CurrentPartic%Geometry%Volume)                              
                                    
                                    Me%Monitor%NbrBoxContaminatedTracers(Box) =                             &
                                    Me%Monitor%NbrBoxContaminatedTracers(Box) + 1
                                End If
                            
                            endif

                            CurrentPartic => CurrentPartic%Next
                        enddo                     
                        
                    endif


                    nProp = nProp + 1
                    CurrentProperty => CurrentProperty%Next
                enddo

                CurrentOrigin => CurrentOrigin%Next
            enddo

            !Calculates Box Mass and Box Mass Fraction By Origin 

                Me%Monitor%NumberOfTracers             (:)   = 0
                Me%Monitor%InstBoxMass                 (:)   = 0.
                Me%Monitor%InstBoxLogMass              (:)   = 0.
                Me%Monitor%InstBoxConc                 (:)   = 0.
                Me%Monitor%InstBoxMassFractionByOrigin (:,:) = 0.
                SumNegativeLogMass                     (:)   = 0.

            do Box = 1, NumberOfBoxes
                      
                CurrentOrigin => Me%FirstOrigin
                do while (associated(CurrentOrigin))

                    Me%Monitor%InstBoxLogMass(Box) =       &
                    Me%Monitor%InstBoxLogMass(Box) +       &
                    Me%Monitor%InstLogMassByOrigin (Box, CurrentOrigin%ID)

                    Me%Monitor%InstBoxMass(Box) =          &
                    Me%Monitor%InstBoxMass(Box) +          &
                    Me%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID)

                    Me%Monitor%NumberOfTracers(Box) =       &
                    Me%Monitor%NumberOfTracers(Box) +       &
                    Me%Monitor%NumberOfTracersFromOrigin (Box, CurrentOrigin%ID)

                    SumNegativeLogMass(Box) =               &
                    SumNegativeLogMass(Box) +               &
                    SumNegativeLogMassByOrigin(Box, CurrentOrigin%ID)
                    
                   CurrentOrigin => CurrentOrigin%Next
                enddo

                ! Computes the Fraction of mass in a box from each origin
                CurrentOrigin => Me%FirstOrigin
                do while (associated(CurrentOrigin))

                    if (Me%Monitor%MassFractionType .EQ. Arithmetic .and. Me%Monitor%InstBoxMass(Box) .GT. 0) then
                        Me%Monitor%InstBoxMassFractionByOrigin (Box, CurrentOrigin%ID) =     &
                        Me%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID) / Me%Monitor%InstBoxMass(Box)                 
                    elseif (Me%Monitor%MassFractionType .EQ. Geometric .and. Me%Monitor%InstBoxLogMass(Box) .GT. 0) then
                        Me%Monitor%InstBoxMassFractionByOrigin (Box, CurrentOrigin%ID) =     &
                        Me%Monitor%InstLogMassByOrigin (Box, CurrentOrigin%ID) / Me%Monitor%InstBoxLogMass(Box)                   
                    end if
 
                   CurrentOrigin => CurrentOrigin%Next
                enddo

                if (Me%Monitor%MassFractionType .EQ. Arithmetic) then
                    Me%Monitor%InstBoxConc (Box) = Me%Monitor%InstBoxMass(Box) / Me%Monitor%InstBoxVolume (Box)

                    Me%Monitor%AverageBoxContaminatedConc (Box) =                            &
                    SumContaminatedTracerConc(Box) / Me%Monitor%VolBoxContaminatedTracers(Box)
                    
                elseif (Me%Monitor%MassFractionType .EQ. Geometric) then
                    if (Me%Monitor%NumberOfTracers(Box) .GT. 0) then
                        Me%Monitor%InstBoxConc (Box) =                                      &
                        10**((Me%Monitor%InstBoxLogMass(Box) - SumNegativeLogMass(Box)) / dble(Me%Monitor%NumberOfTracers(Box))) &
                                                       * dble(Me%Monitor%NumberOfTracers(Box)) / Me%Monitor%InstBoxVolume (Box)
                    end if
                    if (Me%Monitor%NbrBoxContaminatedTracers(Box) .GT. 0) then
                        Me%Monitor%AverageBoxContaminatedConc (Box) =                       &
                        10**(SumContaminatedTracerConc(Box)/dble(Me%Monitor%NbrBoxContaminatedTracers(Box)))     
                    end if                                             
                end if
                If (Me%Monitor%ContaminationDepth .GT. 0.) then
                    Me%Monitor%ContaminationProbability(Box) = Me%Monitor%VolBoxContaminatedTracers(Box) /  &
                                                               Me%Monitor%SurfaceBoxVolume (Box)
                Else
                    Me%Monitor%ContaminationProbability(Box) = Me%Monitor%VolBoxContaminatedTracers(Box) /  &
                                                               Me%Monitor%InstBoxVolume (Box)                    
                End If
            enddo


        End If MMass


        !Integrates the values of the volume contributing from a given origin to the volume
        do Box = 1, NumberOfBoxes
            CurrentOrigin => Me%FirstOrigin
            do while (associated(CurrentOrigin))

                Me%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID) =       &
                    Me%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID) +   &
                    Me%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID)

                CurrentOrigin => CurrentOrigin%Next
            enddo
        enddo


        If (Me%State%MonitorPropMass) then
            allocate (AuxReal (2 + 4*(NumberOfOrigins + 1)))
        Else
            allocate (AuxReal (2*(NumberOfOrigins + 1)))
        End If

        !Writes the Time Serie
        do Box = 1, NumberOfBoxes

            !Instant Volume Values
            AuxReal(1) = Me%Monitor%InstBoxVolume (Box)
            CurrentOrigin => Me%FirstOrigin
            iO = 1
            do while (associated(CurrentOrigin))
                AuxReal(iO+1) = Me%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID)

                iO = iO + 1
                CurrentOrigin => CurrentOrigin%Next
            enddo


            !Integrated Values 
            AuxReal ((NumberOfOrigins + 1) + 1) = Me%Monitor%IntgBoxVolume (Box)
            CurrentOrigin => Me%FirstOrigin
            iO = 1
            do while (associated(CurrentOrigin))
                AuxReal ((NumberOfOrigins + 1) + 1 + iO) = Me%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID)
                iO = iO + 1
                CurrentOrigin => CurrentOrigin%Next
            enddo


            If (Me%State%MonitorPropMass) then

                !Instant Mass Values
                AuxReal(2 * (NumberOfOrigins + 1) + 1) = Me%Monitor%InstBoxMass (Box)
                CurrentOrigin => Me%FirstOrigin
                iO = 1
                do while (associated(CurrentOrigin))
                    AuxReal(2 * (NumberOfOrigins + 1) + 1 + iO) = Me%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID)

                    iO = iO + 1
                    CurrentOrigin => CurrentOrigin%Next
                enddo

                !Instant Geometric Concentration Values and Relative Mass values
                AuxReal(3 * (NumberOfOrigins + 1) + 1) = Me%Monitor%InstBoxConc (Box)
                AuxReal(3 * (NumberOfOrigins + 1) + 2) = Me%Monitor%AverageBoxContaminatedConc (Box)
                AuxReal(3 * (NumberOfOrigins + 1) + 3) = Me%Monitor%ContaminationProbability(Box)
                CurrentOrigin => Me%FirstOrigin
                iO = 1
                do while (associated(CurrentOrigin))
                    AuxReal(3 * (NumberOfOrigins + 1) + 3 + iO) = Me%Monitor%InstBoxMassFractionByOrigin (Box, CurrentOrigin%ID)

                    iO = iO + 1
                    CurrentOrigin => CurrentOrigin%Next
                enddo

            End If

            call WriteTimeSerieLine (Me%Monitor%ObjTimeSerie(Box), DataLine = AuxReal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'MonitorParticle - ModuleLagrangian - ERR02'
        enddo

        deallocate (AuxReal)


        !Unget The Boxes
        call UngetBoxDif(Me%ObjMonBox, MonitorBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MonitorParticle - ModuleLagrangian - ERR03'

    end subroutine MonitorParticle
    
    !--------------------------------------------------------------------------

    subroutine MonitorEulerian
        
        !Local-----------------------------------------------------------------
        real, dimension(:, :, :, :), pointer        :: GridConc
        real, dimension(:, :, :, :), pointer        :: GridMaxTracer
        real, dimension(:, :, :), pointer           :: AuxGrid3D
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: iProp, STAT_CALL, i, j, k
        type (T_Property), pointer                  :: CurrentProperty

        !Begin-----------------------------------------------------------------

        ILB     = Me%ExternalVar%Size%ILB   
        IUB     = Me%ExternalVar%Size%IUB   
        JLB     = Me%ExternalVar%Size%JLB   
        JUB     = Me%ExternalVar%Size%JUB   
        KLB     = Me%ExternalVar%Size%KLB   
        KUB     = Me%ExternalVar%Size%KUB   


        !Allocates auxiliar variable
        allocate (AuxGrid3D (ILB:IUB, JLB:JUB, KLB:KUB))

        !Fills Grid concentration
        call FillGridConcentration  (1, GridConc, GridMaxTracer = GridMaxTracer)

        iProp = 1
        CurrentProperty => Me%FirstOrigin%FirstProperty
        do while (associated(CurrentProperty))
            
            AuxGrid3D(:,:,:) = GridConc(:, :, :, iProp)

            do K = KLB, KUB
            do J = JLB, JUB
            do I = ILB, IUB
                if(Me%ExternalVar%WaterPoints3D(i,j,k) == WaterPoint)then
                    Me%EulerianMonitor%Mass(i,j,k) = AuxGrid3D(i, j, k) * &
                                                     Me%ExternalVar%VolumeZ(i, j, k)
                end if
            end do
            end do
            end do

            call BoxDif(Me%ObjEulMonBox,                                &
                        Me%EulerianMonitor%Mass,                        &
                        "Lagrangian_"//trim(CurrentProperty%Name),      &
                        Me%ExternalVar%WaterPoints3D,                   &
                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'MonitorEulerian - ModuleLagrangian - ERR01'

            CurrentProperty => CurrentProperty%Next
            iProp = iProp + 1
        enddo

        !Deallocates auxiliar variable
        deallocate (AuxGrid3D, GridConc)

        if (Me%OutPut%ConcMaxTracer) deallocate (GridMaxTracer)




    end subroutine MonitorEulerian

    !--------------------------------------------------------------------------

    subroutine ModifyParticStatistic ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        real, dimension(:, :, :, :), pointer     :: GridConc
        real, dimension(:, :, :, :), pointer        :: GridMaxTracer
        type (T_Property), pointer                  :: CurrentProperty
        real, dimension(:, :, :), pointer           :: AuxGrid3D 
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: iProp, STAT_CALL

        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        KLB    = Me%ExternalVar%Size%KLB
        KUB    = Me%ExternalVar%Size%KUB

        !Allocates auxiliar variable
        allocate (AuxGrid3D (ILB:IUB, JLB:JUB, KLB:KUB))

        !Fills Grid concentration
        call FillGridConcentration  (1, GridConc, GridMaxTracer = GridMaxTracer)

        iProp = 1
        CurrentProperty => Me%FirstOrigin%FirstProperty
        do while (associated(CurrentProperty))

            if (CurrentProperty%Statistics) then

                AuxGrid3D(:,:,:) = GridConc(:, :, :, iProp)
                call ModifyStatistic    (CurrentProperty%Statistic1_ID,                  &
                                         Value3D       = AuxGrid3D,                      &
                                         WaterPoints3D = Me%ExternalVar%WaterPoints3D,   &
                                         STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                               &
                    stop 'ModifyParticStatistic - ModuleLagrangian - ERR01'

                if (Me%OutPut%ConcMaxTracer) then

                    AuxGrid3D(:,:,:) = GridMaxTracer(:, :, :, iProp)

                    call ModifyStatistic (CurrentProperty%Statistic2_ID,                 &
                                          Value3D       = AuxGrid3D,                     &
                                          WaterPoints3D = Me%ExternalVar%WaterPoints3D,  &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                           &
                        stop 'ModifyParticStatistic - ModuleLagrangian - ERR02'
                endif 

                if (CurrentProperty%StatisticsLag) then

                    call ComputeStatisticsLag(CurrentProperty%ID, iProp)

                endif

            endif
            CurrentProperty => CurrentProperty%Next
            iProp = iProp + 1
        enddo

        !Deallocates auxiliar variable
        deallocate (AuxGrid3D   )
        deallocate (GridConc     )
        if (Me%OutPut%ConcMaxTracer) deallocate (GridMaxTracer)


    end subroutine ModifyParticStatistic
    !--------------------------------------------------------------------------

    subroutine ComputeStatisticsLag(PropertyID, iProp) 
        !Arguments-------------------------------------------------------------
        integer                                     :: PropertyID, iProp

        !Local-----------------------------------------------------------------
        type (T_Origin  ), pointer                  :: CurrentOrigin
        type (T_Property), pointer                  :: CurrentProperty
        type (T_Partic  ), pointer                  :: CurrentPartic
        real, dimension(:, :, :   ), allocatable    :: VolTracers, MassTracers
        real, dimension(:, :, :, :), allocatable    :: PropTimeStep
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        type (T_Position)                           :: DummYPos    
        integer                                     :: iClass, STAT_CALL   
        real                                        :: DT, Aux, PropSpatial, Concentration

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB
        KLB = Me%ExternalVar%WorkSize%KLB
        KUB = Me%ExternalVar%WorkSize%KUB

        call Search_Property(Me%FirstOrigin, PropertyID, CurrentProperty)

        allocate (VolTracers  (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate (MassTracers (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate (PropTimeStep(ILB:IUB, JLB:JUB, KLB:KUB, 1 : CurrentProperty%nClassesLag))

        VolTracers  (:,:,:  ) = 0.
        MassTracers (:,:,:  ) = 0.
        PropTimeStep(:,:,:,:) = 0.

        !Gets global time step 
        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStatisticsLag - ModuleLagrangian - ERR01'


        call GetStatisticClasses(CurrentProperty%Statistic1_ID,                          &
                                 CurrentProperty%ClassesLag, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ComputeStatisticsLag - ModuleLagrangian - ERR02'


        CurrentOrigin => Me%FirstOrigin

        do while (associated(CurrentOrigin))

            !Compute the total volume tracer by cell 
            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                i  = CurrentPartic%Position%I
                j  = CurrentPartic%Position%J
                k  = CurrentPartic%Position%K
       
WP1:            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then

                    VolTracers (i, j, k) = VolTracers (i, j, k) + CurrentPartic%Geometry%Volume
                    MassTracers(i, j, k) = MassTracers(i, j, k) + CurrentPartic%Mass(iProp)

                endif WP1

                CurrentPartic => CurrentPartic%Next

            enddo


            !Compute the frequency by cell 
            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                i  = CurrentPartic%Position%I
                j  = CurrentPartic%Position%J
                k  = CurrentPartic%Position%K
       
WP2:            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then


doClass:            do iClass = 1, CurrentProperty%nClassesLag
                        if (CurrentPartic%Concentration(iProp) >= CurrentProperty%ClassesLag(iClass, 1) .and. &
                            CurrentPartic%Concentration(iProp)  < CurrentProperty%ClassesLag(iClass, 2)) then

                            PropSpatial = CurrentPartic%Geometry%Volume / Me%ExternalVar%VolumeZ(i, j, k)


                        else
                            PropSpatial = 0
                        endif
                    
                        PropTimeStep(i, j, k, iClass) =  PropTimeStep(i, j, k, iClass) + PropSpatial

                    enddo doClass
                endif WP2

                CurrentPartic => CurrentPartic%Next

            enddo

            CurrentOrigin => CurrentOrigin%Next

        enddo

        do k= KLB, KUB
        do j= JLB, JUB
        do i= ILB, IUB

WP3:        if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then
                    
VT:             if (VolTracers(i, j, k) < Me%ExternalVar%VolumeZ(i, j, k)) then

                    Aux = (Me%ExternalVar%VolumeZ(i, j, k) - VolTracers(i, j, k)) / Me%ExternalVar%VolumeZ(i, j, k)

                    DummyPos%I = i
                    DummyPos%J = j
                    DummyPos%K = k
               
                    call GetAmbientConcentration (CurrentProperty,                       &
                                                  DummyPos,                              &
                                                  Concentration)


                else VT

                    Concentration = MassTracers(i, j, k) / VolTracers(i, j, k)
                    PropTimeStep(i, j, k, :)      = 0
                    Aux                           = 1

                endif VT

               
Class1:         do iClass = 1, CurrentProperty%nClassesLag
                    if (Concentration >= CurrentProperty%ClassesLag(iClass, 1) .and.     &
                        Concentration <  CurrentProperty%ClassesLag(iClass, 2)) then
                        PropSpatial = Aux

                    else
                        PropSpatial = 0
                    endif

                    PropTimeStep(i, j, k, iClass) =  PropTimeStep(i, j, k, iClass) + PropSpatial
        
                    CurrentProperty%FrequencyLag     (i, j, k, iClass) =                 &
                        (CurrentProperty%FrequencyLag(i, j, k, iClass) *                 &
                        (Me%ExternalVar%RunPeriod - DT) +                                &
                         PropTimeStep(i, j, k, iClass)  * DT         ) /                 &
                        (Me%ExternalVar%RunPeriod)

                enddo Class1

                if (sum(CurrentProperty%FrequencyLag(i, j, k, :)) < 0.9) then

                    write(*,*) 'ComputeStatisticsLag - ModuleLagrangian - WARN10'
                    write(*,*) i, j, k, sum(CurrentProperty%FrequencyLag(i, j, k, :)),'<90%'

                endif

                if (sum(CurrentProperty%FrequencyLag(i, j, k, :)) > 1.1) then

                    write(*,*) 'ComputeStatisticsLag - ModuleLagrangian - WARN20'
                    write(*,*) i, j, k, sum(CurrentProperty%FrequencyLag(i, j, k, :)),'>110%'

                endif


            endif WP3

        enddo
        enddo
        enddo


        call UnGetStatistic(CurrentProperty%Statistic1_ID,                               &
                            CurrentProperty%ClassesLag, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeStatisticsLag - ModuleLagrangian - ERR05'


        deallocate (VolTracers )
        deallocate (MassTracers)
        deallocate (PropTimeStep)
        

    end subroutine ComputeStatisticsLag

    !--------------------------------------------------------------------------


    subroutine Convert_XY_CellIJ (Position)

        !Arguments-------------------------------------------------------------
   
        type (T_Position)                           :: Position

        !Local-----------------------------------------------------------------
        integer                                     :: Ipos, Jpos, i, j
        real                                        :: Perc
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: ILower, IUpper, IMiddle
        integer                                     :: JLower, JUpper, JMiddle


        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB

        select case (Me%Grid%GridType)
        
        case (Grid1D)


            JLower = WS_JLB
            JUpper = WS_JUB
            JPos   = JUpper - JLower
            do while (Jpos > 1)
                JMiddle = int((JLower + JUpper)/2)
                if (Position%X >= Me%Grid%ParticX(JMiddle)) then
                    JLower = JMiddle
                else 
                    JUpper = JMiddle
                endif
                JPos = JUpper - JLower
            enddo
            JPos = JLower

            ILower = WS_ILB
            IUpper = WS_IUB
            IPos   = IUpper - ILower
            do while (Ipos > 1)
                IMiddle = int((ILower + IUpper)/2)
                if (Position%Y >= Me%Grid%ParticY(IMiddle)) then
                    ILower = IMiddle
                else 
                    IUpper = IMiddle
                endif
                IPos = IUpper - ILower
            enddo
            IPos = ILower

!doj:        do j = WS_JLB, WS_JUB
!                if (Position%X >= Me%Grid%ParticX(j  ) .and.                  &
!                    Position%X <  Me%Grid%ParticX(j+1)) then
!                    Jpos = j
!                    exit doj
!                endif
!            enddo doj

!doi:        do i = WS_ILB, WS_IUB
!                if (Position%Y >= Me%Grid%ParticY(i  ) .and.                  &
!                    Position%Y <  Me%Grid%ParticY(i+1)) then
!                    Ipos = i
!                    exit doi
!                endif
!            enddo doi

            !CellJ Position
            Perc = (Position%X - Me%Grid%ParticX(Jpos)) /                     &
                   (Me%Grid%ParticX(Jpos+1) - Me%Grid%ParticX(Jpos))
                                                                                       
            Position%CellJ = (Jpos - 1) + Perc                         

            !CellI Position
            Perc = (Position%Y - Me%Grid%ParticY(Ipos)) /                     &
                   (Me%Grid%ParticY(Ipos+1) - Me%Grid%ParticY(Ipos))
                                                                                       
            Position%CellI = (Ipos - 1) + Perc

        case (Grid2D)
        
            Ipos = 1
            Jpos = 1

do1:        do j = WS_JLB, WS_JUB
do2:        do i = WS_ILB, WS_IUB

                if (Position%X >= Me%Grid%ParticXX(i, j  ) .and.              &
                    Position%X <  Me%Grid%ParticXX(i, j+1) .and.              &
                    Position%Y >= Me%Grid%ParticYY(i, j  ) .and.              &
                    Position%Y <  Me%Grid%ParticYY(i+1, j)) then
                    Jpos = j
                    Ipos = i
                    exit do1
                endif   

            enddo do2
            enddo do1


            !CellJ Position
            Perc = (Position%X - Me%Grid%ParticXX(Ipos, Jpos)) /                  &
                   (Me%Grid%ParticXX(Ipos, Jpos+1) - Me%Grid%ParticXX(Ipos, Jpos))
                                                                                       
            Position%CellJ = (Jpos - 1) + Perc                         

            !CellI Position
            Perc = (Position%Y - Me%Grid%ParticYY(Ipos, Jpos)) /                  &
                   (Me%Grid%ParticYY(Ipos+1, Jpos) - Me%Grid%ParticYY(Ipos, Jpos))
                                                                                       
            Position%CellI = (Ipos - 1) + Perc

        end select


    end subroutine Convert_XY_CellIJ

    !--------------------------------------------------------------------------

    subroutine Convert_CellIJ_XY( Position)

        !Arguments-------------------------------------------------------------
   
        type (T_Position)                           :: Position

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        real                                        :: Perc

        !----------------------------------------------------------------------

        !Cell Position
        i           = int(Position%CellI) + 1
        j           = int(Position%CellJ) + 1

        select case (Me%Grid%GridType)
        
        case (Grid1D)

            !Distance
            Perc       = (Position%CellJ - int(Position%CellJ)) *                            &
                         (Me%Grid%ParticX(j+1) - Me%Grid%ParticX(j))
            Position%X = Me%Grid%ParticX(j) + Perc

            !Distance
            Perc       = (Position%CellI - int(Position%CellI)) *                            &
                         (Me%Grid%ParticY(i+1) - Me%Grid%ParticY(i))
            Position%Y = Me%Grid%ParticY(i) + Perc

        case (Grid2D)

            !Distance
            Perc       = (Position%CellJ - int(Position%CellJ)) *                            &
                         (Me%Grid%ParticXX(i, j+1) - Me%Grid%ParticXX(i, j))
            Position%X = Me%Grid%ParticXX(i, j) + Perc

            !Distance
            Perc       = (Position%CellI - int(Position%CellI)) *                            &
                         (Me%Grid%ParticYY(i+1, j) - Me%Grid%ParticYY(i, j))
            Position%Y = Me%Grid%ParticYY(i, j) + Perc

        end select


    end subroutine Convert_CellIJ_XY

    !--------------------------------------------------------------------------

    subroutine Convert_Coord1D_CellIJ (Position, XX, YY)

        !Arguments-------------------------------------------------------------
   
        type (T_Position)                           :: Position
        real,   dimension(:), pointer               :: XX, YY

        !Local-----------------------------------------------------------------
        integer                                     :: Ipos, Jpos
        real                                        :: Perc
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: ILower, IUpper, IMiddle
        integer                                     :: JLower, JUpper, JMiddle


        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB


        JLower = WS_JLB
        JUpper = WS_JUB
        JPos   = JUpper - JLower
        do while (Jpos > 1)
            JMiddle = int((JLower + JUpper)/2)
            if (Position%X >= XX(JMiddle)) then
                JLower = JMiddle
            else 
                JUpper = JMiddle
            endif
            JPos = JUpper - JLower
        enddo
        JPos = JLower

        ILower = WS_ILB
        IUpper = WS_IUB
        IPos   = IUpper - ILower
        do while (Ipos > 1)
            IMiddle = int((ILower + IUpper)/2)
            if (Position%Y >= YY(IMiddle)) then
                ILower = IMiddle
            else 
                IUpper = IMiddle
            endif
            IPos = IUpper - ILower
        enddo
        IPos = ILower


        !CellJ Position
        Perc = (Position%X - XX(Jpos)) /                     &
               (XX(Jpos+1) - XX(Jpos))
                                                                                   
        Position%CellJ = (Jpos - 1) + Perc                         

        !CellI Position
        Perc = (Position%Y - YY(Ipos)) /                     &
               (YY(Ipos+1) - YY(Ipos))
                                                                                   
        Position%CellI = (Ipos - 1) + Perc



    end subroutine Convert_Coord1D_CellIJ

    !--------------------------------------------------------------------------


    subroutine Convert_CellIJ_IJ(Position)

        !Arguments-------------------------------------------------------------
        type (T_Position)                           :: Position

        !----------------------------------------------------------------------

        Position%I = int(Position%CellI) + 1
        Position%J = int(Position%CellJ) + 1

        !----------------------------------------------------------------------

    end subroutine Convert_CellIJ_IJ

    !--------------------------------------------------------------------------

    subroutine Convert_Z_CellK (CurrentOrigin, Position, PositionCorrected)

        !Arguments-------------------------------------------------------------
   
        type (T_Origin),     pointer                :: CurrentOrigin
        type (T_Position)                           :: Position
        logical, optional                           :: PositionCorrected


        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, kFloor
        real                                        :: Perc
        logical                                     :: PositionCorrected_

        !----------------------------------------------------------------------

        i       = Position%I
        j       = Position%J
        k       = Me%ExternalVar%WorkSize%KUB
        kFloor  = Me%ExternalVar%kFloor(i, j)
        
        PositionCorrected_ = .false.

        !If the Particle is located above the surface, its placed closed to the
        !surface
        if (Position%Z < Me%ExternalVar%SZZ(i, j, k)) then
            Position%Z = Me%ExternalVar%SZZ(i, j, k)
            PositionCorrected_ = .true.
        endif

        !If the Particle is located below the bottom, its placed closed to the
        !bottom
        if (Position%Z >= Me%ExternalVar%SZZ(i, j, kFloor-1)) then

            if (CurrentOrigin%State%Deposition) then

                !If the origin emit particles that can be deposited in the bottom and 
                !the particle tries to cross the bottom then is placed at distance from 
                !the bottom lower than the distance beyond each the particle is test 
                !if the shear stress is low enough to consider the particle deposited. 
                Position%Z = Me%ExternalVar%SZZ(I, J, KFloor-1) -  &
                             CurrentOrigin%Deposition%BottomDistance / 2.

            else

                Position%Z = Me%ExternalVar%SZZ(I, J, KFloor) +  &
                             Me%ExternalVar%DWZ(I, J, KFloor) * 0.9

            endif

            PositionCorrected_ = .true.

        endif

        if (present(PositionCorrected)) PositionCorrected = PositionCorrected_


        !Checks the layer
        do while(Me%ExternalVar%SZZ(i, j, k) <= Position%Z)
            k = k - 1
        end do

        Perc = (Me%ExternalVar%SZZ(i, j, k) - Position%Z) /                   &
               (Me%ExternalVar%SZZ(i, j, k) -                                 &
                Me%ExternalVar%SZZ(i, j, k+1))

        !Avoid Rounding erros
        if (Perc >= 0.999) Perc = 0.99000

        Position%CellK = k + Perc

        !----------------------------------------------------------------------

    end subroutine Convert_Z_CellK

    !--------------------------------------------------------------------------

    subroutine Convert_CellK_Z (Position)

        !Arguments-------------------------------------------------------------
   
        type (T_Position)                           :: Position


        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        real                                        :: Perc

        !----------------------------------------------------------------------

        i = Position%I
        j = Position%J
        k = int (Position%CellK) + 1

        Perc       = Position%CellK - int(Position%CellK)

        Position%Z = Me%ExternalVar%SZZ(i, j, k-1) -                          &
                     Perc * Me%ExternalVar%DWZ(i, j, k) 
        
    end subroutine Convert_CellK_Z

    !--------------------------------------------------------------------------

    subroutine Convert_CellK_K (Position)

        !Arguments-------------------------------------------------------------
        type (T_Position)                           :: Position

        Position%K = int(Position%CellK) + 1

    end subroutine Convert_CellK_K

    !--------------------------------------------------------------------------

    subroutine Search_Property(CurrentOrigin, PropertyID, CurrentProperty)

        !Arguments-------------------------------------------------------------
        type(T_Origin),             pointer             :: CurrentOrigin
        type(T_Property),           pointer             :: CurrentProperty
        integer         ,           intent (IN)         :: PropertyID
        !----------------------------------------------------------------------


        CurrentProperty => CurrentOrigin%FirstProperty

        do while (associated(CurrentProperty)) 
            if (CurrentProperty%ID ==PropertyID) then
                exit        
            else
                CurrentProperty => CurrentProperty%Next                 
            end if    
        end do    

        if (.not. associated(CurrentProperty)) then

            call SetError(FATAL_, INTERNAL_, 'Search_Property - ModuleLagrangian - ERR01')

        end if 

 

    end subroutine Search_Property

    !--------------------------------------------------------------------------


    subroutine ParticleOutput ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        real(8), dimension(:), pointer              :: Matrix1D, Matrix1DX, Matrix1DY
        integer                                     :: nP
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: OutPutNumber
        type (T_Time)                               :: Actual
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second
        integer, dimension (6)                      :: TimeAux
        integer                                     :: STAT_CALL
        integer, dimension(:), allocatable          :: TotParticle
        real                                        :: MaximumDepth
        integer                                     :: ParticleBeached
        integer                                     :: iProp, nProp, iGroup, nGroupProp
        type (T_Property), pointer                  :: CurrentProperty, FirstProperty
        type (T_Property), pointer                  :: AgeProperty
        integer                                     :: nPropAge, n
        integer, dimension(:, :, :   ), pointer     :: GridTracerNumber
        real, dimension(:, :, :, :), pointer     :: GridConc
        real,    dimension(:, :, :, :), pointer     :: GridMaxTracer
        real,    dimension(:, :, :   ), pointer     :: GridBottomConc
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        integer, dimension(:, :, :), pointer        :: OpenPoints3D
        real, dimension(:, :, :), pointer           :: SZZ
        real, dimension(:, :   ), pointer           :: OutPutMatrix
        integer                                     :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        character(StringLength)                     :: AuxChar
        real                                        :: Xorig, Yorig, GridAngle, Xaux, Yaux
        integer                                     :: OutPutLines, JetTotalParticles, FirstParticle
        type (T_Position)                           :: Position
        real(8)                                     :: AverageX, AverageY, Stdv, RadiusOfInfluence
        


        !Gets Origin
        call GetGridOrigin  (Me%ObjHorizontalGrid, Xorig, Yorig)

        !Gets GridAngle
        call GetGridAngle   (Me%ObjHorizontalGrid, GridAngle)


        call GetMaximumValue(Me%ObjGridData, MaximumDepth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR10'


        !Shorten Var
        Actual          =  Me%ExternalVar%Now
        OpenPoints3D    => Me%ExternalVar%OpenPoints3D
        SZZ             => Me%ExternalVar%SZZ
        WorkILB         =  Me%ExternalVar%WorkSize%ILB 
        WorkIUB         =  Me%ExternalVar%WorkSize%IUB 
                        
        WorkJLB         =  Me%ExternalVar%WorkSize%JLB 
        WorkJUB         =  Me%ExternalVar%WorkSize%JUB 
                        
        WorkKLB         =  Me%ExternalVar%WorkSize%KLB 
        WorkKUB         =  Me%ExternalVar%WorkSize%KUB 

        if (Me%Output%Write_) then

            OutPutNumber = Me%OutPut%NextOutPut
TOut:       if (Actual >= Me%OutPut%OutTime(OutPutNumber)) then 

                nP = 0                
                CurrentOrigin => Me%FirstOrigin
d1:             do while (associated(CurrentOrigin))

                    !XPosition
                    CurrentPartic   => CurrentOrigin%FirstPartic
                    do while (associated(CurrentPartic))
                        nP = nP + 1
                        CurrentPartic => CurrentPartic%Next
                    enddo
                    CurrentOrigin => CurrentOrigin%Next
                enddo d1

i1:             if (nP>0) then


                    call ExtractDate(Actual, Year = Year, Month  = Month,  Day    = Day,     &
                                             Hour = Hour, Minute = Minute, Second = Second)

                    TimeAux(1) = int(Year  )
                    TimeAux(2) = int(Month )
                    TimeAux(3) = int(Day   )
                    TimeAux(4) = int(Hour  )
                    TimeAux(5) = int(Minute)
                    TimeAux(6) = int(Second)


                    !Writes the Instant - HDF 5
                    call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                        AuxTime(4), AuxTime(5), AuxTime(6))
                    TimePtr => AuxTime
                    call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR20'

                    call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                         Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR30'

                    !Writes SZZ - HDF 5
                    call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,   &
                                         WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR40'

                    call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",  &
                                         "m", Array3D = SZZ, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR50'


                    !Writes OpenPoints - HDF 5
                    call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,            &
                                         WorkJLB, WorkJUB, WorkKLB, WorkKUB,                 &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR60'

                    call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",  &
                                         "-", Array3D = OpenPoints3D, OutputNumber = OutPutNumber, &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR70'


                    allocate (TotParticle(Me%nGroups))
                    TotParticle(:) = 0

                    JetTotalParticles = 0

                    !Writes Data for every origin
                    CurrentOrigin => Me%FirstOrigin
    CurrOr:         do while (associated(CurrentOrigin))

                        allocate   (Matrix1DX(CurrentOrigin%nParticle))
                        allocate   (Matrix1DY(CurrentOrigin%nParticle))
                        Matrix1DX(:) = FillValueReal
                        Matrix1DY(:) = FillValueReal

                        call HDF5SetLimits  (Me%ObjHDF5, 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR80'

                        !XPosition
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            if (Me%Grid%GridType == Grid1D) then
                                Xaux          = CurrentPartic%Position%X
                                Yaux          = CurrentPartic%Position%Y
                                call RodaXY (Xorig, Yorig, GridAngle, Xaux, Yaux)
                                Matrix1DX(nP)  =  Xaux
                                Matrix1DY(nP)  =  Yaux
                            else
                                Matrix1DX(nP)  = CurrentPartic%Position%X
                                Matrix1DY(nP)  = CurrentPartic%Position%Y
                            endif
                            CurrentPartic => CurrentPartic%Next
                        enddo
                        if (nP > 0) then

                            if (CurrentOrigin%AveragePositionON) then
                                AverageX = sum(Matrix1DX(1:nP)) / real(nP)
                                AverageY = sum(Matrix1DY(1:nP)) / real(nP)
                                if (nP > 1) then
                                    Stdv = 0
                                    do n = 1, nP
                                        Stdv = Stdv + ((Matrix1DX(n) - AverageX)**2. +      &
                                                       (Matrix1DY(n) - AverageY)**2.) / (real(nP) - 1)
                                
                                    enddo
                                    RadiusOfInfluence = CurrentOrigin%CoefRadius * sqrt(Stdv)
                                else
                                    RadiusOfInfluence = 0.
                                endif
                            else
                                AverageX          = FillValueReal
                                AverageY          = FillValueReal
                                RadiusOfInfluence = FillValueReal
                            endif 

                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/X Pos", &
                                                "X Position",  "m", Array1D = Matrix1DX,                      &
                                                 Average = AverageX, Radius = RadiusOfInfluence,              &
                                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR90'
                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Y Pos", &
                                                "Y Position",  "m", Array1D = Matrix1DY,                      &
                                                 Average = AverageY, Radius = RadiusOfInfluence,              &
                                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR100'

                            if (Me%OutPut%OriginEnvelope) then
                                call WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY, &
                                                          "X Pos", "Y Pos",  "m", OutputNumber)  
                            endif

                        endif

                      !GeoGraphic Position
                        if (Me%Grid%GridType == Grid2D) then

                            !Longitude, Latitude
                            CurrentPartic   => CurrentOrigin%FirstPartic
                            nP = 0
                            do while (associated(CurrentPartic))
                                nP = nP + 1
                                Matrix1DX(nP)  =  GeographicCoordinates (CurrentPartic%Position, 1)
                                Matrix1DY(nP)  =  GeographicCoordinates (CurrentPartic%Position, 2)
                                CurrentPartic => CurrentPartic%Next
                            enddo

                            if (nP > 0) then

                                if (CurrentOrigin%AveragePositionON) then
                                    AverageX = sum(Matrix1DX(1:nP)) / real(nP)
                                    AverageY = sum(Matrix1DY(1:nP)) / real(nP)

                                    if (nP > 1) then
                                        Stdv = 0.
                                        do n = 1, nP
                                            Stdv = Stdv + ((Matrix1DX(n) - AverageX)**2. +      &
                                                           (Matrix1DY(n) - AverageY)**2.) / (real(nP) - 1)
                                
                                        enddo
                                        RadiusOfInfluence = CurrentOrigin%CoefRadius * sqrt(Stdv)
                                    else
                                        RadiusOfInfluence = 0.
                                    endif

                                else
                                    AverageX          = FillValueReal
                                    AverageY          = FillValueReal
                                    RadiusOfInfluence = FillValueReal
                                endif 

                                call HDF5SetLimits  (Me%ObjHDF5, 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR101'


                                !HDF 5
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Longitude", &
                                                    "Longitude",  "º", Array1D = Matrix1DX,                      &
                                                     Average = AverageX, Radius = RadiusOfInfluence,             &
                                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR110'

                                !HDF 5
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Latitude", &
                                                    "Latitude",  "º", Array1D = Matrix1DY,                       &
                                                     Average = AverageY, Radius = RadiusOfInfluence,             &
                                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR120'

                                if (Me%OutPut%OriginEnvelope) then
                                    call WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY, &
                                                              "Longitude", "Latitude", "º", OutputNumber)  
                                endif


                            endif

                
                        endif
 
                        deallocate   (Matrix1DX)
                        deallocate   (Matrix1DY)


                        allocate    (Matrix1D(CurrentOrigin%nParticle))
                        Matrix1D(:) = FillValueReal

                        call HDF5SetLimits  (Me%ObjHDF5, 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR121'


                        !Real ZPosition
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            Matrix1D(nP)  =  CurrentPartic%Position%Z
                            CurrentPartic => CurrentPartic%Next
                        enddo            
                        if (nP > 0) then
                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Z Pos", &
                                                "Z Position",  "m", Array1D = Matrix1D, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR130'
                        endif

                        !ZPosition for Vertical Cut
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            Matrix1D(nP)  =  MaximumDepth - CurrentPartic%Position%Z
                            CurrentPartic => CurrentPartic%Next
                        enddo            

                        !Volume
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            Matrix1D(nP)  =  CurrentPartic%Geometry%Volume
                            CurrentPartic => CurrentPartic%Next
                        enddo            
                        if (nP > 0) then
                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Volume", &
                                                "Volume",  "m3", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR140'
                        endif

                        !OriginNumber
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            Matrix1D(nP)  =  CurrentOrigin%ID
                            CurrentPartic => CurrentPartic%Next
                        enddo            
                        if (nP > 0) then
                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Origin ID", &
                                                "Origin ID",  "-", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR150'
                        endif

                    
                        if (CurrentOrigin%State%Age) then

                            !Age
                            CurrentPartic   => CurrentOrigin%FirstPartic
                            nP = 0
                            do while (associated(CurrentPartic))
                                nP = nP + 1
                                Matrix1D(nP)  =  CurrentPartic%Age / 86400.
                                CurrentPartic => CurrentPartic%Next
                            enddo            
                            if (nP > 0) then
                                !HDF 5
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Age", &
                                                    "Age",  "days", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR160'
                            endif

                        endif


                        !Oil-Beached Particles
                        if (Me%State%AssociateBeachProb .and. CurrentOrigin%Beaching) then

                            CurrentPartic   => CurrentOrigin%FirstPartic
                            nP = 0
                            do while (associated(CurrentPartic))
                                nP = nP + 1
                                if (.NOT. CurrentPartic%Beached) then 
                                    ParticleBeached = 1
                                else
                                    ParticleBeached = 2
                                end if

                                Matrix1D(nP)  = ParticleBeached
                                CurrentPartic => CurrentPartic%Next
                            enddo            
                            if (nP > 0) then
                                !HDF 5
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Beached",       &
                                                    "Beached",  "[-]", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR170'
                            endif
                    
                        endif

                       if (CurrentOrigin%State%Deposition) then

                            CurrentPartic   => CurrentOrigin%FirstPartic
                            nP = 0
                            do while (associated(CurrentPartic))
                                nP = nP + 1
                                if (CurrentPartic%Deposited) then 

                                    Matrix1D(nP)  = 1 

                                else

                                    Matrix1D(nP)  = 2

                                end if


                                CurrentPartic => CurrentPartic%Next
                            enddo            
                            if (nP > 0) then
                                !HDF 5
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Deposition State", &
                                                    "State",  "ON/OFF", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR180'
                            endif
                    
                        endif


                        CurrentProperty => CurrentOrigin%FirstProperty
                        do while (associated(CurrentProperty))

                            if (CurrentProperty%ID == Fecal_Coliforms_) then

                                CurrentPartic   => CurrentOrigin%FirstPartic
                                nP = 0
                                do while (associated(CurrentPartic))
                                    nP = nP + 1
                                    !T90 in hours
                                    Matrix1D(nP)  = CurrentPartic%T90 / 3600.

                                    CurrentPartic => CurrentPartic%Next
                                enddo            
                                if (nP > 0) then
                                    !HDF 5
                                    call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/T90-Fecal Coli.", &
                                                        "T90",  "hours", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                         STAT = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR190'
                                endif


                            endif

                            CurrentProperty => CurrentProperty%Next
                        enddo



                        !Properties
                        nProp =  1
                        CurrentProperty => CurrentOrigin%FirstProperty
                        do while (associated(CurrentProperty))

                            nullify(AgeProperty)
                            nPropAge = 1
                            if (CurrentProperty%ID == Larvae_) then
                                AgeProperty => CurrentOrigin%FirstProperty
                                do while (AgeProperty%ID /= Age_)
                                    nPropAge = nPropAge + 1
                                    AgeProperty => AgeProperty%Next
                                enddo
                            endif

                            CurrentPartic => CurrentOrigin%FirstPartic
                            nP = 0
                            do while (associated(CurrentPartic))
                                nP = nP + 1
                                Matrix1D(nP)  =  CurrentPartic%Concentration(nProp)
                                CurrentPartic => CurrentPartic%Next
                            enddo            
                            if (nP > 0) then
                                !HDF 5
                                call HDF5WriteData  (Me%ObjHDF5,                  &
                                                     "/Results/"//trim(CurrentOrigin%Name)// &
                                                     "/"//trim(CurrentProperty%Name),        &
                                                     trim(CurrentProperty%Name),             &
                                                     trim(CurrentProperty%Units),            &
                                                     Array1D = Matrix1D,                     &
                                                     OutputNumber = OutPutNumber,            &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR200'
                            endif
                            nProp = nProp + 1
                            CurrentProperty => CurrentProperty%Next
                        enddo

                        deallocate  (Matrix1D)

    iplume:             if (CurrentOrigin%State%ComputePlume) then

                            call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                            OutPutLines = OutPutLines, STAT=STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR210'

                            if (OutPutLines < 1) then

                                call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)

                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR215'

                                call ActualizeJetProperties(CurrentOrigin)

                                call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                                                     OutPutLines = OutPutLines, STAT=STAT_CALL)

                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR220'

                                if (OutPutLines < 1) then
                                    stop 'ParticleOutput - ModuleLagrangian - ERR230'
                                endif

                            endif

                            JetTotalParticles = JetTotalParticles + OutPutLines

                            allocate   (Matrix1DX(OutPutLines))
                            allocate   (Matrix1DY(OutPutLines))

                            Matrix1DX(:) = FillValueReal
                            Matrix1DY(:) = FillValueReal

                            call HDF5SetLimits  (Me%ObjHDF5, 1, OutPutLines, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR240'

                            if (Me%Grid%GridType == Grid1D) then
                                do np = 1, OutPutLines
                                    Xaux          = OutPutMatrix(np,2)
                                    Yaux          = OutPutMatrix(np,3)
                                    call RodaXY (Xorig, Yorig, GridAngle, Xaux, Yaux)
                                    Matrix1DX(np) = Xaux
                                    Matrix1DY(np) = Yaux
                                enddo
                            else
                                do np = 1, OutPutLines
                                    Position%x = OutPutMatrix(np,2)
                                    Position%y = OutPutMatrix(np,3)
                                    call Convert_XY_CellIJ(Position)
                                    call Convert_CellIJ_IJ(Position)
                                    Matrix1DX(np) = GeographicCoordinates (Position, 1)
                                    Matrix1DY(np) = GeographicCoordinates (Position, 2)
                                enddo
                            endif

                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Plume X", &
                                                "Plume X",  "m", Array1D = Matrix1DX, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR250'

                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Plume Y", &
                                                "Plume Y",  "m", Array1D = Matrix1DY, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR260'


                            Matrix1DX(1:OutPutLines) = OutPutMatrix(1:OutPutLines,4)
                            Matrix1DY(1:OutPutLines) = OutPutMatrix(1:OutPutLines,6)


                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Plume Z", &
                                                "Plume Z",  "m", Array1D = Matrix1DX, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR270'

                            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/Plume Conc", &
                                                "Plume Conc",  "a", Array1D = Matrix1DY, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR280'

                            call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR290'

                            deallocate   (Matrix1DX)
                            deallocate   (Matrix1DY)

                        endif iplume


                        !Integrates the sum of the particle which belong to the same 
                        !Group
                        do iGroup = 1, Me%nGroups
                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                TotParticle(iGroup) = TotParticle(iGroup) + CurrentOrigin%nParticle
                            endif
                        enddo


                        CurrentOrigin => CurrentOrigin%Next

                    enddo CurrOr

    !                if (Me%State%WQM .or. Me%State%FCF .or. Me%State%Partition) then
                        do iGroup = 1, Me%nGroups
                            call FillGridConcentration  (iGroup, GridConc,                   &
                                                         GridTracerNumber = GridTracerNumber,&
                                                         GridMaxTracer    = GridMaxTracer,   &
                                                         GridBottomConc   = GridBottomConc)
                            call WriteGridConcentration (iGroup, GridConc,                   &
                                                         GridTracerNumber = GridTracerNumber,&
                                                         GridMaxTracer    = GridMaxTracer,   &
                                                         GridBottomConc   = GridBottomConc)
                        enddo
    !                endif


                    if (Me%State%Oil) then
                        call WriteOilGridThickness 
                        call WriteOilGridConcentration ()
                    endif


                    !Writes 1D data for every group
    AllAtOnes:      if (Me%nOrigins > 1) then

                    do iGroup = 1, Me%nGroups
                
                        if (TotParticle(iGroup) > 1) then
                        allocate    (Matrix1D(TotParticle(iGroup)))
                        Matrix1D(:) = FillValueReal

                        write (AuxChar, fmt='(i3)') Me%GroupIDs(iGroup)

                        call HDF5SetLimits  (Me%ObjHDF5, 1, TotParticle(iGroup),  &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR300'


                        !(XPosition)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    XPos:               do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    if (Me%Grid%GridType == Grid1D) then
                                        Xaux          = CurrentPartic%Position%X
                                        Yaux          = CurrentPartic%Position%Y
                                        call RodaXY (Xorig, Yorig, GridAngle, Xaux, Yaux)
                                        Matrix1D(nP)  =  Xaux
                                    else
                                        Matrix1D(nP)  =  CurrentPartic%Position%X
                                    endif
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo XPos

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5,                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/X Pos",                       &
                                                   "X Position",                             &
                                                   "m",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)

                        !(YPosition)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    YPos:               do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    if (Me%Grid%GridType == Grid1D) then
                                        Xaux          = CurrentPartic%Position%X
                                        Yaux          = CurrentPartic%Position%Y
                                        call RodaXY (Xorig, Yorig, GridAngle, Xaux, Yaux)
                                        Matrix1D(nP)  =  Yaux
                                    else
                                        Matrix1D(nP)  =  CurrentPartic%Position%Y
                                    endif
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo YPos

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5,                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Y Pos",                       &
                                                   "Y Position",                             &
                                                   "m",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)


                        !GeoGraphic Position
                        if (Me%Grid%GridType == Grid2D) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
                            do while (associated(CurrentOrigin))

                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                    CurrentPartic   => CurrentOrigin%FirstPartic
                                    do while (associated(CurrentPartic))
                                        Matrix1D(nP)  =  GeographicCoordinates (CurrentPartic%Position, 1)
                                        CurrentPartic => CurrentPartic%Next
                                        nP = nP + 1
                                    enddo
                                endif
                                CurrentOrigin => CurrentOrigin%Next
                            enddo

                            !HDF 5
                            call HDF5WriteData        (Me%ObjHDF5,                    &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/Longitude",                   &
                                                       "Longitude",                              &
                                                       "º",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)



                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
                            do while (associated(CurrentOrigin))

                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                    CurrentPartic   => CurrentOrigin%FirstPartic
                                    do while (associated(CurrentPartic))
                                        Matrix1D(nP)  =  GeographicCoordinates (CurrentPartic%Position, 2)
                                        CurrentPartic => CurrentPartic%Next
                                        nP = nP + 1
                                    enddo
                                endif
                                CurrentOrigin => CurrentOrigin%Next
                            enddo

                            !HDF 5
                            call HDF5WriteData        (Me%ObjHDF5,                    &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/Latitude",                    &
                                                       "Latitude",                               &
                                                       "º",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)
                        endif


                        !(ZPosition)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    ZPos:               do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    Matrix1D(nP)  = CurrentPartic%Position%Z
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo
                            endif
                                        
                            CurrentOrigin => CurrentOrigin%Next
                        enddo ZPos

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5,                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Z Pos",                       &
                                                   "Z Position",                             &
                                                   "m",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)


                        !(ZPosition for Vertical cut)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    ZPosInv:            do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    Matrix1D(nP)  = MaximumDepth - CurrentPartic%Position%Z
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo            
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo ZPosInv


                        !(Volume)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    Volume:             do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    Matrix1D(nP)  = CurrentPartic%Geometry%Volume
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo            
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo Volume

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5,                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Volume",                      &
                                                   "Volume",                                 &
                                                   "m3",                                     &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)

                        !(OriginNumber)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    Origin:             do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    Matrix1D(nP)  = CurrentOrigin%ID
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo            
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo Origin


                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5,                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Origin ID",                   &
                                                   "Origin ID",                              &
                                                   "-",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)



                        !(Oil-Beached Particles)
                        if (Me%State%AssociateBeachProb) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
    OilState:               do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
    IfBeaching:                     if (CurrentOrigin%Beaching) then
                                        CurrentPartic   => CurrentOrigin%FirstPartic
                                        do while (associated(CurrentPartic))

                                            if (.NOT. CurrentPartic%Beached) then 
                                                ParticleBeached = 1
                                            else
                                                ParticleBeached = 2
                                            end if

                                            Matrix1D(nP)  = ParticleBeached
                                            CurrentPartic => CurrentPartic%Next
                                            nP = nP + 1
                                        enddo
                                    end if IfBeaching
                               endif
                               CurrentOrigin => CurrentOrigin%Next
                            enddo OilState


                            !HDF 5
                            call HDF5WriteData        (Me%ObjHDF5,                               &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                        //"/Data_1D/Beached",                    &
                                                       "Beached",                                &
                                                       "-",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)


                        end if


                        !(Deposited Particles)
                        if (Me%State%Deposition) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
    DepState:               do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then

                                    if (CurrentOrigin%State%Deposition) then

                                        CurrentPartic   => CurrentOrigin%FirstPartic
                                        do while (associated(CurrentPartic))

                                            if (CurrentPartic%Deposited) then 
                                                Matrix1D(nP)  = 1
                                            else
                                                Matrix1D(nP) = 2
                                            end if

                                            CurrentPartic => CurrentPartic%Next
                                            nP = nP + 1
                                        enddo
                                            
                                    endif

                                endif
                
                                CurrentOrigin => CurrentOrigin%Next

                            enddo DepState


                            !HDF 5
                            call HDF5WriteData        (Me%ObjHDF5,                    &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/Deposition State",            &
                                                       "Origin ID",                              &
                                                       "-",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)


                        end if

                        !T90
    iT90:               if (Me%State%T90Variable) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
    dT90:                   do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then

                                    CurrentProperty => CurrentOrigin%FirstProperty
                                    do while (associated(CurrentProperty))

                                        if (CurrentProperty%ID == Fecal_Coliforms_) then

                                            CurrentPartic   => CurrentOrigin%FirstPartic

                                            do while (associated(CurrentPartic))
                                                !T90 in hours
                                                Matrix1D(nP)  = CurrentPartic%T90 / 3600.

                                                CurrentPartic => CurrentPartic%Next
                                                nP = nP + 1
                                            enddo            
                                            exit
                                        endif

                                        CurrentProperty => CurrentProperty%Next
                                   
                                    enddo


                                endif
                            
                                CurrentOrigin => CurrentOrigin%Next

                            enddo dT90

                            !HDF 5
                            call HDF5WriteData        (Me%ObjHDF5,                               &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/T90-Fecal coli.",             &
                                                       "T90",                                    &
                                                       "-",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)


                    
                        endif iT90



                        !Age
                        if (Me%State%Age) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
    Age:                    do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then

                                    if (CurrentOrigin%State%Age) then

                                        CurrentPartic   => CurrentOrigin%FirstPartic
                                        do while (associated(CurrentPartic))

                                            Matrix1D(nP)  =  CurrentPartic%Age / 86400.
                                            CurrentPartic => CurrentPartic%Next
                                            nP = nP + 1
                                        enddo
                                            
                                    endif

                                endif
                
                                CurrentOrigin => CurrentOrigin%Next
                            enddo Age


                            !HDF 5
                            call HDF5WriteData        (Me%ObjHDF5,                    &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/Age",                         &
                                                       "Origin ID",                              &
                                                       "-",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)


                        end if



                        !(Properties)
                        !Just if in all origins are the same amount of properties and in the same
                        !order, writes the properties.
                        CurrentOrigin => Me%FirstOrigin 
    DoCatch:            do while (associated(CurrentOrigin))
                            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                FirstProperty   => CurrentOrigin%FirstProperty
                                nGroupProp      =  CurrentOrigin%nProperties
                                exit DoCatch
                            endif
                            CurrentOrigin => CurrentOrigin%Next
                        enddo DoCatch


                        nProp = 1
                        do while (nProp <= nGroupProp)

                            !Point to the right instant
                            CurrentProperty => FirstProperty
                            do iProp = 1, nProp-1
                                CurrentProperty => CurrentProperty%Next
                            enddo

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
                            do while (associated(CurrentOrigin))

                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                    CurrentPartic     => CurrentOrigin%FirstPartic
                                    do while (associated(CurrentPartic))
                                        Matrix1D(nP)  = CurrentPartic%Concentration(nProp)
                                        CurrentPartic => CurrentPartic%Next
                                        nP = nP + 1
                                    enddo            
                                endif
            
                                CurrentOrigin => CurrentOrigin%Next
                            enddo      

                            !HDF 5
                            call HDF5WriteData        (Me%ObjHDF5,                    &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/"//trim(CurrentProperty%Name), &
                                                       trim(CurrentProperty%Name),               &
                                                       trim(CurrentProperty%Units),              &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)

                            nProp = nProp + 1

                        enddo


    iplume2:            if (Me%State%ComputePlume .and. JetTotalParticles> 0) then
                        
                            allocate   (Matrix1DX(JetTotalParticles))
                            allocate   (Matrix1DY(JetTotalParticles))

                            Matrix1DX(:) = FillValueReal
                            Matrix1DY(:) = FillValueReal

                            call HDF5SetLimits  (Me%ObjHDF5, 1, JetTotalParticles, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR310'

                            FirstParticle = 1
                            CurrentOrigin => Me%FirstOrigin
    iplume3:                do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                
                                    if (CurrentOrigin%State%ComputePlume) then

                                        call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                                                             OutPutLines = OutPutLines, STAT=STAT_CALL)

                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR320'

                                        if (Me%Grid%GridType == Grid1D) then
                                            Matrix1DX(FirstParticle:FirstParticle-1+OutPutLines) = &
                                                      OutPutMatrix(1:OutPutLines,2) + Xorig
                                            Matrix1DY(FirstParticle:FirstParticle-1+OutPutLines) = &
                                                      OutPutMatrix(1:OutPutLines,3) + Yorig
                                        else
                                            do np = 1, OutPutLines
                                                Position%x = OutPutMatrix(np,2)
                                                Position%y = OutPutMatrix(np,3)
                                                call Convert_XY_CellIJ(Position)
                                                call Convert_CellIJ_IJ(Position)
                                                Matrix1DX(FirstParticle-1+np) = GeographicCoordinates (Position, 1)
                                                Matrix1DY(FirstParticle-1+np) = GeographicCoordinates (Position, 2)
                                            enddo
                                        endif

                                        call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)
                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR330'

                                        FirstParticle = FirstParticle + OutPutLines
                                    endif
                            
                                endif

                                CurrentOrigin => CurrentOrigin%Next

                            enddo iplume3
                            

                            call HDF5WriteData  (Me%ObjHDF5, "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume X", "Plume X",  "m", Array1D = Matrix1DX,  &
                                                OutputNumber = OutPutNumber, STAT = STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR340'

                            call HDF5WriteData  (Me%ObjHDF5, "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume Y", "Plume Y",  "m", Array1D = Matrix1DY,  &
                                                OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR350'

                            FirstParticle = 1
                            CurrentOrigin => Me%FirstOrigin
    iplume4:                do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                                
                                    if (CurrentOrigin%State%ComputePlume) then

                                        call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                                                             OutPutLines = OutPutLines, STAT=STAT_CALL)

                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR360'

                                        Matrix1DX(FirstParticle:FirstParticle-1+OutPutLines) = OutPutMatrix(1:OutPutLines,4)
                                        Matrix1DY(FirstParticle:FirstParticle-1+OutPutLines) = OutPutMatrix(1:OutPutLines,6)

                                        call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)
                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR370'

                                        FirstParticle = FirstParticle + OutPutLines
                                    endif
                            
                                endif

                                CurrentOrigin => CurrentOrigin%Next

                            enddo iplume4


                            call HDF5WriteData  (Me%ObjHDF5, "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume Z", "Plume Z",  "m", Array1D = Matrix1DX,  &
                                                OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR380'

                            call HDF5WriteData  (Me%ObjHDF5, "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume Conc", "Plume Conc",  "a",                 &
                                                Array1D = Matrix1DY, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR390'


                            deallocate   (Matrix1DX)
                            deallocate   (Matrix1DY)

                        endif iplume2


                        deallocate  (Matrix1D)

                    endif
                    enddo
                    endif AllAtOnes

                
                    if (Me%State%Monitor) then
                        call WriteMonitorOutput ()
                    endif


                    !Flushes All pending HDF5 commands
                    call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR400'

                endif i1


                !Increments Output number
                Me%OutPut%NextOutPut = OutPutNumber + 1

                !Verifies if all outputs are done (necessary for partic DT smaller then global DT)
                if (Me%OutPut%NextOutPut > size(Me%OutPut%OutTime)) then
                    Me%Output%Write_ = .false.
                endif

            endif  TOut
        endif


    end subroutine ParticleOutput

    !--------------------------------------------------------------------------

    subroutine WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY,                &
                                    StringX, StringY, Units, OutputNumber)

        !Arguments------------------------------------------------------------
        type (T_Origin), pointer            :: CurrentOrigin
        real(8),   dimension(:), pointer    :: Matrix1DX, Matrix1DY
        character(Len=*)                    :: StringX, StringY, Units
        integer                             :: OutputNumber

        !Local-----------------------------------------------------------------
        real,      dimension(:), pointer    :: NodeX, NodeY
        real(8),   dimension(:), pointer    :: Envelope1DX, Envelope1DY
        integer,   dimension(:), pointer    :: BoundaryNodes
        character(Len=StringLength)         :: StringXaux, StringYaux
        real                                :: AverageX, AverageY
        integer                             :: ObjTriangulation = 0, NumberOfBoundaryNodes
        integer                             :: STAT_CALL, NumberOfNodes, i, j


        !Begin-----------------------------------------------------------------

        StringXaux = trim(StringX)//' envelope'
        StringYaux = trim(StringY)//' envelope'

        NumberOfNodes = Size(Matrix1DX)

        if (NumberOfNodes > 3) then

            allocate(NodeX(1:NumberOfNodes), NodeY(1:NumberOfNodes))

            AverageX =sum(Matrix1DX) / real(NumberOfNodes)
            AverageY =sum(Matrix1DY) / real(NumberOfNodes)

            NodeX(:) = Matrix1DX(:) - AverageX
            NodeY(:) = Matrix1DY(:) - AverageY

            call ConstructTriangulation(ObjTriangulation,                              &
                                        NumberOfNodes = NumberOfNodes, NodeX = NodeX,  &
                                        NodeY = NodeY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangian - ERR10'

            call GetNumberOfBoundaryNodes (ObjTriangulation, NumberOfBoundaryNodes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangian - ERR20'

            allocate(Envelope1DX  (1:NumberOfBoundaryNodes),                                &
                     Envelope1DY  (1:NumberOfBoundaryNodes),                                &
                     BoundaryNodes(1:NumberOfBoundaryNodes))


            call GetBoundaryNodes (ObjTriangulation, BoundaryNodes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangian - ERR30'
            !call GetOuterLine

            do i=1, NumberOfBoundaryNodes
                j = BoundaryNodes(i)
                Envelope1DX(i) = NodeX(j) + AverageX
                Envelope1DY(i) = NodeY(j) + AverageY
            enddo

            call HDF5SetLimits  (Me%ObjHDF5, 1, NumberOfBoundaryNodes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangian - ERR40'


            !HDF 5
            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/"//trim(StringXaux), &
                                trim(StringXaux),  trim(Units), Array1D = Envelope1DX,                     &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangian - ERR50'

            !HDF 5
            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(CurrentOrigin%Name)//"/"//trim(StringYaux), &
                                trim(StringYaux),  trim(Units), Array1D = Envelope1DY,                     &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangian - ERR60'

            deallocate(NodeX, NodeY)

            deallocate(Envelope1DX, Envelope1DY, BoundaryNodes)

            call KillTriangulation(ObjTriangulation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangian - ERR70'

        endif

    end subroutine WriteOriginEnvelope

    !--------------------------------------------------------------------------

    subroutine OutputRestartFile
        
        !Local-----------------------------------------------------------------
        real                                :: Year, Month, Day, Hour, Minute, Second
        logical                             :: WriteFinal

        !----------------------------------------------------------------------

        if(Me%ExternalVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then


            call ExtractDate(Me%ExternalVar%Now,                         &
                             Year = Year, Month  = Month,  Day    = Day, &
                             Hour = Hour, Minute = Minute, Second = Second)

            WriteFinal = .true.
!#ifdef _CGI_
            if (Me%RunOnline) then
                WriteFinal = .false.
            endif
!#endif
            if (WriteFinal) call WriteFinalPartic

            Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1

            call SetError(WARNING_, INTERNAL_, "Lagrangian restart file saved          : ", &
                          Year, Month, Day, Hour, Minute, Second)

        end if

    end subroutine OutputRestartFile
    
    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,dimension(:, :, :, :), pointer                  :: GridConc
        real, dimension  (:, :, :), pointer                     :: GridConc3D
        integer                                                 :: STAT_CALL
        type (T_Property), pointer                              :: CurrentProperty
        integer                                                 :: iProp
        integer                                                 :: ILB, IUB, JLB, JUB, KLB, KUB


        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        KLB    = Me%ExternalVar%Size%KLB
        KUB    = Me%ExternalVar%Size%KUB

        !Allocates auxiliar variable
        allocate (GridConc3D (ILB:IUB, JLB:JUB, KLB:KUB         ))

        !Fills Grid concentration
        call FillGridConcentration  (1, GridConc)

        iProp = 0
        CurrentProperty => Me%FirstOrigin%FirstProperty
        do while (associated(CurrentProperty))
            iProp = iProp + 1
            
            GridConc3D = GridConc(:, :, :, iProp)

            call WriteTimeSerie(Me%ObjTimeSerie,                                         &
                                Data3D = GridConc3D,                                     &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangian - ERR01'

            CurrentProperty => CurrentProperty%Next
        enddo

        !Deallocates Temporary Matrixes
        deallocate (GridConc3D)
        deallocate (GridConc  )

    
    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------

    subroutine FillGridConcentration (iGroup, GridConc, GridTracerNumber,                &
                                      GridMaxTracer, GridBottomConc)

        !Arguments-------------------------------------------------------------
        integer                                            :: iGroup
        real,    dimension(:, :, :, :), pointer            :: GridConc
        integer, dimension(:, :, :   ), pointer, optional  :: GridTracerNumber
        real,    dimension(:, :, :, :), pointer, optional  :: GridMaxTracer
        real,    dimension(:, :,    :), pointer, optional  :: GridBottomConc

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                        :: CurrentOrigin
        type (T_Partic), pointer                        :: CurrentPartic
        type (T_Property), pointer                      :: FirstProperty, CurrentProperty
        real(8), dimension(:, :, :   ), pointer         :: GridVolume
        real,    dimension(:, :, :, :), pointer         :: GridMass
        real,    dimension(:, :, :, :), pointer         :: GridGeometricMass      
        real,    dimension(:, :, :   ), pointer         :: GridBottomMass
        real,    dimension(:, :      ), pointer         :: DUX, DVY
        real,    dimension(:), pointer                  :: MeanConc, GeometricMeanConc
        real,    dimension(:), pointer                  :: MinConc, MassVolCel, AmbientConc,AmbientGeoMass
        real,    dimension(:), pointer                  :: MassCoef
        integer, dimension(:), pointer                  :: Geometric
        integer, dimension(:, :, :, :), pointer         :: LocalGridTracerNumber
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                         :: WS_KLB, WS_KUB
        integer                                         :: i, j, k, nProp, status
        integer                                         :: iV, jV, kV, Delta, iAP
        integer                                         :: iInf, iSup
        integer                                         :: jInf, jSup
        integer                                         :: kInf, kSup
        real(8)                                         :: VolumeTotal, Coef
        real(8)                                         :: DiffVolCel
        type (T_Position)                               :: DummYPos    
        logical                                         :: Deposition, FoundSediment
        integer                                         :: Sediment_ID, np
        logical                                         :: GeometricPropertiesPresent


        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        KLB    = Me%ExternalVar%Size%KLB
        KUB    = Me%ExternalVar%Size%KUB
        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB
        WS_KLB = Me%ExternalVar%WorkSize%KLB
        WS_KUB = Me%ExternalVar%WorkSize%KUB
        nProp  = Me%FirstOrigin%nProperties
        
        GeometricPropertiesPresent = .false.

        Deposition = .false.

        CurrentOrigin => Me%FirstOrigin 
Catch:  do while (associated(CurrentOrigin))
            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                FirstProperty   => CurrentOrigin%FirstProperty
                nProp           =  CurrentOrigin%nProperties
                if (CurrentOrigin%State%Deposition) Deposition = .true.
                exit Catch
            endif
            CurrentOrigin => CurrentOrigin%Next
        enddo Catch

        !Allocate GridVolume, GridMass    
        allocate (GridVolume           (ILB:IUB, JLB:JUB, KLB:KUB         ))
        allocate (GridConc             (ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp))
        allocate (GridMass             (ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp))
        allocate (GridGeometricMass    (ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp))
        allocate (GeometricMeanConc    (                           1:nProp))
        allocate (MeanConc             (                           1:nProp))
        allocate (AmbientConc          (                           1:nProp))
        allocate (AmbientGeoMass       (                           1:nProp))
        allocate (MinConc              (                           1:nProp))
        allocate (MassVolCel           (                           1:nProp))
        allocate (MassCOef             (                           1:nProp))
        allocate (Geometric            (                           1:nProp))        
        allocate (LocalGridTracerNumber(ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp))
               
        GridVolume           (:,:,:)   = 0.        
        GridMass             (:,:,:,:) = 0.
        GridGeometricMass    (:,:,:,:) = 0.
        GridConc             (:,:,:,:) = 0.
        MeanConc             (:)       = 0.
        GeometricMeanConc    (:)       = 0.
        AmbientConc          (:)       = 0.
        AmbientGeoMass       (:)       = 0.
        MinConc              (:)       = 0.
        MassVolCel           (:)       = 0.
        Geometric            (:)       = 0
        MassCoef             (:)       = 0.
        LocalGridTracerNumber(:,:,:,:) = 0.        

        if (Deposition .and. present(GridBottomConc)) then
            allocate (GridBottomMass(ILB:IUB, JLB:JUB, 1:nProp))
            allocate (GridBottomConc(ILB:IUB, JLB:JUB, 1:nProp))
            GridBottomMass(:,:,:) = 0.
            GridBottomConc(:,:,:) = 0.
        endif

        if (Me%OutPut%ConcMaxTracer .and. present(GridMaxTracer)) then
            allocate (GridMaxTracer(ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp))
            GridMaxTracer(:,:,:,:) = 0.
        endif

        iAP = 1
        CurrentProperty => FirstProperty
        do while (associated(CurrentProperty))
            Geometric (iAP) = CurrentProperty%Geometric
            if (Geometric (iAP) .gt. 0) GeometricPropertiesPresent = .true.
            iAP = iAP + 1
            CurrentProperty => CurrentProperty%Next
        enddo
       
        if (present(GridTracerNumber)) then
            allocate (GridTracerNumber(ILB:IUB, JLB:JUB, KLB:KUB))
            GridTracerNumber(:,:,:) = 0.
        end if
        

        !Integrates the Volume and the Mass in each GridCell
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

        if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then

            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))
                i = CurrentPartic%Position%I
                j = CurrentPartic%Position%J

                CurrentPartic => CurrentPartic%Next
            enddo


            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

           
                i = CurrentPartic%Position%I
                j = CurrentPartic%Position%J

cd1:            if (.not. CurrentPartic%Deposited) then

                    k = CurrentPartic%Position%K

                    if (present(GridTracerNumber) ) then
                        GridTracerNumber(i, j, k) = GridTracerNumber(i, j, k) + 1
                    endif


                    if (Me%OutPut%ConcMaxTracer .and. present(GridMaxTracer)) then

                        !In this grid is stored the maximum concentration of all tracers present in the cell
                        GridMaxTracer(i, j, k, :) = max(GridMaxTracer(i, j, k, :), CurrentPartic%Concentration(:))

                    endif

                    !Particle fits inside Grid Cell?    
cd2:                if (CurrentPartic%Geometry%Volume <= Me%ExternalVar%VolumeZ(i, j, k)) then

                        GridMass  (i, j, k, :) = GridMass  (i, j, k, :) + CurrentPartic%Mass(:)

                        iAP = 1
                        CurrentProperty => FirstProperty
                        do while (associated(CurrentProperty))
                            if ((Geometric(iAP) == 1) .and. (log10(CurrentPartic%Mass(iAP)) .gt. AllmostZero)) then
                                LocalGridTracerNumber(i, j, k, iAP) = LocalGridTracerNumber(i, j, k, iAP) + 1
                                
                                GridGeometricMass  (i, j, k, iAP) =                         &
                                GridGeometricMass  (i, j, k, iAP) + log10(CurrentPartic%Mass(iAP))
                            end if
                            iAP = iAP + 1
                            CurrentProperty => CurrentProperty%Next
                        enddo

                        GridVolume(i, j, k)    = GridVolume(i, j, k)    + CurrentPartic%Geometry%Volume

                    else  cd2

                        VolumeTotal = Me%ExternalVar%VolumeZ(i, j, k)
                        Delta = 0
                        do while (VolumeTotal < CurrentPartic%Geometry%Volume)

                            Delta = Delta + 1

                            iInf  = max (i-Delta, WS_ILB)
                            jInf  = max (j-Delta, WS_JLB)
                            kInf  = max (k-Delta, WS_KLB)

                            iSup  = min (i+Delta, WS_IUB)
                            jSup  = min (j+Delta, WS_JUB)
                            kSup  = min (k+Delta, WS_KUB)

                            if (iInf == WS_ILB .and. iSup == WS_IUB .and. &
                                jInf == WS_JLB .and. jSup == WS_JUB .and. &
                                kInf == WS_KLB .and. kSup == WS_KUB) then
                        
                                write(*, *)'Particle bigger then domain'
                                stop       'FillGridConcentration - ModuleLagrangian - ERR01'

                            endif

                            VolumeTotal = 0.
                            do iV = iInf, iSup
                            do jV = jInf, jSup
                            do kV = kinf, kSup
                                if (Me%ExternalVar%OpenPoints3D(iV, jV, kV) == OpenPoint) then
                                    VolumeTotal = VolumeTotal +                              &
                                                  Me%ExternalVar%VolumeZ(iV, jV, kV)
                                endif
                            enddo
                            enddo
                            enddo

                        enddo
                
                        do iV = iInf, iSup
                        do jV = jInf, jSup
                        do kV = kinf, kSup
                            if (Me%ExternalVar%OpenPoints3D(iV, jV, kV) == OpenPoint) then
                                Coef        = Me%ExternalVar%VolumeZ(iV, jV, kV) / VolumeTotal
                                MassCoef(:) = CurrentPartic%Mass(:) * Coef

                                GridMass  (iV, jV, kV, :)        = GridMass  (iV, jV, kV, :) + MassCoef(:)

                                iAP = 1
                                CurrentProperty => FirstProperty
                                do while (associated(CurrentProperty))
                                    if ((Geometric(iAP) == 1) .and. (log10(MassCoef(iAP)) .gt. AllmostZero)) then
                                        LocalGridTracerNumber(iV, jV, kV, iAP) = LocalGridTracerNumber(iV, jV, kV, iAP) + 1
                                        
                                        GridGeometricMass(iV, jV, kV, iAP) = GridGeometricMass  (iV, jV, kV, iAP) +   &
                                                                            log10(MassCoef(iAP))

                                    end if
                                    iAP = iAP + 1
                                    CurrentProperty => CurrentProperty%Next
                                enddo
                                                                                                 
                                GridVolume       (iV, jV, kV)    = GridVolume(iV, jV, kV)    +      &
                                                                 CurrentPartic%Geometry%Volume * Coef

                            endif
                        enddo
                        enddo
                        enddo
                
                    endif cd2

                else  cd1 ! The particle is deposited in the bottom
                !In this case no test is made to verify if the the particle occupies more then one cell
                !to maintain the algothim simple.

                      
                    GridBottomMass  (i, j, :) = GridBottomMass  (i, j, :) + CurrentPartic%Mass(:)

                endif cd1

                CurrentPartic => CurrentPartic%Next
            enddo
            
            if (.not. associated(FirstProperty)) FirstProperty => CurrentOrigin%FirstProperty

        endif

        CurrentOrigin => CurrentOrigin%Next
        enddo CurrOr


        if (Deposition .and. present(GridBottomConc)) then !Find if exist the 'sediment' property

            !Ambient Concentration of the place of the particle
            Sediment_ID      = 0
            FoundSediment = .false.
            CurrentProperty => FirstProperty
            do while (associated(CurrentProperty) .and. .not. FoundSediment)
                Sediment_ID = Sediment_ID + 1
                if (CurrentProperty%ID == Sediment) FoundSediment = .true.
                CurrentProperty => CurrentProperty%Next
            enddo
        
        endif



        !Fills up Grid Concentration
        do k = WS_KLB, WS_KUB
        do j = WS_JLB, WS_JUB
        do i = WS_ILB, WS_IUB

            if (Me%ExternalVar%WaterPoints3D (i, j, k) == WaterPoint) then

                !Mean Concentration of the particle
                MeanConc(:) = 0.0
                GeometricMeanConc(:) = 0.0
                if (GridVolume(i, j, k) > 0.0) then

                    iAP = 1
                    CurrentProperty => FirstProperty
                    do while (associated(CurrentProperty))
                        MeanConc(iAP)          = GridMass  (i, j, k, iAP) / GridVolume(i, j, k)

                        if ((Geometric(iAP)==1) .and. (LocalGridTracerNumber (i, j, k, iAP) > 0)) then
                            GeometricMeanConc(iAP) =                                                             &
                            (10.**(GridGeometricMass  (i, j, k, iAP) / real(LocalGridTracerNumber (i, j, k, iAP))))   &
                            * real(LocalGridTracerNumber (i, j, k, iAP)) / GridVolume(i, j, k)
                        end if
                        iAP = iAP + 1
                        CurrentProperty => CurrentProperty%Next
                    enddo
                
                end if

                !Ambient Concentration of the place of the particle
                iAP = 1
                CurrentProperty => FirstProperty
                do while (associated(CurrentProperty))
                
                    DummYPos%I = i
                    DummYPos%J = j
                    DummYPos%K = k
                    
                    call GetAmbientConcentration (CurrentProperty,                          &
                                                  DummYPos,                                 &
                                                  AmbientConc(iAP))

                    MinConc (iAP) = CurrentProperty%Min_concentration

                    iAP = iAP + 1
                    CurrentProperty => CurrentProperty%Next
                enddo


                select case (Me%OutPut%OutPutConcType)

                case (Maximum)


                    !Metodo Max : Admite-se a concentracao igual ao valor maximo de entre 2 valores: 
                    !          - massa total de tracadores presente na celula a dividir pelo volume da celula
                    !          - concentracao media (aritmética ou geométrica) dos tracadores
                                   
                    MassVolCel(:)        = GridMass(i, j, k, :) /                        &
                                           Me%ExternalVar%VolumeZ(i, j, k)
                    
                    where (Geometric(:) == 1)
                        GridConc(i, j, k, :) = max(MassVolCel(:), GeometricMeanConc(:))
                    elsewhere (Geometric(:) == 0)
                        GridConc(i, j, k, :) = max(MassVolCel(:), MeanConc(:)) 
                    end where 
                    
                case (Mean)

                    !Metodo 2 : Se o volume total de tracadores for menor
                    !      que o volume da célula i,j,k, entao a concentracao nesta  e igual a uma media 
                    !      entre a concentracao media (aritmética ou geométrica)dos tracadores e a concentracao ambiente
                    !      caso contrario a concentracao na celula i,j,k e igual a concentracao media 
                    !      dos tracadores

                    if (Me%ExternalVar%VolumeZ(i, j, k) >= GridVolume(i, j, k)) then 

                        DiffVolCel = Me%ExternalVar%VolumeZ(i, j, k) - GridVolume(i, j, k)  
                            where (AmbientConc(:) > 0)            &
                            AmbientGeoMass(:) = log10(DiffVolCel * AmbientConc(:))                        
                
                            where ((Geometric(:) == 1) .and. (LocalGridTracerNumber (i, j, k, :) > 0))                     &
                                GridConc(i, j, k, :) = (    10**                                                        & 
                                                                (                                                       &
                                                                  (                                                     &
                                                                   AmbientGeoMass(:)                                    &
                                                                   +                                                    & 
                                                                   GridGeometricMass  (i, j, k, :)                      &
                                                                   )                                                    &
                                                                   /                                                    &
                                                                   (real(LocalGridTracerNumber (i, j, k, :)+1))            &
                                                                 )                                                      &
                                                        )                                                               &
                                                        *  (real(LocalGridTracerNumber (i, j, k, :)+1))                    &
                                                        /                                                               &
                                                       Me%ExternalVar%VolumeZ(i, j, k)                                  
                           where (Geometric(:) == 0)                                                                    &
                                                                       
                                GridConc(i, j, k, :) = (DiffVolCel * AmbientConc(:) + GridMass  (i, j, k, :)) / &  
                                                        Me%ExternalVar%VolumeZ(i, j, k)
                    else

                        where (Geometric(:) == 0)                                            
                            GridConc(i, j, k, :) = MeanConc(:)
                        elsewhere (Geometric(:) == 1)
                            GridConc(i, j, k, :) = GeometricMeanConc(:)
                        end where
                    endif

                end select

                where (GridConc(i, j, k, :) .lt. MinConc(:))                         &
                    GridConc(i, j, k, :) = AmbientConc(:)

                where (GridConc(i, j, k, :) == 0.)                                   &
                    GridConc(i, j, k, :) = AmbientConc(:)


            else

                GridConc(i, j, k, :) = null_real

            end if

        end do
        end do
        end do
        
        CurrentOrigin => Me%FirstOrigin 
dw1:    do while (associated(CurrentOrigin))
            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)  .and. CurrentOrigin%State%Oil) then
                CurrentProperty  => CurrentOrigin%FirstProperty
                iAP = 1
dw2:            do while (associated(CurrentProperty))
                    if (CurrentProperty%ID == Oil_) then
                        !To be able to do statistics analysis of the oil thickness 
                        GridConc(:, :, WS_KUB, iAP) = CurrentOrigin%GridThickness(:,:)
                        exit dw2                    
                    endif
                    iAP = iAP + 1
                    CurrentProperty => CurrentProperty%Next                                        
                enddo dw2
            endif
            CurrentOrigin => CurrentOrigin%Next
        enddo dw1
        

cd3:    if (Deposition .and. present(GridBottomConc)) then ! fills up the bottom concentration in a simplified way

            !Gets Horizontal Grid
            call GetHorizontalGrid(Me%ObjHorizontalGrid, DUX = DUX, DVY = DVY, STAT = status)
            if (status /= SUCCESS_) call SetError(FATAL_, INTERNAL_, 'FillGridConcentration - ModuleLagrangian - ERR01')

    
            do j = WS_JLB, WS_JUB
            do i = WS_ILB, WS_IUB

                if (Me%ExternalVar%WaterPoints3D (i, j, WS_KUB) == WaterPoint) then

                    if (FoundSediment) then
                        do np = 1, SIZE(GridBottomConc, DIM = 3)

                            if (np /= Sediment_ID) then 
                                !Mass contaminant / Mass sediment
                                if (GridBottomMass(i, j, Sediment_ID) > 0) then
                                    GridBottomConc(i, j, np) = GridBottomMass(i, j, np) / GridBottomMass(i, j, Sediment_ID)
                                else
                                    GridBottomConc(i, j, np) = 0.
                                endif
                            else
                                ! Mass of sediment / m^2
                                GridBottomConc(i, j, Sediment_ID) = GridBottomMass(i, j, Sediment_ID) / DUX(I, j) / DVY(i, j)
                            endif
                       
                        enddo 
                    else !all properties are written in mass / m^2
                        ! Mass / m^2
                        GridBottomConc(i, j, :) = GridBottomMass(i, j, :) / DUX(I, j) / DVY(i, j)

                    endif

                endif

            enddo
            enddo

            !UnGets Horizontal Grid
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DUX,               &
                                   STAT = status)
            if (status /= SUCCESS_) call SetError(FATAL_, INTERNAL_, 'FillGridConcentration - ModuleLagrangian - ERR02')

            !UnGets Horizontal Grid
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, DVY,               &
                                   STAT = status)
            if (status /= SUCCESS_) call SetError(FATAL_, INTERNAL_, 'FillGridConcentration - ModuleLagrangian - ERR03')

            deallocate (GridBottomMass)


        endif cd3

        deallocate (GridVolume        )
        deallocate (GridMass          )
        deallocate (GridGeometricMass )        
        deallocate (GeometricMeanConc )       
        deallocate (MeanConc          )
        deallocate (AmbientConc       )
        deallocate (AmbientGeoMass    )
        deallocate (MinConc           )
        deallocate (MassVolCel        )
        deallocate (MassCoef          )
        deallocate (Geometric         )
        
        deallocate (LocalGridTracerNumber)

        
    end subroutine FillGridConcentration

    !--------------------------------------------------------------------------


    subroutine WriteGridConcentration (iGroup, GridConc, GridTracerNumber,               &
                                       GridMaxTracer, GridBottomConc)

        !Arguments-------------------------------------------------------------
   
        integer                                     :: iGroup
        real, dimension(:, :, :, :), pointer     :: GridConc
        integer, dimension(:, :, :   ), pointer     :: GridTracerNumber
        real,    dimension(:, :, :, :), pointer     :: GridMaxTracer
        real,    dimension(:, :,    :), pointer     :: GridBottomConc

        !Local-----------------------------------------------------------------
        real, dimension(:, :, :), pointer           :: GridConc3D 
        real, dimension(:, :   ), pointer           :: GridConc2D 
        integer                                     :: iAP
        type (T_Origin),   pointer                  :: CurrentOrigin
        type (T_Property), pointer                  :: CurrentProperty, FirstProperty
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB
        character(StringLength)                     :: AuxChar
        logical                                     :: Deposition

        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        KLB    = Me%ExternalVar%Size%KLB
        KUB    = Me%ExternalVar%Size%KUB

        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB
        WS_KLB = Me%ExternalVar%WorkSize%KLB
        WS_KUB = Me%ExternalVar%WorkSize%KUB

        !Allocates auxiliar variable
        allocate (GridConc3D (ILB:IUB, JLB:JUB, KLB:KUB         ))

        !Writes the Group to an auxiliar string
        write (AuxChar, fmt='(i3)') iGroup

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WS_ILB, WS_IUB, WS_JLB, WS_JUB,     &
                              WS_KLB, WS_KUB)

        Deposition = .false.

        CurrentOrigin => Me%FirstOrigin 
Catch:  do while (associated(CurrentOrigin))
            if (CurrentOrigin%GroupID == Me%GroupIDs(iGroup)) then
                FirstProperty   => CurrentOrigin%FirstProperty
                if (CurrentOrigin%State%Deposition) Deposition = .true.
                exit Catch
            endif
            CurrentOrigin => CurrentOrigin%Next
        enddo Catch

        if (Deposition) allocate (GridConc2D (ILB:IUB, JLB:JUB      ))

        iAP = 1
        CurrentProperty => FirstProperty
        do while (associated(CurrentProperty))
        
            GridConc3D(:,:,:) = GridConc(:, :, :, iAP)


            !HDF 5
            call HDF5WriteData        (Me%ObjHDF5,                       &
                                       "/Results/Group_"//trim(adjustl(AuxChar))    &
                                       //"/Data_3D/"//trim(CurrentProperty%Name),   &
                                       trim(CurrentProperty%Name),                  &
                                       trim(CurrentProperty%Units),                 &
                                       Array3D = GridConc3D,                        &
                                       OutputNumber = Me%OutPut%NextOutPut)

            if (Me%OutPut%ConcMaxTracer) then

                GridConc3D(:,:,:) = GridMaxTracer(:, :, :, iAP)


                !HDF 5
                call HDF5WriteData    (Me%ObjHDF5,                                            &
                                       "/Results/Group_"//trim(adjustl(AuxChar))              &
                                       //"/Data_3D_MaxTracer/"//trim(CurrentProperty%Name),   &
                                       trim(CurrentProperty%Name),                            &
                                       trim(CurrentProperty%Units),                           &
                                       Array3D = GridConc3D,                                  &
                                       OutputNumber = Me%OutPut%NextOutPut)

            endif



            if (Deposition) then

                GridConc2D(:,:) = GridBottomConc(:, :, iAP)

                !HDF 5
                call HDF5WriteData        (Me%ObjHDF5,                       &
                                           "/Results/Group_"//trim(adjustl(AuxChar))    &
                                           //"/Bottom/"//trim(CurrentProperty%Name),    &
                                           trim(CurrentProperty%Name),                  &
                                           trim(CurrentProperty%Units),                 &
                                           Array2D = GridConc2D,                        &
                                           OutputNumber = Me%OutPut%NextOutPut)

            endif


            iAP = iAP + 1
            CurrentProperty => CurrentProperty%Next
        enddo

        !HDF 5
        call HDF5WriteData        (Me%ObjHDF5,                                  &
                                   "/Results/Group_"//trim(adjustl(AuxChar))    &
                                   //"/Data_3D/"//"Number",                     &
                                   "Number", "a",                               &
                                   Array3D = GridTracerNumber,                  &
                                   OutputNumber = Me%OutPut%NextOutPut)


        deallocate (GridTracerNumber)
        deallocate (GridConc3D      )
        deallocate (GridConc        )
        
        if (Me%OutPut%ConcMaxTracer) deallocate (GridMaxTracer)

        if (Deposition) then

            deallocate (GridConc2D    )
            deallocate (GridBottomConc)

        endif
       
    end subroutine WriteGridConcentration

    !--------------------------------------------------------------------------

    subroutine WriteFrequencyLag (CurrentProperty)

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:, :, :), pointer           :: GridConc3D 
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB, iClass, STAT_CALL
        character(StringLength)                     :: AuxChar1, AuxChar2, AuxChar


        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        KLB    = Me%ExternalVar%Size%KLB
        KUB    = Me%ExternalVar%Size%KUB

        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB
        WS_KLB = Me%ExternalVar%WorkSize%KLB
        WS_KUB = Me%ExternalVar%WorkSize%KUB

        !Allocates auxiliar variable
        allocate (GridConc3D (ILB:IUB, JLB:JUB, KLB:KUB         ))


        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WS_ILB, WS_IUB, WS_JLB, WS_JUB,                &
                              WS_KLB, WS_KUB)

        call GetStatisticClasses(CurrentProperty%Statistic1_ID,                          &
                                 CurrentProperty%ClassesLag, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'WriteFrequencyLag - ModuleLagrangian - ERR01'

        do iClass = 1, CurrentProperty%nClassesLag


            write(AuxChar1, fmt='(E9.2)')CurrentProperty%ClassesLag(iClass, 1)
            write(AuxChar2, fmt='(E9.2)')CurrentProperty%ClassesLag(iClass, 2)

            AuxChar = trim(adjustl(AuxChar1))//"_"//trim(adjustl(AuxChar2))

        
            GridConc3D(:,:,:) = 100. * CurrentProperty%FrequencyLag(:,:,:,iClass)

                        

            call HDF5WriteData   (Me%ObjHDF5, "/Statistics/Lagrangian/"//trim(CurrentProperty%Name)//"/Classes", &
                                  trim(adjustl(AuxChar)),                                &
                                  "-", Array3D = GridConc3D,                             &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFrequencyLag - ModuleLagrangian - ERR02'

        enddo

        call UnGetStatistic(CurrentProperty%Statistic1_ID,                               &
                            CurrentProperty%ClassesLag, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'WriteFrequencyLag - ModuleLagrangian - ERR03'



        deallocate (GridConc3D )
        
       
    end subroutine WriteFrequencyLag

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine WriteOilGridThickness 

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iGroup
        real, dimension(:, :), pointer              :: OilGridThick2D 
        character(StringLength)                     :: AuxChar
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB



        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB

        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB
        WS_KLB = Me%ExternalVar%WorkSize%KLB
        WS_KUB = Me%ExternalVar%WorkSize%KUB

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WS_ILB, WS_IUB, WS_JLB, WS_JUB,     &
                              WS_KLB, WS_KUB)


        !Allocate GridVolume, GridMass
        allocate (OilGridThick2D (ILB:IUB, JLB:JUB))

Group:  do iGroup = 1, Me%nGroups

            !Writes the Group to an auxiliar string
            write (AuxChar, fmt='(i3)') iGroup

            CurrentOrigin => Me%FirstOrigin
CurrOr:     do while (associated(CurrentOrigin))


                !Just writes the output if there are particle
                if (CurrentOrigin%nParticle > 0) then
                    OilGridThick2D = CurrentOrigin%GridThickness * 1000.0                       &
                                     * (1 - Me%ExternalVar%VWaterContent)


                    !HDF 5
                    call HDF5WriteData        (Me%ObjHDF5,                       &
                                               "/Results/"//trim(CurrentOrigin%Name)        &
                                               //"/Data_2D/Thickness_2D",                   &
                                               "Thickness_2D",                              &
                                               "mm",                                        &
                                               Array2D = OilGridThick2D,                    &
                                               OutputNumber = Me%OutPut%NextOutPut)


                endif

            CurrentOrigin => CurrentOrigin%Next
            enddo CurrOr

        enddo Group

        deallocate (OilGridThick2D )


    end subroutine WriteOilGridThickness


    !--------------------------------------------------------------------------

    subroutine WriteOilGridConcentration 

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: iGroup
        real, dimension(:, :), pointer              :: OilGridConc2D 
        character(StringLength)                     :: AuxChar
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB



        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB

        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB
        WS_KLB = Me%ExternalVar%WorkSize%KLB
        WS_KUB = Me%ExternalVar%WorkSize%KUB

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WS_ILB, WS_IUB, WS_JLB, WS_JUB,     &
                              WS_KLB, WS_KUB)


        !Allocate GridVolume, GridMass
        allocate (OilGridConc2D (ILB:IUB, JLB:JUB))

Group:  do iGroup = 1, Me%nGroups

            !Writes the Group to an auxiliar string
            write (AuxChar, fmt='(i3)') iGroup

            CurrentOrigin => Me%FirstOrigin
CurrOr:     do while (associated(CurrentOrigin))


                !Just writes the output if there are particle
                if (CurrentOrigin%nParticle > 0) then
                    OilGridConc2D = CurrentOrigin%OilGridConcentration

                    !HDF 5
                    call HDF5WriteData        (Me%ObjHDF5,                       &
                                               "/Results/"//trim(CurrentOrigin%Name)        &
                                               //"/Data_2D/OilConcentration_2D",            &
                                               "OilConcentration_2D",                       &
                                               "ppm",                                       &
                                               Array2D = OilGridConc2D,                     &
                                               OutputNumber = Me%OutPut%NextOutPut)


                endif

            CurrentOrigin => CurrentOrigin%Next
            enddo CurrOr

        enddo Group

        deallocate (OilGridConc2D )


    end subroutine WriteOilGridConcentration


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    real(8) function GeographicCoordinates (Position, Direction)

        !Arguments-------------------------------------------------------------
   
        type (T_Position)                           :: Position
        integer                                     :: Direction

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        real                                        :: CellI, CellJ
        real(8)                                     :: Balx, Dx, x1, x2
        !real                                        :: Xorig, Yorig


        !Gets Origin
        !call GetGridOrigin(Me%ObjHorizontalGrid, Xorig, Yorig)

        i       = Position%I
        j       = Position%J
        CellI   = Position%CellI
        CellJ   = Position%CellJ

        if (Direction == 1) then

            Balx  = CellJ - int(CellJ)
            
            x1   = (Me%ExternalVar%XX_IE(i  , j) + Me%ExternalVar%XX_IE(i+1, j  ))/2.
            x2   = (Me%ExternalVar%XX_IE(i, j+1) + Me%ExternalVar%XX_IE(i+1, j+1))/2.

        else

            Balx  = CellI - int(CellI)

            x1   = (Me%ExternalVar%YY_IE(i  , j) + Me%ExternalVar%YY_IE(i  , j+1))/2.
            x2   = (Me%ExternalVar%YY_IE(i+1, j) + Me%ExternalVar%YY_IE(i+1, j+1))/2.

        endif

        Dx   =  x2 - x1
        GeographicCoordinates = x1 + Balx * Dx


    end function GeographicCoordinates 

    !--------------------------------------------------------------------------

    subroutine WriteMonitorOutput 

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer, dimension(:,:,:), pointer          :: MonitorBoxes
        integer                                     :: NumberOfBoxes
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k, Box
        integer                                     :: STAT_CALL
        real, dimension(:, :, :), pointer           :: OutputMatrix
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB

        !Shorten
        ILB    = Me%ExternalVar%Size%ILB
        IUB    = Me%ExternalVar%Size%IUB
        JLB    = Me%ExternalVar%Size%JLB
        JUB    = Me%ExternalVar%Size%JUB
        KLB    = Me%ExternalVar%Size%KLB
        KUB    = Me%ExternalVar%Size%KUB

        WS_ILB = Me%ExternalVar%WorkSize%ILB
        WS_IUB = Me%ExternalVar%WorkSize%IUB
        WS_JLB = Me%ExternalVar%WorkSize%JLB
        WS_JUB = Me%ExternalVar%WorkSize%JUB
        WS_KLB = Me%ExternalVar%WorkSize%KLB
        WS_KUB = Me%ExternalVar%WorkSize%KUB

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WS_ILB, WS_IUB, WS_JLB, WS_JUB,     &
                              WS_KLB, WS_KUB)


        allocate (OutputMatrix (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteMonitorOutput - ModuleLagrangian - ERR00'

        !Get the boxes
        call GetBoxes(Me%ObjMonBox, MonitorBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteMonitorOutput - ModuleLagrangian - ERR01'

        call GetNumberOfBoxes(Me%ObjMonBox, NumberOfBoxes3D = NumberOfBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteMonitorOutput - ModuleLagrangian - ERR01'

        !For every origin writes the instant percentage of the instant volume contributed to 
        !a given monitoring box, refering to the instant volume of the Monitoring Box
        CurrentOrigin => Me%FirstOrigin
        do while (associated(CurrentOrigin))
            OutputMatrix = null_real
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                Box = MonitorBoxes(i, j, k)
                if (Box > 0) then
                    OutputMatrix(i, j, k) = Me%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) / &
                                            Me%Monitor%InstBoxVolume      (Box) * 100.
                endif
            enddo
            enddo
            enddo


            !HDF 5
            call HDF5WriteData         (Me%ObjHDF5,                                &
                                        "/Results/"//trim(CurrentOrigin%Name)//"/InstVolume", &
                                        "InstVolumeContribution",  "%",                       &
                                        Array3D = OutputMatrix,                               &
                                        OutputNumber = Me%OutPut%NextOutPut,       &
                                        STAT = STAT_CALL)

            CurrentOrigin => CurrentOrigin%Next
        enddo


        !For every origin writes the percentage of the integrated volume contributed to 
        !a given monitoring box, refering to the integrated volume of the Monitoring Box
        CurrentOrigin => Me%FirstOrigin
        do while (associated(CurrentOrigin))
            OutputMatrix = null_real
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                Box = MonitorBoxes(i, j, k)
                if (Box > 0) then
                    OutputMatrix(i, j, k) = Me%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID) / &
                                            Me%Monitor%IntgBoxVolume      (Box) * 100.
                endif
            enddo
            enddo
            enddo


            !HDF 5
            call HDF5WriteData         (Me%ObjHDF5,                                &
                                        "/Results/"//trim(CurrentOrigin%Name)//"/IntgVolume", &
                                        "IntgVolumeContribution",  "%",                       &
                                        Array3D = OutputMatrix,                               &
                                        OutputNumber = Me%OutPut%NextOutPut,       &
                                        STAT = STAT_CALL)


            CurrentOrigin => CurrentOrigin%Next
        enddo


        If (Me%State%MonitorPropMass) then

        !For every origin writes the instant percentage of the instant mass contributed to 
        !a given monitoring box, refering to the instant mass of the Monitoring Box
        CurrentOrigin => Me%FirstOrigin
        do while (associated(CurrentOrigin))
            OutputMatrix = null_real
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                Box = MonitorBoxes(i, j, k)
                if (Box > 0) then
                    if (Me%Monitor%InstBoxMass (Box) .EQ. 0. ) then
                        OutputMatrix(i, j, k) = 0
                    else
                        OutputMatrix(i, j, k) = Me%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID) / &
                                                Me%Monitor%InstBoxMass      (Box) * 100.
                    end if
                endif
            enddo
            enddo
            enddo


            !HDF 5
            call HDF5WriteData         (Me%ObjHDF5,                                &
                                        "/Results/"//trim(CurrentOrigin%Name)//"/InstMass", &
                                        "InstMassContribution",  "%",                       &
                                        Array3D = OutputMatrix,                               &
                                        OutputNumber = Me%OutPut%NextOutPut,       &
                                        STAT = STAT_CALL)

            CurrentOrigin => CurrentOrigin%Next
        enddo


        End IF

        !For every monitoring box fills as many cells as the get contributions (in volume)
        !from a given origin. The rest of the cells of the monitoring box is filled with the 
        !value zero
        !Calculates the number of cells in each box
        Me%Monitor%NumberOfCellsPerBox = 0
        do Box = 1,  NumberOfBoxes
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (MonitorBoxes(i, j, k) == Box .and. Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then
                Me%Monitor%NumberOfCellsPerBox(Box) = Me%Monitor%NumberOfCellsPerBox(Box) + 1
            endif
        enddo
        enddo
        enddo
        enddo
        !Number of Cells to fill
        Me%Monitor%NumberOfCellsFromOrigin = 0
        CurrentOrigin => Me%FirstOrigin
        do while (associated(CurrentOrigin))
            do Box = 1,  NumberOfBoxes
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                if (MonitorBoxes(i, j, k) == Box) then
                    Me%Monitor%NumberOfCellsFromOrigin    (Box, CurrentOrigin%ID) = &
                        int(Me%Monitor%NumberOfCellsPerBox(Box) *                   &
                            Me%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) / &
                            Me%Monitor%InstBoxVolume      (Box))
                endif
            enddo
            enddo
            enddo
            enddo
            CurrentOrigin => CurrentOrigin%Next
        enddo
        !Outputs
        OutputMatrix = null_real
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            Box = MonitorBoxes(i, j, k)
            CurrentOrigin => Me%FirstOrigin
            if (Box > 0 .and. Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then

                OutputMatrix(i, j, k) = 0.
                do while (Me%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) == 0)
                    CurrentOrigin => CurrentOrigin%Next
                    if (.not. associated(CurrentOrigin)) exit
                enddo

                if (associated(CurrentOrigin)) then
                    if (Me%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) > 0) then
                        OutputMatrix(i, j, k) = CurrentOrigin%ID
                        Me%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) = &
                            Me%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) - 1
                    endif
                endif
            endif
        enddo
        enddo
        enddo
            


        !HDF 5
        call HDF5WriteData         (Me%ObjHDF5,                          &
                                    "/Results/Monitor",                             &
                                    "VolumeContributedByOrigin",  "%",              &
                                    Array3D = OutputMatrix,                         &
                                    OutputNumber = Me%OutPut%NextOutPut, &
                                    STAT = STAT_CALL)

        !Unget The Boxes
        call UngetBoxDif(Me%ObjMonBox, MonitorBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteMonitorOutput - ModuleLagrangian - ERR02'

        deallocate (OutputMatrix            )

    end subroutine WriteMonitorOutput

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    subroutine KillLagrangian (LagrangianID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: NumberOfBoxes, nB
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL, nUsers
        logical                                     :: WriteFinal

        !----------------------------------------------------------------------
           
        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_)    

cd1:    if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mLAGRANGIAN_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%WritesTimeSerie) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR00'
                endif
                

                if (Me%State%Statistics) then

                    call KillParticleStatistic 

                endif

                !Kills the Light
                if (Me%State%WQM .or. Me%State%T90Variable) then
                    call KillLight 
                endif

                if (Me%OutPut%Write_) then
                    !Closes the transient result file
                    call KillHDF5          (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR01a'
                endif

                WriteFinal = .true.

                if (Me%RunOnline) then
!#ifdef _CGI_
                    WriteFinal = .false.
!#endif
                endif
                !Writes the Final Output
                if (WriteFinal) call WriteFinalPartic(Final = .true.)

                !Kills Monitoring
                if (Me%State%Monitor) then

                    call GetNumberOfBoxes(Me%ObjMonBox, NumberOfBoxes3D = NumberOfBoxes, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                    do nB = 1, NumberOfBoxes
                        call KillTimeSerie (Me%Monitor%ObjTimeSerie(nB), STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR02'
                    enddo

                    deallocate (Me%Monitor%SurfaceBoxVolume           )
                    deallocate (Me%Monitor%InstBoxVolume              )
                    deallocate (Me%Monitor%InstVolumeByOrigin         )
                    deallocate (Me%Monitor%InstBoxMass                )
                    deallocate (Me%Monitor%InstMassByOrigin           )
                    deallocate (Me%Monitor%IntgBoxVolume              )
                    deallocate (Me%Monitor%IntgVolumeByOrigin         )
                    deallocate (Me%Monitor%NumberOfCellsPerBox        )
                    deallocate (Me%Monitor%NumberOfCellsFromOrigin    )
                    deallocate (Me%Monitor%InstBoxLogMass             )         
                    deallocate (Me%Monitor%InstBoxConc                )
                    deallocate (Me%Monitor%NumberOfTracers            )
                    deallocate (Me%Monitor%InstBoxMassFractionByOrigin)
                    deallocate (Me%Monitor%InstLogMassByOrigin        )
                    deallocate (Me%Monitor%NumberOfTracersFromOrigin  )
                    deallocate (Me%Monitor%ContaminationProbability   )    
                    deallocate (Me%Monitor%AverageBoxContaminatedConc )
                    deallocate (Me%Monitor%NbrBoxContaminatedTracers  )
                    deallocate (Me%Monitor%VolBoxContaminatedTracers  )
                    

                    call KillBoxDif (Me%ObjMonBox, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR02b'
                    
                endif

                if (Me%State%EulerianMonitor) then

                    deallocate(Me%EulerianMonitor%Mass)

                    call KillBoxDif (Me%ObjEulMonBox, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR02c'

                end if

                !Kills Oil-Beaching Properties
                if (Me%State%AssociateBeachProb) then
                    if (Me%State%HaveBeachingProbBox) then

                        call KillBoxDif(Me%ObjBeachingProbBox, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR363.'

                    end if

                    deallocate (Me%BeachingProbability        )

                end if
                   

                if (Me%State%Deposition)  then

                    deallocate   (Me%TauErosionGrid )
                    deallocate   (Me%MassSedGrid    )

                endif

                if (Me%State%Filtration)  then

                    deallocate(Me%Filtration%RelativeMassFilter)

                    deallocate(Me%Filtration%MassFiltered)

                endif

                !Kill the OriginList
                call DeallocateOriginList(Me%FirstOrigin, Me%nOrigins)

                !Kill Assimilation
                if (Me%OverLay%Overlay) then

                    call KillAssimilation (Me%ObjAssimilation,  STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR01a'

                    deallocate (Me%OverLay%VelUFinal)
                    deallocate (Me%OverLay%VelVFinal)

                endif

                if (Me%ObjBoxDif /= 0) then
                    call KillBoxDif (Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillLagrangian - ModuleLagrangian - ERR02c'
                endif

                !External Modules
                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime           )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mGRIDDATA_,       Me%ObjGridData       )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap  )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry       )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mMAP_,            Me%ObjMap            )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mHYDRODYNAMIC_,   Me%ObjHydrodynamic   )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mTURBULENCE_,     Me%ObjTurbulence     )
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

                nUsers = DeassociateInstance (mWATERPROPERTIES_,Me%ObjWaterProperties)
                if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'

#ifndef _WAVES_
                if(Me%ObjWaves /= 0)then
                    nUsers = DeassociateInstance (mWAVES_,  Me%ObjWaves)
                    if (nUsers == 0) stop 'KillLagrangian - ModuleLagrangian - ERR01'
                end if
#endif
                call DeallocateInstance

                LagrangianID  = 0
                STAT_         = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_



    end subroutine KillLagrangian

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Lagrangian), pointer                      :: AuxLagrangian      => null()
        type (T_Lagrangian), pointer                      :: PreviousLagrangian => null()

        !Updates pointers
        if (Me%InstanceID == FirstLagrangian%InstanceID) then
            FirstLagrangian => FirstLagrangian%Next
        else
            PreviousLagrangian => FirstLagrangian
            AuxLagrangian      => FirstLagrangian%Next
            do while (AuxLagrangian%InstanceID /= Me%InstanceID)
                PreviousLagrangian => AuxLagrangian
                AuxLagrangian      => AuxLagrangian%Next
            enddo

            !Now update linked list
            PreviousLagrangian%Next => AuxLagrangian%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine DeallocateOriginList (FirstOrigin, nOrigins)

        !Arguments-------------------------------------------------------------
        type (T_origin), pointer                    :: FirstOrigin
        integer                                     :: nOrigins

        !Local-----------------------------------------------------------------
        type (T_origin), pointer                    :: CurrentOrigin => null()

        CurrentOrigin => FirstOrigin
        do while (associated(CurrentOrigin))

            call DeleteOrigin (FirstOrigin, CurrentOrigin, nOrigins)

            CurrentOrigin => FirstOrigin
        enddo


    end subroutine DeallocateOriginList 

    !--------------------------------------------------------------------------

    subroutine WriteFinalPartic(Final)

        !Arguments-------------------------------------------------------------
        logical, optional, intent (IN)               :: Final
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin   => null()
        type (T_Partic), pointer                    :: CurrentPartic   => null()
        type (T_Property), pointer                  :: CurrentProperty => null()
        integer                                     :: UnitID
        integer                                     :: STAT_CALL
        character (Len = Pathlength)                :: filename
        logical                                     :: Final_ = .false.

        if(present(Final)) Final_ = Final

        !Opens File
        call UnitsManager (UnitID, OPEN_FILE, STAT = STAT_CALL)

        !Checks if it's at the end of the run 
        !or !if it's supposed to overwrite the final HDF file
        if (Me%Output%RestartOverwrite .or. Final_) then

            filename = trim(Me%Files%Final)

        else

            filename =  ChangeSuffix(Me%Files%Final,                         &
                               "_"//trim(TimeToString(Me%ExternalVar%Now))//".fin")


        endif

        open   (unit = UnitID, file = trim(filename), status = 'unknown',     &
                form = 'unformatted')
        rewind (unit = UnitID)

        write(UnitID) Me%nOrigins

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            !Writes Origin Information
            write (UnitID) CurrentOrigin%Name
            write (UnitID) CurrentOrigin%ID
            write (UnitID) CurrentOrigin%State%Oil
            write (UnitID) CurrentOrigin%State%Deposition
            write (UnitID) CurrentOrigin%State%Age
            write (UnitID) CurrentOrigin%State%ComputePlume
            write (UnitID) CurrentOrigin%State%PlumeShear
            write (UnitID) CurrentOrigin%State%FarFieldBuoyancy
            write (UnitID) CurrentOrigin%nParticle
            write (UnitID) CurrentOrigin%nProperties

            If (CurrentOrigin%State%Oil)   then
                write (UnitID) CurrentOrigin%GridThickness                      
                write (UnitID) CurrentOrigin%OilGridConcentration
            endif

            !Writes Particle Information
            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                write (UnitID) CurrentPartic%ID

                write (UnitID) CurrentPartic%Position%I
                write (UnitID) CurrentPartic%Position%J
                write (UnitID) CurrentPartic%Position%K
                write (UnitID) CurrentPartic%Position%CellI
                write (UnitID) CurrentPartic%Position%CellJ
                write (UnitID) CurrentPartic%Position%CellK
                write (UnitID) CurrentPartic%Position%X
                write (UnitID) CurrentPartic%Position%Y
                write (UnitID) CurrentPartic%Position%Z

                write (UnitID) CurrentPartic%Geometry%Volume
                write (UnitID) CurrentPartic%Geometry%InitialVolume
                write (UnitID) CurrentPartic%Geometry%VolVar
                write (UnitID) CurrentPartic%Beached

                if (CurrentOrigin%State%Oil) then
                    write (UnitID) CurrentPartic%Geometry%VolumeOil
                endif  

                if (CurrentOrigin%State%Deposition) then
                
                    write (UnitID) CurrentPartic%Deposited
                    write (UnitID) CurrentPartic%TauErosion

                endif

                if (CurrentOrigin%State%Age) then
                
                    write (UnitID) CurrentPartic%Age

                endif

                if (CurrentOrigin%State%PlumeShear) then
                    write (UnitID) CurrentPartic%U
                    write (UnitID) CurrentPartic%V
                endif

                if (CurrentOrigin%State%FarFieldBuoyancy) then
                    write (UnitID) CurrentPartic%W
                endif

                              
                write (UnitID) CurrentPartic%Concentration
                write (UnitID) CurrentPartic%Mass

                write (UnitID) CurrentPartic%TpercursoX
                write (UnitID) CurrentPartic%TpercursoY
                write (UnitID) CurrentPartic%TpercursoZ
                write (UnitID) CurrentPartic%UD_old
                write (UnitID) CurrentPartic%VD_old
                write (UnitID) CurrentPartic%WD_old

                CurrentPartic => CurrentPartic%Next
            enddo

            !Writes Property information
            CurrentProperty => CurrentOrigin%FirstProperty
            do while (associated(CurrentProperty))

                write (UnitID) CurrentProperty%Name

                CurrentProperty => CurrentProperty%Next
            enddo

            !Writes Oil information
            If (CurrentOrigin%State%Oil)                                       &
                call WriteFinalOil(CurrentOrigin%ObjOil, UnitID)

            CurrentOrigin => CurrentOrigin%Next


        enddo CurrOr
            
        !Closes File
        call UnitsManager (UnitID, CLOSE_FILE, STAT = STAT_CALL)


    end subroutine WriteFinalPartic

    !--------------------------------------------------------------------------

    subroutine ReadFinalPartic 

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: NewOrigin    => null()
        type (T_Partic), pointer                    :: NewParticle  => null()
        type (T_Property), pointer                  :: NewProperty  => null()
        integer                                     :: OldOrigins
        integer                                     :: UnitID
        integer                                     :: STAT_CALL
        integer                                     :: nO, nP
        integer                                     :: nParticle, nProperties
        integer                                     :: Dummy

        !Opens File
        call UnitsManager (UnitID, OPEN_FILE, STAT = STAT_CALL)

        open   (unit = UnitID, file = Me%Files%Initial, status = 'old', form = 'unformatted')


        read(UnitID) OldOrigins

        do nO = 1, OldOrigins

            call AllocateNewOrigin (NewOrigin)
            
            NewOrigin%Old   = ON 

            !Reads Origin Information
            read (UnitID) NewOrigin%Name
            read (UnitID) NewOrigin%ID
            read (UnitID) NewOrigin%State%Oil
            read (UnitID) NewOrigin%State%Deposition
            read (UnitID) NewOrigin%State%Age
            read (UnitID) NewOrigin%State%ComputePlume
            read (UnitID) NewOrigin%State%PlumeShear
            read (UnitID) NewOrigin%State%FarFieldBuoyancy
            read (UnitID) nParticle
            read (UnitID) nProperties

            if (NewOrigin%State%Oil) then
        
                !Allocates GridThickness
                allocate (NewOrigin%GridThickness(Me%ExternalVar%Size%ILB:            &
                                                  Me%ExternalVar%Size%IUB,            &
                                                  Me%ExternalVar%Size%JLB:            &
                                                  Me%ExternalVar%Size%JUB),           &
                                                  STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ReadFinalPartic - ModuleLagrangian - ERR04'

                read (UnitID) NewOrigin%GridThickness
            
                !Allocates GridThickness
                allocate (NewOrigin%OilGridConcentration(Me%ExternalVar%Size%ILB:            &
                                                         Me%ExternalVar%Size%IUB,            &
                                                         Me%ExternalVar%Size%JLB:            &
                                                         Me%ExternalVar%Size%JUB),           &
                                                         STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ReadFinalPartic - ModuleLagrangian - ERR04a'

                read (UnitID) NewOrigin%OilGridConcentration

            end if

            !Reads Particle Information
            do nP = 1, nParticle

                call AllocateNewParticle (NewParticle, nProperties, NewOrigin%NextParticID)

                read (UnitID) Dummy !Does not read NewParticle%ID any more
                                    !This is due to an error in attributing 
                                    !identical ID's to different particles
                                    !This error was not completely understood
                                    !but the problem was solved. If you find 
                                    !any problem regarding memory errors and 
                                    !particle ID's please warn Luis or Frank

                read (UnitID) NewParticle%Position%I
                read (UnitID) NewParticle%Position%J
                read (UnitID) NewParticle%Position%K
                read (UnitID) NewParticle%Position%CellI
                read (UnitID) NewParticle%Position%CellJ
                read (UnitID) NewParticle%Position%CellK
                read (UnitID) NewParticle%Position%X
                read (UnitID) NewParticle%Position%Y
                read (UnitID) NewParticle%Position%Z

                read (UnitID) NewParticle%Geometry%Volume
                read (UnitID) NewParticle%Geometry%InitialVolume
                read (UnitID) NewParticle%Geometry%VolVar
                read (UnitID) NewParticle%Beached

                if (NewOrigin%State%Oil) then
                    read (UnitID) NewParticle%Geometry%VolumeOil
                endif  

                if (NewOrigin%State%Deposition) then
                
                    read (UnitID) NewParticle%Deposited
                    read (UnitID) NewParticle%TauErosion

                endif

                if (NewOrigin%State%Age) then
                
                    read (UnitID) NewParticle%Age

                endif

                if (NewOrigin%State%PlumeShear) then
                    read (UnitID) NewParticle%U
                    read (UnitID) NewParticle%V
                endif

                if (NewOrigin%State%FarFieldBuoyancy) then
                    read (UnitID) NewParticle%W
                endif

                read (UnitID) NewParticle%Concentration
                read (UnitID) NewParticle%Mass

                read (UnitID) NewParticle%TpercursoX
                read (UnitID) NewParticle%TpercursoY
                read (UnitID) NewParticle%TpercursoZ
                read (UnitID) NewParticle%UD_old
                read (UnitID) NewParticle%VD_old
                read (UnitID) NewParticle%WD_old

                call InsertParticleToList (NewOrigin, NewParticle, .false.)

            enddo

            !Reads Property Information
            do nP = 1, nProperties

                call AllocateNewProperty (NewProperty)
                read (UnitID) NewProperty%Name

                call InsertPropertyToList(NewOrigin, NewProperty, SetStates = .false.)

            enddo

            !Reads OilInformation
            if (NewOrigin%State%Oil) then
                

                !Starts Oil
                call StartOil(OilID             = NewOrigin%ObjOil,                      &     
                              TimeID            = Me%ObjTime,                            &     
                              EnterDataID       = Me%ObjEnterData,                       &
                              HorizontalGridID  = Me%ObjHorizontalGrid,                  &     
                              GeometryID        = Me%ObjGeometry,                        &     
                              MapID             = Me%ObjMap,                             &     
                              DT                = Me%DT_PARTIC,                          &     
                              ContCalc          = NewOrigin%Old,                         &
                              ExtractType       = FromBlockInBlock,                      &                   
                              STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFinalPartic - ModuleLagrangian - ERR01'

                call ReadFinalOil(NewOrigin%ObjOil, UnitID)

            end if


            call InsertOriginToList (Me%FirstOldOrigin, NewOrigin, Me%nOldOrigins)

        enddo

        !Closes File
        call UnitsManager (UnitID, CLOSE_FILE, STAT = STAT_CALL)


    end subroutine ReadFinalPartic

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine KillLight ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        integer                                     :: iGroup, STAT_CALL

        !Begin-----------------------------------------------------------------

        !Updates light
        do iGroup = 1, Me%nGroups

            call KillLightExtinction(Me%Light%ObjLightExtinction(iGroup), STAT = STAT_CALL)

        enddo

        deallocate (Me%Light%ObjLightExtinction)
        deallocate (Me%Light%TopRadiationCells)

    end subroutine KillLight

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine KillParticleStatistic 

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty  => null()
        integer                                     :: iProp, STAT_CALL

        CurrentProperty => Me%FirstOrigin%FirstProperty
        do while (associated(CurrentProperty))

            iProp = 1
            if (CurrentProperty%Statistics) then

                if (CurrentProperty%StatisticsLag) then
                
                    call WriteFrequencyLag(CurrentProperty)
                
                    deallocate(CurrentProperty%FrequencyLag)

                endif

                call KillStatistic      (CurrentProperty%Statistic1_ID, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillParticleStatistic - ModuleLagrangian - ERR01'

                if (Me%OutPut%ConcMaxTracer) then

                    call KillStatistic      (CurrentProperty%Statistic2_ID, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillParticleStatistic - ModuleLagrangian - ERR02'

                endif 

            endif
            CurrentProperty => CurrentProperty%Next
            iProp = iProp + 1
        enddo

    end subroutine KillParticleStatistic

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (ObjLagrangian_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjLagrangian_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjLagrangian_ID > 0) then
            call LocateObjLagrangian(ObjLagrangian_ID)
            ready_ = VerifyReadLock (mLAGRANGIAN_,  Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjLagrangian (LagrangianID)

        !Arguments-------------------------------------------------------------
   
        integer                                     :: LagrangianID

        !Local-----------------------------------------------------------------

        Me => FirstLagrangian
        do while (associated (Me))
            if (Me%InstanceID == LagrangianID) exit
            Me => Me%Next
        enddo

        if (.not. associated (Me))                                          &
            stop 'ModuleLagrangian - LocateObjLagrangian - ERR01'

    end subroutine LocateObjLagrangian

    !--------------------------------------------------------------------------

    subroutine ReadLockHorizontalGrid

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, GEOG, UTM, MIL_PORT
        integer                                     :: SIMPLE_GEOG, GRID_COORD, NLRD, CoordType 
        
        !----------------------------------------------------------------------
        

        !XX, YY
        call GetHorizontalGrid (Me%ObjHorizontalGrid, XX = Me%ExternalVar%XX,       &
                                YY = Me%ExternalVar%YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangian - ERR20'

        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,         &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,  &
                               NLRD = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(Me%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangian - ERR30'


        if(CoordType == UTM .or. CoordType == MIL_PORT .or.                         &
           CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(Me%ObjHorizontalGrid,                            &
                                   XX_IE = Me%ExternalVar%XX_IE,                    &
                                   YY_IE = Me%ExternalVar%YY_IE,                    &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangian - ERR40'

        else

            call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                     &
                                          GridLatitudeConn  = Me%ExternalVar%YY_IE, &
                                          GridLongitudeConn = Me%ExternalVar%XX_IE, &
                                          STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangian - ERR50'

        end if


    end subroutine ReadLockHorizontalGrid

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------
        
        !Gets Time
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)              
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR10'

        Me%ExternalVar%RunPeriod = Me%ExternalVar%Now - Me%ExternalVar%BeginTime

        call ReadLockHorizontalGrid


i1:     if (.not. Me%RunOnlyMov2D) then

            !Gets Bathymetry
            call GetGridData          (Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT = STAT_CALL)     
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR60'


            !Gets ExteriorPoints 2D
            call GetBoundaries      (Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D, &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR70'


            !WaterColumn
            call GetGeometryWaterColumn(Me%ObjGeometry, Me%ExternalVar%WaterColumn, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR80'


            !SZZ, DWZ    
            call GetGeometryDistances(Me%ObjGeometry,                                   & 
                                      SZZ         = Me%ExternalVar%SZZ,                 &
                                      ZCellCenter = Me%ExternalVar%ZCellCenter,         &
                                      DWZ         = Me%ExternalVar%DWZ,                 &
                                      DWZ_Xgrad   = Me%ExternalVar%DWZ_Xgrad,           &
                                      DWZ_Ygrad   = Me%ExternalVar%DWZ_Ygrad,           &
                                      STAT        = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR90'

            !VolumeZ
            call GetGeometryVolumes(Me%ObjGeometry,                                     &
                                    VolumeZ = Me%ExternalVar%VolumeZ,                   &
                                    STAT    = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR100'

            !kFloorZ
            call GetGeometryKFloor (Me%ObjGeometry,                                     &
                                    Z       = Me%ExternalVar%kFloor,                    &
                                    STAT    = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR110'


            !WaterPoints3D
            call GetWaterPoints3D(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR120'


            !LandPoints3D
            call GetLandPoints3D(Me%ObjMap, Me%ExternalVar%LandPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR130'


            !OpenPoints3D
            call GetOpenPoints3D(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR140'


            !Compute faces
            call GetComputeFaces3D(Me%ObjMap,                                           &
                                   ComputeFacesU3D = Me%ExternalVar%ComputeFaces3D_U,   &
                                   ComputeFacesV3D = Me%ExternalVar%ComputeFaces3D_V,   &
                                   STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR150'


            !Lupward, Ldownward
            call GetMixingLengthVertical(Me%ObjTurbulence,                              &
                                         Lupward   = Me%ExternalVar%Lupward,            &
                                         Ldownward = Me%ExternalVar%Ldownward,          &
                                         STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR160'


            !MixingLengthX, MixingLengthY
            call GetMixingLengthHorizontal(Me%ObjTurbulence,                            &
                                           MixingLengthX = Me%ExternalVar%MixingLengthX,&
                                           MixingLengthY = Me%ExternalVar%MixingLengthY,&
                                           STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR170'

        
            !Velocity_U, Velocity_V
            call GetHorizontalVelocity(Me%ObjHydrodynamic,                              &
                                       Velocity_U = Me%ExternalVar%Velocity_U,          &
                                       Velocity_V = Me%ExternalVar%Velocity_V,          &
                                       STAT       = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR180'

        
            !Velocity_W
            call GetVerticalVelocity(Me%ObjHydrodynamic,                                &
                                     Velocity_W      = Me%ExternalVar%Velocity_W,       &
                                     STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangian - ERR190'

        else i1

        !Allocate 3D matrixes needed

        endif i1


        !----------------------------------------------------------------------

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadLockEulerianDensity

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------

        !Density
        if (Me%State%FarFieldBuoyancy  .or. Me%State%ComputePlume) then
            call GetSigma       (Me%ObjWaterProperties,                            &
                                 Sigma             = Me%ExternalVar%SigmaDensity,  &
                                 CurrentTime       = Me%ExternalVar%Now,           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockEulerianDensity - ModuleLagrangian - ERR10'
        else
            nullify(Me%ExternalVar%SigmaDensity)
        endif
    
    
        if (Me%State%Oil) then
            call GetDensity     (Me%ObjWaterProperties,                       &
                                 Density           = Me%ExternalVar%Density,  &
                                 CurrentTime       = Me%ExternalVar%Now,      &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockEulerianDensity - ModuleLagrangian - ERR20'
        else
            nullify(Me%ExternalVar%Density)
        endif

    end subroutine ReadLockEulerianDensity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine ReadUnLockHorizontalGrid
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------

        !XX, YY
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, Me%ExternalVar%XX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockHorizontalGrid - ModuleLagrangian - ERR10'

        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, Me%ExternalVar%YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockHorizontalGrid - ModuleLagrangian - ERR20'

        !XX_IE, YY_IE
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, Me%ExternalVar%XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR30'

        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, Me%ExternalVar%YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR40'



    end subroutine ReadUnLockHorizontalGrid

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL


        call ReadUnLockHorizontalGrid


i1:     if (.not. Me%RunOnlyMov2D) then

            !Gets Bathymetry
            call UnGetGridData (Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR50'


            !Gets ExteriorPoints 2D
            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D,   &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR60'


            !WaterColumn
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%WaterColumn, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR70'

            !SZZ
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%SZZ, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR80'

            !ZCellCenter
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%ZCellCenter, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR90'

            !DWZ    
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%DWZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR100'

            !VolumeZ
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%VolumeZ, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR110'

            !kFloorZ
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%kFloor, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR120'


            !DWZ_Xgrad
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%DWZ_Xgrad, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR122'

            !DWZ_Ygrad
            call UnGetGeometry (Me%ObjGeometry, Me%ExternalVar%DWZ_Ygrad, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR124'

            !WaterPoints3D
            call UnGetMap      (Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR130'


            !LandPoints3D
            call UnGetMap      (Me%ObjMap, Me%ExternalVar%LandPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR140'


            !OpenPoints3D
            call UnGetMap      (Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR150'

            !Compute faces U
            call UnGetMap      (Me%ObjMap, Me%ExternalVar%ComputeFaces3D_U, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR160'

            !Compute faces V
            call UnGetMap      (Me%ObjMap, Me%ExternalVar%ComputeFaces3D_V, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR170'


            !Lupward
            call UngetTurbulence(Me%ObjTurbulence, Me%ExternalVar%Lupward, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR180'

            !Ldownward
            call UngetTurbulence(Me%ObjTurbulence, Me%ExternalVar%Ldownward, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR190'


            !MixingLengthX
            call UngetTurbulence(Me%ObjTurbulence, Me%ExternalVar%MixingLengthX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR200'


            !MixingLengthY
            call UngetTurbulence(Me%ObjTurbulence, Me%ExternalVar%MixingLengthY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR210'
        

            !Velocity_U
            call UngetHydrodynamic (Me%ObjHydrodynamic, Me%ExternalVar%Velocity_U, STAT  = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR220'


            !Velocity_V
            call UngetHydrodynamic (Me%ObjHydrodynamic, Me%ExternalVar%Velocity_V, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR230'


            !Velocity_W
            call UngetHydrodynamic (Me%ObjHydrodynamic, Me%ExternalVar%Velocity_W, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangian - ERR240'

        else i1

        !Allocate 3D matrixes needed

        endif i1
        !----------------------------------------------------------------------

    end subroutine ReadUnLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockEulerianDensity

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------

        !Sigma Density
        if ((Me%State%FarFieldBuoyancy  .or. Me%State%ComputePlume) .and.               &
            associated(Me%ExternalVar%SigmaDensity)) then

            call UngetWaterProperties (Me%ObjWaterProperties, Me%ExternalVar%SigmaDensity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockEulerianDensity - ModuleLagrangian - ERR10'

        endif

        !Density
        if ((Me%State%Oil) .and.  associated(Me%ExternalVar%Density)) then
            call UngetWaterProperties (Me%ObjWaterProperties, Me%ExternalVar%Density, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockEulerianDensity - ModuleLagrangian - ERR20'
        endif

    end subroutine ReadUnLockEulerianDensity

    !--------------------------------------------------------------------------

end Module ModuleLagrangian

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

