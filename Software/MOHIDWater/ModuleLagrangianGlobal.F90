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

Module ModuleLagrangianGlobal

!BOP
!
! !MODULE: ModuleLagrangianGlobal

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
!   KILL_PART_INSIDE_BOX    : logical                   [0]                     !Kill all particles inside the emission boxes. If  
                                                                                !it is an instantaneous emission.                 

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
                                                                                !a variable concentration

!   NOWQM                   : logical                   [0]                     ! To compute age without running moduleWQM
!   MIN_CONCENTRATION       : real                      []
!   AMBIENT_CONC            : real                      [0.0]                   !Ambient concentration
!   T90_VARIABLE            : logical                   [0]                     !Check if the user wants to compute T90 function 
!                                                                                of ambient properties: salinity,temperature,light
!   T90_VAR_METHOD          : integer                   [1]                     !Fecal decay according to Canteras et al. (1995)
!                                                       [2]                     !Fecal decay according to Chapra (1997)
!   T90                     : real                      [7200]                  !Coliform Decay rate (s)
!   T90_NAME                ; char                      [T90]
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
                                       GetDataOnlineString, SetMatrixValue, TimeToString,   &
                                       ChangeSuffix, ConstructPropertyID
    use ModuleEnterData,        only : ReadFileName, ConstructEnterData, GetData,           &
                                       ExtractBlockFromBuffer, ExtractBlockFromBlock,       &
                                       Block_Unlock, GetOutPutTime, RewindBuffer,           &
                                       GetKeywordFromLine, GetFullBufferLine,               &
                                       ReplaceFullBufferLine, RewindBlock,                  &
                                       GetBlockSize, KillEnterData
    use ModuleWaterQuality,     only : StartWaterQuality, WaterQuality, GetDTWQM,           &
                                       GetWQPropIndex, KillWaterQuality
    use ModuleGridData,         only : GetGridData, GetMaximumValue, UngetGridData
    use ModuleTimeSerie,        only : StartTimeSerie, StartTimeSerieInput, WriteTimeSerie, &
                                       GetNumberOfTimeSeries, GetTimeSerieLocation,         &
                                       CorrectsCellsTimeSerie, GetTimeSerieIntegral,        &
                                       WriteTimeSerieLine, GetTimeSerieValue, KillTimeSerie,&
                                       TryIgnoreTimeSerie
    use ModuleLightExtinction,  only : ConstructLightExtinction, ModifyLightExtinctionField,&
                                       GetLightExtinctionOptions, KillLightExtinction,      &
                                       GetShortWaveExtinctionField, UnGetLightExtinction,   &
                                       GetLongWaveExtinctionCoef, GetRadiationPercentages
    use ModuleHorizontalMap,    only : GetBoundaries, UnGetHorizontalMap
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, WriteHorizontalGrid,              &
                                       UnGetHorizontalGrid, GetGridCoordType, GetCoordTypeList,&
                                       LocateCell, GetDefineCellsMap, GetGridLatitudeLongitude,&
                                       GetXYInsideDomain, GetXYCellZ, GetCellZ_XY,          &
                                       GetLatitudeLongitude, GetGridCellArea,               &
                                       GetGridBorderType
    use ModuleAssimilation,     only : StartAssimilation, GetAssimilationField,             &
                                       UnGetAssimilation, KillAssimilation
    use ModuleGeometry,         only : GetGeometrySize, GetGeometryWaterColumn,             &
                                       GetGeometryDistances, GetGeometryVolumes,            &
                                       GetGeometryKFloor, UnGetGeometry, GetLayer4Level
    use ModuleMap,              only : GetWaterPoints3D, GetLandPoints3D, GetOpenPoints3D,  &
                                       GetComputeFaces3D, UngetMap             
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes, UngetBoxDif,&
                                       BoxDif, CheckIfInsideBox, GetIfBoxInsideDomain, KillBoxDif        
    use ModuleTurbulence,       only : GetMixingLengthVertical, GetMixingLengthHorizontal,  &
                                       UngetTurbulence       
    use ModuleHydrodynamic,     only : StartHydrodynamic,                                   &
                                       GetHorizontalVelocity,                               &
                                       GetVerticalVelocity, UngetHydrodynamic
    use ModuleWaterProperties,  only : Construct_WaterProperties, WaterPropertyExists,      &
                                       GetConcentration, GetDensity, GetSigma,              &
                                       UngetWaterProperties, GetSPM,                        &
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
    public  :: AllocateLagrangianGlobal
    public  :: ConstructLagrangianGlobal
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
    public  :: ModifyLagrangianGlobal
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
    private :: Convert_Z_CellK
    private :: Convert_CellK_Z
    private :: Convert_CellK_K

    private :: Locate_ModelDomain

    private :: Search_Property

    !Selector
    public  :: SetLagrangianShearGlobal
    public  :: SetLagrangianAtmPressureGlobal
    public  :: SetLagrangianWindGlobal
    public  :: GetLagrangianAirOptionsGlobal
    public  :: SetLagSolarRadiationGlobal

    !Destructor
    public  :: DeallocateLagrangianGlobal
    public  :: KillLagrangianGlobal
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
    integer, parameter                          :: FromTimeSerie        = 3


    !Online Emission options 
    integer, parameter                          :: ParticleOne          = 1
    integer, parameter                          :: Particle100          = 2

    character(LEN = StringLength), parameter    :: block_begin          = '<BeginOrigin>'
    character(LEN = StringLength), parameter    :: block_end            = '<EndOrigin>'
    character(LEN = StringLength), parameter    :: statistic_begin      = '<BeginStatistic>'
    character(LEN = StringLength), parameter    :: statistic_end        = '<EndStatistic>'
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
!----------------------------------------------------------------------------
!                                   Eulerian types
!----------------------------------------------------------------------------

    type T_OverLay
        real, dimension(:, :, :), pointer       :: VelUFinal
        real, dimension(:, :, :), pointer       :: VelVFinal
    end type T_OverLay

    type T_ParticleGrid
        integer                                 :: CoordType
        logical                                 :: HaveLatLongGrid = .false.
        logical                                 :: GeoGrid         = .false.
        real                                    :: LatDefault, LongDefault
        real, dimension(:, :), pointer          :: ParticXX
        real, dimension(:, :), pointer          :: ParticYY
    end type T_ParticleGrid

    type T_Light
        integer                                 :: ObjLightExtinction
                          !i,j,k
        real,    dimension(:,:,:), pointer      :: TopRadiationCells, ShortWaveExtinctionField
        logical                                 :: Compute      = OFF
    end type T_Light

    type T_EulerianMonitor
        real(8), dimension(:, :, :), pointer    :: Mass
    end type T_EulerianMonitor

    type T_Monitorization
        real(8), dimension(:),    pointer       :: InstBoxVolume
        real(8), dimension(:),    pointer       :: InstBoxMass
        real(8), dimension(:, :), pointer       :: InstMassByOrigin
        real(8), dimension(:, :), pointer       :: InstVolumeByOrigin
        real(8), dimension(:),    pointer       :: IntgBoxVolume
        real(8), dimension(:, :), pointer       :: IntgVolumeByOrigin
        integer, dimension(:),    pointer       :: NumberOfCellsPerBox
        integer, dimension(:, :), pointer       :: NumberOfCellsFromOrigin
        integer, dimension(:), pointer          :: ObjTimeSerie
                           !i,j,k
        integer, dimension(:,:,:), pointer      :: Boxes
        integer                                 :: NumberOfBoxes
    end type T_Monitorization


    type T_Lag2Euler
                          !i, j, k, p, ig
        real,    dimension(:, :, :, :, :), pointer          :: GridConc
        integer, dimension(:, :, :, :   ), pointer          :: GridTracerNumber
        real,    dimension(:, :, :, :, :), pointer          :: GridMaxTracer
        real,    dimension(:, :,    :, :), pointer          :: GridBottomConc

        real(8), dimension(:, :, :, :   ), pointer          :: GridVolume
        real,    dimension(:, :, :, :, :), pointer          :: GridMass
        real,    dimension(:, :, :, :   ), pointer          :: GridBottomMass
                          !p, ig  
        real,    dimension(:, :),          pointer          :: MeanConc, AmbientConc
        real,    dimension(:, :),          pointer          :: MinConc, MassVolCel

        !Sediments
        real,    dimension(:, :, :),       pointer          :: TauErosionGrid
        real,    dimension(:, :, :),       pointer          :: MassSedGrid

    end type T_Lag2Euler

    type T_PropStatistic
                          !i,j,k,f,ig
        real,    dimension(:,:,:,:,:), pointer :: FrequencyLag
                                  !ig
        integer, dimension(        :), pointer :: Statistic1_ID, Statistic2_ID
    end type T_PropStatistic

    type T_OilSpreading
        !Oil vectors
        real,    pointer, dimension(:,:  )      :: VelocityX
        real,    pointer, dimension(:,:  )      :: VelocityY

        !Oil
                          !i, j 
        real,    dimension(:, :),       pointer          :: GridThickness
        real,    dimension(:, :),       pointer          :: OilGridConcentration

        logical, dimension(:,:), pointer        :: AreaFlag

    end type T_OilSpreading


    type     T_EulerModel

        character(len=StringLength)             :: Name

        !ObjBathymetry
        real,    dimension(:, : ), pointer      :: Bathymetry
        real,    dimension(:    ), pointer      :: XX
        real,    dimension(:    ), pointer      :: YY

        real,    dimension(:, : ), pointer      :: XX_IE
        real,    dimension(:, : ), pointer      :: YY_IE
        real,    dimension(:, : ), pointer      :: DZX, DZY

        real,    dimension(:, : ), pointer      :: GridCellArea

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



        !ObjAssimilation
        real,    pointer, dimension(:,:,:)      :: OverLayU
        real,    pointer, dimension(:,:,:)      :: OverLayV

        !ObjWaterProperties
        real,    pointer, dimension(:,:,:)      :: Density
        real,    pointer, dimension(:,:,:)      :: SigmaDensity
        real,    pointer, dimension(:,:,:)      :: Temperature3D
        real,    pointer, dimension(:,:,:)      :: Salinity3D
        real,    pointer, dimension(:,:,:)      :: FishFood3D
        real,    pointer, dimension(:,:,:)      :: SPM3D


        real   , dimension(:,:,:),    pointer   :: BeachingProbability

        integer                                 :: ObjTimeSerie         = 0
        integer                                 :: ObjWaterProperties   = 0
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


        type(T_ParticleGrid)                    :: Grid
        type(T_EulerianMonitor)                 :: EulerianMonitor

        type(T_Monitorization )                 :: Monitor
                          !i, j, k  
        real,    dimension(:, :, :   ), pointer :: RelativeMassFilter, MassFiltered

        type(T_OverLay        )                 :: Overlay
                                !ig
        type(T_Light), dimension(:),pointer     :: Light
        type(T_Lag2Euler      )                 :: Lag2Euler
                                          !p  
        type(T_PropStatistic  ), dimension(:),pointer :: PropStatistic
                                          !ig  
        type(T_OilSpreading   ), dimension(:),pointer :: OilSpreading

    end type T_EulerModel
!----------------------------------------------------------------------------
!                                   Lagrangian types
!----------------------------------------------------------------------------

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
        logical                                 :: T90Compute           = OFF
        logical                                 :: Filtration           = OFF
        logical                                 :: Boxdif               = OFF
        logical                                 :: OriginDefault        = OFF
        logical                                 :: T90                  = OFF
        logical                                 :: KillPartInsideBox    = OFF
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
        character(PathLength)                   :: Nomfich    
    end type T_Files


    type T_Online
        real,  dimension(:,:), pointer          :: StartDate
        real,  dimension(:),   pointer          :: X, Y 
        real,  dimension(:),   pointer          :: WindCoef
        character(Len=23)                       :: TimeStamp
    end type T_Online

    !Output
    type       T_OutPut
         type (T_Time), pointer, dimension(:)   :: OutTime, RestartOutTime
         integer                                :: NextOutPut, NextRestartOutPut
         integer                                :: OutPutConcType = Maximum
         logical                                :: ConcMaxTracer   
         logical                                :: OriginEnvelope
         logical                                :: Write_, WriteRestartFile = .false. 
         logical                                :: RestartOverwrite
         logical                                :: DummyParticleStartDate
    end type T_OutPut

  
    !Defines a generic Position
    type T_Position
        integer                                 :: I                        = null_int
        integer                                 :: J                        = null_int
        integer                                 :: K                        = null_int
        real                                    :: X                        = null_real
        real                                    :: Y                        = null_real
        real                                    :: CartX                    = null_real
        real                                    :: CartY                    = null_real
        real                                    :: CoordX                   = null_real
        real                                    :: CoordY                   = null_real
        real                                    :: Z                        = null_real
        real                                    :: CellI                    = null_real
        real                                    :: CellJ                    = null_real
        real                                    :: CellK                    = null_real
        integer                                 :: DepthDefinition          = null_int
        logical                                 :: MaintainRelative         = .false.
        real                                    :: Depth                    = null_real
        !Each particle has the eulerian grid index 
        integer                                 :: ModelID                  = null_int
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
        logical                                 :: T90Compute               = OFF
        integer                                 :: T90Var_Method            = null_int
        real                                    :: T90                      = 7200.     !Coliform Bacteria decay time
        character(StringLength)                 :: T90Name
        logical                                 :: T90ON                    = OFF
        character(PathLength)                   :: T90File
        integer                                 :: T90Column                = null_int
        integer                                 :: TimeSerieT90             = 0
        type(T_Partition)                       :: SedimentPartition
        type(T_Partition)                       :: WaterPartition
        logical                                 :: ConcVariable
        integer                                 :: ConcColumn
        logical                                 :: NoWQM                    =.false.
        type (T_Property), pointer              :: Next
        real                                    :: ExtinctionParameter
        logical                                 :: Filtration               = OFF
        logical                                 :: WritesTimeSerie          = OFF
        logical                                 :: WritesPropHDF            = OFF
        real                                    :: MinValue, MaxValue
        logical                                 :: MinON, MaxON
        character(PathLength)                   :: DischargeFile
        integer                                 :: TimeSerieInput           = 0
        logical                                 :: EqualToAmbient           = OFF
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
        logical                                 :: EmissionON
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
        integer                                 :: TimeSerieInputFlow           = 0
        type (T_Time)                           :: InstantEmission
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
        integer                                 :: nPropT90
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
        logical                                 :: Beaching          = OFF
        logical                                 :: Filtration        = OFF
        logical                                 :: Default           = OFF
        logical                                 :: OnlyOnceEmit      = OFF
        logical                                 :: KillPartInsideBox = OFF
        logical                                 :: EstimateMinVol    = OFF
        integer                                 :: MaxPart
    end type T_Origin

    type T_OptionsStat

        type(T_PropertyID)                      :: ID
        integer                                 :: PropOrder
        
        character(len=PathLength)               :: File
        integer                                 :: nClassesLag
        real, dimension(:,:    ), pointer       :: ClassesLag
        logical                                 :: Lag
    end type T_OptionsStat

    type T_Statistic
        integer                                    :: PropNumber
        type(T_OptionsStat), dimension(:), pointer :: OptionsStat
        logical            , dimension(:), pointer :: ON
    end type T_Statistic



    !ExternalVar
    type T_ExternalVar

        !ObjOil
        !real,    pointer, dimension(:,:  )      :: SpreadingVelocityX
        !real,    pointer, dimension(:,:  )      :: SpreadingVelocityY
        real                                    :: DiffVelocity
        real                                    :: VWaterContent
        real                                    :: MWaterContent
        real                                    :: AreaTotal
        real                                    :: OilDensity
        real                                    :: OilViscosity
        real                                    :: FMDispersed
        real                                    :: FMEvaporated
        real                                    :: MDispersed
        integer                                 :: ThicknessGradient, Fay, SpreadingMethod

        !Time - by default is used the time object of the model with higher priority 
        integer                                 :: ObjTime              = 0
        type(T_Time)                            :: Now
        type(T_Time     )                       :: BeginTime
        type(T_Time     )                       :: EndTime
        type(T_Time     )                       :: LastConcCompute

        real                                    :: RunPeriod

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
        logical                                 :: IgnoreON
        integer, dimension(:    ),  pointer     :: GroupIDs
        integer, dimension(:    ),  pointer     :: nOriginsGroup

        character(StringLength)                 :: MonitorProperty
        integer                                 :: MonitorPropertyID           = null_int

        real                                    :: DT_Partic
        real                                    :: BeachingLimit
        real                                    :: DefaultBeachingProbability

        type(T_Statistic)                       :: Statistic

        type(T_Time)                            :: NextCompute
        type(T_Time)                            :: Now

        type (T_Time)                           :: NextFiltration

        type(T_Origin     ), pointer            :: FirstOrigin          => null()
        type(T_Origin     ), pointer            :: FirstOldOrigin       => null()
        type(T_Origin     ), pointer            :: OriginDefault        => null()

        type(T_ExternalVar)                     :: ExternalVar

        type(T_EulerModel ), pointer, dimension(:) :: EulerModel
        integer                                 :: EulerModelNumber

!#ifdef _CGI_
        type(T_Online)                          :: Online    
!#endif      


        logical                                 :: WritesTimeSerie      = .false.
        logical                                 :: RunOnlyMov2D         = .false.
        logical                                 :: Overlay
        logical                                 :: FirstIteration       = .true.

        integer, dimension(:), pointer          :: ObjHDF5              

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

    subroutine AllocateLagrangianGlobal(LagrangianID, STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        integer, optional, intent(OUT)              :: STAT


        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_


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
            
            !Returns ID
            LagrangianID    = Me%InstanceID
            STAT_           = SUCCESS_

        else

            stop 'AllocateLagrangianGlobal - ModuleLagrangianGlobal - ERR10' 

        endif


        if (present(STAT)) STAT = STAT_       
        
    end subroutine AllocateLagrangianGlobal
             
    
    !--------------------------------------------------------------------------

    subroutine ConstructLagrangianGlobal(LagrangianID,                                  &
                                   Nmodels,                                             &
                                   ModelNames,                                          &
                                   FileNomfich,                                         &
                                   LagInstance,                                         &
                                   STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        integer                                     :: Nmodels
        character(len=*), dimension(:),   pointer   :: ModelNames
        character(len=*)                            :: FileNomfich        
        integer,          dimension(:,:), pointer   :: LagInstance
        integer, optional, intent(OUT)              :: STAT


        !Local-----------------------------------------------------------------
        integer, dimension(TotalLagInst_)           :: TimeID,  GridDataID,             &
                                                       HorizontalGridID,                &
                                                       HorizontalMapID,                 &
                                                       GeometryID,                      &
                                                       MapID,                           &
                                                       AssimilationID,                  &
                                                       HydrodynamicID,                  &
                                                       TurbulenceID,                    &
                                                       WavesID,                         &
                                                       WaterPropertiesID

        type (T_EulerModel), pointer                :: EulerModel
        integer                                     :: STAT_, ready_, em, NmodelsFinal
        integer                                     :: STAT_CALL
        integer, dimension(:), pointer              :: IndexMatch

        !Begin--------------------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
em0:        do em =1, Nmodels            
                TimeID           (em) = LagInstance(1,em)
                GridDataID       (em) = LagInstance(2,em)              
                HorizontalGridID (em) = LagInstance(3,em)
                HorizontalMapID  (em) = LagInstance(4,em)
                GeometryID       (em) = LagInstance(5,em)
                MapID            (em) = LagInstance(6,em)
                AssimilationID   (em) = LagInstance(7,em)
                HydrodynamicID   (em) = LagInstance(8,em)  
                TurbulenceID     (em) = LagInstance(9,em) 
                WavesID          (em) = LagInstance(10,em)
                WaterPropertiesID(em) = LagInstance(11,em)
                
            enddo em0
            
            Me%Files%Nomfich  = FileNomfich            

            ! Construct the variable common to all module  
            call ConstructGlobalVariables

            !Construct enter data 
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangianGlobal - ModuleLagrangianGlobal - ERR13'

            allocate(IndexMatch(1:Nmodels))

            call ConstructEulerModelsList(ModelNames, IndexMatch, Nmodels, NmodelsFinal)

            Me%EulerModelNumber = NmodelsFinal
            allocate(Me%EulerModel(Me%EulerModelNumber))

em1:        do em =1, Me%EulerModelNumber

                Me%EulerModel(em)%Name              = ModelNames(IndexMatch(em))

            !External Modules
                Me%EulerModel(em)%ObjTime           = AssociateInstance (mTIME_,           TimeID            (IndexMatch(em)))
                Me%EulerModel(em)%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID        (IndexMatch(em)))
                Me%EulerModel(em)%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID   (IndexMatch(em)))
                Me%EulerModel(em)%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID  (IndexMatch(em)))
                Me%EulerModel(em)%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID        (IndexMatch(em)))
                Me%EulerModel(em)%ObjMap            = AssociateInstance (mMAP_,            MapID             (IndexMatch(em)))
                Me%EulerModel(em)%ObjHydrodynamic   = AssociateInstance (mHYDRODYNAMIC_,   HydrodynamicID    (IndexMatch(em)))
                Me%EulerModel(em)%ObjTurbulence     = AssociateInstance (mTURBULENCE_,     TurbulenceID      (IndexMatch(em)))
                Me%EulerModel(em)%ObjWaterProperties= AssociateInstance (mWATERPROPERTIES_,WaterPropertiesID (IndexMatch(em)))

                if(WavesID(IndexMatch(em)) /= 0)then
                    Me%EulerModel(em)%ObjWaves      = AssociateInstance (mWAVES_,          WavesID(IndexMatch(em)))
                end if

                if (em ==1) then
                    !External Modules ObjTime for the higher priority model
                    Me%ExternalVar%ObjTime          = Me%EulerModel(em)%ObjTime

                    !Gets Time
                    call GetComputeTimeLimits(Me%ExternalVar%ObjTime,                   &
                                              BeginTime = Me%ExternalVar%BeginTime,     &
                                              EndTime   = Me%ExternalVar%EndTime,       &
                                              STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangianGlobal - ModuleLagrangianGlobal - ERR02'

                    ! Actualized the time
                    call GetComputeCurrentTime(Me%ExternalVar%ObjTime,                            &
                                               Me%ExternalVar%Now, STAT = STAT_CALL)                    
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangianGlobal - ModuleLagrangianGlobal - ERR03'
                    
                    call ConstructOuptusFilesNames                    
                
                endif

                !Gets Size
                call GetGeometrySize(Me%EulerModel(em)%ObjGeometry,                                  &
                                     Size     = Me%EulerModel(em)%Size,                  &
                                     WorkSize = Me%EulerModel(em)%WorkSize,              &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangianGlobal - ModuleLagrangianGlobal - ERR04'

            enddo em1

            !Gets Pointer to External modules
            call ReadLockExternalVar()

            call ReadLockEulerianDensity()


em4:        do em =1, Me%EulerModelNumber
                ! Constructs the Particle Grid
                EulerModel => Me%EulerModel(em)
                call ConstructParticleGrid(EulerModel)
                nullify(EulerModel)
            enddo em4 

            !Constructs the origin list 
            call ConstructOrigins

            if (Me%Overlay) then

em2:            do em =1, Me%EulerModelNumber

                    if (AssimilationID(IndexMatch(em)) == 0) then

                        !Associated the Assimilation
                        call StartAssimilation (Me%EulerModel(em)%ObjAssimilation,                       &
                                                Me%EulerModel(em)%ObjTime,                               &
                                                Me%EulerModel(em)%ObjGridData,                           &
                                                Me%EulerModel(em)%ObjHorizontalGrid,                     &
                                                Me%EulerModel(em)%ObjHorizontalMap,                      &
                                                Me%EulerModel(em)%ObjMap,                                &
                                                Me%EulerModel(em)%ObjGeometry,                           &
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangianGlobal - ModuleLagrangianGlobal - ERR08a'
                    
                        AssimilationID = Me%EulerModel(em)%ObjAssimilation

                    else
                    
                        Me%EulerModel(em)%ObjAssimilation = AssociateInstance (mASSIMILATION_, AssimilationID(IndexMatch(em)))

                    endif

                enddo em2
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
            if (Me%State%WQM .or. Me%State%T90Compute) then
                call ConstructParticLightExtinction
            endif

            !Constructs the Time Series
            if (Me%WritesTimeSerie) then
                call ConstructTimeSeries  
            endif

            call ReadParticStatisticOptions   

            !Kills EnterData
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructLagrangianGlobal - ModuleLagrangianGlobal - ERR14'

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
            
            !Allocates all the necessary matrix for integrate lag data in eulerian grids
            call ConstructLag2Euler      

            !Starts the HDF Output
            if (Me%OutPut%Write_) call ConstructHDF5Output      

            !Starts the Statistic
            if (Me%State%Statistics) then
                call NewParticleMass            
                call ConstructParticStatistic   
            endif

            if (Me%State%Oil) call AllocateOil


            !Message to the User
            call ConstructLog             

            deallocate(IndexMatch)


            !Frees Pointer to External modules
            call ReadUnLockExternalVar()

            call ReadUnLockEulerianDensity()

!            nullify(TimeID, GridDataID, HorizontalGridID, HorizontalMapID, GeometryID)
!            nullify(MapID, AssimilationID, HydrodynamicID, TurbulenceID, WavesID, WaterPropertiesID)

            STAT_           = SUCCESS_

        else

            stop 'ConstructLagrangianGlobal - ModuleLagrangianGlobal - ERR99' 

        endif


        if (present(STAT)) STAT = STAT_


    end subroutine ConstructLagrangianGlobal

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

    !------------------------------------------------------------------------------

    subroutine AllocateOil

        !Local---------------------------------------------------------------------
        integer                     :: em, ig, STAT_CALL                          

        !Begin---------------------------------------------------------------------

d1:     do em =1, Me%EulerModelNumber 

            allocate (Me%EulerModel(em)%OilSpreading(1: Me%NGroups), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'AllocateOil - ModuleLagrangianGlobal - ERR10'

d2:         do ig = 1, Me%NGroups

            !Allocates GridThickness
            allocate (Me%EulerModel(em)%OilSpreading(ig)%GridThickness(Me%EulerModel(em)%Size%ILB: &
                                                                       Me%EulerModel(em)%Size%IUB, &
                                                                       Me%EulerModel(em)%Size%JLB: &
                                                                       Me%EulerModel(em)%Size%JUB),&
                                                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'AllocateOil - ModuleLagrangianGlobal - ERR20'


            !Allocates OilGridConcentration
            allocate (Me%EulerModel(em)%OilSpreading(ig)%OilGridConcentration(Me%EulerModel(em)%Size%ILB: &
                                                                              Me%EulerModel(em)%Size%IUB, &
                                                                              Me%EulerModel(em)%Size%JLB: &
                                                                              Me%EulerModel(em)%Size%JUB),&
                                                                              STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'AllocateOil - ModuleLagrangianGlobal - ERR30'

            allocate (Me%EulerModel(em)%OilSpreading(ig)%VelocityX(Me%EulerModel(em)%Size%ILB: &
                                                                   Me%EulerModel(em)%Size%IUB, &
                                                                   Me%EulerModel(em)%Size%JLB: &
                                                                   Me%EulerModel(em)%Size%JUB),&
                                                                   STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'AllocateOil - ModuleLagrangianGlobal - ERR40'

            allocate (Me%EulerModel(em)%OilSpreading(ig)%VelocityY(Me%EulerModel(em)%Size%ILB: &
                                                                   Me%EulerModel(em)%Size%IUB, &
                                                                   Me%EulerModel(em)%Size%JLB: &
                                                                   Me%EulerModel(em)%Size%JUB),&
                                                                   STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'AllocateOil - ModuleLagrangianGlobal - ERR50'


            allocate (Me%EulerModel(em)%OilSpreading(ig)%AreaFlag (Me%EulerModel(em)%Size%ILB: &
                                                                   Me%EulerModel(em)%Size%IUB, &
                                                                   Me%EulerModel(em)%Size%JLB: &
                                                                   Me%EulerModel(em)%Size%JUB),&
                                                                   STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'AllocateOil - ModuleLagrangianGlobal - ERR60'

            Me%EulerModel(em)%OilSpreading(ig)%AreaFlag(:,:) = .true. 


            enddo d2

        enddo d1

    end subroutine AllocateOil

    !------------------------------------------------------------------------------

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
 
        !Input data file
        Message   ='ASCII file used to construct lagrangian module'
        call ReadFileName('PARTIC_DATA', Me%Files%ConstructData,                        &
                           Message = Message, FilesInput = Me%Files%Nomfich,            &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleLagrangianGlobal - ERR05'



cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_  ) then
            write(*,*)  
            write(*,*) 'Initial file not found.'
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleLagrangianGlobal - ERR08'

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then
            write(*,*)  
            write(*,*) 'Keyword for the initial file not found in nomfich.dat.'
            write(*,*) 'ConstructGlobalVariables - ModuleLagrangianGlobal - WRN01'
            write(*,*)  

        else if (STAT_CALL .EQ. SUCCESS_             ) then
            continue
        else
            stop 'ConstructGlobalVariables - ModuleLagrangianGlobal - ERR09'
        end if cd1  

        !----------------------------------------------------------------------

    end subroutine ConstructGlobalVariables

    !--------------------------------------------------------------------------


    !------------------------------------------------------------------------------

    subroutine ConstructOuptusFilesNames

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(len = StringLength)               :: Message

        !----------------------------------------------------------------------

        !Transient HDF File
        Message   ='Instant fields of lagrangian particles in HDF format.'
        call ReadFileName('PARTIC_HDF', Me%Files%TransientHDF,                          &
                           Message = Message,                                           &
                           TIME_END = Me%ExternalVar%EndTime,                           &
                           Extension = 'hdf',                                           &
                           FilesInput = Me%Files%Nomfich,                               & 
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOuptusFilesNames - ModuleLagrangianGlobal - ERR06'


        !Final file
        Message   ='Final Particle Conditions.'
        call ReadFileName('PARTIC_FIN', Me%Files%Final,                                 &
                           Message = Message,                                           &
                           TIME_END = Me%ExternalVar%EndTime,                           &
                           Extension = 'ptf',                                           &
                           FilesInput = Me%Files%Nomfich,                               &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOuptusFilesNames - ModuleLagrangianGlobal - ERR07'


        !Initial Conditions
        Message   ='Initial Particle Conditions.'
        call ReadFileName('PARTIC_INI', Me%Files%Initial,                               &
                           Message = Message, TIME_END = Me%ExternalVar%Now,            &
                           Extension = 'ptf',                                           &
                           FilesInput = Me%Files%Nomfich,                               &  
                           STAT = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_  ) then
            write(*,*)  
            write(*,*) 'Initial file not found.'
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOuptusFilesNames - ModuleLagrangianGlobal - ERR08'

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then
            write(*,*)  
            write(*,*) 'Keyword for the initial file not found in nomfich.dat.'
            write(*,*) 'ConstructOuptusFilesNames - ModuleLagrangianGlobal - WRN01'
            write(*,*)  

        else if (STAT_CALL .EQ. SUCCESS_             ) then
            continue
        else
            stop 'ConstructOuptusFilesNames - ModuleLagrangianGlobal - ERR09'
        end if cd1  

        !----------------------------------------------------------------------

    end subroutine ConstructOuptusFilesNames

    !--------------------------------------------------------------------------
    !------------------------------------------------------------------------------

    subroutine ConstructEulerModelsList(ModelNames, IndexMatch, Nmodels, NmodelsFinal)

        !Arguments-------------------------------------------------------------
        character(len = *), dimension(:), pointer   :: ModelNames
        integer,            dimension(:), pointer   :: IndexMatch
        integer                                     :: Nmodels, NmodelsFinal
        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: AuxChar
        integer                                     :: STAT_CALL, FirstLine, LastLine,  &
                                                       l, em, Naux, ClientNumber, line, iflag
        logical                                     :: BlockFound, MatchModel
        !----------------------------------------------------------------------

 
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                      &
                                    '<BeginModelPriority>', '<EndModelPriority>', BlockFound, &
                                    FirstLine = FirstLine, LastLine = LastLine,         &
                                    STAT = STAT_CALL)
    
        if (STAT_CALL /= SUCCESS_) stop 'ConstructEulerModelsList - ModuleLagrangianGlobal - ERR10'

BF:     if (BlockFound) then

            Naux = LastLine - FirstLine - 1

            if (Naux > Nmodels) stop 'ConstructEulerModelsList - ModuleLagrangianGlobal - ERR20'

            if (Naux < Nmodels) then
                write(*,*) 'The langrangian module will only use ',Naux,' from a total of ', Nmodels
                write(*,*) 'ConstructEulerModelsList - ModuleLagrangianGlobal - WRN20'
            endif

            em = 0
            do line = FirstLine + 1, LastLine - 1

                call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag, &
                             SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ConstructEulerModelsList - ModuleLagrangianGlobal - ERR30'

                em = em + 1
                MatchModel = .false. 
                do l = 1, Nmodels
                    if (trim(AuxChar) == trim(ModelNames(l))) then
                        IndexMatch(l) = em
                        MatchModel = .true.
                        exit
                    endif
                 enddo

                 if (.not.MatchModel) stop 'ConstructEulerModelsList - ModuleLagrangianGlobal - ERR40'

            enddo

        else
        
            do em = 1, Nmodels
                IndexMatch(em) = Nmodels - em + 1
            enddo 
            
            Naux = Nmodels
          
        endif BF

        NmodelsFinal = Naux


        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructEulerModelsList - ModuleLagrangianGlobal - ERR50'



        !----------------------------------------------------------------------

    end subroutine ConstructEulerModelsList
    !--------------------------------------------------------------------------

    
    subroutine ReadParticStatisticOptions

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: StatisticFound
        integer                                     :: STAT_CALL, ClientNumber

        !----------------------------------------------------------------------


        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                   statistic_begin, statistic_end,                  &
                                   StatisticFound, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadParticStatisticOptions - ModuleLagrangianGlobal - ERR10'

if1:    if (StatisticFound) then

            call ReadStatisticProperties(ClientNumber)

        else 
            
            if (Me%State%Statistics) stop 'ReadParticStatisticOptions - ModuleLagrangianGlobal - ERR20'

        endif if1


        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadParticStatisticOptions - ModuleLagrangianGlobal - ERR30'


    end subroutine ReadParticStatisticOptions

    !--------------------------------------------------------------------------

    subroutine ReadStatisticProperties(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                         :: ClientNumber

        !Local-----------------------------------------------------------------
        type(T_OptionsStat), dimension(:), pointer      :: AuxStatistic
        type(T_Property), pointer                       :: CurrentProperty
        type(T_PropertyID)                              :: PropertyID
        character(Len=PathLength)                       :: AuxChar
        logical                                         :: PropertyFound
        integer                                         :: STAT_CALL, nProp, i, j, flag

        !Begin-----------------------------------------------------------------



        allocate(AuxStatistic(Me%OriginDefault%nProperties))

        nProp = 0

DOPROP: do 
 
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                   &
                                       property_begin, property_end,                    &
                                       PropertyFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadStatisticProperties - ModuleLagrangianGlobal - ERR10'
            
i1:         if (PropertyFound) then

                call ConstructPropertyID (PropertyID, Me%ObjEnterData, FromBlockInBlock)

                nProp = nProp + 1

                AuxStatistic(nProp)%ID = PropertyID
 
                Me%State%Statistics = .true.

                !StatisticsFile
                call GetData(AuxChar,                                                   &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='STATISTICS_FILE',                           &
                             ClientModule ='ModuleLagrangianGlobal',                          &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadStatisticProperties - ModuleLagrangianGlobal - ERR40'

                AuxStatistic(nProp)%File = AuxChar

                if (flag /= 1) then
                    write(*,*)'No Statistics file definition for the Property : ',trim(CurrentProperty%Name)
                    write(*,*)'Use the KeyWord                                : STATISTICS_FILE'
                    write(*,*)'OR disable STATISTICS'
                    stop      'ReadStatisticProperties - ModuleLagrangianGlobal - ERR50'
                endif

                !Do frequency analysis considering the tracers concentration
                call GetData(AuxStatistic(nProp)%Lag,                                   &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='STATISTICS_LAG',                            &
                             ClientModule ='ModuleLagrangianGlobal',                          &
                             default      = .false.,                                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadStatisticProperties - ModuleLagrangianGlobal - ERR60'

            else i1

                !all properties were read
                exit 

            endif i1

        enddo DOPROP
        
        Me%Statistic%PropNumber = nProp

        allocate(Me%Statistic%OptionsStat(Me%Statistic%PropNumber))

        Me%Statistic%OptionsStat(1:Me%Statistic%PropNumber) = AuxStatistic(1:Me%Statistic%PropNumber)

        deallocate(AuxStatistic)
        nullify   (AuxStatistic)

        allocate(Me%Statistic%ON(Me%OriginDefault%NProperties))

        Me%Statistic%ON(:) = .false.

        j = 0
        CurrentProperty => Me%OriginDefault%FirstProperty
        do while(associated(CurrentProperty))
            j = j + 1
            do i=1, Me%Statistic%PropNumber
                if (CurrentProperty%ID == Me%Statistic%OptionsStat(i)%ID%IDNumber) then
                    Me%Statistic%ON(j)                    = .true.
                    Me%Statistic%OptionsStat(i)%PropOrder = j
                endif
            enddo
            CurrentProperty => CurrentProperty%Next
        enddo

        nullify(CurrentProperty)


    end subroutine ReadStatisticProperties
                
    !--------------------------------------------------------------------------

    subroutine ConstructLag2Euler

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                         :: WS_KLB, WS_KUB
        integer                                         :: nProp, em
        !Begin------------------------------------------------------------------
                
d1:     do em = 1, Me%EulerModelNumber 

            ILB    = Me%EulerModel(em)%Size%ILB
            IUB    = Me%EulerModel(em)%Size%IUB
            JLB    = Me%EulerModel(em)%Size%JLB
            JUB    = Me%EulerModel(em)%Size%JUB
            KLB    = Me%EulerModel(em)%Size%KLB
            KUB    = Me%EulerModel(em)%Size%KUB
            WS_ILB = Me%EulerModel(em)%WorkSize%ILB
            WS_IUB = Me%EulerModel(em)%WorkSize%IUB
            WS_JLB = Me%EulerModel(em)%WorkSize%JLB
            WS_JUB = Me%EulerModel(em)%WorkSize%JUB
            WS_KLB = Me%EulerModel(em)%WorkSize%KLB
            WS_KUB = Me%EulerModel(em)%WorkSize%KUB


            nProp           =  Me%OriginDefault%nProperties

            !Allocate GridVolume, GridMass    
            allocate (Me%EulerModel(em)%Lag2Euler%GridVolume      (ILB:IUB, JLB:JUB, KLB:KUB         , 1:Me%NGroups))
            allocate (Me%EulerModel(em)%Lag2Euler%GridTracerNumber(ILB:IUB, JLB:JUB, KLB:KUB         , 1:Me%NGroups))

            allocate (Me%EulerModel(em)%Lag2Euler%GridMass        (ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp, 1:Me%NGroups))
            allocate (Me%EulerModel(em)%Lag2Euler%GridConc        (ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp, 1:Me%NGroups))

            allocate (Me%EulerModel(em)%Lag2Euler%MeanConc        (                           1:nProp, 1:Me%NGroups))
            allocate (Me%EulerModel(em)%Lag2Euler%AmbientConc     (                           1:nProp, 1:Me%NGroups))
            allocate (Me%EulerModel(em)%Lag2Euler%MinConc         (                           1:nProp, 1:Me%NGroups))
            allocate (Me%EulerModel(em)%Lag2Euler%MassVolCel      (                           1:nProp, 1:Me%NGroups))

                                                        !i,j,k,p,ig 
            Me%EulerModel(em)%Lag2Euler%GridVolume      (:,:,:,  :) = 0.
            Me%EulerModel(em)%Lag2Euler%GridTracerNumber(:,:,:,  :) = 0.
            Me%EulerModel(em)%Lag2Euler%GridMass        (:,:,:,:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%GridConc        (:,:,:,:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%MeanConc              (:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%AmbientConc           (:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%MinConc               (:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%MassVolCel            (:,:) = 0.


            if (Me%State%Deposition) then
                allocate (Me%EulerModel(em)%Lag2Euler%GridBottomMass(ILB:IUB, JLB:JUB, 1:nProp, 1:Me%NGroups))
                allocate (Me%EulerModel(em)%Lag2Euler%GridBottomConc(ILB:IUB, JLB:JUB, 1:nProp, 1:Me%NGroups))
                                                          !i,j,k,p,ig 
                Me%EulerModel(em)%Lag2Euler%GridBottomMass(:,:,:,  :) = 0.
                Me%EulerModel(em)%Lag2Euler%GridBottomConc(:,:,:,  :) = 0.
            endif

            if (Me%OutPut%ConcMaxTracer) then
                allocate (Me%EulerModel(em)%Lag2Euler%GridMaxTracer(ILB:IUB, JLB:JUB, KLB:KUB, 1:nProp, 1:Me%NGroups))
                                                          !i,j,k,p,ig 
                Me%EulerModel(em)%Lag2Euler%GridMaxTracer (:,:,:,:,:) = 0.
            endif

        enddo d1


d2:     do em =1, Me%EulerModelNumber

            ILB = Me%EulerModel(em)%Size%ILB
            IUB = Me%EulerModel(em)%Size%IUB

            JLB = Me%EulerModel(em)%Size%JLB
            JUB = Me%EulerModel(em)%Size%JUB

            KLB = Me%EulerModel(em)%Size%KLB
            KUB = Me%EulerModel(em)%Size%KUB

            if (Me%State%Deposition)  then
        
                allocate   (Me%EulerModel(em)%Lag2Euler%TauErosionGrid (ILB:IUB, JLB:JUB, 1:Me%nGroups))
                allocate   (Me%EulerModel(em)%Lag2Euler%MassSedGrid    (ILB:IUB, JLB:JUB, 1:Me%nGroups))

                Me%EulerModel(em)%Lag2Euler%TauErosionGrid (:,:,:) = 0.
                Me%EulerModel(em)%Lag2Euler%MassSedGrid    (:,:,:) = 0.

            endif

            if (Me%State%Filtration)  then

                if (Me%nGroups > 1) stop 'ConstructLag2Euler - ModuleLagrangianGlobal - ERR10'

                allocate(Me%EulerModel(em)%RelativeMassFilter(ILB:IUB,JLB:JUB,KLB:KUB))

                allocate(Me%EulerModel(em)%MassFiltered      (ILB:IUB,JLB:JUB,KLB:KUB))

                Me%EulerModel(em)%RelativeMassFilter(:,:,:) = 0. 
                Me%EulerModel(em)%MassFiltered      (:,:,:) = 0. 
            
                Me%NextFiltration = Me%Now

            endif

        enddo d2


        !----------------------------------------------------------------------

    end subroutine ConstructLag2Euler

    !--------------------------------------------------------------------------

    subroutine ConstructParticleGrid(EulerModel)

        !Arguments-------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: i, j, STAT_CALL
        real,    dimension (:, :), pointer          :: XX_IE, YY_IE
        integer                                     :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG
        integer                                     :: GRID_COORD, CoordType, NLRD
        integer, dimension (:, :), pointer          :: DefineCellsMap

        !Begin-----------------------------------------------------------------



        !Shorten
        ILB    = EulerModel%Size%ILB
        IUB    = EulerModel%Size%IUB
        JLB    = EulerModel%Size%JLB
        JUB    = EulerModel%Size%JUB
        WS_ILB = EulerModel%WorkSize%ILB
        WS_IUB = EulerModel%WorkSize%IUB
        WS_JLB = EulerModel%WorkSize%JLB
        WS_JUB = EulerModel%WorkSize%JUB

        !Gets Coordinate Types List
        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(EulerModel%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR10'


        EulerModel%Grid%CoordType = CoordType

        if (CoordType == GEOG .or. CoordType == UTM .or. CoordType == SIMPLE_GEOG)      &
            EulerModel%Grid%HaveLatLongGrid = .true.

        if (CoordType == GEOG .or. CoordType == SIMPLE_GEOG)                            &
            EulerModel%Grid%GeoGrid = .true.


        call GetLatitudeLongitude(EulerModel%ObjHorizontalGrid, Latitude  = EulerModel%Grid%LatDefault,  &
                                                                Longitude = EulerModel%Grid%LongDefault, & 
                                                                STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR20'


        !Gets Horizontal Grid
        call GetHorizontalGrid(EulerModel%ObjHorizontalGrid, XX_IE = XX_IE, YY_IE = YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ER30'

        call GetDefineCellsMap(EulerModel%ObjHorizontalGrid, DefineCellsMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR40'

        allocate(EulerModel%Grid%ParticXX(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR50'

        allocate(EulerModel%Grid%ParticYY(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR60'


        EulerModel%Grid%ParticXX(:, :) = 0.
        EulerModel%Grid%ParticYY(:, :) = 0.

        do j = WS_JLB, WS_JUB
        do i = WS_ILB,   WS_IUB

            if (DefineCellsMap(i, j) == 1) then

                EulerModel%Grid%ParticXX(i  , j  ) = (XX_IE(i,  j  ) + XX_IE(i+1, j  )) / 2.
                EulerModel%Grid%ParticXX(i  , j+1) = (XX_IE(i,  j+1) + XX_IE(i+1, j+1)) / 2.
                EulerModel%Grid%ParticYY(i  , j  ) = (YY_IE(i,  j  ) + YY_IE(i  , j+1)) / 2.
                EulerModel%Grid%ParticYY(i+1, j  ) = (YY_IE(i+1,j  ) + YY_IE(i+1, j+1)) / 2.
            endif
        
        enddo
        enddo


        call UnGetHorizontalGrid(EulerModel%ObjHorizontalGrid, XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR70'

        call UnGetHorizontalGrid(EulerModel%ObjHorizontalGrid, YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR80'

        call UnGetHorizontalGrid(EulerModel%ObjHorizontalGrid, DefineCellsMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticleGrid - ModuleLagrangianGlobal - ERR90'


    end subroutine ConstructParticleGrid    

    !--------------------------------------------------------------------------

    subroutine ConstructOrigins

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: flag
        integer                                     :: STAT_CALL
        integer                                     :: ClientNumber, ClientNumberClone
        real                                        :: DT, DT_PARTIC, MinError
        logical                                     :: BlockFound
        type (T_Origin), pointer                    :: NewOrigin, OriginalOrigin
        logical                                     :: FoundCloneOrigin, ClonesExist
        integer                                     :: Nmax, no, NFirstExtraction
        logical                                     :: FirstExtraction

        !Begin-----------------------------------------------------------------

        nullify(Me%OriginDefault)


        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%ExternalVar%Now,                            &
                           EndTime     = Me%ExternalVar%EndTime,                        &
                           keyword     = 'OUTPUT_TIME',                                 &
                           SearchType  = FromFile,                                      &
                           OutPutsTime = Me%OutPut%OutTime,                             &
                           OutPutsOn   = Me%OutPut%Write_,                              &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR10'

        if (Me%OutPut%Write_) then
            Me%OutPut%NextOutPut = 1
        endif

        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%ExternalVar%Now,                            &
                           EndTime     = Me%ExternalVar%EndTime,                        &
                           keyword     = 'RESTART_FILE_OUTPUT_TIME',                    &
                           SearchType  = FromFile,                                      &
                           OutPutsTime = Me%OutPut%RestartOutTime,                      &
                           OutPutsOn   = Me%OutPut%WriteRestartFile,                    &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR11'

        if(Me%OutPut%WriteRestartFile)then

            Me%OutPut%NextRestartOutput = 1

        end if 

        !Checks wether to overwrite the Restart File OR not
        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='RESTART_FILE_OVERWRITE',                            &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     Default      = .true.,                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR12'

        !Output Concentration Type
        call GetData(Me%OutPut%OutPutConcType,                                          &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_CONC',                                       &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     Default      = Maximum,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR20'

        !Checks if the users wants to output the maximum tracer concentration in each cell
        call GetData(Me%OutPut%ConcMaxTracer,                                           &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_MAX_TRACER',                                 &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR30'


        call GetData(Me%OutPut%OriginEnvelope,                                          & 
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_ORIGIN_ENVELOPE',                            &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR40'

        !Write dummy particle in the start date if start emission > start date 
        call GetData(Me%OutPut%DummyParticleStartDate,                                  &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='WRITE_DUMMY_PART',                                  &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR45'



        call GetComputeTimeStep(Me%ExternalVar%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR70'

        !Get time step for particle from file
        DT_PARTIC = null_real
        call GetData(DT_PARTIC,                                                         &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='DT_PARTIC',                                         &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR80'
                          
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
                    stop      'ConstructOrigins - ModuleLagrangianGlobal - ERR90'
                endif
            elseif (DT > DT_PARTIC) then
                !if (amod(DT, DT_PARTIC) == 0.0) then
                    !OK
                !else
                !    write(*,*)'Particle Time step : ', DT_PARTIC
                !    write(*,*)'Model    Time step : ', DT
                !    write(*,*)'Particle Time step must multiply or submultiply of Model Time step'
                !    stop      'ConstructOrigins - ModuleLagrangianGlobal - ERR100'
                !endif
            
                !The run period must be a multiple of the model DT
                !The abs function is used, to avoid rounding erros
                !The old way was removed, to be able to run with Timesteps lower tehn 1 sec
                !Frank Dec - 2000
                MinError = min (abs(mod (DT, DT_PARTIC)),                               &
                                abs(DT_PARTIC - mod (DT, DT_PARTIC)))
                if (MinError >= 1.e-5) then
                    write(*,*)' Time step error DT_PARTIC - Run period must be a multiple of DT'
                    stop      'ConstructOrigins - ModuleLagrangianGlobal - ERR110.' 
                endif
            endif
        end if

        !Sets time step
        Me%DT_Partic   = DT_PARTIC
        Me%NextCompute = Me%ExternalVar%Now + Me%DT_Partic
        Me%Now         = Me%ExternalVar%Now
        
        call BoxTypeVariablesDefiniton


        !OVERLAY_VELOCITY
        call GetData(Me%Overlay,                                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OVERLAY_VELOCITY',                                  &
                     ClientModule ='ModuleLagrangianGlobal',                            &  
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR240'


        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR250'

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
                     ClientModule ='ModuleLagrangianGlobal',                            &  
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR255'


!#ifdef _CGI_

        if (Me%RunOnline) then

            call ReadQueryString(Nmax)

        endif

!#endif
        !This option allows the definition of discharges in land points without crashing the model
        !When this happens a message is send to the user and the discharge is not considered 
        !This is something done in the framework of COWAMA project
        call GetData(Me%IgnoreON,                                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='IGNORE_ON',                                         &
                     ClientModule ='ModuleLagrangianGlobal',                            &  
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR257'

        !Number of times it read the lagrangian data looking for origins
        !Only when the _CGI_ option is on is able to read several times 
        !the origin blocks
SB:     do no = 1, Nmax

DW:     do  

            if (FoundCloneOrigin) then

                BlockFound      = .true. 

            else 

                call ExtractBlockFromBuffer(Me%ObjEnterData,                            &
                                            ClientNumber,                               &
                                            block_begin, block_end,                     &
                                            BlockFound,                                 &
                                            STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR260'

            endif

BF:         if (BlockFound) then 

                if (FirstExtraction) then 
                    NFirstExtraction = NFirstExtraction + 1
                endif

                !Done for the COWAMA project. Allows not consider discharges define in land points without stopping the model  
                if (.not.(Me%IgnoreON .and. CheckOriginInLandCell())) then
                    
                    !Allocates a new origin
                    call AllocateNewOrigin(NewOrigin)

                    if (.not. FoundCloneOrigin) then
                        nullify(OriginalOrigin)
                        OriginalOrigin => NewOrigin
                    endif


                    call ConstructOneOrigin(NewOrigin, ClientNumber)              

                    if (NewOrigin%Default) then

                        Me%OriginDefault => NewOrigin

                    else

                        !Insert New Origin to the list of origins    
                        call InsertOriginToList (Me%FirstOrigin, NewOrigin, Me%nOrigins)

                    endif

                endif

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

        if (.not. Me%State%OriginDefault) Me%OriginDefault => Me%FirstOrigin


        !Finished reading block -> unlocks block reading
        call Block_Unlock(Me%ObjEnterData,                                              &
                          ClientNumber,                                                 &
                          STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR1600'

        if (Me%ObjEnterDataClone   /= 0)  then
            !Finished reading block -> unlocks block reading
            call Block_Unlock(Me%ObjEnterDataClone,                                     &
                              ClientNumberClone,                                        &
                              STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR1610'

            call KillEnterData(Me%ObjEnterDataClone, STAT = STAT_CALL)  

            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR1620'
        endif

        if (Me%ObjEnterDataOriginal /= 0) then
            call KillEnterData(Me%ObjEnterDataOriginal, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR1630'
        endif

!#ifdef _CGI_

        if (Me%RunOnline) then

            call ChangeOriginOnline(Nmax, NFirstExtraction)
        endif

!#endif

        nullify(OriginalOrigin, NewOrigin)

    end subroutine ConstructOrigins

    !--------------------------------------------------------------------------


    subroutine ConstructOneOrigin(NewOrigin, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: NewOrigin
        integer                                     :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                     :: flag
        integer                                     :: STAT_CALL
        character(3)                                :: Aux
        character(StringLength)                     :: EmissionSpatial
        character(StringLength)                     :: EmissionTemporal
        character(PathLength)                       :: String, String2, RootPath
        real, dimension(:), allocatable             :: Aux2
        real, dimension(1:2)                        :: Position
        real                                        :: Depth , TotalVolume
        logical                                     :: HaveOrigin, PropertyFound, NoDomain
        type (T_Property), pointer                  :: NewProperty
        integer                                     :: i, j, k
        logical                                     :: ret
        integer                                     :: PropertyID, DensityMethod
        logical                                     :: WP_HaveProperty, SedimentDefined = .false.
        logical                                     :: SalOK = .false., TempOK = .false., PressureCorrection
        integer                                     :: em

        !Begin-----------------------------------------------------------------

        call GetData(NewOrigin%Default,                                         &
                     Me%ObjEnterData,                                           &
                     flag,                                                      &
                     SearchType   = FromBlock,                                  &
                     keyword      ='DEFAULT',                                   &
                     ClientModule ='ModuleLagrangianGlobal',                    &
                     Default      = OFF,                                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR250'

        if (NewOrigin%Default) then
            if  (Me%State%OriginDefault) then
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR260'
            else
                Me%State%OriginDefault = .true. 
            endif
        endif

        !Gets its name    
        call GetData(NewOrigin%Name,                                            &
                     Me%ObjEnterData,                                           &
                     flag,                                                      &
                     SearchType   = FromBlock,                                  &
                     keyword      ='ORIGIN_NAME',                               &
                     ClientModule ='ModuleLagrangianGlobal',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR270'

        if (flag == 0) then
            Aux = ' '
            write(Aux,'(i3)') Me%nOrigins + 1
            NewOrigin%Name = trim(adjustl('Origin_'//trim(adjustl(Aux))))
        end if

        !Old Origin    
        call GetData(NewOrigin%Old,                                              &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='OLD',                                        &
                     ClientModule ='ModuleLagrangianGlobal',                     &
                     Default      = .false.,                                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR280'

        if (NewOrigin%Old) then
            Me%State%ContCalc = .true.
        endif

        !Group ID
        call GetData(NewOrigin%GroupID,                                          &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='GROUP_ID',                                   &
                     ClientModule ='ModuleLagrangianGlobal',                     &
                     Default      = 1,                                           &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR290'


        !Gets the spatial emission type
        call GetData(EmissionSpatial,                                            &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='EMISSION_SPATIAL',                           &
                     ClientModule ='ModuleLagrangianGlobal',                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR300'

        !Gets the temporal emission type
        call GetData(EmissionTemporal,                                           &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='EMISSION_TEMPORAL',                          &
                     ClientModule ='ModuleLagrangianGlobal',                     &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR310'



        !Gets the temporal emission type
        call GetData(NewOrigin%EmissionON,                                              &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromBlock,                                          &
                     keyword      ='EMISSION_ON',                                       &
                     default      = ON,                                                 &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR315'


        !Converts Strings to int
        call ConstructEmissionType(NewOrigin, EmissionSpatial, EmissionTemporal)

IT:     if (NewOrigin%EmissionTemporal == Instantaneous_) then

            !Gets the interval between emissions
            call GetData(NewOrigin%InstantEmission,                              &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='INSTANT_PARTIC_EMIT',                    &
                         ClientModule ='ModuleLagrangianGlobal',                 &
                         Default      = Me%ExternalVar%BeginTime,                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR320'
            

        endif IT
        !
ET:     if (NewOrigin%EmissionTemporal == Continuous_) then

            !Gets the interval between emissions
            call GetData(NewOrigin%DT_Emit,                                      &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='DT_EMIT',                                &
                         ClientModule ='ModuleLagrangianGlobal',                 &
                         Default      = Me%DT_Partic,                            &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR325'

            !Gets the interval between emissions
            call GetData(NewOrigin%StartEmission,                                &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='START_PARTIC_EMIT',                      &
                         ClientModule ='ModuleLagrangianGlobal',                 &
                         Default      = Me%ExternalVar%BeginTime,                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR330'

            !Gets the interval between emissions
            call GetData(NewOrigin%StopEmission,                                 &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='STOP_PARTIC_EMIT',                       &
                         ClientModule ='ModuleLagrangianGlobal',                 &
                         Default      = Me%ExternalVar%EndTime,                  &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR340'


            !Gets flow associated to a continuous emission
iP:         if (NewOrigin%EmissionSpatial == Point_) then


                !Flow variable in time
                call GetData(NewOrigin%FlowVariable,                            &
                             Me%ObjEnterData,                                   &
                             flag,                                              &
                             SearchType   = FromBlock,                          &
                             keyword      ='FLOW_VARIABLE',                     &
                             ClientModule ='ModuleLagrangianGlobal',            &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR360'

iFV:            if (NewOrigin%FlowVariable) then

                    !Discharge file name
                    call GetData(NewOrigin%DischargeFile,                               &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlock,                              &
                                 keyword      ='DISCHARGE_FILE',                        &
                                 ClientModule ='ModuleLagrangianGlobal',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR370'

                    !Discharge file name
                    call GetData(NewOrigin%FlowColumn,                          &
                                 Me%ObjEnterData,                               &
                                 flag,                                          &
                                 SearchType   = FromBlock,                      &
                                 keyword      ='FLOW_COLUMN',                   &
                                 ClientModule ='ModuleLagrangianGlobal',        &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR390'

                    call GetData(NewOrigin%EstimateMinVol,                      &
                                 Me%ObjEnterData,                               &
                                 flag,                                          &
                                 SearchType   = FromBlock,                      &
                                 keyword      ='ESTIMATE_MIN_VOL',              &
                                 Default      = .false.,                        &
                                 ClientModule ='ModuleLagrangianGlobal',        &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR392'

iMV:                if (NewOrigin%EstimateMinVol) then

                        call GetData(NewOrigin%MaxPart,                                 &
                                     Me%ObjEnterData,                                   &
                                     flag,                                              &
                                     SearchType   = FromBlock,                          &
                                     keyword      ='MAX_PART',                          &
                                     Default      = 10000,                              &
                                     ClientModule ='ModuleLagrangianGlobal',            &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR393'

                    endif iMV

iDF:                if (.not. NewOrigin%Default) then
                        call StartTimeSerieInput(NewOrigin%TimeSerieInputFlow,          &
                                                 NewOrigin%DischargeFile,               &
                                                 Me%ExternalVar%ObjTime,                &
                                                 STAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR380'

                    endif iDF

                else iFV

                    !If not variable, get flow of the origin
                    call GetData(NewOrigin%Flow,                                &
                                 Me%ObjEnterData,                               &
                                 flag,                                          &
                                 SearchType   = FromBlock,                      &
                                 keyword      ='FLOW',                          &
                                 ClientModule ='ModuleLagrangianGlobal',        &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR394'
                    if (flag /= 1) then
                        write(*,*)'Keyword FLOW not defined at origin :',trim(adjustl(NewOrigin%Name))
                        stop      'ConstructOneOrigin - ModuleLagrangianGlobal - ERR27'
                
                    endif 


                endif iFV



                call GetData(NewOrigin%MovingOrigin,                            &
                             Me%ObjEnterData,                                   &
                             flag,                                              &
                             SearchType   = FromBlock,                          &
                             default      = .false.,                            &
                             keyword      ='MOVING_ORIGIN',                     &
                             ClientModule ='ModuleLagrangianGlobal',            &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR400'

                if (NewOrigin%MovingOrigin) then

                    call GetData(NewOrigin%MovingOriginFile,                     &
                                 Me%ObjEnterData,                     &
                                 flag,                                           &
                                 SearchType   = FromBlock,                       &
                                 keyword      ='MOVING_ORIGIN_FILE',             &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_ .or. flag /= 1) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR410'

                    call GetData(NewOrigin%MovingOriginUnits,                    &
                                 Me%ObjEnterData,                     &
                                 flag,                                           &
                                 default      = 'Cells',                         &
                                 SearchType   = FromBlock,                       &
                                 keyword      ='MOVING_ORIGIN_UNITS',            &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR420'


                    call GetData(NewOrigin%MovingOriginColumnX,                  &
                                 Me%ObjEnterData,                     &
                                 flag,                                           &
                                 SearchType   = FromBlock,                       &
                                 keyword      ='MOVING_ORIGIN_COLUMN_X',         &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_ .or. flag /= 1) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR430'

                    call GetData(NewOrigin%MovingOriginColumnY,                  &
                                 Me%ObjEnterData,                     &
                                 flag,                                           &
                                 SearchType   = FromBlock,                       &
                                 keyword      ='MOVING_ORIGIN_COLUMN_Y',         &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_ .or. flag /= 1) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR440'


                    call StartTimeSerieInput(NewOrigin%ObjTimeSerie,            &
                                             NewOrigin%MovingOriginFile,        &
                                             Me%ExternalVar%ObjTime,            &
                                             STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR450'


                endif

            endif iP

        endif ET

        if (NewOrigin%EmissionSpatial  == Accident_      .or.                           &
           (NewOrigin%EmissionSpatial  == Point_         .and.                          &
           (NewOrigin%EmissionTemporal == Instantaneous_ .or. NewOrigin%FlowVariable))) then
            call GetData(NewOrigin%PointVolume,                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromBlock,                                      &
                         keyword      ='POINT_VOLUME',                                  &
                         ClientModule ='ModuleLagrangianGlobal',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR460'
            if (flag /= 1) then
                write(*,*)'Keyword POINT_VOLUME not defined at origin :',trim(adjustl(NewOrigin%Name))
                stop      'ConstructOneOrigin - ModuleLagrangianGlobal - ERR29'
            endif
        endif

i23:    if (NewOrigin%TimeSerieInputFlow /= 0) then

            if (.not. NewOrigin%Default .and. NewOrigin%EstimateMinVol) then

                TotalVolume = GetTimeSerieIntegral(NewOrigin%TimeSerieInputFlow,        &
                                                   Me%ExternalVar%BeginTime, Me%ExternalVar%EndTime, &
                                                   NewOrigin%FlowColumn, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR384.'

                if (TotalVolume > 0.) then
                    NewOrigin%PointVolume = max (TotalVolume/real(NewOrigin%MaxPart), Me%OriginDefault%PointVolume)
                endif

            endif 
        endif i23

        !Time to double volume
        call GetData(NewOrigin%Movement%TVOL200,                                 &
                     Me%ObjEnterData,                                 &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='TVOL200',                                    &
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR470'

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
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR480'


        !Volume Factor when Particle dies
        call GetData(NewOrigin%Movement%VOLFAC,                                  &
                     Me%ObjEnterData,                                 &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      = 'VOLFAC',                                    &
                     Default      = 10.,                                         &
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR490'

        !How to calculate volume increase
        call GetData(String2,                                                    &
                     Me%ObjEnterData,                                 &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      = 'VOLUME_INCREASE',                           &
                     Default      = Char_Double,                                 &
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR500'

        select case (trim(adjustl(String2)))

        case (Char_Velocity)

            NewOrigin%Movement%TVolType = Velocity_

        case (Char_Double)

            NewOrigin%Movement%TVolType = Double_

        case default

            write(*,*)'Invalid option for keyword VOLUME_INCREASE'
            write(*,*)
            stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR510'

        end select


        !Floating Particle
        call GetData(NewOrigin%Movement%Float,                                   &
                     Me%ObjEnterData,                                 &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='FLOAT',                                      &
                     Default      = OFF,                                         &
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR520'

    
        !Calculates Plume
        call GetData(NewOrigin%State%FarFieldBuoyancy,                                  &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromBlock,                                          &
                     keyword      ='COMPUTE_BUOYANCY',                                  &
                     Default      = OFF,                                                &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR530'

        
        if (NewOrigin%State%FarFieldBuoyancy) then

            Me%State%FarFieldBuoyancy = .true.

            NewOrigin%Movement%InitialVelocityW = 0.

        endif

   
        !Type of Accident
AC:             if (NewOrigin%EmissionSpatial == Accident_) then

            call GetData(NewOrigin%AccidentMethod,                                      &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromBlock,                                      &
                         keyword      ='ACCIDENT_METHOD',                               &
                         ClientModule ='ModuleLagrangianGlobal',                        &
                         Default      = Fay_,                                           &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR670'

            call GetData(NewOrigin%AccidentTime,                                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromBlock,                                      &
                         keyword      ='ACCIDENT_TIME',                                 &
                         ClientModule ='ModuleLagrangianGlobal',                        &
                         Default      = Me%ExternalVar%BeginTime,                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR680'

            if (flag>0) NewOrigin%StartEmission = NewOrigin%AccidentTime
        endif AC

        !Particle Thickness                
        if (NewOrigin%Movement%Float) then
            Me%State%Wind = ON

            if (NewOrigin%EmissionSpatial  == Box_ .or.                                 &
                (NewOrigin%EmissionSpatial == Accident_ .and.                           &
                 NewOrigin%AccidentMethod  == Thickness_)) then                         
                call GetData(NewOrigin%Movement%ThicknessMeters,                        &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlock,                                  &
                             keyword      ='THICKNESS_METERS',                          &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR690'
                if (flag /= 1) then
                    write(*,*)'Keyword THICKNESS_METERS not defined at origin :',trim(adjustl(NewOrigin%Name))
                    stop      'ConstructOneOrigin - ModuleLagrangianGlobal - ERR38'
                endif
            endif

        endif


        !Movement type
        call GetData(String,                                                            &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromBlock,                                          &
                     keyword      ='MOVEMENT',                                          &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR700'

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
                             ClientModule ='ModuleLagrangianGlobal',                   &  
                             Default      = 0.2,                                 &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR710'


                call GetData(NewOrigin%Movement%VarVelH,                         &
                             Me%ObjEnterData,                         &
                             flag,                                               &
                             SearchType   = FromBlock,                           &
                             keyword      ='VARVELH',                            &
                             ClientModule ='ModuleLagrangianGlobal',                   &  
                             Default      = 0.0,                                 &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR720'


                call GetData(String2,                                            &
                             Me%ObjEnterData,                         &
                             flag,                                               &
                             SearchType   = FromBlock,                           &
                             keyword      ='TURB_V',                             &
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR730'

TURB_V:                 if (flag == 1) then   

                    select case(trim(adjustl(String2)))
                    
                    case(Char_Constant)

                        NewOrigin%Movement%StandardDeviationType = VerticalTurbConstant

                        call GetData(NewOrigin%Movement%VarVelVX,                &
                                     Me%ObjEnterData,                 &
                                     flag,                                       &
                                     SearchType   = FromBlock,                   &
                                     keyword      ='VARVELVX',                   &
                                     ClientModule ='ModuleLagrangianGlobal',           &  
                                     Default      = 0.0,                         &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR740'


                        call GetData(NewOrigin%Movement%VarVelV,                 &
                                     Me%ObjEnterData,                 &
                                     flag,                                       &
                                     SearchType   = FromBlock,                   &
                                     keyword      ='VARVELV',                    &
                                     ClientModule ='ModuleLagrangianGlobal',           &  
                                     Default      = 0.0,                         &
                                     STAT         = STAT_CALL)             
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR750'

                    case(Char_Profile )

                        NewOrigin%Movement%StandardDeviationType = VerticalTurb
                        
                        !Will need shear velocity    
                        Me%State%ShearVel = ON
                        
                        !Sets unused variables to dummy
                        NewOrigin%Movement%VarVelVX  = null_real
                        NewOrigin%Movement%VarVelV   = null_real

                    case default

                        write(*,*)'Invalid option for TURB_V'
                        stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR760'

                    end select

                else TURB_V
                
                    write(*,*)'Keyword TURB_V not found'
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR780'

                endif TURB_V

            case(Char_NotRandom    )

                NewOrigin%Movement%MovType = NotRandom_

            case default

                write(*,*)'Invalid horizontal movement keyword'
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR790'
         
            end select

        else MO

            write(*,*)'Keyword MOVEMENT not defined at origin', trim(adjustl(NewOrigin%Name))
            stop      'ConstructOneOrigin - ModuleLagrangianGlobal - ERR800'
        end if MO


        !Advection
        call GetData(NewOrigin%Movement%Advection,                               &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='ADVECTION',                                  &
                     ClientModule ='ModuleLagrangianGlobal',                           &  
                     Default      = ON,                                          &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR810'

        !Kill Land Particles
        call GetData(NewOrigin%Movement%KillLandParticles,                       &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='KILL_LAND_PARTICLES',                        &
                     ClientModule ='ModuleLagrangianGlobal',                           &  
                     Default      = OFF,                                         &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR820'



        !WINDCOEF
        call GetData(NewOrigin%Movement%WindTransferCoef,                        &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='WINDCOEF',                                   &
                     ClientModule ='ModuleLagrangianGlobal',                           &  
                     Default      = 0.03,                                        &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR830'


        allocate(Aux2(2))

        Aux2(1:2) = 0

        !WINDXY
        call GetData(Aux2,                                                       &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='WINDXY',                                     &
                     ClientModule ='ModuleLagrangianGlobal',                           &  
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR840'

        
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
                     ClientModule ='ModuleLagrangianGlobal',                           &  
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR850'

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
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR860'

                NewOrigin%State%Sedimentation = ON

            case(Char_Imposed)
            
                NewOrigin%Movement%SedimentationType = Imposed_

                call GetData(NewOrigin%Movement%SedVel,                          &
                             Me%ObjEnterData,                         &
                             flag,                                               &
                             SearchType   = FromBlock,                           &
                             keyword      ='SED_VELOCITY',                       &
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR870'
                if (flag == 0) then
                    write(*,*)'Sedimentation velocity not defined, keyword SED_VELOCITY'
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR880'
                endif

                NewOrigin%State%Sedimentation     = ON

            case default
            
                write(*,*)'Invalid Sedimentaion type, keyword SEDIMENTATION'
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR890'

            end select

        end if SE

       !Particles with deposition 
        call GetData(NewOrigin%State%Deposition,                                 &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='DEPOSITION',                                 &
                     default      = .false.,                                     & 
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            call SetError(FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR900')

        if (NewOrigin%State%Deposition .and. .not.NewOrigin%State%Sedimentation) &
            call SetError(FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR910')


DE:             if (NewOrigin%State%Deposition) then

            Me%State%Deposition = ON

            call GetData(NewOrigin%Deposition%TauErosion,                        &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='TAU_ERO',                                &
                         default      = 0.2,                                     &                                  
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                           &
                call SetError(FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR920')

            call GetData(NewOrigin%Deposition%TauDeposition,                     &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='TAU_DEP',                                &
                         default      = 0.1,                                     &                                  
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                           &
                call SetError(FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR930')

           call GetData(NewOrigin%Deposition%BottomDistance,                     &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='BOTTOM_DISTANCE',                        &
                         default      = 0.1,                                     &
                         ClientModule ='ModuleLagrangianGlobal',                 &
                         STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                           &
                call SetError(FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR940')

           call GetData(NewOrigin%Deposition%Tdecay,                             &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='TIME_DECAY',                             &
                         !defaul 2 days
                         default      = 172800.,                                 &
                         ClientModule ='ModuleLagrangianGlobal',                 &
                         STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                           &
                call SetError(FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR950')

           call GetData(NewOrigin%Deposition%BottomEmission,                     &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='BOTTOM_EMISSION',                        &
                         !by default the particle are emitted in the water column
                         default      = .false.,                                 &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                           &
                call SetError(FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR960')

           call GetData(NewOrigin%Deposition%ErosionRate,                        &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='EROSION_RATE',                           &
                         !by default the erosion rate in 5e-2 g/m2/s             
                         !This value make sense if the concentration is in mg/l = g/m3
                         default      = 5.e-2,                                   &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                           &
                call SetError(FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR970')

        endif DE


        !Min Sedimentation velocity
        call GetData(NewOrigin%Movement%MinSedVel,                               &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   =  FromBlock,                                  &
                     keyword      = 'MIN_SED_VELOCITY',                          &
                     default      =  0.0,                                        &
                     ClientModule = 'ModuleLagrangianGlobal',                          &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR980'


        !Reads parameter specific to cada Spatial emission type
PA:             if (NewOrigin%EmissionSpatial == Point_ .or.                             &
                    NewOrigin%EmissionSpatial == Accident_ ) then

            !Gets the number of particles to emit
            call GetData(NewOrigin%NbrParticlesIteration,                        &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='NBR_PARTIC',                             &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         Default      = 1,                                       &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR990'


            !Horizontal position in meters
            Position(:) = null_real
            HaveOrigin  = .false.
              
NDF:        if (.not. NewOrigin%Default) then
        
                !Horizontal position in coordinates X, Y
                call GetData(Position,                                  &
                             Me%ObjEnterData,                           &
                             flag,                                      &
                             SearchType   = FromBlock,                  &
                             keyword      ='POSITION_COORDINATES',      &
                             ClientModule ='ModuleLagrangianGlobal',          &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1020'

                if (flag == 2) then
                    HaveOrigin                      = .true.
                endif


                em = Locate_ModelDomain(Position(1), Position(2), NoDomain) 


                if (NoDomain) then
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1023'
                endif

                NewOrigin%Position%CoordX = Position(1)
                NewOrigin%Position%CoordY = Position(2)

                NewOrigin%Position%ModelID =  em

                call Convert_XY_CellIJ(Me%EulerModel(em),NewOrigin%Position, Referential = GridCoord_)

                    
                !Checks if a valid horizontal origin was found
                if (.not. HaveOrigin) then
                    write(*,*)'No Valid Horizontal Location defined for ',trim(adjustl(NewOrigin%Name))
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR66'

                endif

                !Verifies the horizontal location of the origin
                i = NewOrigin%Position%I
                j = NewOrigin%Position%J
                k = Me%EulerModel(em)%WorkSize%KUB

                if (Me%EulerModel(em)%Waterpoints3D(i, j, k) /= WaterPoint) then
                    write(*,*) 'Discharge in a land cell I=',i,' J=',j,'Model name=',trim(Me%EulerModel(em)%Name)

                    write(*,*)'Invalid Location defined for ',trim(adjustl(NewOrigin%Name))
                    write(*,*)'Point [i]:', NewOrigin%Position%I
                    write(*,*)'Point [j]:', NewOrigin%Position%J
                    
                    write(*,*)'Is not a WaterPoint'
                    
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1030'
                endif

                  !Vertical position
                Depth      = null_real
                HaveOrigin = .false.
                call GetData(Depth,                                                 &
                             Me%ObjEnterData,                                       &
                             flag,                                                  &
                             SearchType   = FromBlock,                              &
                             keyword      ='DEPTH_METERS',                          &
                             ClientModule ='ModuleLagrangianGlobal',                &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1040'

                if (flag == 1) then

                    NewOrigin%Position%Z = Depth

                    call Convert_Z_CellK (NewOrigin, Me%EulerModel(em), NewOrigin%Position)
                    call Convert_CellK_K (NewOrigin%Position)

                    NewOrigin%Position%DepthDefinition = Meters

                    HaveOrigin                      = .true.

                else

                    call GetData(Depth,                                              &
                                 Me%ObjEnterData,                                    &
                                 flag,                                               &
                                 SearchType   = FromBlock,                           &
                                 keyword      ='DEPTH_CELLS',                        &
                                 ClientModule ='ModuleLagrangianGlobal',             &
                                 STAT         = STAT_CALL)             
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1050'

                    if (flag == 1) then

                        NewOrigin%Position%CellK = Depth

                        call Convert_CellK_Z (Me%EulerModel(em),NewOrigin%Position)
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
                             ClientModule ='ModuleLagrangianGlobal',                       &
                             STAT         = STAT_CALL)             
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1060'

                if (NewOrigin%Position%MaintainRelative) then
            
                    NewOrigin%Position%Depth = Depth

                endif

                !Initial Position of floating Particle close to the surface
                if (NewOrigin%Movement%Float) then

                    i = NewOrigin%Position%I
                    j = NewOrigin%Position%J
                    k = Me%EulerModel(em)%WorkSize%KUB
                    NewOrigin%Position%Z = Me%EulerModel(em)%SZZ(i, j, k)
            
                    call Convert_Z_CellK  (NewOrigin, Me%EulerModel(em), NewOrigin%Position)
                    call Convert_CellK_K  (NewOrigin%Position)

                    HaveOrigin = .true.

                endif

                !Checks if a valid vertical origin was found
                if (.not. HaveOrigin) then

                    write(*,*)'No Valid Vertical Location defined for ',trim(adjustl(NewOrigin%Name))
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1080'

                endif

                !Verifies the location of the origin
                i = NewOrigin%Position%I
                j = NewOrigin%Position%J
                k = NewOrigin%Position%K
                if (Me%EulerModel(em)%Waterpoints3D(i, j, k) /= WaterPoint) then
                    write(*,*)'Invalid Location defined for ',trim(adjustl(NewOrigin%Name))
                    write(*,*)'Point [i]:', NewOrigin%Position%I
                    write(*,*)'Point [j]:', NewOrigin%Position%J
                    write(*,*)'Point [k]:', NewOrigin%Position%K
                    write(*,*)'Is not a WaterPoint'
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1090'
                endif

            else  NDF

                em = 1

            endif NDF

        endif PA

BX:     if (NewOrigin%EmissionSpatial == Box_) then

            call GetData(NewOrigin%INCRP,                                        &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='INCRP',                                  &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         Default      = 1,                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1100'


            !Volume of the particles inside the
            !If not given OneCell _> One particle
            call GetData(NewOrigin%ParticleBoxVolume,                            &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='BOXVOLINIC',                             &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         Default      = null_real,                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1110'

            !Box Number
            call GetData(NewOrigin%BoxNumber,                                    &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='BOX_NUMBER',                             &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1120'
            if (flag /= 1) then
                write(*,*)'Keyword BOX_NUMBER not defined at origin :', trim(adjustl(NewOrigin%Name))
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1130'
            endif

            if (NewOrigin%EmissionTemporal == Instantaneous_) then

                !kill all particles inside the box previous to the box emission 
                call GetData(NewOrigin%KillPartInsideBox,                               &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlock,                                  &
                             keyword      ='KILL_PART_INSIDE_BOX',                      &
                             Default      = OFF,                                        &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1125'

                if (NewOrigin%KillPartInsideBox) Me%State%KillPartInsideBox = .true. 

            endif

        endif BX

        !Calculates Plume
        call GetData(NewOrigin%State%ComputePlume,                               &
                     Me%ObjEnterData,                                            &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='COMPUTE_PLUME',                              &
                     Default      = OFF,                                         &
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR540'

        !Plume Data
        if (NewOrigin%State%ComputePlume) then


            Me%State%ComputePlume = .true.
            
            if (.not. NewOrigin%State%FarFieldBuoyancy) then
                write(*,*) 'The MOHID JET is ON but'
                write(*,*) 'no buoyancy in the far field is being compute'
                write(*,*) 'ConstructOneOrigin - ModuleLagrangianGlobal - WRN10'
            endif

            !Momentum balance in horizontal direcction
            !Increase of volume tracer due to shear effect
            call GetData(NewOrigin%State%PlumeShear,                             &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='PLUME_SHEAR',                            &
                         Default      = .true.,                                  &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR550'


                                
            !Plume Friction
            call GetData(NewOrigin%Movement%CoefInitialMixing,                   &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='COEF_INITIAL_MIXING',                    &
                         Default      = 1.0,                                     &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR560'


            !Plume Jet Data file
            call GetData(NewOrigin%Movement%JetDataFile,                         &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='JET_DATA_FILE',                          &
                         Default      = "********.***",                          &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR570'
            if (flag == 0)             stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR580'

            !Time interval for the actualization of Plume Jet properties
            call GetData(NewOrigin%Movement%JetDT,                               &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='JET_DT',                                 &
                         Default      = 600.0,                                   &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR590'


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
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR600'
            
            call GetData(NewOrigin%Movement%CorrecPress,                        &
                         Me%ObjEnterData,                                       &
                         flag,                                                  &
                         SearchType   = FromBlock,                              &
                         keyword      ='PRESSURE_CORRECTION',                   &
                         Default      = .true.,                                 &
                         ClientModule ='ModuleLagrangianGlobal',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR610'

            call GetDensityOptions(Me%EulerModel(em)%ObjWaterProperties, DensityMethod,        &
                                   PressureCorrection, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR620'
            
            if (NewOrigin%Movement%DensityMethod/=DensityMethod) then
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR630'
            endif
            
            if (.not. (NewOrigin%Movement%CorrecPress .EQV. PressureCorrection)) then
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR640'
            endif                    
       
        endif

        call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
    
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR650'

       
        !Searches for Properties
DOPROP:         do 
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,            &
                                       property_begin, property_end,             &
                                       PropertyFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1140'
            if (PropertyFound) then

                !Allocates NewProperty
                call AllocateNewProperty(NewProperty)

                !Name of the Property
                call GetData(NewProperty%Name,                                   &
                             Me%ObjEnterData,                                    &
                             flag,                                               &
                             SearchType   = FromBlockInBlock,                    &
                             keyword      ='NAME',                               &
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1150'
                if (flag == 0) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1160'


                !Checks if the property is recognized by the model
                if (.not. CheckPropertyName(trim(NewProperty%Name), NewProperty%ID)) then
                    write (*,*)'Unknown Property : ', trim(NewProperty%Name)
                    stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1170'
                end if

                !Units of the Property
                call GetData(NewProperty%Units,                                         &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='UNITS',                                     &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1180'

                !Concentration of the Property
                call GetData(NewProperty%Concentration,                                 &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='CONCENTRATION',                             &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1190'
                if (flag == 0) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1200'


                !Concentration variable in time
                call GetData(NewProperty%ConcVariable,                                  &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='CONC_VARIABLE',                             &
                             default      = .false.,                                    &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1210'

                !Concentration of the Property equal to the ambient concentration
                call GetData(NewProperty%EqualToAmbient,                                &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='EQUAL_TO_AMBIENT',                          &
                             default      = .false.,                                    &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1202'


                if (NewProperty%EqualToAmbient) NewProperty%ConcVariable = .false. 

                if (NewProperty%ConcVariable) then

                    !Discharge file name
                    call GetData(NewProperty%DischargeFile,                             &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      ='DISCHARGE_FILE',                        &
                                 ClientModule ='ModuleLagrangianGlobal',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1212'

                    if (.not. NewOrigin%Default) then
                        call StartTimeSerieInput(NewProperty%TimeSerieInput,            &
                                                 NewProperty%DischargeFile,             &
                                                 Me%ExternalVar%ObjTime,                &
                                                 STAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1214'
                    endif

                    !Concentration column
                    call GetData(NewProperty%ConcColumn,                                &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      ='CONC_COLUMN',                           &
                                 ClientModule ='ModuleLagrangianGlobal',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1220'
                    if (flag == 0) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1230'

                endif


                !No WQM
                 call GetData(NewProperty%NoWQM,                                        &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      ='NOWQM',                                 &
                                 ClientModule ='ModuleLagrangianGlobal',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1240'


                !MinimunConcentration 
                call GetData(NewProperty%Min_concentration,                             &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='MIN_CONCENTRATION',                         &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             Default      = 0.0,                                        &    
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1250'


                !Time series output
                call GetData(NewProperty%WritesTimeSerie,                               &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             Keyword      = 'TIME_SERIE',                               &
                             ClientModule = 'ModuleLagrangianGlobal',                   &
                             Default      = .false.,                                    &
                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR1255'

                if (NewProperty%WritesTimeSerie) Me%WritesTimeSerie = .true. 


                !Time series output
                call GetData(NewProperty%WritesPropHDF,                                 &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             Keyword      = 'OUTPUT_HDF',                               &
                             ClientModule = 'ModuleLagrangianGlobal',                   &
                             Default      = .false.,                                    &
                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOrigins - ModuleLagrangianGlobal - ERR1257'


                !AmbientConcentration 
                call GetData(NewProperty%AmbientConcentration,                          &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='AMBIENT_CONC',                              &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1260'


                if (flag == 1) then
                    NewProperty%HaveAmbientConcentration = ON
                else
                    !Checks if the WP module have this property
                    WP_HaveProperty = WaterPropertyExists (                             &
                                            Me%EulerModel(em)%ObjWaterProperties,       &
                                            NewProperty%ID, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1270'

                    if (.not. WP_HaveProperty) then
                        write(*,*)'No Ambient Concentration for the Property : ',trim(NewProperty%Name)
                        write(*,*)'Use the KeyWord                           : AMBIENT_CONC'
                        write(*,*)'OR define the property in the Eulerian Module'
                        stop      'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1280'
                    endif
                endif


                !Extra Stuff
                ret = CheckPropertyName (NewProperty%Name, PropertyID)
                select case (PropertyID)

                !T90
                case (Fecal_Coliforms_,E_Coli_ )

                    Me%State%T90      = ON 

                    NewProperty%T90ON = ON

                    call GetData(NewProperty%T90Variable,                               &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      ='T90_VARIABLE',                          &
                                 ClientModule ='ModuleLagrangianGlobal',                &
                                 Default      = .false.,                                &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1290'


                    if (NewProperty%T90Variable) then
                        
                        if (NewProperty%T90Variable) Me%State%T90Variable = .true.
                        call GetData(NewProperty%T90Var_Method,                         &
                                     Me%ObjEnterData,                                   &
                                     flag,                                              &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      ='T90_VAR_METHOD',                    &
                                     ClientModule ='ModuleLagrangianGlobal',            &
                                    STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1300'

                        if (flag /= 1) then
                            write(*,*)'Keyword T90_VAR_METHOD not defined at origin :',trim(adjustl(NewOrigin%Name))
                            stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR94b'
                        endif

                        if (NewProperty%T90Var_Method /= FromTimeSerie) then
                            Me%State%T90Compute     = .true. 
                            NewProperty%T90Compute  = .true. 
                        endif

                        if (NewProperty%T90Var_Method == FromTimeSerie) then
                        
                             

                           !T90 time serie file name
                            call GetData(NewProperty%T90File,                                   &
                                         Me%ObjEnterData,                                       &
                                         flag,                                                  &
                                         SearchType   = FromBlockInBlock,                       &
                                         keyword      ='T90_FILE',                              &
                                         ClientModule ='ModuleLagrangianGlobal',                &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1302'
                            if (flag == 0) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1303'

                            if (NewOrigin%Default) then
                                call StartTimeSerieInput(NewProperty%TimeSerieT90,              &
                                                         NewProperty%T90File,                   &
                                                         Me%ExternalVar%ObjTime,                &
                                                         STAT = STAT_CALL)

                                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1304'
                            endif

                            !Concentration column
                            call GetData(NewProperty%T90Column,                                &
                                         Me%ObjEnterData,                                       &
                                         flag,                                                  &
                                         SearchType   = FromBlockInBlock,                       &
                                         keyword      ='T90_COLUMN',                            &
                                         ClientModule ='ModuleLagrangianGlobal',                &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1305'
                            if (flag == 0) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1306'

                        endif

                    endif

                    call GetData(NewProperty%T90,                                &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='T90',                            &
                                 ClientModule ='ModuleLagrangianGlobal',         &
                                 Default      = 7200.,                           &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1310'

                    call GetData(NewProperty%T90Name,                            &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='T90_NAME',                       &
                                 ClientModule ='ModuleLagrangianGlobal',         &
                                 Default      = 'T90',                           &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1315'

                    if (flag == 0) then

                        if (PropertyID == Fecal_Coliforms_) NewProperty%T90Name = GetPropertyName(T90_       )

                        if (PropertyID == E_Coli_         ) NewProperty%T90Name = GetPropertyName(T90_E_Coli_)
                    
                    endif

                end select

                call GetData(NewProperty%WaterPartition%ON,                      &
                             Me%ObjEnterData,                                    &
                             flag,                                               &
                             SearchType   = FromBlockInBlock,                    &
                             keyword      ='PARTITION_WATER',                    &
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             Default      =.false.,                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                       &
                    call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1320')

WP:                     if (NewProperty%WaterPartition%ON) then

                    Me%State%Partition = ON
                    NewOrigin%State%Partition     = ON

                    call GetData(NewProperty%WaterPartition%Coefficient,         &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='PARTITION_COEF_WATER',           &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 Default      = 0.90,                            &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                   &
                        call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1330')

                    call GetData(NewProperty%WaterPartition%TransferRate,        &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='PARTITION_RATE_WATER',           &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 Default      = 1e-3,                            &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                   &
                        call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1340')


                    call GetData(NewProperty%WaterPartition%CoupleProp,          &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='PARTITION_COUPLE_WATER',         &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 Default      = 0.,                              &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                   &
                        call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1350')

                endif WP

                call GetData(NewProperty%SedimentPartition%ON,                   &
                             Me%ObjEnterData,                                    &
                             flag,                                               &
                             SearchType   = FromBlockInBlock,                    &
                             keyword      ='PARTITION_SED',                      & 
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             Default      =.false.,                              &
                             STAT         = STAT_CALL)                          
                if (STAT_CALL /= SUCCESS_)                                       &
                    call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1360')

                if (NewProperty%SedimentPartition%ON.and..not.NewOrigin%State%Deposition) &
                    call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1370')


SP:             if (NewProperty%SedimentPartition%ON) then

                    Me%State%Partition = ON
                    NewOrigin%State%Partition     = ON

                    call GetData(NewProperty%SedimentPartition%Coefficient,      &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='PARTITION_COEF_SED',             &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 Default      = 0.98,                            &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                   &
                        call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1380')

                    call GetData(NewProperty%SedimentPartition%TransferRate,     &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='PARTITION_RATE_SED',             &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 Default      = 1e-4,                            &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                   &
                        call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1390')


                    call GetData(NewProperty%SedimentPartition%CoupleProp,       &
                                 Me%ObjEnterData,                                &
                                 flag,                                           &
                                 SearchType   = FromBlockInBlock,                &
                                 keyword      ='PARTITION_COUPLE_SED',           &
                                 ClientModule ='ModuleLagrangianGlobal',               &
                                 Default      = 0.,                              &    
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                   &
                        call SetError (FATAL_, KEYWORD_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1400')

                endif SP

                !This property has an extinction parameter. This parameter can be use 
                !to compute the effect of this property in the light extinction
                call GetData(NewProperty%ExtinctionParameter,                    &
                             Me%ObjEnterData,                                    &
                             flag,                                               &
                             SearchType   = FromBlockInBlock,                    &
                             keyword      ='EXTINCTION_PARAMETER',               &
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             default      = 1.,                                  &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1440'


                !This property is being filter from the water column
                call GetData(NewProperty%Filtration,                             &
                             Me%ObjEnterData,                                    &
                             flag,                                               &
                             SearchType   = FromBlockInBlock,                    &
                             keyword      ='FILTRATION',                         &
                             ClientModule ='ModuleLagrangianGlobal',                   &
                             default      = .false.,                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1450'

                if (NewProperty%Filtration) NewOrigin%Filtration = .true. 

                call GetData(NewProperty%MinValue,                               &
                             Me%ObjEnterData,                                    &
                             flag,                                               &
                             SearchType   = FromBlockInBlock,                    &
                             keyword      ='MIN_VALUE',                          &
                             ClientModule ='ModuleLagrangianGlobal',             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1452'

                if (flag > 0) then
                    NewProperty%MinON    = .true. 
                else
                    NewProperty%MinON    = .false. 
                    NewProperty%MinValue = FillValueReal
                endif


                call GetData(NewProperty%MaxValue,                               &
                             Me%ObjEnterData,                                    &
                             flag,                                               &
                             SearchType   = FromBlockInBlock,                    &
                             keyword      ='MAX_VALUE',                          &
                             ClientModule ='ModuleLagrangianGlobal',             &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1452'

                if (flag > 0) then
                    NewProperty%MaxON    = .true. 
                else
                    NewProperty%MaxON    = .false. 
                    NewProperty%MaxValue = -FillValueReal
                endif

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
            call SetError (FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1460')
        endif
        if (NewOrigin%State%Deposition .and. .not. SedimentDefined)              &
            call SetError (FATAL_, INTERNAL_, 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1470')



        !If neceassary Starts The WQM module for this origin
        if (NewOrigin%State%WQM) then
            
            !WQM Data File
            call GetData(NewOrigin%WQM_DataFile,                                 &
                         Me%ObjEnterData,                                        &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='WQM_DATA_FILE',                          &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1480'

            call StartWaterQuality(NewOrigin%WaterQualityID,                     &
                                   NewOrigin%WQM_DataFile,                       &
                                   STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1490'

            NewOrigin%NextWQMCompute = Me%ExternalVar%Now
            call GetDTWQM (NewOrigin%WaterQualityID, DTSecond = NewOrigin%DTWQM, &
                           STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1500'

        endif

        !If neceassary Starts The Oil module for this origin
        if (NewOrigin%State%Oil) then

            call GetData(NewOrigin%UseTheoricArea,                               &
                         Me%ObjEnterData,                             &
                         flag,                                                   &
                         SearchType   = FromBlock,                               &
                         keyword      ='THEORIC_AREA',                           &
                         ClientModule ='ModuleLagrangianGlobal',                       &
                         Default      = .false.,                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1510'

            call ConstructParticOil (NewOrigin, ClientNumber)

        endif

        call GetData(NewOrigin%State%Age,                                        &
                     Me%ObjEnterData,                                 &
                     flag,                                                       &
                     SearchType   = FromBlock,                                   &
                     keyword      ='COMPUTE_AGE',                                &
                     ClientModule ='ModuleLagrangianGlobal',                           &
                     Default      = .false.,                                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1520'

        if (NewOrigin%State%Age) Me%State%Age = .true.


        if (NewOrigin%State%ComputePlume .and. .not. NewOrigin%Default) then

            call Construct_Jet(JetID       = NewOrigin%Movement%ObjJet,          &
                               FileName    = NewOrigin%Movement%JetDataFile,     &
                               PositionX   = NewOrigin%Position%CartX,           & 
                               PositionY   = NewOrigin%Position%CartY,           & 
                               Flow        = NewOrigin%Flow,                     &
                               Salinity    = NewOrigin%Movement%JetSalinity,     &
                               Temperature = NewOrigin%Movement%JetTemperature,  &
                               STAT        = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1530'

            call UnitsManager(NewOrigin%Movement%JetUnit, OPEN_FILE,             &
                              STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1540'
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
                stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1560'
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
                     ClientModule ='ModuleLagrangianGlobal',                          &
                     Default      = OFF,                                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1570'



        call GetData(NewOrigin%AveragePositionON,                               &
                     Me%ObjEnterData,                                           &
                     flag,                                                      &
                     SearchType   = FromBlock,                                  &
                     keyword      ='COMPUTE_AVERAGE_POSITION',                  &
                     ClientModule ='ModuleLagrangianGlobal',                          &
                     Default      = OFF,                                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1580'

        if (NewOrigin%AveragePositionON) then

            call GetData(NewOrigin%CoefRadius,                                  &
                         Me%ObjEnterData,                                       &
                         flag,                                                  &
                         SearchType   = FromBlock,                              &
                         keyword      ='COEF_RADIUS',                           &
                         ClientModule ='ModuleLagrangianGlobal',                      &
                         !90% of the particles inside the circle
                         Default      = 1.645,                                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOneOrigin - ModuleLagrangianGlobal - ERR1590'

        endif
        


    end subroutine ConstructOneOrigin

    !--------------------------------------------------------------------------

    logical function CheckOriginInLandCell ()


        !Local-----------------------------------------------------------------
        type (T_Position)                           :: Position
        character (len = StringLength)              :: OriginName
        real, dimension(1:2)                        :: Aux1D
        character(3)                                :: Aux
        integer                                     :: flag
        integer                                     :: STAT_CALL
        integer                                     :: i, j, k, em
        logical                                     :: Default, HaveOrigin, NoDomain

        !Begin-----------------------------------------------------------------


        CheckOriginInLandCell = .false. 

        !Gets its name    
        call GetData(OriginName,                                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromBlock,                                          &
                     keyword      ='ORIGIN_NAME',                                       &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'CheckOriginInLandCell - ModuleLagrangianGlobal - ERR10'

        if (flag == 0) then
            Aux = ' '
            write(Aux,'(i3)') Me%nOrigins + 1
            OriginName = trim(adjustl('Origin_'//trim(adjustl(Aux))))
        end if

        call GetData(Default,                                                           &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromBlock,                                          &
                     keyword      ='DEFAULT',                                           &
                     ClientModule ='ModuleLagrangianGlobal',                            &
                     Default      = OFF,                                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckOriginInLandCell - ModuleLagrangianGlobal - ERR20'
           
NDF:    if (.not. Default) then
    
            !Horizontal position in coordinates X, Y
            call GetData(Aux1D,                                                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromBlock,                                      &
                         keyword      ='POSITION_COORDINATES',                          &
                         ClientModule ='ModuleLagrangianGlobal',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOriginInLandCell - ModuleLagrangianGlobal - ERR30'

            if (flag == 2) then
                HaveOrigin = .true.
            endif

            em = Locate_ModelDomain(Aux1D(1), Aux1D(2), NoDomain) 

            if (NoDomain) then
                write(*,*) 'Discharge outside the domain - ',trim(OriginName),' - ',trim(Me%EulerModel(Me%EulerModelNumber)%Name)
                write (*,*) 'Origin ',trim(OriginName),' is outside of the outer model domain'
                stop 'CheckOriginInLandCell - ModuleLagrangianGlobal - ERR40'
            endif

            Position%CoordX = Aux1D(1)
            Position%CoordY = Aux1D(2)

            call Convert_XY_CellIJ(Me%EulerModel(em), Position, Referential = GridCoord_)
                
            !Checks if a valid horizontal origin was found
            if (.not. HaveOrigin) then

                write(*,*)'No Valid Horizontal Location defined for ',trim(adjustl(OriginName))
                stop 'CheckOriginInLandCell - ModuleLagrangianGlobal - ERR50'

            endif

            !Verifies the horizontal location of the origin
            i = Position%I
            j = Position%J
            k = Me%EulerModel(em)%WorkSize%KUB

            if (Me%EulerModel(em)%Waterpoints3D(i, j, k) /= WaterPoint) then
                write(*,*) 'Discharge in a land cell I=',i,' J=',j,'Model name=',trim(Me%EulerModel(em)%Name)
                CheckOriginInLandCell = .true.
            endif

            

        endif NDF

    end function CheckOriginInLandCell

    !--------------------------------------------------------------------------


    subroutine BoxTypeVariablesDefiniton

    !Local---------------------------------------------------------------------
    integer                 :: flag, STAT_CALL, em
    !Begin---------------------------------------------------------------------

        !Boxes data file
        call GetData(Me%Files%BoxDataFile,                                              &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='PARTIC_BOX',                                        &
                     ClientModule ='ModuleLagrangianGlobal',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR120'
        if (flag == 1) then


            Me%State%BoxDif = .true.
            
        endif

        !Boxes data file
        call GetData(Me%Files%MonitorBox,                                               &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='MONITOR_BOX',                                       &
                     ClientModule ='ModuleLagrangianGlobal',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR140'
        if (flag == 1) then

            Me%State%Monitor = .true.


            !Name of the Property
            call GetData(Me%MonitorProperty,                                       &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='MONITOR_BOX_PROP_MASS',                         &
                         ClientModule ='ModuleLagrangianGlobal',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR160'

            if (flag == 1) then

                !Checks if the property is recognized by the model       
                if (.not. CheckPropertyName(trim(Me%MonitorProperty), Me%MonitorPropertyID)) then
                    write (*,*)'Unknown Property : ', trim(Me%MonitorProperty)
                    stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR170'
                end if

                Me%State%MonitorPropMass = .true.
            end if

        endif

        !Eulerian boxes data file
        call GetData(Me%Files%EulerianMonitorBox,                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='EULERIAN_MONITOR_BOX',                              &
                     ClientModule ='ModuleLagrangianGlobal',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR180'
        
        if (flag == 1) Me%State%EulerianMonitor = .true.

        !Associates Beaching Probabilities 
        call GetData(Me%State%AssociateBeachProb,                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='ASSOCIATE_BEACH_PROB',                              &
                     Default      = .false.,                                            &
                     ClientModule ='ModuleLagrangianGlobal',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR190'

        if (Me%State%AssociateBeachProb) then

            !Maximum distance between particles and coast for particle beaching 
            call GetData(Me%BeachingLimit,                                              &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='BEACHING_LIMIT',                                &
                         Default      = 5.0,                                            &
                         ClientModule ='ModuleLagrangianGlobal',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR200'

            !Outbox Beaching Probability
            call GetData(Me%DefaultBeachingProbability,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='DEFAULT_BEACHING_PROB',                         &
                         Default      = 0.5,                                            &
                         ClientModule ='ModuleLagrangianGlobal',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR210'


            !Beaching Probabilities Boxes data file
            call GetData(Me%Files%BeachingBoxFileName,                                  &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = FromFile,                                       &
                         keyword      ='BEACHING_BOX_FILENAME',                         &
                         ClientModule ='ModuleLagrangianGlobal',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR220'
            if (flag == 1) then

                Me%State%HaveBeachingProbBox = .true.
        
            endif

        endif

d1:     do em = 1, Me%EulerModelNumber

            if (Me%State%BoxDif) then
                        !Starts BoxDif
                call StartBoxDif(BoxDifID         = Me%EulerModel(em)%ObjBoxDif,        &
                                 TimeID           = Me%EulerModel(em)%ObjTime,          &
                                 HorizontalGridID = Me%EulerModel(em)%ObjHorizontalGrid,&
                                 BoxesFilePath    = Me%Files%BoxDataFile,               &
                                 WaterPoints3D    = Me%EulerModel(em)%Waterpoints3D,    &
                                 Size3D           = Me%EulerModel(em)%Size,             &
                                 WorkSize3D       = Me%EulerModel(em)%WorkSize,         &
                                 STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR130'

            endif

            if (Me%State%Monitor) then

                !Starts BoxDif
                call StartBoxDif(BoxDifID         = Me%EulerModel(em)%ObjMonBox,                           & 
                                 TimeID           = Me%EulerModel(em)%ObjTime,                             &
                                 HorizontalGridID = Me%EulerModel(em)%ObjHorizontalGrid,                   &
                                 BoxesFilePath    = Me%Files%MonitorBox,                    &
                                 WaterPoints3D    = Me%EulerModel(em)%Waterpoints3D,           &
                                 Size3D           = Me%EulerModel(em)%Size,                    &
                                 WorkSize3D       = Me%EulerModel(em)%WorkSize,                &
                                 STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR150'

            endif


            if (Me%State%HaveBeachingProbBox) then

                !Starts BoxDif
                call StartBoxDif(BoxDifID         = Me%EulerModel(em)%ObjBeachingProbBox,              &
                                 TimeID           = Me%EulerModel(em)%ObjTime,                         &
                                 HorizontalGridID = Me%EulerModel(em)%ObjHorizontalGrid,               &
                                 BoxesFilePath    = Me%Files%BeachingBoxFileName,       &
                                 WaterPoints3D    = Me%EulerModel(em)%Waterpoints3D,       &
                                 Size3D           = Me%EulerModel(em)%Size,                &
                                 WorkSize3D       = Me%EulerModel(em)%WorkSize,            &
                                 STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'BoxTypeVariablesDefiniton - ModuleLagrangianGlobal - ERR230'

                Me%State%HaveBeachingProbBox = .true.
        
            endif

        enddo d1

    end subroutine BoxTypeVariablesDefiniton

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
            if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR10'

            call ExtractBlockFromBuffer(Me%ObjEnterDataClone,                           &
                                        ClientNumberClone,                              &
                                        block_begin_Clone, block_end_Clone,             &
                                        BlockFound,                                     &
                                        STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR20'

i2:         if (BlockFound) then
    
                !Construct enter data original again
                call ConstructEnterData(Me%ObjEnterDataOriginal, Me%Files%ConstructData, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR30'

                PreviousStart  = FillValueInt

                call RewindBuffer(Me%ObjEnterDataClone,                                 &
                                  STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR40'


            else  i2

                ClonesExist = .false. 

                !Finished reading block -> unlocks block reading
                call Block_Unlock(Me%ObjEnterDataClone,                                 &
                                  ClientNumberClone,                                    &
                                  STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR50'

                call KillEnterData(Me%ObjEnterDataClone, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR60'

            endif i2

        endif i1

i3:     if (ClonesExist) then

d2:     do
 
            call ExtractBlockFromBuffer(Me%ObjEnterDataClone,                           &
                                        ClientNumberClone,                              &
                                        block_begin_Clone, block_end_Clone,             &
                                        BlockFound,                                     &
                                        STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR70'


i4:         if (BlockFound) then 

                call GetBlockSize(Me%ObjEnterDataClone, ClientNumberClone, StartLine, EndLine, STAT= STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR80'

i5:             if (StartLine > PreviousStart) then

                    call GetData(OriginNameClone,                                       &
                                 Me%ObjEnterDataClone,                                  &
                                 flag,                                                  &
                                 SearchType   = FromBlock,                              &
                                 keyword      ='CLONE',                                 &
                                 ClientModule ='ModuleLagrangianGlobal',                      &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR90'

                    if (flag == 0) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR100'


i6:                if (trim(OriginNameClone) == trim(OriginName)) then
                        FoundCloneOrigin = .true.

                        call ReplaceKeywordsInClone(ClientNumber, ClientNumberClone, StartLine, EndLine)

                        PreviousStart = StartLine

                        exit

                    endif i6

                endif i5

            else i4

                PreviousStart  = FillValueInt

                call RewindBuffer(Me%ObjEnterDataClone, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CheckForOriginClones - ModuleLagrangianGlobal - ERR110'

                FoundCloneOrigin = .false.

                exit

            endif i4

        enddo d2

        endif i3

    end subroutine CheckForOriginClones

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine ReplaceKeywordsInClone(ClientNumber, ClientNumberClone, StartClone, EndClone)


        !Arguments-------------------------------------------------------------
        integer                                     :: ClientNumber, ClientNumberClone, StartClone, EndClone

        !Local-----------------------------------------------------------------
        character(len=PathLength)                   :: FullBufferLine
        character(len=StringLength)                 :: KeywordClone, KeywordOriginal
        logical                                     :: ExistKeywordClone, ExistKeywordOriginal
        logical                                     :: BlockInBlockOn, BlockPropON
        integer                                     :: StartLine, EndLine
        integer                                     :: lineC, lineO
        integer                                     :: STAT_CALL


        !Begin-----------------------------------------------------------------

        call GetBlockSize(Me%ObjEnterData, ClientNumber, StartLine, EndLine, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangianGlobal - ERR10'

        BlockPropON    = .false.
        BlockInBlockOn = .false. 

d1:     do lineC = StartClone, EndClone

            call GetFullBufferLine    (Me%ObjEnterDataClone, lineC, FullBufferLine, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangianGlobal - ERR20'

            if (trim(FullBufferLine) == property_begin) then
                BlockInBlockOn = .true.
                BlockPropON    = .true.
            endif

            if (trim(FullBufferLine) == property_end) then 
                BlockInBlockOn = .false.
                cycle 
            endif

            if (BlockInBlockOn) cycle

            call GetKeywordFromLine(Me%ObjEnterDataClone, lineC, ExistKeywordClone, KeywordClone, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangianGlobal - ERR20'

i1:         if (ExistKeywordClone) then
                             
i2:             if (trim(KeywordClone) /= 'CLONE') then

d2:                 do lineO = StartLine, EndLine                

                        call GetKeywordFromLine(Me%ObjEnterDataOriginal, lineO,         &
                                                ExistKeywordOriginal, KeywordOriginal, STAT = STAT_CALL)  
                        if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangianGlobal - ERR30'
                    
i3:                     if (ExistKeywordOriginal) then

i4:                         if (trim(KeywordClone) == trim(KeywordOriginal)) then

                                call GetFullBufferLine    (Me%ObjEnterDataClone, lineC, FullBufferLine, STAT = STAT_CALL)  
                                if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangianGlobal - ERR40'

                                call ReplaceFullBufferLine(Me%ObjEnterData,      lineO, FullBufferLine, STAT = STAT_CALL)  
                                if (STAT_CALL /= SUCCESS_) stop 'ReplaceKeywordsInClone - ModuleLagrangianGlobal - ERR50'

                                exit

                            endif i4
                        endif i3
                    enddo d2
                endif i2
            endif i1

        enddo d1
        
        if (BlockPropON) call ReplacePropInClone(ClientNumber, ClientNumberClone)


    end subroutine ReplaceKeywordsInClone


    !--------------------------------------------------------------------------

    subroutine ReplacePropInClone(ClientNumber, ClientNumberClone)


        !Arguments-------------------------------------------------------------
        integer                                     :: ClientNumber, ClientNumberClone

        !Local-----------------------------------------------------------------
        character(len=PathLength)                   :: FullBufferLine
        character(len=StringLength)                 :: PropNameClone, PropName, KeywordClone, KeywordOriginal
        logical                                     :: ExistKeywordClone, ExistKeywordOriginal
        logical                                     :: PropertyCloneFound, PropertyFound
        integer                                     :: StartLine, EndLine
        integer                                     :: StartClone, EndClone
        integer                                     :: lineC, lineO, flag
        integer                                     :: STAT_CALL


        !Begin-----------------------------------------------------------------

d1:     do 
 
            call ExtractBlockFromBlock(Me%ObjEnterDataClone, ClientNumberClone,         &
                                       property_begin, property_end,                    &
                                       PropertyCloneFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR10'
            
i1:         if (PropertyCloneFound) then

                call GetData(PropNameClone,                                             &
                             Me%ObjEnterDataClone,                                      &
                             flag,                                                      &
                             SearchType   = FromBlockinBlock,                           &
                             keyword      ='NAME',                                      &
                             ClientModule ='ModuleLagrangianGlobal',                    &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR20'

                if (flag == 0) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR30'


                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
            
                if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR140'

d2:             do 
 
                    call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,           &
                                               property_begin, property_end,            &
                                               PropertyFound, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR30'

                    if (.not. PropertyFound) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR40'

                    call GetData(PropName,                                              &
                                 Me%ObjEnterData,                                       &
                                 flag,                                                  &
                                 SearchType   = FromBlockinBlock,                       &
                                 keyword      ='NAME',                                  &
                                 ClientModule ='ModuleLagrangianGlobal',                &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR50'

                    if (flag == 0) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR60'


i2:                 if (trim(PropNameClone) == trim(PropName)) then

                        call GetBlockSize(Me%ObjEnterData, ClientNumber,                &
                                          StartLine, EndLine, SearchType = FromBlockInBlock, STAT = STAT_CALL)  
                        if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR70'

                        call GetBlockSize(Me%ObjEnterDataClone, ClientNumberClone,      &
                                          StartClone, EndClone, SearchType = FromBlockinBlock, STAT = STAT_CALL)  
                        if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR80'

d3:                     do lineC = StartClone, EndClone

                            call GetFullBufferLine    (Me%ObjEnterDataClone, lineC, FullBufferLine, STAT = STAT_CALL)  
                            if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR90'

                            call GetKeywordFromLine(Me%ObjEnterDataClone, lineC, &
                                                    ExistKeywordClone, KeywordClone, STAT = STAT_CALL)  
                            if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR100'

i3:                         if (ExistKeywordClone) then
                             
i4:                             if (trim(KeywordClone) /= 'NAME') then

d4:                                 do lineO = StartLine, EndLine                

                                        call GetKeywordFromLine(Me%ObjEnterDataOriginal, lineO,         &
                                                                ExistKeywordOriginal, KeywordOriginal, STAT = STAT_CALL)  
                                        if (STAT_CALL /= SUCCESS_) stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR110'
                    
i5:                                     if (ExistKeywordOriginal) then

i6:                                         if (trim(KeywordClone) == trim(KeywordOriginal)) then

                                                call GetFullBufferLine    (Me%ObjEnterDataClone, lineC, &
                                                                           FullBufferLine, STAT = STAT_CALL)
                                                if (STAT_CALL /= SUCCESS_)              &
                                                    stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR120'

                                                call ReplaceFullBufferLine(Me%ObjEnterData,      lineO, &
                                                                           FullBufferLine, STAT = STAT_CALL)  
                                                if (STAT_CALL /= SUCCESS_)              &
                                                    stop 'ReplacePropInClone - ModuleLagrangianGlobal - ERR130'

                                                exit

                                            endif i6
                                        endif i5
                                    enddo d4
                                endif i4
                            endif i3

                        enddo d3
                        
                        exit

                    endif i2
                    
                enddo d2


            else i1

                exit

            endif i1

        enddo d1

        

    end subroutine ReplacePropInClone
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

    subroutine ChangeOriginOnline(Nmax, NFirstExtraction)


        !Arguments-------------------------------------------------------------
        integer                                     :: Nmax, NFirstExtraction

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: NewOrigin
        character(Len=PathLength)                   :: String
        integer                                     :: i, j, em
        logical                                     :: NoDomain


        !Begin-----------------------------------------------------------------

!        if (Me%Online%ON) then

        Me%Files%TransientHDF =trim(Me%Files%TransientHDF)//'_'//trim(Me%Online%TimeStamp)//'.hdf'

        NewOrigin=> Me%FirstOrigin            
        
        do i=1, Nmax
        do j=1, NFirstExtraction

            NewOrigin%Position%ModelID =  Locate_ModelDomain(Me%Online%X(i), Me%Online%Y(i), NoDomain) 
            if (NoDomain) then
                stop 'ChangeOriginOnline - ModuleLagrangianGlobal - ERR10'
            endif

            em = NewOrigin%Position%ModelID

            NewOrigin%Position%CoordX = Me%Online%X(i)
            NewOrigin%Position%CoordY = Me%Online%Y(i)

            call Convert_XY_CellIJ(Me%EulerModel(em), NewOrigin%Position, Referential = GridCoord_)

            !Change Origin name
            write(String,*) i
            NewOrigin%Name = trim(NewOrigin%Name)//'_'//trim(adjustl(String))

            Me%Online%StartDate(i,6) = Me%Online%StartDate(i,6)

            call SetDate (NewOrigin%StartEmission, Me%Online%StartDate(i,1), Me%Online%StartDate(i,2), &
                                                   Me%Online%StartDate(i,3), Me%Online%StartDate(i,4), &
                                                   Me%Online%StartDate(i,5), Me%Online%StartDate(i,6))


            NewOrigin%StopEmission = NewOrigin%StartEmission

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
            stop 'ConstructEmissionType - ModuleLagrangianGlobal - ERR01'

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
            stop 'ConstructEmissionType - ModuleLagrangianGlobal - ERR01'

        end select

        
        !Verifies consistency
        if (NewOrigin%EmissionTemporal == Continuous_ .and. NewOrigin%EmissionSpatial == Accident_) then
            write(*,*)'Cant emit particle in continous as accident'
            write(*,*)'Origin :',trim(adjustl(NewOrigin%Name))
            stop 'ConstructEmissionType - ModuleLagrangianGlobal - ERR02'
        endif

        !If BoxDif isnt associated keyword Partic PARTIC_BOX wasnt found
        if (NewOrigin%EmissionSpatial == Box_ .and. Me%EulerModel(1)%ObjBoxDif == 0) then
            write(*,*)'Cant emit particle in a box, once keyword PARTIC_BOX wasnt given'
            write(*,*)'Origin :',trim(adjustl(NewOrigin%Name))
            stop 'ConstructEmissionType - ModuleLagrangianGlobal - ERR03'
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
            if (CurrentOrigin%EmissionTemporal == Continuous_) then

                call ConstructEmissionTime          (CurrentOrigin)

            endif

            if (CurrentOrigin%EmissionTemporal == Instantaneous_) then

                CurrentOrigin%OnlyOnceEmit = .true.

            endif
    
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

        call ParticleEmission ()

    end subroutine ConstructEmission

    !--------------------------------------------------------------------------

    subroutine ConstructOverlay

        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        integer                                     :: SizeILB, SizeIUB
        integer                                     :: SizeJLB, SizeJUB
        integer                                     :: SizeKLB, SizeKUB, em

        !Begin-----------------------------------------------------------------
        
em1:    do em =1, Me%EulerModelNumber

            !Shorten Size
            SizeILB     = Me%EulerModel(em)%Size%ILB   
            SizeIUB     = Me%EulerModel(em)%Size%IUB   
            SizeJLB     = Me%EulerModel(em)%Size%JLB   
            SizeJUB     = Me%EulerModel(em)%Size%JUB   
            SizeKLB     = Me%EulerModel(em)%Size%KLB   
            SizeKUB     = Me%EulerModel(em)%Size%KUB   
        
            if (Me%OverLay) then

                allocate (Me%EulerModel(em)%OverLay%VelUFinal (SizeILB:SizeIUB, SizeJLB:SizeJUB, SizeKLB:SizeKUB))
                allocate (Me%EulerModel(em)%OverLay%VelVFinal (SizeILB:SizeIUB, SizeJLB:SizeJUB, SizeKLB:SizeKUB))

            else
        
                nullify (Me%EulerModel(em)%OverLay%VelUFinal)
                nullify (Me%EulerModel(em)%OverLay%VelVFinal)

            endif

        enddo em1


    end subroutine ConstructOverlay

    !--------------------------------------------------------------------------

    subroutine ConstructMonitoring 

        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        type (T_EulerModel), pointer                    :: EulerModel
        integer                                         :: nB
        integer                                         :: NumberOfOrigins, iO
        integer                                         :: STAT_CALL
        character(StringLength), dimension(:), pointer  :: PropertyList
        character(StringLength)                         :: AuxChar
        type (T_Property), pointer                      :: CurrentProperty
        integer                                         :: iProp
        integer                                         :: ILB, IUB
        integer                                         :: JLB, JUB
        integer                                         :: KLB, KUB, em, NumberOfBoxes
 
        !Begin-----------------------------------------------------------------

d1:     do em =1, Me%EulerModelNumber
        
            EulerModel => Me%EulerModel(em)
        
            if (Me%State%Monitor) then

                call GetNumberOfBoxes(EulerModel%ObjMonBox, NumberOfBoxes3D = EulerModel%Monitor%NumberOfBoxes, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMonitoring - ModuleLagrangianGlobal - ERR01'

                NumberOfOrigins = Me%nOrigins

                NumberOfBoxes   = EulerModel%Monitor%NumberOfBoxes

                allocate (EulerModel%Monitor%InstBoxVolume           (NumberOfBoxes                  ))
                allocate (EulerModel%Monitor%InstVolumeByOrigin      (NumberOfBoxes, NumberOfOrigins ))
                allocate (EulerModel%Monitor%IntgBoxVolume           (NumberOfBoxes                  ))
                allocate (EulerModel%Monitor%IntgVolumeByOrigin      (NumberOfBoxes, NumberOfOrigins ))
                allocate (EulerModel%Monitor%InstBoxMass             (NumberOfBoxes                  ))
                allocate (EulerModel%Monitor%InstMassByOrigin        (NumberOfBoxes, NumberOfOrigins ))
                allocate (EulerModel%Monitor%NumberOfCellsPerBox     (NumberOfBoxes                  ))
                allocate (EulerModel%Monitor%NumberOfCellsFromOrigin (NumberOfBoxes, NumberOfOrigins ))
                allocate (EulerModel%Monitor%ObjTimeSerie            (NumberOfBoxes                  ))            

                EulerModel%Monitor%IntgBoxVolume      = 0.
                EulerModel%Monitor%IntgVolumeByOrigin = 0.
                EulerModel%Monitor%ObjTimeSerie       = 0
            



                If (Me%State%MonitorPropMass) then
                    allocate (PropertyList (3*(NumberOfOrigins + 1)))
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

                    End If

                    write(AuxChar, fmt=*)nB
                    call StartTimeSerie (EulerModel%Monitor%ObjTimeSerie(nB),                   &
                                         EulerModel%ObjTime,                                    &
                                         Me%Files%ConstructData,                                &
                                         PropertyList,                                          &
                                         Extension   = "LMB",                                   &
                                         ResultFileName = "MonitorBox_"//trim(adjustl(AuxChar)),&
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMonitoring - ModuleLagrangianGlobal - ERR02'

                enddo
                deallocate (PropertyList)

            else

                nullify  (EulerModel%Monitor%InstBoxVolume            )
                nullify  (EulerModel%Monitor%InstVolumeByOrigin       )
                nullify  (EulerModel%Monitor%IntgBoxVolume            )
                nullify  (EulerModel%Monitor%IntgVolumeByOrigin       )
                nullify  (EulerModel%Monitor%InstBoxMass              )
                nullify  (EulerModel%Monitor%InstMassByOrigin         )
                nullify  (EulerModel%Monitor%NumberOfCellsPerBox      )   
                nullify  (EulerModel%Monitor%NumberOfCellsFromOrigin  )
                nullify  (EulerModel%Monitor%ObjTimeSerie             )

            endif

            if (Me%State%EulerianMonitor) then

                ILB     = EulerModel%Size%ILB   
                IUB     = EulerModel%Size%IUB   
                JLB     = EulerModel%Size%JLB   
                JUB     = EulerModel%Size%JUB   
                KLB     = EulerModel%Size%KLB   
                KUB     = EulerModel%Size%KUB   

                !This test is done for simply reason
                if (Me%nGroups > 1) then
                    write(*,*)'Cannot write Eulerian Box Time Series for simulation with more then one Group'
                    stop 'ConstructMonitoring - ModuleLagrangianGlobal - ERR40'
                endif

                !Counts Number of properties
                iProp = 0
                CurrentProperty => Me%OriginDefault%FirstProperty
                do while (associated(CurrentProperty))
                    iProp = iProp + 1    
                    CurrentProperty => CurrentProperty%Next
                enddo

                allocate (PropertyList(1:iProp))

                !Counts Number of properties
                iProp = 0
                CurrentProperty => Me%OriginDefault%FirstProperty
                do while (associated(CurrentProperty))
                    iProp               = iProp + 1
                    PropertyList(iProp) = "Lagrangian_"//trim(CurrentProperty%Name)
                    CurrentProperty => CurrentProperty%Next
                enddo

                allocate(EulerModel%EulerianMonitor%Mass(ILB:IUB, JLB:JUB, KLB:KUB))

                call SetMatrixValue(EulerModel%EulerianMonitor%Mass, EulerModel%Size, &
                                    dble(0.), EulerModel%Waterpoints3D)
            
                !Starts BoxDif
                call StartBoxDif(BoxDifID         = EulerModel%ObjEulMonBox,                        &
                                 TimeID           = EulerModel%ObjTime,                             &
                                 HorizontalGridID = EulerModel%ObjHorizontalGrid,                   &
                                 BoxesFilePath    = Me%Files%EulerianMonitorBox,            &
                                 ScalarOutputList = PropertyList,                           &
                                 WaterPoints3D    = EulerModel%Waterpoints3D,           &
                                 Size3D           = EulerModel%Size,                    &
                                 WorkSize3D       = EulerModel%WorkSize,                &
                                 STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMonitoring - ModuleLagrangianGlobal - ERR050'
            
                deallocate (PropertyList)

            end if

        enddo d1

        nullify(EulerModel)


    end subroutine ConstructMonitoring

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len = StringLength), dimension(:), pointer    :: PropertyList
        integer, pointer, dimension(:, :, :)                    :: WaterPoints3D
        type (T_Property), pointer                              :: CurrentProperty
        integer                                                 :: iProp
        real                                                    :: CoordX, CoordY
        integer                                                 :: STAT_CALL
        integer                                                 :: iflag, em, dn, Id, Jd, TimeSerieNumber
        character(len=PathLength)                               :: TimeSerieLocationFile
        logical                                                 :: CoordON, IgnoreOK


        !This test is done for simply reason
        if (Me%nGroups > 1) then
            write(*,*)'Cannot write Time Series for simulation with more then one Group'
            stop 'ConstructTimeSeries - ModuleLagrangianGlobal - ERR10'
        endif

        !Counts Number of properties
        iProp = 0
        CurrentProperty => Me%OriginDefault%FirstProperty
        do while (associated(CurrentProperty))
            if (CurrentProperty%WritesTimeSerie) then
                iProp = iProp + 1
                if (CurrentProperty%T90ON) iProp = iProp + 1
            endif
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
        CurrentProperty => Me%OriginDefault%FirstProperty
        do while (associated(CurrentProperty))
            if (CurrentProperty%WritesTimeSerie) then
                iProp = iProp + 1    
                PropertyList(iProp) = trim(adjustl(CurrentProperty%Name))
                if (CurrentProperty%T90ON) then
                    iProp = iProp + 1
                    PropertyList(iProp) = trim(adjustl(CurrentProperty%T90Name))
                endif
            endif
            CurrentProperty => CurrentProperty%Next
        enddo

d1:     do em = 1, Me%EulerModelNumber

            !Gets the position of the water points in the Map Module
            call GetWaterPoints3D(Me%EulerModel(em)%ObjMap, WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR20'
        
            call GetData(TimeSerieLocationFile,                                         &
                         Me%ObjEnterData ,iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'TIME_SERIE_LOCATION',                          &
                         ClientModule = 'ModuleLagrangianGlobal',                       &
                         Default      = Me%Files%ConstructData,                         &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR30' 

            !Constructs TimeSerie
            call StartTimeSerie(Me%EulerModel(em)%ObjTimeSerie, Me%EulerModel(em)%ObjTime,&
                                TimeSerieLocationFile,                                  &
                                PropertyList, "srl",                                    &
                                WaterPoints3D = WaterPoints3D,                          &
                                ModelName     = Me%EulerModel(em)%Name,                 &
                                STAT = STAT_CALL)
            if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR40'

            !Ungets the WaterPoints 
            call UnGetMap (Me%EulerModel(em)%ObjMap, WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR50'

            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%EulerModel(em)%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR60'

            do dn = 1, TimeSerieNumber

                call GetTimeSerieLocation(Me%EulerModel(em)%ObjTimeSerie, dn,               &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR70'

                if (CoordON) then
                    call GetXYCellZ(Me%EulerModel(em)%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR80'

                    if (Id < 0 .or. Jd < 0) then
                
                        call TryIgnoreTimeSerie(Me%EulerModel(em)%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR90'

                        if (IgnoreOK) then
                            cycle
                        else
                            stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR100'
                        endif

                    endif


                    call CorrectsCellsTimeSerie(Me%EulerModel(em)%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleLagrangianGlobal - ERR110'
                endif

            enddo


        enddo d1

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
        integer                                     :: iO, i
        logical                                     :: HaveID
          
        !Checks how many different Groups exists
        Me%nGroups = 0
        allocate (AuxGroups    (Me%nOrigins))
        allocate (AuxGroupsIDs (Me%nOrigins))
        CurrentOrigin => Me%FirstOrigin
        iO            =  1
        if (.not. associated(Me%FirstOrigin)) Me%nGroups = 1
         
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


        !Origins per Group
        allocate   (Me%nOriginsGroup  (Me%nGroups))
        Me%nOriginsGroup = 0



        !Group IDs
        allocate   (Me%GroupIDs       (Me%nGroups))
 
        if (associated(Me%FirstOrigin)) then

            do i = 1, Me%nGroups
                Me%GroupIDs (i) = AuxGroupsIDs(i)
            enddo

        else
            Me%GroupIDs (1) = 1
        endif

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
                            stop 'VerifyOriginProperties - ModuleLagrangianGlobal - ERR20'
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
                            stop 'MergeOldWithNewOrigins - ModuleLagrangianGlobal - ERR02'
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
        type (T_EulerModel),  pointer               :: EulerModel
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
        integer                                     :: STAT_CALL, BoxCell, emBox, em
        logical                                     :: BoxInsideDomain

        !Begin------------------------------------------------------------------

        emBox = Me%EulerModelNumber


        do em = 1, Me%EulerModelNumber

            BoxInsideDomain = GetIfBoxInsideDomain(Me%EulerModel(em)%ObjBoxDif, CurrentOrigin%BoxNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'EmissionBox - ModuleLagrangianGlobal - ERR10'

            if (BoxInsideDomain) then
                emBox = em
                exit
            endif

        enddo

        EulerModel => Me%EulerModel(emBox)

        CurrentOrigin%Position%ModelID = emBox

        
        ILB = EulerModel%WorkSize%ILB
        IUB = EulerModel%WorkSize%IUB

        JLB = EulerModel%WorkSize%JLB
        JUB = EulerModel%WorkSize%JUB

        KLB = EulerModel%WorkSize%KLB
        KUB = EulerModel%WorkSize%KUB



        !Get the boxes
        call GetBoxes(EulerModel%ObjBoxDif, Boxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'EmissionBox - ModuleLagrangianGlobal - ERR20'

                       
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

OP1:            if (EulerModel%OpenPoints3D(i, j, KUB) == OpenPoint) then

                if (CurrentOrigin%Movement%Float) then
                
                    BoxCell     = Boxes(i, j, KUB)

                    DepthPartic = EulerModel%SZZ(i, j, KUB)

                else if (CurrentOrigin%Deposition%BottomEmission) then

                    BoxCell     = Boxes(i, j, EulerModel%KFloor(i, j))


                    DepthPartic = EulerModel%SZZ(i, j, EulerModel%KFloor(i, j)-1) +               &
                                               CurrentOrigin%Deposition%BottomDistance / 2.

                endif

BOX2D:          if (BoxCell == CurrentOrigin%BoxNumber) then

                if (CurrentOrigin%ParticleBoxVolume > 0.0) then

                    VOLCEL = EulerModel%GridCellArea(i, j) *                         &
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

                        call Convert_CellIJ_XY(EulerModel, NewParticle%Position)

                        !Z position
                        NewParticle%Position%Z = DepthPartic

                        call Convert_Z_CellK  (CurrentOrigin, EulerModel, NewParticle%Position)
                        call Convert_CellK_K  (               NewParticle%Position)

                        NewParticle%Geometry%VolVar = 0.0

                        
                        if (CurrentOrigin%ParticleBoxVolume > 0.0) then
                            NewParticle%Geometry%Volume     = CurrentOrigin%ParticleBoxVolume
                        else
                            NewParticle%Geometry%Volume =                                                   &
                                (EulerModel%Grid%ParticXX(i, j+1) - EulerModel%Grid%ParticXX(i, j)) * &
                                (EulerModel%Grid%ParticYY(i+1, j) - EulerModel%Grid%ParticYY(i, j)) * &
                                 BoxThickness
                        endif

                        !Stores initial Volume
                        NewParticle%Geometry%InitialVolume = NewParticle%Geometry%Volume
                        
                        if (CurrentOrigin%State%Oil)                                    &
                            NewParticle%Geometry%VolumeOil     = NewParticle%Geometry%Volume

                        NewParticle%Position%ModelID = emBox
                       
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

OP:         if ((EulerModel%OpenPoints3D(i, j, k) == OpenPoint) .and. &
                (Boxes(i, j, k) == CurrentOrigin%BoxNumber)) then
 
                if (CurrentOrigin%ParticleBoxVolume > 0.0) then                                 

                    VOLCEL = EulerModel%VolumeZ(I,J,K) + FACTOR    
                                                       
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

                        call Convert_CellIJ_XY(EulerModel, NewParticle%Position)

                        !Random Z position
                        call RANDOM_NUMBER(auxZ)
                        NewParticle%Position%CellK = ZI + DPPZ * auxZ 
                        call Convert_CellK_Z  (EulerModel,NewParticle%Position)
                        call Convert_CellK_K  (NewParticle%Position)

                        NewParticle%Geometry%VolVar     = 0.0
                        
                        if (CurrentOrigin%ParticleBoxVolume > 0.0) then
                            NewParticle%Geometry%Volume        = CurrentOrigin%ParticleBoxVolume
                        else
                            NewParticle%Geometry%Volume        =                             &
                                   EulerModel%VolumeZ(NewParticle%Position%I, &
                                                                     NewParticle%Position%J, &
                                                                     NewParticle%Position%K)
                        endif

                        !Stores initial Volume
                        NewParticle%Geometry%InitialVolume = NewParticle%Geometry%Volume

                        if (CurrentOrigin%State%Oil)                                    &
                            NewParticle%Geometry%VolumeOil     = NewParticle%Geometry%Volume

                        NewParticle%Position%ModelID = emBox

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
        call UngetBoxDif(EulerModel%ObjBoxDif, Boxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'EmissionBox - ModuleLagrangianGlobal - ERR30'

        nullify(EulerModel)


    end subroutine EmissionBox

    !--------------------------------------------------------------------------

    subroutine EmissionPoint (CurrentOrigin)
        
        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Local-----------------------------------------------------------------
        integer                                     :: nP
        type (T_Partic), pointer                    :: NewParticle
        !logical                                     :: CanEmit
        real                                        :: ParticleVolume

        !Begin-----------------------------------------------------------------

!        CanEmit = .true. 

        !Actualizes Variable Flow
        if (CurrentOrigin%EmissionTemporal == Continuous_)  then
            call ActualizeOrigin (CurrentOrigin, ParticleVolume)
        endif

!        if (CurrentOrigin%EmissionTemporal == Continuous_ .and. CurrentOrigin%Flow<=0)  then
            !CanEmit = .false.
!            CurrentOrigin%Flow = 1e-15
!        endif

!ic:     if (CanEmit) then

        do nP = 1, CurrentOrigin%NbrParticlesIteration

            !Allocates a new Particle
            call AllocateNewParticle (NewParticle, CurrentOrigin%nProperties,            &
                                      CurrentOrigin%NextParticID)
            
            !Sets Initial Position
            NewParticle%Position = CurrentOrigin%Position


            if (CurrentOrigin%Position%MaintainRelative .and. CurrentOrigin%Position%DepthDefinition == Cells) then

                NewParticle%Position%CellK = NewParticle%Position%Depth

                call Convert_CellK_Z (Me%EulerModel(CurrentOrigin%Position%ModelID), NewParticle%Position)
                call Convert_CellK_K (NewParticle%Position)


            endif
            

            NewParticle%Geometry%VolVar     = 0.0

            !Volume
            select case (CurrentOrigin%EmissionTemporal)

            case (Continuous_)


                if (CurrentOrigin%FlowVariable) then
                    !Volume of the Seed Particle
                    !See subroutine ActualizeOrigin 
                    NewParticle%Geometry%Volume = ParticleVolume
                else

                    NewParticle%Geometry%Volume = CurrentOrigin%Flow     *              &
                                                  CurrentOrigin%DT_EMIT  /              &
                                                  CurrentOrigin%NbrParticlesIteration

                endif

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

!        endif ic


    end subroutine EmissionPoint
    
    !--------------------------------------------------------------------------

    subroutine ActualizeOrigin (CurrentOrigin, ParticleVolume)

        !Arguments-------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        real                                        :: ParticleVolume

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        type (T_Time)                               :: Time1, Time2, StartTime, EndTime
        real                                        :: Value1, Value2, TotalVolume
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (CurrentOrigin%FlowVariable) then
            
            !Gets Flow values that limited the current instant
            call GetTimeSerieValue(CurrentOrigin%TimeSerieInputFlow,                        &
                                   Me%Now,                                                  &
                                   CurrentOrigin%FlowColumn,                                &
                                   Time1, Value1, Time2, Value2, TimeCycle,                 &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ActualizeOrigin - ModuleLagrangianGlobal - ERR10.'

            if (TimeCycle) then
                CurrentOrigin%Flow = Value1
                ParticleVolume = CurrentOrigin%Flow     *                                   &
                                 CurrentOrigin%DT_EMIT  /                                   &
                                 CurrentOrigin%NbrParticlesIteration

            else
                !Interpolates Value for current instant
    !            call InterpolateValueInTime(Me%ExternalVar%Now, Time1,               &
    !                                        Value1, Time2, Value2, CurrentOrigin%Flow)

                StartTime = Me%Now 
                EndTime   = Me%Now + CurrentOrigin%DT_EMIT
                if (EndTime > Me%ExternalVar%EndTime) then
                    CurrentOrigin%NbrParticlesIteration = 0
                else

                    TotalVolume = GetTimeSerieIntegral(CurrentOrigin%TimeSerieInputFlow,        &
                                                       StartTime, EndTime,                      &
                                                       CurrentOrigin%FlowColumn, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizeOrigin - ModuleLagrangianGlobal - ERR20.'

                    if (TotalVolume == 0.) then
                        CurrentOrigin%NbrParticlesIteration = 0
                    else
                        CurrentOrigin%NbrParticlesIteration = int (TotalVolume / CurrentOrigin%PointVolume) + 1
                        ParticleVolume                      = TotalVolume / real(CurrentOrigin%NbrParticlesIteration)
                    endif
                endif
            endif

        endif

        CurrentProperty => CurrentOrigin%FirstProperty
        do while (associated(CurrentProperty))

            if (CurrentProperty%ConcVariable) then
                call ActualizeConcentration (CurrentProperty)
            endif

            CurrentProperty => CurrentProperty%Next
        enddo

        nullify(CurrentProperty)


    end subroutine ActualizeOrigin 

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine ActualizeOriginDefault 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty


        !Begin-----------------------------------------------------------------
        if (associated(Me%OriginDefault)) then
            CurrentProperty => Me%OriginDefault%FirstProperty
            do while (associated(CurrentProperty))

                if (CurrentProperty%T90Variable .and. CurrentProperty%T90Var_Method == FromTimeSerie) then
                    call ActualizeT90           (CurrentProperty)
                endif

                CurrentProperty => CurrentProperty%Next
            enddo

            nullify(CurrentProperty)
        endif

    end subroutine ActualizeOriginDefault  

    !--------------------------------------------------------------------------
    
    subroutine ActualizeConcentration (CurrentProperty)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Gets Flow values that limited the current instant
        call GetTimeSerieValue(CurrentProperty%TimeSerieInput,                          &
                               Me%Now,                                      &
                               CurrentProperty%ConcColumn,                              &
                               Time1, Value1, Time2, Value2, TimeCycle,                 &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeOrigin - ModuleLagrangianGlobal - ERR10.'

        if (TimeCycle) then
            CurrentProperty%Concentration = Value1
        else
            !Interpolates Value for current instant
            call InterpolateValueInTime(Me%Now, Time1,                      &
                                        Value1, Time2, Value2,                          &
                                        CurrentProperty%Concentration)
        endif

    end subroutine ActualizeConcentration        

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine ActualizeT90 (CurrentProperty)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: Time1, Time2
        real                                        :: Value1, Value2
        logical                                     :: TimeCycle
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Gets Flow values that limited the current instant
        call GetTimeSerieValue(CurrentProperty%TimeSerieT90,                            &
                               Me%Now,                                      &
                               CurrentProperty%T90Column,                               &
                               Time1, Value1, Time2, Value2, TimeCycle,                 &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeT90 - ModuleLagrangianGlobal - ERR10.'

        if (TimeCycle) then
            CurrentProperty%T90 = Value1
        else
            !Interpolates Value for current instant
            call InterpolateValueInTime(Me%Now, Time1,                      &
                                        Value1, Time2, Value2,                          &
                                        CurrentProperty%T90)
        endif

    end subroutine ActualizeT90        

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
        integer                                     :: STAT_CALL, I, J, OutPutNumber, em
        logical                                     :: OutPutJet
        !Begin-----------------------------------------------------------------

        I = CurrentOrigin%Position%I
        J = CurrentOrigin%Position%J

        em = CurrentOrigin%Position%ModelID

        
        !Gets the temperature and the Salinity from the Eulerian model
        call GetConcentration(Me%EulerModel(em)%ObjWaterProperties,                                    &
                              ConcentrationX    = Temperature3D,                        &
                              PropertyXIDNumber = Temperature_,                         &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR01'

        call GetConcentration(Me%EulerModel(em)%ObjWaterProperties,                                    &
                              ConcentrationX    = Salinity3D,                           &
                              PropertyXIDNumber = Salinity_,                            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR02'

        OutPutNumber = Me%OutPut%NextOutPut

        if (Me%Output%Write_) then 
            if (Me%Now >= Me%OutPut%OutTime(OutPutNumber)) then 
                OutPutJet = .true.
            else
                OutPutJet = .false.
            endif

        else
            OutPutJet = .false.
        endif
            
        call ModifyJet(JetID          = CurrentOrigin%Movement%ObjJet,                  &
                       Salinity       = Salinity3D,                                     &
                       Temperature    = Temperature3D,                                  &
                       VelU           = Me%EulerModel(em)%Velocity_U,                   &
                       VelV           = Me%EulerModel(em)%Velocity_V,                   &
                       VelW           = Me%EulerModel(em)%Velocity_W,                   &
                       SZZ            = Me%EulerModel(em)%SZZ,                          &
                       I              = I,                                              &
                       J              = J,                                              &
                       BottomLayer    = Me%EulerModel(em)%kFloor(I, J),                 &
                       SurfaceLayer   = Me%EulerModel(em)%WorkSize%KUB,                 &
                       OutPutOK       = .true.,                                         &
                       STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR03'

        !UnGets the temperature and the Salinity from the Eulerian model
        call UnGetWaterProperties(Me%EulerModel(em)%ObjWaterProperties,                                &
                                  Temperature3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR04'

        call UnGetWaterProperties(Me%EulerModel(em)%ObjWaterProperties,                                &
                                  Salinity3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR05'


        call ExtractDate(Me%Now, Year, Month, Day, Hour, Minutes, Seconds)


        call GetPlumeDilution(CurrentOrigin%Movement%ObjJet, PlumeDilution, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR07'

        call GetPlumeLocation(CurrentOrigin%Movement%ObjJet, PlumeX, PlumeY, PlumeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR08'

        call GetPlumeDensity(CurrentOrigin%Movement%ObjJet, PlumeDensity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR09'

        call GetPlumeTemperature(CurrentOrigin%Movement%ObjJet, PlumeTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR10'

        call GetPlumeSalinity(CurrentOrigin%Movement%ObjJet, PlumeSalinity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR11'

        call GetPlumeVelocity(CurrentOrigin%Movement%ObjJet, PlumeU, PlumeV, PlumeW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR12'

        call GetPlumeMixingHorLength(CurrentOrigin%Movement%ObjJet, PlumeMixingHorLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR13'

        call GetPlumeThickness(CurrentOrigin%Movement%ObjJet, PlumeThickness, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeJetProperties - ModuleLagrangianGlobal - ERR14'


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
        integer                                     :: iP, em
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        em = NewParticle%Position%ModelID


        call GetPlumeVelocity(JetID     = CurrentOrigin%Movement%ObjJet,                 &
                              PlumeVelU = CurrentOrigin%Movement%InitialVelocityU,       &
                              PlumeVelV = CurrentOrigin%Movement%InitialVelocityV,       &
                              PlumeVelW = CurrentOrigin%Movement%InitialVelocityW,       &
                              STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangianGlobal - ERR01'

        NewParticle%U = CurrentOrigin%Movement%InitialVelocityU
        NewParticle%V = CurrentOrigin%Movement%InitialVelocityV
        NewParticle%W = CurrentOrigin%Movement%InitialVelocityW

        call GetPlumeDilution(JetID         = CurrentOrigin%Movement%ObjJet,             &
                              PlumeDilution = CurrentOrigin%Movement%PlumeDilution,      &
                              STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangianGlobal - ERR02'

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

                if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangianGlobal - ERR03'
               


            else if (CurrentProperty%ID == Temperature_) then
                
                call GetPlumeTemperature (JetID            = CurrentOrigin%Movement%ObjJet,&
                                          PlumeTemperature = NewParticle%Concentration(iP),&
                                          STAT             = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangianGlobal - ERR04'

            else

                !Properties dilution
                call GetAmbientConcentration (CurrentProperty,                          &
                                              em,                                       &
                                              NewParticle%Position,                     &
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
                              PlumeX = NewParticle%Position%CartX,                       &
                              PlumeY = NewParticle%Position%CartY,                       &
                              PlumeZ = NewParticle%Position%Z,                           &
                              STAT   = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'GiveJetPropertiesToParticle - ModuleLagrangianGlobal - ERR05'

        call Convert_Z_CellK   (CurrentOrigin, Me%EulerModel(em), NewParticle%Position)
        call Convert_CellK_K   (NewParticle%Position)

        call Convert_XY_CellIJ (Me%EulerModel(em), NewParticle%Position, Referential = Cartesian_)


        if (Me%EulerModel(em)%OpenPoints3D(NewParticle%Position%i,                      &
                                                   NewParticle%Position%j,              &
                                                   NewParticle%Position%k) /= OpenPoint) then
            
            !Sets Horizontal Position equal to origin
            NewParticle%Position%X = CurrentOrigin%Position%X
            NewParticle%Position%Y = CurrentOrigin%Position%Y

            call Convert_XY_CellIJ (Me%EulerModel(em), NewParticle%Position, Referential = AlongGrid_)

            write(*,*) 'The jet of Origin = ', trim(CurrentOrigin%Name), ' intersected land' 

        endif

    end subroutine GiveJetPropertiesToParticle

    !--------------------------------------------------------------------------

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
        integer                                     :: i, j, k, KUB, em
        logical                                     :: Emited
        type (T_Position)                           :: NewPosition


        !Cell of the Accident
        JJ = CurrentOrigin%position%J
        II = CurrentOrigin%position%I         

        em = CurrentOrigin%Position%ModelID

        
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
            if (STAT_CALL /= SUCCESS_) stop 'EmissionAccident - ModuleLagrangianGlobal - ERR01'
            
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

                !call LocateEulerModel(NewParticle)

                em = NewParticle%Position%ModelID

                !Converts Horizontal Position
                call Convert_XY_CellIJ (Me%EulerModel(em), NewPosition, Referential = AlongGrid_)

                i   = NewPosition%I
                j   = NewPosition%J
                KUB = Me%EulerModel(em)%WorkSize%KUB
                
                !Stores the Particle
                if (Me%EulerModel(em)%OpenPoints3D(i, j, KUB) == OpenPoint) then


                    NewPosition%ModelID = Locate_ModelDomain(NewPosition%CoordX, NewPosition%CoordY) 


                    !If the vertical position of the Particle isnt a waterpoint, put it 
                    !close to the suface
                    call Convert_Z_CellK (CurrentOrigin, Me%EulerModel(em), NewPosition)
                    call Convert_CellK_K (               NewPosition)


                    k = NewPosition%K
                    if (Me%EulerModel(em)%OpenPoints3D(i, j, k) == OpenPoint) then
                        NewParticle%Position = NewPosition
                    else
                        NewPosition%Z = Me%EulerModel(em)%SZZ(i, j, KUB)
                        call Convert_Z_CellK (CurrentOrigin, Me%EulerModel(em), NewPosition)
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
            CurrentOrigin%NextEmission = Me%Now +                 &
                                         CurrentOrigin%DT_Emit
        else
            CurrentOrigin%NextEmission = Me%Now
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
        integer                                     :: STAT_CALL, em


        !Extracts the Oil Section from the Origin
        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                       &
                                   oil_begin, oil_end,                                  &
                                   OilSectionFound, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructParticOil - ModuleLagrangianGlobal - ERR10'


i1:     if (OilSectionFound .and. .not. NewOrigin%Old) then

            em = NewOrigin%Position%ModelID
            

            !Starts Oil
            call StartOil(OilID             = NewOrigin%ObjOil,                     &     
                          TimeID            = Me%EulerModel(em)%ObjTime,            &     
                          EnterDataID       = Me%ObjEnterData,                      &
                          HorizontalGridID  = Me%EulerModel(em)%ObjHorizontalGrid,  &     
                          GeometryID        = Me%EulerModel(em)%ObjGeometry,        &     
                          MapID             = Me%EulerModel(em)%ObjMap,             &     
                          DT                = Me%DT_PARTIC,                         &     
                          ContCalc          = NewOrigin%Old,                        &
                          ExtractType       = FromBlockInBlock,                     &
                          STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticOil - ModuleLagrangianGlobal - ERR10'

        elseif (.not. OilSectionFound) then i1

            write (*,*)'Oil Section not defined for Origin : ', trim(NewOrigin%Name)
            stop       'ConstructParticOil - ModuleLagrangianGlobal - ERR40'

        endif i1
  
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
        integer                                     :: i, j, k, Box, em
        integer                                     :: STAT_CALL

        !-----------------------------------------------------------------------
        !Get the boxes only
        em = 1
        call GetNumberOfBoxes(Me%EulerModel(em)%ObjBeachingProbBox, NumberOfBoxes3D = NumberOfBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangianGlobal - ERR03'

        allocate (BoxesBeachingProbability(NumberOfBoxes))

        call GetData(BoxesBeachingProbability,                                  &
                     Me%ObjEnterData,                                           &
                     flag,                                                      &
                     SearchType = FromFile,                                     &
                     keyword    = 'BOXES_BEACHING_PROB',                        &
                     ClientModule = 'ModuleLagrangianGlobal',                         &
                     STAT       = STAT_CALL)
        if       (STAT_CALL .EQ. SIZE_ERR_)  then
            write(*,*) 
            write(*,*) 'Error calling GetData.  '
            write(*,*) 'Number of box values is incorrect:'
            write(*,*) '    NumberOfBoxes =', NumberOfBoxes
            write(*,*) '    BoxesData   =', flag
            stop       'Subroutine VerifyBeachingProbabilities; Module ModuleLagrangianGlobal. ERR04.'

        else if ((STAT_CALL .NE. SIZE_ERR_) .AND.                               &
                 (STAT_CALL .NE. SUCCESS_)) then                                                                        
            stop 'Subroutine VerifyBeachingProbabilities; Module ModuleLagrangianGlobal. ERR05.'
        end if               

        if (flag==0) then
            write(*,*) 
            write(*,*) 'Error do not have the box beaching probability.'           
            stop       'Subroutine VerifyBeachingProbabilities; Module ModuleLagrangianGlobal. ERR06.'
        end if



d1:     do em = 1, Me%EulerModelNumber
        
            allocate (Me%EulerModel(em)%BeachingProbability(Me%EulerModel(em)%Size%ILB: &
                                                            Me%EulerModel(em)%Size%IUB, &
                                                            Me%EulerModel(em)%Size%JLB: &
                                                            Me%EulerModel(em)%Size%JUB, &
                                                            Me%EulerModel(em)%Size%KLB: &
                                                            Me%EulerModel(em)%Size%KUB),&
                                                            STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangianGlobal - ERR01'
        
            if ( .NOT. Me%State%HaveBeachingProbBox) then

                Me%EulerModel(em)%BeachingProbability = Me%DefaultBeachingProbability

            else if (Me%State%HaveBeachingProbBox) then

                !Shorten
                ILB    = Me%EulerModel(em)%WorkSize%ILB
                IUB    = Me%EulerModel(em)%WorkSize%IUB
                JLB    = Me%EulerModel(em)%WorkSize%JLB
                JUB    = Me%EulerModel(em)%WorkSize%JUB
                KLB    = Me%EulerModel(em)%WorkSize%KLB
                KUB    = Me%EulerModel(em)%WorkSize%KUB
        

                call GetBoxes(Me%EulerModel(em)%ObjBeachingProbBox, BeachingProbBoxes, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangianGlobal - ERR02'

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    Box = BeachingProbBoxes(i, j, k)
            
                    if (Box .GE. 0 ) then
            
                        Me%EulerModel(em)%BeachingProbability (i,j,k) = BoxesBeachingProbability(Box)
                     
                    else if (Box .LT. 0 ) then
            
                        Me%EulerModel(em)%BeachingProbability (i,j,k) = Me%DefaultBeachingProbability
            
                    endif
                enddo
                enddo
                enddo

                deallocate (BoxesBeachingProbability)

                !Unget The Boxes
                call UngetBoxDif(Me%EulerModel(em)%ObjBeachingProbBox, BeachingProbBoxes, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyBeachingProbabilities - ModuleLagrangianGlobal - ERR07'

            end if

        enddo d1

    end subroutine VerifyBeachingProbabilities

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Output

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel
        real, pointer, dimension(:, :)              :: Bathymetry
        integer, pointer, dimension(:, :, :)        :: WaterPoints3D
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: HDF5_CREATE, em, ic

 
        !Begin-----------------------------------------------------------------

        allocate(Me%ObjHDF5(Me%EulerModelNumber))

        Me%ObjHDF5(:) = 0

em1:    do em =1, Me%EulerModelNumber

            EulerModel => Me%EulerModel(em)

            !Bounds
            ILB = EulerModel%WorkSize%ILB
            IUB = EulerModel%WorkSize%IUB

            JLB = EulerModel%WorkSize%JLB
            JUB = EulerModel%WorkSize%JUB

            KLB = EulerModel%WorkSize%KLB
            KUB = EulerModel%WorkSize%KUB

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

            ic = len(trim(Me%Files%TransientHDF))

            !Opens HDF File
            call ConstructHDF5      (Me%ObjHDF5(em),                                    &
                                     trim(Me%Files%TransientHDF(1:ic-4))//"_"//         &
                                     trim(Me%EulerModel(em)%Name)//".hdf5",             &
                                     HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR05'


            if (.not. Me%RunOnline) then
                !Write the Horizontal Grid
                call WriteHorizontalGrid(EulerModel%ObjHorizontalGrid, Me%ObjHDF5(em),  &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR06'

                !Gets a pointer to Bathymetry
                call GetGridData        (EulerModel%ObjGridData, Bathymetry, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR01'

                !Gets WaterPoints3D
                call GetWaterPoints3D   (EulerModel%ObjMap, WaterPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR03'


                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5(em), ILB, IUB, JLB, JUB, KLB, KUB,       &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR17'


                !Writes the Grid
                call HDF5WriteData   (Me%ObjHDF5(em), "/Grid", "Bathymetry", "m",         &
                                      Array2D = Bathymetry,                                      &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR18'

                call HDF5WriteData   (Me%ObjHDF5(em), "/Grid", "WaterPoints3D", "-",      &
                                      Array3D = WaterPoints3D,                                   &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR19'

                !Flushes All pending HDF5 commands
                call HDF5FlushMemory (Me%ObjHDF5(em), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR23'

                !Ungets the Bathymetry
                call UngetGridData (EulerModel%ObjGridData, Bathymetry, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR24'

                !Ungets the Waterpoints3D
                call UnGetMap        (EulerModel%ObjMap, WaterPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleLagrangianGlobal - ERR25'

            endif

        enddo em1

        nullify(EulerModel)

    end subroutine ConstructHDF5Output
  
    !--------------------------------------------------------------------------

    subroutine ConstructParticStatistic

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character (StringLength)                        :: GroupName
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                         :: WS_KLB, WS_KUB
        integer                                         :: nProp, em, nP, ig, STAT_CALL
        !Begin------------------------------------------------------------------
                
d1:     do em = 1, Me%EulerModelNumber 

            ILB    = Me%EulerModel(em)%Size%ILB
            IUB    = Me%EulerModel(em)%Size%IUB
            JLB    = Me%EulerModel(em)%Size%JLB
            JUB    = Me%EulerModel(em)%Size%JUB
            KLB    = Me%EulerModel(em)%Size%KLB
            KUB    = Me%EulerModel(em)%Size%KUB
            WS_ILB = Me%EulerModel(em)%WorkSize%ILB
            WS_IUB = Me%EulerModel(em)%WorkSize%IUB
            WS_JLB = Me%EulerModel(em)%WorkSize%JLB
            WS_JUB = Me%EulerModel(em)%WorkSize%JUB
            WS_KLB = Me%EulerModel(em)%WorkSize%KLB
            WS_KUB = Me%EulerModel(em)%WorkSize%KUB

            nProp           =  Me%Statistic%PropNumber

        !Fills Grid concentration
        !if (Me%Now > Me%ExternalVar%LastConcCompute)  call FillGridConcentration  

            allocate(Me%EulerModel(em)%PropStatistic(1:nProp))

            do nP = 1, nProp

            allocate(Me%EulerModel(em)%PropStatistic(nP)%Statistic1_ID(1:Me%NGroups))

            Me%EulerModel(em)%PropStatistic(nP)%Statistic1_ID(1:Me%NGroups) = 0

            if (Me%OutPut%ConcMaxTracer) then

                allocate(Me%EulerModel(em)%PropStatistic(nP)%Statistic2_ID(1:Me%NGroups))

                Me%EulerModel(em)%PropStatistic(nP)%Statistic2_ID(1:Me%NGroups) = 0

            endif

            do ig = 1, Me%NGroups

                GroupName = trim(Me%EulerModel(em)%Name)//"_"

                call ConstructStatistic (Me%EulerModel(em)%PropStatistic(nP)%Statistic1_ID(ig),     &
                                         ObjTime          = Me%EulerModel(em)%ObjTime,              &
                                         ObjHDF5          = Me%ObjHDF5(em),                             &
                                         Size             = Me%EulerModel(em)%Size,                 &
                                         WorkSize         = Me%EulerModel(em)%WorkSize,             &
                                         DataFile         = Me%Statistic%OptionsStat(nP)%File,      &
                                         Name             = Me%Statistic%OptionsStat(nP)%ID%Name,   &
                                         GroupName        = trim(GroupName),                        &
                                         STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ConstructParticStatistic - ModuleLagrangianGlobal - ERR02'

                if (Me%OutPut%ConcMaxTracer) then

                    GroupName = trim(GroupName)//"MaxTracer_"

                    call ConstructStatistic (Me%EulerModel(em)%PropStatistic(nP)%Statistic2_ID(ig), &
                                             ObjTime      = Me%EulerModel(em)%ObjTime,              &
                                             ObjHDF5      = Me%ObjHDF5(em),                             &
                                             Size         = Me%EulerModel(em)%Size,                 &
                                             WorkSize     = Me%EulerModel(em)%WorkSize,             &
                                             DataFile     = Me%Statistic%OptionsStat(nP)%File,      &
                                             Name         = Me%Statistic%OptionsStat(nP)%ID%Name,   &
                                             GroupName    = GroupName,                              &
                                             STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'ConstructParticStatistic - ModuleLagrangianGlobal - ERR03'

                endif

            enddo
            enddo

            do nP = 1, nProp
            do ig = 1, Me%NGroups

                if (Me%Statistic%OptionsStat(nP)%Lag) then
                    
                    call GetStatisticClassesNumber(Me%EulerModel(em)%PropStatistic(nP)%Statistic1_ID(ig),&
                                                   Me%Statistic%OptionsStat(nP)%nClassesLag, STAT= STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'ConstructParticStatistic - ModuleLagrangianGlobal - ERR04'

                    allocate(Me%EulerModel(em)%PropStatistic(nP)%FrequencyLag           &
                                   (ILB:IUB,JLB:JUB,KLB:KUB,                            &
                                    1:Me%Statistic%OptionsStat(nP)%nClassesLag,         &
                                    1:Me%NGroups))

                    Me%EulerModel(em)%PropStatistic(nP)%FrequencyLag(:,:,:,:,:) = 0.
                
                endif
            enddo
            enddo
        enddo d1


    end subroutine ConstructParticStatistic

    !--------------------------------------------------------------------------

    subroutine ConstructParticLightExtinction
        
        !Local-----------------------------------------------------------------
        integer                                     :: ig, STAT_CALL  
        integer                                     :: LightExtinctionID, em
          
        !Begin-----------------------------------------------------------------

d1:     do em = 1, Me%EulerModelNumber

            allocate (Me%EulerModel(em)%Light(Me%nGroups))

            do ig = 1, Me%nGroups

                Me%EulerModel(em)%Light(ig)%ObjLightExtinction = 0

                allocate(Me%EulerModel(em)%Light(ig)%TopRadiationCells(Me%EulerModel(em)%Size%ILB: &
                                                                       Me%EulerModel(em)%Size%IUB, &
                                                                       Me%EulerModel(em)%Size%JLB: &
                                                                       Me%EulerModel(em)%Size%JUB, &
                                                                       Me%EulerModel(em)%Size%KLB: &
                                                                       Me%EulerModel(em)%Size%KUB))

            enddo

            call GetWaterPropertiesSubModulesID(WaterPropertiesID      = Me%EulerModel(em)%ObjWaterProperties, &
                                                LightExtinctionID      = LightExtinctionID,                    &
                                                STAT                   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructParticLightExtinction - Lagrangian - ERR00'


            if(LightExtinctionID == 0)then

                Me%EulerModel(em)%Light(:)%Compute = ON
          
                do ig = 1, Me%nGroups

                    call ConstructLightExtinction (Me%EulerModel(em)%Light(ig)%ObjLightExtinction, &
                                                   Me%EulerModel(em)%ObjTime, Me%ObjEnterData,         &
                                                   Me%EulerModel(em)%Size, Me%EulerModel(em)%WorkSize, &
                                                   STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructParticLightExtinction - Lagrangian - ERR01'
                enddo

            else

                Me%EulerModel(em)%Light(:)%Compute = OFF

                do ig = 1, Me%nGroups
                    Me%EulerModel(em)%Light(ig)%ObjLightExtinction = AssociateInstance(mLIGHTEXTINCTION_,LightExtinctionID)
                enddo

            end if

        enddo d1




    end subroutine ConstructParticLightExtinction

    !--------------------------------------------------------------------------

    subroutine ConstructLog

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: ig


        write(*, *)"----------------------- LAGRANGIAN -----------------------"
        write(*, *)
        write(*, *)"Number of Origins : ", Me%nOrigins
        write(*, *)"Number of Groups  : ", Me%nGroups
        write(*, *)

        do ig = 1, Me%nGroups
            write(*, *)"GroupID           : ", Me%GroupIDs(ig)
            write(*, *)"---Number of Orig.: ", Me%nOriginsGroup(ig)
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
        if (STAT_CALL /= SUCCESS_) stop 'AllocateNewOrigin - ModuleLagrangianGlobal - ERR01'
        nullify  (NewOrigin%Next)
        nullify  (NewOrigin%FirstPartic)
        nullify  (NewOrigin%FirstProperty)

        NewOrigin%nProperties  = 0
        NewOrigin%nPropT90     = 0
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
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangianGlobal - ERR10'
                endif

                !Kill WaterQuality
                if (OriginToDelete%State%WQM) then
                    call KillWaterQuality (OriginToDelete%WaterQualityID, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangianGlobal - ERR20'
                endif

                !Kill TimeSerie
                if (OriginToDelete%MovingOrigin) then
                    call KillTimeSerie(OriginToDelete%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangianGlobal - ERR30'
                endif

                !Kills ObjJet
                if (OriginToDelete%Movement%ObjJet /= 0) then
                   

                    write(OriginToDelete%Movement%JetUnit,*) '<EndTimeSerie>'

                    call UnitsManager(OriginToDelete%Movement%JetUnit, CLOSE_FILE, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangianGlobal - ERR40'

                    call KillJet (JetID = OriginToDelete%Movement%ObjJet,  STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangianGlobal - ERR50'
                endif

                
                if (OriginToDelete%FlowVariable) then
                    call KillTimeSerie(OriginToDelete%TimeSerieInputFlow, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangianGlobal - ERR60'
                endif



                !Kill PropertyList
                CurrentProperty => OriginToDelete%FirstProperty
                do while (associated(CurrentProperty))
                    if (CurrentProperty%ConcVariable) then
                        call KillTimeSerie(CurrentProperty%TimeSerieInput, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'DeleteOrigin - ModuleLagrangianGlobal - ERR70'
                    endif
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
        if (STAT_CALL /= SUCCESS_) stop 'AllocateNewProperty - ModuleLagrangianGlobal - ERR01'
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

        if(NewProperty%T90ON) then
            Origin%nPropT90 = Origin%nPropT90 + 1
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
        integer                                     :: iP, emp
        real                                        :: Concentration

        !Initializes Mass and Concentration of the new particle
        if (InitConcentration) then

            iP = 1
            CurrentProperty => Origin%FirstProperty
            do while (associated(CurrentProperty))
                

                if (CurrentProperty%EqualToAmbient) then
                    emp = NewParticle%Position%ModelID

                    call GetAmbientConcentration (CurrentProperty,                      &
                                                  emp,                                  & 
                                                  NewParticle%Position,                 &
                                                  Concentration)

                else

                    Concentration = CurrentProperty%Concentration
                endif


                NewParticle%Concentration(iP) = Concentration
                NewParticle%Mass         (iP) = Concentration * NewParticle%Geometry%Volume
        
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
            stop 'DeleteParticleDeleted - ModuleLagrangianGlobal - ERR10'
        endif

    end subroutine DeleteParticle
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !----------------------------------------------------------------------------
    
    subroutine SetLagrangianAtmPressureGlobal(LagrangianID, ModelName, AtmPressure, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: LagrangianID
        character(len=*)                            :: ModelName
        real, pointer, dimension(:,:)               :: AtmPressure
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_, ModelID          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            ModelID = ReturnModelIndex (ModelName)

            Me%EulerModel(ModelID)%AtmPressure => AtmPressure
                        
            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetLagrangianAtmPressureGlobal

    !----------------------------------------------------------------------------
    
    subroutine SetLagSolarRadiationGlobal(LagrangianID, ModelName, SolarRadiation, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: LagrangianID
        character(len=*)                            :: ModelName
        real, pointer, dimension(:,:)               :: SolarRadiation
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_, ModelID          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            ModelID = ReturnModelIndex (ModelName)

            Me%EulerModel(ModelID)%SurfaceRadiation => SolarRadiation
                        
            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetLagSolarRadiationGlobal

    !----------------------------------------------------------------------------
    
    subroutine SetLagrangianWindGlobal(LagrangianID, ModelName, WindX, WindY, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: LagrangianID
        character(len=*)                            :: ModelName
        real, pointer, dimension(:,:)               :: WindX, WindY
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_, ModelID          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            if(Me%State%Wind)then

                ModelID = ReturnModelIndex (ModelName)

                Me%EulerModel(ModelID)%WindX => WindX
                Me%EulerModel(ModelID)%WindY => WindY
            end if


            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetLagrangianWindGlobal


    !----------------------------------------------------------------------

    subroutine SetLagrangianShearGlobal(LagrangianID, ModelName, ShearStress, ShearVelocity, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: LagrangianID
        character(len=*)                :: ModelName
        real, dimension(:,:), pointer   :: ShearStress, ShearVelocity
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_, ModelID             
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            ModelID = ReturnModelIndex (ModelName)

            if(Me%State%Deposition) Me%EulerModel(ModelID)%BottomStress  => ShearStress            
            if(Me%State%ShearVel  ) Me%EulerModel(ModelID)%ShearVelocity => ShearVelocity 
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetLagrangianShearGlobal


    !----------------------------------------------------------------------

    subroutine GetLagrangianAirOptionsGlobal(LagrangianID, Oil, Wind, WaterQuality, T90Variable, STAT)

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


    end subroutine GetLagrangianAirOptionsGlobal  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine ModifyLagrangianGlobal(LagrangianID, STAT)

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

                if (Me%OverLay) then
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

                if (Me%State%WQM .or. Me%State%T90Compute) then
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

                !Check the concentration limits
                call CheckConcLimits        ()

                !Eliminates Particles
                call PurgeParticles         ()

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

                Me%FirstIteration = .false. 
            enddo
            !endif

            call ReadUnLockExternalVar  ()

            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLagrangianGlobal

    !--------------------------------------------------------------------------
    subroutine OverLayProcess ()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: WorkSizeILB, WorkSizeIUB
        integer                                     :: WorkSizeJLB, WorkSizeJUB
        integer                                     :: WorkSizeKLB, WorkSizeKUB
        integer                                     :: i, j, k
        integer                                     :: STAT_CALL, em

        do em = 1, Me%EulerModelNumber

            !Shorten
            WorkSizeILB = Me%EulerModel(em)%WorkSize%ILB
            WorkSizeIUB = Me%EulerModel(em)%WorkSize%IUB
            WorkSizeJLB = Me%EulerModel(em)%WorkSize%JLB
            WorkSizeJUB = Me%EulerModel(em)%WorkSize%JUB
            WorkSizeKLB = Me%EulerModel(em)%WorkSize%KLB
            WorkSizeKUB = Me%EulerModel(em)%WorkSize%KUB


            call GetAssimilationField(Me%EulerModel(em)%ObjAssimilation,                 &
                                      ID          = VelocityU_,                      &
                                      Field3D     = Me%EulerModel(em)%OverLayU,  &
                                      STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangianGlobal - ERR01' 

            call GetAssimilationField(Me%EulerModel(em)%ObjAssimilation,                 &
                                      ID          = VelocityV_,                      &
                                      Field3D     = Me%EulerModel(em)%OverLayV,  &
                                      STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangianGlobal - ERR02' 

            Me%EulerModel(em)%OverLay%VelUFinal = 0.
            Me%EulerModel(em)%OverLay%VelVFinal = 0.

            do k = WorkSizeKLB, WorkSizeKUB
            do j = WorkSizeJLB, WorkSizeJUB
            do i = WorkSizeILB, WorkSizeIUB

                Me%EulerModel(em)%OverLay%VelUFinal(i, j, k) = (Me%EulerModel(em)%OverLayU   (i, j  , k)     +   &
                                                            Me%EulerModel(em)%OverLayU   (i, j+1, k))/2. +   &
                                                            Me%EulerModel(em)%Velocity_U (i, j  , k)

                Me%EulerModel(em)%OverLay%VelVFinal(i, j, k) = (Me%EulerModel(em)%OverLayV   (i  , j, k)     +   &
                                                            Me%EulerModel(em)%OverLayV   (i+1, j, k))/2. +   &
                                                            Me%EulerModel(em)%Velocity_V (i  , j, k)
            enddo
            enddo
            enddo

            call UngetAssimilation (Me%EulerModel(em)%ObjAssimilation, Me%EulerModel(em)%OverLayU, &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangianGlobal - ERR03' 

            call UngetAssimilation (Me%EulerModel(em)%ObjAssimilation, Me%EulerModel(em)%OverLayV, &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OverLayProcess - ModuleLagrangianGlobal - ERR04' 

        enddo


    end subroutine OverLayProcess

    !--------------------------------------------------------------------------

    subroutine ParticleEmission ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin

        !Begin-----------------------------------------------------------------

        if (Me%State%KillPartInsideBox) then
            call KillParticInsideBoxes ()

            call PurgeParticles        ()
        endif

        call ActualizeOriginDefault


        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

iON:        if (CurrentOrigin%EmissionON) then

                !Constructs the Temporal emission 
                select case (CurrentOrigin%EmissionTemporal)

                case (Continuous_)


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


                case (Instantaneous_)

                    if (.not. CurrentOrigin%Old) then

                        if (Me%Now >= CurrentOrigin%InstantEmission .and. CurrentOrigin%OnlyOnceEmit) then

                            !Realizes only one emission
                            select case (CurrentOrigin%EmissionSpatial)

                            case (Box_)

                                call EmissionBox           (CurrentOrigin)

                            case (Point_)

                                call EmissionPoint         (CurrentOrigin)

                            case (Accident_)

                                call EmissionAccident      (CurrentOrigin)

                            end select

                            CurrentOrigin%OnlyOnceEmit = .false.

                        endif

                    endif

                end select

            endif iON

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr


    end subroutine ParticleEmission

   !--------------------------------------------------------------------------

    subroutine FillGridThickness ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, KUB
        integer                                     :: i, j, em, ig

        !Begin-----------------------------------------------------------------
                !Fills Grid concentration
        if (Me%Now > Me%ExternalVar%LastConcCompute) call FillGridConcentration 

d1:     do em = 1, Me%EulerModelNumber

            ILB = Me%EulerModel(em)%WorkSize%ILB
            JLB = Me%EulerModel(em)%WorkSize%JLB
            IUB = Me%EulerModel(em)%WorkSize%IUB
            JUB = Me%EulerModel(em)%WorkSize%JUB
            KUB = Me%EulerModel(em)%WorkSize%KUB

d2:         do ig = 1, Me%NGroups 
d3:         do j  = JLB, JUB
d4:         do i  = ILB, IUB
                             
i1:             if (Me%EulerModel(em)%WaterPoints3D(i, j, KUB) == WaterPoint) then 

                !Calculates the GridThickness
                    Me%EulerModel(em)%OilSpreading(ig)%GridThickness(i, j) =            &
                        Me%EulerModel(em)%Lag2Euler%GridVolume(i, j, KUB, ig) /         &
                        Me%EulerModel(em)%GridCellArea(i, j)

                endif i1

            enddo d4
            enddo d3
            enddo d2

        enddo d1

    end subroutine FillGridThickness

    !--------------------------------------------------------------------------

    subroutine OilGridConcentration (CurrentOrigin, WaveHeight, WaterDensity)

        !Arguments-------------------------------------------------------------
        type(T_Origin), pointer                     :: CurrentOrigin
        real                                        :: WaveHeight
        real                                        :: WaterDensity

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, KUB
        integer                                     :: i, j, em, ig

        !Begin-----------------------------------------------------------------
                !Fills Grid concentration
i1:     if (Me%State%Oil) then

            em = CurrentOrigin%Position%ModelID
            ig = CurrentOrigin%GroupID

            ILB = Me%EulerModel(em)%WorkSize%ILB
            JLB = Me%EulerModel(em)%WorkSize%JLB
            IUB = Me%EulerModel(em)%WorkSize%IUB
            JUB = Me%EulerModel(em)%WorkSize%JUB
            KUB = Me%EulerModel(em)%WorkSize%KUB

d3:         do j  = JLB, JUB
d4:         do i  = ILB, IUB
                         
i2:             if (Me%EulerModel(em)%WaterPoints3D(i, j, KUB) == WaterPoint) then 

                !Calculates the Concentration
                ! CurrentOrigin%OilGridConcentration = (Specific_MSurface + Specific_MDispersed) / (MixingDepth * WaterDensity)
                    Me%EulerModel(em)%OilSpreading(ig)%OilGridConcentration(i, j) =                     &
                                  1.e6 * ((Me%EulerModel(em)%OilSpreading(ig)%GridThickness(i, j) * &
                                  Me%ExternalVar%OilDensity) + & 
                                  (Me%ExternalVar%MDispersed / Me%EulerModel(em)%GridCellArea(i, j))) / & 
                                  (1.5 * WaveHeight * WaterDensity)                      

                endif i2

            enddo d4
            enddo d3

        endif i1


    end subroutine OilGridConcentration

    !--------------------------------------------------------------------------

     subroutine ParticleDensity () 
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin),   pointer                  :: CurrentOrigin
        type (T_Partic),   pointer                  :: CurrentPartic
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: TemperatureID, SalinityID, emp, iProp
        real                                        :: T, S, Depth



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
                        stop      'ParticleDensity - ModuleLagrangianGlobal - ERR10'
                    endif
                    CurrentPartic => CurrentPartic%Next
                enddo
                
                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))
                    T = CurrentPartic%Concentration (TemperatureID)
                    S = CurrentPartic%Concentration (SalinityID   )

                    emp = CurrentPartic%Position%ModelID
                    
                    Depth = - Me%EulerModel(emp)%ZCellCenter(CurrentPartic%Position%I, &
                                                            CurrentPartic%Position%J, &
                                                            CurrentPartic%Position%K)

                    if (CurrentOrigin%Movement%DensityMethod == UNESCOState_) then
                        CurrentPartic%SigmaDensity = SigmaUNESCOPressureCorrection (T, S, Depth,CurrentPartic%SigmaDensity)
           
                    elseif (CurrentOrigin%Movement%DensityMethod == LeendertseState_) then
                        write(*,*)'Invalid Density Method'
                        stop      'ParticleDensity - ModuleLagrangianGlobal - ERR20'

                    else 
                        write(*,*)'Invalid Density Method'
                        stop      'ParticleDensity - ModuleLagrangianGlobal - ERR30'
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
        real,   dimension(:,:,:), pointer           :: MassSedGrid, TauErosionGrid       
        integer                                     :: i, j, ILB, IUB, JLB, JUB, &
                                                       Sediment_ID, em, ig, emp
        logical                                     :: FoundSediment

        !Begin-----------------------------------------------------------------

d1:     do em = 1, Me%EulerModelNumber        

            Me%EulerModel(em)%Lag2Euler%TauErosionGrid(:, :, :) = 0.
            Me%EulerModel(em)%Lag2Euler%MassSedGrid   (:, :, :) = 0.

        enddo d1

d2:     do ig = 1, Me%NGroups 

            CurrentOrigin => Me%FirstOrigin
CurrOr:         do while (associated(CurrentOrigin))

                if (CurrentOrigin%GroupID /= ig) then
                    CurrentOrigin => CurrentOrigin%Next
                    cycle
                endif
    
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

                    i  = CurrentPartic%Position%I
                    j  = CurrentPartic%Position%J

                    emp = CurrentPartic%Position%ModelID

                    if (CurrentPartic%Deposited) then

                        Me%EulerModel(emp)%Lag2Euler%TauErosionGrid (i, j, ig) =     &
                            Me%EulerModel(emp)%Lag2Euler%TauErosionGrid (i, j, ig) + &
                            CurrentPartic%Mass(Sediment_ID) * CurrentOrigin%Deposition%TauErosion

                        Me%EulerModel(emp)%Lag2Euler%MassSedGrid(i, j, ig) =                                  &
                            Me%EulerModel(emp)%Lag2Euler%MassSedGrid(i, j, ig) + CurrentPartic%Mass(Sediment_ID) 

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

                    emp= CurrentPartic%Position%ModelID

                    if (CurrentPartic%Deposited .and. Me%EulerModel(emp)%Lag2Euler%MassSedGrid(i, j, ig) > 1e-12) then
                    
                        CurrentPartic%ErosionRateProbability =                      &
                            min (1., CurrentOrigin%Deposition%ErosionRate *         &
                                     Me%DT_Partic              /                    &
                                     Me%EulerModel(emp)%Lag2Euler%MassSedGrid(i, j, ig))

                    else 
                        CurrentPartic%ErosionRateProbability = 0. 
                    endif

                    CurrentPartic => CurrentPartic%Next
                enddo

                CurrentOrigin => CurrentOrigin%Next

            enddo CurrOr1
            

        enddo d2

        
d3:     do em = 1, Me%EulerModelNumber
            
            ILB = Me%EulerModel(em)%WorkSize%ILB
            JLB = Me%EulerModel(em)%WorkSize%JLB
            IUB = Me%EulerModel(em)%WorkSize%IUB
            JUB = Me%EulerModel(em)%WorkSize%JUB

            MassSedGrid    => Me%EulerModel(em)%Lag2Euler%MassSedGrid
            TauErosionGrid => Me%EulerModel(em)%Lag2Euler%TauErosionGrid 


            do ig =   1, Me%nGroups
            do j  = JLB, JUB
            do i  = ILB, IUB
                
                if (MassSedGrid(i, j, ig) > 1e-12) then

                    TauErosionGrid (i, j, ig) = TauErosionGrid (i, j, ig) / MassSedGrid(i, j, ig)

                else

                    TauErosionGrid (i, j, ig) = 0.

                endif

            enddo
            enddo
            enddo

            nullify(MassSedGrid   )
            nullify(TauErosionGrid)

        enddo d3


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
        integer                                     :: FilterProp_ID, STAT_CALL, em, emp

        !Begin-----------------------------------------------------------------

i1:     if (Me%Now >= Me%NextFiltration) then

            !Begin-----------------------------------------------------------------
d1:         do em = 1, Me%EulerModelNumber

                ILB = Me%EulerModel(em)%WorkSize%ILB
                JLB = Me%EulerModel(em)%WorkSize%JLB
                IUB = Me%EulerModel(em)%WorkSize%IUB
                JUB = Me%EulerModel(em)%WorkSize%JUB
                KLB = Me%EulerModel(em)%WorkSize%KLB
                KUB = Me%EulerModel(em)%WorkSize%KUB

                Me%EulerModel(em)%RelativeMassFilter(:,:,:) = 0. 
                Me%EulerModel(em)%MassFiltered      (:,:,:) = 0. 

                call GetFiltrationRate (Me%EulerModel(em)%ObjWaterProperties, FiltrationRateX,             &
                                        DTEulerian, Fecal_Coliforms_, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop "Lagrangian - ActualizesMassFilterGrid - ERR10"

                if (DTEulerian < Me%DT_Partic) then
                    stop "Lagrangian - ActualizesMassFilterGrid - ERR20"
                endif


                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

                    if (Me%EulerModel(em)%Waterpoints3D(i, j, k) == WaterPoint .and.           &
                        FiltrationRateX(i, j, k) > 0.) then
                        ![ ] = [1/T] * [T] / [L3]
                        Me%EulerModel(em)%RelativeMassFilter(i,j,k) = FiltrationRateX(i, j, k) * &
                                                                  DTEulerian               

                        if (Me%EulerModel(em)%RelativeMassFilter(i,j,k) > 1.)                   &
                            stop "Lagrangian - ActualizesMassFilterGrid - ERR30"
                    endif

                enddo
                enddo
                enddo

                call UnGetWaterProperties (Me%EulerModel(em)%ObjWaterProperties, FiltrationRateX, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop "Lagrangian - ActualizesMassFilterGrid - ERR50"

            enddo d1

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

                    emp = CurrentPartic%Position%ModelID

                    if (Me%EulerModel(emp)%RelativeMassFilter(i,j,k) > 0.) then

                        AuxFilter = Me%EulerModel(emp)%RelativeMassFilter(i,j,k) * CurrentPartic%Mass(FilterProp_ID)

                        Me%EulerModel(emp)%MassFiltered(i, j, k) = Me%EulerModel(emp)%MassFiltered(i, j, k) + AuxFilter

                        CurrentPartic%Mass(FilterProp_ID) = CurrentPartic%Mass(FilterProp_ID) - AuxFilter
                    endif

                    CurrentPartic => CurrentPartic%Next
                enddo

                CurrentOrigin => CurrentOrigin%Next

            enddo CurrOr

d2:         do em = 1, Me%EulerModelNumber

                call SetLagrangianSinksSources(Me%EulerModel(em)%ObjWaterProperties, Fecal_Coliforms_,     &
                                               Me%EulerModel(em)%MassFiltered, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop "Lagrangian - ActualizesMassFilterGrid - ERR60"

            enddo d2

            Me%NextFiltration = Me%NextFiltration + DTEulerian

        endif i1


    end subroutine ActualizesMassFilterGrid

    !--------------------------------------------------------------------------

    subroutine VerifyParticleBeaching()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: i, j, k, emp
        real                                        :: CellI, CellJ, CellK
        real                                        :: BalX, BalY
        real                                        :: Rand1

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))


IfBeaching: if (CurrentOrigin%Beaching) then
                    
                CurrentPartic => CurrentOrigin%FirstPartic
CurrPart:       do while (associated(CurrentPartic))

                    !Grid Cell of the particle
                    i     = CurrentPartic%Position%I
                    j     = CurrentPartic%Position%J
                    k     = CurrentPartic%Position%K

                    !Cell Position
                    CellI = CurrentPartic%Position%CellI
                    CellJ = CurrentPartic%Position%CellJ
                    CellK = CurrentPartic%Position%CellK

                    emp    = CurrentPartic%Position%ModelID

                    !Fraction of the cell
                    BALX  = CellJ - int(CellJ)
                    BALY  = CellI - int(CellI)

IfParticNotBeached: if (.NOT. CurrentPartic%Beached) then

                        call RANDOM_NUMBER(Rand1)

                        if ((Me%EulerModel(emp)%OpenPoints3D(i, j+1, k).NE. OpenPoint)  .AND. &
                              ((1-BALX) *  (Me%EulerModel(emp)%Grid%ParticXX(i, j+1) - &
                               Me%EulerModel(emp)%Grid%ParticXX(i, j)) .LT. Me%BeachingLimit)) then

                            if ((Me%EulerModel(emp)%BeachingProbability(i,j+1,k) .GT. Rand1)   &
                                .OR. (Me%EulerModel(emp)%BeachingProbability(i,j+1,k) .EQ. 1)) &  
                               CurrentPartic%Beached = ON
       
                        else if ((Me%EulerModel(emp)%OpenPoints3D(i, j-1, k).NE. OpenPoint)  .AND.  &
                                   (BALX *  (Me%EulerModel(emp)%Grid%ParticXX(i, j+1) - &
                                    Me%EulerModel(emp)%Grid%ParticXX(i, j)).LT. Me%BeachingLimit)) then


                            if ((Me%EulerModel(emp)%BeachingProbability(i,j-1,k) .GT. Rand1)   &
                                .OR. (Me%EulerModel(emp)%BeachingProbability(i,j-1,k) .EQ. 1)) &  
                               CurrentPartic%Beached = ON
     
                        else if ((Me%EulerModel(emp)%OpenPoints3D(i+1, j, k).NE. OpenPoint)  .AND. &
                                   ((1-BALY) *  (Me%EulerModel(emp)%Grid%ParticYY(i+1, j) - &
                                   Me%EulerModel(emp)%Grid%ParticYY(i, j)).LT. Me%BeachingLimit)) then

                            if ((Me%EulerModel(emp)%BeachingProbability(i+1,j,k) .GT. Rand1) &
                                .OR. (Me%EulerModel(emp)%BeachingProbability(i+1,j,k) .EQ. 1)) &  
                               CurrentPartic%Beached = ON
       
                        else if ((Me%EulerModel(emp)%OpenPoints3D(i-1, j, k).NE. OpenPoint)  .AND. &
                                   (BALY *  (Me%EulerModel(emp)%Grid%ParticYY(i+1, j) -            &
                                   Me%EulerModel(emp)%Grid%ParticYY(i, j)) .LT. Me%BeachingLimit)) then

                            if ((Me%EulerModel(emp)%BeachingProbability(i-1,j,k) .GT. Rand1)   &
                                .OR. (Me%EulerModel(emp)%BeachingProbability(i-1,j,k) .EQ. 1)) &  
                               CurrentPartic%Beached = ON
                
                        end if

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
        integer                                     :: STAT_CALL, emp


        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%MovingOrigin) then

                !Gets Value arround current instant     
                call GetTimeSerieValue(CurrentOrigin%ObjTimeSerie,                       &
                                       Me%Now,                    &
                                       CurrentOrigin%MovingOriginColumnX,                &
                                       Time1, Value1, Time2, Value2, TimeCycle,          &
                                       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MoveOrigin - ModuleLagrangianGlobal - ERR01'

                if (TimeCycle) then
                    NewValueX = Value1
                else
                    !Interpolates Value for current instant
                    call InterpolateValueInTime(Me%Now, Time1,    &
                                                Value1, Time2, Value2, NewValueX)
                endif
                                

                !Gets Value arround current instant    
                call GetTimeSerieValue(CurrentOrigin%ObjTimeSerie,                       &
                                       Me%Now,                    &
                                       CurrentOrigin%MovingOriginColumnY,                &
                                       Time1, Value1, Time2, Value2, TimeCycle,          &
                                       STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MoveOrigin - ModuleLagrangianGlobal - ERR02'

                if (TimeCycle) then
                    NewValueY = Value1
                else
                    !Interpolates Value for current instant
                    call InterpolateValueInTime(Me%Now, Time1,    &
                                                Value1, Time2, Value2, NewValueY)
                endif

                emp = CurrentOrigin%Position%ModelID
                 
                if (trim(CurrentOrigin%MovingOriginUnits) == trim(Char_Cells)) then

                    CurrentOrigin%Position%CellI = NewValueY
                    CurrentOrigin%Position%CellJ = NewValueX

                    !Convert Coordinates
                    call Convert_CellIJ_XY(Me%EulerModel(emp), CurrentOrigin%Position)


                else
                
                    CurrentOrigin%Position%CoordX = NewValueX
                    CurrentOrigin%Position%CoordY = NewValueY

                    !Convert Coordinates
                    call Convert_XY_CellIJ(Me%EulerModel(emp),CurrentOrigin%Position, Referential = GridCoord_)


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
        integer                                     :: i, j, k, em, ig
        integer                                     :: ThicknessGradient, Fay
        integer                                     :: SpreadingMethod

        !Begin-----------------------------------------------------------------



        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))


            !Calls Oil Active Processes to calculate Spreading Velocity
            if (CurrentOrigin%State%Oil .and. CurrentOrigin%nParticle > 0) then

                i            = CurrentOrigin%Position%I
                j            = CurrentOrigin%Position%J
                k            = CurrentOrigin%Position%K
                
                em           = CurrentOrigin%Position%ModelID
                ig           = CurrentOrigin%GroupID

                !Gets the temperature and the Density from the Eulerian model
                call GetConcentration(Me%EulerModel(em)%ObjWaterProperties,             &
                                      ConcentrationX    = Temperature3D,                &
                                      PropertyXIDNumber = Temperature_,                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    write (*,*)'The Lagrangian Module needs the Water Temperature from the eulerian module'
                    stop 'MovePartic - ModuleLagrangianGlobal - ERR10'
                endif


                TemperatureX = Temperature3D (i, j, k)

                call UngetWaterProperties (Me%EulerModel(em)%ObjWaterProperties, Temperature3D,  &
                                           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangianGlobal - ERR70'


                DensityX     = Me%EulerModel(em)%Density (i, j, k)


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
                call OilActiveProcesses(OilID              = CurrentOrigin%ObjOil,                            &
                                        GridThickness      = Me%EulerModel(em)%OilSpreading(ig)%GridThickness,&
                                        WaterTemperature   = TemperatureX,              &
                                        WaterDensity       = DensityX,                  &
                                        VolInic            = VolInic,                   &
                                        DT                 = Me%DT_PARTIC,              &
                                        AreaTotal          = CurrentOrigin%AreaTotal,   &
                                        STAT               = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangianGlobal - ERR20'

                !Gets Oil Spreading Parameters
                call GetOilSpreadingList(ThicknessGradient = ThicknessGradient,          &
                                         Fay               = Fay)
                call GetOilSpreading(CurrentOrigin%ObjOil,                               &    
                                     SpreadingMethod       = SpreadingMethod,            &
                                     STAT                  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangianGlobal - ERR30'

                
                if     (SpreadingMethod == ThicknessGradient) then

                    !call GetOilSpreadingVelocity(CurrentOrigin%ObjOil,                      &
                    !     SpreadingVelocityX = Me%EulerModel(em)%SpreadingVelocityX, &
                    !     SpreadingVelocityY = Me%EulerModel(em)%SpreadingVelocityY, &
                    !     STAT               = STAT_CALL)
                    !if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangianGlobal - ERR03b'

                elseif (SpreadingMethod == Fay              ) then

                    call GetOilSpreadingVelocity(CurrentOrigin%ObjOil,                  &
                         DiffVelocity   = Me%ExternalVar%DiffVelocity,                  &
                         STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangianGlobal - ERR40'
                endif

            endif

            !Moves the particles horizontaly
            call MoveParticHorizontal(CurrentOrigin, ThicknessGradient, Fay, SpreadingMethod)

            !Ungets Spreading Velocity
            if (CurrentOrigin%State%Oil .and. CurrentOrigin%nParticle > 0) then

                if     (SpreadingMethod == ThicknessGradient) then

                    !call UnGetOil (CurrentOrigin%ObjOil,                                 &
                    !               Me%EulerModel(em)%SpreadingVelocityX,         &
                    !               STAT = STAT_CALL)
                    !if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangianGlobal - ERR50'

                    !call UnGetOil (CurrentOrigin%ObjOil,                                 &
                    !               Me%EulerModel(em)%SpreadingVelocityY,         &
                    !               STAT = STAT_CALL)
                    !if (STAT_CALL /= SUCCESS_) stop 'MovePartic - ModuleLagrangianGlobal - ERR60'

                endif

            endif

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr


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
        type (T_Position)                           :: NewPosition, OldPosition
        integer                                     :: NewI, NewJ, KUB
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        logical                                     :: PositionCorrected
        real                                        :: TauAux
        integer                                     :: ID_Group, emp, em, STAT_CALL
        real                                        :: Radius, Area
        real                                        :: VolOld, VolNew, dVol
        real                                        :: Modulus, WindX, WindY
        real                                        :: GradDWx, GradDWy, Aux
        logical                                     :: NoIntU, NoIntV, ComputeTrajectory, NotChangeDomain, HaveDomain, MovePartic





        CurrentPartic => CurrentOrigin%FirstPartic
CP:     do while (associated (CurrentPartic))

            ComputeTrajectory = .true. 

CT:         do while (ComputeTrajectory) 

                emp    = CurrentPartic%Position%ModelID

                WS_ILB = Me%EulerModel(emp)%WorkSize%ILB
                WS_IUB = Me%EulerModel(emp)%WorkSize%IUB
                WS_JLB = Me%EulerModel(emp)%WorkSize%JLB
                WS_JUB = Me%EulerModel(emp)%WorkSize%JUB
            
                !Grid Cell of the particle
                i        = CurrentPartic%Position%I
                j        = CurrentPartic%Position%J
                k        = CurrentPartic%Position%K

                ID_Group = CurrentOrigin%GroupID

                emp       = CurrentPartic%Position%ModelID

                !Shorten Var
                if (Me%OverLay) then
                    Velocity_U => Me%EulerModel(emp)%OverLay%VelUFinal
                    Velocity_V => Me%EulerModel(emp)%OverLay%VelVFinal
                else
                    Velocity_U => Me%EulerModel(emp)%Velocity_U
                    Velocity_V => Me%EulerModel(emp)%Velocity_V
                endif

BD:             if (CurrentPartic%Beached .or. CurrentPartic%Deposited) then

                    if (CurrentPartic%Deposited) then

                        TauAux = CurrentPartic%TauErosion + (Me%EulerModel(emp)%Lag2Euler%TauErosionGrid( i, j, ID_Group) - &
                                 CurrentPartic%TauErosion) * Me%DT_Partic / CurrentOrigin%Deposition%Tdecay

                        CurrentPartic%TauErosion = max (CurrentOrigin%Deposition%TauErosion, TauAux)


                        if (VerifyRessuspension(CurrentPartic)) then
                            CurrentPartic%Deposited = .false. 
                            !The effect of other sediment classes ends when the particle is ressuspended
                            CurrentPartic%TauErosion = CurrentOrigin%Deposition%TauErosion
                        endif

                    endif

                    ComputeTrajectory = .false.
                    
                    exit
            
                endif BD
                
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


                !Spagnol et al. (Mar. Ecol. Prog. Ser., 235, 299-302, 2002).                        
                !Linear Interpolation to obtain the thickness gradient
                GradDWx = LinearInterpolation(Me%EulerModel(emp)%DWZ_Xgrad(i,  j  ,k),  &
                                              Me%EulerModel(emp)%DWZ_Xgrad(i,  j+1,k), Balx)
                GradDWy = LinearInterpolation(Me%EulerModel(emp)%DWZ_Ygrad(i,  j  ,k),  &
                                              Me%EulerModel(emp)%DWZ_Ygrad(i+1,j  ,k), Baly)

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

                        WindX = Me%EulerModel(emp)%WindX(i, j)
                        WindY = Me%EulerModel(emp)%WindY(i, j)

                    endif


                    UWind = CurrentOrigin%Movement%WindTransferCoef * WindX
                        
                    VWind = CurrentOrigin%Movement%WindTransferCoef * WindY

                    !Plume velocity
                    UPlume = 0.
                    VPlume = 0.

                    if (CurrentOrigin%State%Oil) then

                        if      (SpreadingMethod == ThicknessGradient) then
                
                            UOil = LinearInterpolation(                                      &
                                    Me%EulerModel(emp)%OilSpreading(ID_Group)%VelocityX(i, j  ),    &
                                    Me%EulerModel(emp)%OilSpreading(ID_Group)%VelocityX(i, j+1),    &
                                    Balx)

                            VOil = LinearInterpolation(                                      &
                                    Me%EulerModel(emp)%OilSpreading(ID_Group)%VelocityY(i, j  ),    &
                                    Me%EulerModel(emp)%OilSpreading(ID_Group)%VelocityY(i+1, j),    &
                                    Baly)
                                            
                        elseif  (SpreadingMethod == Fay              ) then

                            call RANDOM_NUMBER(R1)
                            call RANDOM_NUMBER(R2)

                            UOil = R1 * cos(2 * PI * R2) * Me%ExternalVar%DiffVelocity
                            VOil = R1 * sin(2 * PI * R2) * Me%ExternalVar%DiffVelocity

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
                    Layers = Me%EulerModel(emp)%WorkSize%KUB -                        &
                             Me%EulerModel(emp)%WorkSize%KLB + 1 

                    if (Layers > 1) then !More than 1 layer

                        EspSup = Me%EulerModel(emp)%DWZ(i, j, k+1)
                        Esp    = Me%EulerModel(emp)%DWZ(i, j, k  )
                        EspInf = Me%EulerModel(emp)%DWZ(i, j, k-1)

                        NoIntU = .false.
                        NoIntV = .false. 

                        if (Me%EulerModel(emp)%ComputeFaces3D_U(i, j,   k) /= Covered .or.  &
                            Me%EulerModel(emp)%ComputeFaces3D_U(i, j+1, k) /= Covered ) then
                           !No vertical interpolation
                            NoIntU = .true.
                        endif

                        if (Me%EulerModel(emp)%ComputeFaces3D_V(i, j, k  ) /= Covered .or.  &
                            Me%EulerModel(emp)%ComputeFaces3D_V(i+1, j, k) /= Covered ) then
                           !No vertical interpolation
                            NoIntV = .true.
                        endif


                        !Not Close to the bottom
                        if (k /= Me%EulerModel(emp)%kFloor(i, j) .and. BALZ < 0.5) then
                    
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
                        else if ((K /= Me%EulerModel(emp)%WorkSize%KUB) .and. (BALZ > 0.5 )) then

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
                        TlagrangeX = Me%EulerModel(emp)%MixingLengthX(i, j, k) /      &
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
                       
                        !Spagnol et al. (Mar. Ecol. Prog. Ser., 235, 299-302, 2002).                        
                        Aux                      = (1.0 - 2.0 * RAND)
                        if (Aux >= GradDWx)   UD = - UD                                  
                    else
                        UD                       = CurrentPartic%UD_old                               
                        CurrentPartic%TpercursoX = CurrentPartic%TpercursoX + Me%DT_Partic
                    end if


                    if (VStandardDeviation > 0.0) then
                        TlagrangeY = Me%EulerModel(emp)%MixingLengthY(i, j, k) /      &
                                     VStandardDeviation
                    else
                        TlagrangeY = 0.0
                    endif

                    if (CurrentPartic%TpercursoY >= TlagrangeY) then
                        call random_number(RAND)
                        !VD                       = (1.0 - 2.0 * RAND) * 1.732050808 *        &
                        !                           VStandardDeviation  

                        !Spagnol et al. (Mar. Ecol. Prog. Ser., 235, 299-302, 2002).                        
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

                OldPosition   = CurrentPartic%Position
        
                !New Position
                NewPosition%X = CurrentPartic%Position%X + DX
                NewPosition%Y = CurrentPartic%Position%Y + DY

                NotChangeDomain = GetXYInsideDomain(Me%EulerModel(emp)%ObjHorizontalGrid, &
                NewPosition%X, NewPosition%Y, Referential= AlongGrid_, STAT = STAT_CALL) 

                if (STAT_CALL /= SUCCESS_) stop 'Locate_ModelDomain - ModuleLagrangianGlobal - ERR10'


iNCD:           if (NotChangeDomain) then

                    ComputeTrajectory   = .false.
                    HaveDomain          = .true.
                    NewPosition%ModelID = emp

                else iNCD

                    HaveDomain = .false. 

d1:                 do em = emp+1, Me%EulerModelNumber
                    
                        HaveDomain = GetXYInsideDomain(Me%EulerModel(em)%ObjHorizontalGrid, &
                                       CurrentPartic%Position%CoordX,                       &
                                       CurrentPartic%Position%CoordY,                       &
                                       Referential= GridCoord_, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Locate_ModelDomain - ModuleLagrangianGlobal - ERR10'

                    
iHD:                    if (HaveDomain) then

                            !Convert Coordinates
                            call Convert_XY_CellIJ(Me%EulerModel(em),CurrentPartic%Position, Referential = GridCoord_)

iOpen2D:                    if  (Me%EulerModel(em)%OpenPoints3D(CurrentPartic%Position%i, CurrentPartic%Position%j, &
                                                                Me%EulerModel(em)%WorkSize%KUB) == OpenPoint) then


                                call Convert_Z_CellK (CurrentOrigin, Me%EulerModel(em), CurrentPartic%Position, PositionCorrected)
                                call Convert_CellK_K (                                  CurrentPartic%Position)

                                CurrentPartic%Position%ModelID = em

                                ComputeTrajectory = .true.

                            else iOpen2D

                                !If a particle doesnt move the freazed state is ON
                                CurrentPartic%Freazed    = ON

                                CurrentPartic%TpercursoX = abs(null_real)
                                CurrentPartic%TpercursoY = abs(null_real)

                                ComputeTrajectory = .false.

                            endif iOpen2D

                            exit
                        endif iHD

                    enddo d1

                    if (.not. HaveDomain) then                            
                        ComputeTrajectory        = .false. 
                        CurrentPartic%KillPartic = ON 
                    endif
            
                endif iNCD

            enddo CT

            MovePartic = .true.

            if (CurrentPartic%Beached           .or.                                    &
                CurrentPartic%Deposited         .or.                                    &
                CurrentPartic%Freazed           .or.                                    &
                CurrentPartic%KillPartic) then
                MovePartic = .false.
            endif 
            
            if (MovePartic .and. NewPosition%ModelID < 0) then
                MovePartic = .false.
            endif
            

iFKP:       if (MovePartic) then

                !Convert Coordinates
                call Convert_XY_CellIJ(Me%EulerModel(NewPosition%ModelID),NewPosition, Referential = AlongGrid_)

                !Verifies new position
                NewI = NewPosition%i
                NewJ = NewPosition%j

                KUB  = Me%EulerModel(NewPosition%ModelID)%WorkSize%KUB                   

                !If it isnt a OpenPoint, donesnt move, reset TPercurso
ie:             if  (Me%EulerModel(NewPosition%ModelID)%OpenPoints3D(NewI, NewJ, KUB) /= OpenPoint) then

                    !If a particle doesnt move the freazed state is ON
                    CurrentPartic%Freazed    = ON

                    CurrentPartic%TpercursoX = abs(null_real)
                    CurrentPartic%TpercursoY = abs(null_real)

                !Moves it 
                else ie

                    call MoveParticVertical  (CurrentOrigin, CurrentPartic, NewPosition, VelModH)

                    call Convert_Z_CellK (CurrentOrigin, Me%EulerModel(NewPosition%ModelID), &
                                                         NewPosition, PositionCorrected)
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

                endif ie

            endif iFKP

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
        integer                                     :: i, j, emp

        i = CurrentPartic%Position%I
        j = CurrentPartic%Position%J
        emp= CurrentPartic%Position%ModelID

        VerifyRessuspension = OFF

        VerifyShearStress   = OFF
        VerifyErosionRate   = OFF

        !Partheniades, E., 1965. Erosion and deposition of cohesive soils.
        !J. Hydr. Div., ASCE 91 HY1 , 105–139.
     
        !The ressuspension of a tracer is function of the bottom shear stress
        !and of the erosion rate. The erosion rate (Er) is quantifiied in the form of
        !a probability that is equal to  min (1, Er * dt / MassSedTotal(i,j))
cd1:    if (Me%EulerModel(emp)%BottomStress(CurrentPartic%Position%I,             &
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
        integer                                     :: i, j, kbottom, emp

        i = CurrentPartic%Position%I
        j = CurrentPartic%Position%J
        !k = CurrentPartic%Position%K
        emp= CurrentPartic%Position%ModelID

        kbottom = Me%EulerModel(emp)%KFloor(i, j)

        VerifyDeposition = OFF

cd1:    if (CurrentOrigin%Deposition%BottomDistance    >=                               &
            (Me%EulerModel(emp)%SZZ(i, j, kbottom -1)- CurrentPartic%Position%Z)) then

!Odd, N.V.M., Owen, M.W., 1972. A two-layer model of mud transport in the Thames estuary. 
!In: Proceedings. Institution of Civil Engineers, London, pp. 195-202.

cd2:        if (Me%EulerModel(emp)%BottomStress(i,j) <                                      &
                CurrentOrigin%Deposition%TauDeposition) then
                            
                call random_number(Ranval)

                Aux = (CurrentOrigin%Deposition%TauDeposition -                         &
                       Me%EulerModel(emp)%BottomStress(i,j)) /                              &
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
        integer                                     :: i, j, k, KUB, emp
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

        emp      = CurrentPartic%Position%ModelID

        CellK   = CurrentPartic%Position%CellK
        KUB     = Me%EulerModel(emp)%WorkSize%KUB
        BALZ    = CellK - int(CellK)


MF:     if (CurrentOrigin%Movement%Float) then

            NewPosition%Z = Me%EulerModel(emp)%SZZ(i, j, KUB)

        else MF

SA:         if (CurrentOrigin%Movement%MovType == SullivanAllen_) then

                CompZ_Up   = Me%EulerModel(emp)%Lupward  (i, j, k)
                CompZ_Down = Me%EulerModel(emp)%Ldownward(i, j, k)

                if ((CompZ_Up < 0.0) .OR. (CompZ_Down < 0.0)) then
                    CompZ_Up   = 0.0
                    CompZ_Down = 0.0
                end if
    
                if ((Me%EulerModel(emp)%WorkSize%KUB -                                     &
                     Me%EulerModel(emp)%WorkSize%KLB + 1) > 1) then 

                    EspSup = Me%EulerModel(emp)%DWZ(i, j, k+1)
                    Esp    = Me%EulerModel(emp)%DWZ(i, j, k  )
                    EspInf = Me%EulerModel(emp)%DWZ(i, j, k-1)

                    if ((K /= Me%EulerModel(emp)%WorkSize%KLB) .AND.              &
                        (BALZ < 0.5)) then       !Close to the bottom

                        CompZ1_Up   = Me%EulerModel(emp)%Lupward  (i, j, k-1)
                        CompZ1_Down = Me%EulerModel(emp)%Ldownward(i, j, k-1)

                        AuxCompMisturaZ_Up   = 2.0 * (CompZ1_Up * (0.5 - BALZ) * Esp +   &
                                               CompZ_Up * BALZ * Esp + CompZ_Up * 0.5 *  &
                                               EspInf) / (Esp + EspInf)

                        AuxCompMisturaZ_Down = 2.0 * (CompZ1_Down * (0.5 - BALZ) * Esp + &
                                               CompZ_Down * BALZ * Esp + CompZ_Down *    &
                                               0.5 * EspInf) / (Esp + EspInf)
                 
                    else if ((K /= Me%EulerModel(emp)%WorkSize%KUB) .AND.                  &
                             (BALZ .GT. 0.5 )) then !Close to the surface

                        CompZ1_Up   = Me%EulerModel(emp)%Lupward  (i, j, k+1)
                        CompZ1_Down = Me%EulerModel(emp)%Ldownward(i, j, k+1)

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

                    WStandardDeviation = 1.0975 * Me%EulerModel(emp)%ShearVelocity(i, j) 

                    WD = WD_(CurrentPartic, WStandardDeviation, AuxCompMisturaZ_Up,      &
                             AuxCompMisturaZ_Down, Me%DT_Partic)

                end select

            else if (CurrentOrigin%Movement%MovType .EQ. NotRandom_    ) then MT

                WD = 0.0

            end if MT


            if (CurrentOrigin%Movement%Advection) then

                if (Me%EulerModel(emp)%WorkSize%KUB -                             &
                    Me%EulerModel(emp)%WorkSize%KLB + 1 == 1) then

                    W = 0.

                else

                    W = Me%EulerModel(emp)%Velocity_W(i, j, k  ) * (1.0 - BALZ) + &
                        Me%EulerModel(emp)%Velocity_W(i, j, k+1) * BALZ 

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
                DeltaD = Me%EulerModel(emp)%SigmaDensity(i, j, k) - CurrentPartic%SigmaDensity

                !Buoyancy is consider null for low density gradients 
                if (abs(DeltaD) > 0.5) then
                    ai     = Gravity * DeltaD / (Me%EulerModel(emp)%SigmaDensity(i, j, k) + &
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

                dwt = Me%EulerModel(emp)%DWZ(i, j, k) / Me%DT_Partic

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
        integer                                     :: i, j, k, emp
        real                                        :: CellI, CellJ
        real                                        :: Balx, Baly
        real                                        :: U, V, VelQ
        real                                        :: FKK
        real                                        :: VolOld
        real                                        :: PlumeCoef, RelativeVel, VarianceTurb
        real, dimension(:, :, :), pointer           :: Velocity_U, Velocity_V

      

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

                        emp    = CurrentPartic%Position%ModelID

                        !Shorten Var
                        if (Me%OverLay) then
                            Velocity_U => Me%EulerModel(emp)%OverLay%VelUFinal
                            Velocity_V => Me%EulerModel(emp)%OverLay%VelVFinal
                        else
                            Velocity_U => Me%EulerModel(emp)%Velocity_U
                            Velocity_V => Me%EulerModel(emp)%Velocity_V
                        endif

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
                                      (Me%EulerModel(emp)%WaterColumn(i, j)/10.0) *       &  !(H/HREF)
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
        integer                                     :: nProp, emp
        real                                        :: OldVolume

        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%State%VariableGeom) then

                nProp = 1
                CurrentProperty => CurrentOrigin%FirstProperty
dw2:            do while (associated(CurrentProperty))

                    CurrentPartic => CurrentOrigin%FirstPartic
dw3:                do while (associated(CurrentPartic))

DB:                     if (.not. CurrentPartic%Deposited .and.                         &
                            .not. CurrentPartic%Beached   .and.                         &
                            .not. CurrentPartic%Freazed        ) then

                            emp = CurrentPartic%Position%ModelID

                            call GetAmbientConcentration (CurrentProperty,              &
                                                          emp,                          & 
                                                          CurrentPartic%Position,       &
                                                          Concentration)

                            if (CurrentPartic%Geometry%VolVar < 0) then

                                stop 'Dilution - ModuleLagrangianGlobal - ERR01'

                            endif

                            OldVolume = CurrentPartic%Geometry%Volume -                 &
                                        CurrentPartic%Geometry%VolVar

                            !if (nprop <1 .or. nprop > 4 .or. CurrentPartic%Geometry%Volume <0) then
                            !    write(*,*) CurrentOrigin%name
                            !    write(*,*) CurrentProperty%ID
                            !    write(*,*) CurrentPartic%ID
                            !    write(*,*) nprop
                            !    write(*,*) CurrentPartic%Geometry%Volume 
                            !endif

                            !Calculates New Concentration
                            CurrentPartic%Concentration(nProp) =                        &
                                (CurrentPartic%Concentration(nProp) * OldVolume +       &
                                 Concentration * CurrentPartic%Geometry%VolVar) /       &
                                 CurrentPartic%Geometry%Volume

                        endif DB
                    
                        CurrentPartic => CurrentPartic%Next
                    
                    enddo dw3
        
                    nProp = nProp + 1
                    CurrentProperty => CurrentProperty%Next
                    
                enddo dw2

            endif
            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr


    end subroutine Dilution

    !--------------------------------------------------------------------------

    subroutine GetAmbientConcentration (Property, ModelID, Position, Concentration)

        !Arguments-------------------------------------------------------------
   
        type (T_Property), pointer                  :: Property
        integer                                     :: ModelID
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

            call GetConcentration(Me%EulerModel(ModelID)%ObjWaterProperties,            &
                                  ConcentrationX    = ConcentrationX,                   &
                                  PropertyXIDNumber = Property%ID,                      &
                                  STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                
                write(*,*)'AmbientConcentration not given for : ', trim(Property%Name)
                stop 'GetAmbientConcentration - ModuleLagrangianGlobal - ERR10'

            endif

            Concentration = ConcentrationX(Position%I, Position%J, Position%K)

            call UngetWaterProperties(Me%EulerModel(ModelID)%ObjWaterProperties,        &
                                      ConcentrationX,                                   &
                                      STAT              = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                
                stop 'GetAmbientConcentration - ModuleLagrangianGlobal - ERR20'

            endif

        endif


    end subroutine GetAmbientConcentration

    !--------------------------------------------------------------------------

    subroutine VerifyLandParticles ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        integer                                     :: i, j, k, emp


        !Set KILL state of all particle which are in intertidal areas to ON

        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))
                
            if (CurrentOrigin%Movement%KillLandParticles) then

                CurrentPartic  => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J
                    k = CurrentPartic%Position%K
                    emp= CurrentPartic%Position%ModelID

                    if (Me%EulerModel(emp)%OpenPoints3D(i, j, k)  /= OpenPoint .and. &
                        Me%EulerModel(emp)%Waterpoints3D(i, j, k) == WaterPoint) then

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

    subroutine KillParticInsideBoxes ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: BoxOrigin, CurrentOrigin   
        type (T_Partic), pointer                    :: CurrentPartic 
        logical                                     :: Inside
        integer                                     :: em, STAT_CALL   
        !Begin-----------------------------------------------------------------


        BoxOrigin => Me%FirstOrigin
BoxOr:  do while (associated(BoxOrigin))

            if (BoxOrigin%KillPartInsideBox) then

                if (Me%Now >= BoxOrigin%InstantEmission .and. BoxOrigin%OnlyOnceEmit) then

                    CurrentOrigin => Me%FirstOrigin
CurrOr:             do while (associated(CurrentOrigin))                    
                
                        CurrentPartic  => CurrentOrigin%FirstPartic
                        do while (associated(CurrentPartic))

                            em = CurrentPartic%Position%ModelID

                            Inside = CheckIfInsideBox (Me%EulerModel(em)%ObjBoxDif,     &
                                                       BoxOrigin%BoxNumber,             &
                                                       CurrentPartic%Position%CoordX,   &
                                                       CurrentPartic%Position%CoordY,   &
                                                       STAT = STAT_CALL)
                
                            if (STAT_CALL /= SUCCESS_) stop 'KillParticInsideBoxes - ModuleLagrangianGlobal - ERR10'

                            !Deletes all Particle with KillPartic State on
                            if (Inside) CurrentPartic%KillPartic = ON

                            CurrentPartic => CurrentPartic%Next


                        enddo

                        CurrentOrigin => CurrentOrigin%Next

                    enddo CurrOr

                endif

            endif

            BoxOrigin => BoxOrigin%Next

        enddo BoxOr

    end subroutine KillParticInsideBoxes

       
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine LightEvolution ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentParticle
        type (T_Property), pointer                  :: CurrentProperty

        real,    dimension(:, :, :   ), pointer     :: Concentration
        real                                        :: UnitsCoef, dh1, CenterRadiation, SWPercentage
        integer                                     :: ig, STAT_CALL, iProp, i, j, k, em
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB, Kbottom
        logical                                     :: NeedsParameters, NeedsConcentrations, &
                                                       NeedsPhyto, NeedsSPM, FoundPhyto, FoundSed, ValidProp
        !Begin-----------------------------------------------------------------


        !Updates light
d1:     do em     = 1, Me%EulerModelNumber

            !Shorten
            ILB    = Me%EulerModel(em)%WorkSize%ILB
            IUB    = Me%EulerModel(em)%WorkSize%IUB
            JLB    = Me%EulerModel(em)%WorkSize%JLB
            JUB    = Me%EulerModel(em)%WorkSize%JUB
            KLB    = Me%EulerModel(em)%WorkSize%KLB
            KUB    = Me%EulerModel(em)%WorkSize%KUB


            !Allocates auxiliar variable
            allocate (Concentration (Me%EulerModel(em)%Size%ILB:Me%EulerModel(em)%Size%IUB, &
                                     Me%EulerModel(em)%Size%JLB:Me%EulerModel(em)%Size%JUB, &
                                     Me%EulerModel(em)%Size%KLB:Me%EulerModel(em)%Size%KUB))
            
d2:         do ig = 1, Me%nGroups
        
i0:             if(Me%EulerModel(em)%Light(ig)%Compute)then

                    call GetLightExtinctionOptions(LightExtinctionID    = Me%EulerModel(em)%Light(ig)%ObjLightExtinction,& 
                                                   NeedsParameters      = NeedsParameters,       &
                                                   NeedsConcentrations  = NeedsConcentrations,   &
                                                   NeedsPhyto           = NeedsPhyto,            &
                                                   NeedsSPM             = NeedsSPM,              &
                                                   STAT                 = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'LightEvolution - ModuleLagrangianGlobal - ERR10'

                    call GetRadiationPercentages (Me%EulerModel(em)%Light(ig)%ObjLightExtinction,           &
                                                  SWPercentage = SWPercentage, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)stop 'LightEvolution - ModuleLagrangianGlobal - ERR20'


    i1:             if(NeedsConcentrations)then

                        if (Me%Now > Me%ExternalVar%LastConcCompute) call FillGridConcentration

                        if (associated(Me%FirstOrigin)) then
                            CurrentOrigin => Me%FirstOrigin 
                        else
                            CurrentOrigin => Me%OriginDefault
                        endif

    Catch:              do while (associated(CurrentOrigin))
                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
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

                            if (CurrentProperty%ID == Phytoplankton_ .or. CurrentProperty%ID == Cohesive_Sediment_) then
                                ValidProp = .true.
                            else
                                ValidProp = .false.
                            endif

    i2:                     if ((NeedsPhyto.and.FoundPhyto).or.(NeedsSPM.and.FoundSed).and.ValidProp) then

                                UnitsCoef = 1e-3
                         
                                Concentration(:,:,:) = Me%EulerModel(em)%Lag2Euler%GridConc(:,:,:,iProp,ig)

    i3:                         if(NeedsParameters)then

                            
                                    call ModifyLightExtinctionField                            &
                                        (LightExtinctionID   = Me%EulerModel(em)%Light(ig)%ObjLightExtinction,&
                                         WaterPoints3D       = Me%EulerModel(em)%Waterpoints3D,&
                                         CurrentTime         = Me%Now,             &
                                         PropertyID          = CurrentProperty%ID,             &
                                         Concentration       = Concentration,                  &
                                         UnitsCoef           = UnitsCoef,                      &
                                         ExtinctionParameter = CurrentProperty%ExtinctionParameter,&
                                         STAT                = STAT_CALL)
                                    if (STAT_CALL/= SUCCESS_) stop 'LightEvolution - ModuleLagrangianGlobal - ERR30'

                                else i3

                                    call ModifyLightExtinctionField                            &
                                        (LightExtinctionID   = Me%EulerModel(em)%Light(ig)%ObjLightExtinction,&
                                         WaterPoints3D       = Me%EulerModel(em)%Waterpoints3D,&
                                         CurrentTime         = Me%Now,             &
                                         PropertyID          = CurrentProperty%ID,             &
                                         Concentration       = Concentration,                  &
                                         UnitsCoef           = UnitsCoef,                      &
                                         STAT                = STAT_CALL)
                                    if (STAT_CALL/= SUCCESS_) stop 'LightEvolution - ModuleLagrangianGlobal - ERR40'

                                end if i3


                            endif i2

                            CurrentProperty => CurrentProperty%Next
                            iProp = iProp + 1

                        enddo dw1

                        if (NeedsPhyto.and..not.FoundPhyto) then

                            stop 'LightEvolution - ModuleLagrangianGlobal - ERR50'

                        endif
                
                        if (NeedsSPM.and..not. FoundSed) then

                            stop 'LightEvolution - ModuleLagrangianGlobal - ERR60'

                        endif

                    else i1

                        call ModifyLightExtinctionField(LightExtinctionID   = Me%EulerModel(em)%Light(ig)%ObjLightExtinction,&
                                                        WaterPoints3D       = Me%EulerModel(em)%Waterpoints3D,   &
                                                        CurrentTime         = Me%Now,             &
                                                        STAT                = STAT_CALL)
                        if (STAT_CALL/= SUCCESS_) stop 'LightEvolution - ModuleLagrangianGlobal - ERR70'

                    end if i1

                end if i0



            call GetShortWaveExtinctionField(Me%EulerModel(em)%Light(ig)%ObjLightExtinction, &
                   Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'LightEvolution - ModuleLagrangianGlobal - ERR80'

d3:         do j = JLB, JUB
d4:         do i = ILB, IUB

                
i5:             if (Me%EulerModel(em)%Waterpoints3D(i, j, KUB) == WaterPoint) then 
                    
                    kbottom = Me%EulerModel(em)%KFloor(i, j)

                    Me%EulerModel(em)%Light(ig)%TopRadiationCells(i, j, KUB) =          &
                        Me%EulerModel(em)%SurfaceRadiation(i, j) * SWPercentage
                   
d5:                 do k= KUB-1, kbottom - 1, -1

                        Me%EulerModel(em)%Light(ig)%TopRadiationCells(i, j, k) =        &
                            Me%EulerModel(em)%Light(ig)%TopRadiationCells   (i,j,k+1) * &
                            exp(-Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField(i,j,k+1)* &
                            Me%EulerModel(em)%DWZ(i, j, k + 1))

                        Me%EulerModel(em)%Light(ig)%TopRadiationCells(i, j, k)=         &
                            max(Me%EulerModel(em)%Light(ig)%TopRadiationCells(i, j, k), 0.001)

                    end do d5

                end if i5

            end do d4
            end do d3

        enddo d2

        deallocate(Concentration)

        enddo d1


        CurrentOrigin => Me%FirstOrigin 
dw2:    do while (associated(CurrentOrigin))
                        
            ig = CurrentOrigin%GroupID

            !It is admited that the tracers thickness is always smaller than the cells thickness
            CurrentParticle => CurrentOrigin%FirstPartic
dw3:            do while (associated(CurrentParticle))

                i  = CurrentParticle%Position%I
                j  = CurrentParticle%Position%J
                k  = CurrentParticle%Position%K

                em = CurrentParticle%Position%ModelID


                dh1 = CurrentParticle%Position%Z - Me%EulerModel(em)%SZZ(i, j, k)

                !Radiation in the center of the particle (I2)
                !I2 = I1*exp(-K * Thickness/2) <=> I1 = I2 / exp(-K * Thickness/2)
                !I1 = top particle face
                !I3 = bottom particle face
                !I3 = I2*exp(-K * Thickness/2)
                CenterRadiation =  Me%EulerModel(em)%Light(ig)%TopRadiationCells(i, j, k) *    &
                                   exp( - Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField(i,j,k) * dh1)

                !dh2 = CurrentParticle%Geometry%Thickness

                !Iaverage = (I1 - I3) / (K * Thickness) = I2 / (K * Thickness) * 
                !           (1/exp(-K * Thickness/2) - exp(-K * Thickness/2))
                !Iaverage = I2 / (K * Thickness) * (exp(K * Thickness/2) - exp(-K * Thickness/2)) 
                !AverageRadiation = CenterRadiation / Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField(i,j,k) / dh2 * &
                !                   (exp(  Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField(i,j,k) * dh2/2)         - &
                !                    exp(- Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField(i,j,k) + dh2/2)) 

                CurrentParticle%Radiation    = CenterRadiation  
                CurrentParticle%ShortWaveExt = Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField(i, j, k)

                CurrentParticle => CurrentParticle%Next
            enddo dw3


            CurrentOrigin => CurrentOrigin%Next
        enddo dw2


d10:    do em     = 1, Me%EulerModelNumber
        do ig = 1, Me%NGroups
            call UnGetLightExtinction(Me%EulerModel(em)%Light(ig)%ObjLightExtinction, &
                                      Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'LightEvolution - ModuleLagrangianGlobal - ERR100'
        enddo
        enddo d10

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
        integer                                     :: nParticle
        integer                                     :: i, j, k, nProp, em, emp
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
        Actual = Me%Now

d1:     do em =1, Me%EulerModelNumber

            !Gets the temperature and the Salinity from the Eulerian model
            call GetConcentration(Me%EulerModel(em)%ObjWaterProperties,                 &
                                  ConcentrationX    = Me%EulerModel(em)%Temperature3D,  &
                                  PropertyXIDNumber = Temperature_,                     &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR10'

            call GetConcentration(Me%EulerModel(em)%ObjWaterProperties,                 &
                                  ConcentrationX    = Me%EulerModel(em)%Salinity3D,     &
                                  PropertyXIDNumber = Salinity_,                        &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR03'


            !Gets the FishFood from the Eulerian model
            if (Me%State%Larvae) then
                call GetConcentration(Me%EulerModel(em)%ObjWaterProperties,                        &
                                      ConcentrationX    = Me%EulerModel(em)%FishFood3D,                         &
                                      PropertyXIDNumber = FishFood_,                          &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR003'
            end if

        enddo d1

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
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR01'


                !Allocates the temporary matrixes
                allocate (TemperatureX        (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR05'

                allocate (SalinityX           (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR06'

                allocate (SWRadiationX        (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR07'

                allocate (LightExtCoefX       (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR07a'

                allocate (Thickness           (CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR07b'

                allocate (WQMMass(CurrentOrigin%nProperties, CurrentOrigin%nParticle), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR08'

                if (CurrentOrigin%State%Larvae)  then            
                    allocate (FishFoodX       (CurrentOrigin%nParticle), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR008'
                end if


                CurrentPartic  => CurrentOrigin%FirstPartic
                nParticle      = 1
                do while (associated(CurrentPartic))

                    i = CurrentPartic%Position%I
                    j = CurrentPartic%Position%J
                    k = CurrentPartic%Position%K

                    emp= CurrentPartic%Position%ModelID

                    !Stores Temperature
                    TemperatureX (nParticle) = Me%EulerModel(emp)%Temperature3D (i, j, k)

                    !Stores Salinity
                    SalinityX    (nParticle) = Me%EulerModel(emp)%Salinity3D    (i, j, k)

                   !Stores FishFood
                    if (CurrentOrigin%State%Larvae)                                        &            
                        FishFoodX (nParticle) = Me%EulerModel(emp)%FishFood3D (i, j, k)
                    

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
                            WQMMass(PhytoWQM,            nParticle) = &
                            CurrentPartic%Concentration (nProp)
                                                            
                        case (Zooplankton_                    )
                            WQMMass(ZooWQM,              nParticle) = &
                            CurrentPartic%Concentration (nProp)
                                                            
                        case (Larvae_                         )
                            WQMMass(LarvaeWQM,           nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (Age_                         )
                            WQMMass(AgeWQM,              nParticle) = &
                            CurrentPartic%Concentration (nProp)                                                           
 
                        case (PON_                            )
                            WQMMass(PartOrganicNitrogenWQM, nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (DONRefractory_                  )
                            WQMMass(DONRefractoryWQM,    nParticle) = &
                            CurrentPartic%Concentration (nProp)
                           
                        case (DONNon_Refractory_              )
                            WQMMass(DONNonRefractoryWQM, nParticle) = &
                            CurrentPartic%Concentration (nProp)
                           
                        case (Ammonia_                        )
                            WQMMass(AmmoniaWQM,          nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (Nitrate_                        )
                            WQMMass(NitrateWQM,          nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (Nitrite_                        )
                            WQMMass(NitriteWQM,          nParticle) = &
                            CurrentPartic%Concentration (nProp)
     
                        case (BOD_                            )
                            WQMMass(BODWQM,              nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (Oxygen_                         )
                            WQMMass(OxygenWQM,           nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (Ciliate_                        )
                            WQMMass(CiliateWQM,          nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (Bacteria_                       )
                            WQMMass(BacteriaWQM,         nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (POP_                            )
                            WQMMass(PartOrganicPhosphorusWQM,nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (DOPRefractory_                  )
                            WQMMass(DOPRefractoryWQM,    nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (DOPNon_Refractory_              )
                            WQMMass(DOPNonRefractoryWQM, nParticle) = &
                            CurrentPartic%Concentration (nProp)

                        case (Inorganic_Phosphorus_           )
                            WQMMass(InorganicPhosphorusWQM, nParticle) = &
                            CurrentPartic%Concentration (nProp)
    
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
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR09a'

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
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR09b'

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
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR10'
                nullify (TemperatureX) 

                deallocate (SalinityX,              STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR11'
                nullify (SalinityX)

                deallocate (SWRadiationX,           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR12'
                nullify (SWRadiationX)

                deallocate (LightExtCoefX,          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR12'
                nullify (LightExtCoefX)

                deallocate (Thickness,               STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR12'
                nullify (Thickness)

                deallocate (WQMMass,                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR13'
                nullify (WQMMass)

                if (CurrentOrigin%State%Larvae)  then            
                    deallocate (FishFoodX,           STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR113'
                    nullify (FishFoodX) 
                end if
                

                CurrentOrigin%NextWQMCompute = CurrentOrigin%NextWQMCompute + CurrentOrigin%DTWQM

            end if

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

d2:     do em =1, Me%EulerModelNumber

            !Ungets Concentration from the eulerian module
            call UngetWaterProperties (Me%EulerModel(em)%ObjWaterProperties,            &
                                       Me%EulerModel(em)%Temperature3D,                 &
                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR14'

            call UngetWaterProperties (Me%EulerModel(em)%ObjWaterProperties,            &
                                       Me%EulerModel(em)%Salinity3D,                    &
                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR15'
       
            if (Me%State%Larvae)  then            
                call UngetWaterProperties (Me%EulerModel(em)%ObjWaterProperties,        &
                                           Me%EulerModel(em)%FishFood3D,                &
                                           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PropertiesEvolution - ModuleLagrangianGlobal - ERR115'
            end if

        enddo d2

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

P1:             if (Me%State%T90Compute) then

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
                        stop 'ColiformDecay - ModuleLagrangianGlobal - ERR01'

                    endif

                    if (iTemp == -99) then

                        write(*,*) 'The origin ='//trim(CurrentOrigin%Name)// &
                                    ' must have a temperature property to compute a variable T90'
                        stop 'ColiformDecay - ModuleLagrangianGlobal - ERR02'

                    endif


                endif P1


                nProp = 1
                CurrentProperty => CurrentOrigin%FirstProperty
CurrProp:       do while (associated(CurrentProperty))

                    if (CurrentProperty%ID == Fecal_Coliforms_ .or.CurrentProperty%ID == E_Coli_) then

                        CurrentPartic => CurrentOrigin%FirstPartic
                        do while (associated(CurrentPartic))

                            if (CurrentProperty%T90Compute) then

  
                                CurrentPartic%T90 = ComputeT90(CurrentPartic%Concentration(iSal),&
                                                               CurrentPartic%Concentration(iTemp),&
                                                               CurrentPartic%Radiation, CurrentProperty%T90Var_Method)
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

    real function ComputeT90 (Sal, Temp, Radiation, Method)

        !Arguments-------------------------------------------------------------
        real                                        :: Sal, Temp, Radiation
        integer                                     :: Method
        !Local-----------------------------------------------------------------
        real                                        :: Convert, Light

        !Begin-----------------------------------------------------------------

 


!______Mortality model selection

        If (Method .eq. Canteras) then

            ComputeT90 = ComputeT90_Canteras (Temp,Sal,Radiation)
            

        elseif (Method .eq. Chapra) then

            !Converts W in ly/hr
            Convert    = 0.086325
            Light      = Convert * Radiation

            ComputeT90 = ComputeT90_Chapra (Temp, Sal, Light)  

        else

            write (*,*) 'T90 calculation method unknown'
            stop 'ComputeT90 - ModuleLagrangianGlobal - ERR1'
        
        endif


    end function ComputeT90

    !--------------------------------------------------------------------------
    subroutine ComputeT90Matrix (em, ig, CurrentProperty, GridConc3D)

        !Arguments----------------------------------------------------------------------
        integer                                     :: em, ig
        type (T_Property), pointer                  :: CurrentProperty
        real,   dimension(:,:,:), pointer           :: GridConc3D
        !Local--------------------------------------------------------------------------
        real                                        :: Sal, Temp, Radiation, Thickness
        integer                                     :: iSal, iTemp, nProp, i, j, k, STAT_CALL
        type (T_Property), pointer                  :: AuxProperty

        !Begin--------------------------------------------------------------------------

        if (CurrentProperty%T90Compute) then

            AuxProperty => Me%OriginDefault%FirstProperty
            nProp =   1
            iSal  = -99
            iTemp = -99
            do while (associated(AuxProperty))

                if (AuxProperty%ID == Salinity_   ) iSal  = nProp
                if (AuxProperty%ID == Temperature_) iTemp = nProp

                nProp       =  nProp + 1
                AuxProperty => AuxProperty%Next

            enddo

            call GetShortWaveExtinctionField(Me%EulerModel(em)%Light(ig)%ObjLightExtinction, &
                                             Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeT90Matrix - ModuleLagrangianGlobal - ERR10'

        endif
         
        do k = Me%EulerModel(em)%WorkSize%KLB, Me%EulerModel(em)%WorkSize%KUB
        do j = Me%EulerModel(em)%WorkSize%JLB, Me%EulerModel(em)%WorkSize%JUB
        do i = Me%EulerModel(em)%WorkSize%ILB, Me%EulerModel(em)%WorkSize%IUB
            
i1:         if (Me%EulerModel(em)%Waterpoints3D (i, j, k) == WaterPoint) then
                                
i2:             if (CurrentProperty%T90Compute) then

                    Sal       = Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, iSal,  ig)

                    Temp      = Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, iTemp, ig)

                    Thickness = Me%EulerModel(em)%DWZ(i, j, k)

                    Radiation = (Me%EulerModel(em)%Light(ig)%TopRadiationCells(i, j, k) - &
                                 Me%EulerModel(em)%Light(ig)%TopRadiationCells(i, j, k-1))/&
                                (Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField(i, j, k) * Thickness)
                
                    GridConc3D(i, j, k) = ComputeT90 (Sal, Temp, Radiation, CurrentProperty%T90Var_Method)

   
                else i2

                    GridConc3D(i, j, k) =  CurrentProperty%T90       
    
                endif i2
                
            else i1
            
                    GridConc3D(i, j, k) =  FillValueReal
            
            endif i1

        enddo
        enddo
        enddo                    


        if (CurrentProperty%T90Compute) then

            call UnGetLightExtinction(Me%EulerModel(em)%Light(ig)%ObjLightExtinction, &
                                      Me%EulerModel(em)%Light(ig)%ShortWaveExtinctionField, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeT90Matrix - ModuleLagrangianGlobal - ERR20'

        endif

    end subroutine ComputeT90Matrix
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
        integer                                     :: i, j, em, ig, emp
        real                                        :: VolInic
        real                                        :: API

        do em = 1, Me%EulerModelNumber
        do ig = 1, Me%NGroups
            Me%EulerModel(em)%OilSpreading(ig)%AreaFlag(:,:) = .true. 
        enddo
        enddo

        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            CurrentOrigin%VolTotOilBeached = 0.
            CurrentOrigin%VolTotBeached    = 0.
            CurrentOrigin%AreaTotal        = 0.
            CurrentOrigin%VolumeTotal      = 0.
            CurrentOrigin%VolumeOilTotal   = 0.

            ig = CurrentOrigin%GroupID

            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                i = CurrentPartic%Position%I
                j = CurrentPartic%Position%J
                emp= CurrentPartic%Position%ModelID

                if (CurrentPartic%Beached)    then

                    CurrentOrigin%VolTotBeached     = CurrentOrigin%VolTotBeached +         &
                                                      CurrentPartic%Geometry%Volume

                    CurrentOrigin%VolTotOilBeached  = CurrentOrigin%VolTotOilBeached +      &
                                                      CurrentPartic%Geometry%VolumeOil


                else if (.NOT. CurrentPartic%Beached) then

                    if (.not. CurrentOrigin%UseTheoricArea) then
                        if (Me%EulerModel(emp)%OilSpreading(ig)%AreaFlag(i, j)) then
                            CurrentOrigin%AreaTotal = CurrentOrigin%AreaTotal +              &
                                                      Me%EulerModel(emp)%GridCellArea (i, j)
                            Me%EulerModel(emp)%OilSpreading(ig)%AreaFlag(i, j) = .false.
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


    end subroutine ComputeAreaVolume

    !--------------------------------------------------------------------------

    subroutine InternalParticOil ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic    
        real, dimension(:, :, :), pointer           :: Temperature3D
        real, dimension(:, :, :), pointer           :: SPM3D
        real, dimension(:, :   ), pointer           :: WavePeriod2D, WaveHeight2D
        integer                                     :: i, j, k, emp 
        real                                        :: WaterTemperature, WaterDensity, SPM
        real                                        :: UWIND, VWIND, Wind
        real                                        :: AtmPressure
        real                                        :: WaveHeight, WavePeriod
        real                                        :: Factor
        real                                        :: VWaterContent, MWaterContent
        real                                        :: MDispersed
        real                                        :: OilDensity
        real                                        :: OilViscosity
        real                                        :: FMEvaporated
        real                                        :: FMDispersed     
        real                                        :: AreaTotalOUT
        real                                        :: VolumeTotalOUT, VolOld
        integer                                     :: STAT_CALL

 
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%State%Oil .and. CurrentOrigin%nParticle > 0) then

                i       = CurrentOrigin%Position%I
                j       = CurrentOrigin%Position%J
                k       = CurrentOrigin%Position%K

                emp      = CurrentOrigin%Position%ModelID


                !Gets the temperature, the Density and the SPM from the Eulerian model
                call GetConcentration(Me%EulerModel(emp)%ObjWaterProperties,             &
                                      ConcentrationX    = Temperature3D,                &
                                      PropertyXIDNumber = Temperature_,                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR01'

                call GetSPM          (Me%EulerModel(emp)%ObjWaterProperties,             &
                                      SPM               = SPM3D,                        &
                                      STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR02'

#ifndef _WAVES_
                call GetWaves (WavesID    = Me%EulerModel(emp)%ObjWaves,                 &
                               WavePeriod = WavePeriod2D,                               &
                               WaveHeight = WaveHeight2D,                               &
                               STAT       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR03'
#endif
        
                UWIND   = Me%EulerModel(emp)%WindX (i, j)
                VWIND   = Me%EulerModel(emp)%WindY (i, j)
                Wind    = abs(cmplx(UWIND, VWIND))

                AtmPressure         = Me%EulerModel(emp)%AtmPressure(i, j)
                
                WaveHeight          = WaveHeight2D (i, j)

                WavePeriod          = WavePeriod2D (i, j)

                WaterTemperature    = Temperature3D             (i, j, k)
                WaterDensity        = Me%EulerModel(emp)%Density (i, j, k)
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
                                          MWaterContent         = MWaterContent,                  &
                                          MDispersed            = MDispersed,                     &
                                          OilDensity            = OilDensity,                     &
                                          OilViscosity          = OilViscosity,                   &
                                          FMEvaporated          = FMEvaporated,                   &
                                          FMDispersed           = FMDispersed,                    &
                                          VolTotOilBeached      = CurrentOrigin%VolTotOilBeached, &
                                          VolTotBeached         = CurrentOrigin%VolTotBeached,    &
                                          VolumeTotalIN         = CurrentOrigin%VolumeOilTotal,   &   
                                          VolumeTotalOUT        = VolumeTotalOUT,                 &     
                                          AreaTotal             = CurrentOrigin%AreaTotal,        & 
                                          AreaTotalOUT          = AreaTotalOUT,                   & 
                                          WaveHeight            = WaveHeight,                     &
                                          WavePeriod            = WavePeriod,                     &
                                          STAT                  = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR07'


                Me%ExternalVar%VWaterContent = VWaterContent
                Me%ExternalVar%MWaterContent = MWaterContent
                Me%ExternalVar%MDispersed    = MDispersed
                Me%ExternalVar%OilDensity    = OilDensity
                Me%ExternalVar%OilViscosity  = OilViscosity
                Me%ExternalVar%FMEvaporated  = FMEvaporated
                Me%ExternalVar%FMDispersed   = FMDispersed
                Me%ExternalVar%AreaTotal     = AreaTotalOUT

                call OilGridConcentration  (CurrentOrigin, WaveHeight, WaterDensity)       

                !Modifies OilVolume
                Factor                                  =  VolumeTotalOUT / CurrentOrigin%VolumeTotal
                CurrentPartic                           => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    if (.NOT. CurrentPartic%Beached ) then
                        
                        !Old Volume
                        VolOld  = CurrentPartic%Geometry%Volume
                    
                        !New Volume
                        CurrentPartic%Geometry%Volume    = CurrentPartic%Geometry%Volume * Factor
                        CurrentPartic%Geometry%VolumeOil = CurrentPartic%Geometry%Volume *                 &
                                                           (1 - Me%ExternalVar%VWaterContent)

                        !Volume Variation
                        CurrentPartic%Geometry%VolVar = CurrentPartic%Geometry%Volume - VolOld


                    else if (CurrentPartic%Beached) then

                        CurrentPartic%Geometry%VolVar = 0.0
         
                    end if 
         
                    CurrentPartic => CurrentPartic%Next

                enddo

            endif    


#ifndef _WAVES_
            call UnGetWaves(Me%EulerModel(emp)%ObjWaves, WavePeriod2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR08'

            call UnGetWaves(Me%EulerModel(emp)%ObjWaves, WaveHeight2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR09'
#endif
            !Ungets Concentration from the eulerian module
            call UngetWaterProperties (Me%EulerModel(emp)%ObjWaterProperties, Temperature3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR10'

            call UngetWaterProperties (Me%EulerModel(emp)%ObjWaterProperties,  SPM3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InternalParticOil - ModuleLagrangianGlobal - ERR11'

            CurrentOrigin => CurrentOrigin%Next
        enddo CurrOr

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

    subroutine CheckConcLimits ()

        !Arguments-------------------------------------------------------------
   

        !Local------------------------------------------------------------------
        type (T_Origin),   pointer                    :: CurrentOrigin
        type (T_Property), pointer                    :: CurrentProperty
        type (T_Partic),   pointer                    :: CurrentPartic
        integer                                       :: np

        !Begin------------------------------------------------------------------
        
        CurrentOrigin => Me%FirstOrigin
CurrOr: do while (associated(CurrentOrigin))

            CurrentProperty => CurrentOrigin%FirstProperty

            do np =1, CurrentOrigin%nProperties
              
                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    if (CurrentProperty%MinON) then
                        if (CurrentPartic%Concentration(np) < CurrentProperty%MinValue) then
                            CurrentPartic%KillPartic = ON
                        endif
                    endif


                    if (CurrentProperty%MaxON) then
                        if (CurrentPartic%Concentration(np) > CurrentProperty%MaxValue) then
                            CurrentPartic%KillPartic = ON
                        endif
                    endif

                    CurrentPartic => CurrentPartic%Next
                enddo

                CurrentProperty => CurrentProperty%Next

            enddo

            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

    end subroutine CheckConcLimits

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
        integer                                     :: NumberOfOrigins
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k, Box, iO
        integer                                     :: STAT_CALL
        real, dimension(:), pointer                 :: AuxReal
        integer                                     :: nProp, em 

        !Begin--------------------------------------------------------------------------

        NumberOfOrigins = Me%nOrigins

d1:     do em =1, Me%EulerModelNumber

            !Shorten
            ILB    = Me%EulerModel(em)%WorkSize%ILB
            IUB    = Me%EulerModel(em)%WorkSize%IUB
            JLB    = Me%EulerModel(em)%WorkSize%JLB
            JUB    = Me%EulerModel(em)%WorkSize%JUB
            KLB    = Me%EulerModel(em)%WorkSize%KLB
            KUB    = Me%EulerModel(em)%WorkSize%KUB
        
            !Get the boxes
            call GetBoxes(Me%EulerModel(em)%ObjMonBox, Me%EulerModel(em)%Monitor%Boxes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'MonitorParticle - ModuleLagrangianGlobal - ERR01'

            !Calculates the volume of each monitoring boxes
            Me%EulerModel(em)%Monitor%InstBoxVolume = 0.
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                Box = Me%EulerModel(em)%Monitor%Boxes(i, j, k)
                if (Box > 0 .and. Me%EulerModel(em)%OpenPoints3D(i, j, k) == OpenPoint) then
                    Me%EulerModel(em)%Monitor%InstBoxVolume (Box) = Me%EulerModel(em)%Monitor%InstBoxVolume (Box) + &
                                                                     Me%EulerModel(em)%VolumeZ   (i, j, k)
                endif
            enddo
            enddo
            enddo

            !Integrates the volume of each monitoring box
            do Box = 1, Me%EulerModel(em)%Monitor%NumberOfBoxes
                Me%EulerModel(em)%Monitor%IntgBoxVolume (Box) = Me%EulerModel(em)%Monitor%IntgBoxVolume (Box) + &
                                                                Me%EulerModel(em)%Monitor%InstBoxVolume (Box)
            enddo

        enddo d1


d2:     do em =1, Me%EulerModelNumber

            !Calculates the volume contributed from a given origin to the volume of a 
            !monitoring cell
            Me%EulerModel(em)%Monitor%InstVolumeByOrigin = 0.
            CurrentOrigin => Me%FirstOrigin
            do while (associated(CurrentOrigin))

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i  = CurrentPartic%Position%I
                    j  = CurrentPartic%Position%J
                    k  = CurrentPartic%Position%K
                    if (em /= CurrentPartic%Position%ModelID) then
                        CurrentPartic => CurrentPartic%Next
                        cycle
                    endif
               
                    !Number of box in which Particle is located
                    Box = Me%EulerModel(em)%Monitor%Boxes(i, j, k)

                    if (Box > 0 .and. Me%EulerModel(em)%OpenPoints3D(i, j, k) == OpenPoint) then
                        Me%EulerModel(em)%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) =     &
                            Me%EulerModel(em)%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) + &
                            dble(CurrentPartic%Geometry%Volume)
                    endif

                    CurrentPartic => CurrentPartic%Next
                enddo

                CurrentOrigin => CurrentOrigin%Next

            enddo

        enddo d2

MMass:  If (Me%State%MonitorPropMass) then

d3:         do em =1, Me%EulerModelNumber

                ! Calculates the mass contributed from a given origin to the volume of a 
                ! monitoring group of cells

                Me%EulerModel(em)%Monitor%InstMassByOrigin = 0.
                CurrentOrigin => Me%FirstOrigin
                do while (associated(CurrentOrigin))

                    nProp = 1
                    CurrentProperty => CurrentOrigin%FirstProperty
                    do while (associated(CurrentProperty))

                        if (CurrentProperty%Name == Me%MonitorProperty) then
                            CurrentPartic => CurrentOrigin%FirstPartic
                            do while (associated(CurrentPartic))
                                i  = CurrentPartic%Position%I
                                j  = CurrentPartic%Position%J
                                k  = CurrentPartic%Position%K

                                if (em /= CurrentPartic%Position%ModelID) then
                                    CurrentPartic => CurrentPartic%Next
                                    cycle
                                endif
           
                                !Number of box in which Particle is located
                                Box = Me%EulerModel(em)%Monitor%Boxes(i, j, k)

                                if (Box > 0 .and. Me%EulerModel(em)%OpenPoints3D(i, j, k) == OpenPoint) then
                                    Me%EulerModel(em)%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID) =     &
                                    Me%EulerModel(em)%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID) + &
                                    (CurrentPartic%Concentration(nProp) * dble(CurrentPartic%Geometry%Volume)) 
                                endif

                                CurrentPartic => CurrentPartic%Next
                            enddo

                        endif

                        nProp = nProp + 1
                        CurrentProperty => CurrentProperty%Next
                    enddo

                    CurrentOrigin => CurrentOrigin%Next
                enddo

                !Calculates Box Mass
                Me%EulerModel(em)%Monitor%InstBoxMass (:) = 0.
                do Box = 1, Me%EulerModel(em)%Monitor%NumberOfBoxes
                    CurrentOrigin => Me%FirstOrigin
                    do while (associated(CurrentOrigin))

                        if (em /= CurrentOrigin%Position%ModelID) then
                            CurrentOrigin => CurrentOrigin%Next
                            cycle
                        endif

                        Me%EulerModel(em)%Monitor%InstBoxMass(Box) =       &
                        Me%EulerModel(em)%Monitor%InstBoxMass(Box) +   &
                        Me%EulerModel(em)%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID)

                        CurrentOrigin => CurrentOrigin%Next
                    enddo
                enddo

            enddo d3

        End If MMass

d4:     do em =1, Me%EulerModelNumber

            !Integrates the values of the volume contributing from a given origin to the volume
            do Box = 1, Me%EulerModel(em)%Monitor%NumberOfBoxes
                CurrentOrigin => Me%FirstOrigin
                do while (associated(CurrentOrigin))

                    if (em /= CurrentOrigin%Position%ModelID) then
                        CurrentOrigin => CurrentOrigin%Next
                        cycle
                    endif

                    Me%EulerModel(em)%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID) =       &
                        Me%EulerModel(em)%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID) +   &
                        Me%EulerModel(em)%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID)

                    CurrentOrigin => CurrentOrigin%Next
                enddo
            enddo

        enddo d4


        If (Me%State%MonitorPropMass) then
            allocate (AuxReal (3*(NumberOfOrigins + 1)))
        Else
            allocate (AuxReal (2*(NumberOfOrigins + 1)))
        End If


d5:     do em =1, Me%EulerModelNumber

            !Writes the Time Serie
            do Box = 1, Me%EulerModel(em)%Monitor%NumberOfBoxes

                !Instant Volume Values
                AuxReal(1) = Me%EulerModel(em)%Monitor%InstBoxVolume (Box)
                CurrentOrigin => Me%FirstOrigin
                iO = 1
                do while (associated(CurrentOrigin))
                    AuxReal(iO+1) = Me%EulerModel(em)%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID)

                    iO = iO + 1
                    CurrentOrigin => CurrentOrigin%Next
                enddo


                !Integrated Values 
                AuxReal ((NumberOfOrigins + 1) + 1) = Me%EulerModel(em)%Monitor%IntgBoxVolume (Box)
                CurrentOrigin => Me%FirstOrigin
                iO = 1
                do while (associated(CurrentOrigin))
                    AuxReal ((NumberOfOrigins + 1) + 1 + iO) = Me%EulerModel(em)%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID)
                    iO = iO + 1
                    CurrentOrigin => CurrentOrigin%Next
                enddo


                If (Me%State%MonitorPropMass) then

                !Instant Mass Values
                AuxReal(2 * (NumberOfOrigins + 1) + 1) = Me%EulerModel(em)%Monitor%InstBoxMass (Box)
                CurrentOrigin => Me%FirstOrigin
                iO = 1
                do while (associated(CurrentOrigin))
                    AuxReal(2 * (NumberOfOrigins + 1) + 1 + iO) = Me%EulerModel(em)%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID)

                    iO = iO + 1
                    CurrentOrigin => CurrentOrigin%Next
                enddo

                End If

                call WriteTimeSerieLine (Me%EulerModel(em)%Monitor%ObjTimeSerie(Box), DataLine = AuxReal, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MonitorParticle - ModuleLagrangianGlobal - ERR02'
            enddo

            deallocate (AuxReal)


            !Unget The Boxes
            call UngetBoxDif(Me%EulerModel(em)%ObjMonBox, Me%EulerModel(em)%Monitor%Boxes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'MonitorParticle - ModuleLagrangianGlobal - ERR03'

        enddo d5

    end subroutine MonitorParticle
    
    !--------------------------------------------------------------------------

    subroutine MonitorEulerian
        
        !Local-----------------------------------------------------------------
        !real, dimension(:, :, :, :), pointer        :: GridMaxTracer
        real, dimension(:, :, :), pointer           :: AuxGrid3D 
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: iProp, STAT_CALL, i, j, k, em
        type (T_Property), pointer                  :: CurrentProperty

        !Begin-----------------------------------------------------------------

        !Fills Grid concentration
        if (Me%Now > Me%ExternalVar%LastConcCompute) call FillGridConcentration 

!d1:     do em =1, Me%EulerModelNumber

        ! In this version only writes and computes the concentration for the higher priority eulerian model

            em = Me%EulerModelNumber

            ILB     = Me%EulerModel(em)%Size%ILB   
            IUB     = Me%EulerModel(em)%Size%IUB   
            JLB     = Me%EulerModel(em)%Size%JLB   
            JUB     = Me%EulerModel(em)%Size%JUB   
            KLB     = Me%EulerModel(em)%Size%KLB   
            KUB     = Me%EulerModel(em)%Size%KUB   


            !Allocates auxiliar variable
            allocate (AuxGrid3D (ILB:IUB, JLB:JUB, KLB:KUB))

            iProp = 1
            CurrentProperty => Me%OriginDefault%FirstProperty
            do while (associated(CurrentProperty))

                AuxGrid3D(:,:,:) = Me%EulerModel(em)%Lag2Euler%GridConc(:, :, :,iProp, 1)

                do K = KLB, KUB
                do J = JLB, JUB
                do I = ILB, IUB
                    if(Me%EulerModel(em)%Waterpoints3D(i,j,k) == WaterPoint) then
                        Me%EulerModel(em)%EulerianMonitor%Mass(i,j,k) = AuxGrid3D(i, j, k) * &
                                                         Me%EulerModel(em)%VolumeZ(i, j, k)
                    end if
                end do
                end do
                end do

                call BoxDif(Me%EulerModel(em)%ObjEulMonBox,                             &
                            Me%EulerModel(em)%EulerianMonitor%Mass,                     &
                            "Lagrangian_"//trim(CurrentProperty%Name),                  &
                            Me%EulerModel(em)%Waterpoints3D,                            &
                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MonitorEulerian - ModuleLagrangianGlobal - ERR01'

                CurrentProperty => CurrentProperty%Next
                iProp = iProp + 1
            enddo

            !Deallocates auxiliar variable
            deallocate (AuxGrid3D)

            !if (Me%OutPut%ConcMaxTracer) deallocate (GridMaxTracer)

        !enddo d1



    end subroutine MonitorEulerian

    !--------------------------------------------------------------------------

    subroutine ModifyParticStatistic ()

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        real, dimension(:, :, :), pointer           :: AuxGrid3D 
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: iProp, STAT_CALL, em, p, ig

        !Begin--------------------------------------------------------------------------


        !Fills Grid concentration
        if (Me%Now > Me%ExternalVar%LastConcCompute) call FillGridConcentration 

d1:     do em =1, Me%EulerModelNumber

               
d2:     do ig = 1, Me%NGroups

            !Shorten
            ILB    = Me%EulerModel(em)%Size%ILB
            IUB    = Me%EulerModel(em)%Size%IUB
            JLB    = Me%EulerModel(em)%Size%JLB
            JUB    = Me%EulerModel(em)%Size%JUB
            KLB    = Me%EulerModel(em)%Size%KLB
            KUB    = Me%EulerModel(em)%Size%KUB

            !Allocates auxiliar variable
            allocate (AuxGrid3D (ILB:IUB, JLB:JUB, KLB:KUB))

            iProp = 1
            CurrentProperty => Me%FirstOrigin%FirstProperty
            do while (associated(CurrentProperty))

                if (Me%Statistic%ON(iProp)) then

                    do p=1, Me%Statistic%PropNumber
                    
                        if (CurrentProperty%ID == Me%Statistic%OptionsStat(p)%ID%IDNumber) exit

                    enddo

                    AuxGrid3D(:,:,:) = Me%EulerModel(em)%Lag2Euler%GridConc(:, :, :, iProp, ig)
                    call ModifyStatistic    (Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID(ig), &
                                             Value3D       = AuxGrid3D,                         &
                                             WaterPoints3D = Me%EulerModel(em)%Waterpoints3D,   &
                                             STAT          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                               &
                        stop 'ModifyParticStatistic - ModuleLagrangianGlobal - ERR01'

                    if (Me%OutPut%ConcMaxTracer) then

                        AuxGrid3D(:,:,:) = Me%EulerModel(em)%Lag2Euler%GridMaxTracer(:, :, :, iProp, ig)

                        call ModifyStatistic (Me%EulerModel(em)%PropStatistic(p)%Statistic2_ID(ig),                 &
                                              Value3D       = AuxGrid3D,                     &
                                              WaterPoints3D = Me%EulerModel(em)%Waterpoints3D,  &
                                              STAT          = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                           &
                            stop 'ModifyParticStatistic - ModuleLagrangianGlobal - ERR02'
                    endif 

                    if (Me%Statistic%OptionsStat(p)%Lag) then

                        call ComputeStatisticsLag(CurrentProperty%ID, p, ig, iProp)

                    endif

                endif
                CurrentProperty => CurrentProperty%Next
                iProp = iProp + 1
            enddo

            !Deallocates auxiliar variable
            deallocate (AuxGrid3D   )
            !if (Me%OutPut%ConcMaxTracer) deallocate (GridMaxTracer)

        enddo d2
        enddo d1


    end subroutine ModifyParticStatistic
    !--------------------------------------------------------------------------

    subroutine ComputeStatisticsLag(PropertyID, p, ig, iProp) 
        !Arguments-------------------------------------------------------------
        integer                                     :: PropertyID, p, ig, iProp

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
        integer                                     :: iClass, STAT_CALL, em
        real                                        :: DT, Aux, PropSpatial, Concentration

        !Begin--------------------------------------------------------------------------

d1:     do em =1, Me%EulerModelNumber

        ! In this version only writes and computes the statistic of tracers located in the higher priority eulerian model

!            em = Me%EulerModelNumber

            !Shorten
            ILB = Me%EulerModel(em)%WorkSize%ILB
            IUB = Me%EulerModel(em)%WorkSize%IUB
            JLB = Me%EulerModel(em)%WorkSize%JLB
            JUB = Me%EulerModel(em)%WorkSize%JUB
            KLB = Me%EulerModel(em)%WorkSize%KLB
            KUB = Me%EulerModel(em)%WorkSize%KUB

            call Search_Property(Me%OriginDefault, PropertyID, CurrentProperty)

            allocate (VolTracers  (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate (MassTracers (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate (PropTimeStep(ILB:IUB, JLB:JUB, KLB:KUB, 1 : Me%Statistic%OptionsStat(p)%nClassesLag))

            VolTracers  (:,:,:  ) = 0.
            MassTracers (:,:,:  ) = 0.
            PropTimeStep(:,:,:,:) = 0.

            !Gets global time step 
            call GetComputeTimeStep(Me%ExternalVar%ObjTime, DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStatisticsLag - ModuleLagrangianGlobal - ERR01'


            call GetStatisticClasses(Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID(ig),                          &
                                     Me%Statistic%OptionsStat(p)%ClassesLag, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ComputeStatisticsLag - ModuleLagrangianGlobal - ERR02'


            CurrentOrigin => Me%FirstOrigin

            do while (associated(CurrentOrigin))

                !Compute the total volume tracer by cell 
                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

                    i  = CurrentPartic%Position%I
                    j  = CurrentPartic%Position%J
                    k  = CurrentPartic%Position%K

                    if (em /= CurrentPartic%Position%ModelID) then
                        CurrentPartic => CurrentPartic%Next
                        cycle
                    endif
       
    WP1:            if (Me%EulerModel(em)%Waterpoints3D(i, j, k) == WaterPoint) then

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

                    if (em /= CurrentPartic%Position%ModelID) then
                        CurrentPartic => CurrentPartic%Next
                        cycle
                    endif
       
    WP2:            if (Me%EulerModel(em)%Waterpoints3D(i, j, k) == WaterPoint) then


    doClass:            do iClass = 1, Me%Statistic%OptionsStat(p)%nClassesLag
                            if (CurrentPartic%Concentration(iProp) >= Me%Statistic%OptionsStat(p)%ClassesLag(iClass, 1) .and. &
                                CurrentPartic%Concentration(iProp)  < Me%Statistic%OptionsStat(p)%ClassesLag(iClass, 2)) then

                                PropSpatial = CurrentPartic%Geometry%Volume / Me%EulerModel(em)%VolumeZ(i, j, k)


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

    WP3:        if (Me%EulerModel(em)%Waterpoints3D(i, j, k) == WaterPoint) then
                    
    VT:             if (VolTracers(i, j, k) < Me%EulerModel(em)%VolumeZ(i, j, k)) then

                        Aux = (Me%EulerModel(em)%VolumeZ(i, j, k) - VolTracers(i, j, k)) / Me%EulerModel(em)%VolumeZ(i, j, k)

                        DummyPos%I = i
                        DummyPos%J = j
                        DummyPos%K = k
               
                        call GetAmbientConcentration (CurrentProperty,                      &
                                                      em,                                   &
                                                      DummyPos,                             &
                                                      Concentration)


                    else VT

                        Concentration = MassTracers(i, j, k) / VolTracers(i, j, k)
                        PropTimeStep(i, j, k, :)      = 0
                        Aux                           = 1

                    endif VT

               
    Class1:         do iClass = 1, Me%Statistic%OptionsStat(p)%nClassesLag
                        if (Concentration >= Me%Statistic%OptionsStat(p)%ClassesLag(iClass, 1) .and.     &
                            Concentration <  Me%Statistic%OptionsStat(p)%ClassesLag(iClass, 2)) then
                            PropSpatial = Aux

                        else
                            PropSpatial = 0
                        endif

                        PropTimeStep(i, j, k, iClass) =  PropTimeStep(i, j, k, iClass) + PropSpatial
        
                        Me%EulerModel(em)%PropStatistic(p)%FrequencyLag     (i, j, k, iClass, ig) = &
                            (Me%EulerModel(em)%PropStatistic(p)%FrequencyLag(i, j, k, iClass, ig) * &
                            (Me%ExternalVar%RunPeriod - DT) +                                       &
                             PropTimeStep(i, j, k, iClass)  * DT         ) /                        &
                            (Me%ExternalVar%RunPeriod)

                    enddo Class1

                    if (sum(Me%EulerModel(em)%PropStatistic(p)%FrequencyLag(i, j, k, :, ig)) < 0.90) then

                        write(*,*) 'ComputeStatisticsLag - ModuleLagrangianGlobal - WARN10'
                        write(*,*) i, j, k, sum(Me%EulerModel(em)%PropStatistic(p)%FrequencyLag(i, j, k, :, ig)),'< 90%'

                    endif

                    if (sum(Me%EulerModel(em)%PropStatistic(p)%FrequencyLag(i, j, k, :, ig)) > 1.1) then

                        write(*,*) 'ComputeStatisticsLag - ModuleLagrangianGlobal - WARN20'
                        write(*,*) i, j, k, sum(Me%EulerModel(em)%PropStatistic(p)%FrequencyLag(i, j, k, :, ig)),'> 110%'

                    endif


                endif WP3

            enddo
            enddo
            enddo


            call UnGetStatistic(Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID(ig),                               &
                                Me%Statistic%OptionsStat(p)%ClassesLag, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeStatisticsLag - ModuleLagrangianGlobal - ERR05'


            deallocate (VolTracers )
            deallocate (MassTracers)
            deallocate (PropTimeStep)

        enddo d1
        

    end subroutine ComputeStatisticsLag

    !--------------------------------------------------------------------------


    subroutine Convert_XY_CellIJ (EulerModel,Position, Referential)

        !Arguments-------------------------------------------------------------
        type (T_EulerModel)                         :: EulerModel   
        type (T_Position)                           :: Position
        integer,    optional, intent(IN)            :: Referential

        !Local-----------------------------------------------------------------
        integer                                     :: Ipos, Jpos, STAT_CALL, GridBorderType
        real                                        :: PercI, PercJ


        GridBorderType = GetGridBorderType(EulerModel%ObjHorizontalGrid, Referential = AlongGrid_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Convert_XY_CellIJ - ModuleLagrangianGlobal - ERR10'


        if   (Referential == GridCoord_) then

            call GetXYCellZ(EulerModel%ObjHorizontalGrid, Position%CoordX,              &
                            Position%CoordY, Ipos, Jpos, PercI, PercJ,  Referential = GridCoord_, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'Convert_XY_CellIJ - ModuleLagrangianGlobal - ERR20'

            Position%CartX   = Position%CoordX
            Position%CartY   = Position%CoordY

            Position%X       = Position%CoordX
            Position%Y       = Position%CoordY

            if (EulerModel%Grid%GeoGrid) then

                call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Ipos, Jpos, PercI, PercJ,&
                                 Position%CartX, Position%CartY, Referential = Cartesian_, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'Convert_XY_CellIJ - ModuleLagrangianGlobal - ERR30'

                Position%X       = Position%CartX
                Position%Y       = Position%CartY

            endif

            if (GridBorderType /= Rectang_) then

                call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Ipos, Jpos, PercI, PercJ,&
                                 Position%X, Position%Y, Referential = AlongGrid_, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Convert_XY_CellIJ - ModuleLagrangianGlobal - ERR40'
            
            endif

        else if (Referential == AlongGrid_) then

            call GetXYCellZ(EulerModel%ObjHorizontalGrid, Position%X, Position%Y,       &
                            Ipos, Jpos, PercI, PercJ, Referential = AlongGrid_, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'Convert_XY_CellIJ - ModuleLagrangianGlobal - ERR50'

            Position%CartX   = Position%X
            Position%CartY   = Position%Y

            Position%CoordX  = Position%X
            Position%CoordY  = Position%Y

            if (EulerModel%Grid%GeoGrid) then

                call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Ipos, Jpos, PercI, PercJ, Position%CoordX, &
                                 Position%CoordY, Referential = GridCoord_, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'Convert_CellIJ_XY - ModuleLagrangianGlobal - ERR60'

            endif

            if (GridBorderType /= Rectang_) then
                call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Ipos, Jpos, PercI, PercJ, Position%CartX, &
                                 Position%CartY, Referential = Cartesian_, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'Convert_CellIJ_XY - ModuleLagrangianGlobal - ERR70'

                if (.not. EulerModel%Grid%GeoGrid) then

                    Position%CoordX = Position%CartX
                    Position%CoordY = Position%CartY
                    
                endif
            endif

       else if (Referential == Cartesian_) then


            call GetXYCellZ(EulerModel%ObjHorizontalGrid, Position%CartX, Position%CartY,       &
                            Ipos, Jpos, PercI, PercJ, Referential = Cartesian_, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'Convert_XY_CellIJ - ModuleLagrangianGlobal - ERR90'


            Position%X       = Position%CartX
            Position%Y       = Position%CartY

            Position%CoordX  = Position%CartX
            Position%CoordY  = Position%CartY

            if (EulerModel%Grid%GeoGrid) then

                call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Ipos, Jpos, PercI, PercJ, Position%CoordX, &
                                 Position%CoordY, Referential = GridCoord_, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'Convert_CellIJ_XY - ModuleLagrangianGlobal - ERR100'

            endif

            if (GridBorderType /= Rectang_) then

                call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Ipos, Jpos, PercI, PercJ, Position%X, &
                                 Position%Y, Referential = AlongGrid_, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'Convert_CellIJ_XY - ModuleLagrangianGlobal - ERR110'
 
            endif


        endif

        !CellJ Position
        Position%CellJ = (Jpos - 1) + PercJ                         
        Position%CellI = (Ipos - 1) + PercI


        Position%I = Ipos
        Position%J = Jpos


    end subroutine Convert_XY_CellIJ

    !--------------------------------------------------------------------------

    subroutine Convert_CellIJ_XY(EulerModel,Position)

        !Arguments-------------------------------------------------------------
        type (T_EulerModel)                         :: EulerModel   
        type (T_Position)                           :: Position


        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real                                        :: PercI, PercJ

        !----------------------------------------------------------------------

        !Cell Position
        Position%i           = int(Position%CellI) + 1
        Position%j           = int(Position%CellJ) + 1

        PercI       =  Position%CellI - int(Position%CellI)
        PercJ       =  Position%CellJ - int(Position%CellJ)

        call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Position%i, Position%j, PercI, PercJ, Position%CoordX, &
                         Position%CoordY, Referential = GridCoord_, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'Convert_CellIJ_XY - ModuleLagrangianGlobal - ERR10'


        call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Position%i, Position%j, PercI, PercJ, Position%X,      &
                         Position%Y, Referential = AlongGrid_, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'Convert_CellIJ_XY - ModuleLagrangianGlobal - ERR20'

        call GetCellZ_XY(EulerModel%ObjHorizontalGrid, Position%i, Position%j, PercI, PercJ, Position%CartX,  &
                         Position%CartY, Referential = Cartesian_, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'Convert_CellIJ_XY - ModuleLagrangianGlobal - ERR30'


    end subroutine Convert_CellIJ_XY

    !--------------------------------------------------------------------------


    subroutine Convert_Z_CellK (CurrentOrigin, EulerModel, Position, PositionCorrected)

        !Arguments-------------------------------------------------------------
   
        type (T_Origin),     pointer                :: CurrentOrigin
        type (T_EulerModel)                         :: EulerModel   
        type (T_Position)                           :: Position
        logical, optional                           :: PositionCorrected


        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, kFloor
        real                                        :: Perc
        logical                                     :: PositionCorrected_

        !----------------------------------------------------------------------

        i       = Position%I
        j       = Position%J
        k       = EulerModel%WorkSize%KUB
        kFloor  = EulerModel%kFloor(i, j)
        
        PositionCorrected_ = .false.

        !If the Particle is located above the surface, its placed closed to the
        !surface
        if (Position%Z < EulerModel%SZZ(i, j, k)) then
            Position%Z = EulerModel%SZZ(i, j, k)
            PositionCorrected_ = .true.
        endif

        !If the Particle is located below the bottom, its placed closed to the
        !bottom
        if (kFloor <= 0) then
            write(*,*) 'Kfloor =',KFloor,' i= ',i,' j= ',j
            stop 'LagrangianGlobal - Convert_Z_CellK - ERR10'
        endif
        if (Position%Z >= EulerModel%SZZ(i, j, kFloor-1)) then

            if (CurrentOrigin%State%Deposition) then

                !If the origin emit particles that can be deposited in the bottom and 
                !the particle tries to cross the bottom then is placed at distance from 
                !the bottom lower than the distance beyond each the particle is test 
                !if the shear stress is low enough to consider the particle deposited. 
                Position%Z = EulerModel%SZZ(I, J, KFloor-1) -  &
                             CurrentOrigin%Deposition%BottomDistance / 2.

            else

                Position%Z = EulerModel%SZZ(I, J, KFloor) +  &
                             EulerModel%DWZ(I, J, KFloor) * 0.9

            endif

            PositionCorrected_ = .true.

        endif

        if (present(PositionCorrected)) PositionCorrected = PositionCorrected_


        !Checks the layer
        do while(EulerModel%SZZ(i, j, k) <= Position%Z)
            k = k - 1
        end do

        Perc = (EulerModel%SZZ(i, j, k) - Position%Z) / (EulerModel%SZZ(i, j, k) - EulerModel%SZZ(i, j, k+1))

        !Avoid Rounding erros
        if (Perc >= 0.999) Perc = 0.99000

        Position%CellK = k + Perc

        !----------------------------------------------------------------------

    end subroutine Convert_Z_CellK

    !--------------------------------------------------------------------------

    subroutine Convert_CellK_Z (EulerModel,Position)

        !Arguments-------------------------------------------------------------
        type (T_EulerModel)                         :: EulerModel   
        type (T_Position)                           :: Position

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        real                                        :: Perc

        !----------------------------------------------------------------------

        i = Position%I
        j = Position%J
        k = int (Position%CellK) + 1

        Perc       = Position%CellK - int(Position%CellK)

        Position%Z = EulerModel%SZZ(i, j, k-1) - Perc * EulerModel%DWZ(i, j, k) 
        
    end subroutine Convert_CellK_Z

    !--------------------------------------------------------------------------

    subroutine Convert_CellK_K (Position)

        !Arguments-------------------------------------------------------------
        type (T_Position)                           :: Position

        Position%K = int(Position%CellK) + 1

    end subroutine Convert_CellK_K

    
    !--------------------------------------------------------------------------

    integer function Locate_ModelDomain (X,Y, NoDomain, Referential, emDif)

        !Arguments-------------------------------------------------------------
        real                                        :: X, Y
        logical, optional, intent(OUT)              :: NoDomain
        integer, optional, intent(IN)               :: Referential,  emDif      
        !Local-----------------------------------------------------------------
        integer                                     :: Referential_, em, STAT_CALL, emDif_
        logical                                     :: NoDomain_
        !----------------------------------------------------------------------

        if (present(Referential)) then

            Referential_ = Referential

        else

            Referential_ = GridCoord_

        endif

        if (present(emDif)) then

            emDif_ = emDif

        else

            emDif_ = FillValueInt

        endif

        NoDomain_ = .true.

d1:     do em =1, Me%EulerModelNumber

            if (GetXYInsideDomain(Me%EulerModel(em)%ObjHorizontalGrid, X, Y, Referential= Referential_, STAT = STAT_CALL) &
                .and. em /= emDif_) then
                if (STAT_CALL /= SUCCESS_) stop 'Locate_ModelDomain - ModuleLagrangianGlobal - ERR10'
                Locate_ModelDomain = em
                NoDomain_ = .false.
                exit
            endif


        enddo d1

        if (present(NoDomain)) then

            NoDomain = NoDomain_

        endif


 
        
    end function Locate_ModelDomain

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    integer function ReturnModelIndex (ModelName)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: ModelName

        !Local-----------------------------------------------------------------
        integer                                     :: em
        logical                                     :: NameFound = .false.

        !----------------------------------------------------------------------


d1:     do em =1, Me%EulerModelNumber

            if (trim(Me%EulerModel(em)%Name)==trim(ModelName)) then
                NameFound = .true.
                exit     
            endif

        enddo d1
        
        if (.not. NameFound) stop 'ReturnModelIndex - ModuleLagrangianGlobal - ERR10'

        ReturnModelIndex = em
        
    end function ReturnModelIndex

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

            call SetError(FATAL_, INTERNAL_, 'Search_Property - ModuleLagrangianGlobal - ERR01')

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
        integer                                     :: ParticleBeached
        integer                                     :: iProp, nProp, ig, nGroupProp
        type (T_Property), pointer                  :: CurrentProperty, FirstProperty
        type (T_Property), pointer                  :: AgeProperty
        integer                                     :: nPropAge, n
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        integer, dimension(:, :, :), pointer        :: OpenPoints3D
        real, dimension(:, :, :), pointer           :: SZZ
        real, dimension(:, :   ), pointer           :: OutPutMatrix
        integer                                     :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        character(StringLength)                     :: AuxChar
        real, dimension(:),       pointer           :: MaximumDepth
        integer                                     :: OutPutLines, JetTotalParticles, FirstParticle, em, em1, emMax, emp
        type (T_Position)                           :: Position
        real(8)                                     :: AverageX, AverageY, Stdv, RadiusOfInfluence
        
        !Begin--------------------------------------------------------------------------

        if (Me%Output%Write_) then

            OutPutNumber = Me%OutPut%NextOutPut
            !Shorten Var
            Actual       =  Me%Now

TOut:       if (Actual >= Me%OutPut%OutTime(OutPutNumber)) then 

                emMax = Me%EulerModelNumber 

                allocate(MaximumDepth(1:emMax))

                do em1 = 1, Me%EulerModelNumber 

                    call GetMaximumValue(Me%EulerModel(em1)%ObjGridData, MaximumDepth(em1), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR10'

                enddo

de:             do em = 1, Me%EulerModelNumber

                OpenPoints3D    => Me%EulerModel(em)%OpenPoints3D
                SZZ             => Me%EulerModel(em)%SZZ
                WorkILB         =  Me%EulerModel(em)%WorkSize%ILB 
                WorkIUB         =  Me%EulerModel(em)%WorkSize%IUB 
                        
                WorkJLB         =  Me%EulerModel(em)%WorkSize%JLB 
                WorkJUB         =  Me%EulerModel(em)%WorkSize%JUB 
                        
                WorkKLB         =  Me%EulerModel(em)%WorkSize%KLB 
                WorkKUB         =  Me%EulerModel(em)%WorkSize%KUB 

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


i0:             if (Me%RunOnline .and. em == emMax) then 

                    if (nP == 0) then
                        call DummyParticleStartDate(em)
                    else
                        call WriteRunOnline(em)
                    endif 

                    cycle

                endif i0

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
                call HDF5SetLimits  (Me%ObjHDF5(em), 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR20'

                call HDF5WriteData  (Me%ObjHDF5(em), "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR30'


                !Writes SZZ - HDF 5
                call HDF5SetLimits  (Me%ObjHDF5(em), WorkILB, WorkIUB, WorkJLB,   &
                                     WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR40'

                call HDF5WriteData  (Me%ObjHDF5(em), "/Grid/VerticalZ", "Vertical",  &
                                     "m", Array3D = SZZ, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR50'


                !Writes OpenPoints - HDF 5
                call HDF5SetLimits  (Me%ObjHDF5(em), WorkILB, WorkIUB,            &
                                     WorkJLB, WorkJUB, WorkKLB, WorkKUB,                 &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR60'

                call HDF5WriteData  (Me%ObjHDF5(em), "/Grid/OpenPoints", "OpenPoints",  &
                                     "-", Array3D = OpenPoints3D, OutputNumber = OutPutNumber, &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR70'

                if (Me%Now > Me%ExternalVar%LastConcCompute .and. em == 1) call FillGridConcentration

                call WriteGridConcentration(em) 

                !Flushes All pending HDF5 commands
                call HDF5FlushMemory (Me%ObjHDF5(em), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR75'

i1:             if (nP>0) then


                    if (Me%State%Oil) then
                        call WriteOilGridThickness     (em) 
                        call WriteOilGridConcentration (em)
                    endif

                    if (em /= emMax) cycle

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

                        call HDF5SetLimits  (Me%ObjHDF5(em), 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR80'

                        !XPosition
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            emp = CurrentPartic%Position%ModelID
                            Matrix1DX(nP)  = CurrentPartic%Position%CartX
                            Matrix1DY(nP)  = CurrentPartic%Position%CartY
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
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/X Pos", &
                                                "X Position",  "m", Array1D = Matrix1DX,                      &
                                                 Average = AverageX, Radius = RadiusOfInfluence,              &
                                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR90'
                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Y Pos", &
                                                "Y Position",  "m", Array1D = Matrix1DY,                      &
                                                 Average = AverageY, Radius = RadiusOfInfluence,              &
                                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR100'

                            if (Me%OutPut%OriginEnvelope) then
                                call WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY, &
                                                          "X Pos", "Y Pos",  "m", OutputNumber, em)  
                            endif

                        endif

                      !GeoGraphic Position

                        !Longitude, Latitude
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            Matrix1DX(nP)  =  GeographicCoordinates (CurrentPartic%Position%ModelID, CurrentPartic%Position, 1)
                            Matrix1DY(nP)  =  GeographicCoordinates (CurrentPartic%Position%ModelID, CurrentPartic%Position, 2)
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


                            call HDF5SetLimits  (Me%ObjHDF5(em), 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR80'


                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Longitude", &
                                                "Longitude",  "º", Array1D = Matrix1DX,                      &
                                                 Average = AverageX, Radius = RadiusOfInfluence,             &
                                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR110'

                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Latitude", &
                                                "Latitude",  "º", Array1D = Matrix1DY,                       &
                                                 Average = AverageY, Radius = RadiusOfInfluence,             &
                                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR120'

                            if (Me%OutPut%OriginEnvelope) then
                                call WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY, &
                                                          "Longitude", "Latitude", "º", OutputNumber, em)  
                            endif

                
                        endif
 
                        deallocate   (Matrix1DX)
                        deallocate   (Matrix1DY)


                        allocate    (Matrix1D(CurrentOrigin%nParticle))
                        Matrix1D(:) = FillValueReal


                        call HDF5SetLimits  (Me%ObjHDF5(em), 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR80'


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
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Z Pos", &
                                                "Z Position",  "m", Array1D = Matrix1D, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR130'
                        endif

                        !ZPosition for Vertical Cut
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            emp = CurrentPartic%Position%ModelID
                            Matrix1D(nP)  =  MaximumDepth(emp) - CurrentPartic%Position%Z
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
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Volume", &
                                                "Volume",  "m3", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR140'
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
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Origin ID", &
                                                "Origin ID",  "-", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR150'
                        endif

                        !Model ID
                        CurrentPartic   => CurrentOrigin%FirstPartic
                        nP = 0
                        do while (associated(CurrentPartic))
                            nP = nP + 1
                            Matrix1D(nP)  =  CurrentPartic%Position%ModelID
                            CurrentPartic => CurrentPartic%Next
                        enddo            
                        if (nP > 0) then
                            !HDF 5
                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Model ID", &
                                                "Model ID",  "-", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR150'
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
                                call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Age", &
                                                    "Age",  "days", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR160'
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
                                call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Beach Prop.", &
                                                    "Volume",  "m3", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR170'
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
                                call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Deposition State", &
                                                    "State",  "ON/OFF", Array1D = Matrix1D, OutputNumber = OutPutNumber,     &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR180'
                            endif
                    
                        endif


                        CurrentProperty => CurrentOrigin%FirstProperty
                        do while (associated(CurrentProperty))

                            if (CurrentProperty%ID == Fecal_Coliforms_ .or.CurrentProperty%ID == E_Coli_) then

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
                                    call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)// &
                                                         "/"//trim(CurrentProperty%T90Name),                     &
                                                         trim(CurrentProperty%T90Name),  "hours",                &
                                                         Array1D = Matrix1D, OutputNumber = OutPutNumber,        &
                                                         STAT = STAT_CALL)
                                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR190'
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
                                call HDF5WriteData  (Me%ObjHDF5(em),                  &
                                                     "/Results/"//trim(CurrentOrigin%Name)// &
                                                     "/"//trim(CurrentProperty%Name),        &
                                                     trim(CurrentProperty%Name),             &
                                                     trim(CurrentProperty%Units),            &
                                                     Array1D = Matrix1D,                     &
                                                     OutputNumber = OutPutNumber,            &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR200'
                            endif
                            nProp = nProp + 1
                            CurrentProperty => CurrentProperty%Next
                        enddo

                        deallocate  (Matrix1D)

    iplume:             if (CurrentOrigin%State%ComputePlume) then

                            emp = CurrentOrigin%Position%ModelID

                            call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                                                 OutPutLines = OutPutLines, STAT=STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR210'

                            if (OutPutLines < 1) then

                                call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)

                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR215'

                                call ActualizeJetProperties(CurrentOrigin)

                                call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                                                     OutPutLines = OutPutLines, STAT=STAT_CALL)

                                if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR220'

                                if (OutPutLines < 1) then
                                    stop 'ParticleOutput - ModuleLagrangianGlobal - ERR230'
                                endif

                            endif

                            JetTotalParticles = JetTotalParticles + OutPutLines

                            allocate   (Matrix1DX(OutPutLines))
                            allocate   (Matrix1DY(OutPutLines))

                            Matrix1DX(:) = FillValueReal
                            Matrix1DY(:) = FillValueReal

                            call HDF5SetLimits  (Me%ObjHDF5(em), 1, OutPutLines, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR240'

                            do np = 1, OutPutLines
                                Position%CartX = OutPutMatrix(np,2)
                                Position%CartY = OutPutMatrix(np,3)
                                call Convert_XY_CellIJ(Me%EulerModel(emp), Position, Referential = Cartesian_)
                                Matrix1DX(np) = GeographicCoordinates (emp, Position, 1)
                                Matrix1DY(np) = GeographicCoordinates (emp, Position, 2)
                            enddo

                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Plume X", &
                                                "Plume X",  "m", Array1D = Matrix1DX, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR250'

                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Plume Y", &
                                                "Plume Y",  "m", Array1D = Matrix1DY, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR260'


                            Matrix1DX(1:OutPutLines) = OutPutMatrix(1:OutPutLines,4)
                            Matrix1DY(1:OutPutLines) = OutPutMatrix(1:OutPutLines,6)


                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Plume Z", &
                                                "Plume Z",  "m", Array1D = Matrix1DX, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR270'

                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Plume Conc", &
                                                "Plume Conc",  "a", Array1D = Matrix1DY, OutputNumber = OutPutNumber,    &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR280'

                            call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR290'

                            deallocate   (Matrix1DX)
                            deallocate   (Matrix1DY)

                        endif iplume


                        !Integrates the sum of the particle which belong to the same 
                        !Group
                        do ig = 1, Me%nGroups
                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                TotParticle(ig) = TotParticle(ig) + CurrentOrigin%nParticle
                            endif
                        enddo


                        CurrentOrigin => CurrentOrigin%Next

                    enddo CurrOr


                    !Writes 1D data for every group
AllAtOnes:          if (Me%nOrigins > 1) then

dig:                do ig = 1, Me%nGroups
iTP:                if (TotParticle(ig) > 1) then
                        allocate    (Matrix1D(TotParticle(ig)))
                        Matrix1D(:) = FillValueReal

                        write (AuxChar, fmt='(i3)') Me%GroupIDs(ig)

                        call HDF5SetLimits  (Me%ObjHDF5(em), 1, TotParticle(ig),  &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR300'


                        !(XPosition)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    XPos:               do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    Matrix1D(nP)  =  CurrentPartic%Position%CartX
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo XPos

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
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

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    emp  = CurrentPartic%Position%ModelID
                                    Matrix1D(nP)  =  CurrentPartic%Position%CartY
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo YPos

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Y Pos",                       &
                                                   "Y Position",                             &
                                                   "m",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)


                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
                        do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    emp = CurrentPartic%Position%ModelID
                                    Matrix1D(nP)  =  GeographicCoordinates (emp, CurrentPartic%Position, 1)
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo
                            endif
                            CurrentOrigin => CurrentOrigin%Next
                        enddo

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
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

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                
                                    emp = CurrentPartic%Position%ModelID
                                    Matrix1D(nP)  =  GeographicCoordinates (emp, CurrentPartic%Position, 2)
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1

                                enddo
                            endif
                            CurrentOrigin => CurrentOrigin%Next
                        enddo

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Latitude",                    &
                                                   "Latitude",                               &
                                                   "º",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)


                        !(ZPosition)
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
    ZPos:               do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
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
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
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

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    emp = CurrentPartic%Position%ModelID
                                    Matrix1D(nP)  = MaximumDepth(emp) - CurrentPartic%Position%Z
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

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
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
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
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

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
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
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Origin ID",                   &
                                                   "Origin ID",                              &
                                                   "-",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)
                        !Model ID
                        nP = 1
                        CurrentOrigin => Me%FirstOrigin
ModelID:                do while (associated(CurrentOrigin))

                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                CurrentPartic   => CurrentOrigin%FirstPartic
                                do while (associated(CurrentPartic))
                                    Matrix1D(nP)  = CurrentPartic%Position%ModelID
                                    CurrentPartic => CurrentPartic%Next
                                    nP = nP + 1
                                enddo            
                            endif
                
                            CurrentOrigin => CurrentOrigin%Next
                        enddo ModelID


                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                    &
                                                   "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                   //"/Data_1D/Model ID",                    &
                                                   "Model ID",                               &
                                                   "-",                                      &
                                                   Array1D = Matrix1D,                       &
                                                   OutputNumber = OutPutNumber,              &
                                                   STAT = STAT_CALL)

                        !(Oil-Beached Particles)
iobp:                   if (Me%State%AssociateBeachProb) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
    OilState:               do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
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
                            call HDF5WriteData        (Me%ObjHDF5(em),                    &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/Beach Prop.",                 &
                                                       "Origin ID",                              &
                                                       "-",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)


                        end if iobp


                        !(Deposited Particles)
idp:                    if (Me%State%Deposition) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
    DepState:               do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then

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
                            call HDF5WriteData        (Me%ObjHDF5(em),                    &
                                                       "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                       //"/Data_1D/Deposition State",            &
                                                       "Origin ID",                              &
                                                       "-",                                      &
                                                       Array1D = Matrix1D,                       &
                                                       OutputNumber = OutPutNumber,              &
                                                       STAT = STAT_CALL)


                        end if idp

                        !T90
    iT90:               if (Me%State%T90Variable) then

                            nP = 1
                            CurrentOrigin => Me%FirstOrigin
    dT90:                   do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then

                                    CurrentProperty => CurrentOrigin%FirstProperty
                                    do while (associated(CurrentProperty))

                                        if (CurrentProperty%ID == Fecal_Coliforms_ .or.CurrentProperty%ID == E_Coli_) then

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
                            call HDF5WriteData        (Me%ObjHDF5(em),                               &
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
                                if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then

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
                            call HDF5WriteData        (Me%ObjHDF5(em),                    &
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
                            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
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

                                if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
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
                            call HDF5WriteData        (Me%ObjHDF5(em),                    &
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

                            call HDF5SetLimits  (Me%ObjHDF5(em), 1, JetTotalParticles, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR310'

                            FirstParticle = 1
                            CurrentOrigin => Me%FirstOrigin
    iplume3:                do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                
                                    if (CurrentOrigin%State%ComputePlume) then

                                        call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                                                             OutPutLines = OutPutLines, STAT=STAT_CALL)

                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR320'

                                        emp = CurrentOrigin%Position%ModelID

                                        do np = 1, OutPutLines
                                            Position%CartX = OutPutMatrix(np,2)
                                            Position%CartY = OutPutMatrix(np,3)
                                            call Convert_XY_CellIJ(Me%EulerModel(emp), Position, Referential = Cartesian_)
                                            Matrix1DX(FirstParticle-1+np) = GeographicCoordinates (emp, Position, 1)
                                            Matrix1DY(FirstParticle-1+np) = GeographicCoordinates (emp, Position, 2)
                                        enddo

                                        call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)
                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR330'

                                        FirstParticle = FirstParticle + OutPutLines
                                    endif
                            
                                endif

                                CurrentOrigin => CurrentOrigin%Next

                            enddo iplume3
                            

                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume X", "Plume X",  "m", Array1D = Matrix1DX,  &
                                                OutputNumber = OutPutNumber, STAT = STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR340'

                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume Y", "Plume Y",  "m", Array1D = Matrix1DY,  &
                                                OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR350'

                            FirstParticle = 1
                            CurrentOrigin => Me%FirstOrigin
    iplume4:                do while (associated(CurrentOrigin))
                                if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then
                                
                                    if (CurrentOrigin%State%ComputePlume) then

                                        call GetOutPutMatrix(CurrentOrigin%Movement%ObjJet, OutPutMatrix, &
                                                             OutPutLines = OutPutLines, STAT=STAT_CALL)

                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR360'

                                        Matrix1DX(FirstParticle:FirstParticle-1+OutPutLines) = OutPutMatrix(1:OutPutLines,4)
                                        Matrix1DY(FirstParticle:FirstParticle-1+OutPutLines) = OutPutMatrix(1:OutPutLines,6)

                                        call UnGetJet(CurrentOrigin%Movement%ObjJet, OutPutMatrix, STAT=STAT_CALL)
                                        if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR370'

                                        FirstParticle = FirstParticle + OutPutLines
                                    endif
                            
                                endif

                                CurrentOrigin => CurrentOrigin%Next

                            enddo iplume4


                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume Z", "Plume Z",  "m", Array1D = Matrix1DX,  &
                                                OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR380'

                            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/Group_"//trim(adjustl(AuxChar)) &
                                                //"/Data_1D/"//"/Plume Conc", "Plume Conc",  "a",                 &
                                                Array1D = Matrix1DY, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR390'


                            deallocate   (Matrix1DX)
                            deallocate   (Matrix1DY)

                        endif iplume2


                        deallocate  (Matrix1D)

                    endif iTP
                    enddo dig
                    endif AllAtOnes

                
                    if (Me%State%Monitor) then
                        call WriteMonitorOutput ()
                    endif


                    !Flushes All pending HDF5 commands
                    call HDF5FlushMemory (Me%ObjHDF5(em), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ParticleOutput - ModuleLagrangianGlobal - ERR400'

                endif i1

                enddo de


                !Increments Output number
                Me%OutPut%NextOutPut = OutPutNumber + 1

                !Verifies if all outputs are done (necessary for partic DT smaller then global DT)
                if (Me%OutPut%NextOutPut > size(Me%OutPut%OutTime)) then
                    Me%Output%Write_ = .false.
                endif

                deallocate(MaximumDepth)

            endif  TOut
        endif


    end subroutine ParticleOutput

    !--------------------------------------------------------------------------

!--------------------------------------------------------------------------


    subroutine DummyParticleStartDate(em)

        !Arguments-------------------------------------------------------------
        integer                                     :: em
        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        real(8), dimension(:), pointer              :: Matrix1D, Matrix1DX, Matrix1DY
        integer                                     :: STAT_CALL
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        real(8)                                     :: AverageX, AverageY, RadiusOfInfluence
        
        !Begin--------------------------------------------------------------------------


        !Writes the Instant - HDF 5
        call ExtractDate   (Me%Now, AuxTime(1), AuxTime(2), AuxTime(3),                 &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime
        call HDF5SetLimits  (Me%ObjHDF5(em), 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR10'

        call HDF5WriteData  (Me%ObjHDF5(em), "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                             Array1D = TimePtr, OutputNumber = Me%OutPut%NextOutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR20'


        !Writes Data for every origin
        CurrentOrigin => Me%FirstOrigin
CurrOr:     do while (associated(CurrentOrigin))

            allocate   (Matrix1DX(1))
            allocate   (Matrix1DY(1))

            call HDF5SetLimits  (Me%ObjHDF5(em), 1, 1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR70'


          !GeoGraphic Position

            !Longitude, Latitude
            Matrix1DX(1)  =  CurrentOrigin%Position%CoordX
            Matrix1DY(1)  =  CurrentOrigin%Position%CoordY

            if (CurrentOrigin%AveragePositionON) then
                AverageX = CurrentOrigin%Position%CoordX
                AverageY = CurrentOrigin%Position%CoordY
                RadiusOfInfluence = 0.
            else
                AverageX          = FillValueReal
                AverageY          = FillValueReal
                RadiusOfInfluence = FillValueReal
            endif 


            call HDF5SetLimits  (Me%ObjHDF5(em), 1, 1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR80'


            !HDF 5
            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Longitude", &
                                "Longitude",  "º", Array1D = Matrix1DX,                 &
                                 Average = AverageX, Radius = RadiusOfInfluence,        &
                                 OutputNumber = Me%OutPut%NextOutPut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR90'

            !HDF 5
            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Latitude", &
                                "Latitude",  "º", Array1D = Matrix1DY,                  &
                                 Average = AverageY, Radius = RadiusOfInfluence,        &
                                 OutputNumber = Me%OutPut%NextOutPut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR100'

            if (Me%OutPut%OriginEnvelope) then
                call WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY,           &
                                          "Longitude", "Latitude", "º", Me%OutPut%NextOutPut, em)  
            endif



            deallocate   (Matrix1DX)
            deallocate   (Matrix1DY)

            allocate    (Matrix1D(1))
            Matrix1D(:) = 0.

            call HDF5SetLimits  (Me%ObjHDF5(em), 1, 1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR110'


         !HDF 5
            call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Z Pos", &
                                "Z Position",  "m", Array1D = Matrix1D, OutputNumber = Me%OutPut%NextOutPut,    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR120'

            deallocate  (Matrix1D)


            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5(em), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DummyParticleStartDate - ModuleLagrangianGlobal - ERR130'



    end subroutine DummyParticleStartDate

    !--------------------------------------------------------------------------

!--------------------------------------------------------------------------


    subroutine WriteRunOnline(em)

        !Arguments-------------------------------------------------------------
        integer                                     :: em
        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        type (T_Partic), pointer                    :: CurrentPartic
        real(8), dimension(:), pointer              :: Matrix1D, Matrix1DX, Matrix1DY
        integer                                     :: STAT_CALL
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        real(8)                                     :: AverageX, AverageY, Stdv, RadiusOfInfluence
        integer                                     :: nP, n

        
        !Begin--------------------------------------------------------------------------


        !Writes the Instant - HDF 5
        call ExtractDate   (Me%Now, AuxTime(1), AuxTime(2), AuxTime(3),                 &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime
        call HDF5SetLimits  (Me%ObjHDF5(em), 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR10'

        call HDF5WriteData  (Me%ObjHDF5(em), "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                             Array1D = TimePtr, OutputNumber = Me%OutPut%NextOutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR20'


        !Writes Data for every origin
        CurrentOrigin => Me%FirstOrigin
CurrOr:     do while (associated(CurrentOrigin))

                allocate   (Matrix1DX(CurrentOrigin%nParticle))
                allocate   (Matrix1DY(CurrentOrigin%nParticle))
                Matrix1DX(:) = FillValueReal
                Matrix1DY(:) = FillValueReal

                call HDF5SetLimits  (Me%ObjHDF5(em), 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangian - ERR30'


              !GeoGraphic Position

                !Longitude, Latitude
                CurrentPartic   => CurrentOrigin%FirstPartic
                nP = 0
                do while (associated(CurrentPartic))
                    nP = nP + 1
                    Matrix1DX(nP)  =  GeographicCoordinates (CurrentPartic%Position%ModelID, CurrentPartic%Position, 1)
                    Matrix1DY(nP)  =  GeographicCoordinates (CurrentPartic%Position%ModelID, CurrentPartic%Position, 2)
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


                    call HDF5SetLimits  (Me%ObjHDF5(em), 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR40'


                    !HDF 5
                    call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Longitude", &
                                        "Longitude",  "º", Array1D = Matrix1DX,                      &
                                         Average = AverageX, Radius = RadiusOfInfluence,             &
                                         OutputNumber = Me%OutPut%NextOutPut, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR50'

                    !HDF 5
                    call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Latitude", &
                                        "Latitude",  "º", Array1D = Matrix1DY,                       &
                                         Average = AverageY, Radius = RadiusOfInfluence,             &
                                         OutputNumber = Me%OutPut%NextOutPut, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR60'

                    if (Me%OutPut%OriginEnvelope) then
                        call WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY, &
                                                  "Longitude", "Latitude", "º", Me%OutPut%NextOutPut, em)  
                    endif

        
                endif

                deallocate   (Matrix1DX)
                deallocate   (Matrix1DY)


                allocate    (Matrix1D(CurrentOrigin%nParticle))
                Matrix1D(:) = FillValueReal


                call HDF5SetLimits  (Me%ObjHDF5(em), 1, CurrentOrigin%nParticle, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR70'


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
                    call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/Z Pos", &
                                        "Z Position",  "m", Array1D = Matrix1D, OutputNumber = Me%OutPut%NextOutPut,    &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR80'
                        endif

            deallocate  (Matrix1D)


            CurrentOrigin => CurrentOrigin%Next

        enddo CurrOr

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5(em), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteRunOnline - ModuleLagrangianGlobal - ERR90'



    end subroutine WriteRunOnline

    !--------------------------------------------------------------------------

    subroutine WriteOriginEnvelope(CurrentOrigin, Matrix1DX, Matrix1DY,                &
                                    StringX, StringY, Units, OutputNumber, em)

        !Arguments------------------------------------------------------------
        type (T_Origin), pointer            :: CurrentOrigin
        real(8),   dimension(:), pointer    :: Matrix1DX, Matrix1DY
        character(Len=*)                    :: StringX, StringY, Units
        integer                             :: OutputNumber, em

        !Local-----------------------------------------------------------------
        real,      dimension(:), pointer    :: NodeX, NodeY
        real(8),   dimension(:), pointer    :: Envelope1DX, Envelope1DY
        integer,   dimension(:), pointer    :: BoundaryNodes
        character(Len=StringLength)         :: StringXaux, StringYaux
        real                                :: AverageX, AverageY
        integer                             :: ObjTriangulation, NumberOfBoundaryNodes
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
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangianGlobal - ERR10'

            call GetNumberOfBoundaryNodes (ObjTriangulation, NumberOfBoundaryNodes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangianGlobal - ERR20'

            allocate(Envelope1DX  (1:NumberOfBoundaryNodes),                                &
                     Envelope1DY  (1:NumberOfBoundaryNodes),                                &
                     BoundaryNodes(1:NumberOfBoundaryNodes))


            call GetBoundaryNodes (ObjTriangulation, BoundaryNodes, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangianGlobal - ERR30'
            !call GetOuterLine

            do i=1, NumberOfBoundaryNodes
                j = BoundaryNodes(i)
                Envelope1DX(i) = NodeX(j) + AverageX
                Envelope1DY(i) = NodeY(j) + AverageY
            enddo

            call KillTriangulation(ObjTriangulation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangianGlobal - ERR70'

            deallocate(NodeX, NodeY)
            deallocate(BoundaryNodes)

        else if (NumberOfNodes == 1) then
            NumberOfBoundaryNodes = 3

            allocate(Envelope1DX  (1:NumberOfBoundaryNodes),                                &
                     Envelope1DY  (1:NumberOfBoundaryNodes))
            
            Envelope1DX  (1:1) = Matrix1DX(1) - 1.e-5
            Envelope1DY  (1:1) = Matrix1DY(1) - 1.e-5
            Envelope1DX  (2:2) = Matrix1DX(1)
            Envelope1DY  (2:2) = Matrix1DY(1)
            Envelope1DX  (3:3) = Matrix1DX(1) + 1.e-5
            Envelope1DY  (3:3) = Matrix1DY(1) - 1.e-5

        endif

        call HDF5SetLimits  (Me%ObjHDF5(em), 1, NumberOfBoundaryNodes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangianGlobal - ERR40'


        !HDF 5
        call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/"//trim(StringXaux), &
                            trim(StringXaux),  trim(Units), Array1D = Envelope1DX,                     &
                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangianGlobal - ERR50'

        !HDF 5
        call HDF5WriteData  (Me%ObjHDF5(em), "/Results/"//trim(CurrentOrigin%Name)//"/"//trim(StringYaux), &
                            trim(StringYaux),  trim(Units), Array1D = Envelope1DY,                     &
                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteOriginEnvelope - ModuleLagrangianGlobal - ERR60'


        deallocate(Envelope1DX, Envelope1DY)


    end subroutine WriteOriginEnvelope

    !--------------------------------------------------------------------------

    subroutine OutputRestartFile
        
        !Local-----------------------------------------------------------------
        real                                :: Year, Month, Day, Hour, Minute, Second
        logical                             :: WriteFinal

        !----------------------------------------------------------------------

        if(Me%NextCompute >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then


            call ExtractDate(Me%Now,                         &
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
        real, dimension(:, :, :), pointer               :: GridConc3D
        integer                                         :: STAT_CALL
        type (T_Property), pointer                      :: CurrentProperty
        integer                                         :: iProp, em, ig
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: TimeSerieNumber, dn, id, jd, kd
        logical                                         :: DepthON, IgnoreOK
        real                                            :: DepthLevel

        !Begin-----------------------------------------------------------------
        !Fills Grid concentration
        if (Me%Now > Me%ExternalVar%LastConcCompute) call FillGridConcentration 

d1:     do em =1, Me%EulerModelNumber

        ! In this version only writes times for the first origin group
!d2:    do ig =1, Me%NGroups

        !em = Me%EulerModelNumber
        ig = 1

        !Shorten
        ILB    = Me%EulerModel(em)%Size%ILB
        IUB    = Me%EulerModel(em)%Size%IUB
        JLB    = Me%EulerModel(em)%Size%JLB
        JUB    = Me%EulerModel(em)%Size%JUB
        KLB    = Me%EulerModel(em)%Size%KLB
        KUB    = Me%EulerModel(em)%Size%KUB

        !Allocates auxiliar variable
        allocate (GridConc3D (ILB:IUB, JLB:JUB, KLB:KUB         ))

        !Corrects if necessary the cell of the time serie based in the time serie depth
        call GetNumberOfTimeSeries(Me%EulerModel(em)%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR10'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%EulerModel(em)%ObjTimeSerie, dn,               &  
                                      LocalizationI = id,                               &
                                      LocalizationJ = jd,                               &
                                      DepthLevel    = DepthLevel,                       &
                                      DepthON       = DepthON,                          & 
                                      STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR20'

            if (DepthON) then

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%EulerModel(em)%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR30'

                    if (IgnoreOK) then
                        cycle
                    else
                        stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR40'
                    endif

                endif

                kd = GetLayer4Level(Me%EulerModel(em)%ObjGeometry, id, jd, DepthLevel, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR50'

                call CorrectsCellsTimeSerie(Me%EulerModel(em)%ObjTimeSerie, dn,  k = kd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR60'

                
                if (Me%EulerModel(em)%WaterPoints3D(id, jd, kd) /= WaterPoint .and. Me%FirstIteration) then

                    write(*,*) 'Time serie station I=',Id, 'J=',Jd,'K=',kd,'is located in land' 
                    write(*,*) 'OutPut_TimeSeries - ModuleWaterProperties - WRN100'

                endif
            endif


        enddo

        iProp = 0
        CurrentProperty => Me%OriginDefault%FirstProperty
        do while (associated(CurrentProperty))

            if (CurrentProperty%WritesTimeSerie) then
                iProp = iProp + 1
            
                GridConc3D = Me%EulerModel(em)%Lag2Euler%GridConc(:, :, :, iProp, ig)

                call WriteTimeSerie(Me%EulerModel(em)%ObjTimeSerie,                         &
                                    Data3D = GridConc3D,                                    &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR70'

                if (CurrentProperty%T90ON) then

                    call ComputeT90Matrix(em, ig, CurrentProperty, GridConc3D)

                    call WriteTimeSerie(Me%EulerModel(em)%ObjTimeSerie,                     &
                                        Data3D = GridConc3D,                                &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleLagrangianGlobal - ERR80'

                endif
            endif

            CurrentProperty => CurrentProperty%Next
        enddo

        !Deallocates Temporary Matrixes
        deallocate (GridConc3D)


!        enddo d2
        enddo d1

    
    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------

    subroutine FillGridConcentration

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                        :: CurrentOrigin
        type (T_Partic), pointer                        :: CurrentPartic
        type (T_Property), pointer                      :: FirstProperty, CurrentProperty
        real,    dimension(:, :      ), pointer         :: DUX, DVY
        integer                                         :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                         :: WS_KLB, WS_KUB
        integer                                         :: i, j, k, nProp, status
        integer                                         :: iV, jV, kV, Delta, iAP
        integer                                         :: iInf, iSup
        integer                                         :: jInf, jSup
        integer                                         :: kInf, kSup
        real(8)                                         :: VolumeTotal, Coef
        real                                            :: DiffVolCel
        type (T_Position)                               :: DummYPos    
        logical                                         :: FoundSediment, PartInside
        integer                                         :: Sediment_ID, np, em, ig, STAT_CALL

        !Begin----------------------------------------------------------------------


        Me%ExternalVar%LastConcCompute = Me%Now

d1:     do em = 1, Me%EulerModelNumber 


            WS_ILB = Me%EulerModel(em)%WorkSize%ILB
            WS_IUB = Me%EulerModel(em)%WorkSize%IUB
            WS_JLB = Me%EulerModel(em)%WorkSize%JLB
            WS_JUB = Me%EulerModel(em)%WorkSize%JUB
            WS_KLB = Me%EulerModel(em)%WorkSize%KLB
            WS_KUB = Me%EulerModel(em)%WorkSize%KUB


            nProp           =  Me%OriginDefault%nProperties 
                                                        !i,j,k,p,ig 
            Me%EulerModel(em)%Lag2Euler%GridVolume      (:,:,:,  :) = 0.
            Me%EulerModel(em)%Lag2Euler%GridTracerNumber(:,:,:,  :) = 0.
            Me%EulerModel(em)%Lag2Euler%GridMass        (:,:,:,:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%GridConc        (:,:,:,:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%MeanConc              (:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%AmbientConc           (:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%MinConc               (:,:) = 0.
            Me%EulerModel(em)%Lag2Euler%MassVolCel            (:,:) = 0.


            if (Me%State%Deposition) then
                                                          !i,j,k,p,ig 
                Me%EulerModel(em)%Lag2Euler%GridBottomMass(:,:,:,  :) = 0.
                Me%EulerModel(em)%Lag2Euler%GridBottomConc(:,:,:,  :) = 0.
            endif

            if (Me%OutPut%ConcMaxTracer) then
                Me%EulerModel(em)%Lag2Euler%GridMaxTracer (:,:,:,:,:) = 0.
            endif

dg:     do ig = 1, Me%NGroups 

            !Integrates the Volume and the Mass in each GridCell
            CurrentOrigin => Me%FirstOrigin
    CurrOr: do while (associated(CurrentOrigin))

            if (CurrentOrigin%GroupID == Me%GroupIDs(ig)) then

                CurrentPartic => CurrentOrigin%FirstPartic
                do while (associated(CurrentPartic))

       
                    if (em == CurrentPartic%Position%ModelID) then

                        i = CurrentPartic%Position%I
                        j = CurrentPartic%Position%J
                    
                    else

                        PartInside = GetXYInsideDomain(Me%EulerModel(em)%ObjHorizontalGrid, CurrentPartic%Position%CoordX, &
                                                       CurrentPartic%Position%CoordY,Referential = GridCoord_, STAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) stop 'FillGridConcentration - ModuleLagrangianGlobal - ERR05'

                        if (PartInside) then

                            call GetXYCellZ(Me%EulerModel(em)%ObjHorizontalGrid, CurrentPartic%Position%CoordX, &
                                            CurrentPartic%Position%CoordY, i, j, Referential = GridCoord_, STAT = STAT_CALL)

                            if (STAT_CALL /= SUCCESS_) stop 'FillGridConcentration - ModuleLagrangianGlobal - ERR10'

                        else
                            CurrentPartic => CurrentPartic%Next
                            cycle
                        endif
                    
                    endif



    cd1:            if (.not. CurrentPartic%Deposited) then


                        if (em == CurrentPartic%Position%ModelID) then
                            k = CurrentPartic%Position%K
                        else
                            k = GetLayer4Level(Me%EulerModel(em)%ObjGeometry, i, j, CurrentPartic%Position%Z, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'FillGridConcentration - ModuleLagrangianGlobal - ERR20'
                        endif


                        Me%EulerModel(em)%Lag2Euler%GridTracerNumber(i, j, k, ig) = &
                            Me%EulerModel(em)%Lag2Euler%GridTracerNumber(i, j, k, ig) + 1

                        if (Me%OutPut%ConcMaxTracer) then

                            !In this grid is stored the maximum concentration of all tracers present in the cell
                            Me%EulerModel(em)%Lag2Euler%GridMaxTracer(i, j, k, :, ig) = &
                                max(Me%EulerModel(em)%Lag2Euler%GridMaxTracer(i, j, k, :, ig), &
                                CurrentPartic%Concentration(:))

                        endif

                        !Particle fits inside Grid Cell?    
    cd2:                if (CurrentPartic%Geometry%Volume <= Me%EulerModel(em)%VolumeZ(i, j, k)) then

                            Me%EulerModel(em)%Lag2Euler%GridMass  (i, j, k, :, ig) = &
                                Me%EulerModel(em)%Lag2Euler%GridMass  (i, j, k, :, ig) + CurrentPartic%Mass(:)
                            Me%EulerModel(em)%Lag2Euler%GridVolume (i, j, k, ig)   = &
                                Me%EulerModel(em)%Lag2Euler%GridVolume(i, j, k, ig)    + CurrentPartic%Geometry%Volume

                        else  cd2

                            WS_ILB = Me%EulerModel(em)%WorkSize%ILB
                            WS_IUB = Me%EulerModel(em)%WorkSize%IUB
                            WS_JLB = Me%EulerModel(em)%WorkSize%JLB
                            WS_JUB = Me%EulerModel(em)%WorkSize%JUB
                            WS_KLB = Me%EulerModel(em)%WorkSize%KLB
                            WS_KUB = Me%EulerModel(em)%WorkSize%KUB

                            VolumeTotal = Me%EulerModel(em)%VolumeZ(i, j, k)
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
                                    stop       'FillGridConcentration - ModuleLagrangianGlobal - ERR30'

                                endif

                                VolumeTotal = 0.
                                do iV = iInf, iSup
                                do jV = jInf, jSup
                                do kV = kinf, kSup
                                    if (Me%EulerModel(em)%OpenPoints3D(iV, jV, kV) == OpenPoint) then
                                        VolumeTotal = VolumeTotal +                              &
                                                      Me%EulerModel(em)%VolumeZ(iV, jV, kV)
                                    endif
                                enddo
                                enddo
                                enddo

                            enddo
            
                            do iV = iInf, iSup
                            do jV = jInf, jSup
                            do kV = kinf, kSup
                                if (Me%EulerModel(em)%OpenPoints3D(iV, jV, kV) == OpenPoint) then
                    
                                    Coef = Me%EulerModel(em)%VolumeZ(iV, jV, kV) / VolumeTotal
                    
                                    Me%EulerModel(em)%Lag2Euler%GridMass  (iV, jV, kV, :, ig) = &
                                        Me%EulerModel(em)%Lag2Euler%GridMass  (iV, jV, kV, :, ig) + &
                                        CurrentPartic%Mass(:) * Coef

                                    Me%EulerModel(em)%Lag2Euler%GridVolume(iV, jV, kV, ig) = &
                                        Me%EulerModel(em)%Lag2Euler%GridVolume(iV, jV, kV, ig)    + &
                                        CurrentPartic%Geometry%Volume * Coef

                                endif
                            enddo
                            enddo
                            enddo
            
                        endif cd2

                    else  cd1 ! The particle is deposited in the bottom
                    !In this case no test is made to verify if the the particle occupies more then one cell
                    !to maintain the algothim simple.

                  
                        Me%EulerModel(em)%Lag2Euler%GridBottomMass (i, j, :, ig) =  &
                            Me%EulerModel(em)%Lag2Euler%GridBottomMass (i, j, :, ig) + CurrentPartic%Mass(:)

                    endif cd1

                    CurrentPartic => CurrentPartic%Next
                enddo

                if (.not. associated(FirstProperty)) FirstProperty => CurrentOrigin%FirstProperty

            endif

            CurrentOrigin => CurrentOrigin%Next
            enddo CurrOr


            if (Me%State%Deposition) then !Find if exist the 'sediment' property

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

        enddo dg

        enddo d1

d2:     do em = 1, Me%EulerModelNumber 

            WS_ILB = Me%EulerModel(em)%WorkSize%ILB
            WS_IUB = Me%EulerModel(em)%WorkSize%IUB
            WS_JLB = Me%EulerModel(em)%WorkSize%JLB
            WS_JUB = Me%EulerModel(em)%WorkSize%JUB
            WS_KLB = Me%EulerModel(em)%WorkSize%KLB
            WS_KUB = Me%EulerModel(em)%WorkSize%KUB

            !Fills up Grid Concentration
g2:         do ig = 1, Me%NGroups
             
            do k = WS_KLB, WS_KUB
            do j = WS_JLB, WS_JUB
            do i = WS_ILB, WS_IUB

                if (Me%EulerModel(em)%Waterpoints3D (i, j, k) == WaterPoint) then

                    !Mean Concentration of the particle
                    Me%EulerModel(em)%Lag2Euler%MeanConc(:,ig) = 0.0
                    if (Me%EulerModel(em)%Lag2Euler%GridVolume(i, j, k, ig) > 0.0) then
                        Me%EulerModel(em)%Lag2Euler%MeanConc  (:, ig) = &
                            Me%EulerModel(em)%Lag2Euler%GridMass  (i, j, k, :, ig) / &
                            Me%EulerModel(em)%Lag2Euler%GridVolume(i, j, k, ig)
                    endif

                    !Ambient Concentration of the place of the particle
                    iAP = 1
!                    CurrentProperty => FirstProperty

                    CurrentProperty => Me%OriginDefault%FirstProperty

                    do while (associated(CurrentProperty))
            
                        DummYPos%I = i
                        DummYPos%J = j
                        DummYPos%K = k
                
                        call GetAmbientConcentration (CurrentProperty,                  &
                                                      em,                               &
                                                      DummYPos,                         &
                                                      Me%EulerModel(em)%Lag2Euler%AmbientConc(iAP, ig))

                        Me%EulerModel(em)%Lag2Euler%MinConc (iAP,ig) = CurrentProperty%Min_concentration

                        iAP = iAP + 1
                        CurrentProperty => CurrentProperty%Next
                    enddo


                    select case (Me%OutPut%OutPutConcType)

                    case (Maximum)


                        !Metodo Max : Admite-se a concentracao igual ao valor maximo de entre 2 valores: 
                        !          - massa total de tracadores presente na celula a dividir pelo volume da celula
                        !          - concentracao media dos tracadores
                               
                        Me%EulerModel(em)%Lag2Euler%MassVolCel(:, ig) = &
                            Me%EulerModel(em)%Lag2Euler%GridMass(i, j, k, :, ig) / &
                            Me%EulerModel(em)%VolumeZ(i, j, k)
                        
                        Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig) = &
                            max(Me%EulerModel(em)%Lag2Euler%MassVolCel(:, ig), &
                            Me%EulerModel(em)%Lag2Euler%MeanConc(:, ig)) 

                    case (Mean)

                        !Metodo 2 : Se o volume total de tracadores for menor
                        !      que o volume da célula i,j,k, entao a concentracao nesta  e igual a uma media 
                        !      entre a concentracao media dos tracadores e a concentracao ambiente
                        !      caso contrario a concentracao na celula i,j,k e igual a concentracao media 
                        !      dos tracadores

                        if (Me%EulerModel(em)%VolumeZ(i, j, k) >=                       &
                            Me%EulerModel(em)%Lag2Euler%GridVolume(i, j, k, ig)) then 
                            DiffVolCel = Me%EulerModel(em)%VolumeZ(i, j, k) -           &
                            Me%EulerModel(em)%Lag2Euler%GridVolume(i, j, k, ig)  
            
                            Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig) =  &
                            (DiffVolCel * Me%EulerModel(em)%Lag2Euler%AmbientConc(:,ig) + &
                                          Me%EulerModel(em)%Lag2Euler%GridMass  (i, j, k, :,ig)) / &  
                                          Me%EulerModel(em)%VolumeZ(i, j, k)

                        else

                            Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig) =  &
                                Me%EulerModel(em)%Lag2Euler%MeanConc(:, ig)

                        endif

                    end select

                    where (Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig).lt.     &
                           Me%EulerModel(em)%Lag2Euler%MinConc(:, ig))                  &
                           Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig) =       &
                           Me%EulerModel(em)%Lag2Euler%AmbientConc(:,ig)

                    where (Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig)  == 0.)                                   &
                           Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig)  =      &
                           Me%EulerModel(em)%Lag2Euler%AmbientConc(:,ig)


                else

                     Me%EulerModel(em)%Lag2Euler%GridConc(i, j, k, :, ig) = null_real

                end if

            end do
            end do
            end do

            end do g2

    cd3:    if (Me%State%Deposition) then ! fills up the bottom concentration in a simplified way

                !Gets Horizontal Grid
                call GetHorizontalGrid(Me%EulerModel(em)%ObjHorizontalGrid, DUX = DUX, DVY = DVY, STAT = status)
                if (status /= SUCCESS_) call SetError(FATAL_, INTERNAL_, 'FillGridConcentration - ModuleLagrangianGlobal - ERR40')

g3:             do ig = 1, Me%NGroups

                    do j = WS_JLB, WS_JUB
                    do i = WS_ILB, WS_IUB

                        if (Me%EulerModel(em)%Waterpoints3D (i, j, WS_KUB) == WaterPoint) then

                            if (FoundSediment) then
                                do np = 1, SIZE(Me%EulerModel(em)%Lag2Euler%GridBottomConc, DIM = 3)

                                    if (np /= Sediment_ID) then 
                                        !Mass contaminant / Mass sediment
                                        if (Me%EulerModel(em)%Lag2Euler%GridBottomMass(i, j, Sediment_ID, ig) > 1e-12) then
                                            Me%EulerModel(em)%Lag2Euler%GridBottomConc(i, j, np, ig) = &
                                            Me%EulerModel(em)%Lag2Euler%GridBottomMass(i, j, np, ig) / &
                                            Me%EulerModel(em)%Lag2Euler%GridBottomMass(i, j, Sediment_ID, ig)
                                        else
                                            Me%EulerModel(em)%Lag2Euler%GridBottomConc(i, j, np, ig) = 0.
                                        endif
                                    else
                                        ! Mass of sediment / m^2
                                        Me%EulerModel(em)%Lag2Euler%GridBottomConc(i, j, Sediment_ID, ig) = &
                                        Me%EulerModel(em)%Lag2Euler%GridBottomMass(i, j, Sediment_ID, ig) / DUX(I, j) / DVY(i, j)
                                    endif
                   
                                enddo 
                            else !all properties are written in mass / m^2
                                ! Mass / m^2
                                Me%EulerModel(em)%Lag2Euler%GridBottomConc(i, j, :, ig) = &
                                    Me%EulerModel(em)%Lag2Euler%GridBottomMass(i, j, :, ig) / DUX(I, j) / DVY(i, j)

                            endif

                        endif

                    enddo
                    enddo

                enddo g3

                !UnGets Horizontal Grid
                call UnGetHorizontalGrid(Me%EulerModel(em)%ObjHorizontalGrid, DUX,               &
                                       STAT = status)
                if (status /= SUCCESS_) call SetError(FATAL_, INTERNAL_, 'FillGridConcentration - ModuleLagrangianGlobal - ERR50')

                !UnGets Horizontal Grid
                call UnGetHorizontalGrid(Me%EulerModel(em)%ObjHorizontalGrid, DVY,               &
                                       STAT = status)
                if (status /= SUCCESS_) call SetError(FATAL_, INTERNAL_, 'FillGridConcentration - ModuleLagrangianGlobal - ERR60')

            endif cd3

        enddo d2


    end subroutine FillGridConcentration

    !--------------------------------------------------------------------------


    subroutine WriteGridConcentration(em)

        !Arguments-------------------------------------------------------------
        integer                                     :: em

        !Local-----------------------------------------------------------------
        integer, dimension(:, :, :   ), pointer     :: WaterPoints3D
        real, dimension(:, :, :), pointer           :: GridConc3D
        real, dimension(:, :   ), pointer           :: GridConc2D 
        integer                                     :: ig, p
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB
        character(StringLength)                     :: AuxChar, AuxChar2, AuxChar3

        !Shorten

        ig = 1

d2:     do ig = 1, Me%nGroups

            ILB    = Me%EulerModel(em)%Size%ILB
            IUB    = Me%EulerModel(em)%Size%IUB
            JLB    = Me%EulerModel(em)%Size%JLB
            JUB    = Me%EulerModel(em)%Size%JUB
            KLB    = Me%EulerModel(em)%Size%KLB
            KUB    = Me%EulerModel(em)%Size%KUB

            WS_ILB = Me%EulerModel(em)%WorkSize%ILB
            WS_IUB = Me%EulerModel(em)%WorkSize%IUB
            WS_JLB = Me%EulerModel(em)%WorkSize%JLB
            WS_JUB = Me%EulerModel(em)%WorkSize%JUB
            WS_KLB = Me%EulerModel(em)%WorkSize%KLB
            WS_KUB = Me%EulerModel(em)%WorkSize%KUB

            WaterPoints3D    => Me%EulerModel(em)%WaterPoints3D


            !Allocates auxiliar variable
            allocate (GridConc3D (ILB:IUB, JLB:JUB, KLB:KUB         ))

            if (Me%State%Deposition) allocate (GridConc2D (ILB:IUB, JLB:JUB      ))

            !Writes the Group to an auxiliar string
            write (AuxChar, fmt='(i3)') ig

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5(em), WS_ILB, WS_IUB, WS_JLB, WS_JUB,           &
                                  WS_KLB, WS_KUB)

            if (Me%nGroups == 1) then

                AuxChar2 = "/Results/"

            else

                AuxChar2 = "/Results/Group_"//trim(adjustl(AuxChar))//"/Data_3D/"

            endif

            CurrentProperty => Me%OriginDefault%FirstProperty


dp:         do p = 1, Me%OriginDefault%NProperties

ih:             if (CurrentProperty%WritesPropHDF) then

                    GridConc3D(:,:,:) = Me%EulerModel(em)%Lag2Euler%GridConc(:, :, :, p, ig)

                    AuxChar3 = trim(AuxChar2)//trim(CurrentProperty%Name)

                    !HDF 5
                    call HDF5WriteData        (Me%ObjHDF5(em),                                  &
                                               trim(AuxChar3),                              &
                                               trim(CurrentProperty%Name),                  &
                                               trim(CurrentProperty%Units),                 &
                                               Array3D = GridConc3D,                        &
                                               OutputNumber = Me%OutPut%NextOutPut)

                    if (CurrentProperty%T90ON) then

                        call ComputeT90Matrix(em, ig, CurrentProperty, GridConc3D)


                        AuxChar3 = trim(AuxChar2)//trim(CurrentProperty%T90Name)

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                          &
                                                   trim(AuxChar3),                          &
                                                   trim(CurrentProperty%T90Name),           &
                                                   's',                                     &
                                                   Array3D = GridConc3D,                    &
                                                   OutputNumber = Me%OutPut%NextOutPut)


                    endif

                    if (Me%OutPut%ConcMaxTracer) then

                        GridConc3D(:,:,:) = Me%EulerModel(em)%Lag2Euler%GridMaxTracer(:, :, :, p, ig)

                        AuxChar3 = "/Results/Group_"//trim(adjustl(AuxChar))//&
                                   "/Data_3D_MaxTracer/"//trim(CurrentProperty%Name)

                        !HDF 5
                        call HDF5WriteData    (Me%ObjHDF5(em),                                  &
                                               trim(AuxChar3),                              &
                                               trim(CurrentProperty%Name),                  &
                                               trim(CurrentProperty%Units),                 &
                                               Array3D = GridConc3D,                        &
                                               OutputNumber = Me%OutPut%NextOutPut)

                    endif

                    if (Me%State%Deposition) then

                        GridConc2D(:,:) = Me%EulerModel(em)%Lag2Euler%GridBottomConc(:, :, p, ig)

                        AuxChar3 = "/Results/Group_"//trim(adjustl(AuxChar))//&
                                   "/Bottom/"//trim(CurrentProperty%Name)

                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                              &
                                                   trim(AuxChar3),                          &
                                                   trim(CurrentProperty%Name),              &
                                                   trim(CurrentProperty%Units),             &
                                                   Array2D = GridConc2D,                    &
                                                   OutputNumber = Me%OutPut%NextOutPut)

                    endif

                endif ih

                CurrentProperty => CurrentProperty%Next


            enddo dp

            nullify(CurrentProperty)

            AuxChar3 = trim(AuxChar2)//"Number"


            GridConc3D(:,:,:) = Me%EulerModel(em)%Lag2Euler%GridTracerNumber(:, :, :, ig)

            !HDF 5
            call HDF5WriteData        (Me%ObjHDF5(em),                                      &
                                       trim(AuxChar3),                                  &
                                       "Number", "a",                                   &
                                       Array3D = GridConc3D,                            &
                                       OutputNumber = Me%OutPut%NextOutPut)

        enddo d2
!        enddo d1

        nullify(GridConc3D  )
        nullify(WaterPoints3D)
        nullify(GridConc2D  )  


    end subroutine WriteGridConcentration

    !--------------------------------------------------------------------------

    subroutine WriteFrequencyLag (em, p, ig)

        !Arguments-------------------------------------------------------------
        integer                                     :: em, p, ig
        !Local-----------------------------------------------------------------
        real, dimension(:, :, :), pointer           :: GridConc3D 
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB, iClass, STAT_CALL
        character(StringLength)                     :: AuxChar1, AuxChar2, AuxChar


        !Shorten
        ILB    = Me%EulerModel(em)%Size%ILB
        IUB    = Me%EulerModel(em)%Size%IUB
        JLB    = Me%EulerModel(em)%Size%JLB
        JUB    = Me%EulerModel(em)%Size%JUB
        KLB    = Me%EulerModel(em)%Size%KLB
        KUB    = Me%EulerModel(em)%Size%KUB

        WS_ILB = Me%EulerModel(em)%WorkSize%ILB
        WS_IUB = Me%EulerModel(em)%WorkSize%IUB
        WS_JLB = Me%EulerModel(em)%WorkSize%JLB
        WS_JUB = Me%EulerModel(em)%WorkSize%JUB
        WS_KLB = Me%EulerModel(em)%WorkSize%KLB
        WS_KUB = Me%EulerModel(em)%WorkSize%KUB

        !Allocates auxiliar variable
        allocate (GridConc3D (ILB:IUB, JLB:JUB, KLB:KUB         ))


        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5(em), WS_ILB, WS_IUB, WS_JLB, WS_JUB,               &
                              WS_KLB, WS_KUB)


i1:     if (Me%Statistic%OptionsStat(p)%Lag) then
            
            call GetStatisticClassesNumber(Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID(ig),&
                                           Me%Statistic%OptionsStat(p)%nClassesLag, STAT= STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'WriteFrequencyLag - ModuleLagrangianGlobal - ERR10'

            call GetStatisticClasses(Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID(ig),                          &
                                     Me%Statistic%OptionsStat(p)%ClassesLag, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'WriteFrequencyLag - ModuleLagrangianGlobal - ERR20'


            do iClass = 1, Me%Statistic%OptionsStat(p)%nClassesLag

                write(AuxChar1, fmt='(E9.2)')Me%Statistic%OptionsStat(p)%ClassesLag(iClass, 1)
                write(AuxChar2, fmt='(E9.2)')Me%Statistic%OptionsStat(p)%ClassesLag(iClass, 2)

                AuxChar = trim(adjustl(AuxChar1))//"_"//trim(adjustl(AuxChar2))

        
                GridConc3D(:,:,:) = 100. * Me%EulerModel(em)%PropStatistic(p)%FrequencyLag(:,:,:,iClass, ig)

                        

                call HDF5WriteData   (Me%ObjHDF5(em), "/Statistics/Lagrangian/"         &
                                      //trim(Me%Statistic%OptionsStat(p)%ID%Name)       &
                                      //"/Classes", trim(adjustl(AuxChar)),             &
                                      "-", Array3D = GridConc3D,                        &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFrequencyLag - ModuleLagrangianGlobal - ERR10'

            enddo

            call UnGetStatistic(Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID(ig),       &
                                Me%Statistic%OptionsStat(p)%ClassesLag, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'WriteFrequencyLag - ModuleLagrangianGlobal - ERR20'

            deallocate (GridConc3D )

        endif i1
        
       
    end subroutine WriteFrequencyLag

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine WriteOilGridThickness(em) 

        !Arguments-------------------------------------------------------------
        integer                                     :: em   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: ig
        real, dimension(:, :), pointer              :: OilGridThick2D 
        character(StringLength)                     :: AuxChar
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB


!em1:    do em =1, Me%EulerModelNumber

            !Shorten
            ILB    = Me%EulerModel(em)%Size%ILB
            IUB    = Me%EulerModel(em)%Size%IUB
            JLB    = Me%EulerModel(em)%Size%JLB
            JUB    = Me%EulerModel(em)%Size%JUB

            WS_ILB = Me%EulerModel(em)%WorkSize%ILB
            WS_IUB = Me%EulerModel(em)%WorkSize%IUB
            WS_JLB = Me%EulerModel(em)%WorkSize%JLB
            WS_JUB = Me%EulerModel(em)%WorkSize%JUB
            WS_KLB = Me%EulerModel(em)%WorkSize%KLB
            WS_KUB = Me%EulerModel(em)%WorkSize%KUB

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5(em), WS_ILB, WS_IUB, WS_JLB, WS_JUB,     &
                                  WS_KLB, WS_KUB)


            !Allocate GridVolume, GridMass
            allocate (OilGridThick2D (ILB:IUB, JLB:JUB))

Group:      do ig = 1, Me%nGroups

                !Writes the Group to an auxiliar string
                write (AuxChar, fmt='(i3)') ig

                CurrentOrigin => Me%FirstOrigin
CurrOr:         do while (associated(CurrentOrigin))


                    !Just writes the output if there are particle
                    if (CurrentOrigin%nParticle > 0) then
                        OilGridThick2D = Me%EulerModel(em)%OilSpreading(ig)%GridThickness * 1000.0 &
                                         * CurrentOrigin%VolumeOilTotal / CurrentOrigin%VolumeTotal 


                        !HDF 5
                        call HDF5WriteData        (Me%ObjHDF5(em),                       &
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

        !enddo em1

    end subroutine WriteOilGridThickness


    !--------------------------------------------------------------------------

    subroutine WriteOilGridConcentration(em) 

        !Arguments-------------------------------------------------------------
        integer                                     :: em

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: ig
        real, dimension(:, :), pointer              :: OilGridConc2D 
        character(StringLength)                     :: AuxChar
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB


        !em = 1

        !Shorten
        ILB    = Me%EulerModel(em)%Size%ILB
        IUB    = Me%EulerModel(em)%Size%IUB
        JLB    = Me%EulerModel(em)%Size%JLB
        JUB    = Me%EulerModel(em)%Size%JUB

        WS_ILB = Me%EulerModel(em)%WorkSize%ILB
        WS_IUB = Me%EulerModel(em)%WorkSize%IUB
        WS_JLB = Me%EulerModel(em)%WorkSize%JLB
        WS_JUB = Me%EulerModel(em)%WorkSize%JUB
        WS_KLB = Me%EulerModel(em)%WorkSize%KLB
        WS_KUB = Me%EulerModel(em)%WorkSize%KUB

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5(em), WS_ILB, WS_IUB, WS_JLB, WS_JUB,     &
                              WS_KLB, WS_KUB)


        !Allocate GridVolume, GridMass
        allocate (OilGridConc2D (ILB:IUB, JLB:JUB))

Group:  do ig = 1, Me%nGroups

            !Writes the Group to an auxiliar string
            write (AuxChar, fmt='(i3)') ig

            CurrentOrigin => Me%FirstOrigin
CurrOr:     do while (associated(CurrentOrigin))


                !Just writes the output if there are particle
                if (CurrentOrigin%nParticle > 0) then
                    !OilGridConc2D = CurrentOrigin%OilGridConcentration
                    OilGridConc2D = Me%EulerModel(em)%OilSpreading(ig)%OilGridConcentration

                    !HDF 5
                    call HDF5WriteData        (Me%ObjHDF5(em),                       &
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

    real(8) function GeographicCoordinates (em, Position, Direction)

        !Arguments-------------------------------------------------------------
   
        integer                                     :: em
        type (T_Position)                           :: Position
        integer                                     :: Direction

        !Local-----------------------------------------------------------------



        if (Me%EulerModel(em)%Grid%HaveLatLongGrid) then

             if (Direction == 1) then
                GeographicCoordinates = Position%CoordX
            else
                GeographicCoordinates = Position%CoordY
            endif

        else
            if (Direction == 1) then
                GeographicCoordinates = Me%EulerModel(em)%Grid%LongDefault
            else
                GeographicCoordinates = Me%EulerModel(em)%Grid%LatDefault
            endif
        endif


    end function GeographicCoordinates 

    !--------------------------------------------------------------------------

    subroutine WriteMonitorOutput 

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        type (T_Origin), pointer                    :: CurrentOrigin
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k, Box
        integer                                     :: STAT_CALL
        real,    dimension(:, :, :), pointer        :: OutputMatrix
        integer, dimension(:, :, :), pointer        :: MonitorBoxes
        integer                                     :: WS_ILB, WS_IUB, WS_JLB, WS_JUB
        integer                                     :: WS_KLB, WS_KUB, em

        em = 1

        !Shorten
        ILB    = Me%EulerModel(em)%Size%ILB
        IUB    = Me%EulerModel(em)%Size%IUB
        JLB    = Me%EulerModel(em)%Size%JLB
        JUB    = Me%EulerModel(em)%Size%JUB
        KLB    = Me%EulerModel(em)%Size%KLB
        KUB    = Me%EulerModel(em)%Size%KUB

        WS_ILB = Me%EulerModel(em)%WorkSize%ILB
        WS_IUB = Me%EulerModel(em)%WorkSize%IUB
        WS_JLB = Me%EulerModel(em)%WorkSize%JLB
        WS_JUB = Me%EulerModel(em)%WorkSize%JUB
        WS_KLB = Me%EulerModel(em)%WorkSize%KLB
        WS_KUB = Me%EulerModel(em)%WorkSize%KUB

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5(em), WS_ILB, WS_IUB, WS_JLB, WS_JUB,     &
                              WS_KLB, WS_KUB)


        allocate (OutputMatrix (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteMonitorOutput - ModuleLagrangianGlobal - ERR00'

        !Get the boxes
        call GetBoxes(Me%EulerModel(em)%ObjMonBox, MonitorBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteMonitorOutput - ModuleLagrangianGlobal - ERR01'

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
                    OutputMatrix(i, j, k) = Me%EulerModel(em)%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) / &
                                            Me%EulerModel(em)%Monitor%InstBoxVolume      (Box) * 100.
                endif
            enddo
            enddo
            enddo


            !HDF 5
            call HDF5WriteData         (Me%ObjHDF5(em),                                &
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
                    OutputMatrix(i, j, k) = Me%EulerModel(em)%Monitor%IntgVolumeByOrigin (Box, CurrentOrigin%ID) / &
                                            Me%EulerModel(em)%Monitor%IntgBoxVolume      (Box) * 100.
                endif
            enddo
            enddo
            enddo


            !HDF 5
            call HDF5WriteData         (Me%ObjHDF5(em),                                &
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
                    OutputMatrix(i, j, k) = Me%EulerModel(em)%Monitor%InstMassByOrigin (Box, CurrentOrigin%ID) / &
                                            Me%EulerModel(em)%Monitor%InstBoxMass      (Box) * 100.
                endif
            enddo
            enddo
            enddo


            !HDF 5
            call HDF5WriteData         (Me%ObjHDF5(em),                                &
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
        Me%EulerModel(em)%Monitor%NumberOfCellsPerBox = 0
        do Box = 1,  Me%EulerModel(em)%Monitor%NumberOfBoxes
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (MonitorBoxes(i, j, k) == Box .and. Me%EulerModel(em)%OpenPoints3D(i, j, k) == OpenPoint) then
                Me%EulerModel(em)%Monitor%NumberOfCellsPerBox(Box) = Me%EulerModel(em)%Monitor%NumberOfCellsPerBox(Box) + 1
            endif
        enddo
        enddo
        enddo
        enddo
        !Number of Cells to fill
        Me%EulerModel(em)%Monitor%NumberOfCellsFromOrigin = 0
        CurrentOrigin => Me%FirstOrigin
        do while (associated(CurrentOrigin))
            do Box = 1,  Me%EulerModel(em)%Monitor%NumberOfBoxes
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                if (MonitorBoxes(i, j, k) == Box) then
                    Me%EulerModel(em)%Monitor%NumberOfCellsFromOrigin    (Box, CurrentOrigin%ID) = &
                        int(Me%EulerModel(em)%Monitor%NumberOfCellsPerBox(Box) *                   &
                            Me%EulerModel(em)%Monitor%InstVolumeByOrigin (Box, CurrentOrigin%ID) / &
                            Me%EulerModel(em)%Monitor%InstBoxVolume      (Box))
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
            if (Box > 0 .and. Me%EulerModel(em)%OpenPoints3D(i, j, k) == OpenPoint) then

                OutputMatrix(i, j, k) = 0.
                do while (Me%EulerModel(em)%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) == 0)
                    CurrentOrigin => CurrentOrigin%Next
                    if (.not. associated(CurrentOrigin)) exit
                enddo

                if (associated(CurrentOrigin)) then
                    if (Me%EulerModel(em)%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) > 0) then
                        OutputMatrix(i, j, k) = CurrentOrigin%ID
                        Me%EulerModel(em)%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) = &
                            Me%EulerModel(em)%Monitor%NumberOfCellsFromOrigin(Box, CurrentOrigin%ID) - 1
                    endif
                endif
            endif
        enddo
        enddo
        enddo
            


        !HDF 5
        call HDF5WriteData         (Me%ObjHDF5(em),                          &
                                    "/Results/Monitor",                             &
                                    "VolumeContributedByOrigin",  "%",              &
                                    Array3D = OutputMatrix,                         &
                                    OutputNumber = Me%OutPut%NextOutPut, &
                                    STAT = STAT_CALL)

        !Unget The Boxes
        call UngetBoxDif(Me%EulerModel(em)%ObjMonBox, MonitorBoxes, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteMonitorOutput - ModuleLagrangianGlobal - ERR02'

        deallocate (OutputMatrix            )

    end subroutine WriteMonitorOutput

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    subroutine KillLagrangianGlobal (LagrangianID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: nUsers

        !----------------------------------------------------------------------
           
        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_)    

cd1:    if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mLAGRANGIAN_,  Me%InstanceID)

            if (nUsers == 0) then

                call DeallocateInstance

                LagrangianID  = 0
                STAT_         = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_



    end subroutine KillLagrangianGlobal

    !--------------------------------------------------------------------------
    subroutine DeallocateLagrangianGlobal (LagrangianID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LagrangianID
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        type (T_Property),  pointer                 :: CurrentProperty
        integer                                     :: nB
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL, nUsers, em
        logical                                     :: WriteFinal

        !----------------------------------------------------------------------
           
        STAT_ = UNKNOWN_

        call Ready(LagrangianID, ready_)    

cd1:    if (ready_ .NE. OFF_ERR_) then

            if (Me%WritesTimeSerie) then
d3:             do em = 1, Me%EulerModelNumber
                    call KillTimeSerie(Me%EulerModel(em)%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR10'
                enddo d3
            endif
            

            if (Me%State%Statistics) then

                call KillParticleStatistic 

            endif


            if (Me%OutPut%Write_) then
d4:             do em = 1, Me%EulerModelNumber
                    !Closes the transient result file
                    call KillHDF5          (Me%ObjHDF5(em), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR20'
                enddo d4
            endif

            WriteFinal = .true.

            if (Me%RunOnline) then
!#ifdef _CGI_
                WriteFinal = .false.
!#endif
            endif
            !Writes the Final Output
            if (WriteFinal) call WriteFinalPartic(Final = .true.)

d1:         do em = 1, Me%EulerModelNumber

                !Kills the Light
                if (Me%State%WQM .or. Me%State%T90Compute) then
                    call KillLight(em) 
                endif

                !Kills Monitoring
                if (Me%State%Monitor) then

                    do nB = 1, Me%EulerModel(em)%Monitor%NumberOfBoxes
                        call KillTimeSerie (Me%EulerModel(em)%Monitor%ObjTimeSerie(nB), STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR30'
                    enddo

                    deallocate (Me%EulerModel(em)%Monitor%InstBoxVolume      )
                    deallocate (Me%EulerModel(em)%Monitor%InstVolumeByOrigin )
                    deallocate (Me%EulerModel(em)%Monitor%InstBoxMass        )
                    deallocate (Me%EulerModel(em)%Monitor%InstMassByOrigin   )
                    deallocate (Me%EulerModel(em)%Monitor%IntgBoxVolume      )
                    deallocate (Me%EulerModel(em)%Monitor%IntgVolumeByOrigin )

                    call KillBoxDif (Me%EulerModel(em)%ObjMonBox, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR40'
                
                endif

                if (Me%State%EulerianMonitor) then

                    deallocate(Me%EulerModel(em)%EulerianMonitor%Mass)

                    call KillBoxDif (Me%EulerModel(em)%ObjEulMonBox, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR50'

                end if

                !Kills Oil-Beaching Properties
                if (Me%State%AssociateBeachProb) then
                    if (Me%State%HaveBeachingProbBox) then

                        call KillBoxDif(Me%EulerModel(em)%ObjBeachingProbBox, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR60'

                    end if

                    deallocate (Me%EulerModel(em)%BeachingProbability)

                end if
               

                if (Me%State%Deposition)  then

                    deallocate   (Me%EulerModel(em)%Lag2Euler%TauErosionGrid)
                    deallocate   (Me%EulerModel(em)%Lag2Euler%MassSedGrid   )

                endif

                if (Me%State%Filtration)  then

                    deallocate(Me%EulerModel(em)%RelativeMassFilter)

                    deallocate(Me%EulerModel(em)%MassFiltered)

                endif

            enddo d1

            if (Me%State%Oil)  then
                call DeAllocateOil
            endif

            !Kill the OriginList
            call DeallocateOriginList(Me%FirstOrigin, Me%nOrigins)

            if (associated(Me%OriginDefault)) then
                CurrentProperty => Me%OriginDefault%FirstProperty
                do while (associated(CurrentProperty))

                    if (CurrentProperty%T90Variable .and. CurrentProperty%T90Var_Method == FromTimeSerie) then
                        call KillTimeSerie(CurrentProperty%TimeSerieT90, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR70'
                    endif

                    CurrentProperty => CurrentProperty%Next
                enddo

                nullify(CurrentProperty)

            endif

d2:         do em = 1, Me%EulerModelNumber

                !Kill Assimilation
                if (Me%Overlay) then

                    call KillAssimilation (Me%EulerModel(em)%ObjAssimilation,  STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR70'

                    deallocate (Me%EulerModel(em)%OverLay%VelUFinal)
                    deallocate (Me%EulerModel(em)%OverLay%VelVFinal)

                endif

                if (Me%EulerModel(em)%ObjBoxDif /= 0) then
                    call KillBoxDif (Me%EulerModel(em)%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR80'
                endif

                !External Modules
                nUsers = DeassociateInstance (mTIME_,           Me%EulerModel(em)%ObjTime           )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR90'

                nUsers = DeassociateInstance (mGRIDDATA_,       Me%EulerModel(em)%ObjGridData       )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR100'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%EulerModel(em)%ObjHorizontalMap  )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR110'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%EulerModel(em)%ObjHorizontalGrid )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR120'

                nUsers = DeassociateInstance (mGEOMETRY_,       Me%EulerModel(em)%ObjGeometry       )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR130'

                nUsers = DeassociateInstance (mMAP_,            Me%EulerModel(em)%ObjMap            )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR140'

                nUsers = DeassociateInstance (mHYDRODYNAMIC_,   Me%EulerModel(em)%ObjHydrodynamic   )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR150'

                nUsers = DeassociateInstance (mTURBULENCE_,     Me%EulerModel(em)%ObjTurbulence     )
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR160'

                nUsers = DeassociateInstance (mWATERPROPERTIES_,Me%EulerModel(em)%ObjWaterProperties)
                if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR170'

#ifndef _WAVES_
                if(Me%EulerModel(em)%ObjWaves /= 0)then
                    nUsers = DeassociateInstance (mWAVES_,  Me%EulerModel(em)%ObjWaves)
                    if (nUsers == 0) stop 'DeallocateLagrangianGlobal - ModuleLagrangianGlobal - ERR180'
                end if
#endif
            enddo d2

            STAT_         = SUCCESS_

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_



    end subroutine DeallocateLagrangianGlobal

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

        deallocate(Me%EulerModel)

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
                               "_"//trim(TimeToString(Me%Now))//".fin")


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
            write (UnitID) CurrentOrigin%Position%ModelID
            write (UnitID) CurrentOrigin%State%Oil
            write (UnitID) CurrentOrigin%State%Deposition
            write (UnitID) CurrentOrigin%State%Age
            write (UnitID) CurrentOrigin%State%ComputePlume
            write (UnitID) CurrentOrigin%State%PlumeShear
            write (UnitID) CurrentOrigin%State%FarFieldBuoyancy
            write (UnitID) CurrentOrigin%nParticle
            write (UnitID) CurrentOrigin%nProperties

!            If (CurrentOrigin%State%Oil)   then
!                write (UnitID) CurrentOrigin%OilSpreading%GridThickness                      
!                write (UnitID) CurrentOrigin%OilGridConcentration
!            endif

            !Writes Particle Information
            CurrentPartic => CurrentOrigin%FirstPartic
            do while (associated(CurrentPartic))

                write (UnitID) CurrentPartic%ID

                write (UnitID) CurrentPartic%Position%ModelID

                write (UnitID) CurrentPartic%Position%I
                write (UnitID) CurrentPartic%Position%J
                write (UnitID) CurrentPartic%Position%K
                write (UnitID) CurrentPartic%Position%CellI
                write (UnitID) CurrentPartic%Position%CellJ
                write (UnitID) CurrentPartic%Position%CellK
                write (UnitID) CurrentPartic%Position%CoordX
                write (UnitID) CurrentPartic%Position%CoordY
                write (UnitID) CurrentPartic%Position%CartX
                write (UnitID) CurrentPartic%Position%CartY
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
        integer                                     :: Dummy, em

        !Opens File
        call UnitsManager (UnitID, OPEN_FILE, STAT = STAT_CALL)

        open   (unit = UnitID, file = Me%Files%Initial, status = 'old', form = 'unformatted')


        read(UnitID) OldOrigins

d1:     do nO = 1, OldOrigins

            call AllocateNewOrigin (NewOrigin)
            
            NewOrigin%Old   = ON 

            !Reads Origin Information
            read (UnitID) NewOrigin%Name
            read (UnitID) NewOrigin%ID
            read (UnitID) NewOrigin%Position%ModelID
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
!                allocate (NewOrigin%OilSpreading%GridThickness(EulerModel%Size%ILB:            &
!                                                  EulerModel%Size%IUB,            &
!                                                  EulerModel%Size%JLB:            &
!                                                  EulerModel%Size%JUB),           &
!                                                  STAT = STAT_CALL)

!                if (STAT_CALL /= SUCCESS_) stop 'ReadFinalPartic - ModuleLagrangianGlobal - ERR04'

!                read (UnitID) NewOrigin%OilSpreading%GridThickness
            
!                !Allocates GridThickness
!                allocate (NewOrigin%OilGridConcentration(EulerModel%Size%ILB:            &
!                                                         EulerModel%Size%IUB,            &
!                                                         EulerModel%Size%JLB:            &
!                                                         EulerModel%Size%JUB),           &
!                                                         STAT = STAT_CALL)

!                if (STAT_CALL /= SUCCESS_) stop 'ReadFinalPartic - ModuleLagrangianGlobal - ERR04a'

!                read (UnitID) NewOrigin%OilGridConcentration

            end if

            !Reads Particle Information
d2:         do nP = 1, nParticle

                call AllocateNewParticle (NewParticle, nProperties, NewOrigin%NextParticID)

                read (UnitID) Dummy !Does not read NewParticle%ID any more
                                    !This is due to an error in attributing 
                                    !identical ID's to different particles
                                    !This error was not completely understood
                                    !but the problem was solved. If you find 
                                    !any problem regarding memory errors and 
                                    !particle ID's please warn Luis or Frank

                read (UnitID) NewParticle%Position%ModelID

                read (UnitID) NewParticle%Position%I
                read (UnitID) NewParticle%Position%J
                read (UnitID) NewParticle%Position%K
                read (UnitID) NewParticle%Position%CellI
                read (UnitID) NewParticle%Position%CellJ
                read (UnitID) NewParticle%Position%CellK
                read (UnitID) NewParticle%Position%CoordX
                read (UnitID) NewParticle%Position%CoordY
                read (UnitID) NewParticle%Position%CartX
                read (UnitID) NewParticle%Position%CartY
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

            enddo d2

            !Reads Property Information
d3:         do nP = 1, nProperties

                call AllocateNewProperty (NewProperty)
                read (UnitID) NewProperty%Name

                call InsertPropertyToList(NewOrigin, NewProperty, SetStates = .false.)

            enddo d3

            !Reads OilInformation
i3:         if (NewOrigin%State%Oil) then
                
                em = NewOrigin%Position%ModelID

                !Starts Oil
                call StartOil(OilID             = NewOrigin%ObjOil,                     &
                              TimeID            = Me%EulerModel(em)%ObjTime,            &
                              EnterDataID       = Me%ObjEnterData,                      &
                              HorizontalGridID  = Me%EulerModel(em)%ObjHorizontalGrid,  &
                              GeometryID        = Me%EulerModel(em)%ObjGeometry,        &
                              MapID             = Me%EulerModel(em)%ObjMap,             &
                              DT                = Me%DT_PARTIC,                         &     
                              ContCalc          = NewOrigin%Old,                        &
                              ExtractType       = FromBlockInBlock,                     &                   
                              STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFinalPartic - ModuleLagrangianGlobal - ERR01'

                call ReadFinalOil(NewOrigin%ObjOil, UnitID) 

            end if i3


            call InsertOriginToList (Me%FirstOldOrigin, NewOrigin, Me%nOldOrigins)

        enddo d1

        !Closes File
        call UnitsManager (UnitID, CLOSE_FILE, STAT = STAT_CALL)


    end subroutine ReadFinalPartic

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine KillLight (em)

        !Arguments-------------------------------------------------------------
        integer                                     :: em

        !Local-----------------------------------------------------------------
        integer                                     :: ig, STAT_CALL

        !Begin-----------------------------------------------------------------

        !Updates light
        do ig = 1, Me%nGroups

            call KillLightExtinction(Me%EulerModel(em)%Light(ig)%ObjLightExtinction, STAT = STAT_CALL)

            deallocate (Me%EulerModel(em)%Light(ig)%TopRadiationCells)

        enddo

    end subroutine KillLight

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine KillParticleStatistic 

        !Arguments-------------------------------------------------------------
   

        !Local-----------------------------------------------------------------
        integer                                     :: nProp, em, STAT_CALL, p, ig

        !Begin-----------------------------------------------------------------

d1:     do em = 1, Me%EulerModelNumber 

            nProp   =  Me%Statistic%PropNumber
            
d2:         do p = 1, nProp
                                
d3:             do ig = 1, Me%NGroups

                    call WriteFrequencyLag(em, p, ig)
                
                    call KillStatistic    (Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID(ig), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillParticleStatistic - ModuleLagrangianGlobal - ERR10'

                    if (Me%OutPut%ConcMaxTracer) then

                        call KillStatistic(Me%EulerModel(em)%PropStatistic(p)%Statistic2_ID(ig), STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillParticleStatistic - ModuleLagrangianGlobal - ERR20'

                    endif
                    
                enddo d3

                deallocate(Me%EulerModel(em)%PropStatistic(p)%Statistic1_ID)

                if (Me%OutPut%ConcMaxTracer) then

                    deallocate(Me%EulerModel(em)%PropStatistic(p)%Statistic2_ID)

                endif
                    


                if (Me%Statistic%OptionsStat(p)%Lag) then
                    deallocate(Me%EulerModel(em)%PropStatistic(p)%FrequencyLag)
                endif

            enddo d2

            deallocate(Me%EulerModel(em)%PropStatistic)

        enddo d1

    end subroutine KillParticleStatistic

    !--------------------------------------------------------------------------

    subroutine KillLag2Euler

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer             :: em
        !Begin------------------------------------------------------------------
                
d1:     do em = 1, Me%EulerModelNumber 

            !deAllocate 
            deallocate (Me%EulerModel(em)%Lag2Euler%GridVolume        )
            deallocate (Me%EulerModel(em)%Lag2Euler%GridTracerNumber  )

            deallocate (Me%EulerModel(em)%Lag2Euler%GridMass          )
            deallocate (Me%EulerModel(em)%Lag2Euler%GridConc          )

            deallocate (Me%EulerModel(em)%Lag2Euler%MeanConc          )
            deallocate (Me%EulerModel(em)%Lag2Euler%AmbientConc       )
            deallocate (Me%EulerModel(em)%Lag2Euler%MinConc           )
            deallocate (Me%EulerModel(em)%Lag2Euler%MassVolCel        )

            if (Me%State%Deposition) then
                deallocate (Me%EulerModel(em)%Lag2Euler%GridBottomMass)
                deallocate (Me%EulerModel(em)%Lag2Euler%GridBottomConc)
            endif

            if (Me%OutPut%ConcMaxTracer) then
                deallocate (Me%EulerModel(em)%Lag2Euler%GridMaxTracer)
            endif

        enddo d1

        !----------------------------------------------------------------------

    end subroutine KillLag2Euler

    !--------------------------------------------------------------------------


    !------------------------------------------------------------------------------

    subroutine DeAllocateOil

        !Local---------------------------------------------------------------------
        integer                     :: em, ig, STAT_CALL                          

        !Begin---------------------------------------------------------------------

d1:     do em =1, Me%EulerModelNumber 

            deallocate (Me%EulerModel(em)%OilSpreading, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeAllocateOil - ModuleLagrangianGlobal - ERR10'

d2:         do ig = 1, Me%NGroups

            !Allocates GridThickness
            deallocate (Me%EulerModel(em)%OilSpreading(ig)%GridThickness, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeAllocateOil - ModuleLagrangianGlobal - ERR20'


            !Allocates OilGridConcentration
            deallocate (Me%EulerModel(em)%OilSpreading(ig)%OilGridConcentration, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'DeAllocateOil - ModuleLagrangianGlobal - ERR30'

            deallocate (Me%EulerModel(em)%OilSpreading(ig)%VelocityX, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'DeAllocateOil - ModuleLagrangianGlobal - ERR40'

            deallocate (Me%EulerModel(em)%OilSpreading(ig)%VelocityY, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'DeAllocateOil - ModuleLagrangianGlobal - ERR50'

            deallocate (Me%EulerModel(em)%OilSpreading(ig)%AreaFlag, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'DeAllocateOil - ModuleLagrangianGlobal - ERR60'

            enddo d2

        enddo d1

    end subroutine DeAllocateOil

    !------------------------------------------------------------------------------

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
            stop 'ModuleLagrangianGlobal - LocateObjLagrangian - ERR01'

    end subroutine LocateObjLagrangian

    !--------------------------------------------------------------------------

    subroutine ReadLockHorizontalGrid(EulerModel)
        
        !Arguments-------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, GEOG, UTM, MIL_PORT
        integer                                     :: SIMPLE_GEOG, GRID_COORD, NLRD, CoordType 
        
        !----------------------------------------------------------------------
        

        !XX, YY
        call GetHorizontalGrid (EulerModel%ObjHorizontalGrid, XX = EulerModel%XX,       &
                                YY = EulerModel%YY, DZX = EulerModel%DZX,               &
                                DZY = EulerModel%DZY, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangianGlobal - ERR20'

        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(EulerModel%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangianGlobal - ERR30'


        if(CoordType == UTM .or. CoordType == MIL_PORT .or.                             &
           CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(EulerModel%ObjHorizontalGrid,                        &
                                   XX_IE = EulerModel%XX_IE,                            &
                                   YY_IE = EulerModel%YY_IE,                            &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangianGlobal - ERR40'

        else

            call GetGridLatitudeLongitude(EulerModel%ObjHorizontalGrid,                 &
                                          GridLatitudeConn  = EulerModel%YY_IE,         &
                                          GridLongitudeConn = EulerModel%XX_IE,         &
                                          STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangianGlobal - ERR50'

        end if

        call GetGridCellArea (EulerModel%ObjHorizontalGrid, EulerModel%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockHorizontalGrid - ModuleLagrangianGlobal - ERR60'


    end subroutine ReadLockHorizontalGrid

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar

        !Local-----------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel
        integer                                     :: STAT_CALL, em
        
        !----------------------------------------------------------------------
        
        !Gets Time
        call GetComputeCurrentTime(Me%ExternalVar%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)              
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR10'

        Me%ExternalVar%RunPeriod = Me%ExternalVar%Now - Me%ExternalVar%BeginTime

em1:    do em =1, Me%EulerModelNumber

            EulerModel => Me%EulerModel(em)

            call ReadLockHorizontalGrid(EulerModel)


    i1:     if (.not. Me%RunOnlyMov2D) then

                !Gets Bathymetry
                call GetGridData          (EulerModel%ObjGridData, EulerModel%Bathymetry, STAT = STAT_CALL)     
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR20'


                !Gets ExteriorPoints 2D
                call GetBoundaries      (EulerModel%ObjHorizontalMap, EulerModel%BoundaryPoints2D, &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR30'


                !WaterColumn
                call GetGeometryWaterColumn(EulerModel%ObjGeometry, EulerModel%WaterColumn, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR40'


                !SZZ, DWZ    
                call GetGeometryDistances(EulerModel%ObjGeometry,                       & 
                                          SZZ         = EulerModel%SZZ,                 &
                                          ZCellCenter = EulerModel%ZCellCenter,         &
                                          DWZ         = EulerModel%DWZ,                 &
                                          DWZ_Xgrad   = EulerModel%DWZ_Xgrad,           &
                                          DWZ_Ygrad   = EulerModel%DWZ_Ygrad,           &
                                          STAT        = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR50'

                !VolumeZ
                call GetGeometryVolumes(EulerModel%ObjGeometry,                         &
                                        VolumeZ = EulerModel%VolumeZ,                   &
                                        STAT    = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR60'

                !kFloorZ
                call GetGeometryKFloor (EulerModel%ObjGeometry,                         &
                                        Z       = EulerModel%kFloor,                    &
                                        STAT    = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR70'


                !WaterPoints3D
                call GetWaterPoints3D(EulerModel%ObjMap, EulerModel%Waterpoints3D, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR80'


                !LandPoints3D
                call GetLandPoints3D(EulerModel%ObjMap, EulerModel%LandPoints3D, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR90'


                !OpenPoints3D
                call GetOpenPoints3D(EulerModel%ObjMap, EulerModel%OpenPoints3D, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR100'


                !Compute faces
                call GetComputeFaces3D(EulerModel%ObjMap,                               &
                                       ComputeFacesU3D = EulerModel%ComputeFaces3D_U,   &
                                       ComputeFacesV3D = EulerModel%ComputeFaces3D_V,   &
                                       STAT= STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR110'


                !Lupward, Ldownward
                call GetMixingLengthVertical(EulerModel%ObjTurbulence,                  &
                                             Lupward   = EulerModel%Lupward,            &
                                             Ldownward = EulerModel%Ldownward,          &
                                             STAT      = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR120'


                !MixingLengthX, MixingLengthY
                call GetMixingLengthHorizontal(EulerModel%ObjTurbulence,                &
                                               MixingLengthX = EulerModel%MixingLengthX,&
                                               MixingLengthY = EulerModel%MixingLengthY,&
                                               STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR130'

        
                !Velocity_U, Velocity_V
                call GetHorizontalVelocity(EulerModel%ObjHydrodynamic,                  &
                                           Velocity_U = EulerModel%Velocity_U,          &
                                           Velocity_V = EulerModel%Velocity_V,          &
                                           STAT       = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR140'

        
                !Velocity_W
                call GetVerticalVelocity(EulerModel%ObjHydrodynamic,                    &
                                         Velocity_W      = EulerModel%Velocity_W,       &
                                         STAT            = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleLagrangianGlobal - ERR150'

            else i1

            !Allocate 3D matrixes needed

            endif i1

        enddo em1

        nullify(EulerModel)


        !----------------------------------------------------------------------

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadLockEulerianDensity

        !Local-----------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel
        integer                                     :: STAT_CALL, em
        
        !----------------------------------------------------------------------

em1:    do em =1, Me%EulerModelNumber

            EulerModel => Me%EulerModel(em)

            !Density
            if (Me%State%FarFieldBuoyancy  .or. Me%State%ComputePlume) then
                call GetSigma       (EulerModel%ObjWaterProperties,                     &
                                     Sigma             = EulerModel%SigmaDensity,       &
                                     CurrentTime       = Me%Now,            &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockEulerianDensity - ModuleLagrangianGlobal - ERR10'
            else
                nullify(EulerModel%SigmaDensity)
            endif
    
    
            if (Me%State%Oil) then
                call GetDensity     (EulerModel%ObjWaterProperties,                     &
                                     Density           = EulerModel%Density,            &
                                     CurrentTime       = Me%Now,            &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockEulerianDensity - ModuleLagrangianGlobal - ERR20'
            else
                nullify(EulerModel%Density)
            endif

        enddo em1

        nullify(EulerModel)

    end subroutine ReadLockEulerianDensity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine ReadUnLockHorizontalGrid(EulerModel)
    
        !Arguments-------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------

        !XX, YY
        call UnGetHorizontalGrid (EulerModel%ObjHorizontalGrid, EulerModel%XX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockHorizontalGrid - ModuleLagrangianGlobal - ERR10'

        call UnGetHorizontalGrid (EulerModel%ObjHorizontalGrid, EulerModel%YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockHorizontalGrid - ModuleLagrangianGlobal - ERR20'

        !XX_IE, YY_IE
        call UnGetHorizontalGrid (EulerModel%ObjHorizontalGrid, EulerModel%XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR30'

        call UnGetHorizontalGrid (EulerModel%ObjHorizontalGrid, EulerModel%YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR40'

        call UnGetHorizontalGrid (EulerModel%ObjHorizontalGrid, EulerModel%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR50'

        call UnGetHorizontalGrid (EulerModel%ObjHorizontalGrid, EulerModel%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR60'

        call UnGetHorizontalGrid (EulerModel%ObjHorizontalGrid, EulerModel%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR70'

    end subroutine ReadUnLockHorizontalGrid

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel
        integer                                     :: STAT_CALL, em

        !Begin-----------------------------------------------------------------

em1:    do em =1, Me%EulerModelNumber

            EulerModel => Me%EulerModel(em)

            call ReadUnLockHorizontalGrid(EulerModel)

i1:         if (.not. Me%RunOnlyMov2D) then

                !Gets Bathymetry
                call UnGetGridData (EulerModel%ObjGridData, EulerModel%Bathymetry, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR50'


                !Gets ExteriorPoints 2D
                call UngetHorizontalMap (EulerModel%ObjHorizontalMap, EulerModel%BoundaryPoints2D,   &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR60'


                !WaterColumn
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%WaterColumn, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR70'

                !SZZ
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%SZZ, STAT  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR80'

                !ZCellCenter
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%ZCellCenter, STAT  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR90'

                !DWZ    
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%DWZ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR100'

                !VolumeZ
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%VolumeZ, STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR110'

                !kFloorZ
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%kFloor, STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR120'

                !DWZ_Xgrad
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%DWZ_Xgrad, STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR122'

                !DWZ_Ygrad
                call UnGetGeometry (EulerModel%ObjGeometry, EulerModel%DWZ_Ygrad, STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR124'


                !WaterPoints3D
                call UnGetMap      (EulerModel%ObjMap, EulerModel%Waterpoints3D, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR130'


                !LandPoints3D
                call UnGetMap      (EulerModel%ObjMap, EulerModel%LandPoints3D, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR140'


                !OpenPoints3D
                call UnGetMap      (EulerModel%ObjMap, EulerModel%OpenPoints3D, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR150'

                !Compute faces U
                call UnGetMap      (EulerModel%ObjMap, EulerModel%ComputeFaces3D_U, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR160'

                !Compute faces V
                call UnGetMap      (EulerModel%ObjMap, EulerModel%ComputeFaces3D_V, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR170'


                !Lupward
                call UngetTurbulence(EulerModel%ObjTurbulence, EulerModel%Lupward, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR180'

                !Ldownward
                call UngetTurbulence(EulerModel%ObjTurbulence, EulerModel%Ldownward, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR190'


                !MixingLengthX
                call UngetTurbulence(EulerModel%ObjTurbulence, EulerModel%MixingLengthX, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR200'


                !MixingLengthY
                call UngetTurbulence(EulerModel%ObjTurbulence, EulerModel%MixingLengthY, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR210'
        

                !Velocity_U
                call UngetHydrodynamic (EulerModel%ObjHydrodynamic, EulerModel%Velocity_U, STAT  = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR220'


                !Velocity_V
                call UngetHydrodynamic (EulerModel%ObjHydrodynamic, EulerModel%Velocity_V, STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR230'


                !Velocity_W
                call UngetHydrodynamic (EulerModel%ObjHydrodynamic, EulerModel%Velocity_W, STAT = STAT_CALL)                    
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleLagrangianGlobal - ERR240'

            else i1

            !Allocate 3D matrixes needed

            endif i1

        enddo em1

        nullify(EulerModel)
        !----------------------------------------------------------------------

    end subroutine ReadUnLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockEulerianDensity

        !Local-----------------------------------------------------------------
        type (T_EulerModel), pointer                :: EulerModel
        integer                                     :: STAT_CALL, em

        !Begin-----------------------------------------------------------------

em1:    do em =1, Me%EulerModelNumber

            EulerModel => Me%EulerModel(em)

            !Sigma Density
            if ((Me%State%FarFieldBuoyancy  .or. Me%State%ComputePlume) .and.           &
                associated(EulerModel%SigmaDensity)) then

                call UngetWaterProperties (EulerModel%ObjWaterProperties, EulerModel%SigmaDensity, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockEulerianDensity - ModuleLagrangianGlobal - ERR10'

            endif

            !Density
            if ((Me%State%Oil) .and.  associated(EulerModel%Density)) then
                call UngetWaterProperties (EulerModel%ObjWaterProperties, EulerModel%Density, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockEulerianDensity - ModuleLagrangianGlobal - ERR20'
            endif

        enddo em1

        nullify(EulerModel)
        !----------------------------------------------------------------------

    end subroutine ReadUnLockEulerianDensity

    !--------------------------------------------------------------------------

end Module ModuleLagrangianGlobal

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
