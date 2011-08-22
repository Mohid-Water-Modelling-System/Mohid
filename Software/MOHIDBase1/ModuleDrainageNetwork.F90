!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : DrainageNetwork
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig / Rosa Trancoso
! DESCRIPTION   : Module which simulates a 1D Drainage Network System
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
!
! Keywords read in the Data File
!
! Keyword                   : Data Type         Default     !Comment
!
! NETWORK_FILE              : char              -           !Path to drainage network file
! CHECK_NODES               : 0/1               [1]         !Ckeck nodes consistency
! CHECK_REACHES             : 0/1               [1]         !Check reaches consistency

! DISCHARGES                : 0/1               [0]         !Use module discharges (WWTP, etc)
! HYDRODYNAMIC_APROX        : int               [1]         !1 - KinematicWave, 2 - DiffusionWave, 3 - DynamicWave
! NUMERICAL_SCHEME          : int               [0]         !0 - ExplicitScheme, 1 - ImplicitScheme
! If ImplicitScheme ------------------------------------------------------------
! TIME_WEIGHT_FACTOR        : real              [0.7]       !Factor de ponderacao do peso dos termos explicitos e implicitos
! RELAXATION_FACTOR

! MASS_ERR                  : real(8)           [0.001]     !Max error in mass conservation
! GLOBAL_MANNING            : real              -           !Rugosity in Channels
! MIN_WATER_DEPTH           : real              [0.001]     !Min water depth in nodes (For h < MIN_WATER_DEPTH water stops flowing)
! MIN_WATER_DEPTH_PROCESS   : real              [0.01]      !Water Quality Process / Surface Fluxes shutdown
! INITIAL_WATER_DEPTH       : real              [0.0]       !Initial water depth
! TRANSMISSION_LOSSES       : 0/1               [0]         !If user wants to use transmission losses
! HYDRAULIC_CONDUCTIVITY    : real              -           !Hydraulic Conductivity to calculate transmission losses
! REMOVE_OVERTOP            : 0/1               [0]         !Removes Water if channels are overtoped
! STORM_WATER_MODEL_LINK    : 0/1               [0]         !If linked to a StormWaterModel
! MINIMUM_SLOPE             : real              [0.0]       !Minimum Slope for Kinematic Wave
! STABILIZE                 : 0/1               [0]         !Restart time iteration if high volume gradients
! STABILIZE_FACTOR          : real              [0.1]       !max gradient in time steps as fraction of old volume
! MAX_ITERATIONS            : int               [100]       !Max iterations for stabilized check
! DT_FACTOR                 : real              [0.8]       !Factor for DT Prediction
! MAX_DT_FLOOD              : real              [10.0]      !Max DT if channel water level exceeds full bank 
! AERATION_METHOD           : int               [-]         !1 - PoolAndRifle, 2 - ChannelControled_
! T90_DECAY_MODEL           : 0                 [1]         !0 - Constant, 1 - Canteras, 2 - Chapra
! T90                       : real              [7200.]     !if T90_DECAY_MODEL = Constant
! SHADING_FACTOR            : real              [1.]        !0-1 fraction of riparian shading
! FRACTION_SEDIMENT         : 0/1               [0]
! GLOBAL_TOXICITY           : char              ['SUM']     !Global Toxicity Computation Method : SUM,MAX,RISKRATIO
! GEO_CONVERSATION_FACTOR   : real              [1.]        !Lat to Meters rough estimation
! OUTPUT_TIME               : int int...        [-]
! DOWNSTREAM_BOUNDARY       : int               [1]         !0 - Dam, 1 - ZDG, 2 - CD, 3 - ImposedWaterLevel, 4 - ImposedVelocity
!   If ImposedWaterLevel--------------------------------------------------------
!       FILE_IN_TIME        : char              [NONE]      !NONE, TIMESERIE
!       DEFAULT_VALUE       : real              -           !Default value for water level at downstream boundary
!       If FILE_IN_TIME = TIMESERIE---------------------------------------------
!           FILENAME        : char              -           !Name of timeserie file for the downstream boundary
!           DATA_COLUMN     : int               -           !Number of column with data
!
! TIME_SERIE_LOCATION       : char              -           !Path to time serie especification nodes
! MAX_BUFFER_SIZE           : 1000
! COMPUTE_RESIDUAL          : 1 
! DT_OUTPUT_TIME            : 1200
! TIME_SERIE_BY_NODES       : 0/1               [0]         !Keyword to see if the user wants the time series to be written by 
                                                            !nodes, i.e.,
                                                            !One file per node, with all variables in the headers list
                                                            !if FALSE, its one file per variable with nodes in the headers.

!<BeginNodeTimeSerie>  / <EndNodeTimeSerie>


!<beginproperty>
!   NAME                    : cohesive sediment 
!   UNITS                   : mg/L 
!   DESCRIPTION             : cohesive sediment
!   DEFAULTVALUE            : 100.00
!   MIN_VALUE               : 0.0
!   ADVECTION_DIFUSION      : 1
!       ADVECTION_SCHEME    : 1       !Upwind
!       DIFFUSION_SCHEME    : 5       !CentralDif
!       DIFFUSIVITY         : 1E-8    !m2/s
!       VIRTUAL_COEF        : 0.01
!   WATER_QUALITY           : 0
!   BENTHOS                 : 0
!   DECAY_T90               : 0       !uses T90 decay model for fecal coliforms
!   DECAY_GENERIC           : 0       !uses generic decay (for now 1st order)
!                                       [2] -                                   
!   TIME_SERIE              : 1
!<endproperty>

!
!Network file ##################################################################
!
!<BeginNode>
!   ID                      : int               -           !Node ID number
!   COORDINATES             : real real         -           !Node coordinates
!   GRID_I                  : int               -           !I position of node, if grid
!   GRID_J                  : int               -           !J position of node, if grid
!   TERRAIN_LEVEL            : real              -           !Bottom level of cross section
!   MANNING_CHANNEL         : real         GLOBAL_MANNING   !Node rugosity
!   WATER_DEPTH             : real      INITIAL_WATER_DEPTH !Node initial water depth
!   CROSS_SECTION_TYPE      : int               [1]         !1 - Trapezoidal, 2 - TrapezoidalFlood, 3 - Tabular
!   1 - Trapezoidal, 2 - TrapezoidalFlood
!       BOTTOM_WIDTH        : real              -           !Bottom width of cross section
!       TOP_WIDTH           : real              -           !Top width of cross section
!       HEIGHT              : real              -           !Max height of cross section
!   2 - TrapezoidalFlood
!       MIDDLE_WIDTH        : real              -           !Middle width of cross section
!       MIDDLE_HEIGHT       : real              -           !Middle height of cross section
!   3 - Tabular
!       N_STATIONS          : integer           -           !number os stations that define the cross section
!       STATION             : real real ...     -           !station values
!       ELEVATION/LEVEL     : real real ...     -           !elevation values        
!<EndNode>    
!<BeginReach>
!   ID                      : int               -           !Reach ID Number
!   DOWNSTREAM_NODE         : int               -           !Downstream node ID
!   UPSTREAM_NODE           : int               -           !Upstream node ID
!   ACTIVE                  : boolean           -           !Active Reach. If Inactive, no flow is calculated
!<EndReach>
!
! EcoToxicity model ################################################################
! 
! Every toxic property must be discharged.
! Its concentration in the river network is set to 0.0.
! Discharge concentration must be equal to 1, because we are measuring the dilution
! D = 1 - C_new / C_ini
! the variable property%toxicity%concentration represents C/c_ini so it starts by being 1.
! This is not even close to a final version. 
! For more details, or sugestions/corrections, contact Rosa.
   

Module ModuleDrainageNetwork

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleHDF5
    use ModuleFunctions            , only: InterpolateValueInTime, ConstructPropertyID, ComputeT90_Chapra,  &
                                           ComputeT90_Canteras, LongWaveDownward, LongWaveUpward,           &
                                           LatentHeat, SensibleHeat, OxygenSaturation,                      &
                                           OxygenSaturationHenry, OxygenSaturationCeQualW2, AerationFlux,   &
                                           TimeToString, ChangeSuffix, DistanceBetweenTwoGPSPoints
    use ModuleTimeSerie            , only: StartTimeSerie, StartTimeSerieInput, WriteTimeSerieLine,         &
                                           GetTimeSerieValue, KillTimeSerie
    use ModuleStopWatch            , only: StartWatch, StopWatch
    use ModuleDischarges           , only: Construct_Discharges, GetDischargesNumber, GetDischargesNodeID,  &
                                           GetDischargeWaterFlow, GetDischargeConcentration, Kill_Discharges
    use ModuleLightExtinction      , only: ConstructLightExtinction, GetLightExtinctionOptions,             &
                                           GetRadiationPercentages, GetShortWaveExtinctionField,            &
                                           ModifyLightExtinctionField, UnGetLightExtinction,                &
                                           KillLightExtinction
    use ModuleInterface            , only: ConstructInterface, Modify_Interface, KillInterface, GetWQRatio

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructDrainageNetwork
    private ::      AllocateInstance
    private ::      ReadDataFile
    private ::      ConstructDownstreamBoundary
        
    private ::      ConstructNetwork
    private ::          ConstructNodeList
    private ::              CountTotalNodes
    private ::              ConstructNode 
    private ::                  InitializeTabularCrossSection
    private ::                      ComputeExtraArea
    private ::                          TrapezoidGeometry
    private ::              CheckNodesConsistency        
    private ::          ConstructReachList
    private ::              CountTotalReaches    
    private ::              ConstructReach    
    private ::              CheckReachesConsistency
    private ::              CalculateReaches
    private ::          ConnectNetwork
    private ::          OrderNodes
    private ::          WriteOrderedNodes
    private ::          CountOutlets

    private ::      ConstructPropertyList
    private ::          ConstructProperty
    private ::              ConstructPropertyValues
    private ::              InitializeProperty
    private ::          Add_Property
    private ::      CheckSelectedProp
    
    private ::      InitializeVariables
    private ::          ReadInitialFile
    private ::          InitializeNodes
    private ::              ComputeXSFromWaterDepth
    private ::              TabularGeometry
    private ::          InitializeReaches

    private ::      ConstructSubModules
    private ::          CoupleLightExtinction
    private ::          CoupleWaterQuality
    private ::          CoupleBenthos

    private ::      ConstructOutput
    private ::          ReadTimeSerieNodeList
    private ::          ConstructTimeSerieList
    private ::          ConstructTimeSeries    
    private ::              FillPropNameVector
    private ::          ConstructHDF5Output
    private ::      ConstructLog

    private ::      FindNodePosition
    private ::      FindReachPosition


    !Selector
    public  :: GetDrainageSize    
    public  :: GetChannelsID
    public  :: GetChannelsStationName
    public  :: GetChannelsWaterLevel
    public  :: GetChannelsBottomLevel
    public  :: GetChannelsSurfaceWidth
    public  :: GetChannelsBottomWidth
    public  :: GetChannelsBankSlope
    public  :: GetChannelsNodeLength
    public  :: GetChannelsOpenProcess
    public  :: GetHasProperties
    public  :: GetDNnProperties
    public  :: GetDNPropertiesIDByIdx   
    public  :: GetHasToxicity
    public  :: GetPropHasBottomFluxes
    public  :: GetNeedsRadiation
    public  :: GetNeedsAtmosphere
    public  :: GetNextDrainageNetDT
    public  :: GetVolumes
    public  :: GetDNConcentration
    public  :: GetDNMassBalance             !To Basin get the property mass balance values
    public  :: CheckDNProperty             
    public  :: UnGetDrainageNetwork
    public  :: SetAtmosphereDrainageNet     !To be called from MOHID Land
    public  :: SetAtmosphereRiverNet        !To be called from River Network
    public  :: SetPMPConcDN                 !DrainageNetwork gets the conc from Porous Media Properties
    public  :: SetRPConcDN                  !DrainageNetwork gets the conc from Runoff Properties
    public  :: SetGWFlowLayersToDN          !DrainageNetwork gets the Porous Media layers limits for GWFlow (faster process)
    private :: SearchProperty
        
    !Modifier
    public  :: FillOutPutMatrix
    public  :: ModifyDrainageNetwork
    private ::      ModifyDrainageNetLocal
    private ::          StoreInitialValues
    private ::          ModifyWaterDischarges
    private ::          ModifyWaterExchange
    private ::          ModifyTransmissionLosses
    private ::          UpdateCrossSections
    private ::              ComputeCrossSection
    private ::                  TrapezoidWaterHeight
    private ::                  TabularWaterLevel
    private ::                  ModifyDownstreamTimeSerie       
    private ::              UpdateReachCrossSection
    private ::          UpdateComputeFaces          
    private ::          UpdateOpenPoints    

    private ::          ModifyHydrodynamics
    private ::              ModifyReach
    private ::                  ComputeCriticalFlow
    private ::                  ComputeKinematicWave
    private ::                  ComputeStVenant  
    private ::                      HydroAdvection
    private ::              ModifyNode  
    private ::                  ComputeNodeInFlow
    private ::                  ComputeNodeOutFlow 
    private ::              VerifyMinimumVolume
    private ::              Cascade    
    private ::          ResetToInitialValues 
    private ::          TransportProperties
    private ::              Advection_Diffusion
    private ::                  ComputeAdvection
    private ::                  ComputeDiffusion
    private ::              SetMinimumConcentration
    private ::          ModifyTopRadiation
    private ::          ComputeSurfaceFluxes
    private ::          ColiformDecay
    private ::          ModifyToxicity
    private ::              ComputeToxicityForEachEffluent
    private ::          ModifyWaterQuality
    private ::          ModifyBenthos
    private ::          ComputeBottomFluxes
    private ::              ModifyShearStress
    private ::              ComputeErosionFluxes
    private ::              ComputeDepositionFluxes
    private ::                  SettlingVelocity
    private ::          UpdateChannelsDynamicMatrix
    private ::          ComputeNextDT
    private ::          WriteTimeSeries
    private ::              WriteTimeSeriesByNodes
    private ::              WriteTimeSeriesByProp
    private ::          HDF5Output
    private ::          MaxStationValues
    private ::          CalculateLoad
    private ::          CalculateTSS
    private ::          CalculateVSS
       
    !Destructor
    public  :: KillDrainageNetwork
    private ::      MaxStationValuesOutput
    private ::      WriteFinalFile
    private ::      Write_Errors_Messages                                                    
    

    !Management
    private ::      Ready
    private ::          LocateObjDrainageNetwork 
    
    !Interfaces----------------------------------------------------------------

    interface ModifyDrainageNetwork
        module procedure ModifyDrainageNetWithGrid
        module procedure ModifyDrainageNetWithoutGrid
    end interface

    interface UnGetDrainageNetwork
        module procedure UnGetDrainageNetworkR4
        module procedure UnGetDrainageNetworkI4
        module procedure UnGetDrainageNetworkA4
        module procedure UnGetDrainageNetwork1DR4
    end interface


    !Parameters------------------------------------------------------------------------------------
    character(StringLength), parameter              :: BeginNode            = '<BeginNode>'
    character(StringLength), parameter              :: EndNode              = '<EndNode>'
    character(StringLength), parameter              :: BeginReach           = '<BeginReach>'
    character(StringLength), parameter              :: EndReach             = '<EndReach>'
    character(StringLength), parameter              :: BeginNodeTimeSerie   = '<BeginNodeTimeSerie>'
    character(StringLength), parameter              :: EndNodeTimeSerie     = '<EndNodeTimeSerie>'
    character(StringLength), parameter              :: BeginReachTimeSerie  = '<BeginReachTimeSerie>'
    character(StringLength), parameter              :: EndReachTimeSerie    = '<EndReachTimeSerie>'
    
    character(LEN = StringLength), parameter        :: block_begin          = '<beginproperty>'
    character(LEN = StringLength), parameter        :: block_end            = '<endproperty>'
    
    !CrossSections
    integer, parameter                              :: Trapezoidal          = 1
    integer, parameter                              :: TrapezoidalFlood     = 2
    integer, parameter                              :: Tabular              = 3
    
    !DownstreamBoundary
    integer, parameter                              :: Dam                  = 0
    integer, parameter                              :: ZeroDepthGradient    = 1
    integer, parameter                              :: CriticalDepth        = 2
    integer, parameter                              :: ImposedWaterLevel    = 3
    integer, parameter                              :: ImposedVelocity      = 4
    

    !HydrodynamicApproximation
    integer, parameter                              :: KinematicWave        = 1 !Manning com declive de fundo
    integer, parameter                              :: DiffusionWave        = 2 !Manning com declive de superficie
    integer, parameter                              :: DynamicWave          = 3 !Todos os termos de StVenant

    !Toxicity Function types
    integer, parameter                              :: Saturation           = 1
    integer, parameter                              :: Linear               = 2
    integer, parameter                              :: RiskRatio            = 3
    
    !Variable downstream boundary
    integer, parameter                              :: None                 = 1
    integer, parameter                              :: ReadTimeSerie        = 2
    integer, parameter                              :: OpenMI               = 3

    !TimeSerie hydrodynamic properties
    integer, parameter                              :: pWaterDepth          = 1 
    integer, parameter                              :: pWaterLevel          = 2 
    integer, parameter                              :: pPercentageMaxVolume = 3 
    integer, parameter                              :: pVerticalArea        = 4
    integer, parameter                              :: pFlowToChannels      = 5 
    integer, parameter                              :: pVolume              = 6    
    integer, parameter                              :: pFlow                = 7 
    integer, parameter                              :: pVelocity            = 8 
    integer, parameter                              :: pGWFlowToChannels    = 9
    integer, parameter                              :: pPoolDepth           = 10
    integer, parameter                              :: pDT                  = 11
    integer, parameter                              :: pDTLocal             = 12

    integer, parameter                              :: BaseTimeSeries       = 12

    !OutputHydro
    integer, parameter                              :: pHydroTimeGradient   = 13
    integer, parameter                              :: pHydroAdvection      = 14
    integer, parameter                              :: pHydroPressure       = 15
    integer, parameter                              :: pHydroGravity        = 16
    integer, parameter                              :: pHydroFriction       = 17
 
    !T90 Calc Method
    integer, parameter                              :: Constant             = 0
    integer, parameter                              :: Canteras             = 1
    integer, parameter                              :: Chapra               = 2

    !O2 Aeration Method
    integer, parameter                              :: PoolAndRifle_        = 1
    integer, parameter                              :: ChannelControled_    = 2

    !TimeSerie hydrodynamic properties
    character(StringLength), parameter              :: Char_WaterDepth           = trim(adjustl('channel water depth'))
    character(StringLength), parameter              :: Char_WaterLevel           = trim(adjustl('channel water level'))
    character(StringLength), parameter              :: Char_PercentageMaxVolume  = trim(adjustl('percentage max volume'))
    character(StringLength), parameter              :: Char_VerticalArea         = trim(adjustl('vertical area'))
    character(StringLength), parameter              :: Char_FlowToChannels       = trim(adjustl('flow to channels'))
    character(StringLength), parameter              :: Char_Volume               = trim(adjustl('volume'))
    character(StringLength), parameter              :: Char_Flow                 = trim(adjustl('channel flow'))
    character(StringLength), parameter              :: Char_Velocity             = trim(adjustl('velocity'))
    character(StringLength), parameter              :: Char_GWFlowToChannels     = trim(adjustl('GW flow to channels'))
    character(StringLength), parameter              :: Char_PoolDepth            = trim(adjustl('pool water depth'))
    character(StringLength), parameter              :: Char_DT                   = trim(adjustl('DT'))
    character(StringLength), parameter              :: Char_DTLocal              = trim(adjustl('Local DT'))

    character(StringLength), parameter              :: Char_HydroTimeGradient    = trim(adjustl('hydro time gradient'))
    character(StringLength), parameter              :: Char_HydroAdvection       = trim(adjustl('hydro advection'))
    character(StringLength), parameter              :: Char_HydroPressure        = trim(adjustl('hydro pressure'))
    character(StringLength), parameter              :: Char_HydroGravity         = trim(adjustl('hydro gravity'))
    character(StringLength), parameter              :: Char_HydroFriction        = trim(adjustl('hydro friction'))

    integer, parameter                              :: UnitMax              = 80

    !Types---------------------------------------------------------------------

    type T_FlowFrequency
        type (T_Time)                               :: StartDate
        type (T_Time)                               :: StopDate
        real                                        :: MinimumFlow
    end type T_FlowFrequency

    type T_OutPut
        type (T_Time), dimension(:), pointer        :: OutTime
        type (T_Time), dimension(:), pointer        :: RestartOutTime
        integer                                     :: NextOutPut
        logical                                     :: Yes = .false.
        logical                                     :: WriteRestartFile     = .false.
        logical                                     :: RestartOverwrite     = .false.
        integer                                     :: NextRestartOutput    = 1
        logical                                     :: ComputeFlowFrequency = .false.
        type (T_FlowFrequency)                      :: FlowFrequency        
    end type T_OutPut

    type T_Files 
         character(PathLength)                      :: InputData
         character(PathLength)                      :: FinalFile
         character(PathLength)                      :: HDFFile
         character(PathLength)                      :: Initial
         character(PathLength)                      :: Network
         integer                                    :: ObjEnterDataNetwork   = 0
         integer                                    :: ObjEnterDataInitial   = 0
    end type T_Files

    type T_CrossSection
        integer                                     :: Form                     = null_int        
        real                                        :: BottomWidth              = null_real        
        real                                        :: TopWidth                 = null_real
        real                                        :: Slope                    = null_real     
        real                                        :: Height                   = null_real     !Total: from bottomlevel to surface       
        real                                        :: TerrainLevel             = null_real     !dado input da net
        real                                        :: BottomLevel              = null_real     !isto passa a ser calculado
        real                                        :: ManningCH                = null_real
        real                                        :: PoolDepth                = null_real    
        real                                        :: MiddleWidth              = null_real             
        real                                        :: MiddleHeight             = null_real
        real                                        :: SlopeTop                 = null_real
        !Tabular
        integer                                     :: IBottom                  = 0
        integer                                     :: NStations                = 0
        integer                                     :: NLevels                  = 0
        real, dimension(:), pointer                 :: Station                  !length NStations
        real, dimension(:), pointer                 :: Elevation                !length NStations
        real, dimension(:), pointer                 :: BankSlope                !length NStations

        real, dimension(:), pointer                 :: Level                    !length NLevels
        real, dimension(:), pointer                 :: LevelSlopeLeft           !length NLevels
        real, dimension(:), pointer                 :: LevelSlopeRight          !length NLevels
        real, dimension(:), pointer                 :: LevelBottomWidth         !length NLevels
        real, dimension(:), pointer                 :: LevelVerticalArea        !length NLevels
        real, dimension(:), pointer                 :: LevelWetPerimeter        !length NLevels       
        real, dimension(:), pointer                 :: LevelSurfaceWidth        !length NLevels       
       
    end type   T_CrossSection    

    type T_MaxValues
        real                                        :: Depth                 = null_real
        real                                        :: Flow                  = null_real
        real                                        :: Vel                   = null_real
        character(len=StringLength)                 :: Time                  = null_str
    end type

    type T_Node
        integer                                     :: ID                       = null_int
        real                                        :: X                        = null_real
        real                                        :: Y                        = null_real
        real                                        :: VerticalArea             = null_real
        real                                        :: WaterDepth               = null_real !cotas (inclui z bottom)
        real                                        :: InitialWaterDepth        = null_real
        real                                        :: WaterLevel               = null_real
        real(8)                                     :: VolumeNew                = null_real
        real(8)                                     :: VolumeOld                = null_real        
        real(8)                                     :: InitialVolumeOld         = null_real
        real(8)                                     :: InitialVolumeNew         = null_real
        real                                        :: VolumeMax                = null_real
        real                                        :: VolumeMaxTrapez1         = null_real
        real                                        :: VolumeMin                = null_real
        real                                        :: WetPerimeter             = null_real
        real                                        :: Length                   = null_real        
        real                                        :: SurfaceArea              = null_real
        real                                        :: SurfaceWidth             = null_real
        integer                                     :: GridI                    = null_int
        integer                                     :: GridJ                    = null_int
        integer                                     :: nUpstreamReaches         = 0
        integer                                     :: nDownstreamReaches       = 0
        integer                                     :: Order                    = null_int
        integer, dimension (:), pointer             :: UpstreamReaches
        integer, dimension (:), pointer             :: DownstreamReaches
        logical                                     :: TimeSerie                = .FALSE.
        logical                                     :: Discharges               = .FALSE.
        type (T_CrossSection)                       :: CrossSection                             
        character(len=StringLength)                 :: StationName              = ''
        real                                        :: SingCoef                 = 1.0
        type(T_MaxValues)                           :: Max
        real                                        :: EVTP                     = null_real !m/s evapotranspiration in pools
    end type  T_Node
    
    type T_Reach
        private
        integer                                     :: ID                       = null_int
        logical                                     :: Active                   = .true.
        character(LEN = StringLength)               :: Name                     
        integer                                     :: UpstreamNode             = null_int
        integer                                     :: DownstreamNode           = null_int
        real                                        :: Length                   = null_real    
        real                                        :: Slope                    = null_real
        real                                        :: FlowNew                  = 0.0        
        real                                        :: FlowOld                  = 0.0
        real                                        :: InitialFlowOld           = 0.0
        real                                        :: InitialFlowNew           = 0.0
        real                                        :: Velocity                 = 0.0        
        real                                        :: VerticalArea             = 0.0
        real                                        :: PoolVerticalArea         = 0.0
        real                                        :: HydraulicRadius          = 0.0
        real                                        :: Manning                  = 0.0
        logical                                     :: TimeSerie

        real                                        :: HydroTimeGradient        = 0.0
        real                                        :: HydroAdvection           = 0.0
        real                                        :: HydroPressure            = 0.0
        real                                        :: HydroGravity             = 0.0
        real                                        :: HydroFriction            = 0.0
        
        !Flow accumulation analisys
        real                                        :: InitialFlowAccTime       = 0.0
        real                                        :: FlowAccTime              = 0.0
        real                                        :: FlowAccPerc              = 0.0

    end type   T_Reach

    type T_TimeSerie
        integer                                                 :: ObjEnterData = 0
        logical                                                 :: ByNodes      = .false.
        character(PathLength)                                   :: Location        
        integer                                                 :: nNodes       = 0
        integer                                                 :: nProp        = 0
        integer                 , dimension (:), pointer        :: ObjTimeSerie        
        character(StringLength) , dimension (:), pointer        :: Name
        real                    , dimension (:), pointer        :: X
        real                    , dimension (:), pointer        :: Y        
        real, dimension(:), pointer                             :: DataLine        
    end type T_TimeSerie

    type       T_ExtVar
        real                                        :: DT
        logical                                     :: CoupledPMP = .false.
        logical                                     :: CoupledRP = .false.
    end type T_ExtVar

    type T_Downstream
        integer                                     :: Boundary
        integer                                     :: Evolution
        real                                        :: DefaultValue !WaterColumn in meters
        character(PathLength)                       :: FileName
        integer                                     :: DataColumn
        integer                                     :: ObjTimeSerie     = 0
     end type T_Downstream

    type T_Toxicity
        integer                                     :: Evolution        
        real                                        :: Slope
        real                                        :: EC50  !Concentration that causes 50% of effect (Tox = 0.5)        
        real, dimension (:), pointer                :: Field                    => null()
    end type    T_Toxicity

    type T_ComputeOptions
        logical                                     :: TimeSerie                = .false.    
        logical                                     :: Discharges               = .false.
        logical                                     :: Toxicity                 = .false.
        logical                                     :: T90_Decay                = .false.
        logical                                     :: Generic_Decay            = .false.
        logical                                     :: SurfaceFluxes            = .false.
        logical                                     :: BottomFluxes             = .false.
        logical                                     :: Erosion                  = .false.
        logical                                     :: Deposition               = .false.
        logical                                     :: AdvectionDiffusion       = .false.
        logical                                     :: WaterQuality             = .false.
        logical                                     :: Benthos                  = .false.
        logical                                     :: CeQualW2                 = .false.
        logical                                     :: Life                     = .false.
        logical                                     :: MinConcentration         = .false.
        logical                                     :: TopRadiation             = .false.
        logical                                     :: LightExtinction          = .false.
        logical                                     :: TransmissionLosses       = .false.
        logical                                     :: RemoveOverTop            = .false.
        logical                                     :: SumTotalConc             = .false.
        logical                                     :: ComputeLoad              = .false.
        logical                                     :: CalcFractionSediment     = .false.
        logical                                     :: EVTPFromReach            = .false.
        logical                                     :: StormWaterModelLink      = .false.
    end type T_ComputeOptions

    type       T_Coupling
         type(T_Time)                               :: NextCompute
         real                                       :: DT_Compute               = FillValueReal
         logical                                    :: Yes                      = .false.
         integer                                    :: NumberOfProperties       = 0
    end type T_Coupling  

    type       T_Coupled
         type(T_Coupling)                           :: WQM
         type(T_Coupling)                           :: CEQUALW2
         type(T_Coupling)                           :: Life
         type(T_Coupling)                           :: Benthos
    end type T_Coupled

    type T_MassBalance
        real(8)                                     :: TotalStoredMass
        real(8)                                     :: TotalDischargeMass
        real(8)                                     :: TotalOutFlowMass
    end type T_MassBalance

    type        T_Property
        type (T_PropertyID)                         :: ID
        type (T_ComputeOptions)                     :: ComputeOptions

        !Concentrations
        real, dimension (:), pointer                :: Concentration            => null()
        real, dimension (:), pointer                :: ConcentrationOld         => null()
        real, dimension (:), pointer                :: InitialConcentration     => null()
        real, dimension (:), pointer                :: InitialConcentrationOld  => null()
        real, dimension (:), pointer                :: MassCreated              => null()   !kg
        real, dimension (:), pointer                :: OverLandConc             => null()
        real, dimension (:), pointer                :: GWaterConc               => null()
        real, dimension (:, :, :), pointer          :: GWaterConcLayers         => null()   !for computation by layers
        real, dimension (:), pointer                :: DWaterConc               => null()
        real, dimension (:), pointer                :: BottomConc               => null()   !kg m-2
        real, dimension (:), pointer                :: MassInKg                 => null()   !kg (run with Benthos)
        real, dimension (:), pointer                :: Load                     => null()
        real, dimension (:), pointer                :: TotalConc                => null()   !**WASSIM 16/11/2005 
        real, dimension (:), pointer                :: ErosionRate              => null()   !kg m-2 s-1
        real, dimension (:), pointer                :: DepositionRate           => null()   !kg m-3 s-1
        real, dimension (:), pointer                :: Ws                       => null()   !m s-1 (vertical velocity)
                                                                                            !positive direction is downswards  
        real                                        :: MinValue     
        real                                        :: InitialValue 
        real                                        :: BottomMinConc                        !kg m-2
        real                                        :: BoundaryConcentration
        
        !Advection Diffusion
        real                                        :: Diffusivity 
        integer                                     :: Advection_Scheme
        integer                                     :: Diffusion_Scheme

        !Toxicity
        type (T_Toxicity)                           :: Toxicity
        
        !Decay
        real                                        :: DecayRate

        
        type (T_MassBalance)                        :: MB

        real                                        :: IScoefficient
        real                                        :: ExtinctionCoefficient
        real                                        :: ErosionCriticalShear        
        real                                        :: DepositionCriticalShear
        real                                        :: ErosionCoefficient
        real                                        :: CHS
        integer                                     :: Ws_Type
        real                                        :: Ws_Value
        real                                        :: KL
        real                                        :: KL1
        real                                        :: ML
        real                                        :: M

        character(PathLength)                       :: OutputName
        type (T_Property), pointer                  :: Next, Prev
      end type    T_Property

    type T_StormWaterModelLink
        integer                                     :: nOutflowNodes        = 0  !Nº of nodes where water flows from here to SWMM
        integer                                     :: nInflowNodes         = 0  !Nº of nodes where water flows from SWMM to here
        integer, dimension(:), allocatable          :: OutflowIDs           
        integer, dimension(:), allocatable          :: InflowIDs            
        real, dimension(:), allocatable             :: Outflow              
        real, dimension(:), allocatable             :: Inflow               
    end type T_StormWaterModelLink

    type T_DrainageNetwork
        integer                                     :: InstanceID            = 0
        character(len=StringLength)                 :: ModelName        
        integer                                     :: ObjEnterData          = 0  
        integer                                     :: ObjDischarges         = 0
        integer                                     :: ObjTime               = 0 
        integer                                     :: ObjInterface          = 0
        integer                                     :: ObjBenthicInterface   = 0
        integer                                     :: ObjLightExtinction    = 0
        integer                                     :: ObjHDF5               = 0
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime
        type (T_Time)                               :: CurrentTime
        type (T_Node) , dimension(:), pointer       :: Nodes                 => null()
        type (T_Reach), dimension(:), pointer       :: Reaches               => null()
        integer       , dimension(:), pointer       :: ComputeFaces          => null()
        integer       , dimension(:), pointer       :: OpenPointsFlow        => null()
        integer       , dimension(:), pointer       :: OpenPointsProcess     => null()
        integer       , dimension(:), pointer       :: RiverPoints           => null()        
        integer                                     :: TotalNodes            = 0        
        integer                                     :: TotalReaches          = 0
        integer                                     :: TotalOutlets          = 0
        integer                                     :: OutletReachPos
        integer                                     :: OutletNodePos
        integer                                     :: HighestOrder          = 0
        logical                                     :: CheckNodes
        logical                                     :: CheckReaches
        integer                                     :: XSCalc
        logical                                     :: HasGrid
        integer                                     :: CoordType
        type (T_OutPut)                             :: OutPut
        type (T_ComputeOptions)                     :: ComputeOptions
        type (T_TimeSerie)                          :: TimeSerie
        type (T_Files )                             :: Files
        type (T_Coupled  )                          :: Coupled
        type (T_StormWaterModelLink)                :: StormWaterModelLink
        logical                                     :: Continuous
        logical                                     :: StopOnWrongDate
        type (T_Property), pointer                  :: FirstProperty         => null()
        type (T_Property), pointer                  :: LastProperty          => null()
        integer                                     :: PropertiesNumber
        logical                                     :: HasProperties         = .false.
        
        
        real                                        :: GlobalManning
        logical                                     :: AllowBackwardWater    = .false.
        real                                        :: MinimumSlope
        real                                        :: InitialWaterDepth
        real                                        :: MinimumWaterDepth
        real                                        :: MinimumWaterDepthProcess
        
        integer                                     :: HydrodynamicApproximation
        real                                        :: NumericalScheme
        real,    dimension(:)  , pointer            :: RunOffVector         => null()
        real,    dimension(:)  , pointer            :: GroundVector         => null()
        real,    dimension(:,:,:), pointer          :: GroundVectorLayers   => null()
        real,    dimension(:)  , pointer            :: DiffuseVector        => null()
        real,    dimension(:)  , pointer            :: TransmissionFlow     => null()
        
        logical                                     :: GWFlowByLayers
        integer, dimension(:), pointer              :: GWFlowBottomLayer
        integer, dimension(:), pointer              :: GWFlowTopLayer

        real,    dimension(:,:), pointer            :: ChannelsWaterLevel   => null()
        real,    dimension(:,:), pointer            :: ChannelsBottomLevel  => null()
        real,    dimension(:,:), pointer            :: ChannelsBottomWidth  => null()
        real,    dimension(:,:), pointer            :: ChannelsSurfaceWidth => null()
        real,    dimension(:,:), pointer            :: ChannelsBankSlope    => null()
        real,    dimension(:,:), pointer            :: ChannelsNodeLength   => null()
        integer, dimension(:,:), pointer            :: ChannelsOpenProcess  => null()               
        real,    dimension(:)  , pointer            :: ShortWaveExtinction  => null()
        real,    dimension(:)  , pointer            :: ShortWaveField       => null()
        real,    dimension(:)  , pointer            :: LongWaveField        => null()
        real,    dimension(:)  , pointer            :: NodesDWZ             => null()
        real,    dimension(:)  , pointer            :: TopRadiation         => null()
        real,    dimension(:)  , pointer            :: AirTemperature       => null()
        real,    dimension(:)  , pointer            :: CloudCover           => null()
        real,    dimension(:)  , pointer            :: RelativeHumidity     => null()
        real,    dimension(:)  , pointer            :: WindSpeed            => null()
        real,    dimension(:)  , pointer            :: SedimentTemperature  => null()
        integer, dimension(:)  , pointer            :: DischargesLink       => null()
        real,    dimension(:)  , pointer            :: DischargesFlow       => null()
        real,    dimension(:,:), pointer            :: DischargesConc       => null()

        integer                    , dimension(:,:), pointer :: ChannelsID          => null()   
        character(len=StringLength), dimension(:,:), pointer :: ChannelsStationName => null()   

        logical                                     :: Discharges           = OFF   

        type (T_Downstream)                         :: Downstream
        type (T_Size2D)                             :: Size
        type (T_ExtVar)                             :: ExtVar        
        real                                        :: NextDT               = null_real
        integer                                     :: LastGoodNiter        = 1
        integer                                     :: InternalTimeStepSplit= 5

        real, dimension (:), pointer                :: GlobalToxicity       => null()
        integer                                     :: nToxicProp           = 0
        character(len=StringLength)                 :: GlobalToxicityEvolution
       

        !MassBalance
        logical                                     :: CheckMass 
        real(8)                                     :: TotalStoredVolume            = 0.0
        real(8)                                     :: TotalOutputVolume            = 0.0
        real(8)                                     :: TotalFlowVolume              = 0.0 !TotalOutput - EVTP
        real(8)                                     :: TotalInputVolume             = 0.0 !by discharges
        real(8)                                     :: TotalOverTopVolume           = 0.0 !OverTopping
        real(8)                                     :: TotalStormWaterOutput        = 0.0 !Total outflow to the Storm Water System
        real(8)                                     :: TotalStormWaterInput         = 0.0
        real(8)                                     :: InitialTotalOutputVolume     = 0.0
        real(8)                                     :: InitialTotalFlowVolume       = 0.0
        real(8)                                     :: InitialTotalInputVolume      = 0.0 !by discharges

        logical                                     :: Stabilize                    = .true.
        real                                        :: StabilizeFactor
        real                                        :: StabilizeCoefficient
        integer                                     :: MaxIterations
        real                                        :: DTFactor
        real                                        :: MaxDTFlood

        integer                                     :: nPropWithDischarges          = 0   !Performance
        !T90    
        integer                                     :: T90Var_Method                = null_int
        real                                        :: T90

        !Ripirian Shading
        real                                        :: ShadingFactor

        !Transmission Losses
        real                                        :: HydraulicConductivity                 

        integer                                     :: AerationEquation
        
      
        real,    dimension(:)  , pointer            :: ShearStress                  => null()        

        logical                                     :: WriteMaxStationValues        = .false.

        logical                                     :: OutputHydro                  = .false.
        
        !Evapotranspirate in reach pools
        real                                        :: EVTPMaximumDepth
        real                                        :: EVTPCropCoefficient

        type (T_DrainageNetwork), pointer           :: Next                         => null()
    end type  T_DrainageNetwork


    !Global Module Variables
    type (T_DrainageNetwork), pointer               :: FirstDrainageNetwork => null()
    type (T_DrainageNetwork), pointer               :: Me                   => null()

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !----------------------------------------------------------------------------

    subroutine ConstructDrainageNetwork(ModelName, DrainageNetworkID, TimeID, Size,   &
                                        CheckMass, CoupledPMP, CoupledRP, STAT)

        !Arguments---------------------------------------------------------------
        character(len=*)                                :: ModelName
        integer                                         :: DrainageNetworkID  
        integer                                         :: TimeID
        type (T_Size2D), optional                       :: Size
        logical, optional                               :: CheckMass
        logical, optional                               :: CoupledPMP, CoupledRP     
        integer, optional, intent(OUT)                  :: STAT     

        !Local-------------------------------------------------------------------
        integer                                         :: nDischarges, iDis, NodePos, NodeID
        logical                                         :: Found
        integer                                         :: STAT_CALL
        integer                                         :: ready_         
        type (T_Node), pointer                          :: CurrNode
        type(T_Property), pointer                       :: Property
        real                                            :: BottomMass

        !------------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mDRAINAGENETWORK_)) then
            nullify (FirstDrainageNetwork)
            call RegisterModule (mDrainageNetwork_) 
        endif
 
        call Ready(DrainageNetworkID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
         
            call AllocateInstance

            Me%ModelName        = ModelName

            !Associates module Time
            Me%ObjTime = AssociateInstance   (mTIME_, TimeID)
            
            if (present(CoupledPMP)) then
                Me%ExtVar%CoupledPMP = CoupledPMP
            endif
            if (present(CoupledRP)) then
                Me%ExtVar%CoupledRP  = CoupledRP
            endif

            !Gets Current Compute Time
            call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR01'
            
            !Gets Compute Time Limits
            call GetComputeTimeLimits (Me%ObjTime, BeginTime = Me%BeginTime,           &
                                       EndTime = Me%EndTime, STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR01a'
            

            !Verifies if Drainage Network Runs coupled to a grid or not
            if (present(Size)) then
                Me%HasGrid = .true.
                Me%Size    = Size
            else
                Me%HasGrid = .false.
            endif

            if (present(CheckMass)) then
                Me%CheckMass = CheckMass
            else
                Me%CheckMass = .false.
            end if
                       
            !Reads main user options
            call ReadDataFile

            call ConstructDownstreamBoundary

            !Connects nodes / reaches
            call ConstructNetwork

            !Set up properties to be transported
            call ConstructPropertyList

            !Verifies Global consistence of properties
            call CheckSelectedProp
            
            !Initial all variables
            call InitializeVariables

            if (Me%ComputeOptions%Discharges) then
                call Construct_Discharges(Me%ObjDischarges,                              &
                                          Me%ObjTime,                                    &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR02'  
                
                !Build Discharge NodeID / NodePos link
                !Gets the number of discharges
                call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR03'

                allocate(Me%DischargesLink(nDischarges))
                allocate(Me%DischargesFlow(nDischarges))
                allocate(Me%DischargesConc(nDischarges, Me%nPropWithDischarges))

                do iDis = 1, nDischarges

                    call GetDischargesNodeID  (Me%ObjDischarges, iDis, NodeID, STAT = STAT_CALL)
                    if (STAT_CALL/=SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR04'

                    call FindNodePosition   (NodeID, NodePos, Found)

                    CurrNode => Me%Nodes(NodePos)
                    CurrNode%Discharges = .true.
                    
                    if (Found) then
                        Me%DischargesLink(iDis) = NodePos
                    else
                        write (*,*) 'Discharge Node not found'
                        write (*,*) 'Node ID = ', NodeID            
                        stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR05'
                    end if
                    
                end do
                
            endif        
            
            !Link to StormWaterModel
            if (Me%ComputeOptions%StormWaterModelLink) then
                call ConstructStormWaterModelLink
            endif
            
            !Couples other modules
            call ConstructSubModules

            !Opens Output files
            call ConstructOutput
            
            !First HDF Output
            call HDF5Output

            !User Feed-Back
            call ConstructLog
            
            if (Me%CheckMass) then
                Me%TotalStoredVolume = 0.0
                Property => Me%FirstProperty
                do while (associated(Property))
                    Property%MB%TotalStoredMass    = 0.0
                    Property => Property%Next
                enddo
                
                do NodeID = 1, Me%TotalNodes
                    if (Me%Nodes(NodeID)%nDownStreamReaches /= 0) then
                        Me%TotalStoredVolume = Me%TotalStoredVolume + Me%Nodes(NodeID)%VolumeNew
                    endif
                
                    Property => Me%FirstProperty
                    do while (associated(Property))
                        
                        CurrNode => Me%Nodes(NodeID)
                        BottomMass = 0.0
                        if (Check_Particulate_Property(Property%ID%IDNumber)) then
                            ![kg] = [kg/m2] * [m2]
                            BottomMass = Property%BottomConc(NodeID) * CurrNode%CrossSection%BottomWidth * CurrNode%Length
                        endif
                        
                        ![kg] = [kg] + [kg] + [g/m3] * [m3] * [1e-3kg/g]
                        Property%MB%TotalStoredMass = Property%MB%TotalStoredMass  +  BottomMass  &
                                                      + Property%Concentration (NodeID)           &
                                                      * Property%ISCoefficient                    &
                                                      * Me%Nodes(NodeID)%VolumeNew
                        
                        Property => Property%Next
                    enddo
                enddo                
                
            end if
            

            !Close input data file
            call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR03'

            !Returns ID
            DrainageNetworkID = Me%InstanceID

            STAT_CALL = SUCCESS_
            

        else cd0
            
            stop 'ModuleDrainageNetwork - ConstructDrainageNetwork - ERR04'

        end if cd0

        if (present(STAT)) STAT = STAT_CALL

        !-----------------------------------------------------------------------

    end subroutine ConstructDrainageNetwork
 
    !---------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments--------------------------------------------------------------
                                                    
        !Local------------------------------------------------------------------
        type (T_DrainageNetwork), pointer           :: NewObjDrainageNetwork
        type (T_DrainageNetwork), pointer           :: PreviousObjDrainageNetwork


        !Allocates new instance
        allocate (NewObjDrainageNetwork)
        nullify  (NewObjDrainageNetwork%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstDrainageNetwork)) then
            FirstDrainageNetwork            => NewObjDrainageNetwork
            Me                              => NewObjDrainageNetwork
        else
            PreviousObjDrainageNetwork      => FirstDrainageNetwork
            Me                              => FirstDrainageNetwork%Next
            do while (associated(Me))
                PreviousObjDrainageNetwork  => Me
                Me                          => Me%Next
            enddo
            Me                              => NewObjDrainageNetwork
            PreviousObjDrainageNetwork%Next => NewObjDrainageNetwork
        endif

        Me%InstanceID = RegisterNewInstance (mDrainageNetwork_)


    end subroutine AllocateInstance

    !---------------------------------------------------------------------------

    subroutine ReadDataFile

        !Local------------------------------------------------------------------
        integer                                     :: flag, STAT_CALL
        integer                                     :: GeoConversationFactor
        character(len=StringLength)                 :: AuxString

        !Reads name of the data file from nomfich.dat
        call ReadFileName('DRAINAGE_NETWORK', Me%Files%InputData,               &
                           Message = "Drainage Network Data File",              &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR01' 

        call ReadFileName('DRAINAGE_NETWORK_FIN', Me%Files%FinalFile,           &
                           Message = "Drainage Network Final File",             &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR02' 

        call ReadFileName('DRAINAGE_NETWORK_HDF', Me%Files%HDFFile,             &
                           Message = "Drainage Network HDF File",               &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR03' 
        
        call ConstructEnterData (Me%ObjEnterData, Me%Files%InputData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR04' 

        call GetData(Me%Files%Network,                                          &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'NETWORK_FILE',                             &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR05'        
        
        call GetData(Me%CheckNodes,                                             &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'CHECK_NODES',                              &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = .true.,                                     &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR06'        
        

        call GetData(Me%CheckReaches,                                           &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'CHECK_REACHES',                            &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = .true.,                                     &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR07'        

        call GetData(Me%ComputeOptions%Discharges,                              &
                     Me%ObjEnterData, flag,                                     &
                     keyword    = 'DISCHARGES',                                 & 
                     default    = .false.,                                      &
                     SearchType = FromFile,                                     &                     
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR08'        

        call GetData(Me%HydrodynamicApproximation,                              &
                      Me%ObjEnterData, flag,                                    &  
                      keyword      = 'HYDRODYNAMIC_APROX',                      &
                      ClientModule = 'DrainageNetwork',                         &
                      SearchType   = FromFile,                                  &
                      Default      = KinematicWave,                            &
                      STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR09'        

        call GetData(Me%NumericalScheme,                                        &
                      Me%ObjEnterData, flag,                                    &  
                      keyword      = 'NUMERICAL_SCHEME',                        &
                      ClientModule = 'DrainageNetwork',                         &
                      SearchType   = FromFile,                                  &
                      Default      = ExplicitScheme,                            &
                      STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR010'        

        if (Me%NumericalScheme /= ExplicitScheme .and. Me%NumericalScheme /= ImplicitScheme) &
        stop 'ModuleDrainageNetwork - ReadDataFile - ERR09b'        

        call GetData(Me%GlobalManning,                                          &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'GLOBAL_MANNING',                           &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = null_real,                                  &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR11'        

        call GetData(Me%AllowBackwardWater,                                     &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'ALLOW_BACKWATER',                          &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = .false.,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR11a'        
            

        call GetData(Me%MinimumSlope,                                           &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'MINIMUM_SLOPE',                            &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = 0.0,                                        &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR12'        


        call GetData(Me%MinimumWaterDepth,                                      &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'MIN_WATER_DEPTH',                          &   
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = 0.001,                                      &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR13'        

        if (Me%MinimumWaterDepth.LT.0.0) then
            write (*,*)'Invalid Number of Minimum Water Level [MIN_WATER_DEPTH]'
            stop 'ModuleDrainageNetwork - ReadDataFile - ERR14'
        end if

        call GetData(Me%MinimumWaterDepthProcess,                               &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'MIN_WATER_DEPTH_PROCESS',                  &   
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = 0.01,                                       &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR15'        

        if (Me%MinimumWaterDepthProcess.LT.0.0) then
            write (*,*)'Invalid Number of Minimum Water Level [MIN_WATER_DEPTH_PROCESS]'
            stop 'ModuleDrainageNetwork - ReadDataFile - ERR16'
        end if


        call GetData(Me%Continuous,                                             &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'CONTINUOUS',                               &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = OFF,                                        &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR17'        

        if (Me%Continuous) then
            call ReadFileName('DRAINAGE_NETWORK_INI', Me%Files%Initial,         &
                              Message = "Drainage Network Initial File",        &
                              STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR18'

            call GetData(Me%StopOnWrongDate,                                        &
                         Me%ObjEnterData, flag,                                     &  
                         keyword      = 'STOP_ON_WRONG_DATE',                       &
                         ClientModule = 'DrainageNetwork',                          &
                         SearchType   = FromFile,                                   &
                         Default      = .true.,                                     &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR19'        


        else        

            call GetData(Me%InitialWaterDepth,                                      &
                         Me%ObjEnterData, flag,                                     &  
                         keyword      = 'INITIAL_WATER_DEPTH',                      &
                         ClientModule = 'DrainageNetwork',                          &
                         SearchType   = FromFile,                                   &
                         Default      = 0.0,                                        &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR20'        

            if (Me%InitialWaterDepth.LT.0.0) then
                write (*,*)'Invalid Number of Initial Water Level [INITIAL_WATER_DEPTH]'
                stop 'ModuleDrainageNetwork - ReadDataFile - ERR21'
            end if

        end if

        call GetData(Me%Stabilize,                                              &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'STABILIZE',                                &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = .true.,                                     &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR22'        

        if (Me%Stabilize) then

            call GetData(Me%StabilizeFactor,                                    &
                         Me%ObjEnterData, flag,                                 &  
                         keyword      = 'STABILIZE_FACTOR',                     &
                         ClientModule = 'DrainageNetwork',                      &
                         SearchType   = FromFile,                               &
                         Default      = 0.1,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR23'        

            call GetData(Me%StabilizeCoefficient,                               &
                         Me%ObjEnterData, flag,                                 &  
                         keyword      = 'STABILIZE_COEFFICIENT',                &
                         ClientModule = 'DrainageNetwork',                      &
                         SearchType   = FromFile,                               &
                         Default      = 0.05,                                   &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR23'        

            call GetData(Me%MaxIterations,                                      &
                         Me%ObjEnterData, flag,                                 &  
                         keyword      = 'MAX_ITERATIONS',                       &
                         ClientModule = 'DrainageNetwork',                      &
                         SearchType   = FromFile,                               &
                         Default      = 100,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR24'        

        end if

        !Factor for DT Prediction
        call GetData(Me%DTFactor,                                           &
                     Me%ObjEnterData, flag,                                 &  
                     keyword      = 'DT_FACTOR',                            &
                     ClientModule = 'DrainageNetwork',                      &
                     SearchType   = FromFile,                               &
                     Default      = 1.05,                                   &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR25'        

        if (Me%DTFactor <= 1.0) then
            write (*,*)'Invalid DT Factor [DT_FACTOR]'
            write (*,*)'Value must be greater then 1.0'
            stop 'ModuleDrainageNetwork - ReadDataFile - ERR28a'              
        endif

        !Max DT if channel water level exceeds full bank 
        call GetData(Me%MaxDTFlood,                                         &
                     Me%ObjEnterData, flag,                                 &  
                     keyword      = 'MAX_DT_FLOOD',                         &
                     ClientModule = 'DrainageNetwork',                      &
                     SearchType   = FromFile,                               &
                     Default      = 10.0,                                   &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR26'        

        call GetData(Me%AerationEquation,                                   &
                     Me%ObjEnterData, flag,                                 &
                     Keyword      ='AERATION_METHOD',                       &
                     SearchType   = FromFile,                               &
                     ClientModule = 'DrainageNetwork',                      &
                     Default      = PoolAndRifle_,                          &
                    STAT          = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR27'  

        if (Me%AerationEquation /= PoolAndRifle_ .and. Me%AerationEquation /= ChannelControled_) then
            write (*,*)'Invalid O2 Aeration Method'
            stop 'ModuleDrainageNetwork - ReadDataFile - ERR28'  
        endif


        call GetData(Me%T90Var_Method,                                      &   
                     Me%ObjEnterData, flag,                                 &
                     Keyword        = 'T90_DECAY_MODEL',                    &
                     ClientModule   = 'ModuleDrainageNetwork',              &
                     SearchType     = FromFile,                             &
                     Default        = Canteras,                             &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR29' 

        if (Me%T90Var_Method == Constant) then
            call GetData(Me%T90,                                            &   
                         Me%ObjEnterData, flag,                             &
                         Keyword        = 'T90',                            &
                         ClientModule   = 'ModuleDrainageNetwork',          &
                         SearchType     = FromFile,                         &
                         Default        = 7200.,                            &
                         STAT           = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR30' 
        endif

        call GetData(Me%ShadingFactor,                                      &
                    Me%ObjEnterData,                                        &
                    flag,                                                   &
                    SearchType   = FromFile,                                &
                    keyword      = 'SHADING_FACTOR',                        &
                    Default      = 1.0,                                     &
                    ClientModule = 'DrainageNetwork',                       &
                    STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR31' 


        call GetData(Me%ComputeOptions%TransmissionLosses,                  &
                     Me%ObjEnterData,                                       &
                     flag,                                                  &
                     SearchType   = FromFile,                               &
                     keyword      = 'TRANSMISSION_LOSSES',                  &
                     Default      = .false.,                                &
                     ClientModule = 'DrainageNetwork',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR32' 


        call GetData(Me%ComputeOptions%RemoveOverTop,                       &
                     Me%ObjEnterData,                                       &
                     flag,                                                  &
                     SearchType   = FromFile,                               &
                     keyword      = 'REMOVE_OVERTOP',                       &
                     Default      = .false.,                                &
                     ClientModule = 'DrainageNetwork',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR33' 

        call GetData(Me%ComputeOptions%CalcFractionSediment,                &
                     Me%ObjEnterData,                                       &
                     flag,                                                  &
                     SearchType   = FromFile,                               &
                     keyword      = 'FRACTION_SEDIMENT',                    &
                     Default      = .false.,                                &
                     ClientModule = 'DrainageNetwork',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR33_Wassim_FractionSediment' 



        if (Me%ComputeOptions%TransmissionLosses) then

            call GetData(Me%HydraulicConductivity,                          &
                         Me%ObjEnterData,                                   &
                         flag,                                              &
                         SearchType   = FromFile,                           &
                         keyword      = 'HYDRAULIC_CONDUCTIVITY',           &
                         Default      = 1.e-5,                              &
                         ClientModule = 'DrainageNetwork',                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR34' 

        endif

        !Reads Global Toxicity Computation Method
        call GetData(AuxString, Me%ObjEnterData,  flag,                     &
                     keyword        = 'GLOBAL_TOXICITY',                    &                         
                     ClientModule   = 'DrainageNetwork',                    &
                     SearchType     = FromFile,                             &
                     Default        = 'SUM',                                &                         
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR35'

        select case (trim(adjustl(AuxString)))
            case ("Max",       "MAX", "max")
                Me%GlobalToxicityEvolution = 'MAX'
            case ("Sum",       "SUM", "sum")
                Me%GlobalToxicityEvolution = 'SUM'            
            case ("Riskratio", "RiskRatio", "RISKRATIO", "riskratio")
                Me%GlobalToxicityEvolution = 'RISKRATIO'
            case default
                write(*,*)'Invalid option for keyword GLOBAL_TOXICITY'
                stop 'ModuleDrainageNetwork - ReadDataFile - ERR36' 
        end select


        !Reads Global GeoConversation Factor (Lat/ to Meters) rough estimation
        call GetData(GeoConversationFactor, Me%ObjEnterData,  flag,                      &
                     keyword        = 'GEO_CONVERSATION_FACTOR',                            &                         
                     ClientModule   = 'DrainageNetwork',                                    &
                     SearchType     = FromFile,                                             &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR37'
        
        if (flag == 1) then
            call SetError(WARNING_, INTERNAL_, 'The keyword GEO_CONVERSATION_FACTOR is obselete and not used any more', ON)
        endif


        !Output Hydrodynamic properties
        call GetData(Me%OutputHydro, Me%ObjEnterData, flag,                             &
                     keyword        = 'OUTPUT_HYDRO',                                   &                         
                     ClientModule   = 'DrainageNetwork',                                &
                     SearchType     = FromFile,                                         &
                     Default        = .FALSE.,                                          &                         
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR38' 


        !Sets Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                              &
                           CurrentTime = Me%CurrentTime,                                 &
                           EndTime     = Me%EndTime,                                     &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%Yes,                                  &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadDataFile - ERR40' 

        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%CurrentTime,                               &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleDrainageNetwork - ERR41'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleBasin',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR42'

        call GetData(Me%Output%ComputeFlowFrequency,                                    &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OUTPUT_FLOW_FREQUENCY',                            &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBasin',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR50'
        
        if (Me%Output%ComputeFlowFrequency) then
            !Reads Begin Time for frequency analisys
            call GetData(Me%Output%FlowFrequency%StartDate,                                 &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = FromFile,                                           &
                         keyword      = 'FLOW_FREQUENCY_STARTDATE',                         &
                         Default      = Me%BeginTime,                                       &
                         ClientModule = 'ModuleBasin',                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR55'            


            call GetData(Me%Output%FlowFrequency%StopDate,                                  &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = FromFile,                                           &
                         keyword      = 'FLOW_FREQUENCY_ENDDATE',                           &
                         Default      = Me%EndTime,                                         &
                         ClientModule = 'ModuleBasin',                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR60'            
            
            call GetData(Me%Output%FlowFrequency%MinimumFlow,                               &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = FromFile,                                           &
                         keyword      = 'FLOW_FREQUENCY_MINIMUMFLOW',                       &
                         Default      = 0.0,                                                &
                         ClientModule = 'ModuleBasin',                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR65'            

        
        endif
        
        !to evapotrnaspirate from reach - in drying pools where vegetation accumulates and removes water
        call GetData(Me%ComputeOptions%EVTPFromReach,                                   &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'EVTP_FROM_REACH',                                  &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBasin',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR70'
        
        if (Me%ComputeOptions%EVTPFromReach) then
            
            !maximum depth to happen evtp (vegetation only installs in low flow conditions)
            call GetData(Me%EVTPMaximumDepth,                                               &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = FromFile,                                           &
                         keyword      = 'EVTP_MAXIMUM_DEPTH',                               &
                         Default      = 0.1,                                                &
                         ClientModule = 'ModuleBasin',                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR80'


            !crop coefficient - multiply by potential evapotransp.
            call GetData(Me%EVTPCropCoefficient,                                            &
                         Me%ObjEnterData,                                                   &
                         flag,                                                              &
                         SearchType   = FromFile,                                           &
                         keyword      = 'EVTP_CROP_COEF',                                   &
                         Default      = 1.0,                                                &
                         ClientModule = 'ModuleBasin',                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR90'
            
        endif

        !If linked to a StormWaterModel
        call GetData(Me%ComputeOptions%StormWaterModelLink,                             &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'STORM_WATER_MODEL_LINK',                           &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleDrainageNetwork',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleDrainageNetwork - ERR100'


    end subroutine ReadDataFile
    
    !---------------------------------------------------------------------------

    subroutine ConstructDownstreamBoundary
        
        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                         :: flag, STAT_CALL
        character(len=StringLength)                     :: AuxString

        call GetData(Me%Downstream%Boundary,                                    &
                     Me%ObjEnterData, flag,                                     &  
                     keyword      = 'DOWNSTREAM_BOUNDARY',                      &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromFile,                                   &
                     Default      = ZeroDepthGradient,                          &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR01'

if1:    if (Me%Downstream%Boundary == ImposedWaterLevel .or. Me%Downstream%Boundary == ImposedVelocity) then

            !Reads Time Evolution
            call GetData(AuxString, Me%ObjEnterData,  flag,                     &
                         keyword        = 'FILE_IN_TIME',                       &                         
                         ClientModule   = 'DrainageNetwork',                    &
                         SearchType     = FromFile,                             &
                         Default        = 'None',                               &                         
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR03'

            select case (trim(adjustl(AuxString)))

                case ("None",       "NONE", "none")
                    
                    Me%Downstream%Evolution    = None

                case ("Timeserie",  "TIMESERIE",    "timeserie",    "TimeSerie")
                    
                    Me%Downstream%Evolution    = ReadTimeSerie
                    
                    if (Me%Downstream%Boundary == ImposedVelocity)    &
                        stop 'not ready - ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR04a'
                        
                case ("OpenMI", "OPENMI", "openmi", "OpenMi")
                
                    Me%Downstream%Evolution    = OpenMI

                case default

                    write(*,*)'Invalid option for keyword FILE_IN_TIME'
                    stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR04'

            end select

            call GetData(Me%Downstream%DefaultValue, Me%ObjEnterData,  flag,    &
                         keyword        = 'DEFAULT_VALUE',                       &                         
                         ClientModule   = 'DrainageNetwork',                    &
                         SearchType     = FromFile,                             &
                         Default        = FillValueReal,                        &                         
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR05'


            if (flag == 0) then
                write(*,*)'Please define default value for downstream boundary'
                stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR06'
            end if
           
if2:        if (Me%Downstream%Evolution == ReadTimeSerie) then
                            
                call GetData(Me%Downstream%FileName,                            &
                             Me%ObjEnterData , flag,                            &
                             SearchType   = FromFile,                           &
                             keyword      = 'FILENAME',                         &
                             ClientModule = 'DrainageNetwork',                  &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR07'

                if (flag==0)then
                    write(*,*)'Time Serie File Name not given'
                    stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR08'
                endif

                call GetData(Me%Downstream%DataColumn,                           &
                             Me%ObjEnterData , flag,                            &
                             SearchType   = FromFile,                           &
                             keyword      = 'DATA_COLUMN',                      &
                             ClientModule = 'DrainageNetwork',                  &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR09'

                if (flag==0)then
                    write(*,*)'Data Column not given'
                    stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR10'
                endif

                !Starts Time Serie
                call StartTimeSerieInput(Me%Downstream%ObjTimeSerie,            &
                                         Me%Downstream%FileName,                &
                                         Me%ObjTime,                            &
                                         STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructDownstreamBoundary - ERR11'
                
            end if if2


        endif if1

    end subroutine ConstructDownstreamBoundary
    
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ConstructNetwork


        !Local------------------------------------------------------------------
        integer                                     :: flag, STAT_CALL

        call ConstructEnterData (Me%Files%ObjEnterDataNetwork, Me%Files%Network, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNetwork - ERR01'

        !Checks for the COORD_TIP
        call GetData(Me%CoordType,                                                       &
                     Me%Files%ObjEnterDataNetwork, flag,                                 & 
                     keyword      = 'COORDINATE_TYPE',                                   &
                     ClientModule = 'DrainageNetwork',                                   &
                     SearchType   = FromFile,                                            &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNetwork - ERR02'

        if (flag == 0 .or. (Me%CoordType /= 1 .and. Me%CoordType /= 2)) then
            write(*,*)'The Drainage Network does not contain a valid specification for the coordinate system.'
            write(*,*)'Please set the keyword COORD_TIP to a valid option'
            write(*,*)'Allowed options are:'
            write(*,*)'COORDINATE_TYPE      : 1 ! Geographic Coordinates'
            write(*,*)'COORDINATE_TYPE      : 2 ! Projected Coordinates'
            call SetError (FATAL_, INTERNAL_, "Invalid Coordinates")
        endif
        
        !Rewinds buffer for subsequent readings
        call RewindBuffer(Me%Files%ObjEnterDataNetwork, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNetwork - ERR03'

        call ConstructNodeList

        call ConstructReachList 

        call ConnectNetwork

        if (Me%NumericalScheme == ImplicitScheme) then
            call OrderNodes
            call ReconnectNetwork
            call WriteOrderedNodes
        end if
        
        !Checks consistency and finds outlet Node / Reach Position
        call CountOutlets ()

                        
                        
        call KillEnterData (Me%Files%ObjEnterDataNetwork, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNetwork - ERR04'

    end subroutine ConstructNetwork

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ConstructNodeList

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine            
        integer                                     :: STAT_CALL, NodePos
               
        !-----------------------------------------------------------------------

        call CountTotalNodes
 
        nullify  (Me%Nodes)
        allocate (Me%Nodes (1:Me%TotalNodes))
        
        NodePos = 0

do1:    do 

            call ExtractBlockFromBuffer(Me%Files%ObjEnterDataNetwork, ClientNumber,     &
                                        BeginNode, EndNode, BlockFound,                  &
                                        FirstLine, LastLine, STAT_CALL) 

if1:        if (STAT_CALL .EQ. SUCCESS_) then    

if2:              if (BlockFound) then                 
                    
                    NodePos = NodePos + 1

                    call ConstructNode (NodePos)
                                      
                  else if2

                    if (NodePos /= Me%TotalNodes) stop 'ModuleDrainageNetwork - ConstructNodeList - ERR01'

                    call Block_Unlock(Me%Files%ObjEnterDataNetwork, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'ModuleDrainageNetwork - ConstructNodeList - ERR02'

                    exit do1    !No more blocks

                  end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1

                    stop 'ModuleDrainageNetwork - ConstructNodeList - ERR02.'

            end if if1

        end do do1

        if (Me%CheckNodes) call CheckNodesConsistency

    end subroutine ConstructNodeList 

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine CountTotalNodes

        !This subroutine counts the total number of nodes and checks the
        !existence of valid and repeated NodeIDs
        !Local------------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine            
        integer                                     :: STAT_CALL
        integer                                     :: NodeID, OldNodeID
        !integer                                     :: MaxNodeID, MinNodeID 
        integer                                     :: flag

        
        Me%TotalNodes = 0
        !MinNodeID = - null_int
        !MaxNodeID = null_int
        OldNodeID = null_int


do1:    do 

            call ExtractBlockFromBuffer(Me%Files%ObjEnterDataNetwork, ClientNumber,     &
                                        BeginNode, EndNode, BlockFound,                 &
                                        FirstLine, LastLine, STAT_CALL) 

if1:        if (STAT_CALL .EQ. SUCCESS_) then    

if2:            if (BlockFound) then                 
                    
                    !Gets ID
                    call GetData(NodeID,                                        &
                                 Me%Files%ObjEnterDataNetwork, flag,            & 
                                 keyword      = 'ID',                           &
                                 ClientModule = 'DrainageNetwork',              &
                                 SearchType   = FromBlock,                      &
                                 STAT         = STAT_CALL)                                  
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - CountTotalNodes - ERR01'

                    if (flag /= 1) then
                        write (*,*)'Invalid Node ID [ID]'
                        stop 'ModuleDrainageNetwork - CountTotalNodes - ERR02'
                    endif

                    Me%TotalNodes = Me%TotalNodes + 1
                                        
                    !if (NodeID .LT. MinNodeID ) MinNodeID = NodeID
                    !if (NodeID .GT. MaxNodeID ) MaxNodeID = NodeID
                   
                    if (NodeID .EQ. OldNodeID ) then
                        write (*,*) 'Repeated Node ID = ', NodeID
                        stop 'ModuleDrainageNetwork - CountTotalNodes - ERR03'
                    else
                        OldNodeID = NodeID
                    end if

                else if2
                    
                    !if (MinNodeID.NE. 1) then
                    !    write (*,*) 'Inconsistency in Node IDs - Missing NodeID = 1'
                    !    stop 'ModuleDrainageNetwork - CountTotalNodes - ERR04'
                    !else if (MaxNodeID.NE. Me%TotalNodes) then
                    !    write (*,*) 'Inconsistency in Node IDs - Missing NodeID =', Me%TotalNodes
                    !    stop 'ModuleDrainageNetwork - CountTotalNodes - ERR05'
                    !end if


                    call Block_Unlock(Me%Files%ObjEnterDataNetwork, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - CountTotalNodes - ERR01'

                    call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - CountTotalNodes - ERR02'
                    
                    exit do1    !No more blocks

                  end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1

                    stop 'ModuleDrainageNetwork - CountTotalNodes - ERR03.'

            end if if1

        end do do1
        
    end subroutine CountTotalNodes
       
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ConstructNode (NodePos)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                         :: NodePos
        !External---------------------------------------------------------------
        type (T_Node), pointer                      :: NewNode
        integer                                     :: STAT_CALL
        integer                                     :: flag, NStations
        real, dimension (2)                         :: AuxCoord 
        logical                                     :: ComputeElevation

        !Local------------------------------------------------------------------
     
        nullify (NewNode)
        NewNode => Me%Nodes (NodePos)

        call GetData(NewNode%ID,                                                &
                     Me%Files%ObjEnterDataNetwork, flag,                        & 
                     keyword      = 'ID',                                       &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR01'

        if (flag /= 1) then
            write (*,*)'Invalid Node ID [ID]'
            stop 'ModuleDrainageNetwork - ConstructNode - ERR02'
        endif

        call GetData(NewNode%StationName,                                       &
                     Me%Files%ObjEnterDataNetwork, flag,                        & 
                     keyword        = 'ASSOCIATEDSTATION_NAME',                &
                     default        = null_str,                                 &
                     ClientModule   = 'DrainageNetwork',                        &
                     SearchType     = FromBlock,                                &
                     STAT           = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR01a'

        if (NewNode%StationName /= null_str) Me%WriteMaxStationValues = .TRUE.
        !Gets Location        
        if (NewNode%X.EQ.null_real.AND. NewNode%Y.EQ.null_real) then
            call GetData(AuxCoord,                                              &
                         Me%Files%ObjEnterDataNetwork, flag,                    &
                         keyword      = 'COORDINATES',                          &
                         ClientModule = 'DrainageNetwork',                      &
                         SearchType   = FromBlock,                              &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR03'
                                                                                
            if (flag .EQ. 2) then
                NewNode%X = AuxCoord(1)
                NewNode%Y = AuxCoord(2)                       
            else
                write(*,*) 'Invalid Node Coordenates [COORDINATES]'
                stop 'ModuleDrainageNetwork - ConstructNode - ERR04'
            end if 

        else
            write (*,*) 'Repeated Node = ', NewNode%ID
            stop 'ModuleDrainageNetwork - ConstructNode - ERR05'
        end if

        !Gets associated Grid Point I
        call GetData(NewNode%GridI,                                             &
                     Me%Files%ObjEnterDataNetwork, flag,                        &
                     keyword      = 'GRID_I',                                   &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     default      = null_int,                                   &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR06'
                                                                                
        !Gets associated Grid Point J
        call GetData(NewNode%GridJ,                                             &
                     Me%Files%ObjEnterDataNetwork, flag,                        &
                     keyword      = 'GRID_J',                                   &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     default      = null_int,                                   &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR07'

        !TerrainLevel (before InitializeTabularCrossSection)
        call GetData(NewNode%CrossSection%TerrainLevel,                         &
                     Me%Files%ObjEnterDataNetwork, flag,                        &  
                     keyword      = 'TERRAIN_LEVEL',                            &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR08'
        if (flag /= 1) then
            write (*,*)'Invalid Node Terrain Level [TERRAIN_LEVEL]'
            stop 'ModuleDrainageNetwork - ConstructNode - ERR21'
        endif

        !Singularity Coef - % available vertical area = 1 - % reduction Av by singularity
        call GetData(NewNode%SingCoef,                                          &
                     Me%Files%ObjEnterDataNetwork, flag,                        &  
                     keyword      = 'SING_COEF',                                &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     default      = 1.0,                                        &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR08a'

        if (NewNode%SingCoef <= AlmostZero) then 
            write (*,*)'Invalid Singularity Coefficient [SING_COEF]'
            stop 'ModuleDrainageNetwork - ConstructNode - ERR22'
        endif

        !Cross Section Type        
        call GetData(NewNode%CrossSection%Form,                                 &
                     Me%Files%ObjEnterDataNetwork, flag,                        &  
                     keyword      = 'CROSS_SECTION_TYPE',                       &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     default      = Trapezoidal,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR09'   


ifXS:   if (NewNode%CrossSection%Form == Trapezoidal .or.                       &
            NewNode%CrossSection%Form == TrapezoidalFlood) then

            !Bottom Width
            call GetData(NewNode%CrossSection%BottomWidth,                      &
                         Me%Files%ObjEnterDataNetwork, flag,                    &  
                         keyword      = 'BOTTOM_WIDTH',                         &
                         ClientModule = 'DrainageNetwork',                      &
                         SearchType   = FromBlock,                              &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR10'
            if (flag /= 1) then
                write (*,*)'Invalid Node Bottom Width [BOTTOM_WIDTH]'
                stop 'ModuleDrainageNetwork - ConstructNode - ERR10a'
            endif

            if ( NewNode%CrossSection%BottomWidth == 0.0) NewNode%CrossSection%BottomWidth = AllmostZero

                 
            !Top Width
            call GetData(NewNode%CrossSection%TopWidth,                         &
                         Me%Files%ObjEnterDataNetwork, flag,                    &  
                         keyword      = 'TOP_WIDTH',                            &
                         ClientModule = 'DrainageNetwork',                      &
                         SearchType   = FromBlock,                              &
                         Default      = null_real,                              &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR11'
            if (flag /= 1 ) then
                write (*,*)'Invalid Node Top Width [TOP_WIDTH]'
                stop 'ModuleDrainageNetwork - ConstructNode - ERR11a'
            endif


            !Height
            call GetData(NewNode%CrossSection%Height,                           &
                         Me%Files%ObjEnterDataNetwork, flag,                    &  
                         keyword      = 'HEIGHT',                               &
                         ClientModule = 'DrainageNetwork',                      &
                         SearchType   = FromBlock,                              &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR12'
            if (flag /= 1) then
                write (*,*)'Invalid Node Height [HEIGHT]'
                stop 'ModuleDrainageNetwork - ConstructNode - ERR12a'
            endif

            NewNode%CrossSection%Slope = (( NewNode%CrossSection%TopWidth        &
                                          - NewNode%CrossSection%BottomWidth )   &
                                          / 2 ) / NewNode%CrossSection%Height

            NewNode%CrossSection%BottomLevel = NewNode%CrossSection%TerrainLevel - NewNode%CrossSection%Height


            if (NewNode%CrossSection%Form == TrapezoidalFlood) then

                !MiddleWidth
                call GetData(NewNode%CrossSection%MiddleWidth,                      &
                             Me%Files%ObjEnterDataNetwork, flag,                    &  
                             keyword      = 'MIDDLE_WIDTH',                         &
                             ClientModule = 'DrainageNetwork',                      &
                             SearchType   = FromBlock,                              &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR13'
                if (flag /= 1) then
                    write (*,*)'Invalid Node Middle Width [MIDDLE_WIDTH]'
                    stop 'ModuleDrainageNetwork - ConstructNode - ERR13a'
                endif


                !MiddleHeight
                call GetData(NewNode%CrossSection%MiddleHeight,                     &
                             Me%Files%ObjEnterDataNetwork, flag,                    &  
                             keyword      = 'MIDDLE_HEIGHT',                        &
                             ClientModule = 'DrainageNetwork',                      &
                             SearchType   = FromBlock,                              &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR14'
                if (flag /= 1) then
                    write (*,*)'Invalid Node Middle Height [MIDDLE_HEIGHT]'
                    stop 'ModuleDrainageNetwork - ConstructNode - ERR14a'
                endif

                if (NewNode%CrossSection%MiddleHeight >= NewNode%CrossSection%Height) then
                    write (*,*)'Node Middle Height must be <= than Height'
                    stop 'ModuleDrainageNetwork - ConstructNode - ERR14b'

                endif

                NewNode%CrossSection%Slope = (( NewNode%CrossSection%MiddleWidth     &
                                              - NewNode%CrossSection%BottomWidth )   &
                                              / 2 ) / NewNode%CrossSection%MiddleHeight

                NewNode%CrossSection%SlopeTop = (( NewNode%CrossSection%TopWidth     &
                                              - NewNode%CrossSection%MiddleWidth )   &
                                              / 2 )                                  &
                                              / (NewNode%CrossSection%Height -       &
                                                 NewNode%CrossSection%MiddleHeight)

            endif


        elseif (NewNode%CrossSection%Form == Tabular) then !ifXS


            call GetData(NewNode%CrossSection%NStations,                            &
                         Me%Files%ObjEnterDataNetwork, flag,                        &
                         keyword      = 'N_STATIONS',                               &
                         ClientModule = 'DrainageNetwork',                          &
                         SearchType   = FromBlock,                                  &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR15'
            if (flag /= 1 .or. NewNode%CrossSection%NStations <= 3) then
                write (*,*)'Minimum mumber of station points for Tabular Cross Section is 3 [N_STATIONS].'
                write (*,*)'in node ', NewNode%ID
                stop 'ModuleDrainageNetwork - ConstructNode - ERR15a'
            endif
        
            NStations = NewNode%CrossSection%NStations

            allocate(NewNode%CrossSection%Station      (NStations))
            allocate(NewNode%CrossSection%Elevation    (NStations))
            allocate(NewNode%CrossSection%BankSlope    (NStations))

            NewNode%CrossSection%Station      = null_real
            NewNode%CrossSection%Elevation    = null_real
            NewNode%CrossSection%BankSlope    = null_real

            call GetData(NewNode%CrossSection%Station,                              &
                         Me%Files%ObjEnterDataNetwork, flag,                        &
                         keyword      = 'STATION',                                  &
                         ClientModule = 'DrainageNetwork',                          &
                         SearchType   = FromBlock,                                  &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR16'

            if (flag /= NStations) then
                write(*,*) 'Invalid Node Station data [STATION]'
                stop 'ModuleDrainageNetwork - ConstructNode - ERR16a'
            end if 

            call GetData(NewNode%CrossSection%Elevation,                            &
                         Me%Files%ObjEnterDataNetwork, flag,                        &
                         keyword      = 'ELEVATION',                                &
                         ClientModule = 'DrainageNetwork',                          &
                         SearchType   = FromBlock,                                  &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR17'

            if (flag == 0) then
                
                ComputeElevation = .true.

                call GetData(NewNode%CrossSection%Elevation,                            &
                             Me%Files%ObjEnterDataNetwork, flag,                        &
                             keyword      = 'LEVEL',                                    &
                             ClientModule = 'DrainageNetwork',                          &
                             SearchType   = FromBlock,                                  &
                             STAT         = STAT_CALL)                                  
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR17'

                if (flag /= NStations) then
                    write(*,*) 'Invalid Node LEVEL data [LEVEL]'
                    stop 'ModuleDrainageNetwork - ConstructNode - ERR17a'
                end if 

            else
            
                ComputeElevation = .false.

                if (flag /= NStations) then
                    write(*,*) 'Invalid Node Elevation data [ELEVATION]'
                    stop 'ModuleDrainageNetwork - ConstructNode - ERR17a'
                end if 

            endif

            call InitializeTabularCrossSection(NewNode, ComputeElevation)                                   

        else !ifXS
                    
            write (*,*)'Invalid Cross Section Form'
            stop 'ModuleDrainageNetwork - ConstructNode - ERR60'

        end if ifXS
           
        !Pool Depth
        call GetData(NewNode%CrossSection%PoolDepth,                        &
                     Me%Files%ObjEnterDataNetwork, flag,                    &  
                     keyword      = 'POOL_DEPTH',                           &
                     ClientModule = 'DrainageNetwork',                      &
                     SearchType   = FromBlock,                              &
                     Default      = 0.0,                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructNode - ERR61'


        !Variaveis para os calculos hidrodinamicos
  
        call GetData(NewNode%CrossSection%ManningCH,                            &
                 Me%Files%ObjEnterDataNetwork, flag,                            &  
                 keyword      = 'MANNING_CHANNEL',                              &
                 ClientModule = 'DrainageNetwork',                              &
                 SearchType   = FromBlock,                                      &
                 Default      = Me%GlobalManning,                               &
                 STAT         = STAT_CALL)                                      
                                                                                
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructNode - ERR62' 
                
        if (NewNode%CrossSection%ManningCH.LT.0.0) then
            write (*,*)'Invalid Number of Manning Coeficient [MANNING_CHANNEL]'
            stop 'ModuleDrainageNetwork - ConstructNode - ERR62a'
        endif

        if (.not. Me%Continuous) then                            
            call GetData(NewNode%WaterDepth,                                    &
                     Me%Files%ObjEnterDataNetwork, flag,                        &  
                     keyword      = 'WATER_DEPTH',                              &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     default      = Me%InitialWaterDepth,                       &
                     STAT         = STAT_CALL)                                      
                                                                                
            if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ModuleDrainageNetwork - ConstructNode - ERR63'   

            if (NewNode%WaterDepth.LT.0.0) then
                write (*,*)'Invalid Number of Water Level [WATER_DEPTH]'
                stop 'ModuleDrainageNetwork - ConstructNode - ERR63a'
            endif
        end if
                                                                                  
    end subroutine ConstructNode
 
     !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine InitializeTabularCrossSection (CurrNode, ComputeElevation)

        !Arguments-----------------------------------------------------------------
        type (T_Node), pointer                          :: CurrNode
        logical                                         :: ComputeElevation

        !Local---------------------------------------------------------------------
        integer                                         :: N, i, IBL, NLevels, ilev, itab, itabZ
        real                                            :: Z, ZLow, aux
        real                                            :: dH, Av,Pw, Sw, Bw, Sl, Sr        
        real, dimension(:), pointer                     :: Slope, Station, Elevation
        real, dimension(:), allocatable                 :: Elev2
        integer, dimension(1)                           :: k        
        character(len=PathLength)                       :: AuxString
        logical                                         :: mustStop = .false.

        !Calcula area vertical para a seccao toda, com base nos trapezios
        !Depois é so calcular o extra
        !Slope = dx / dy
        !Slope = 0 <=> talude vertical
        !Slope = Inf <=> talude horizontal (nao ha problema porque dH=0)
        
        N = CurrNode%CrossSection%NStations

        Slope     => CurrNode%CrossSection%BankSlope
        Station   => CurrNode%CrossSection%Station
        Elevation => CurrNode%CrossSection%Elevation

        do i = 1,N-1
            aux = Elevation(i+1) - Elevation(i)
            if (abs(aux) > AllmostZero) Slope(i) = (Station(i+1) - Station(i)) / aux
        enddo

        k = minloc(Elevation)
        IBL = k(1)
        CurrNode%CrossSection%IBottom = IBL

        !---------------------------------------------
        !CHECKS
        !---------------------------------------------

        !Check that stations increase

        do i = 1,N-1
            if (Station(i+1)-Station(i) < 0) then
                write(*,*) 'Stations must increase in Node ', CurrNode%ID
                mustStop = .true.
            endif
        enddo    
        if (mustStop) stop 'InitializeTabularCrossSection - ModuleDrainageNetwork - ERR02'

        !Check that Elevations decrease up to bottom level and increase after

        do i=1,IBL-1
            if (Elevation(i+1)-Elevation(i) > 0) then
                write(*,*) 'Left bank elevations must decrease in Node ', CurrNode%ID
                stop 'InitializeTabularCrossSection - ModuleDrainageNetwork - ERR03'
            endif
        enddo

        do i=IBL,N-1
            if (Elevation(i+1)-Elevation(i) < 0) then
                write(*,*) 'Right bank elevations must increase in Node ', CurrNode%ID
                stop 'InitializeTabularCrossSection - ModuleDrainageNetwork - ERR04'
            endif
        enddo

        !Check that Elevations start and end at the same value
        !Change station value to do this

        if (Elevation(1) /= Elevation(N)) then

            Z = min(Elevation(1), Elevation(N))
            if (Elevation(1) /= Z) then

                aux = Station(1) + (Z - Elevation(1)) * Slope(1)

                !write(*,*) 'Changed Station 1 in Node ', CurrNode%ID
                !write(*,*) 'Old : ', Station(1), Elevation(1)
                !write(*,*) 'New : ', aux, Z

                write(AuxString,*) 'Changed station 1 in Node ', CurrNode%ID, ' to ', aux
                call SetError(WARNING_, INTERNAL_, AuxString, ON)

                Station(1) = aux
                Elevation(1) = Z
            
             else

                aux = Station(N-1) + (Z - Elevation(N-1)) * Slope(N-1)

                !write(*,*) 'Changed Station N in Node ', CurrNode%ID
                !write(*,*) 'Old : ', Station(N), Elevation(N)
                !write(*,*) 'New : ', aux, Z

                write(AuxString,*) 'Changed station N in Node ', CurrNode%ID, ' to ', aux
                call SetError(WARNING_, INTERNAL_, AuxString, ON)

                Station(N) = aux
                Elevation(N) = Z 

            endif

        endif

        !ComputeBottomLevel     
        if (ComputeElevation) then

            CurrNode%CrossSection%BottomLevel = CurrNode%CrossSection%TerrainLevel - maxval(Elevation)
        
            do i = 1, N
                CurrNode%CrossSection%Elevation(i) = CurrNode%CrossSection%Elevation(i) + CurrNode%CrossSection%BottomLevel
            enddo
        
        else
            CurrNode%CrossSection%BottomLevel = minval(Elevation)
        endif    

        !Compute Number of Levels
        allocate(Elev2(CurrNode%CrossSection%NStations))
        Elev2   = Elevation       
        NLevels = 1        
        Z       = minval(Elev2)

        do while(Z /= Elevation(1))

            k = minloc(Elev2)   !Lower level
            itab = k(1)
            ZLow = Elev2(itab)

            Elev2(itab) = - null_real
            Z = minval(Elev2)
            k = minloc(Elev2)
            itabZ = k(1)
               
            dH = Z - ZLow

            if (dH > AlmostZero) NLevels = NLevels + 1            

        enddo
    
        CurrNode%CrossSection%NLevels = NLevels

        allocate(CurrNode%CrossSection%Level             (NLevels))
        allocate(CurrNode%CrossSection%LevelSlopeLeft    (NLevels))
        allocate(CurrNode%CrossSection%LevelSlopeRight   (NLevels))
        allocate(CurrNode%CrossSection%LevelBottomWidth  (NLevels))
        allocate(CurrNode%CrossSection%LevelVerticalArea (NLevels))
        allocate(CurrNode%CrossSection%LevelWetPerimeter (NLevels))
        allocate(CurrNode%CrossSection%LevelSurfaceWidth (NLevels))
       
        CurrNode%CrossSection%Level             = null_real
        CurrNode%CrossSection%LevelSlopeLeft    = null_real
        CurrNode%CrossSection%LevelSlopeRight   = null_real
        CurrNode%CrossSection%LevelBottomWidth  = null_real
        CurrNode%CrossSection%LevelVerticalArea = null_real
        CurrNode%CrossSection%LevelWetPerimeter = null_real
        CurrNode%CrossSection%LevelSurfaceWidth = null_real

        !--------------------------------------------------------------
        ! COMPUTE LEVEL PROPERTIES 
        !--------------------------------------------------------------

        !Start at BottomLevel
        CurrNode%CrossSection%Level(1)        = Elevation(IBL)
        CurrNode%CrossSection%LevelVerticalArea(1) = 0.0   
        CurrNode%CrossSection%LevelWetPerimeter(1) = 0.0 
        CurrNode%CrossSection%LevelSurfaceWidth(1) = 0.0 

        !allocate(Elev2(CurrNode%CrossSection%NStations))
        Elev2 = Elevation       
        ilev  = 2        
        Z     = CurrNode%CrossSection%Level(1)

        do while(Z /= Elevation(1))

            k = minloc(Elev2)   !Lower level
            itab = k(1)
            ZLow = Elev2(itab)

            Elev2(itab) = - null_real
            Z = minval(Elev2)
            k = minloc(Elev2)
            itabZ = k(1)
               
            dH = Z - ZLow

            if (dH > AlmostZero) then
                call ComputeExtraArea (CurrNode%CrossSection,       &
                                         itab, dH,      &
                                         SlopeLeft    = Sl,         &
                                         SlopeRight   = Sr,         &
                                         BottomWidth  = Bw,         &
                                         VerticalArea = Av,         &
                                         WetPerimeter = Pw,         &
                                         SurfaceWidth = Sw)


                CurrNode%CrossSection%Level            (ilev)   = Z
                CurrNode%CrossSection%LevelVerticalArea(ilev)   = CurrNode%CrossSection%LevelVerticalArea(ilev-1) + Av
                CurrNode%CrossSection%LevelWetPerimeter(ilev)   = CurrNode%CrossSection%LevelWetPerimeter(ilev-1) + Pw
                CurrNode%CrossSection%LevelSurfaceWidth(ilev)   = Sw
                
                CurrNode%CrossSection%LevelSlopeLeft   (ilev-1) = Sl
                CurrNode%CrossSection%LevelSlopeRight  (ilev-1) = Sr            
                CurrNode%CrossSection%LevelBottomWidth (ilev-1) = Bw

                if (ilev > 2)      &
                    CurrNode%CrossSection%LevelWetPerimeter(ilev) = CurrNode%CrossSection%LevelWetPerimeter(ilev) - Bw

                ilev = ilev + 1
            
            endif

        enddo

        deallocate(Elev2)

        ilev = ilev -1 
        if (ilev /= NLevels) then
            write(*,*) 'NLevels wrong in Node ', CurrNode%ID
            stop 'InitializeTabularCrossSection - ModuleDrainageNetwork - ERR05'
        endif

        !write(*,*) 'Level, VerticalArea, WetPerimeter, SurfaceWidth, SlopeLeft, SlopeRight, BottomWidth'
        !do i = 1,NLevels
        !    write(*,*) CurrNode%CrossSection%Level(i),              &
        !               CurrNode%CrossSection%LevelVerticalArea(i),  &
        !               CurrNode%CrossSection%LevelWetPerimeter(i),  &
        !               CurrNode%CrossSection%LevelSurfaceWidth(i),  &
        !               CurrNode%CrossSection%LevelSlopeLeft   (i),  &
        !               CurrNode%CrossSection%LevelSlopeRight  (i),  &
        !               CurrNode%CrossSection%LevelBottomWidth (i)            
        !enddo

        CurrNode%CrossSection%BottomWidth = CurrNode%CrossSection%LevelBottomWidth(1)

        !for when the vol > vol max (it happens but then restarts)
        CurrNode%CrossSection%LevelSlopeLeft   (NLevels) = CurrNode%CrossSection%LevelSlopeLeft    (NLevels-1)
        CurrNode%CrossSection%LevelSlopeRight  (NLevels) = CurrNode%CrossSection%LevelSlopeRight   (NLevels-1)         
        CurrNode%CrossSection%LevelBottomWidth (NLevels) = CurrNode%CrossSection%LevelSurfaceWidth (NLevels)


    end subroutine InitializeTabularCrossSection

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ComputeExtraArea (CrossSection, ITab, dH, SlopeLeft, SlopeRight, &
                                 BottomWidth, VerticalArea, WetPerimeter, SurfaceWidth)

        !Arguments ------------------------------------------------------------
        type(T_CrossSection)                :: CrossSection
        integer , intent(in)                :: ITab
        real    , intent(in)                :: dH        
        real    , intent(out)               :: SlopeLeft, SlopeRight, BottomWidth
        real    , intent(out)               :: VerticalArea, WetPerimeter, SurfaceWidth

        !Locals----------------------------------------------------------------
        integer                             :: N, I, IBL
        real                                :: OtherSideStation, Z
        real, dimension(:), pointer         :: Slope, Station, Elevation

        N         = CrossSection%NStations
        IBL       = CrossSection%IBottom
        Slope     => CrossSection%BankSlope
        Station   => CrossSection%Station
        Elevation => CrossSection%Elevation

        Z = Elevation(ITab)

        do i = IBL,N
            if (Elevation(i) <= Z) SlopeRight = Slope(i)
        enddo
        
        do i = 1,IBL
            if (Elevation(i) > Z) SlopeLeft = Slope(i)
        enddo
             
        if (ITab <= IBL) then   ! pertence ao left bank
        ! procurar a station no right bank
        ! para calcular a bottom width deste trapezio

            do i = IBL,N
                if (Elevation(i) == Z) then
                    OtherSideStation = Station(i)
                    !exit - aqui nao porque quero a station mais a direita
                elseif (Elevation(i) < Z) then
                    OtherSideStation = Station(i) + Slope(i) * (Z - Elevation(i))
                endif
            enddo

        else    !pertence ao right bank
        ! procurar a station no left bank
        ! para calcular a bottom width deste trapezio

            do i = 1,IBL
                if (Elevation(i) == Z) then
                    OtherSideStation = Station(i)
                    exit !aqui poe-se exit porque quero a station mais a esquerda
                elseif (Elevation(i) > Z) then
                    OtherSideStation = Station(i) + Slope(i) * (Z - Elevation(i))                    
                endif
            enddo

        endif

        BottomWidth  = abs(Station(ITab) - OtherSideStation)

        call TrapezoidGeometry (b  = BottomWidth,   &
                                mL = SlopeLeft,     &
                                mR = SlopeRight,    &
                                h  = dH,            &
                                Av = VerticalArea,  &
                                P  = WetPerimeter,  &
                                Sw = SurfaceWidth) 


        !VerticaArea  = BottomWidth * dH + abs(SlopeLeft) * dH * dH * 0.5 + SlopeRight * dH * dH * 0.5
        !WetPerimeter = ( sqrt( 1. + SlopeLeft**2 ) + sqrt( 1. + SlopeRight**2 ) ) * dH
        !SurfaceWidth = BottomWidth + SlopeLeft * dH + SlopeRight * dH

!100 format(I2, 1x, 2F6.2,1x, F8.2, 1x, F8.2, 1x, 2F6.1,2F8.3)
!        write(*,100) ITab, Z, dH, SlopeLeft, SlopeRight, Station(ITab), OtherSideStation, dH, VerticalArea

        nullify(Slope, Station, Elevation)    

    end subroutine ComputeExtraArea

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine CheckNodesConsistency

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Node), pointer                          :: CurrNode, NextNode
        integer                                         :: NodeID, NextNodeID

        do NodeID = 1, Me%TotalNodes

            CurrNode => Me%Nodes (NodeID)
   
            do NextNodeID = 1, Me%TotalNodes            

                NextNode => Me%Nodes (NextNodeID)
                
                if (NextNodeID /= NodeID) then
                    !Verifies if there are two nodes in the same location
                    if (abs(NextNode%X - CurrNode%X) < AllmostZero) then
                        if (abs(NextNode%Y - CurrNode%Y) < AllmostZero) then
                            write (*,*)'Two nodes at the same location'
                            write (*,*)'Node ID', NodeID
                            stop 'CheckNodesConsistency - ModuleDrainageNetwork - ERR01'
                        endif
                    endif
                endif

            end do
        end do

    end subroutine CheckNodesConsistency

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    
    subroutine ConstructReachList

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine            
        integer                                     :: STAT_CALL, ReachPos


        call CountTotalReaches

        nullify  (Me%Reaches)
        allocate (Me%Reaches (1:Me%TotalReaches)) 
        
        ReachPos = 0       

do1:    do 

            call ExtractBlockFromBuffer(Me%Files%ObjEnterDataNetwork, ClientNumber,      &
                                        BeginReach, EndReach, BlockFound,                &
                                        FirstLine, LastLine, STAT_CALL) 

if1:        if (STAT_CALL .EQ. SUCCESS_) then    

if2:              if (BlockFound) then                 

                    ReachPos = ReachPos + 1   
                                                                                          
                    call ConstructReach (ReachPos)
                                      
                  else if2

                    if (ReachPos /= Me%TotalReaches) stop 'ModuleDrainageNetwork - ConstructReachList - ERR01'

                    call Block_Unlock(Me%Files%ObjEnterDataNetwork, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)                                         &
                        stop 'ModuleDrainageNetwork - ConstructReachList - ERR02'

                    exit do1    !No more blocks

                  end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1

                    stop 'ModuleDrainageNetwork - ConstructReachList - ERR03.'

            end if if1

        end do do1

        !Checks Consistency
        if (Me%CheckReaches) call CheckReachesConsistency
        
        !Calculates Length / Slope
        call CalculateReaches

    end subroutine ConstructReachList 

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine CountTotalReaches
      
        !Local------------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine            
        integer                                     :: STAT_CALL
        integer                                     :: ReachID, OldReachID
        !integer                                     :: MaxReachID, MinReachID
        integer                                     :: flag


        Me%TotalReaches =  0
        !MinReachID = - null_int
        !MaxReachID = null_int
        OldReachID = null_int


do1:    do 

            call ExtractBlockFromBuffer(Me%Files%ObjEnterDataNetwork, ClientNumber,      &
                                        BeginReach, EndReach, BlockFound,                &
                                        FirstLine, LastLine, STAT_CALL) 

if1:        if (STAT_CALL .EQ. SUCCESS_) then    

if2:              if (BlockFound) then                 

                    !Gets ID
                    call GetData(ReachID,                                       &
                                 Me%Files%ObjEnterDataNetwork, flag,            & 
                                 keyword      = 'ID',                           &
                                 ClientModule = 'DrainageNetwork',              &
                                 SearchType   = FromBlock,                      &
                                 STAT         = STAT_CALL)                                  
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - CountTotalReaches - ERR01'

                    if (flag /= 1) then
                        write (*,*)'Invalid Reach ID [ID]'
                        stop 'ModuleDrainageNetwork - CountTotalReaches - ERR02'
                    endif

                    Me%TotalReaches = Me%TotalReaches + 1
                                        
                    !if (ReachID .LT. MinReachID ) MinReachID = ReachID
                    !if (ReachID .GT. MaxReachID ) MaxReachID = ReachID
                   
                    if (ReachID .EQ. OldReachID ) then
                        write (*,*) 'Repeated Reach ID = ', ReachID
                        stop 'ModuleDrainageNetwork - CountTotalReaches - ERR03'
                    else
                        OldReachID = ReachID
                    end if
                    
                                      
                  else if2

                    !if (MinReachID.NE. 1) then
                    !    write (*,*) 'Inconsistency in Reach IDs - Missing ReachID = 1'
                    !    stop 'ModuleDrainageNetwork - CountTotalReaches - ERR04'
                    !else if (MaxReachID.NE. Me%TotalReaches) then
                    !    write (*,*) 'Inconsistency in Reach IDs - Missing ReachID =', Me%TotalReaches
                    !    stop 'ModuleDrainageNetwork - CountTotalReaches - ERR05'
                    !end if

                    call Block_Unlock(Me%Files%ObjEnterDataNetwork, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - CountTotalReaches - ERR06'

                    call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - CountTotalReaches - ERR07'


                    exit do1    !No more blocks

                  end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1

                    stop 'ModuleDrainageNetwork - CountTotalReaches - ERR08.'

            end if if1

        end do do1
        
    end subroutine CountTotalReaches

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ConstructReach (ReachPos)
               
        !Arguments--------------------------------------------------------------
        integer, intent(IN)                         :: ReachPos
        !Local------------------------------------------------------------------
        integer                                     :: DownNodeID, UpNodeID
        type (T_Reach), pointer                     :: NewReach
        integer                                     :: flag, STAT_CALL
        character(LEN = StringLength)               :: str_UpNode, str_DownNode
        logical                                     :: Found

        !------------------------------------------------------------------------

        nullify (NewReach)
        NewReach => Me%Reaches (ReachPos)

        !Gets ID
        call GetData(NewReach%ID,                                               &
                     Me%Files%ObjEnterDataNetwork, flag,                        &
                     keyword      = 'ID',                                       &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructReach - ERR01'
        if (flag /= 1) then
            write (*,*)'Invalid Reach ID [ID]'
            stop 'ModuleDrainageNetwork - ConstructReach - ERR01'
        endif

        !Gets Active Flag
        call GetData(NewReach%Active,                                           &
                     Me%Files%ObjEnterDataNetwork, flag,                        &
                     keyword      = 'ACTIVE',                                   &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     Default        = .TRUE.,                                   &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructReach - ERR01a'
       

        !Gets Downstream Node
        call GetData(DownNodeID,                                                &
                     Me%Files%ObjEnterDataNetwork, flag,                        &  
                     keyword      = 'DOWNSTREAM_NODE',                          &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructReach - ERR02'
        if (flag /= 1) then
            write (*,*)'Invalid Downstream Node [DOWNSTREAM_NODE]'
            stop 'ModuleDrainageNetwork - ConstructReach - ERR02'
        endif

        call FindNodePosition (DownNodeID, NewReach%DownstreamNode, Found)     

        if (.NOT.Found) then
            write (*,*) 'Downstream Node not found'
            write (*,*) 'Node ID = ', DownNodeID
            write (*,*) 'ReachID = ', NewReach%ID
            stop 'ModuleDrainageNetwork - ConstructReach - ERR03'
        end if


        !Gets Upstream Node
        call GetData(UpNodeID,                                                  &
                     Me%Files%ObjEnterDataNetwork, flag,                        &  
                     keyword      = 'UPSTREAM_NODE',                            &
                     ClientModule = 'DrainageNetwork',                          &
                     SearchType   = FromBlock,                                  &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructReach - ERR05'
        if (flag /= 1) then
            write (*,*)'Invalid Upstream Node [UPSTREAM_NODE]'
            stop 'ModuleDrainageNetwork - ConstructReach - ERR04'
        endif

        call FindNodePosition (UpNodeID, NewReach%UpstreamNode, Found)

         if (.NOT.Found) then
            write (*,*) 'Upsream Node not found'
            write (*,*) 'Node ID = ', UpNodeID
            write (*,*) 'ReachID = ', NewReach%ID
            stop 'ModuleDrainageNetwork - ConstructReach - ERR05'
        end if
        
        str_UpNode   =''
        str_DownNode =''
        write(str_UpNode  , '(i10)') UpNodeID
        write(str_DownNode, '(i10)') DownNodeID

        NewReach%Name ='Reach_'//trim(adjustl(adjustr(str_UpNode)))//           &
                       '_'//trim(adjustl(adjustr(str_DownNode)))

    end subroutine ConstructReach
    

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine CheckReachesConsistency      

        !Local------------------------------------------------------------------
        integer                                         :: ReachID, NodeID
        integer                                         :: NextReachID
        type (T_Reach), pointer                         :: CurrReach, NextReach
        logical                                         :: DownstreamNodeExists
        logical                                         :: UpstreamNodeExists

        do ReachID = 1, Me%TotalReaches

        CurrReach => Me%Reaches (ReachID)

            !Verifies if Reach DownStreamNode is not equal UpStreamNode
            if (CurrReach%DownstreamNode ==  CurrReach%UpstreamNode) then
                write (*,*)'Downstream Node must be different from Upstream Node'
                write (*,*)'Reach ID', CurrReach%ID
                write (*,*)'Reach Position', ReachID
                stop 'CheckReachesConsistency - ModuleDrainageNetwork - ERR02'
            endif

            do NextReachID = 1, Me%TotalReaches
    
            NextReach => Me%Reaches (NextReachID)            
                
                if (NextReachID /= ReachID) then
                    !Verifies if there are two reaches with the same end nodes
                    if (CurrReach%UpstreamNode   == NextReach%UpstreamNode .and.               &
                        CurrReach%DownstreamNode == NextReach%DownstreamNode) then
                        write (*,*)'Two reaches with the same end nodes ID'
                        write (*,*)'Reach ID', ReachID
                        stop 'CheckReachesConsistency - ModuleDrainageNetwork - ERR04'
                    endif
                end if 
                               
            end do
            
        end do


        do ReachID = 1, Me%TotalReaches
        
        CurrReach => Me%Reaches (ReachID)
    
            DownstreamNodeExists = .false.
            UpStreamNodeExists   = .false.
            
            do NodeID = 1, Me%TotalNodes
                       
                if (NodeID == CurrReach%DownstreamNode) then
                    DownstreamNodeExists = .true.
                endif
                if (NodeID == CurrReach%UpstreamNode) then
                    UpStreamNodeExists   = .true.
                endif
                
            enddo

            if (.not. DownstreamNodeExists) then
                write (*,*)'Downstream Node does not exists for Reach'
                write (*,*)'Reach Pos', ReachID
                stop 'CheckReachesConsistency - ModuleDrainageNetwork - ERR05'
            endif
            
            if (.not. UpStreamNodeExists) then
                write (*,*)'Upstream Node does not exists for Reach'
                write (*,*)'Reach Pos', ReachID
                stop 'CheckReachesConsistency - ModuleDrainageNetwork - ERR06'
            endif

        enddo    
                
    end subroutine CheckReachesConsistency

    !---------------------------------------------------------------------------

    subroutine CalculateReaches

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: ReachID
        type (T_Reach), pointer                     :: CurrReach
        type (T_Node), pointer                      :: DownstreamNode
        type (T_Node), pointer                      :: UpstreamNode
        character(len=256)                          :: AuxString

        do ReachID = 1, Me%TotalReaches
        
            CurrReach => Me%Reaches (ReachID)
        
            nullify (DownstreamNode)
            nullify (UpstreamNode  )

            DownstreamNode => Me%Nodes (CurrReach%DownstreamNode)
            UpstreamNode   => Me%Nodes (CurrReach%UpstreamNode  )

            if (Me%CoordType == 1) then
                CurrReach%Length = DistanceBetweenTwoGPSPoints(DownstreamNode%X, DownstreamNode%Y, UpstreamNode%X, UpstreamNode%Y)
            else
                CurrReach%Length = sqrt( (DownstreamNode%X - UpstreamNode%X) ** 2. +         &
                                         (DownstreamNode%Y - UpstreamNode%Y) ** 2.)         
            endif
            
            CurrReach%Slope  = (UpstreamNode%CrossSection%BottomLevel -                  &
                                DownstreamNode%CrossSection%BottomLevel) /               &
                                CurrReach%Length

            !Nao sei se isto é boa ideia...                                
            if (CurrReach%Slope <= - Me%MinimumSlope ) then

                write(AuxString,*) 'Negative slope in reach',CurrReach%ID
                call SetError(WARNING_, INTERNAL_, AuxString, ON)                            

!                if (Me%HydrodynamicApproximation == KinematicWave) then
!                    Me%HydrodynamicApproximation = DiffusionWave                
!                    call SetError(WARNING_, INTERNAL_, 'Changing HydrodynamicApproximation to DiffusionWave', ON)            
!                endif

            endif 

            !Sets minimum slope so water doesnt not stop in reaches with little drainage
            !Handle with care...
            !CurrReach%Slope = max(CurrReach%Slope, Me%MinimumSlope)
            if (abs(CurrReach%Slope) < Me%MinimumSlope) then
                CurrReach%Slope = sign(Me%MinimumSlope, CurrReach%Slope)
            endif
            
        enddo    

    end subroutine CalculateReaches

        

    !---------------------------------------------------------------------------

    subroutine ConnectNetwork
    
        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: ReachID, NodeID
        type (T_Node), pointer                      :: CurrNode
        type (T_Reach), pointer                     :: CurrReach
        integer                                     :: iUp, iDown
        
        do NodeID = 1, Me%TotalNodes

            CurrNode => Me%Nodes (NodeID)        

            !Check Upstream /Downstream Reaches to Node
            do ReachID = 1, Me%TotalReaches 
            
                CurrReach => Me%Reaches (ReachID)           

                if (CurrReach%UpstreamNode == NodeID) then
                    CurrNode%nDownstreamReaches = CurrNode%nDownstreamReaches + 1
                endif

                if (CurrReach%DownstreamNode == NodeID) then
                    CurrNode%nUpstreamReaches = CurrNode%nUpstreamReaches     + 1
                endif
            
            enddo

            !Tests consistency
            if (CurrNode%nDownstreamReaches > 1) then
                write (*,*)'Node with more then one Downstream reaches [ID]', NodeID
            endif

            !Allocates DownstreamReaches / UpstreamReaches
            if (CurrNode%nUpstreamReaches > 0) then
                allocate (CurrNode%UpstreamReaches  (CurrNode%nUpstreamReaches))
            endif

            if (CurrNode%nDownstreamReaches > 0) then
                allocate (CurrNode%DownstreamReaches(CurrNode%nDownstreamReaches))
            
            endif
            
            !Fills DownstreamReaches / UpstreamReaches to Node
            iUp   = 0
            iDown = 0

            do ReachID = 1, Me%TotalReaches
                
                CurrReach => Me%Reaches (ReachID)           

                if (CurrReach%UpstreamNode == NodeID) then
                    iDown = iDown + 1
                    CurrNode%DownstreamReaches(iDown) = ReachID
                endif

                if (CurrReach%DownstreamNode == NodeID) then
                    iUp = iUp + 1
                    CurrNode%UpstreamReaches  (iUp  ) = ReachID
                endif

            enddo

       enddo


    end subroutine ConnectNetwork

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    
    subroutine OrderNodes
    
        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: NodePos
        type (T_Node), pointer                      :: CurrNode, UpNode
        type (T_Reach), pointer                     :: UpReach
        integer                                     :: iOrder, i
        logical                                     :: Done

        

        do NodePos = 1, Me%TotalNodes
            Me%Nodes(NodePos)%Order = null_int
        end do

        !Initialize Heads
        do NodePos = 1, Me%TotalNodes            

            if (Me%Nodes(NodePos)%nUpstreamReaches == 0)                          &
                Me%Nodes(NodePos)%Order = 1
        end do

        Done = .false.
do1:    do while (.not.Done)
            
            !Assign other nodes orders        
            do NodePos = 1, Me%TotalNodes            
            
                CurrNode => Me%Nodes (NodePos)            
                
                iOrder = 0
                do i = 1, CurrNode%nUpstreamReaches
                
                    UpReach => Me%Reaches (CurrNode%UpstreamReaches(i))
                    UpNode  => Me%Nodes   (UpReach%UpstreamNode)
                
                    if (UpNode%Order .GT. 0) then
                        iOrder = iOrder + UpNode%Order
                    else
                        iOrder = null_int
                        exit
                    end if

                end do

                CurrNode%Order = 1 + iOrder

            end do

            Done = .true.
            do NodePos = 1, Me%TotalNodes
                if (Me%Nodes(NodePos)%Order .LT. 0.0) then
                    Done = .false.
                    exit
                end if        
            end do

        end do do1

        !With this ordering method, the highest order is equal to the sum of all nodes
        Me%HighestOrder = Me%TotalNodes

    end subroutine OrderNodes

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ReconnectNetwork

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: NodePos, iPos, iOrder, i
        type (T_Node) , dimension(:), pointer       :: NewNodes
        type (T_Node) , pointer                     :: CurrNode
        type (T_Reach), pointer                     :: CurrReach

        nullify (NewNodes)       
        allocate (NewNodes (1:Me%TotalNodes))        

        iPos = 1
        do iOrder  = 1, Me%HighestOrder
        do NodePos = 1, Me%TotalNodes

            CurrNode => Me%Nodes (NodePos)
            
            if (CurrNode%Order == iOrder) then
                
                NewNodes (iPos) = CurrNode
                
                if (CurrNode%nDownstreamReaches == 1) then
                    CurrReach => Me%Reaches (CurrNode%DownstreamReaches (1))
                    CurrReach%UpstreamNode = iPos
                end if

                do i = 1, CurrNode%nUpstreamReaches
                    CurrReach => Me%Reaches (CurrNode%UpstreamReaches (i))
                    CurrReach%DownstreamNode = iPos
                end do
                
                iPos = iPos + 1

            end if

        end do        
        end do

        Me%Nodes = NewNodes
        deallocate (NewNodes)

        call CheckReachesConsistency


    end subroutine ReconnectNetwork
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    
    subroutine WriteOrderedNodes

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: Unit, NodePos
        type (T_Node), pointer                      :: CurrNode       

        call UnitsManager (Unit, OPEN_FILE)
        open (unit = unit, file = trim('nodes ordered.xyz'), status = 'unknown')

        write (unit, *)"<begin_xyz>"

        do NodePos = 1, Me%TotalNodes

            CurrNode => Me%Nodes (NodePos)                        
            write (unit, *) CurrNode%X, CurrNode%Y, CurrNode%Order

        enddo

        write (unit, *)"<end_xyz>"           

        call UnitsManager (Unit, CLOSE_FILE)

    end subroutine WriteOrderedNodes

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine CountOutlets ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                 :: NodePos, ReachPos
        type (T_Node ), pointer                 :: CurrNode

        Me%TotalOutlets = 0

        do NodePos = 1, Me%TotalNodes

            CurrNode => Me%Nodes(NodePos)
            if (CurrNode%nDownstreamReaches .EQ. 0) then
                Me%TotalOutlets = Me%TotalOutlets + 1
                if (CurrNode%nUpstreamReaches .NE. 1) then
                    write (*,*)
                    write (*,*) 'Outlet node must have a single upstream reach'
                    write (*,*) 'Check Node ID = ', CurrNode%ID
                    stop 'CountOutlets - ModuleDrainageNetwork - ERR01'
                end if

                ReachPos = CurrNode%UpstreamReaches (1)                
                Me%OutletReachPos = ReachPos
                Me%OutletNodePos = NodePos
            end if

        end do

        if (Me%TotalOutlets /= 1) stop 'ModuleDrainageNetwork - CountOutlets - ERR01'

    end subroutine CountOutlets

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
        
    subroutine ConstructPropertyList

        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: NewProperty
        integer                                     :: ClientNumber
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound


        ! Initialize the properties number   
        Me%PropertiesNumber = 0

        ! Initialize the properties list   
        nullify (Me%FirstProperty)        

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                        block_begin, block_end, BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                                        
                    call ConstructProperty(NewProperty)
                    
                    call Add_Property(NewProperty)

                else
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructPropertyList - ModuleDrainageNetwork - ERR01'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructPropertyList - ModuleDrainageNetwork - ERR02'
            end if cd1
        end do do1

        if (Me%PropertiesNumber .GT. 0) then
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
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProperty - ModuleDrainageNetwork - ERR01'

        allocate (NewProperty%Concentration            (1:Me%TotalNodes))
        allocate (NewProperty%ConcentrationOld         (1:Me%TotalNodes))
        allocate (NewProperty%InitialConcentration     (1:Me%TotalNodes))
        allocate (NewProperty%InitialConcentrationOld  (1:Me%TotalNodes))
        allocate (NewProperty%MassCreated              (1:Me%TotalNodes))
        allocate (NewProperty%OverLandConc             (1:Me%TotalNodes))
        allocate (NewProperty%GWaterConc               (1:Me%TotalNodes))
        allocate (NewProperty%DWaterConc               (1:Me%TotalNodes))
        allocate (NewProperty%TotalConc                (1:Me%TotalNodes))
        allocate (NewProperty%Load                     (1:Me%TotalNodes))
        allocate (NewProperty%MassInKg                 (1:Me%TotalNodes))
            
       
        NewProperty%Concentration           = 0.0
        NewProperty%ConcentrationOld        = 0.0
        NewProperty%InitialConcentration    = 0.0
        NewProperty%InitialConcentrationOld = 0.0
        NewProperty%MassCreated             = 0.0
        NewProperty%OverLandConc            = 0.0
        NewProperty%GWaterConc              = 0.0
        NewProperty%TotalConc               = 0.0
        NewProperty%Load                    = 0.0
        NewProperty%MassInKg                = 0.0

        call ConstructPropertyID     (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call ConstructPropertyValues (NewProperty)
        
        if (NewProperty%ComputeOptions%Toxicity) then

            nullify(NewProperty%Toxicity%Field)            
            allocate (NewProperty%Toxicity%Field (1:Me%TotalNodes))
            NewProperty%Toxicity%Field = 0.0

        end if

    end subroutine ConstructProperty

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyValues (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),   pointer                 :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        real                                        :: OverLandConcentration
        real                                        :: GWaterConcentration
        real                                        :: DWaterConcentration
        real                                        :: BottomInitialConc
        logical                                     :: Aux
        
        !Begin-----------------------------------------------------------------

        call GetData(NewProperty%ComputeOptions%ComputeLoad,                    &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'COMPUTE_LOAD',                           &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     Default        = .FALSE.,                                  &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR01'
            
        if (NewProperty%ComputeOptions%ComputeLoad) then
            Me%ComputeOptions%ComputeLoad = .true.
        endif

        call GetData(NewProperty%InitialValue,                                  &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'DEFAULT_VALUE',                          &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     STAT           = STAT_CALL)              
        if (iflag /= 0) then
            write(*,*)'The keyword DEFAULT_VALUE in Drainage Network file'
            write(*,*)'is obsolete. Use DEFAULTVALUE instead'
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR01.5'         
        endif
        
        call GetData(NewProperty%InitialValue,                                  &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'DEFAULTVALUE',                           &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     Default        = 0.0,                                      &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR02' 


        call GetData(NewProperty%BoundaryConcentration,                         &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'DEFAULTBOUNDARY',                        &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     Default        = NewProperty%InitialValue,                 &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR02aa' 


        call GetData(OverLandConcentration,                                     &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'OVERLAND_CONCENTRATION',                 &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     Default        = 0.0,                                      &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR02a' 
        
        if (Me%ExtVar%CoupledRP) then
            if (iflag .ne. 0) then
                write(*,*)'Using Module RunoffProperties for overland concentration'
                write(*,*)'keyword OVERLAND_CONCENTRATION in each property is redundant'
                write(*,*)'and not consistent, please remove it.'
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR02b'
            endif
        else
            NewProperty%OverLandConc = OverLandConcentration
        endif
        
        call GetData(GWaterConcentration,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'GROUNDWATER_CONCENTRATION',              &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     Default        = 0.0,                                      &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR02c' 
        
        if (Me%ExtVar%CoupledPMP) then
            if (iflag .ne. 0) then
                write(*,*)'Using Module PorousMediaProperties for groundwater concentration.'
                write(*,*)'keyword GROUNDWATER_CONCENTRATION in each property is redundant'
                write(*,*)'and not consistent, please remove it.'
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR02d'
            endif
        else        
            NewProperty%GWaterConc = GWaterConcentration
        endif

        call GetData(DWaterConcentration,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'DIFFUSEWATER_CONCENTRATION',             &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     Default        = OverLandConcentration,                    &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR02e' 
        NewProperty%DWaterConc = DWaterConcentration


        call GetData(NewProperty%MinValue,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     Keyword        = 'MIN_VALUE',                              &
                     ClientModule   = 'ModuleDrainageNetwork',                  &
                     SearchType     = FromBlock,                                &
                     Default        = 0.0,                                      &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR03'
        
        if (iflag==1)  then
            NewProperty%ComputeOptions%MinConcentration = .true.
        else
            NewProperty%ComputeOptions%MinConcentration = .false.
        endif

        call GetData(NewProperty%ComputeOptions%AdvectionDiffusion,                 &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'ADVECTION_DIFUSION',                         &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = ON,                                           &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR04'

if1:    if (NewProperty%ComputeOptions%AdvectionDiffusion) then
            
            Me%ComputeOptions%AdvectionDiffusion = .true.
                
            !Numerical Discretization of Advection
            call GetData(NewProperty%Advection_Scheme,                             &
                         Me%ObjEnterData, iflag,                                   &
                         Keyword        = 'ADVECTION_SCHEME',                      &
                         ClientModule   = 'ModuleDrainageNetwork',                 &
                         SearchType     = FromBlock,                               &
                         Default        = UpwindOrder1,                            &
                         STAT           = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_)                                           &
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR05'

            if (NewProperty%Advection_Scheme /= UpwindOrder1) then
                write (*,*) 'Invalid option for keyword [ADVECTION_SCHEME]'
                stop 'ConstructPropertyValues - ModuleDrainageNetwork - ERR06'
            end if
!            if (NewProperty%Advection_Scheme /= UpwindOrder1 .AND.                 &
!                NewProperty%Advection_Scheme /= CentralDif) then
!                write (*,*) 'Invalid keyword [ADVECTION_SCHEME]'
!                stop 'ConstructPropertyValues - ModuleDrainageNetwork - ERR06'
!            end if


             call GetData(NewProperty%Diffusion_Scheme,                             &
                          Me%ObjEnterData, iflag,                                   &
                          Keyword        = 'DIFFUSION_SCHEME',                      &
                          ClientModule   = 'ModuleDrainageNetwork',                 &
                          SearchType     = FromBlock,                               &
                          Default        = CentralDif,                              &
                          STAT           = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR07'
 
             
            if (NewProperty%Diffusion_Scheme /= CentralDif) then
                write (*,*) 'Invalid keyword [DIFFUSION_SCHEME]'
                stop 'ConstructPropertyValues - ModuleDrainageNetwork - ERR08'
            end if
 
            !Molecular diffusivity of property in m2/s
            call GetData(NewProperty%Diffusivity,                                   &
                         Me%ObjEnterData, iflag,                                    &
                         Keyword        = 'DIFFUSIVITY',                            &
                         ClientModule   = 'ModuleDrainageNetwork',                  &
                         SearchType     = FromBlock,                                &
                         Default        = 1e-8,                                     &
                         STAT           = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR09' 
        
        
        end if if1

        call GetData(NewProperty%ComputeOptions%Discharges,                         &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'DISCHARGES',                                 &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR10' 

        if (NewProperty%ComputeOptions%Discharges .and. .not. Me%ComputeOptions%Discharges) then
            call SetError(WARNING_, INTERNAL_, 'Missing keyword [DISCHARGES]', ON)            
            Me%ComputeOptions%Discharges = ON
        end if
        
        if (NewProperty%ComputeOptions%Discharges) Me%nPropWithDischarges = Me%nPropWithDischarges + 1

        call GetData(NewProperty%ComputeOptions%Toxicity,                           &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'TOXICITY',                                   &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR11' 

        if (NewProperty%ComputeOptions%Toxicity .and. .not. NewProperty%ComputeOptions%Discharges) then
            write (*,*) 'Toxic property NOT discharged',                        &
                         trim(adjustl(adjustr(NewProperty%ID%Name)))
            !stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR12' 
        end if
         
        if (NewProperty%ComputeOptions%Toxicity .and. .not. Me%ComputeOptions%Toxicity) then            

            Me%ComputeOptions%Toxicity = ON
            nullify     (Me%GlobalToxicity)                       
            allocate    (Me%GlobalToxicity (1:Me%TotalNodes))
            Me%GlobalToxicity = 0.0

        end if

        if (NewProperty%ComputeOptions%Toxicity) Me%nToxicProp = Me%nToxicProp + 1


ifTox:  if (NewProperty%ComputeOptions%Toxicity) then

            call GetData(NewProperty%Toxicity%Evolution,                            &   
                         Me%ObjEnterData, iflag,                                    &
                         Keyword        = 'TOX_EVOLUTION',                          &
                         ClientModule   = 'ModuleDrainageNetwork',                  &
                         SearchType     = FromBlock,                                &
                         Default        = RiskRatio,                                &
                         STAT           = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR13' 

            if (NewProperty%Toxicity%Evolution /= Saturation .and.                  &
                NewProperty%Toxicity%Evolution /= Linear     .and.                  &
                NewProperty%Toxicity%Evolution /= RiskRatio ) then
                write (*,*) 'Toxicity Evolution badly defined',                     &
                trim(adjustl(adjustr(NewProperty%ID%Name)))
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR13a' 
            end if


if2:        if (NewProperty%Toxicity%Evolution == Saturation .OR.                   &
                NewProperty%Toxicity%Evolution == RiskRatio        )  then

                !WARNING - this is in fraction of initial concentration units [%]
                !See header of this module for more information on toxic model
                call GetData(NewProperty%Toxicity%EC50,                             &   
                             Me%ObjEnterData, iflag,                                &
                             Keyword        = 'EC50',                               &
                             ClientModule   = 'ModuleDrainageNetwork',              &
                             SearchType     = FromBlock,                            &
                             Default        = 0.5,                                  &
                             STAT           = STAT_CALL)              
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR13b'
            else !if2
            
                call GetData(NewProperty%Toxicity%Slope,                            &   
                             Me%ObjEnterData, iflag,                                &
                             Keyword        = 'SLOPE',                              &
                             ClientModule   = 'ModuleDrainageNetwork',              &
                             SearchType     = FromBlock,                            &
                             Default        = 1.0,                                  &
                             STAT           = STAT_CALL)              
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR13c'
           
           end if if2


        end if ifTox


        call GetData(Aux,                                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DECAY',                                             &
                     ClientModule = 'ModuleDrainageNetwork',                             &
                     STAT         = STAT_CALL)
        if (iflag == 1) then
            write(*,*)
            write(*,*) 'ERROR: '
            write(*,*) 'DECAY keyword is now obsolete in drainage properties'
            write(*,*) 'if want to simulate Coliform decay please use DECAY_T90'
            write(*,*) 'if want to simulate generic decay please use DECAY_GENERIC'
            stop 'ConstructPropertyValues - ModuleDrainageNetwork - ERR60'
        endif
        
        !Dacay Time
        call GetData(NewProperty%ComputeOptions%T90_Decay,                          &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'DECAY_T90',                                  &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR11' 
        
        if (NewProperty%ComputeOptions%T90_Decay .and. .not. NewProperty%ComputeOptions%Discharges) then
            write (*,*) 'Decaying properties must be discharged',                   &
                         trim(adjustl(adjustr(NewProperty%ID%Name)))
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR12' 
        end if
         
        if (NewProperty%ComputeOptions%T90_Decay ) then
            Me%ComputeOptions%T90_Decay = ON
        endif

        !Generic decay
        call GetData(NewProperty%ComputeOptions%Generic_Decay,                           &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DECAY_GENERIC',                                     &
                     ClientModule = 'ModuleDrainageNetwork',                             &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ConstructPropertyValues - ModuleDrainageNetwork - ERR50'
        if (NewProperty%ComputeOptions%Generic_Decay) then
            Me%ComputeOptions%Generic_Decay       = .true.
        endif

        if (NewProperty%ComputeOptions%Generic_Decay) then
            
            !Decay rate k (s-1) in P = Po*exp(-kt)
            call GetData(NewProperty%DecayRate,                                              &
                         Me%ObjEnterData,iflag,                                              &
                         SearchType   = FromBlock,                                           &
                         keyword      = 'DECAY_RATE',                                        &
                         ClientModule = 'ModuleDrainageNetwork',                             &
                         default      = 0.,                                                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                     &
                stop 'ConstructPropertyValues - ModuleDrainageNetwork - ERR70'            
            
        endif

        !Checks for Surface Fluxes
        call GetData(NewProperty%ComputeOptions%SurfaceFluxes,                      &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'SURFACE_FLUXES',                               &
                     Default      = .false.,                                        & 
                     ClientModule = 'ModuleDrainageNetwork',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

        if (NewProperty%ComputeOptions%SurfaceFluxes)                               &
            Me%ComputeOptions%SurfaceFluxes = .true.

        !Checks for Bottom Fluxes
        call GetData(NewProperty%ComputeOptions%BottomFluxes,                       &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'BOTTOM_FLUXES',                                &
                     Default      = .false.,                                        & 
                     ClientModule = 'ModuleDrainageNetwork',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

        if (NewProperty%ComputeOptions%BottomFluxes) then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 
            end if            
        endif

ifB:    if (NewProperty%ComputeOptions%BottomFluxes) then
            Me%ComputeOptions%BottomFluxes = .true.

            allocate (NewProperty%BottomConc     (1:Me%TotalNodes))
            allocate (NewProperty%ErosionRate    (1:Me%TotalNodes))
            allocate (NewProperty%DepositionRate (1:Me%TotalNodes))
            allocate (NewProperty%Ws             (1:Me%TotalNodes))

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
                         ClientModule = 'ModuleDrainageNetwork',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 
            
            NewProperty%BottomConc = BottomInitialConc

            !Bottom Initial Concentration
            call GetData(NewProperty%BottomMinConc,                                     &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOTTOM_MIN_CONC',                              &
                         Default      =  0.0,                                           & 
                         ClientModule = 'ModuleDrainageNetwork',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

            !Compute erosion fluxes
            call GetData(NewProperty%ComputeOptions%Erosion,                            &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'EROSION',                                      &
                         Default      =  .true.,                                        & 
                         ClientModule = 'ModuleDrainageNetwork',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

            if (NewProperty%ComputeOptions%Erosion) then
                !Critial Erosion Shear Stress [Pa]
                call GetData(NewProperty%ErosionCriticalShear,                              &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'CRIT_SS_EROSION',                              &
                             Default      =  0.2,                                           & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

                !Erosion Coefficient [kg m-2 s-1]
                call GetData(NewProperty%ErosionCoefficient,                                &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'EROSION_COEF',                                 &
                             Default      =  5.0E-4,                                        & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 
            
            end if


            !Compute deposition fluxes
            call GetData(NewProperty%ComputeOptions%Deposition,                         &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'DEPOSITION',                                   &
                         Default      =  .true.,                                        & 
                         ClientModule = 'ModuleDrainageNetwork',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

            if (NewProperty%ComputeOptions%Deposition) then
                !Critial Deposition Shear Stress [Pa]
                call GetData(NewProperty%DepositionCriticalShear,                           &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'CRIT_SS_DEPOSITION',                           &
                             Default      =  0.1,                                           & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

                if (NewProperty%ComputeOptions%Erosion .and. &
                    NewProperty%DepositionCriticalShear >= NewProperty%ErosionCriticalShear) then
                    write (*,*) '[CRIT_SS_EROSION] must be higher than [CRIT_SS_DEPOSITION]'
                    stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 
                end if
    

                !See ModuleFreeVerticalMovement - Hindered settling  - CHS - kg m-3
                call GetData(NewProperty%CHS,                                               &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'CHS',                                          &
                             Default      =  4.0,                                           & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

                !See ModuleFreeVerticalMovement - Hindered settling  - CHS
                !Settling type: WSConstant = 1, SPMFunction = 2
                call GetData(NewProperty%Ws_Type,                                           &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'WS_TYPE',                                      &
                             Default      =  WSConstant,                                    & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

                !See ModuleFreeVerticalMovement - Constant settling velocity [m s-1]
                call GetData(NewProperty%Ws_Value,                                          &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'WS_VALUE',                                     &
                             Default      =  0.0001,                                        & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 


                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%KL,                                                &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'KL',                                           &
                             Default      =  0.1,                                           & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%KL1,                                               &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'KL1',                                          &
                             Default      =  0.1,                                           & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%ML,                                                &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'ML',                                           &
                             Default      =  4.62,                                          & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 


                !See ModuleFreeVerticalMovement
                call GetData(NewProperty%M,                                                 &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'M',                                            &
                             Default      =  1.0,                                           & 
                             ClientModule = 'ModuleDrainageNetwork',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR14' 

            end if

        end if ifB
        
        !Checks for Water Quality model
        call GetData(NewProperty%ComputeOptions%WaterQuality,                       &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'WATER_QUALITY',                              &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR15' 
        
        if (NewProperty%ComputeOptions%WaterQuality) then
            Me%ComputeOptions%WaterQuality = .true.
        end if


        !Checks for Benthos model
        call GetData(NewProperty%ComputeOptions%Benthos,                            &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'BENTHOS',                                    &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR15a' 
        
        if (NewProperty%ComputeOptions%Benthos) then
            
            Me%ComputeOptions%Benthos = .true.

        end if


        !Checks for CEQUAL_W2 model
        call GetData(NewProperty%ComputeOptions%CeQualW2,                           &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'CEQUALW2',                                   &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR16' 
        
        if (NewProperty%ComputeOptions%CeQualW2) then
            Me%ComputeOptions%CeQualW2 = .true.
        end if

        !Checks for Life model
        call GetData(NewProperty%ComputeOptions%Life,                               &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'LIFE',                                       &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR16' 
        
        if (NewProperty%ComputeOptions%Life) then
            Me%ComputeOptions%Life = .true.
        end if

        !Checks if user wants to calculate total Concentration (Column + Bottom)
        call GetData(NewProperty%ComputeOptions%SumTotalConc,                       &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'SUMTOTALCONC',                                 &
                     Default      =  .FALSE.,                                       & 
                     ClientModule = 'ModuleDrainageNetwork',                        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR16a'

        if (NewProperty%ComputeOptions%SumTotalConc) then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' is not'
                write(*,*) 'recognised as PARTICULATE and does not have Bottom_ or total_Conc'
                stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR16b' 
            end if            
            
            Me%ComputeOptions%SumTotalConc = .true.
        endif


        !IS Coeficient
        call GetData(NewProperty%IScoefficient,                                     &
                     Me%ObjEnterData, iflag,                                        &
                     KeyWord        = 'IS_COEF',                                    &
                     Default        = 1.e-3,                                        &      
                     SearchType     = FromBlock,                                    &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     STAT           = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR17' 

        !IS ExtinctionCoef
        call GetData(NewProperty%ExtinctionCoefficient,                             &
                     Me%ObjEnterData, iflag,                                        &
                     KeyWord        = 'EXTINCTION_PARAMETER',                       &
                     Default        = 1.0,                                          &      
                     SearchType     = FromBlock,                                    &
                     ClientModule   = 'ModuleDrainageNetwork',                      &
                     STAT           = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR18' 

        call GetData(NewProperty%ComputeOptions%TimeSerie,                          &
             Me%ObjEnterData, iflag,                                                &
             Keyword        = 'TIME_SERIE',                                         &
             ClientModule   = 'ModuleWaterProperties',                              &
             SearchType     = FromBlock,                                            &
             Default        = .false.,                                              &
             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR99'

        call GetData(NewProperty%OutputName,                                        &
             Me%ObjEnterData, iflag,                                                &
             Keyword        = 'OUTPUT_NAME',                                        &
             ClientModule   = 'ModuleWaterProperties',                              &
             SearchType     = FromBlock,                                            &
             Default        = 'NAME',                                               &
             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR99a'

        if (NewProperty%OutputName /= 'NAME' .and. NewProperty%OutputName /= 'DESCRIPTION') &
            stop 'ModuleDrainageNetwork - ConstructPropertyValues - ERR99b'

    end subroutine ConstructPropertyValues    
    
    !---------------------------------------------------------------------------
    
    subroutine InitializeProperty
    
        !Arguments--------------------------------------------------------------
        
        !Local------------------------------------------------------------------
        type(T_Property), pointer                   :: Property

        Property => Me%FirstProperty
        do while (associated(Property))
               
            if (.not. Me%Continuous) then
                Property%Concentration = Property%InitialValue
                
                Property%Concentration(Me%OutletNodePos) = Property%BoundaryConcentration
                
            endif
                
                
            call ComputeToxicityForEachEffluent
                       
            Property => Property%Next
        end do
       
    end subroutine InitializeProperty
    
    !--------------------------------------------------------------------------- 
    
    subroutine Add_Property(NewProperty)

        !Arguments--------------------------------------------------------------
        type(T_Property),           pointer     :: NewProperty


        if (.not.associated(Me%FirstProperty)) then
            Me%PropertiesNumber  = 1
            Me%FirstProperty     => NewProperty
            Me%LastProperty      => NewProperty
        else
            NewProperty%Prev     => Me%LastProperty
            Me%LastProperty%Next => NewProperty
            Me%LastProperty      => NewProperty
            Me%PropertiesNumber  = Me%PropertiesNumber + 1
        end if 

    end subroutine Add_Property 

    !---------------------------------------------------------------------------   

  
    subroutine CheckSelectedProp  ()  

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------

        !When Fecal_Coliforms are simulated, Temperature and sal are neeed,  too
        if (PropertyExists(Fecal_Coliforms_)) then            

            if (.not. PropertyExists(Temperature_) .or. .not. PropertyExists(Salinity_))  then 
                write (*,*) 'For coliform decay Temperature and salinity are needed'
                stop 'CheckSelectedProp - ModuleDrainageNetwork - ERR01' 
            endif
        
        endif
        

    end subroutine CheckSelectedProp    

    !---------------------------------------------------------------------------   

    subroutine InitializeVariables

        !Local------------------------------------------------------------------
        type (T_Node), pointer                      :: CurrNode        
        integer                                     :: NodeID, iUp, iDown, NLevels
        type (T_Reach), pointer                     :: UpReach, DownReach
        real                                        :: Av, P, Sw, m, TopH

        allocate (Me%RunOffVector       (Me%TotalNodes))
        allocate (Me%GroundVector       (Me%TotalNodes))
        allocate (Me%GWFlowBottomLayer  (Me%TotalNodes))
        allocate (Me%GWFlowTopLayer     (Me%TotalNodes))
        allocate (Me%DiffuseVector      (Me%TotalNodes))
        allocate (Me%ComputeFaces       (Me%TotalNodes))
        allocate (Me%OpenPointsFlow     (Me%TotalNodes))
        allocate (Me%OpenPointsProcess  (Me%TotalNodes))
        allocate (Me%RiverPoints        (Me%TotalNodes))
        
        if (Me%ComputeOptions%TransmissionLosses) then
            allocate (Me%TransmissionFlow (Me%TotalNodes))
            Me%TransmissionFlow = 0.0
        endif

        if  (Me%ComputeOptions%SurfaceFluxes .or. Me%ComputeOptions%WaterQuality .or. &
             Me%ComputeOptions%CeQualW2 .or. Me%ComputeOptions%Life)  then
            Me%ComputeOptions%TopRadiation = .true.
        endif

        Me%RunOffVector     = 0.0  
        Me%GroundVector     = 0.0
        Me%GWFlowBottomLayer = null_int
        Me%GWFlowTopLayer    = null_int
        Me%DiffuseVector    = 0.0
        Me%ComputeFaces     = 0  
        Me%OpenPointsFlow   = 0  
        Me%OpenPointsProcess= 0  
        Me%RiverPoints      = 1

        nullify (CurrNode)

        do NodeID = 1, Me%TotalNodes
        
            CurrNode => Me%Nodes (NodeID)               

            !Length of reaches in a node control volume
            CurrNode%Length = 0.0
        
            do iUp=1, CurrNode%nUpstreamReaches
                nullify (UpReach)
                UpReach => Me%Reaches (CurrNode%UpstreamReaches (iUp))
                CurrNode%Length = CurrNode%Length + UpReach%Length / 2.0
            end do

            do iDown=1, CurrNode%nDownstreamReaches
                nullify (DownReach) 
                DownReach => Me%Reaches (CurrNode%DownstreamReaches  (iDown))
                CurrNode%Length = CurrNode%Length + DownReach%Length / 2.0
            end do
            
            !if (CurrNode%nUpstreamReaches == 0 .or. CurrNode%nDownstreamReaches == 0) then
            !    CurrNode%Length = 2.* CurrNode%Length
            !endif

        end do            
                            
        !if continuous, read water depth from Initial File
        !if not continuous, water depth initialized in subroutine ReadDataFile            
        
        if (Me%Continuous) call ReadInitialFile          

        call InitializeNodes
        call InitializeReaches
        call InitializeProperty

if1:    if (Me%HasGrid) then

            !Channels WaterLevel
            allocate(Me%ChannelsWaterLevel(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%ChannelsWaterLevel   = 0.0
            call UpdateChannelsDynamicMatrix

            allocate(Me%ChannelsBottomLevel   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%ChannelsBottomWidth   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%ChannelsSurfaceWidth  (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))            
            allocate(Me%ChannelsBankSlope     (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%ChannelsNodeLength    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%ChannelsOpenProcess   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%ChannelsID            (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%ChannelsStationName   (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%ChannelsBottomLevel  = null_real
            Me%ChannelsBottomWidth  = null_real
            Me%ChannelsSurfaceWidth = null_real
            Me%ChannelsBankSlope    = null_real
            Me%ChannelsNodeLength   = null_real
            Me%ChannelsOpenProcess  = null_int
            Me%ChannelsID           = null_int
            Me%ChannelsStationName  = ''

            do NodeID = 1, Me%TotalNodes

                CurrNode => Me%Nodes (NodeID)         

                Me%ChannelsID          (CurrNode%GridI, CurrNode%GridJ) = CurrNode%ID
                Me%ChannelsStationName (CurrNode%GridI, CurrNode%GridJ) = CurrNode%StationName
                       
                Me%ChannelsBottomLevel  (CurrNode%GridI, CurrNode%GridJ) = CurrNode%CrossSection%BottomLevel
                Me%ChannelsNodeLength   (CurrNode%GridI, CurrNode%GridJ) = CurrNode%Length

                if (CurrNode%CrossSection%Form == Trapezoidal) then

                    call TrapezoidGeometry (b  = CurrNode%CrossSection%BottomWidth, &
                                            mR = CurrNode%CrossSection%Slope,       &
                                            mL = CurrNode%CrossSection%Slope,       &
                                            h  = CurrNode%CrossSection%Height,      &
                                            Av = Av,                                &
                                            P  = P,                                 &
                                            Sw = Sw)

                    Me%ChannelsBottomWidth  (CurrNode%GridI, CurrNode%GridJ) = CurrNode%CrossSection%BottomWidth
                    Me%ChannelsSurfaceWidth (CurrNode%GridI, CurrNode%GridJ) = Sw
                    Me%ChannelsBankSlope    (CurrNode%GridI, CurrNode%GridJ) = CurrNode%CrossSection%Slope

                elseif (CurrNode%CrossSection%Form == TrapezoidalFlood) then

                    TopH = CurrNode%CrossSection%Height - CurrNode%CrossSection%MiddleHeight

                    call TrapezoidGeometry (b  = CurrNode%CrossSection%MiddleWidth, &
                                            mR = CurrNode%CrossSection%SlopeTop,    &
                                            mL = CurrNode%CrossSection%SlopeTop,    &
                                            h  = TopH,                              &
                                            Av = Av,                                &
                                            P  = P,                                 &
                                            Sw = Sw)

                    Me%ChannelsBottomWidth  (CurrNode%GridI, CurrNode%GridJ) = CurrNode%CrossSection%BottomWidth
                    Me%ChannelsSurfaceWidth (CurrNode%GridI, CurrNode%GridJ) = Sw
                    Me%ChannelsBankSlope    (CurrNode%GridI, CurrNode%GridJ) = CurrNode%CrossSection%SlopeTop

                elseif (CurrNode%CrossSection%Form == Tabular) then
                    
                    NLevels = CurrNode%CrossSection%NLevels    

                    m = 0.5 * ( abs(CurrNode%CrossSection%LevelSlopeLeft(NLevels))  &
                            +  CurrNode%CrossSection%LevelSlopeRight(NLevels) )

                    Me%ChannelsBottomWidth  (CurrNode%GridI, CurrNode%GridJ) = CurrNode%CrossSection%LevelWetPerimeter(NLevels)
                    Me%ChannelsSurfaceWidth (CurrNode%GridI, CurrNode%GridJ) = CurrNode%CrossSection%LevelSurfaceWidth(NLevels)
                    Me%ChannelsBankSlope    (CurrNode%GridI, CurrNode%GridJ) = m
        
                endif

            enddo

        endif if1
                            
    end subroutine InitializeVariables

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ReadInitialFile


        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------

        !Local------------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: InitialFile
        type (T_Time)                               :: EndTimeFile
        real                                        :: DT_error
        integer                                     :: STAT_CALL
        type(T_Property), pointer                   :: Property        

        !-----------------------------------------------------------------------

        call UnitsManager(InitialFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadInitialFile - ERR01'

        open(Unit = InitialFile, File = Me%Files%Initial, Form = 'UNFORMATTED',     &
             status = 'OLD', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadInitialFile - ERR02'

        !Reads Date
        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

       
        DT_error = EndTimeFile - Me%BeginTime

        !Avoid rounding errors - Frank 08-2001
        if (abs(DT_error) >= 0.01) then
            
            write(*,*) 'The end time of the previous run is different from the start time of this run'
            write(*,*) 'Date in the file'
            write(*,*) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
            write(*,*) 'DT_error', DT_error
            if (Me%StopOnWrongDate) stop 'ModuleDrainageNetwork - ReadInitialFile - ERR04'   

        endif

        read(InitialFile)Me%Nodes%WaterLevel
        read(InitialFile)Me%Reaches%FlowNew

        Property => Me%FirstProperty
        do while (associated(Property))
            read (InitialFile) Property%Concentration
            Property => Property%Next
        end do

        read(InitialFile)Me%LastGoodNiter

        Property => Me%FirstProperty
        do while (associated(Property))
            read (InitialFile, err=10) Property%BottomConc
            Property => Property%Next
        end do

    10  continue

        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - ReadInitialFile - ERR04'

    end subroutine ReadInitialFile

    !---------------------------------------------------------------------------

    subroutine InitializeNodes
        
        !Local------------------------------------------------------------------
        type (T_Node), pointer                          :: CurrNode
        integer                                         :: NodeID
        real                                            :: TopH, AvTrapez1, AvTrapez2
        real                                            :: Av, Pw, Sw
        real                                            :: DownStreamLevel
                    
        if (Me%Continuous) then
            
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes(NodeID)
                CurrNode%WaterDepth = CurrNode%WaterLevel - CurrNode%CrossSection%BottomLevel
            end do

        else

            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes(NodeID)
                CurrNode%WaterLevel = CurrNode%WaterDepth + CurrNode%CrossSection%BottomLevel
            end do

        end if

        
        do NodeID = 1, Me%TotalNodes    
            
            CurrNode => Me%Nodes(NodeID)

            call ComputeXSFromWaterDepth (CurrNode)
                                           
if1:        if (CurrNode%nDownstreamReaches /= 0) then
            
                if (CurrNode%CrossSection%Form .EQ. Trapezoidal) then            
                               
                    CurrNode%VolumeMax = (( CurrNode%CrossSection%BottomWidth   &
                                         + CurrNode%CrossSection%Slope          &
                                         * CurrNode%CrossSection%Height )       &
                                         * CurrNode%CrossSection%Height )       &
                                         * CurrNode%Length

                    CurrNode%VolumeMin = (( CurrNode%CrossSection%BottomWidth   &
                                         + CurrNode%CrossSection%Slope          &
                                         * Me%MinimumWaterDepth)                &
                                         * Me%MinimumWaterDepth)                &
                                         * CurrNode%Length
                                         
                else if (CurrNode%CrossSection%Form .EQ. TrapezoidalFlood) then

                    AvTrapez1 = ( CurrNode%CrossSection%BottomWidth             &
                                 + CurrNode%CrossSection%Slope                  &
                                 * CurrNode%CrossSection%MiddleHeight )         &
                                 * CurrNode%CrossSection%MiddleHeight

                    TopH      = CurrNode%CrossSection%Height - CurrNode%CrossSection%MiddleHeight

                    AvTrapez2 = ( CurrNode%CrossSection%MiddleWidth             &
                                 + CurrNode%CrossSection%SlopeTop               &
                                 * TopH ) * TopH

                    CurrNode%VolumeMaxTrapez1 = AvTrapez1 * CurrNode%Length
                    CurrNode%VolumeMax        = (AvTrapez1 + AvTrapez2) * CurrNode%Length

                    if (Me%MinimumWaterDepth > CurrNode%CrossSection%MiddleHeight) then

                        TopH = Me%MinimumWaterDepth - CurrNode%CrossSection%MiddleHeight

                        AvTrapez2 = ( CurrNode%CrossSection%MiddleWidth         &
                                     + CurrNode%CrossSection%SlopeTop           &
                                     * TopH ) * TopH
                    else 

                        AvTrapez1 = ( CurrNode%CrossSection%BottomWidth         &
                                     + CurrNode%CrossSection%Slope              &
                                     * Me%MinimumWaterDepth )                   &
                                     * Me%MinimumWaterDepth

                        AvTrapez2 = 0.0

                    endif

                    CurrNode%VolumeMin = (AvTrapez1 + AvTrapez2) * CurrNode%Length

                elseif (CurrNode%CrossSection%Form .EQ. Tabular) then

                    Av = CurrNode%CrossSection%LevelVerticalArea (CurrNode%CrossSection%NLevels)
                    CurrNode%VolumeMax = Av * CurrNode%Length

                    TopH = CurrNode%CrossSection%BottomLevel + Me%MinimumWaterDepth
                    call TabularGeometry (CurrNode%CrossSection, TopH, Av, Pw, Sw)                    
                    CurrNode%VolumeMin = Av * CurrNode%Length

                else                 
                    stop 'Invalid cross section form - InitializeNodes - ModuleDrainageNetwork - ERR01'
                end if
                                
                CurrNode%VolumeOld = CurrNode%VolumeNew     
                    
                CurrNode%VolumeMax = CurrNode%SingCoef * CurrNode%VolumeMax
            
            else 
            
                if (Me%Downstream%Boundary == ImposedWaterLevel) then !if1
                
                    if (Me%Downstream%Evolution  == None .or. Me%DownStream%Evolution == OpenMI) then
                       CurrNode%WaterLevel = Me%Downstream%DefaultValue
                    else if (Me%Downstream%Evolution == ReadTimeSerie) then
                       call ModifyDownstreamTimeSerie (CurrNode%WaterLevel)
                    end if

                    CurrNode%WaterDepth = CurrNode%WaterLevel - CurrNode%CrossSection%BottomLevel
                          
                    DownStreamLevel = CurrNode%WaterLevel
                    
                else

                    DownStreamLevel = CurrNode%WaterLevel
                    
                endif
            
            end if if1

        end do

        !Downstream level, so inicial water level makes sense
         if (.not. Me%Continuous) then
         
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes(NodeID)
                
                if (CurrNode%WaterLevel < DownStreamLevel) then
                
                    CurrNode%WaterLevel = DownStreamLevel
                    
                    CurrNode%WaterDepth = CurrNode%WaterLevel - CurrNode%CrossSection%BottomLevel
                
                    call ComputeXSFromWaterDepth (CurrNode)
                    
                endif
                
            enddo         
         endif


        !if (Me%CheckMass) then
        !    Me%TotalStoredVolume = 0.0
        !    do NodeID = 1, Me%TotalNodes
        !        if (Me%Nodes(NodeID)%nDownStreamReaches /= 0) then
        !            Me%TotalStoredVolume = Me%TotalStoredVolume + Me%Nodes(NodeID)%VolumeNew
        !        endif
        !    end do
        !end if


    end subroutine InitializeNodes
    
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ComputeXSFromWaterDepth (CurrNode)

        !Arguments--------------------------------------------------------------
        type (T_Node), pointer                          :: CurrNode
        
        !Local------------------------------------------------------------------        
        real(8)                                         :: PoolVolume
        real                                            :: TopH, AvTrapez1, AvTrapez2
        real                                            :: PTrapez1, aux


        if (CurrNode%nDownstreamReaches /= 0) then

            if (CurrNode%CrossSection%Form == Trapezoidal .OR.            &
               (CurrNode%CrossSection%Form == TrapezoidalFlood .AND.      &
                CurrNode%WaterDepth <= CurrNode%CrossSection%MiddleHeight))  then

                call TrapezoidGeometry (b  = CurrNode%CrossSection%BottomWidth, &
                                        mR = CurrNode%CrossSection%Slope,       &
                                        mL = CurrNode%CrossSection%Slope,       &
                                        h  = CurrNode%WaterDepth,               &
                                        Av = CurrNode%VerticalArea,             &
                                        P  = CurrNode%WetPerimeter,             &
                                        Sw = CurrNode%SurfaceWidth)

            elseif (CurrNode%CrossSection%Form == TrapezoidalFlood) then
            ! from the previous if
            ! we already know that CurrNode%WaterDepth > CurrNode%CrossSection%MiddleHeigh

                call TrapezoidGeometry (b  = CurrNode%CrossSection%BottomWidth, &
                                        mR = CurrNode%CrossSection%Slope,       &
                                        mL = CurrNode%CrossSection%Slope,       &
                                        h  = CurrNode%CrossSection%MiddleHeight,&
                                        Av = AvTrapez1,                         &
                                        P  = PTrapez1,                          &
                                        Sw = aux)

                TopH = CurrNode%WaterDepth - CurrNode%CrossSection%MiddleHeight

                call TrapezoidGeometry (b  = CurrNode%CrossSection%MiddleWidth, &
                                        mR = CurrNode%CrossSection%SlopeTop,    &
                                        mL = CurrNode%CrossSection%SlopeTop,    &
                                        h  = TopH,                              &
                                        Av = AvTrapez2,                         &
                                        P  = aux,                               &
                                        Sw = CurrNode%SurfaceWidth)

                CurrNode%VerticalArea = AvTrapez1 + AvTrapez2
                CurrNode%WetPerimeter = PTrapez1 + 2. * TopH * sqrt (1. + CurrNode%CrossSection%SlopeTop**2.)


            elseif (CurrNode%CrossSection%Form == Tabular) then
                
                call TabularGeometry (CurrNode%CrossSection,   &
                                      CurrNode%WaterLevel,     &
                                      CurrNode%VerticalArea,   & 
                                      CurrNode%WetPerimeter,   &
                                      CurrNode%SurfaceWidth)

            else                 
                stop 'Invalid cross section form - ComputeXSFromWaterDepth - ModuleDrainageNetwork - ERR01'
            end if


            CurrNode%SurfaceArea  = CurrNode%SurfaceWidth * CurrNode%Length

            !Correction of Surface Area for low water depth - invented by Frank
            if (CurrNode%CrossSection%PoolDepth < AllmostZero) then
                if     (CurrNode%WaterDepth < 0.01) then
                    CurrNode%SurfaceArea = CurrNode%SurfaceArea  * 0.1
                    CurrNode%SurfaceWidth= CurrNode%SurfaceWidth * 0.1
                elseif (CurrNode%WaterDepth < 0.05) then
                    CurrNode%SurfaceArea = CurrNode%SurfaceArea  * 0.5
                    CurrNode%SurfaceWidth= CurrNode%SurfaceWidth * 0.5
                endif
            endif    

            if (CurrNode%VerticalArea .LT. AllmostZero) then
                !Pools are always half in the beginning
                PoolVolume          = 0.5 * CurrNode%CrossSection%PoolDepth * CurrNode%CrossSection%BottomWidth * CurrNode%Length
                CurrNode%VolumeNew  = 0.0 + PoolVolume
            else
                !Pools are always full in the beginning
                PoolVolume            = 1.0 * CurrNode%CrossSection%PoolDepth * CurrNode%CrossSection%BottomWidth * CurrNode%Length 
                CurrNode%VerticalArea = CurrNode%SingCoef * CurrNode%VerticalArea
                !CurrNode%WetPerimeter = CurrNode%SingCoef * CurrNode%WetPerimeter
                CurrNode%VolumeNew    = CurrNode%VerticalArea * CurrNode%Length + PoolVolume
            endif

                
        end if

    end subroutine ComputeXSFromWaterDepth
    
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine TrapezoidGeometry (b, mL, mR, h, Av, P, Sw) 


        !Arguments--------------------------------------------------------------
        real, intent(in)                        :: b, mL, mR, h
        real, intent(out)                       :: Av, P, Sw
        real                                    :: m
        
        m = 0.5 * ( abs(mL) + mR )
        Av = (b + m  * h) * h
        P  =  b + h  * ( sqrt(1. + mL**2.) + sqrt(1. + mR**2.) )
        Sw =  b + 2  * m * h

    end subroutine TrapezoidGeometry

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine TabularGeometry (CrossSection, WaterLevel, Av, P, Sw)

        !Arguments--------------------------------------------------------------
        type(T_CrossSection)                :: CrossSection
        real, intent(in)                    :: WaterLevel
        real, intent(out)                   :: Av, P, Sw

        !Locals----------------------------------------------------------------
        integer                             :: i, ilev
        real                                :: dH, b, m

        if (WaterLevel < CrossSection%BottomLevel)  stop 'WaterLevel lower than bottom level'
        !if (WaterLevel > CrossSection%Elevation(1)) write(*,*) 'WaterLevel higher than max elevation'
        
        Av = 0.
        P  = 0.
        Sw = 0.
        dH = null_real

        do i=  1, CrossSection%NLevels
            if (CrossSection%Level(i) <= WaterLevel) then
                dH = WaterLevel - CrossSection%Level(i)
                ilev = i                
                !exit nao porque quero o lowest level mais aproximado
            endif
        enddo   

        if (dH <= 1e-6) then
            Av = CrossSection%LevelVerticalArea(ilev)
            P  = CrossSection%LevelWetPerimeter(ilev)
            Sw = CrossSection%LevelSurfaceWidth(ilev)
        else        

            b = CrossSection%LevelBottomWidth(ilev)
            m = 0.5 * ( abs(CrossSection%LevelSlopeLeft(ilev)) + CrossSection%LevelSlopeRight(ilev) )
        
            call TrapezoidGeometry (b  = CrossSection%LevelBottomWidth(ilev),   &
                                    mL = CrossSection%LevelSlopeLeft(ilev),     &
                                    mR = CrossSection%LevelSlopeRight(ilev),    &
                                    h  = dH,                                    &
                                    Av = Av,                                    &
                                    P  = P,                                     &
                                    Sw = Sw)
            Av = Av + CrossSection%LevelVerticalArea(ilev)

            if (ilev >= 2) &
                    P  = P + CrossSection%LevelWetPerimeter(ilev) - CrossSection%LevelBottomWidth(ilev)
            
        endif

    end subroutine TabularGeometry

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine InitializeReaches

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Reach), pointer                 :: CurrReach      
        integer                                 :: ReachID

        do ReachID = 1, Me%TotalReaches
            CurrReach => Me%Reaches(ReachID)                
            call UpdateReachCrossSection (CurrReach)
        end do
        
        if (.not. Me%Continuous) then
    
!            if (Me%HydrodynamicApproximation == DiffusionWave) then

!                do ReachID = 1, Me%TotalReaches
!                    CurrReach => Me%Reaches(ReachID)                
!
!                    !Update Slope based on water level
!                    CurrReach%Slope = (Me%Nodes(CurrReach%UpstreamNode)%Waterlevel -            &
!                                       Me%Nodes(CurrReach%DownstreamNode)%Waterlevel) /         &
!                                       CurrReach%Length
!                                       
!                    !Don't allow negative slopes and impose a minimum slope..
!                    if (.not. Me%AllowBackwardWater) then
!                        CurrReach%Slope = max(CurrReach%Slope, Me%MinimumSlope)
!                    endif
!
!                    call ComputeKinematicWave (CurrReach) 
!                end do
!
!            else

                do ReachID = 1, Me%TotalReaches
                    Me%Reaches(ReachID)%FlowNew = 0.0                
                end do
        
!            end if

        end if

    end subroutine InitializeReaches

    !---------------------------------------------------------------------------

    subroutine ConstructStormWaterModelLink

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: NodeID, iout, iin, iUp
        type (T_Reach), pointer                     :: DownStreamReach
        type (T_Node), pointer                      :: CurrNode
        logical                                     :: upStreamActive 

        !Counts nodes
        Me%StormWaterModelLink%nOutflowNodes = 0
        Me%StormWaterModelLink%nInflowNodes  = 0
        

        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            if (CurrNode%nDownStreamReaches == 1) then
                DownStreamReach => Me%Reaches (CurrNode%DownstreamReaches (1))
                
                upStreamActive = .false.
                do iUp = 1, CurrNode%nUpStreamReaches
                    if (Me%Reaches (CurrNode%UpstreamReaches(iUp))%Active) then
                        upStreamActive = .true.
                    endif
                enddo

                !Flow to Storm Water System -> In nodes where downstream node is inactive and upstream nodes active
                if (.not. DownStreamReach%Active .and. upStreamActive) then
                    Me%StormWaterModelLink%nOutflowNodes = Me%StormWaterModelLink%nOutflowNodes + 1
                endif
                
                !Flow from Storm Water System -> In nodes where upstream nodes are inactive and downstream nodes are active
                if (DownStreamReach%Active .and. .not. upStreamActive .and. CurrNode%nUpStreamReaches > 0) then
                    Me%StormWaterModelLink%nInflowNodes = Me%StormWaterModelLink%nInflowNodes + 1
                endif
               
            endif
        enddo
        
        !Allocates Matrixes
        if (Me%StormWaterModelLink%nOutflowNodes > 0) then
            allocate (Me%StormWaterModelLink%OutflowIDs (Me%StormWaterModelLink%nOutflowNodes))
            allocate (Me%StormWaterModelLink%Outflow    (Me%StormWaterModelLink%nOutflowNodes))
        endif
        
        if (Me%StormWaterModelLink%nInflowNodes > 0) then
            allocate (Me%StormWaterModelLink%InflowIDs  (Me%StormWaterModelLink%nInflowNodes))
            allocate (Me%StormWaterModelLink%Inflow     (Me%StormWaterModelLink%nInflowNodes))
        endif
        
        !Fills Matrixes
        iout = 1
        iin  = 1
        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            if (CurrNode%nDownStreamReaches == 1) then
                DownStreamReach => Me%Reaches (CurrNode%DownstreamReaches (1))
                
                upStreamActive = .false.
                do iUp = 1, CurrNode%nUpStreamReaches
                    if (Me%Reaches (CurrNode%UpstreamReaches(iUp))%Active) then
                        upStreamActive = .true.
                    endif
                enddo

                !Flow to Storm Water System -> In nodes where downstream node is inactive and upstream nodes active
                if (.not. DownStreamReach%Active .and. upStreamActive) then
                    Me%StormWaterModelLink%OutflowIDs(iout) = NodeID
                    iout = iout + 1
                endif
                
                !Flow from Storm Water System -> In nodes where upstream nodes are inactive and downstream nodes are active
                if (DownStreamReach%Active .and. .not. upStreamActive .and. CurrNode%nUpStreamReaches > 0) then
                    Me%StormWaterModelLink%InflowIDs(iout) = NodeID
                    iin = iin + 1
                endif
               
            endif
        enddo

    end subroutine ConstructStormWaterModelLink

    !---------------------------------------------------------------------------

    subroutine ConstructSubModules

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: nWaterQualityModels
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: STAT_CALL

        !Begin------------------------------------------------------------------

        !If needed allocate TopRadiation
        if (Me%ComputeOptions%TopRadiation .or. Me%ComputeOptions%T90_Decay)  then
            allocate (Me%TopRadiation       (Me%TotalNodes))
            Me%TopRadiation         = null_real
        endif

        !If Needed allocate AirTemperature and SedimentTemperture
        call SearchProperty (PropertyX, PropertyXIDNumber = Temperature_, STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            allocate (Me%AirTemperature     (Me%TotalNodes))
            allocate (Me%SedimentTemperature(Me%TotalNodes))
            Me%AirTemperature      = null_real
            Me%SedimentTemperature = PropertyX%Concentration
        else
            if (Me%ComputeOptions%TopRadiation .or. Me%ComputeOptions%T90_Decay) then
                write (*,*)'Drainage Network needs Property Temperature'
                stop 'ConstructSubModules - ModuleDrainageNetwork - ERR00'
            endif
        endif


        !If neeed constructs Light Extinctions
        if  (Me%ComputeOptions%TopRadiation)  then
            call CoupleLightExtinction
        endif

      
        !Verifies number of water quality models
        nWaterQualityModels = 0
        if(Me%ComputeOptions%WaterQuality) nWaterQualityModels = nWaterQualityModels + 1
        if(Me%ComputeOptions%CeQualW2    ) nWaterQualityModels = nWaterQualityModels + 1
        if(Me%ComputeOptions%Life        ) nWaterQualityModels = nWaterQualityModels + 1

        if (nWaterQualityModels > 1) then
            write(*,*)'Cannot run more then one Water Quality model in the same simulation'
            stop      'ConstructSubModules - ModuleDrainageNetwork - ERR01'
        end if

        if (Me%ComputeOptions%WaterQuality) then
            call CoupleWaterQuality
        endif

        if (Me%ComputeOptions%Benthos) then
            call CoupleBenthos
        endif

        if (Me%ComputeOptions%BottomFluxes) then
            call SearchProperty(PropertyX, PropertyXIDNumber = TSS_, STAT = STAT_CALL)
            !give a warning, if TSS is not chosen as property
            !cohesive sediment concentration must be taken in this case instead!!!
            if (STAT_CALL /= SUCCESS_) then
              write(*,*) 'Bottom Fluxes are activated, but TSS is not chosen as property'
              write(*,*) 'Cohesive sediment will be taken to calculate erosion rates!'
            endif             
            if (Me%ComputeOptions%CalcFractionSediment)then
                call SearchProperty(PropertyX, PropertyXIDNumber = COHSED_FINE_, STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    allocate (Me%ShearStress (Me%TotalReaches))            
                    Me%ShearStress = 0.0
                else
                    call SearchProperty(PropertyX, PropertyXIDNumber = COHSED_MEDIUM_, STAT = STAT_CALL)
                    if (STAT_CALL == SUCCESS_) then
                        allocate (Me%ShearStress (Me%TotalReaches))            
                        Me%ShearStress = 0.0
                    else
                        call SearchProperty(PropertyX, PropertyXIDNumber = COHSED_COARSE_, STAT = STAT_CALL)
                        if (STAT_CALL == SUCCESS_) then
                            allocate (Me%ShearStress (Me%TotalReaches))            
                            Me%ShearStress = 0.0
                        else
                            write (*,*)
                            write (*,*) 'Bottom Fluxes needs at least one Cohesive Sediment Fraction'
                            stop        'ConstructSubModules - ModuleDrainageNetwork - ERR02_Wassim'
                        end if
                    end if
                end if
            else
                call SearchProperty(PropertyX, PropertyXIDNumber = Cohesive_Sediment_, STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    allocate (Me%ShearStress (Me%TotalReaches))            
                    Me%ShearStress = 0.0
                else
                    write (*,*)
                    write (*,*) 'Bottom Fluxes needs Cohesive_Sediment_'
                    stop        'ConstructSubModules - ModuleDrainageNetwork - ERR02'
                end if
            end if  
        end if

               

    end subroutine ConstructSubModules

    !---------------------------------------------------------------------------
    
    subroutine CoupleLightExtinction

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Size1D)                             :: Size1D
        logical                                     :: NeedsPhyto      = .false.
        logical                                     :: NeedsSPM        = .false.
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: STAT_CALL

        Size1D%ILB = 1
        Size1D%IUB = Me%TotalReaches
        call ConstructLightExtinction(LightExtinctionID = Me%ObjLightExtinction,        &
                                      TimeID            = Me%ObjTime,                   &
                                      EnterDataID       = Me%ObjEnterData,              &
                                      Size1D            = Size1D,                       &
                                      STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleLightExtinction - ModuleDrainageNetwork - ERR01'


        !Allocates Variables - Local Var
        allocate (Me%ShortWaveExtinction(Me%TotalNodes))
        allocate (Me%ShortWaveField     (Me%TotalNodes))
        allocate (Me%LongWaveField      (Me%TotalNodes))
        allocate (Me%NodesDWZ           (Me%TotalNodes))

        !Allocates Variables - External 
        allocate (Me%CloudCover         (Me%TotalNodes))
        allocate (Me%RelativeHumidity   (Me%TotalNodes))
        allocate (Me%WindSpeed          (Me%TotalNodes))

        Me%ShortWaveExtinction  = null_real
        Me%ShortWaveField       = null_real
        Me%LongWaveField        = null_real
        Me%NodesDWZ             = null_real

        Me%CloudCover           = null_real
        Me%RelativeHumidity     = null_real
        Me%WindSpeed            = null_real

        !
        call GetLightExtinctionOptions(LightExtinctionID = Me%ObjLightExtinction,       &
                                       NeedsPhyto        = NeedsPhyto,                  &
                                       NeedsSPM          = NeedsSPM,                    &
                                       STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'CoupleLightExtinction - ModuleDrainageNetwork - ERR02'

        if (NeedsPhyto) then

            call SearchProperty(PropertyX, PropertyXIDNumber = Phytoplankton_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CoupleLightExtinction - ModuleDrainageNetwork - ERR03'
        
            PropertyX%ComputeOptions%LightExtinction = .true.

        end if

        if (NeedsSPM) then

            call SearchProperty(PropertyX, PropertyXIDNumber = Cohesive_Sediment_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CoupleLightExtinction - ModuleDrainageNetwork - ERR04'
        
            PropertyX%ComputeOptions%LightExtinction = .true.

        end if

     
    end subroutine CoupleLightExtinction

    !---------------------------------------------------------------------------

    subroutine CoupleWaterQuality

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        integer, pointer, dimension(:)                      :: WaterQualityPropertyList
        integer                                             :: STAT_CALL
        real                                                :: WaterQualityDT
        integer                                             :: nProp = 0 
        type (T_Size1D)                                     :: Size1D

        !Begin------------------------------------------------------------------

        !Counts the number of Properties which has WaterQuality option set to true
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%ComputeOptions%WaterQuality) then
                nProp = nProp + 1
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates Array wto hold IDs
        allocate (WaterQualityPropertyList(1:nProp))

        !Fills Array
        PropertyX => Me%FirstProperty
        nProp = 0
        do while (associated(PropertyX))
            if (PropertyX%ComputeOptions%WaterQuality) then
                nProp = nProp + 1
                WaterQualityPropertyList(nProp) = PropertyX%ID%IDNumber
            endif
            PropertyX => PropertyX%Next
        enddo

        !Start Interface
        Size1D%ILB = 1
        Size1D%IUB = Me%TotalReaches
        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = WaterQualityModel,             &
                                DT                  = WaterQualityDT,                &
                                PropertiesList      = WaterQualityPropertyList,      &
                                RiverPoints1D       = Me%RiverPoints,                &
                                Size1D              = Size1D,                        &
                                STAT = STAT_CALL)

        deallocate (WaterQualityPropertyList)

        Me%Coupled%WQM%DT_Compute  = WaterQualityDT 
        Me%Coupled%WQM%NextCompute = Me%CurrentTime

    end subroutine CoupleWaterQuality
    
    !---------------------------------------------------------------------------
    
    subroutine CoupleBenthos

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        integer, pointer, dimension(:)                      :: BenthosPropertyList
        integer                                             :: STAT_CALL
        real                                                :: BenthosDT
        integer                                             :: nProp = 0 
        type (T_Size1D)                                     :: Size1D

        !Begin------------------------------------------------------------------

        !Counts the number of Properties which has WaterQuality option set to true
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%ComputeOptions%Benthos) then
                nProp = nProp + 1
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates Array wto hold IDs
        allocate (BenthosPropertyList(1:nProp))

        !Fills Array
        PropertyX => Me%FirstProperty
        nProp = 0
        do while (associated(PropertyX))
            if (PropertyX%ComputeOptions%Benthos) then
                nProp = nProp + 1
                BenthosPropertyList(nProp) = PropertyX%ID%IDNumber
            endif
            PropertyX => PropertyX%Next
        enddo

        !Start Interface
        Size1D%ILB = 1
        Size1D%IUB = Me%TotalReaches
        call ConstructInterface(InterfaceID         = Me%ObjBenthicInterface,           &
                                TimeID              = Me%ObjTime,                       &
                                SinksSourcesModel   = BenthosModel,                     &
                                DT                  = BenthosDT,                        &
                                PropertiesList      = BenthosPropertyList,              &
                                RiverPoints1D       = Me%RiverPoints,                   &
                                Size1D              = Size1D,                           &
                                STAT = STAT_CALL)

        deallocate (BenthosPropertyList)

        Me%Coupled%Benthos%DT_Compute  = BenthosDT 
        Me%Coupled%Benthos%NextCompute = Me%CurrentTime

    end subroutine CoupleBenthos
        
    !---------------------------------------------------------------------------

    subroutine ConstructOutput

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------

        !Verifies in whichs nodes to write time series 
        call ReadTimeSerieNodeList
       
        !Basic Time Series Data (Hydrodynamic Properties) & Transported Properties
        call ConstructTimeSerieList

        !Opens all Time Series Data Files
        call ConstructTimeSeries
        
        !Opens HDF5 File
        if (Me%Output%Yes) then
            call ConstructHDF5Output 
        endif
        
       
    end subroutine ConstructOutput

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ReadTimeSerieNodeList

        !This subroutine reads selected NodeID and selects also the downstream reach.
        !Local------------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound, Found
        integer                                     :: FirstLine, LastLine            
        integer                                     :: STAT_CALL, flag
        integer                                     :: line, NodeID
        integer                                     :: NodePos, DownReachPos
        type (T_Node), pointer                      :: CurrNode
        character (Len = StringLength)              :: aux1, aux2               
        
 
        !-----------------------------------------------------------------------

        !Keyword to see if the user wants the time series to be written by nodes, i.e.,
        !One file per node, with all variables in the headers list
        !if FALSE, its one file per variable with nodes in the headers.
        call GetData(Me%TimeSerie%ByNodes,                                  &
                     Me%ObjEnterData, flag,                                 &
                     SearchType   = FromFile,                               &
                     keyword      = 'TIME_SERIE_BY_NODES',                  &
                     default      = .false.,                                &
                     ClientModule = 'ModuleDrainageNetwork',                &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'ReadTimeSerieNodeList - ModuleDrainageNetwork - ERR01'

    
        call GetData(Me%TimeSerie%Location,                                 &
                     Me%ObjEnterData, flag,                                 &
                     SearchType   = FromFile,                               &
                     keyword      = 'TIME_SERIE_LOCATION',                  &
                     ClientModule = 'ModuleDrainageNetwork',                &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'ReadTimeSerieNodeList - ModuleDrainageNetwork - ERR02'
            
if1:    if (flag==1) then   
     
            Me%TimeSerie%nNodes =  0

            call ConstructEnterData (Me%TimeSerie%ObjEnterData,             &
                                     Me%TimeSerie%Location, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                      &
                stop 'ReadTimeSerieNodeList - ModuleDrainageNetwork - ERR03'
        
            !Nodes--------------------------------------------------------------
            call ExtractBlockFromBuffer(Me%TimeSerie%ObjEnterData, ClientNumber,    &
                                        BeginNodeTimeSerie, EndNodeTimeSerie,       &
                                        BlockFound, FirstLine, LastLine, STAT_CALL) 

            if (STAT_CALL .EQ. SUCCESS_) then    

                if (BlockFound) then                 
                    
                    line = FirstLine
                    do while (line < LastLine - 1)

                        line = line + 1                    

                        call GetData(NodeID, Me%TimeSerie%ObjEnterData,             &
                                     flag, Buffer_Line  = line,                     &
                                     STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                  &
                            stop 'ReadTimeSerieNodeList - ModuleDrainageNetwork - ERR04'

                        call FindNodePosition (NodeID, NodePos, Found)

                        if (.NOT.Found) then
                            write (*,*) 'Node not found'
                            write (*,*) 'Node ID = ', NodeID                           
                            stop 'ReadTimeSerieNodeList - ModuleDrainageNetwork - ERR05'
                        end if

                        
                        CurrNode => Me%Nodes(NodePos)
                        
                        if (CurrNode%TimeSerie) then
                            write (*,*) 'Repeated node in time series: ', CurrNode%ID
                            stop 'ReadTimeSerieNodeList - ModuleDrainageNetwork - ERR05a'
                        end if
                            
                        CurrNode%TimeSerie = .TRUE.
                        Me%TimeSerie%nNodes = Me%TimeSerie%nNodes + 1

                        if (CurrNode%nDownstreamReaches /= 0) then
                            DownReachPos = CurrNode%DownstreamReaches (1)
                            Me%Reaches(DownReachPos)%TimeSerie = .TRUE.
                        else
                            write(aux1,*) NodeID
                            aux2 = 'Requested Node Time Series of Ghost Node ID = '// trim(adjustl(adjustr(aux1)))
                            call SetError(FATAL_, KEYWORD_, aux2, ON)
                        end if


                    end do
                                      
                else 

                    call Block_Unlock(Me%TimeSerie%ObjEnterData, ClientNumber)
                    call RewindBlock (Me%TimeSerie%ObjEnterData, ClientNumber)

                end if

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then 

                stop 'ReadTimeSerieNodeList - ModuleDrainageNetwork - ERR06'

            end if
        
        end if if1


    end subroutine ReadTimeSerieNodeList

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ConstructTimeSerieList

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer                          :: Property
        integer                                             :: i, NodePos
        character(LEN = StringLength)                       :: aux

        
             
        Me%TimeSerie%nProp = BaseTimeSeries

        if (Me%HasProperties) then
            Property => Me%FirstProperty
            do while (associated (Property))            
                if (Property%ComputeOptions%TimeSerie   ) then
                    Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1
                    if (Property%ComputeOptions%Toxicity    ) Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1
                    !deposited
                    !Rosa - erosion_rate + deposition_rate are not yet implemented because they have to be integrated 
                    !over de time serie DT interval
                    if (Property%ComputeOptions%BottomFluxes) Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1
                    if (Property%ComputeOptions%SumTotalConc) Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1
                    if (Property%ComputeOptions%ComputeLoad) Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1
                end if
                Property => Property%Next    
            end do
        end if

        !Shear Stress
        if (Me%ComputeOptions%BottomFluxes  ) Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1

        !Global toxicity
        if (Me%ComputeOptions%Toxicity      ) Me%TimeSerie%nProp = Me%TimeSerie%nProp + 1

        !Output hydrodynamic properties
        if (Me%OutputHydro                  ) Me%TimeSerie%nProp = Me%TimeSerie%nProp + 5


if0:    if (Me%TimeSerie%ByNodes) then


            allocate (Me%TimeSerie%ObjTimeSerie (1:Me%TimeSerie%nNodes))
            allocate (Me%TimeSerie%Name         (1:Me%TimeSerie%nNodes))
            allocate (Me%TimeSerie%X            (1:Me%TimeSerie%nNodes))
            allocate (Me%TimeSerie%Y            (1:Me%TimeSerie%nNodes))

            Me%TimeSerie%ObjTimeSerie = 0
            i = 0

            do NodePos = 1, Me%TotalNodes
                                
                if (Me%Nodes(NodePos)%TimeSerie) then
                    
                    i = i + 1                    
                    aux = ''
                    write(aux,'(i10)') Me%Nodes(NodePos)%ID                    
                    Me%TimeSerie%Name(i)  = 'Node_'//trim(adjustl(adjustr(aux)))
                    Me%TimeSerie%X(i)     = Me%Nodes(NodePos)%X
                    Me%TimeSerie%Y(i)     = Me%Nodes(NodePos)%Y
                end if

            end do


        else if0

            
            allocate (Me%TimeSerie%ObjTimeSerie (1:Me%TimeSerie%nProp))
            allocate (Me%TimeSerie%Name         (1:Me%TimeSerie%nProp))

            Me%TimeSerie%ObjTimeSerie = 0

            call FillPropNameVector (Me%TimeSerie%Name)

        end if if0

    end subroutine ConstructTimeSerieList

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ConstructTimeSeries

        !Local------------------------------------------------------------------
        integer                                              :: STAT_CALL
        integer                                              :: NodePos, ReachPos        
        integer                                              :: nNodes, nReaches, i
        character(LEN = StringLength)                        :: aux
        character(LEN = StringLength), dimension(:), pointer :: NodeHeaderList 
        character(LEN = StringLength), dimension(:), pointer :: ReachHeaderList 
        character(LEN = StringLength), dimension(:), pointer :: PropHeaderList 
       
       
        
if1:    if (Me%TimeSerie%nNodes > 0) then

if2:        if (Me%TimeSerie%ByNodes) then
            
            
                allocate (PropHeaderList           (1:Me%TimeSerie%nProp ))
                allocate (Me%TimeSerie%DataLine    (1:Me%TimeSerie%nProp ))            

                call FillPropNameVector (PropHeaderList)

                do i = 1, Me%TimeSerie%nNodes

                    call StartTimeSerie(Me%TimeSerie%ObjTimeSerie(i), Me%ObjTime,       &
                                        TimeSerieDataFile = trim(Me%TimeSerie%Location),&
                                        PropertyList      = PropHeaderList,             &
                                        Extension         = "srn",                      &
                                        ResultFileName    = Me%TimeSerie%Name (i),      &
                                        ModelName         = Me%ModelName,               &
                                        CoordX            = Me%TimeSerie%X    (i),      &
                                        CoordY            = Me%TimeSerie%Y    (i),      &
                                        STAT              = STAT_CALL)
                    if (STAT_CALL /= 0) stop 'ConstructTimeSeries - ModuleDrainageNetwork - ERR01'
                
                end do

                deallocate (PropHeaderList)
            
            else !if2
                        
                allocate (NodeHeaderList         (1:Me%TimeSerie%nNodes ))
                allocate (ReachHeaderList        (1:Me%TimeSerie%nNodes ))            
                allocate (Me%TimeSerie%DataLine  (1:Me%TimeSerie%nNodes ))            

                !NodeHeaderList-------------------------------------------------
                nNodes = 1
                do NodePos = 1, Me%TotalNodes
                
                    if (Me%Nodes(NodePos)%TimeSerie) then
                        aux = ''
                        write(aux,'(i10)') Me%Nodes(NodePos)%ID                    
                        NodeHeaderList(nNodes) = 'Node_'//trim(adjustl(adjustr(aux)))
                        nNodes = nNodes + 1 
                    endif                
                enddo
            
                !ReachHeaderList------------------------------------------------
                nReaches = 1
                do ReachPos = 1, Me%TotalReaches
                
                    if (Me%Reaches(ReachPos)%TimeSerie) then
                        ReachHeaderList(nReaches) = Me%Reaches(ReachPos)%Name
                        nReaches = nReaches + 1 
                    endif                

                enddo

                !BaseTimeSeries + PropertiesTimeSeries--------------------------
                do i = 1, Me%TimeSerie%nProp


                    if (Me%TimeSerie%Name (i) == Char_Flow              .or. &
                        Me%TimeSerie%Name (i) == Char_Velocity          .or. &
                        Me%TimeSerie%Name (i) == Char_HydroTimeGradient .or. &
                        Me%TimeSerie%Name (i) == Char_HydroAdvection    .or. &
                        Me%TimeSerie%Name (i) == Char_HydroPressure     .or. &
                        Me%TimeSerie%Name (i) == Char_HydroGravity      .or. &
                        Me%TimeSerie%Name (i) == Char_HydroFriction ) then


                        call StartTimeSerie(Me%TimeSerie%ObjTimeSerie(i), Me%ObjTime,           &
                                            TimeSerieDataFile = trim(Me%TimeSerie%Location),    &
                                            PropertyList      = ReachHeaderList,                &
                                            Extension         = "srn",                          &
                                            ResultFileName    = Me%TimeSerie%Name (i),          &
                                            ModelName         = Me%ModelName,                   &                
                                            STAT              = STAT_CALL)
                        if (STAT_CALL /= 0) stop 'ConstructTimeSeries - ModuleDrainageNetwork - ERR02'

                    else

                        call StartTimeSerie(Me%TimeSerie%ObjTimeSerie(i), Me%ObjTime,           &
                                            TimeSerieDataFile = trim(Me%TimeSerie%Location),    &
                                            PropertyList      = NodeHeaderList,                 &
                                            Extension         = "srn",                          &
                                            ResultFileName    = Me%TimeSerie%Name (i),          &
                                            ModelName         = Me%ModelName,                   &                
                                            STAT              = STAT_CALL)
                        if (STAT_CALL /= 0) stop 'ConstructTimeSeries - ModuleDrainageNetwork - ERR03'

                    end if

                end do    
                                
                deallocate (NodeHeaderList  )
                deallocate (ReachHeaderList )

            end if if2

        end if if1


    end subroutine ConstructTimeSeries
    
    !---------------------------------------------------------------------------
    
    subroutine FillPropNameVector (PropVector)

        !Arguments--------------------------------------------------------------
        character(StringLength) , dimension (:), pointer    :: PropVector  

        !Local------------------------------------------------------------------
        type (T_Property), pointer                          :: Property
        integer                                             :: i


        PropVector (pWaterDepth         ) = Char_WaterDepth
        PropVector (pWaterLevel         ) = Char_WaterLevel
        PropVector (pVolume             ) = Char_Volume
        PropVector (pPercentageMaxVolume) = Char_PercentageMaxVolume
        PropVector (pFlow               ) = Char_Flow
        PropVector (pVelocity           ) = Char_Velocity
        PropVector (pVerticalArea       ) = Char_VerticalArea
        PropVector (pFlowToChannels     ) = Char_FlowToChannels
        PropVector (pGWFlowToChannels   ) = Char_GWFlowToChannels
        PropVector (pPoolDepth          ) = Char_PoolDepth
        PropVector (pDT                 ) = Char_DT
        PropVector (pDTLocal            ) = Char_DTLocal
        i = BaseTimeSeries

        if (Me%OutputHydro) then

            PropVector (pHydroTimeGradient  ) = Char_HydroTimeGradient
            PropVector (pHydroAdvection     ) = Char_HydroAdvection
            PropVector (pHydroPressure      ) = Char_HydroPressure
            PropVector (pHydroGravity       ) = Char_HydroGravity
            PropVector (pHydroFriction      ) = Char_HydroFriction
            i = i + 5

        end if

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

                    if (Property%ComputeOptions%BottomFluxes) then
                        if (Property%OutputName == 'NAME') then
                            PropVector (i) = 'Bottom_'//trim(adjustl(adjustr(Property%ID%Name)))
                        else if (Property%OutputName == 'DESCRIPTION') then
                            PropVector (i) = 'Bottom_'//trim(adjustl(adjustr(Property%ID%Description)))
                        end if
                        i = i + 1
                        
                        if (Property%ComputeOptions%SumTotalConc) then
                            if (Property%OutputName == 'NAME') then
                                PropVector (i) = 'Total_'//trim(adjustl(adjustr(Property%ID%Name)))
                            else if (Property%OutputName == 'DESCRIPTION') then
                                PropVector (i) = 'Total_'//trim(adjustl(adjustr(Property%ID%Description)))
                            end if
                            i = i + 1
                        end if
                        
                    end if

                    if (Property%ComputeOptions%Toxicity) then
                        if (Property%OutputName == 'NAME') then
                            PropVector (i) = trim(adjustl(adjustr(Property%ID%Name)))//'_toxicity'
                        else if (Property%OutputName == 'DESCRIPTION') then
                            PropVector (i) = trim(adjustl(adjustr(Property%ID%Description)))//'_toxicity'
                        end if
                        i = i + 1              
                    end if

                    if (Property%ComputeOptions%ComputeLoad) then
                        if (Property%OutputName == 'NAME') then
                            PropVector (i) = 'Load_'//trim(adjustl(adjustr(Property%ID%Name)))
                        else if (Property%OutputName == 'DESCRIPTION') then
                            PropVector (i) = 'Load_'//trim(adjustl(adjustr(Property%ID%Description)))
                        end if
                        i = i + 1              
                    end if

                end if                        
                    
                Property => Property%Next
            end do

            if (Me%ComputeOptions%BottomFluxes) then
                PropVector (i) = 'shear_stress'   
                i = i + 1
            end if
             
            if (Me%ComputeOptions%Toxicity    ) then
                PropVector (i) = 'global_toxicity'    
                i = i + 1
            end if
        
        end if if0
        
        nullify (Property) 
        

    end subroutine FillPropNameVector

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Output

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                             :: HDF5_CREATE, STAT_CALL
        integer                                             :: iNode, ReachID
        real, dimension(:), pointer                         :: NodeX, NodeY
        integer, dimension(:), pointer                      :: NodeID, ReachIDs, UpNode, DownNode
        
        
        Me%OutPut%NextOutPut = 1  
        
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%HDFFile)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleDrainageNetwork - ERR01'
        
        !Writes information about the Nodes / Reaches
        !Writes the Nodes X / Y
        allocate(NodeID(1: Me%TotalNodes))
        allocate(NodeX (1: Me%TotalNodes))
        allocate(NodeY (1: Me%TotalNodes))
        
        do iNode = 1, Me%TotalNodes
            NodeID(iNode) = Me%Nodes(iNode)%ID
            NodeX(iNode)  = Me%Nodes(iNode)%X
            NodeY(iNode)  = Me%Nodes(iNode)%Y
        enddo

        allocate(UpNode  (1: Me%TotalReaches))
        allocate(DownNode(1: Me%TotalReaches))
        allocate(ReachIDs(1: Me%TotalReaches))
        
        do ReachID = 1, Me%TotalReaches
            ReachIDs(ReachID) = ReachID
            UpNode  (ReachID) = Me%Nodes(Me%Reaches(ReachID)%UpstreamNode)%ID
            DownNode(ReachID) = Me%Nodes(Me%Reaches(ReachID)%DownstreamNode)%ID
        enddo

        !Nodes
        call HDF5SetLimits  (Me%ObjHDF5, 1, Me%TotalNodes, STAT = STAT_CALL)

        call HDF5WriteData   (Me%ObjHDF5, "/Nodes", "ID", "m",                          &
                              Array1D = NodeID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleDrainageNetwork - ERR02'
        
        call HDF5WriteData   (Me%ObjHDF5, "/Nodes", "X", "m",                           &
                              Array1D = NodeX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleDrainageNetwork - ERR02'

        call HDF5WriteData   (Me%ObjHDF5, "/Nodes", "Y", "m",                           &
                              Array1D = NodeY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleDrainageNetwork - ERR03'

        !Reaches
        call HDF5SetLimits  (Me%ObjHDF5, 1, Me%TotalReaches, STAT = STAT_CALL)
        
        call HDF5WriteData   (Me%ObjHDF5, "/Reaches", "ID", "-",                         &
                              Array1D = ReachIDs, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleDrainageNetwork - ERR04'


        call HDF5WriteData   (Me%ObjHDF5, "/Reaches", "Up", "-",                         &
                              Array1D = UpNode, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleDrainageNetwork - ERR04'

        call HDF5WriteData   (Me%ObjHDF5, "/Reaches", "Down", "-",                       &
                              Array1D = DownNode, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleDrainageNetwork - ERR05'


        deallocate(NodeX, NodeY)
        deallocate(DownNode, UpNode)
        
    end subroutine ConstructHDF5Output

    !---------------------------------------------------------------------------

    subroutine ConstructLog

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: Property


#ifndef _OPENMI_
        write(*, *)
        write(*, *)"-------------------- DRAINAGE NETWORK --------------------"         
        write(*, *)

        !Writes General Stats
        write(*, *)"Number of Nodes      : ", Me%TotalNodes
        write(*, *)"Number of Reaches    : ", Me%TotalReaches
        write(*, *)"Number of Outlets    : ", Me%TotalOutlets

        !Writes Hydrodynamic Approximation
        select case (Me%HydrodynamicApproximation)
            case (KinematicWave)
                write(*,*) 'Hydrodynamic Approx. : Kinematic wave'
            case (DiffusionWave)
                write(*,*) 'Hydrodynamic Approx. : Diffuse wave'
            case (DynamicWave)
                write(*,*) 'Hydrodynamic Approx. : Dynamic Wave'
            case default 
                write(*,*)'Hydrodynamic Approx. : INVALID'
                stop 'ModuleDrainageNetwork - ConstructLog - ERR01'
        end select

        !Writes NumericalScheme
        if (Me%NumericalScheme == ExplicitScheme) then
            write(*,*) 'Numerical Scheme     : ExplicitScheme'
        else 
            write(*,*) 'Numerical Scheme     : ImplicitScheme'
        endif 

        
        !Writes Downstream Boundary
        select case (Me%Downstream%Boundary)
            case (Dam)
                write(*,*) 'Downstream Boundary  : Dam'
            case(ZeroDepthGradient)
                write(*,*) 'Downstream Boundary  : Zero depth gradient'
            case(CriticalDepth)
                write(*,*) 'Downstream Boundary  : Critical depth'
            case(ImposedWaterLevel)
                write(*,*) 'Downstream Boundary  : Imposed water level'
            case(ImposedVelocity)
                write(*,*) 'Downstream Boundary  : Imposed velocity'
            case default
                write(*, *)"Downstream Boundary  : INVALID"
                stop 'ModuleDrainageNetwork - ConstructLog - ERR02'
        end select

        write(*, *)

        if (Me%HasProperties) then

            write(*, *)"--------------- DRAINAGE NETWORK PROPERTIES --------------"      
            write(*, *)
            write(*, *)"Num of Properties : ", Me%PropertiesNumber
            write(*, *)

            Property => Me%FirstProperty
            do while (associated(Property))
                
                if (Property%OutputName == 'NAME') then
                    write(*, *)"Property                : ", trim(Property%ID%Name)                
                else if (Property%OutputName == 'DESCRIPTION') then
                    write(*, *)"Property                : ", trim(Property%ID%Description)                                
                end if

                    write(*, *)"---Advection Diffusion  : ", Property%ComputeOptions%AdvectionDiffusion
                    write(*, *)"---Discharges           : ", Property%ComputeOptions%Discharges
                    write(*, *)"---Toxicity             : ", Property%ComputeOptions%Toxicity
                    write(*, *)"---SurfaceFluxes        : ", Property%ComputeOptions%SurfaceFluxes
                    write(*, *)"---BottomFluxes         : ", Property%ComputeOptions%BottomFluxes
                    write(*, *)"---WaterQuality         : ", Property%ComputeOptions%WaterQuality                    
                    write(*, *)"---Benthos              : ", Property%ComputeOptions%Benthos   
                    write(*, *)"---GenericDecay         : ", Property%ComputeOptions%Generic_Decay                 
                    write(*, *)"---T90Decay             : ", Property%ComputeOptions%T90_Decay
                    write(*, *)
                
                                
                Property=>Property%Next
            enddo
        end if
#endif  

    end subroutine ConstructLog

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine FindNodePosition (NodeID, NodePos, Found)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                             :: NodeID
        integer, intent(OUT)                            :: NodePos
        logical, intent(OUT)                            :: Found
        !Local------------------------------------------------------------------
        type (T_Node), pointer                          :: CurrNode

        Found = .FALSE.

        do NodePos = 1, Me%TotalNodes
        
            nullify(CurrNode)
            CurrNode => Me%Nodes (NodePos)
            if (CurrNode%ID == NodeID) then
                Found = .TRUE.
                exit
            end if

        end do

    end subroutine FindNodePosition

    !---------------------------------------------------------------------------

    subroutine FindReachPosition (ReachID, ReachPos, Found)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                             :: ReachID
        integer, intent(OUT)                            :: ReachPos
        logical, intent(OUT)                            :: Found
        !Local------------------------------------------------------------------
        type (T_Reach), pointer                         :: CurrReach
        
        Found = .FALSE.

        do ReachPos = 1, Me%TotalReaches
        
            nullify(CurrReach)
            CurrReach => Me%Reaches (ReachPos)
            if (CurrReach%ID == ReachID) then
                Found = .TRUE.
                exit
            end if

        end do

    end subroutine FindReachPosition

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SEL

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !---------------------------------------------------------------------------

    subroutine GetDrainageSize (DrainageNetworkID, nNodes, nReaches, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        integer, intent(OUT), optional                  :: nNodes
        integer, intent(OUT), optional                  :: nReaches
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present (nNodes   )) nNodes    = Me%TotalNodes
            if (present (nReaches )) nReaches  = Me%TotalReaches

            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetDrainageSize

    !---------------------------------------------------------------------------        

    subroutine GetChannelsID (DrainageNetworkID, ChannelsID, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        integer, dimension (:,:), pointer               :: ChannelsID
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsID => Me%ChannelsID                             

            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsID

    !---------------------------------------------------------------------------        

    subroutine GetChannelsStationName (DrainageNetworkID, ChannelsStationName, STAT)

        !Arguments--------------------------------------------------------------
        integer                                                 :: DrainageNetworkID
        character(len=StringLength), dimension (:,:), pointer   :: ChannelsStationName
        integer, intent(OUT), optional                          :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsStationName => Me%ChannelsStationName                             

            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsStationName

    !---------------------------------------------------------------------------

    subroutine GetChannelsWaterLevel (DrainageNetworkID, ChannelsWaterLevel, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension (:,:), pointer                  :: ChannelsWaterLevel
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsWaterLevel => Me%ChannelsWaterLevel
                    
          

            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsWaterLevel

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetChannelsBottomLevel (DrainageNetworkID, ChannelsBottomLevel, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension (:,:), pointer                  :: ChannelsBottomLevel
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsBottomLevel => Me%ChannelsBottomLevel
                    
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsBottomLevel

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetChannelsBankSlope(DrainageNetworkID, ChannelsBankSlope, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension (:,:), pointer                  :: ChannelsBankSlope
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsBankSlope => Me%ChannelsBankSlope
                    
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsBankSlope

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetChannelsSurfaceWidth (DrainageNetworkID, ChannelsSurfaceWidth, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension (:,:), pointer                  :: ChannelsSurfaceWidth
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsSurfaceWidth => Me%ChannelsSurfaceWidth
                    
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsSurfaceWidth

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetChannelsBottomWidth (DrainageNetworkID, ChannelsBottomWidth, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension (:,:), pointer                  :: ChannelsBottomWidth
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsBottomWidth => Me%ChannelsBottomWidth
                    
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsBottomWidth

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetChannelsNodeLength (DrainageNetworkID, ChannelsNodeLength, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension (:,:), pointer                  :: ChannelsNodeLength
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsNodeLength => Me%ChannelsNodeLength
                    
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsNodeLength

    !---------------------------------------------------------------------------

    subroutine GetChannelsOpenProcess (DrainageNetworkID, ChannelsOpenProcess, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        integer, dimension (:,:), pointer               :: ChannelsOpenProcess
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

           
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID)
            ChannelsOpenProcess => Me%ChannelsOpenProcess
                    
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetChannelsOpenProcess

    !---------------------------------------------------------------------------

    subroutine GetHasProperties (DrainageNetworkID, HasProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        logical                                         :: HasProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            HasProperties       = Me%HasProperties
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetHasProperties

    !---------------------------------------------------------------------------

    subroutine GetDNConcentration(DrainageNetworkID, ConcentrationX, PropertyXIDNumber, &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: DrainageNetworkID
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

        call Ready(DrainageNetworkID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mDRAINAGENETWORK_, Me%InstanceID) 

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
                write(*,*) 'Looking for Property in Drainage Network', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'GetDNConcentration - ModuleDrainageNetwork - ERR010'
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetDNConcentration

    !--------------------------------------------------------------------------------

    subroutine GetDNnProperties (DrainageNetworkID, nProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        integer                                         :: nProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            nProperties       = Me%PropertiesNumber
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetDNnProperties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetDNPropertiesIDByIdx (DrainageNetworkID, Idx, ID,PropAdvDiff, Particulate, OutputName, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        integer, intent(IN)                             :: Idx
        integer, intent(OUT)                            :: ID
        logical, intent(OUT)                            :: PropAdvDiff
        logical, intent(OUT), optional                  :: Particulate
        integer, intent(OUT), optional                  :: STAT
        character (Len = StringLength), optional        :: OutputName

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, i
        type (T_Property), pointer                      :: CurrProp

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrProp => Me%FirstProperty
            do i = 1, idx - 1
                CurrProp => CurrProp%Next
            enddo

            ID          = CurrProp%ID%IDNumber
            PropAdvDiff = CurrProp%ComputeOptions%AdvectionDiffusion
            
            if (present (Particulate)) then
                Particulate = Check_Particulate_Property(CurrProp%ID%IDNumber)
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

    end subroutine GetDNPropertiesIDByIdx

    !---------------------------------------------------------------------------
 
     subroutine CheckDNProperty (DrainangeNetworkID,                        &
                                  PropertyID,                               &
                                  STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainangeNetworkID
        integer                                         :: PropertyID
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type (T_Property), pointer                      :: PropertyX
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainangeNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then 
                if (.not. PropertyX%ComputeOptions%AdvectionDiffusion) then
                    write(*,*)
                    write(*,*)'Property', GetPropertyName(PropertyID)
                    write(*,*)'has advection diffusion inactive in Drainage Network Module'
                    write(*,*)'and it is unconsistent with activation in other Modules'
                    stop 'CheckDNProperty - ModuleDraianageNetwork - ERR01' 
                else               
                    STAT_CALL = SUCCESS_
                endif
            else
                write(*,*)
                write(*,*)'Could not find property', GetPropertyName(PropertyID)
                write(*,*)'in DraianageNetwork Module'
                stop 'CheckDNProperty - ModuleDraianageNetwork - ERR010'
            endif
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine CheckDNProperty

    !---------------------------------------------------------------------------
 
    subroutine GetHasToxicity (DrainageNetworkID, HasToxicity, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        logical                                         :: HasToxicity
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            HasToxicity       = Me%ComputeOptions%Toxicity
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetHasToxicity

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetPropHasBottomFluxes (ObjDrainageNetworkID, PropIDNumber, HasBottomFluxes, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjDrainageNetworkID        
        integer, intent (IN )                           :: PropIDNumber
        logical, intent (OUT)                           :: HasBottomFluxes
        integer, intent (OUT), optional                 :: STAT        

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        type (T_Property), pointer                      :: Property

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ObjDrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call SearchProperty (Property, PropIDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - GetPropHasBottomFluxes - ERR01'

            HasBottomFluxes = Property%ComputeOptions%BottomFluxes

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetPropHasBottomFluxes

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    
    subroutine GetNeedsRadiation (DrainageNetworkID, NeedRadiation, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        logical                                         :: NeedRadiation
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_
        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            NeedRadiation       = Me%ComputeOptions%SurfaceFluxes .or. Me%ComputeOptions%WaterQuality .or. &
                                  Me%ComputeOptions%CeQualW2      .or. Me%ComputeOptions%Life         .or. &
                                  Me%ComputeOptions%T90_Decay
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNeedsRadiation
    
    !---------------------------------------------------------------------------

    subroutine GetNeedsAtmosphere (DrainageNetworkID, NeedAtmosphere, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        logical                                         :: NeedAtmosphere
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_
        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            NeedAtmosphere = Me%ComputeOptions%SurfaceFluxes
            STAT_CALL = SUCCESS_

        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNeedsAtmosphere

    !---------------------------------------------------------------------------

    subroutine GetNextDrainageNetDT (DrainageNetworkID, DT, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, intent(OUT)                               :: DT
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            DT        = Me%NextDT

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNextDrainageNetDT

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetVolumes(DrainageNetworkID, TotalInputVolume,                           &
                          TotalOutputVolume, TotalStoredVolume, TotalFlowVolume,         &
                          TotalOvertopVolume, TotalStormWaterOutput,                     &
                          TotalStormWaterInput, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real(8), intent(OUT), optional                  :: TotalInputVolume
        real(8), intent(OUT), optional                  :: TotalOutputVolume
        real(8), intent(OUT), optional                  :: TotalStoredVolume
        real(8), intent(OUT), optional                  :: TotalFlowVolume
        real(8), intent(OUT), optional                  :: TotalOvertopVolume
        real(8), intent(OUT), optional                  :: TotalStormWaterOutput
        real(8), intent(OUT), optional                  :: TotalStormWaterInput
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            if (present (TotalInputVolume      )) TotalInputVolume      = Me%TotalInputVolume
            if (present (TotalOutputVolume     )) TotalOutputVolume     = Me%TotalOutputVolume
            if (present (TotalStoredVolume     )) TotalStoredVolume     = Me%TotalStoredVolume
            if (present (TotalFlowVolume       )) TotalFlowVolume       = Me%TotalFlowVolume
            if (present (TotalOvertopVolume    )) TotalOvertopVolume    = Me%TotalOvertopVolume
            if (present (TotalStormWaterOutput )) TotalStormWaterOutput = Me%TotalStormWaterOutput
            if (present (TotalStormWaterInput  )) TotalStormWaterInput  = Me%TotalStormWaterInput
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL
    
    end subroutine GetVolumes

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetDNMassBalance(DrainageNetworkID, PropertyID, TotalDischargeMass,       &
                          TotalOutFlowMass, TotalStoredMass, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        integer                                         :: PropertyID
        real(8), intent(OUT), optional                  :: TotalDischargeMass
        real(8), intent(OUT), optional                  :: TotalOutFlowMass
        real(8), intent(OUT), optional                  :: TotalStoredMass
!        real(8), intent(OUT), optional                  :: TotalOverTopMass
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, STAT_, ready_
        type (T_Property), pointer                      :: PropertyX
        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then 
                if (present (TotalDischargeMass   )) TotalDischargeMass   = PropertyX%MB%TotalDischargeMass
                if (present (TotalOutFlowMass     )) TotalOutFlowMass     = PropertyX%MB%TotalOutFlowMass
                if (present (TotalStoredMass      )) TotalStoredMass      = PropertyX%MB%TotalStoredMass
!                if (present (TotalOverTopMass     )) TotalOverTopMass     = Property%MB%TotalOverTopMass
                STAT_CALL = SUCCESS_
            else
                stop 'GetDNMassBalance - ModuleRunoffProperties - ERR01'            
            endif
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL
    
    end subroutine GetDNMassBalance

    !---------------------------------------------------------------------------

    subroutine UnGetDrainageNetworkR4 (ObjDrainageNetworkID, Array, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjDrainageNetworkID
        real, dimension(:, :), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDrainageNetworkID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mDRAINAGENETWORK_, Me%InstanceID, "UnGetDrainageNetwork")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetDrainageNetworkR4

    !---------------------------------------------------------------------------

    subroutine UnGetDrainageNetwork1DR4 (ObjDrainageNetworkID, Array, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjDrainageNetworkID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDrainageNetworkID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mDRAINAGENETWORK_, Me%InstanceID, "UnGetDrainageNetwork")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetDrainageNetwork1DR4

    !---------------------------------------------------------------------------    

    subroutine UnGetDrainageNetworkI4 (ObjDrainageNetworkID, Array, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjDrainageNetworkID
        integer, dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDrainageNetworkID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mDRAINAGENETWORK_, Me%InstanceID, "UnGetDrainageNetwork")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetDrainageNetworkI4
    !---------------------------------------------------------------------------

    subroutine UnGetDrainageNetworkA4 (ObjDrainageNetworkID, Array, STAT)

        !Arguments--------------------------------------------------------------
        integer                                                 :: ObjDrainageNetworkID
        character(len=StringLength), dimension(:, :), pointer   :: Array
        integer, intent(OUT), optional                          :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDrainageNetworkID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mDRAINAGENETWORK_, Me%InstanceID, "UnGetDrainageNetwork")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetDrainageNetworkA4

    !---------------------------------------------------------------------------

    subroutine SetAtmosphereRiverNet  (ObjDrainageNetworkID, TopRadiation, AirTemperature,      &
                                           CloudCover, RelativeHumidity, WindSpeed, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjDrainageNetworkID
        real, optional                                  :: TopRadiation
        real, optional                                  :: AirTemperature
        real, optional                                  :: CloudCover
        real, optional                                  :: RelativeHumidity
        real, optional                                  :: WindSpeed
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: NodeID
        type (T_Node), pointer                          :: CurrNode
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDrainageNetworkID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes(NodeID)
                
                if (present(TopRadiation))      Me%TopRadiation     (NodeID) = TopRadiation     
                if (present(AirTemperature))    Me%AirTemperature   (NodeID) = AirTemperature   
                if (present(CloudCover))        Me%CloudCover       (NodeID) = CloudCover           
                if (present(RelativeHumidity))  Me%RelativeHumidity (NodeID) = RelativeHumidity 
                if (present(WindSpeed))         Me%WindSpeed        (NodeID) = WindSpeed        
            enddo

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_

    end subroutine SetAtmosphereRiverNet

    !---------------------------------------------------------------------------

    subroutine SetAtmosphereDrainageNet   (ObjDrainageNetworkID, TopRadiation, AirTemperature,      &
                                           CloudCover, RelativeHumidity, WindSpeed, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjDrainageNetworkID
        real, dimension(:, :), pointer, optional        :: TopRadiation
        real, dimension(:, :), pointer, optional        :: AirTemperature
        real, dimension(:, :), pointer, optional        :: CloudCover
        real, dimension(:, :), pointer, optional        :: RelativeHumidity
        real, dimension(:, :), pointer, optional        :: WindSpeed
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: NodeID
        type (T_Node), pointer                          :: CurrNode
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjDrainageNetworkID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes(NodeID)
                
                if (present(TopRadiation))      Me%TopRadiation     (NodeID) = TopRadiation     (CurrNode%GridI, CurrNode%GridJ)
                if (present(AirTemperature))    Me%AirTemperature   (NodeID) = AirTemperature   (CurrNode%GridI, CurrNode%GridJ)
                if (present(CloudCover))        Me%CloudCover       (NodeID) = CloudCover       (CurrNode%GridI, CurrNode%GridJ)    
                if (present(RelativeHumidity))  Me%RelativeHumidity (NodeID) = RelativeHumidity (CurrNode%GridI, CurrNode%GridJ)
                if (present(WindSpeed))         Me%WindSpeed        (NodeID) = WindSpeed        (CurrNode%GridI, CurrNode%GridJ)
            enddo

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_

    end subroutine SetAtmosphereDrainageNet 

    !---------------------------------------------------------------------------
    
    subroutine SetPMPConcDN   (DrainageNetworkID, ConcentrationX2D, ConcentrationX3D, &
                                                PropertyXIDNumber, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension(:, :), pointer, optional        :: ConcentrationX2D
        real, dimension(:, :, :), pointer, optional     :: ConcentrationX3D
        integer                                         :: PropertyXIDNumber
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: NodeID
        type (T_Node), pointer                          :: CurrNode
        integer                                         :: STAT_, ready_
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
          
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_)

            if (STAT_ == SUCCESS_) then
                
                if (present(ConcentrationX2D)) then
                    do NodeID = 1, Me%TotalNodes
                        CurrNode => Me%Nodes(NodeID)
                        
                        PropertyX%GWaterConc (NodeID) = ConcentrationX2D     (CurrNode%GridI, CurrNode%GridJ)

                    enddo
                elseif (present(ConcentrationX3D)) then
                    PropertyX%GWaterConcLayers => ConcentrationX3D
                endif               

            else
                write(*,*) 'Looking for Porous Media Property in Drainage Network', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetPMPConcDN - ModuleDrainageNetwork - ERR010'
            end if

        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetPMPConcDN

    !---------------------------------------------------------------------------

    subroutine SetRPConcDN   (DrainageNetworkID, ConcentrationX,      &
                                                PropertyXIDNumber, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        real, dimension(:, :), pointer                  :: ConcentrationX
        integer                                         :: PropertyXIDNumber
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: NodeID
        type (T_Node), pointer                          :: CurrNode
        integer                                         :: STAT_, ready_
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_)

            if (STAT_ == SUCCESS_) then
            
                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes(NodeID)
                    
 !                   PropertyX%ConcentrationRP (NodeID) = ConcentrationX     (CurrNode%GridI, CurrNode%GridJ)
                    PropertyX%OverLandConc (NodeID) = ConcentrationX     (CurrNode%GridI, CurrNode%GridJ)

                enddo                

            else
                write(*,*) 'Looking for Runoff Property in Drainage Network', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetPMPConcDN - ModuleDrainageNetwork - ERR010'
            end if

        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetRPConcDN

    !---------------------------------------------------------------------------

    subroutine SetGWFlowLayersToDN   (DrainageNetworkID, GWFlowBottomLayer,        &
                                      GWFlowTopLayer, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: DrainageNetworkID
        integer, dimension(:, :), pointer               :: GWFlowBottomLayer
        integer, dimension(:, :), pointer               :: GWFlowTopLayer
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: NodeID
        type (T_Node), pointer                          :: CurrNode
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
                
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes(NodeID)
                
                Me%GWFlowTopLayer (NodeID)    = GWFlowTopLayer     (CurrNode%GridI, CurrNode%GridJ)
                Me%GWFlowBottomLayer (NodeID) = GWFlowBottomLayer  (CurrNode%GridI, CurrNode%GridJ)

            enddo
            STAT_ = SUCCESS_ 
        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetGWFlowLayersToDN

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
                if (PrintWarning) write (*,*)'Property Not Found in Module DrainageNetwork ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if

        if (present(STAT)) STAT = STAT_

        !-----------------------------------------------------------------------

    end subroutine SearchProperty

   !----------------------------------------------------------------------------

    logical function PropertyExists (PropertyXIDNumber)

        !Arguments--------------------------------------------------------------
        integer, intent (IN)    :: PropertyXIDNumber

        !Local------------------------------------------------------------------
        type(T_Property), pointer         :: PropertyX


        PropertyX => Me%FirstProperty

do2 :   do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber==PropertyXIDNumber) then
                PropertyExists = .true.
                return
            else
                PropertyX => PropertyX%Next                 
            end if
        end do do2

        PropertyExists = .false.
        return

        !-----------------------------------------------------------------------

    end function PropertyExists
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !---------------------------------------------------------------------------

    subroutine ModifyDrainageNetWithGrid(DrainageNetworkID, OLFlowToChannels,            & 
                                         GWFlowToChannels, GWFlowToChannelsLayer, DiffuseFlow, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: DrainageNetworkID          
        real, dimension(:, :), pointer              :: OLFlowToChannels
        real, dimension(:, :), pointer              :: GWFlowToChannels
        real, dimension(:, :, :), pointer, optional :: GWFlowToChannelsLayer
        real, dimension(:, :), pointer              :: DiffuseFlow
        integer, intent(OUT)          , optional    :: STAT

        !Local------------------------------------------------------------------
        type (T_Node), pointer                      :: CurrNode
        integer                                     :: STAT_CALL, ready_
        integer                                     :: NodeID
        
        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)
        
        if (ready_ .EQ. IDLE_ERR_) then
            
            !Update RunOff and GW Fluxes
            do NodeID = 1, Me%TotalNodes            
                CurrNode => Me%Nodes (NodeID)
                Me%RunOffVector(NodeID) = OLFlowToChannels(CurrNode%GridI, CurrNode%GridJ)
                
                if (associated(GWFlowToChannels)) then             
                    Me%GroundVector(NodeID) = GWFlowToChannels(CurrNode%GridI, CurrNode%GridJ)
                endif
                
                if (present(GWFlowToChannelsLayer)) then 
                    Me%GWFlowByLayers     = .true.            
                    Me%GroundVectorLayers => GWFlowToChannelsLayer
                endif 
                
                if (associated(DiffuseFlow)) then
                    Me%DiffuseVector(NodeID) = DiffuseFlow(CurrNode%GridI, CurrNode%GridJ)
                endif
            enddo
            
            
            call ModifyDrainageNetLocal
            
                     
            STAT_CALL = SUCCESS_
        else               
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine ModifyDrainageNetWithGrid

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
            
    subroutine ModifyDrainageNetWithoutGrid(DrainageNetworkID, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: DrainageNetworkID  
        integer, intent(OUT)        , optional      :: STAT

        !Local------------------------------------------------------------------
        integer                                     :: STAT_CALL, ready_                
        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)
        
        if (ready_ .EQ. IDLE_ERR_) then
       
            Me%HasGrid = .false.

!            !UpdateRunOffFluxes - use this in the future????
!            do NodeID = 1, Me%TotalNodes
!                Me%RunOffVector (NodeID) = 0.0
!            end do 

            call ModifyDrainageNetLocal

            STAT_CALL = SUCCESS_
        else               
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine ModifyDrainageNetWithoutGrid

    !---------------------------------------------------------------------------

    subroutine ModifyDrainageNetLocal

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        real                                        :: SumDT, LocalDT, gA
        integer                                     :: iter, Niter, NodeID, ReachID
        integer                                     :: STAT_CALL
        logical                                     :: Restart       
        type (T_Property), pointer                  :: Property
        type(T_Reach), pointer                      :: CurrReach
        type(T_Node), pointer                       :: UpNode, DownNode, CurrNode
        real                                        :: BottomMass
        !Begin-------------------------------------------------------------------
       
        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ModifyDrainageNet")

        !Mass Balance
        Me%TotalOutputVolume       = 0.0
        Me%TotalFlowVolume        = 0.0
        Me%TotalInputVolume       = 0.0
        Me%TotalOverTopVolume     = 0.0
        Me%TotalStormWaterOutput  = 0.0
        
        if (Me%CheckMass) then
            Property => Me%FirstProperty
            do while (associated(Property))
                Property%MB%TotalStoredMass    = 0.0
                Property%MB%TotalDischargeMass = 0.0
                Property%MB%TotalOutFlowMass   = 0.0
!                Property%MB%TotalOverTopMass   = 0.0
                
                Property => Property%Next
            enddo
        endif
        
        
        !Stores Initial Value for the case that a Volume gets negative,
        !So Drainage Network can start all over
        call StoreInitialValues

        !Gets Current Time
        call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyDrainageNetLocal - ModuleDrainageNetwork - ERR01'
        
        call GetComputeTimeStep (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyDrainageNetLocal - ModuleDrainageNetwork - ERR02'

        if (Me%OutputHydro) then
            !Compute St. Venant terms - except time gradient
            do ReachID = 1, Me%TotalReaches
                CurrReach => Me%Reaches(ReachID)
                UpNode    => Me%Nodes(CurrReach%UpstreamNode)
                DownNode  => Me%Nodes(CurrReach%DownstreamNode)
                gA = Gravity * UpNode%VerticalArea
                        
                CurrReach%HydroPressure = gA * (UpNode%WaterLevel - DownNode%WaterLevel) / CurrReach%Length
                CurrReach%HydroGravity  = - gA * CurrReach%Slope
                if (CurrReach%VerticalArea .LE. AlmostZero) then
                    CurrReach%HydroFriction = 0.0
                else                        
                    CurrReach%HydroFriction = CurrReach%Manning * CurrReach%FlowNew / (CurrReach%VerticalArea        &
                                              * CurrReach%HydraulicRadius**(2./3.) )
                    CurrReach%HydroFriction = gA * CurrReach%HydroFriction**2
                endif
                                        
            enddo
        endif

        SumDT       = 0.0
        Restart     = .false.
        Niter       = Me%InternalTimeStepSplit
        !Niter       = max(Me%LastGoodNiter, 1)
        LocalDT     = Me%ExtVar%DT / Niter
        iter        = 1
        do while (iter <= Niter)

            !Transmission Losses - Should be used when running MOHID River Network only
            if (.not. Me%HasGrid .and. Me%ComputeOptions%TransmissionLosses) then
                call ModifyTransmissionLosses   (LocalDT)
            endif

            !Inputs Water from discharges
            if (Me%ComputeOptions%Discharges) then
                call ModifyWaterDischarges  (LocalDT, iter)                
            endif
            
            !Inputs Water from StormWaterModel
            if (Me%ComputeOptions%Discharges) then
                call FlowFromStormWater     (LocalDT)
            endif
            
            !Water From OverLandFlow / GW Exchange
            if (Me%HasGrid) then
                call ModifyWaterExchange    (LocalDT)
            endif

            
            call UpdateCrossSections

            !Actualizes Mapping
            call UpdateComputeFaces
            call UpdateOpenPoints
          
            !Runs Hydrodynamic
            call ModifyHydrodynamics        (LocalDT, Restart, Niter)

            !If Hydrodynamic return Restart as true, Restart with initial Solution
            if (Restart) then
                !Niter   = Niter + 1
                Niter = Niter + Me%InternalTimeStepSplit
                LocalDT = Me%ExtVar%DT / Niter
                call WriteDTLog ('ModuleDrainageNetwork', Niter, LocalDT)
                if (LocalDT < Me%ExtVar%DT / 10000.) then
                    write(*,*)'LocalDT below limit'
                    write(*,*)'Module Drainage Network Decrease Global DT'
                    stop 'ModifyDrainageNetLocal - ModuleDrainageNetwork - ERR03'
                endif
                call ResetToInitialValues ()
                call UpdateCrossSections  ()
                SumDT   = 0.0
                Restart = .false.
                iter    = 1

            else
                
                !Runs Advection / Diffusion
                if (Me%ComputeOptions%AdvectionDiffusion) call TransportProperties      (LocalDT)

                SumDT = SumDT + LocalDT
                iter  = iter  + 1

                if (Me%Output%ComputeFlowFrequency) then
                    call ComputeFlowFrequency (LocalDT, SumDT)
                endif                

                
            endif

        enddo

        Me%LastGoodNiter = Niter

        !So far, total removed volume is the flow volume
        Me%TotalOutputVolume = Me%TotalFlowVolume

        !Removes OverTop - 
        if (.not. Me%HasGrid .and. Me%ComputeOptions%RemoveOverTop) then
            call ModifyOverToped
        endif
        
        !Storm Water Inflow
        if (Me%ComputeOptions%StormWaterModelLink) then
            call FlowToStormWater
        endif

        !Top Radiation
        if (Me%ComputeOptions%TopRadiation)       call ModifyTopRadiation       ()   
            
        !Exchages heat with surface
        if (Me%ComputeOptions%SurfaceFluxes     ) call ComputeSurfaceFluxes     ()      

        !Transpiration by vegetation inside river - in drying pools
        if (Me%ComputeOptions%EVTPFromReach     ) call ComputeEVTPFromReach     ()   
        
        !Exponential decay
        if (Me%ComputeOptions%Generic_Decay     ) call GenericDecay             ()

        !Coliform decay
        if (Me%ComputeOptions%T90_Decay         ) call ColiformDecay            ()

        !Toxicity processes
        if (Me%ComputeOptions%Toxicity          ) call ModifyToxicity           ()

        !WaterQuality
        if (Me%ComputeOptions%WaterQuality      ) call ModifyWaterQuality       ()

        !Benthos
        if (Me%ComputeOptions%Benthos           ) call ModifyBenthos            ()

        !Bottom fluxes
        if (Me%ComputeOptions%BottomFluxes      ) call ComputeBottomFluxes      ()
        
        !Load Integration
        if (Me%ComputeOptions%ComputeLoad       ) call CalculateLoad            ()
        
        if (Me%HasGrid) then
            call UpdateChannelsDynamicMatrix
        endif

        if (Me%CheckMass) then
            Me%TotalStoredVolume = 0.0
            do NodeID = 1, Me%TotalNodes
                if (Me%Nodes(NodeID)%nDownStreamReaches /= 0) then
                    Me%TotalStoredVolume = Me%TotalStoredVolume + Me%Nodes(NodeID)%VolumeNew
                endif
                Property => Me%FirstProperty
                do while (associated(Property))
                    
                    CurrNode => Me%Nodes(NodeID)
                    BottomMass = 0.0
                    if (Check_Particulate_Property(Property%ID%IDNumber)) then
                        ![kg] = [kg/m2] * [m2]
                        BottomMass = Property%BottomConc(NodeID) * CurrNode%CrossSection%BottomWidth * CurrNode%Length
                    endif                
                    
                    ![kg] = [kg] + [g/m3] * [m3] * [1e-3kg/g]
                    Property%MB%TotalStoredMass = Property%MB%TotalStoredMass + BottomMass &
                                                  + Property%Concentration (NodeID)        &
                                                  * Property%ISCoefficient                 &
                                                  * Me%Nodes(NodeID)%VolumeNew
                    Property => Property%Next
                enddo                                    
!                endif
            end do
        end if

        if (Me%OutputHydro) then
            !Compute St. Venant terms - time gradient + advection
            !advection here because HydroAdvection function uses FlowOld
            do ReachID = 1, Me%TotalReaches
                CurrReach => Me%Reaches(ReachID)
                CurrReach%HydroTimeGradient = (CurrReach%FlowNew - CurrReach%FlowOld) / Me%ExtVar%DT
                CurrReach%HydroAdvection    = HydroAdvection(CurrReach, Me%ExtVar%DT)
            enddo
        endif

        if (Niter == Me%InternalTimeStepSplit) then
            Me%NextDT = Me%ExtVar%DT * Me%DTFactor
        else
            Me%NextDT = Me%ExtVar%DT / Niter
        endif
        
!        if (Niter == 1) then
!            Me%NextDT = Me%ExtVar%DT * Me%DTFactor
!        else
!            Me%NextDT = Me%ExtVar%DT / Me%LastGoodNiter
!            Me%LastGoodNiter = max(Me%LastGoodNiter -1, 1)
!        endif

!        else if (Niter <= 10) then
!            Me%NextDT = Me%ExtVar%DT * 1.0
!        else
!            Me%NextDT = Me%ExtVar%DT / 2.0
!        endif
        
        
        !call ComputeNextDT

        Property => Me%FirstProperty
        do while (associated (Property))
            if (Property%ID%IDNumber == VSS_) then
                call CalculateVSS(Property)
            end if
            if (Property%ID%IDNumber == TSS_) then
                call CalculateTSS(Property)

            end if
            Property => Property%Next    
        end do

        if (Me%TimeSerie%nNodes .GT.0) call WriteTimeSeries (LocalDT)
        
        if (Me%Output%Yes) call HDF5Output

        if (Me%WriteMaxStationValues) call MaxStationValues 

        !Restart Output
        if (Me%Output%WriteRestartFile .and. .not. (Me%CurrentTime == Me%EndTime)) then
            if(Me%CurrentTime >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                call WriteFinalFile
                Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
            endif
        endif


        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ModifyDrainageNet")


    end subroutine ModifyDrainageNetLocal

    !---------------------------------------------------------------------------

    subroutine MaxStationValues

        !Local-----------------------------------------------------------------
        integer                                 :: NodeID
        type(T_Node) , pointer                  :: CurrNode
        type(T_Reach), pointer                  :: CurrReach
        character(len=StringLength)             :: AuxString
        real                                    :: Year, Month, Day, Hour, Minute, Second
        character(len=4)                        :: CYear
        character(len=2)                        :: CMonth, CDay, CHour, CMinute, CSecond


        do NodeID = 1, Me%TotalNodes

            CurrNode => Me%Nodes(NodeID)

            

            if (CurrNode%StationName /= null_str .and. CurrNode%nDownstreamReaches /= 0) then

                CurrReach => Me%Reaches(CurrNode%DownstreamReaches(1))

                if (CurrReach%FlowNew > CurrNode%Max%Flow) then

                    CurrNode%Max%Flow  = CurrReach%FlowNew                  
                    CurrNode%Max%Vel   = CurrReach%Velocity
                    CurrNode%Max%Depth = CurrNode%WaterDepth

                    call ExtractDate(Me%CurrentTime, Year = Year, Month = Month,   &
                                 Day = Day, Hour = Hour, Minute = Minute, Second = Second)


                    write(CYear  , '(I4)') int(Year)
                    write(CMonth , '(I2)') int(Month)   ; if (Month < 10)   CMonth (1:1) = '0' 
                    write(CDay   , '(I2)') int(Day)     ; if (Day < 10)     CDay   (1:1) = '0'
                    write(CHour  , '(I2)') int(Hour)    ; if (Hour < 10)    CHour  (1:1) = '0'    
                    write(CMinute, '(I2)') int(Minute)  ; if (Minute < 10)  CMinute(1:1) = '0'
                    write(CSecond, '(I2)') int(Second)  ; if (Second < 10)  CSecond(1:1) = '0'

                    write(AuxString,*) CYear, '-', CMonth, '-', CDay, ' ', CHour, ':', CMinute, ':', CSecond

                    CurrNode%Max%Time = AuxString

                endif

                nullify(CurrReach)

            endif

        enddo

        nullify(CurrNode)

    end subroutine MaxStationValues

    !---------------------------------------------------------------------------

    subroutine ComputeFlowFrequency(LocalDT, SumDT)
        
        !Argument--------------------------------------------------------------
        real                                    :: LocalDT, SumDT
        !Local-----------------------------------------------------------------
        integer                                 :: ReachID
        type(T_Reach), pointer                  :: CurrReach
        real                                    :: TimeWindow
        !Begin-----------------------------------------------------------------
        
        if (Me%CurrentTime .gt. Me%Output%FlowFrequency%StartDate .and. Me%CurrentTime .lt. Me%Output%FlowFrequency%StopDate) then
            
            !Accumulate time with local DT because it converged and result was approved
            do ReachID = 1, Me%TotalReaches

                CurrReach => Me%Reaches(ReachID)

                if (CurrReach%FlowNew > Me%Output%FlowFrequency%MinimumFlow) then
                    
                    !s
                    CurrReach%FlowAccTime = CurrReach%FlowAccTime + LocalDT


                endif

                nullify(CurrReach)

            enddo
            
            
            !End of time step - update flow percentage
            if (Abs(Me%ExtVar%DT - SumDT) .lt. AllmostZero) then

                do ReachID = 1, Me%TotalReaches

                    CurrReach => Me%Reaches(ReachID)

                    !s
                    TimeWindow            = Me%CurrentTime - Me%Output%FlowFrequency%StartDate
                    CurrReach%FlowAccPerc = CurrReach%FlowAccTime / TimeWindow
                
                    nullify(CurrReach)
                    
                    
                enddo
                
            endif
            
            
        endif
        
    end subroutine ComputeFlowFrequency

    !---------------------------------------------------------------------------

    subroutine CalculateLoad()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer              :: Property
        integer                                 :: NodeID, i       
        type (T_Node), pointer                  :: CurrNode
        type (T_Reach), pointer                 :: DownReach
        real(8)                                 :: Flow

        nullify (Property)
        Property => Me%FirstProperty

        do while (associated (Property))

            if (Property%ComputeOptions%ComputeLoad) then 
            
                do NodeID = 1, Me%TotalNodes
                
                    !if (Me%OpenpointsFlow(NodeID) == OpenPoint) then

                        CurrNode  => Me%Nodes (NodeID)

                        Flow = 0.0

                        !Adds Inflow due to channel flow
                        do i = 1, CurrNode%nDownStreamReaches
                            DownReach => Me%Reaches (CurrNode%DownStreamReaches (i))
                            
                            if (Me%ComputeFaces(DownReach%ID) == OpenPoint) then
                                Flow = Flow +  dble(DownReach%FlowNew)
                            endif
                            
                        enddo
                
                        Property%Load(NodeID) = Property%Concentration(NodeID) * Flow
                        
                    !endif
                    
                enddo

            end if
        
            Property => Property%Next

        enddo

    end subroutine CalculateLoad

    !---------------------------------------------------------------------------

    subroutine CalculateTSS(TSSProperty)

        !Arguments--------------------------------------------------------------
        type (T_Property), pointer              :: TSSProperty

        !Local------------------------------------------------------------------
        type (T_Property), pointer              :: Property
        type (T_Property) , pointer             :: PropertySedF !Cohesive Sediment Fine
        type (T_Property) , pointer             :: PropertySedM !Cohesive Sediment Medium
        type (T_Property) , pointer             :: PropertySedC !Cohesive Sediment Coarse
        type (T_Property) , pointer             :: PropertyVSS  

        integer                                 :: NodeID      
        !-----------------------------------------------------------------------

        !Sum Porperties in TSS

        if(Me%ComputeOptions%CalcFractionSediment)then
            call SearchProperty(PropertySedF, PropertyXIDNumber = COHSED_FINE_)
            call SearchProperty(PropertySedM, PropertyXIDNumber = COHSED_MEDIUM_)
            call SearchProperty(PropertySedC, PropertyXIDNumber = COHSED_COARSE_)
            call SearchProperty(PropertyVSS,  PropertyXIDNumber = VSS_)
            do NodeID = 1, Me%TotalNodes               
               !TSS in mgTS/l
               TSSProperty%Concentration(NodeID) = PropertySedF%Concentration(NodeID)                   &
                                                 + PropertySedM%Concentration(NodeID)                   &  
                                                 + PropertySedC%Concentration(NodeID)                   &  
                                                 + PropertyVSS%Concentration(NodeID) 

               TSSProperty%BottomConc(NodeID)    = PropertySedF%BottomConc(NodeID)                   &
                                                 + PropertySedM%BottomConc(NodeID)                   &  
                                                 + PropertySedC%BottomConc(NodeID)                   &  
                                                 + PropertyVSS%BottomConc(NodeID)                    
            enddo
        else
            call SearchProperty(Property,     PropertyXIDNumber = Cohesive_Sediment_)
            call SearchProperty(PropertyVSS,  PropertyXIDNumber = VSS_)
            do NodeID = 1, Me%TotalNodes               
               !TSS in mgTS/l
               TSSProperty%Concentration(NodeID) = Property%Concentration(NodeID)                       &
                                                 + PropertyVSS%Concentration(NodeID)                   

               TSSProperty%BottomConc(NodeID)    = Property%BottomConc(NodeID)                       &
                                                 + PropertyVSS%BottomConc(NodeID)                   
            enddo
        end if
        
        !-----------------------------------------------------------------------

    end subroutine CalculateTSS

    !To sum up the concetrations of all VSSproperties in VSS property

    subroutine CalculateVSS(VSSProperty)

        !Arguments--------------------------------------------------------------
        type (T_Property), pointer              :: VSSProperty

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local------------------------------------------------------------------
        type (T_Node    ), pointer              :: CurrNode
        type (T_Property), pointer              :: Property
        integer                                 :: NodeID
        real                                    :: Ratio

        !-----------------------------------------------------------------------

        !Sum Porperties in VSS
                
                nullify (CurrNode)
                do NodeID = 1, Me%TotalNodes

            VSSProperty%Concentration(NodeID) = 0.0
            VSSProperty%BottomConc(NodeID)    = 0.0
                    
                    Property => Me%FirstProperty
                    do while (associated (Property))

                if (isVSS(Property%ID%IDNumber)) then

                    call GetWQRatio(InterfaceID = Me%ObjInterface,          &
                                    PropertyID  = Property%ID%IDNumber,     &
                                    Ratio       = Ratio,                    &
                                    STAT        = STAT_CALL)
                   
                   !VSS in mgTS/l
                   VSSProperty%Concentration(NodeID) = VSSProperty%Concentration(NodeID) + Ratio * Property%Concentration(NodeID)
                   if(Property%ComputeOptions%BottomFluxes) then
                      !VSS in kgTS/m2
                      VSSProperty%BottomConc(NodeID)    = VSSProperty%BottomConc(NodeID) + Ratio * Property%BottomConc(NodeID)
                   endif
                        end if
                        Property => Property%Next    
                    end do
                enddo

        !-----------------------------------------------------------------------

    end subroutine CalculateVSS

    subroutine ComputeNextDT

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: NodeID, i
        real                                        :: MinTime, TimeToEmpty
        real                                        :: OutFlow, InFlow
        type (T_Node), pointer                      :: CurrNode
        type (T_Reach), pointer                     :: DownReach, UpReach
        integer                                     :: STAT_CALL
        character(len=StringLength)                 :: AuxString
        logical                                     :: VariableDT


        call GetVariableDT(Me%ObjTime, VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModuleDrainageNetwork -  ERR00'
        
        if (VariableDT) then

            call GetMaxComputeTimeStep(Me%ObjTime, MinTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModuleDrainageNetwork -  ERR01'

            if (MinTime == 0.0) MinTime = Me%ExtVar%DT

            do NodeID = 1, Me%TotalNodes

                if (Me%OpenpointsFlow(NodeID) == OpenPoint) then

                    OutFlow   = 0.0
                    InFlow    = 0.0

                    CurrNode  => Me%Nodes (NodeID)

                    if (.not. CurrNode%Discharges) then

                        !Adds Outflow due to channel flow
                        do i = 1, CurrNode%nDownstreamReaches
                            DownReach => Me%Reaches (CurrNode%DownstreamReaches (i))
                            OutFlow = OutFlow + dble(DownReach%FlowNew)
                        enddo

                        !Adds Inflow due to channel flow
                        do i = 1, CurrNode%nUpstreamReaches
                            UpReach => Me%Reaches (CurrNode%UpstreamReaches (i))
                            InFlow  = InFlow + dble(UpReach%FlowNew)
                        enddo


                        !Adds Outflow due to exchange with land
                        if (Me%RunOffVector (NodeID) < 0.0) then
                            OutFlow = Outflow - Me%RunOffVector (NodeID)
                        else
                            InFlow  = InFlow  + Me%RunOffVector (NodeID)
                        endif

                        !Adds Outflow due to GW exchange
                        if (Me%GroundVector(NodeID) < 0.0) then
                            OutFlow = Outflow - Me%GroundVector (NodeID)
                        else
                            InFlow  = Inflow  + Me%GroundVector (NodeID)
                        endif

                        if (Inflow < OutFlow) then

                            
                            if (Me%HasGrid) then
                                TimeToEmpty = CurrNode%VolumeNew / (OutFlow - Inflow) * Me%DTFactor * max(Me%LastGoodNiter - 1, 1)
                            else
                                !For simple river net application it does not make sense to get a big difference between the 
                                !internal DT of Drainage Network & the main program...
                                TimeToEmpty = CurrNode%VolumeNew / (OutFlow - Inflow) * Me%DTFactor
                            endif

                            MinTime = Min(MinTime, TimeToEmpty)
                        endif

                        !Checks for flooding
                        if (CurrNode%WaterDepth > CurrNode%CrossSection%Height) then
                            write(AuxString, fmt=*)CurrNode%ID
                            call SetError(WARNING_, INTERNAL_, 'Flood at '//trim(AuxString), OFF)
                            MinTime = Min(MinTime, Me%MaxDTFlood)
                        endif
                        
                    endif

                endif

            end do

            Me%NextDT = MinTime
            
        else
        
            Me%NextDT = Me%ExtVar%DT
        
        endif

    end subroutine ComputeNextDT

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ModifyWaterDischarges (LocalDT, iter)

        !Arguments--------------------------------------------------------------
        real                                    :: LocalDT
        integer                                 :: iter

        !Local------------------------------------------------------------------
        type (T_Property), pointer              :: Property
        type (T_Node), pointer                  :: CurrNode
        integer                                 :: iDis, nDischarges, NodeID, NodePos        
        real(8)                                 :: VolumeNew    
        integer                                 :: STAT_CALL
        integer                                 :: iProp

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ModifyWaterDischarges")

        !Actualize Volumes
        do NodeID = 1, Me%TotalNodes
            Me%Nodes(NodeID)%VolumeOld = Me%Nodes(NodeID)%VolumeNew
        end do

        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleDrainageNetwork - ModifyWaterDischarges - ERR01'

        do iDis = 1, nDischarges

            NodePos = Me%DischargesLink(iDis)

            CurrNode => Me%Nodes(NodePos)

            if (iter == 1) then
                call GetDischargeWaterFlow(Me%ObjDischarges,                                &
                                        Me%CurrentTime, iDis,                               &
                                        Me%Nodes(NodePos)%WaterDepth,                       &
                                        Me%DischargesFlow(iDis), STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_) stop 'ModuleDrainageNetwork - ModifyWaterDischarges - ERR04'
            endif
            
            VolumeNew = CurrNode%VolumeNew + Me%DischargesFlow(iDis) * LocalDT

            if (Me%CheckMass) Me%TotalInputVolume = Me%TotalInputVolume + Me%DischargesFlow(iDis) * LocalDT

            nullify (Property)
            Property => Me%FirstProperty
            iProp = 0
            do while (associated (Property))
                
                if (Property%ComputeOptions%Discharges) then 
                    
                    iProp = iProp + 1

                    !Gets Discharge Concentration for this cycle of iter
                    if (iter == 1) then
                        call GetDischargeConcentration (Me%ObjDischarges,                           &
                                                        Me%CurrentTime,                             &
                                                        iDis, Me%DischargesConc(iDis, iProp),       &
                                                        Property%ID%IDNumber,                       &
                                                        STAT = STAT_CALL)
                        if (STAT_CALL/=SUCCESS_) then
                            if (STAT_CALL == NOT_FOUND_ERR_) then 
                                !When a property is not found associated to a discharge
                                !by default is consider that the concentration is zero
                                Me%DischargesConc(iDis, iProp) = 0.
                            else
                                stop 'ModuleDrainageNetwork - ModifyWaterDischarges - ERR05'
                            endif
                        endif
                    endif

                    call DischargeProperty (Me%DischargesFlow(iDis), Me%DischargesConc(iDis, iProp),        &
                                            CurrNode%VolumeNew, Property%Concentration(NodePos),            &
                                            Property%IScoefficient, LocalDT, .false.)
                                            
                    if (Me%CheckMass) then
                        !kg = kg + m3/s * s * g/m3 * 1e-3kg/g
                        Property%MB%TotalDischargeMass = Property%MB%TotalDischargeMass + (Me%DischargesFlow(iDis)           &
                                                         * LocalDT * Me%DischargesConc(iDis, iProp) * Property%IScoefficient)
                    endif

                end if
                Property => Property%Next

            enddo
        
            CurrNode%VolumeNew = VolumeNew

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ModifyWaterDischarges")

    end subroutine ModifyWaterDischarges
    
    !---------------------------------------------------------------------------

    subroutine ModifyWaterExchange (LocalDT)

        !Arguments--------------------------------------------------------------
        real                                        :: LocalDT

        !Local------------------------------------------------------------------
        integer                                     :: NodeID, K
        type (T_Node), pointer                      :: CurrNode
        type (T_Property), pointer                  :: Property
        real                                        :: GWConc

        !-----------------------------------------------------------------------

        !Actualize Volumes
        do NodeID = 1, Me%TotalNodes
            Me%Nodes(NodeID)%VolumeOld = Me%Nodes(NodeID)%VolumeNew
        end do


        !Actualizes Concentrations
        nullify (Property)
        Property => Me%FirstProperty
        do while (associated (Property))

            !Temperature
            if (Property%ID%IDNumber == Temperature_) then
                Property%OverLandConc(:) = Me%AirTemperature        (:)
                Property%GWaterConc  (:) = Me%SedimentTemperature   (:)
                Property%DWaterConc  (:) = Me%AirTemperature        (:)
            endif

            !All other properties are given in the data file

            Property => Property%Next
        enddo

        !Discharges OverLandFlow
        nullify (Property)
        Property => Me%FirstProperty
        do while (associated (Property))
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                call DischargeProperty (Me%RunOffVector (NodeID), Property%OverLandConc(NodeID),&
                                        CurrNode%VolumeNew,   Property%Concentration(NodeID),   &
                                        Property%IScoefficient, LocalDT, .false.)
            enddo
            Property => Property%Next
        enddo

        !Actualize VolumeNew
        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            CurrNode%VolumeNew = CurrNode%VolumeNew + (Me%RunOffVector (NodeID) * LocalDT)
        enddo        

        !Discharges GroundWaterFlow - Particulate property will not exit
        nullify (Property)
        Property => Me%FirstProperty
        do while (associated (Property))
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                !not compute for phantom node - it has not flow associated in porous media
                !and no limits for K defined
                if (CurrNode%nDownstreamReaches .gt. 0) then
                    if (.not. Me%GWFlowByLayers) then
                        
                        !if property particulate and flow going to river, conc matrix value is zero (not changed
                        !since the allocation because PMP particulate properties are not linked to DN)
                        call DischargeProperty (Me%GroundVector (NodeID), Property%GWaterConc(NodeID),  &
                                                CurrNode%VolumeNew,   Property%Concentration(NodeID),   &
                                                Property%IScoefficient, LocalDT,                        &
                                                Check_Particulate_Property(Property%ID%IDNumber))
                        
                    else
                        do K = Me%GWFlowBottomLayer(NodeID), Me%GWFlowTopLayer(NodeID)
                            
                            !if property particulate and flow going to river, conc matrix value should be zero 
                            !but this matrix is not allocated in DN because is 3D (only a pointer and exists for dissolved).
                            if (Check_Particulate_Property(Property%ID%IDNumber)) then
!                            if ((Check_Particulate_Property(Property%ID%IDNumber)) .and.                         &
!                                (Me%GroundVectorLayers (CurrNode%GridI, CurrNode%GridJ, k) .gt. 0.0)) then
                                GWConc = 0.0
                            else
                                GWConc = Property%GWaterConcLayers(CurrNode%GridI, CurrNode%GridJ, k)
                            endif
                            
                            call DischargeProperty (Me%GroundVectorLayers (CurrNode%GridI, CurrNode%GridJ, k),    &
                                                    GWConc,                                                       &
                                                    CurrNode%VolumeNew,   Property%Concentration(NodeID),         &
                                                    Property%IScoefficient, LocalDT,                              &
                                                    Check_Particulate_Property(Property%ID%IDNumber)) 
                                                    
                        enddo
                    endif
                endif
            enddo
            Property => Property%Next
        enddo

        !Actualize VolumeNew
        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            if (.not. Me%GWFlowByLayers) then
                CurrNode%VolumeNew = CurrNode%VolumeNew + (Me%GroundVector (NodeID) * LocalDT)
            else
                do K = Me%GWFlowBottomLayer(NodeID), Me%GWFlowTopLayer(NodeID)
                
                    CurrNode%VolumeNew = CurrNode%VolumeNew                                               &
                                         + (Me%GroundVectorLayers (CurrNode%GridI, CurrNode%GridJ, k)     &
                                         * LocalDT)             
                enddo
            endif
        enddo        
    
        !Discharges DiffuseFlow
        nullify (Property)
        Property => Me%FirstProperty
        do while (associated (Property))
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                call DischargeProperty (Me%DiffuseVector (NodeID), Property%DWaterConc(NodeID), &
                                        CurrNode%VolumeNew,   Property%Concentration(NodeID),   &
                                        Property%IScoefficient, LocalDT, .false.)
            enddo
            Property => Property%Next
        enddo

        !Actualize VolumeNew
        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            CurrNode%VolumeNew = CurrNode%VolumeNew + (Me%DiffuseVector (NodeID) * LocalDT)
        enddo        
    
               
    end subroutine ModifyWaterExchange 

    !---------------------------------------------------------------------------

    subroutine ModifyOverToped 

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: NodeID
        type (T_Node), pointer                      :: CurrNode 

        !Actualize VolumeNew
        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            if (CurrNode%WaterDepth > CurrNode%CrossSection%Height) then
                
                Me%TotalOverTopVolume   = Me%TotalOverTopVolume + (CurrNode%VolumeNew - CurrNode%VolumeMax)
                
!                Property => Me%FirstProperty
!                do while (associated(Property))
!                    !kg = kg + (m3 * g/m3 * 1e-3 kg/g)
!                    Property%MB%TotalOverTopMass = Property%MB%TotalOverTopMass + ((CurrNode%VolumeNew - CurrNode%VolumeMax)  &
!                                                   * Property%Concentration(CurrNode) * Property%ISCoefficient)
!                    Property => Property%Next
!                enddo

                CurrNode%WaterDepth     = CurrNode%CrossSection%Height
                CurrNode%VolumeNew      = CurrNode%VolumeMax

                call ComputeXSFromWaterDepth (CurrNode)

!                CurrNode%WetPerimeter   = CurrNode%CrossSection%BottomWidth         &
!                                        + 2. * CurrNode%WaterDepth                  &
!                                        * sqrt (1. + CurrNode%CrossSection%Slope**2.)

!                CurrNode%SurfaceWidth = (CurrNode%CrossSection%BottomWidth  & 
!                                      + 2. * CurrNode%CrossSection%Slope *  &
!                                        CurrNode%WaterDepth)

!                CurrNode%SurfaceArea  = CurrNode%SurfaceWidth *              &
!                                        CurrNode%Length
            
            endif
        enddo        


    end subroutine ModifyOverToped

    !---------------------------------------------------------------------------

    subroutine ModifyTransmissionLosses (LocalDT)

        !Arguments--------------------------------------------------------------
        real                                        :: LocalDT

        !Local------------------------------------------------------------------
        integer                                     :: NodeID
        type (T_Node), pointer                      :: CurrNode
        type (T_Property), pointer                  :: Property

        !Actualize Volumes
        do NodeID = 1, Me%TotalNodes
            Me%Nodes(NodeID)%VolumeOld = Me%Nodes(NodeID)%VolumeNew
        end do

        !Calculates Transmission Losses
        do NodeID = 1, Me%TotalNodes 
            CurrNode => Me%Nodes (NodeID)      

            if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
            !if (CurrNode%WaterDepth > Me%MinimumWaterDepth) then
                ![m3/s]                      =       [m/s]                    * [m]                               *    [m]
                Me%TransmissionFlow (NodeID) = -1.0* Me%HydraulicConductivity * CurrNode%CrossSection%BottomWidth *  &
                                               CurrNode%Length

                !Don't permite the volume to get negative
!                DeadVolume                   = CurrNode%CrossSection%BottomWidth * CurrNode%Length * Me%MinimumWaterDepth
!                Me%TransmissionFlow (NodeID) = min(Me%TransmissionFlow (NodeID), (CurrNode%VolumeNew - DeadVolume) / LocalDT)
            else
                Me%TransmissionFlow (NodeID) = 0.0
            endif                
        enddo

        !Discharge Properties
        nullify (Property)
        Property => Me%FirstProperty
        do while (associated (Property))
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
                !if (CurrNode%WaterDepth > Me%MinimumWaterDepth) then
                    call DischargeProperty (Me%TransmissionFlow (NodeID), Property%Concentration(NodeID),   &
                                            CurrNode%VolumeNew,   Property%Concentration(NodeID),           &
                                            Property%IScoefficient, LocalDT,                                &  
                                            Check_Particulate_Property(Property%ID%IDNumber))
                endif
            enddo
            Property => Property%Next
        enddo

        !Actualize VolumeNew
        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            CurrNode%VolumeNew = CurrNode%VolumeNew + (Me%TransmissionFlow (NodeID) * LocalDT)
        enddo        


    end subroutine ModifyTransmissionLosses

    !---------------------------------------------------------------------------

    subroutine DischargeProperty (DischargeFlow, DischargeConc, Volume,         &
                                  Concentration, ISCoef, LocalDT, Accumulate)

        !Arguments--------------------------------------------------------------
        real                                        :: DischargeFlow, DischargeConc
        real(8)                                     :: Volume
        real                                        :: Concentration, ISCoef
        real                                        :: LocalDT
        logical                                     :: Accumulate

        !Local------------------------------------------------------------------
        real(8)                                     :: DischargeVolume
        real(8)                                     :: OldMass, NewMass
        real                                        :: ISDischargeConc, ISConcentration

        ISDischargeConc = DischargeConc * ISCoef
        ISConcentration = Concentration * ISCoef

        if (abs(DischargeFlow) > AllmostZero) then
        
            if      (DischargeFlow > 0.0) then

                !Explicit discharges input 
                DischargeVolume  = dble(LocalDT)*dble(DischargeFlow)

                OldMass          = dble(ISConcentration) * Volume
                NewMass          = OldMass + DischargeVolume * dble(ISDischargeConc)                                       

                ISConcentration = NewMass / (Volume + DischargeFlow * LocalDT)

            elseif (DischargeFlow < 0.0) then
                    
                !If the discharge flow is negative (Output) then the concentration
                !to consider is the concentration of the NOD ID where the discharge
                !is located

                !If the property acculumlates in the water column 
                !(e.g particulate properties during infiltration) then the concentration will increase

                !Implicit discharges output
                DischargeVolume  = dble(LocalDT)*dble(DischargeFlow)

                OldMass          = dble(ISConcentration) * Volume

                if (Accumulate) then
                    !Concentration    = OldMass / (Volume + DischargeFlow * LocalDT)
                    NewMass          = OldMass
                else
                    NewMass          = OldMass * (1.0 + DischargeVolume / Volume)
                endif

                ISConcentration    = NewMass / (Volume + DischargeFlow * LocalDT)

            endif

        else
            
            !Do Nothing            

        endif
                    
        
        Concentration = ISConcentration / ISCoef


    end subroutine DischargeProperty

    !---------------------------------------------------------------------------

    subroutine FlowFromStormWater(LocalDT)

        !Arguments--------------------------------------------------------------
        real                                        :: LocalDT

        !Local------------------------------------------------------------------
        integer                                     :: NodeID
        type (T_Node), pointer                      :: CurrNode

        !Actualize VolumeNew
        do NodeID = 1, Me%StormWaterModelLink%nInflowNodes
            CurrNode                                => Me%Nodes (Me%StormWaterModelLink%InflowIDs(NodeID))
            Me%TotalStormWaterInput                 = Me%TotalStormWaterInput + Me%StormWaterModelLink%Inflow(NodeID)
            CurrNode%VolumeNew                      = CurrNode%VolumeNew + Me%StormWaterModelLink%Inflow(NodeID) * LocalDT
        enddo        

    end subroutine FlowFromStormWater

    !---------------------------------------------------------------------------

    subroutine FlowToStormWater

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: NodeID
        type (T_Node), pointer                      :: CurrNode

        !Actualize VolumeNew
        do NodeID = 1, Me%StormWaterModelLink%nOutflowNodes
            CurrNode                                => Me%Nodes (Me%StormWaterModelLink%OutflowIDs(NodeID))
            Me%TotalStormWaterOutput                = Me%TotalStormWaterOutput + CurrNode%VolumeNew
            Me%StormWaterModelLink%Outflow(NodeID)  = CurrNode%VolumeNew / Me%ExtVar%DT
            CurrNode%VolumeNew                      = 0.0
        enddo        


    end subroutine FlowToStormWater

    !---------------------------------------------------------------------------

    subroutine ModifyHydrodynamics (LocalDT, Restart, Niter)

        !Arguments--------------------------------------------------------------
        real                                        :: LocalDT
        logical                                     :: Restart
        integer                                     :: Niter

        !Local------------------------------------------------------------------
        integer                                     :: NodeID, ReachID
        type (T_Reach), pointer                     :: CurrReach
        !$ integer                                  :: CHUNK

        !$OMP PARALLEL PRIVATE(NodeID,ReachID)

        !$ CHUNK = 10

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ModifyHydrodynamics")

        !Actualize Volumes
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do NodeID = 1, Me%TotalNodes
            Me%Nodes(NodeID)%VolumeOld = Me%Nodes(NodeID)%VolumeNew
        end do
        !$OMP END DO

        !Actualize Reaches
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do ReachID = 1, Me%TotalReaches
            Me%Reaches(ReachID)%FlowOld = Me%Reaches(ReachID)%FlowNew
        end do
        !$OMP END DO

        !$OMP END PARALLEL

        if (Me%NumericalScheme == ExplicitScheme) then

            do ReachID = 1, Me%TotalReaches                        
                call ModifyReach (ReachID, LocalDT)
            end do

            do NodeID = 1, Me%TotalNodes
                call ModifyNode          (NodeID, LocalDT)
                call VerifyMinimumVolume (NodeID, Restart, Niter)
                if (Restart) exit
            end do
                  
            if (Me%CheckMass .and. .not. Restart) then               
                CurrReach => Me%Reaches (Me%OutletReachPos)
                Me%TotalFlowVolume   = Me%TotalFlowVolume + CurrReach%FlowNew * LocalDT
            end if

        else if (Me%NumericalScheme == ImplicitScheme) then

            call Cascade (LocalDT, Restart, Niter)

            if (Me%CheckMass .and. .not. Restart) then               
                CurrReach => Me%Reaches (Me%OutletReachPos)                    
                Me%TotalFlowVolume = Me%TotalFlowVolume + (CurrReach%FlowNew + CurrReach%FlowOld) / 2. * LocalDT
            end if
               
        end if
               
        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ModifyHydrodynamics")

    end subroutine ModifyHydrodynamics            

    !---------------------------------------------------------------------------            
    !---------------------------------------------------------------------------            

    subroutine ModifyReach (ReachID, DT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ReachID
        real                                            :: DT

        !Local-----------------------------------------------------------------          
        type(T_Reach), pointer                          :: CurrReach !, UpReach
        type (T_Node), pointer                          :: DownNode  !, UpNode
        !integer                                         :: iReach     
        !----------------------------------------------------------------------        
         

        CurrReach => Me%Reaches (ReachID)                 

if0:    if (Me%ComputeFaces (ReachID) == Compute) then
            
            DownNode => Me%Nodes (CurrReach%DownstreamNode)

if1:        if (DownNode%nDownstreamReaches .EQ. 0) then
            
                !Outlet
                select case(Me%Downstream%Boundary)
                                
                    case (Dam)
        
                        CurrReach%FlowNew   = 0.0
                        CurrReach%Velocity  = 0.0
        
                    case (ZeroDepthGradient)

                        call ComputeKinematicWave (CurrReach) 

                    case (CriticalDepth)

                        call ComputeCriticalFlow (CurrReach)

                    case(ImposedWaterLevel)

                        call ComputeStVenant (CurrReach, DT)
                                        
                    case (ImposedVelocity)

                        CurrReach%Velocity = Me%Downstream%DefaultValue
                        CurrReach%FlowNew  = CurrReach%Velocity * CurrReach%VerticalArea

                    case default

                        write(*,*) 'Invalid downstream boundary'
                        stop 'ModifyReach - ModuleDrainageNetwork - ERR02'

                end select


            else if (Me%HydrodynamicApproximation == KinematicWave) then !if1                      
                    
                call ComputeKinematicWave (CurrReach)
                
            else if (Me%HydrodynamicApproximation == DiffusionWave) then
            
                !Update Slope based on water level
                CurrReach%Slope = (Me%Nodes(CurrReach%UpstreamNode)%Waterlevel -            &
                                   Me%Nodes(CurrReach%DownstreamNode)%Waterlevel) /         &
                                   CurrReach%Length
                
                
!                !Don't allow negative slopes and impose a minimum slope..
!                if (.not. Me%AllowBackwardWater) then
!                    CurrReach%Slope = max(CurrReach%Slope, Me%MinimumSlope)
!                endif
                
                call ComputeKinematicWave (CurrReach)

            else !if1

                call ComputeStVenant (CurrReach, DT)                
        
            end if if1
                        
        else !if0

            CurrReach%FlowNew = 0.0
            CurrReach%Velocity  = 0.0

        end if if0

    end subroutine ModifyReach

    !---------------------------------------------------------------------------            

    subroutine ComputeCrossSection (CurrNode)

        !Arguments--------------------------------------------------------------
        type (T_Node), pointer                      :: CurrNode

        !Local------------------------------------------------------------------
        real                                        :: Av_New, h_New, AvTrapez2, TopH
        real(8)                                     :: PoolVolume, VolNewAux
        type(T_Reach), pointer                      :: UpReach
        type(T_Node), pointer                       :: UpNode

if1:    if (CurrNode%nDownstreamReaches /= 0) then
            
            PoolVolume = CurrNode%CrossSection%PoolDepth * CurrNode%Length * CurrNode%CrossSection%BottomWidth

            !Volume greater then volume in pools
if2:        if (CurrNode%VolumeNew > PoolVolume) then
            
                !--------------------------------------------------------------
                !COMPUTE VerticalArea AND WaterDepth --------------------------
                !--------------------------------------------------------------
    
                VolNewAux = CurrNode%VolumeNew / CurrNode%SingCoef
                Av_New    = (VolNewAux - PoolVolume) / CurrNode%Length
                
                if (CurrNode%CrossSection%Form == Trapezoidal .OR.                    &
                   (CurrNode%CrossSection%Form == TrapezoidalFlood .AND.              &
                    VolNewAux <= CurrNode%VolumeMaxTrapez1)) then

                    h_New = TrapezoidWaterHeight (b  = CurrNode%CrossSection%BottomWidth,  &
                                                  m  = CurrNode%CrossSection%Slope,        &
                                                  Av = Av_New)

                elseif (CurrNode%CrossSection%Form == TrapezoidalFlood) then
                ! from the previous if
                ! we already know that CurrNode%WaterDepth > CurrNode%CrossSection%MiddleHeigh


                    AvTrapez2 = (VolNewAux - CurrNode%VolumeMaxTrapez1) / CurrNode%Length

                    TopH =  TrapezoidWaterHeight (b  = CurrNode%CrossSection%MiddleWidth,  &
                                                  m  = CurrNode%CrossSection%SlopeTop,     &
                                                  Av = AvTrapez2)

                    h_New = CurrNode%CrossSection%MiddleHeight + TopH

                elseif (CurrNode%CrossSection%Form == Tabular) then

                    call TabularWaterLevel (CurrNode%CrossSection, Av_New, CurrNode%WaterLevel)
                    h_New = CurrNode%WaterLevel - CurrNode%CrossSection%BottomLevel

                else
                    
                    stop 'Invalid cross section form - ComputeCrossSection - ModuleDrainageNetwork - ERR01'
                end if

                !If Height exceeds channel height, just consider volume of full
                !bank width
                !if (h_New > CurrNode%CrossSection%Height) then
                !    h_New  = CurrNode%CrossSection%Height + (CurrNode%VolumeNew - PoolVolume - CurrNode%VolumeMax) &
                !             / CurrNode%Length / CurrNode%CrossSection%TopWidth
                !endif
  
                CurrNode%WaterDepth     = h_New

                !--------------------------------------------------------------
                !COMPUTE OTHER CROSS SECTION PROPERTIES -----------------------
                !--------------------------------------------------------------

                call ComputeXSFromWaterDepth (CurrNode)

                !Substarcts minumum area (Stability resons)                                                                 
                !CurrNode%VerticalArea   = max(CurrNode%VerticalArea - Me%MinimumWaterDepth * CurrNode%CrossSection%BottomWidth, 0.0)

            else !if2

                CurrNode%VerticalArea   = 0.0
                CurrNode%WaterDepth     = 0.0
                CurrNode%WetPerimeter   = 0.0

                if (CurrNode%CrossSection%PoolDepth < AllmostZero) then
                    CurrNode%SurfaceWidth   = 0.0
                    CurrNode%SurfaceArea    = 0.0
                else
                    CurrNode%SurfaceWidth   = CurrNode%CrossSection%BottomWidth
                    CurrNode%SurfaceArea    = CurrNode%SurfaceWidth * CurrNode%Length
                endif
                
            endif if2

        else !if1
        
            if (Me%Downstream%Boundary == ImposedWaterLevel) then
           
                    if (Me%Downstream%Evolution  == None .or. Me%Downstream%Evolution == OpenMI) then
                       CurrNode%WaterLevel = Me%Downstream%DefaultValue
                    else if (Me%Downstream%Evolution == ReadTimeSerie) then
                       call ModifyDownstreamTimeSerie (CurrNode%WaterLevel)
                    end if

                    CurrNode%WaterDepth = CurrNode%WaterLevel - CurrNode%CrossSection%BottomLevel
                
            else
            
                !Assumes constant slope in the last reach
                UpReach             => Me%Reaches (CurrNode%UpstreamReaches (1))
                UpNode              => Me%Nodes   (UpReach%UpstreamNode) 
                CurrNode%WaterLevel = UpNode%WaterLevel - UpReach%Slope * UpReach%Length
                CurrNode%WaterDepth = max(CurrNode%WaterLevel - CurrNode%CrossSection%BottomLevel, 0.0)
                
            endif
                  
        end if if1

        CurrNode%WaterLevel = CurrNode%WaterDepth + CurrNode%CrossSection%BottomLevel
                                                               
    end subroutine ComputeCrossSection

    !---------------------------------------------------------------------------            
    !---------------------------------------------------------------------------            

    subroutine TabularWaterLevel (CrossSection, Av, WaterLevel)

        !Arguments--------------------------------------------------------------
        type(T_CrossSection)                :: CrossSection
        real, intent(in)                    :: Av
        real, intent(out)                   :: WaterLevel

        !Locals----------------------------------------------------------------
        integer                             :: i, ilev
        real                                :: dAv, b, m, dH

        
        !if (Av > CrossSection%LevelVerticalArea(CrossSection%NLevels)) write(*,*) 'Av higher than total vertical area'

        dAv = 0.0

        do i=  1, CrossSection%NLevels
            if (CrossSection%LevelVerticalArea(i) <= Av) then
                dAv  = Av - CrossSection%LevelVerticalArea(i)
                ilev = i            
                !exit nao porque quero o lowest level mais aproximado
            endif
        enddo   

        if (dAv <= 1e-6) then
            WaterLevel = CrossSection%Level(ilev)
        else        
            
            b = CrossSection%LevelBottomWidth(ilev)
            m = 0.5 * ( abs(CrossSection%LevelSlopeLeft(ilev)) + CrossSection%LevelSlopeRight(ilev) )

            dH = TrapezoidWaterHeight (b, m, dAv)

            WaterLevel = CrossSection%Level(ilev) + dH
            
        endif

    end subroutine TabularWaterLevel


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    real function TrapezoidWaterHeight (b, m, Av)
    
        !Arguments-------------------------------------------------------------
        real                                :: b,m, Av

        !Locals----------------------------------------------------------------
        real                                :: binomio, sqrt_binomio
        real                                :: h

        if (m .LE. AllmostZero) then                   
            !Rectangular
            h  = Av  / b
        else
            binomio = b*b + 4. * m * Av
            if (binomio .LT. 0.0) then
                stop 'TrapezoidWaterHeight - ModuleDrainageNetwork - ERR01'
            else if (binomio .LE. AllmostZero) then
                sqrt_binomio = 0.
            else
                sqrt_binomio = sqrt (binomio)
            endif

            h = (- b + sqrt_binomio) / (2. * m)

        endif
        
        TrapezoidWaterHeight = h

    end function TrapezoidWaterHeight

    !---------------------------------------------------------------------------            
    !---------------------------------------------------------------------------
        
    subroutine UpdateReachCrossSection (CurrReach)
    
        !Arguments--------------------------------------------------------------
        !Local------------------------------------------------------------------
        type (T_Reach), pointer                 :: CurrReach
        type (T_Node ), pointer                 :: UpNode, DownNode    
        real                                    :: PoolDepth, PoolVolume   

        UpNode    => Me%Nodes(CurrReach%UpstreamNode)
        DownNode  => Me%Nodes(CurrReach%DownstreamNode)

        !Volume greater then volume in pools
        PoolVolume = UpNode%CrossSection%PoolDepth * UpNode%Length * UpNode%CrossSection%BottomWidth
        if (UpNode%VolumeNew > PoolVolume) then
            PoolDepth = UpNode%CrossSection%PoolDepth
        else
            PoolDepth = UpNode%VolumeNew / (UpNode%Length * UpNode%CrossSection%BottomWidth)
        endif

        if (Me%AllowBackwardWater .and. UpNode%WaterLevel < DownNode%WaterLevel .and. DownNode%ID /= Me%OutletNodePos) then
            CurrReach%VerticalArea    = (UpNode%VerticalArea + DownNode%VerticalArea) / 2.0
            CurrReach%HydraulicRadius = CurrReach%VerticalArea / ((UpNode%WetPerimeter + DownNode%WetPerimeter) / 2.0)
        else
            CurrReach%VerticalArea    = UpNode%VerticalArea
            if (UpNode%WetPerimeter > 0.0) then
                CurrReach%HydraulicRadius = UpNode%VerticalArea / UpNode%WetPerimeter
            else            
                CurrReach%HydraulicRadius = 0.0
            end if
        endif
        
        
        CurrReach%PoolVerticalArea= PoolDepth * UpNode%CrossSection%BottomWidth
        CurrReach%Manning         = UpNode%CrossSection%ManningCH


        
    end subroutine UpdateReachCrossSection

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ComputeKinematicWave (CurrReach)

        !Arguments--------------------------------------------------------------        
        type (T_Reach), pointer                 :: CurrReach
        real                                    :: sign

        if (abs(CurrReach%Slope) >  AllmostZero) then
        
            if (CurrReach%Slope > 0.0) then
                sign = 1.0
            else
                sign = -1.0
            endif


            CurrReach%FlowNew  = sign * CurrReach%VerticalArea * CurrReach%HydraulicRadius **(2.0/3.0) &
                                 * sqrt(abs(CurrReach%Slope))  / CurrReach%Manning
                                 

            CurrReach%Velocity = CurrReach%FlowNew / (CurrReach%VerticalArea + CurrReach%PoolVerticalArea)

        else
        
            CurrReach%FlowNew  = 0.0
            CurrReach%Velocity = 0.0

        endif
    
    end subroutine ComputeKinematicWave

    !---------------------------------------------------------------------------

    subroutine ComputeStVenant (CurrReach, DT) 
     

        !Arguments--------------------------------------------------------------        
        type (T_Reach), pointer                 :: CurrReach        
        real, intent (in)                       :: DT
        
        !Internal---------------------------------------------------------------
        type (T_Node ), pointer                 :: UpNode, DownNode 
        real                                    :: LevelSlope, Pressure
        real                                    :: Friction, Advection !, AdvectionUp, AdvectionDown

        UpNode   => Me%Nodes (CurrReach%UpstreamNode  )
        DownNode => Me%Nodes (CurrReach%DownstreamNode)
                            
        !PRESSURE - explicit ----------------------------------------------------
        
        !m/m              =   m                                       / m      
        LevelSlope        = (UpNode%WaterLevel - DownNode%WaterLevel) / CurrReach%Length 

        if (abs(LevelSlope) > AllmostZero) then              
            !m3/s             = s  * m/s2    * m2                     * m/m
            Pressure          = DT * Gravity * CurrReach%VerticalArea * LevelSlope
        else
            Pressure          = 0.0
        endif

        !FRICTION - semi-implicit -----------------------------------------------
!        if (UpNode%WaterDepth > Me%MinimumWaterDepth) then
            !        =  s * m/s2    * m3/s * s/m(1/3) / (m * m.... units are correct, I check on the paper.. FB
            Friction = DT * Gravity * abs(CurrReach%FlowOld) * CurrReach%Manning ** 2. &
                     / ( CurrReach%VerticalArea * CurrReach%HydraulicRadius ** (4./3.) ) 
!        else
!            Friction = 0.0
!        endif
        

        !ADVECTION - upwind (in - out)-------------------------------------------
        !positive direction is downstream
        Advection = HydroAdvection(CurrReach, DT)

        !EQUACAO ---------------------------------------------------------------

        CurrReach%FlowNew = ( CurrReach%FlowOld + Advection + Pressure )    &
                          / ( 1. + Friction )
       
!        CurrReach%FlowNew = 0.0;
!        if (abs(CurrReach%FlowNew) > 0.0) then
!            CurrReach%FlowNew = 0.0;
!        endif
       
        CurrReach%Velocity = CurrReach%FlowNew / (CurrReach%VerticalArea + CurrReach%PoolVerticalArea)

    
    end subroutine ComputeStVenant

    !---------------------------------------------------------------------------

    real function HydroAdvection (CurrReach, DT)
    
        !Arguments--------------------------------------------------------------        
        type (T_Reach), pointer                 :: CurrReach        
        real, intent (in)                       :: DT
        
        !Internal---------------------------------------------------------------
        type (T_Node ), pointer                 :: UpNode, DownNode 
        type (T_Reach), pointer                 :: UpReach, DownReach       
        real                                    :: AdvectionUp, AdvectionDown
        real                                    :: DownFlux, UpFlux
        real                                    :: DownProp, UpProp
        integer                                 :: i 

        !ADVECTION - upwind (in - out)-------------------------------------------
        !positive direction is downstream

        UpNode   => Me%Nodes (CurrReach%UpstreamNode  )
        DownNode => Me%Nodes (CurrReach%DownstreamNode)


        HydroAdvection = 0.0

       !Down Flux
        AdvectionDown = 0.0
        do i = 1, DownNode%nDownstreamReaches
            DownReach   => Me%Reaches (DownNode%DownstreamReaches (i))
            DownFlux = (CurrReach%FlowOld + DownReach%FlowOld) / 2.

            if (DownFlux.GE.0.0) then
                DownProp = CurrReach%Velocity
            else
                DownProp = DownReach%Velocity 
            end if
                                                                        
            !m4/s2        = m4/s2           m3/s     * m/s
            AdvectionDown = AdvectionDown + DownFlux * DownProp
        end do
       
        !UpFlux
        AdvectionUp = 0.0
        do i = 1, UpNode%nUpstreamReaches                    
           UpReach   => Me%Reaches (UpNode%UpstreamReaches (i))
           UpFlux = ( CurrReach%FlowOld + UpReach%FlowOld) / 2.

            if (UpFlux.GE.0.0) then
                UpProp = UpReach%Velocity
            else
                UpProp = CurrReach%Velocity
            end if

           AdvectionUp    =  AdvectionUp + UpFlux * UpProp
        end do
   
        !m3/s          =  m4/s2                        *  s / m
        HydroAdvection = (AdvectionUp - AdvectionDown) * DT / CurrReach%Length

        nullify(UpNode)
        nullify(DownNode)
                       
    end function HydroAdvection

    !---------------------------------------------------------------------------

   
    subroutine ComputeCriticalFlow (CurrReach)
    
        !Arguments--------------------------------------------------------------        
        type (T_Reach), pointer                 :: CurrReach        
    
        !Local------------------------------------------------------------------
        type (T_Node), pointer                  :: CurrNode
        real                                    :: h      !hydraulic mean depth


            nullify(CurrNode)
            CurrNode => Me%Nodes(CurrReach%UpstreamNode)

            if (CurrNode%VerticalArea .LE. AlmostZero) then
                CurrReach%FlowNew  = 0.0
                CurrReach%Velocity = 0.0
            else
                h = CurrNode%VerticalArea / CurrNode%SurfaceWidth
                CurrReach%FlowNew = CurrNode%VerticalArea * sqrt(Gravity*h)
                CurrReach%Velocity = CurrReach%FlowNew / CurrNode%VerticalArea
            endif

    end subroutine ComputeCriticalFlow

    !---------------------------------------------------------------------------
                
    subroutine ModifyNode (NodeID, DT)

        !Arguments-------------------------------------------------------------
        integer                                     :: NodeID
        real                                        :: DT
         
        !Local-----------------------------------------------------------------          
        type (T_Node ), pointer                     :: CurrNode       
        real(8)                                     :: InFlow, OutFlow
                                       
        CurrNode => Me%Nodes(NodeID)                                    

        !if (Me%OpenPointsFlow (NodeID) == OpenPoint) then
                    
            call ComputeNodeInFlow  (CurrNode, InFlow)
            call ComputeNodeOutFlow (CurrNode, OutFlow)
     
            CurrNode%VolumeNew = CurrNode%VolumeOld + ( InFlow - OutFlow ) * DT
            
        !else
        
            !CurrNode%VolumeNew = CurrNode%VolumeOld
                        
        !end if        
        
    end subroutine ModifyNode

    !---------------------------------------------------------------------------            
    subroutine VerifyMinimumVolume (NodeID, Restart, Niter)

        !Arguments--------------------------------------------------------------
        integer                                     :: NodeID
        logical                                     :: Restart
        integer                                     :: Niter

        !Local-----------------------------------------------------------------
        type (T_Node), pointer                      :: CurrNode
        
        CurrNode => Me%Nodes (NodeID)

        if (Me%Stabilize .and. Niter < Me%MaxIterations) then
            if (CurrNode%nDownstreamReaches /= 0 .and. CurrNode%VolumeOld > Me%StabilizeCoefficient * CurrNode%VolumeMax) then ! CurrNode%VolumeMin) then   
                if (abs(CurrNode%VolumeNew - CurrNode%VolumeOld) > Me%StabilizeFactor * CurrNode%VolumeMax) then
                    Restart = .true.
                    return
                endif
                if (CurrNode%VolumeNew < 0.0) then
                    Restart = .true.
                    return
                endif
            endif
        end if

    end subroutine VerifyMinimumVolume

    !---------------------------------------------------------------------------            

    subroutine Cascade (DT, Restart, Niter)

        !Arguments--------------------------------------------------------------
        real                                        :: DT
        logical                                     :: Restart
        integer                                     :: Niter

        !Local-----------------------------------------------------------------
        integer                                     :: NodeID, ReachID
        type (T_Node ), pointer                     :: CurrNode !, DownNode
        type (T_Reach), pointer                     :: CurrReach
        real(8)                                     :: InFlow, OutFlow
        real(8)                                     :: OutFlowNew, OutFlowOld, Vol
        integer                                     :: iter, MaxIter
        logical                                     :: Iterate
        real                                        :: Error, Tolerance        

        Tolerance = 0.001
        MaxIter   = 100

        
do1:    do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
        
if1:        if (Me%OpenPointsFlow (NodeID) == OpenPoint) then

                ReachID = CurrNode%DownstreamReaches (1)
                CurrReach => Me%Reaches (ReachID)
                             
                call ComputeNodeInFlow  (CurrNode, InFlow)                                    
                call ComputeNodeOutFlow (CurrNode, OutFlowOld)
               
                Iterate = .true.
                iter = 0            
do2:            do while (Iterate)

                    iter = iter + 1 
                    
                    call ComputeCrossSection     (CurrNode)
                    call UpdateReachCrossSection (CurrReach)
                    call ModifyReach             (ReachID, DT)

                    call ComputeNodeOutFlow (CurrNode, OutFlowNew)

                    OutFlow = (OutFlowOld + OutFlowNew) / 2.
                    
                    Vol = CurrNode%VolumeOld + DT * (InFlow - OutFlow)

                    Error = abs(Vol - CurrNode%VolumeNew)
                    
                    CurrNode%VolumeNew = Vol
                    
                    if (Error > Tolerance) then
                    
                        Iterate = .true.
                        CurrNode%VolumeNew = Vol    

                        if (iter >= MaxIter) then
                            write(*,*) 'Max number of iterations exceeded in Cascade'
                            stop 'ModuleDrainageNetwrok - Cascade - ERR01'
                        end if
                
                    else

                        Iterate = .false.      
                        call VerifyMinimumVolume (NodeID, Restart, Niter)
                        if (restart) exit                    

                    end if

                end do do2
                
                if (Restart) exit

            else !if1

                

                if (CurrNode%nDownstreamReaches /= 0) then
                    
                    CurrNode%VolumeNew = CurrNode%VolumeOld
                    ReachID = CurrNode%DownstreamReaches (1)
                    CurrReach => Me%Reaches (ReachID)
                    CurrReach%FlowNew  = 0.0
                    CurrReach%Velocity = 0.0

                 end if

                    
            end if if1             

        end do do1

    end subroutine Cascade

    !---------------------------------------------------------------------------            
    !---------------------------------------------------------------------------            

    subroutine UpdateCrossSections ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Node) , pointer                     :: CurrNode
        type (T_Reach), pointer                     :: CurrReach
        integer                                     :: NodeID, ReachID
        
        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "UpdateCrossSections")

        do NodeID = 1, Me%TotalNodes
            CurrNode => Me%Nodes (NodeID)
            call ComputeCrossSection (CurrNode)
        end do

        do ReachID = 1, Me%TotalReaches
            CurrReach => Me%Reaches (ReachID)
            call UpdateReachCrossSection (CurrReach)
        end do

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "UpdateCrossSections")

    end subroutine UpdateCrossSections
    
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine UpdateComputeFaces

        !Local------------------------------------------------------------------
        type (T_Reach), pointer                     :: CurrReach
        type (T_Node) , pointer                     :: UpNode, DownNode
        integer                                     :: ReachID
        integer                                     :: ComputeFaceUpDown, ComputeFaceDownUp
        real                                        :: Min_Level_Up, Min_Level_Down
        real                                        :: Level_Up, Level_Down
        
        !-----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "UpdateComputeFaces")
       
        !Based on ModuleHorizontalMap
        Me%ComputeFaces = 0

        do ReachID = 1, Me%TotalReaches

            ComputeFaceUpDown = 0
            ComputeFaceDownUp = 0
                      
            CurrReach => Me%Reaches (ReachID)

            if (CurrReach%Active) then

                UpNode    => Me%Nodes (CurrReach%UpstreamNode)
                DownNode  => Me%Nodes (CurrReach%DownstreamNode)

                Min_Level_Up   = UpNode%CrossSection%BottomLevel   + Me%MinimumWaterDepth
                Min_Level_Down = DownNode%CrossSection%BottomLevel + Me%MinimumWaterDepth

                Level_Up   = UpNode%WaterLevel
                Level_Down = DownNode%WaterLevel

                !V1----------------------
                !Alterar tambem ComputeKinematic
                if (Level_Up > Min_Level_Up .and. Level_Up > Min_Level_Down)        &
                    ComputeFaceUpDown = 1

                !V2------------------------
                !Alterar tambem ComputeKinematic
                !if (Level_Up > Min_Level_Up .and. Level_Up > Level_Down)                &
                !    ComputeFaceUpDown = 1

                !Open the face in Down Up Direction if Hydrodynamic aprox. allows backwater
                if (Me%AllowBackwardWater) then
                    if (Level_Down > Min_Level_Down .and. Level_Down > Min_Level_Up)    &
                        ComputeFaceDownUp = 1
                endif

                if (ComputeFaceUpDown + ComputeFaceDownUp > 0) then
                    Me%ComputeFaces (ReachID) = 1
                endif
                
            else
            
                !In active Reach
                Me%ComputeFaces (ReachID) = 0 
            
            endif

        end do

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "UpdateComputeFaces")

    end subroutine UpdateComputeFaces

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine UpdateOpenPoints

        !Local-----------------------------------------------------------------        
        type (T_Node) , pointer                     :: CurrNode
        integer                                     :: NodeID, i
        integer                                     :: UpReachID, DownReachID
        integer                                     :: Sum
        real                                        :: hInPool

        Me%OpenPointsFlow    = 0
        Me%OpenPointsProcess = 0

        if (Me%HasGrid) then
            Me%ChannelsOpenProcess = 0
        endif

        do NodeID = 1, Me%TotalNodes

            Sum = 0.0
            CurrNode => Me%Nodes (NodeID)

            if (CurrNode%nDownstreamReaches /= 0) then

                do i = 1, CurrNode%nUpstreamReaches
                    UpReachID = CurrNode%UpstreamReaches (i)
                    Sum = Sum + Me%ComputeFaces (UpReachID)
                end do

                do i = 1, CurrNode%nDownstreamReaches
                    DownReachID = CurrNode%DownstreamReaches (i)
                    Sum = Sum + Me%ComputeFaces (DownReachID)
                end do

                if (Sum > 0) Me%OpenPointsFlow (NodeID) = OpenPoint

                if (CurrNode%CrossSection%PoolDepth > 0.0) then
                    hInPool = CurrNode%VolumeNew / (CurrNode%CrossSection%BottomWidth * CurrNode%Length)
                    if (hInPool > Me%MinimumWaterDepthProcess) then
                        Me%OpenPointsProcess(NodeID) = OpenPoint
                        if (Me%HasGrid) then
                            Me%ChannelsOpenProcess(CurrNode%GridI, CurrNode%GridJ) = OpenPoint
                        endif
                    endif
                else
                    if (CurrNode%WaterDepth > Me%MinimumWaterDepthProcess) then
                        Me%OpenPointsProcess(NodeID) = OpenPoint
                        if (Me%HasGrid) then
                            Me%ChannelsOpenProcess(CurrNode%GridI, CurrNode%GridJ) = OpenPoint
                        endif
                    endif
                endif
            
            end if

        end do


    end subroutine UpdateOpenPoints

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine UpdateChannelsDynamicMatrix

        !Local-----------------------------------------------------------------
        integer                                     :: NodeID
        type (T_Node), pointer                      :: CurrNode

        !----------------------------------------------------------------------
       
        do NodeID = 1, Me%TotalNodes
            
            CurrNode => Me%Nodes (NodeID)                
            
            Me%ChannelsWaterLevel       (CurrNode%GridI, CurrNode%GridJ) = CurrNode%WaterLevel
                    
        enddo

    end subroutine UpdateChannelsDynamicMatrix 

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    
    subroutine ModifyDownstreamTimeSerie (NewValue)

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        type (T_Time)                                   :: Time1, Time2
        real                                            :: Value1, Value2, NewValue
        logical                                         :: TimeCycle
        
        !Begin----------------------------------------------------------------
        
        !Gets Value for current Time
        call GetTimeSerieValue (Me%Downstream%ObjTimeSerie, Me%CurrentTime,             &
                                Me%Downstream%DataColumn,                               &
                                Time1, Value1, Time2, Value2, TimeCycle,                &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ModifyDownstreamTimeSerie - ModuleDrainageNetwork - ERR02'

        
            if (TimeCycle) then
                NewValue = Value1
            else
                !Interpolates Value for current instant
                call InterpolateValueInTime(Me%CurrentTime,                             &
                                            Time1, Value1,                              &
                                            Time2, Value2, NewValue)
            endif

    end subroutine ModifyDownstreamTimeSerie

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ComputeNodeInFlow (CurrNode, InFlow)

        !Arguments---------------------------------------------------------------
        type (T_Node), pointer                      :: CurrNode 
        real(8), intent(OUT)                        :: InFlow
        
        !Local-------------------------------------------------------------------
        integer                                     :: i
        type (T_Reach), pointer                     :: UpReach
        !------------------------------------------------------------------------
        
        InFlow = 0.0        
        do i = 1, CurrNode%nUpstreamReaches

            nullify (UpReach)
            UpReach => Me%Reaches (CurrNode%UpstreamReaches (i))            

            if (Me%ComputeFaces(UpReach%ID) == OpenPoint) then

                if (Me%NumericalScheme == ExplicitScheme) then
                    InFlow = InFlow + UpReach%FlowNew
                else 
                    InFlow = InFlow + (UpReach%FlowNew + UpReach%FlowOld) / 2. 
                end if
            
            endif

            !InFlow = InFlow + (1. - Me%NumericalScheme) * UpReach%FlowOld &
            !       + Me%NumericalScheme * UpReach%FlowNew            

        end do

    end subroutine ComputeNodeInFlow

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ComputeNodeOutFlow (CurrNode, OutFlow)

        !Arguments---------------------------------------------------------------
        type (T_Node), pointer                      :: CurrNode 
        !integer                                     :: NodeID
        real(8), intent(OUT)                        :: OutFlow
        
        !Local-------------------------------------------------------------------
        integer                                     :: i
        type (T_Reach), pointer                     :: DownReach       
        !------------------------------------------------------------------------
        
        OutFlow = 0.0       
     
        do i = 1, CurrNode%nDownstreamReaches

            nullify (DownReach)
            DownReach => Me%Reaches (CurrNode%DownstreamReaches (i))
            
            if (Me%ComputeFaces(DownReach%ID) == OpenPoint) then
                OutFlow = OutFlow + dble(DownReach%FlowNew)
            endif
        end do

    end subroutine ComputeNodeOutFlow

    !---------------------------------------------------------------------------
    
    subroutine TransportProperties (LocalDT)

        !Arguments-------------------------------------------------------------
        real                                        :: LocalDT

        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: Property
!        type (T_Reach), pointer                     :: CurrReach
!        integer                                     :: CurrNode
        !Begin------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "TransportProperties")

       
        !OLD NEW Stuff
        nullify (Property)
        Property => Me%FirstProperty                                                    
        do while (associated (Property))
            Property%ConcentrationOld = Property%Concentration
            Property => Property%Next
        enddo
        
        !Transports Properties
        call Advection_Diffusion      (LocalDT)
        
        !Set MinimumConcentration of Properties - This will create Mass
        call SetMinimumConcentration 
        
        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "TransportProperties")
               
                
    end subroutine TransportProperties 
    
    !---------------------------------------------------------------------------
    
    subroutine Advection_Diffusion (LocalDT) 

        !Arguments--------------------------------------------------------------
        real                                    :: LocalDT

        !Local------------------------------------------------------------------
        type (T_Property), pointer              :: Property
        type (T_Node    ), pointer              :: CurrNode
        real                                    :: Advection, Diffusion
        real                                    :: AdvOutFlow, DifOutFlow
        integer                                 :: NodeID

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "Advection_Diffusion")


        nullify (Property)
        Property => Me%FirstProperty                                                    

        do while (associated (Property))

            if (Me%CheckMass) then
                !Mass Outflow
                AdvOutFlow = 0.0
                DifOutFlow = 0.0
            endif
            
            !FB - 12/08/2011
            !1. I check a changed some signs, which during backwater efects, where wrong.
            !2. This routine should be optimized. Allocate one diffusion matrix / advection matrix in the beginning of
            !   the simulation. Then calculate first the fluxes (ones for each reach)
            !3. Then update the concenctration in all points based on the vector
            !4. Mass checking routines should be placed in the end into one single if block
            !
            if (Property%ComputeOptions%AdvectionDiffusion) then
                
                do NodeID = 1, Me%TotalNodes

                    if (Me%OpenPointsFlow (NodeID) == OpenPoint) then
                       
                        CurrNode => Me%Nodes (NodeID)
                
                        !Sai - Entra
                        call ComputeAdvection  (Advection, Property, NodeID, AdvOutFlow)  !gX/s
                        call ComputeDiffusion  (Diffusion, Property, NodeID, DifOutFlow)  !gX/s               

                        
                        !New Concentration    
                        Property%Concentration (NodeID) = (Property%ConcentrationOld (NodeID) * CurrNode%VolumeOld + &
                                                           LocalDT * (Advection + Diffusion )) / CurrNode%VolumeNew
                                                           
                   end if
                end do
                
            endif
            
            
            if (Me%CheckMass) then
                !kg = g/s * s * 1-3kg/g
                Property%MB%TotalOutFlowMass = Property%MB%TotalOutFlowMass + (AdvOutFlow + DifOutFlow) * LocalDT * 1e-3
            endif
            
            Property => Property%Next
        
        end do

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "Advection_Diffusion")


    end subroutine Advection_Diffusion

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ComputeAdvection  (AdvectionFlux, Property, NodePos, AdvOutFlow)  

        !Arguments--------------------------------------------------------------
        real                                    :: AdvectionFlux  !kg/s     
        type (T_Property), pointer              :: Property
        integer                                 :: NodePos
        real                                    :: AdvOutflow
        !Local------------------------------------------------------------------
        type (T_Node    ), pointer              :: CurrNode              
        type (T_Reach   ), pointer              :: DownReach, UpReach, OutletReach
        integer                                 :: i
        real                                    :: DownProp, UpProp
        real                                    :: DownFlux, UpFlux
        !Begin------------------------------------------------------------------

!        !AdvectionFlux = (Conc.Q)_Sai - (Conc.Q)_Entra 
        !AdvectionFlux = (Conc.Q)_Down - (Conc.Q)_Up !Revision 6/4/2010 David
        
        CurrNode => Me%Nodes(NodePos)


        !Down Flux - Positive if flow is downstream. Reduces mass in nodes
        DownFlux = 0.0
        do i = 1, CurrNode%nDownstreamReaches
            DownReach   => Me%Reaches (CurrNode%DownstreamReaches (i))

            if (Me%ComputeFaces(DownReach%ID) == OpenPoint) then

                if (DownReach%FlowNew.GE.0.0) then
                    DownProp = Property%ConcentrationOld (NodePos    )
                else
                    DownProp = Property%ConcentrationOld (DownReach%DownstreamNode)
                end if
                                                                                
                DownFlux    =  DownFlux + DownReach%FlowNew * DownProp

            endif

        end do
        
        if (Me%CheckMass) then
            OutletReach => Me%Reaches (Me%OutletReachPos)               
            if (NodePos == OutletReach%UpstreamNode) then
                !Flow exiting
                if (DownFlux .gt. 0.0) then
                    !g/s
                    AdvOutflow = DownFlux
                endif
            endif
        endif
           
        !UpFlux
        UpFlux = 0.0
        do i = 1, CurrNode%nUpstreamReaches                    
           UpReach   => Me%Reaches (CurrNode%UpstreamReaches (i))
           
            if (Me%ComputeFaces(UpReach%ID) == OpenPoint) then
           
                if (UpReach%FlowNew.GE.0.0) then
                    UpProp = Property%ConcentrationOld (UpReach%UpstreamNode)
                else
                    UpProp = Property%ConcentrationOld (NodePos  )
                end if

                UpFlux    =  UpFlux + UpReach%FlowNew * UpProp
            
            endif
            
        end do
         
        !Advection flux is the entering flux to the node, less the exiting flux
        AdvectionFlux = UpFlux - DownFlux
        

    end subroutine ComputeAdvection

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    
    subroutine  ComputeDiffusion  (DiffusionFlux, Property, NodePos, DifOutFlow)

        !Arguments--------------------------------------------------------------
        real                                    :: DiffusionFlux  ![]/s     
        type (T_Property), pointer              :: Property
        integer                                 :: NodePos
        real                                    :: DifOutFlow
        !Local------------------------------------------------------------------
        type (T_Node    ), pointer              :: CurrNode, DownNode, UpNode              
        type (T_Reach   ), pointer              :: DownReach, UpReach, OutletReach
        integer                                 :: DownNodePos, UpNodePos, i
        real                                    :: DownFlux, UpFlux
        real                                    :: GradProp


       
        CurrNode => Me%Nodes(NodePos)

if1:    if (Property%Diffusion_Scheme == CentralDif) then

            DownFlux = 0.0
            do i = 1, CurrNode%nDownstreamReaches
                DownReach   => Me%Reaches (CurrNode%DownstreamReaches (i))
                DownNodePos =  DownReach%DownstreamNode
                DownNode    => Me%Nodes (DownNodePos)

                if (Me%OpenPointsFlow(DownNodePos) == OpenPoint) then

                    GradProp    = (Property%ConcentrationOld (DownNodePos) - Property%ConcentrationOld (NodePos)) &
                                / DownReach%Length

                    DownFlux    =  DownFlux + Property%Diffusivity * GradProp * DownReach%VerticalArea

                endif

            end do
 
             if (Me%CheckMass) then
                OutletReach => Me%Reaches (Me%OutletReachPos)               
                if (NodePos == OutletReach%UpstreamNode) then
                    !Diffusion exiting
                    if (DownFlux .lt. 0.0) then
                        !g/s
                        DifOutflow = - DownFlux
                    endif
                endif
            endif
 
            UpFlux = 0.0
            do i = 1, CurrNode%nUpstreamReaches
                UpReach   => Me%Reaches (CurrNode%UpstreamReaches (i))
                UpNodePos =  UpReach%UpstreamNode
                UpNode    => Me%Nodes (UpNodePos)

                if (Me%OpenPointsFlow(UpNodePos) == OpenPoint) then

                    GradProp  = (Property%ConcentrationOld (NodePos) - Property%ConcentrationOld (UpNodePos)) &
                                 / UpReach%Length

                    UpFlux    =  UpFlux + Property%Diffusivity * GradProp * UpReach%VerticalArea

                endif
            end do

        else if (Property%Diffusion_Scheme == UpwindOrder1) then

            write (*,*) 'Upwind discretization for diffusion of properties'
            write (*,*) 'Not yet implemented'
            stop 'ComputeDiffusion - ModuleDrainageNetwork - ERR01'
                
        end if if1        
        
        DiffusionFlux = UpFlux - DownFlux

    end subroutine  ComputeDiffusion
    
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    
    subroutine SetMinimumConcentration 

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer              :: Property
        type (T_Node), pointer                  :: CurrNode
        integer                                 :: NodeID

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "SetMinimumConcentration")


        nullify (Property)
        Property => Me%FirstProperty                                                    

        do while (associated (Property))

            if (Property%ComputeOptions%MinConcentration) then

                do NodeID = 1, Me%TotalNodes

                    CurrNode => Me%Nodes (NodeID)

                    if (Property%Concentration (NodeID) .LT. Property%MinValue) then
                        
                        if (CurrNode%VolumeNew .GT. AllMostZero) then

                            Property%MassCreated (NodeID) = Property%MassCreated (NodeID)       &
                                                          + (Property%MinValue                  &
                                                          -  Property%Concentration (NodeID)    &
                                                          *  Property%ISCoefficient)            &
                                                          * CurrNode%VolumeNew 
                
                            Property%Concentration (NodeID) = Property%MinValue
                    
                        endif

                    endif

                enddo
                
            endif

            Property => Property%Next
        
        end do

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "SetMinimumConcentration")

     
    end subroutine SetMinimumConcentration
    
    !---------------------------------------------------------------------------

    subroutine ModifyTopRadiation

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: i
        real                                        :: SWPercentage, LWPercentage
        logical                                     :: NeedsParameters      = .false.
        logical                                     :: NeedsConcentrations  = .false.
        type(T_Property), pointer                   :: PropertyX

        !----------------------------------------------------------------------

        call GetRadiationPercentages(Me%ObjLightExtinction, SWPercentage, LWPercentage, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyTopRadiation - ModuleDrainageNetwork - ERR01'

        do i = 1, Me%TotalReaches
            Me%ShortWaveField(i) = SWPercentage * Me%TopRadiation(i)
            Me%LongWaveField (i) = LWPercentage * Me%TopRadiation(i)
        enddo

        !Updates Light Extinction Coefs
        call GetLightExtinctionOptions(LightExtinctionID    = Me%ObjLightExtinction,        & 
                                       NeedsParameters      = NeedsParameters,              &
                                       NeedsConcentrations  = NeedsConcentrations,          &
                                       STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Compute_SWExtCoefField - ModuleDrainageNetwork - ERR02'
        
        if (NeedsConcentrations) then

            PropertyX => Me%FirstProperty
            do while (associated(PropertyX))

                if (PropertyX%ComputeOptions%LightExtinction) then

                    if(NeedsParameters)then

                        call ModifyLightExtinctionField(LightExtinctionID   = Me%ObjLightExtinction,          &
                                                        RiverPoints1D       = Me%RiverPoints,                 &
                                                        CurrentTime         = Me%CurrentTime,                 &
                                                        PropertyID          = PropertyX%ID%IDNumber,          &
                                                        Concentration       = PropertyX%Concentration,        &
                                                        UnitsCoef           = PropertyX%IScoefficient,        &
                                                        ExtinctionParameter = PropertyX%ExtinctionCoefficient,&
                                                        STAT                = STAT_CALL)
                        if (STAT_CALL/= SUCCESS_) stop 'Compute_SWExtCoefField - ModuleDrainageNetwork - ERR03'

                    else

                        call ModifyLightExtinctionField(LightExtinctionID   = Me%ObjLightExtinction,          &
                                                        RiverPoints1D       = Me%RiverPoints,                 &
                                                        CurrentTime         = Me%CurrentTime,                 &
                                                        PropertyID          = PropertyX%ID%IDNumber,          &
                                                        Concentration       = PropertyX%Concentration,        &
                                                        UnitsCoef           = PropertyX%IScoefficient,        &
                                                        STAT                = STAT_CALL)
                        if (STAT_CALL/= SUCCESS_) stop 'Compute_SWExtCoefField - ModuleDrainageNetwork - ERR04'

                    end if


                endif

                PropertyX=>PropertyX%Next

            enddo

            nullify(PropertyX)

        else

            call ModifyLightExtinctionField(LightExtinctionID   = Me%ObjLightExtinction,          &
                                            RiverPoints1D       = Me%RiverPoints,                 &
                                            CurrentTime         = Me%CurrentTime,                 &
                                            STAT                = STAT_CALL)
            if (STAT_CALL/= SUCCESS_) stop 'Compute_SWExtCoefField - ModuleDrainageNetwork - ERR05'

        end if

    end subroutine ModifyTopRadiation

    !---------------------------------------------------------------------------

    subroutine ComputeSurfaceFluxes ()
            
        !Arguments--------------------------------------------------------------
        
        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: Temperature, Oxygen
        type(T_Property), pointer                   :: Property
        type (T_Node    ) , pointer                 :: CurrNode
        integer                                     :: NodeID, STAT
        real                                        :: TotalHeatFlux
        real                                        :: InfraRed_, LatentHeat_
        real                                        :: SensibleHeat_, GroundHeatFlux_
        real                                        :: DOSAT, Palt, BottomSolarFlux
        real,    parameter                          :: LatentHeatOfVaporization = 2.5e6         ![J/kg]
        real,    parameter                          :: ReferenceDensity         = 1000.         ![kg/m3]
        real                                        :: Evaporation, Ka, Kl
        real                                        :: Flow, Vel, Slope

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ComputeSurfaceFluxes")

        !Surface Heat Fluxes             
        call SearchProperty (Temperature, Temperature_, .true., STAT)

        if (STAT == SUCCESS_ .and. Temperature%ComputeOptions%SurfaceFluxes) then

            do NodeID = 1, Me%TotalNodes

                if (Me%OpenPointsProcess(NodeID) == OpenPoint) then
            
           
                    !Downwelling Longwave Radiation (Wunderlich et. al. 1968)
                    InfraRed_       = LongWaveDownward (Me%CloudCover    (NodeID),              &
                                                        Me%AirTemperature(NodeID)) +            &
                                      LongWaveUpward   (Temperature%Concentration(NodeID))

                    !Latent Heat
                    LatentHeat_     = LatentHeat (ReferenceDensity,                             &
                                                  Temperature%Concentration(NodeID),            &
                                                  Me%AirTemperature  (NodeID),                  &
                                                  Me%RelativeHumidity(NodeID),                  &
                                                  Me%WindSpeed       (NodeID))

                    !SensibleHeat
                    SensibleHeat_   = SensibleHeat(ReferenceDensity,                            &
                                                   Temperature%Concentration(NodeID),           &
                                                   Me%AirTemperature  (NodeID),                 &
                                                   Me%WindSpeed       (NodeID))


                    !Exchange With Ground - CE-QUAL-W2 manual, pag. 231. Temperture for sediment = 1 day
                    Me%SedimentTemperature(NodeID) = (Temperature%Concentration(NodeID)     *   &
                                                      Me%ExtVar%DT                          +   &
                                                      Me%SedimentTemperature   (NodeID)     *   &
                                                      (86400.0 - Me%ExtVar%DT))             /   &
                                                      86400.0

                    GroundHeatFlux_ = 2.2 * (Temperature%Concentration(NodeID) - Me%SedimentTemperature(NodeID))
                
                    !Net Solar (in - out)
                    BottomSolarFlux = Me%ShadingFactor * Me%TopRadiation (NodeID) * exp(-1./20. *  &
                                      Me%Nodes(NodeID)%WaterDepth)

                    !Sum all
                    TotalHeatFlux   = Me%ShadingFactor * Me%TopRadiation (NodeID) +                  &
                                      LatentHeat_ + SensibleHeat_ +  InfraRed_ - GroundHeatFlux_ -   &
                                      BottomSolarFlux

                    ![Celsius]      =  [Celsius]  + [Joules / m^2 / s] * [s] / [kg/m^3] / [Joules/kg/Celsius] / [m]
                    Temperature%Concentration(NodeID) = Temperature%Concentration(NodeID)   +   &
                                                        TotalHeatFlux * Me%ExtVar%DT        *   &
                                                        Me%Nodes(NodeID)%SurfaceArea        /   &
                                                        ReferenceDensity                    /   &
                                                        SpecificHeatDefault                 /   & 
                                                        Me%Nodes(NodeID)%VolumeNew


                    !Calculates Evaporation
                    ![m3/s]                    = [J/m2/s] / [J/kg] / [kg/m3] * [m] * [m]
                    Evaporation     = LatentHeat_                                           /   &
                                      LatentHeatOfVaporization                              /   &
                                      ReferenceDensity                                      *   &
                                      Me%Nodes(NodeID)%SurfaceArea
                    
                    !Just considers loss of water
                    if (Evaporation < 0.0) then
                        
                        !Discharge Properties
                        nullify (Property)
                        Property => Me%FirstProperty
                        do while (associated (Property))
                            CurrNode => Me%Nodes (NodeID)                            
                            call DischargeProperty (Evaporation, Property%Concentration(NodeID),            &
                                                    CurrNode%VolumeNew,   Property%Concentration(NodeID),   &
                                                    Property%IScoefficient, Me%ExtVar%DT, ON)
                            Evaporation = Evaporation
                            Property => Property%Next
                        enddo
                            
                        
                        
                        !Update volume
                        Me%Nodes(NodeID)%VolumeNew = Me%Nodes(NodeID)%VolumeNew + Evaporation * Me%ExtVar%DT

                        if (Me%CheckMass) then
                            Me%TotalOutputVolume = Me%TotalOutputVolume - Evaporation * Me%ExtVar%DT
                        endif
                    endif

                endif

            end do

            !Updates Cross Sections    
            call UpdateCrossSections

            !Actualizes Mapping
            call UpdateComputeFaces
            call UpdateOpenPoints


        endif

        !Oxygen Surface Flux
        call SearchProperty (Oxygen, Oxygen_, .false., STAT)

        if (STAT == SUCCESS_) then
            if (Oxygen%ComputeOptions%SurfaceFluxes) then

                call SearchProperty (Temperature, Temperature_, .true., STAT)
                if (STAT /= SUCCESS_) stop 'ComputeSurfaceFluxes - ModuleDrainageNetwork - ERR10' 

                !
                ! Eq. taken from CE-QUAL-W2 user manual (Version 3.1)
                ! p. 275
                ! Referencia Melching and Flores (1999)
                !  
                !

                do NodeID = 1, Me%TotalNodes

                    if (Me%OpenPointsProcess(NodeID) == OpenPoint) then
                    
                        !Points to Node
                        CurrNode => Me%Nodes (NodeID)

                        !Mortimer Altitude Correction
                        Palt  = (1.0 - CurrNode%CrossSection%BottomLevel / 1000.0 / 44.3) ** 5.25
                        DOSAT = Palt * exp(7.7117 - 1.31403*(log(Temperature%Concentration(NodeID) + 45.93)))
                

                        if (CurrNode%nDownstreamReaches .NE. 0) then

                            Flow  = Me%Reaches(CurrNode%DownstreamReaches(1))%FlowNew
                            Vel   = Me%Reaches(CurrNode%DownstreamReaches(1))%Velocity
                            Slope = Me%Reaches(CurrNode%DownstreamReaches(1))%Slope

                            select case (Me%AerationEquation)

                            case (PoolAndRifle_)

                                !Ka [1/day]
                                if (Flow > AllmostZero) then
                                    if (Flow <= 0.556) then
                                        Ka = 517. * (Vel * Slope) ** 0.524 * Flow ** (-0.242)
                                    else
                                        Ka = 596. * (Vel * Slope) ** 0.528 * Flow ** (-0.136)
                                    endif
                                else
                                    Ka = 0.0  !So it will be set to the minimum value
                                endif

                            case (ChannelControled_)

                                !Ka [1/day]
                                if (Flow > AllmostZero) then
                                    if (Flow <= 0.556) then
                                        Ka = 88.0 * (Vel * Slope) ** 0.313 * CurrNode%WaterDepth ** (-0.353)
                                    else
                                        Ka = 142. * (Vel * Slope) ** 0.333 * CurrNode%WaterDepth ** (-0.660) * &
                                             CurrNode%SurfaceArea ** (-0.243)
                                    endif
                                else
                                    Ka = 0.0  !So it will be set to the minimum value
                                endif
                        
                            end select                            
                        
                            !Kl [m/day]
                            Kl = Ka * CurrNode%WaterDepth

                            !minimum value of KL in m/day
                            if (KL <= 0.6) KL = 0.6

                            !temperature correction
                            KL = KL * 1.024**(Temperature%Concentration(NodeID) - 20.0)

                            !convert from m/d to m/s
                            KL = KL / 86400.             


                            !New Concentration
                            Oxygen%Concentration(NodeID) = Oxygen%Concentration(NodeID)                     +   &
                                                           KL                                               *   &
                                                           CurrNode%SurfaceArea                             *   &
                                                           (DOSAT - Oxygen%Concentration(NodeID))           *   &
                                                           Me%ExtVar%DT                                     /   &
                                                           CurrNode%VolumeNew        


                        endif

                    endif

                enddo


            endif

        endif

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ComputeSurfaceFluxes")


    end subroutine ComputeSurfaceFluxes

    !---------------------------------------------------------------------------


    subroutine ComputeEVTPFromReach ()
        
        !if this process connected it may be removing water from evapotranspiration and from
        !evaporation by latent heat in routine surface fluxes - check this
            
        !Arguments--------------------------------------------------------------
        
        !Local------------------------------------------------------------------
        type(T_Property), pointer                   :: Property
        type (T_Node    ) , pointer                 :: CurrNode
        integer                                     :: NodeID
        real                                        :: EVTPFlux

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ComputeEVTPFromReach")


        do NodeID = 1, Me%TotalNodes
            
            CurrNode => Me%Nodes (NodeID)
            
            !Evaporate in pools with water depth lower than maximum
            if ((Me%OpenPointsProcess(NodeID) == OpenPoint)                                       &
                 .and. (CurrNode%CrossSection%PoolDepth .gt. 0.0)                                 &
                 .and. (CurrNode%WaterDepth .gt. 0.0)                                             &
                 .and. (CurrNode%WaterDepth .le. Me%EVTPMaximumDepth)) then
        
                !m3/s = m/s * m2 * -
                EVTPFlux = CurrNode%EVTP * CurrNode%Length * CurrNode%CrossSection%BottomWidth * Me%EVTPCropCoefficient
                    
                !Discharge Properties
                nullify (Property)
                Property => Me%FirstProperty
                do while (associated (Property))
                    CurrNode => Me%Nodes (NodeID)                            
                    call DischargeProperty (EVTPFlux, Property%Concentration(NodeID),             &
                                            CurrNode%VolumeNew,   Property%Concentration(NodeID), &
                                            Property%IScoefficient, Me%ExtVar%DT, ON)
                    Property => Property%Next
                enddo
                        
                    
                !Update volume
                Me%Nodes(NodeID)%VolumeNew = Me%Nodes(NodeID)%VolumeNew + EVTPFlux * Me%ExtVar%DT

                if (Me%CheckMass) then
                    Me%TotalOutputVolume = Me%TotalOutputVolume - EVTPFlux * Me%ExtVar%DT
                endif

            endif

        end do

        !Updates Cross Sections    
        call UpdateCrossSections

        !Actualizes Mapping
        call UpdateComputeFaces
        call UpdateOpenPoints


        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ComputeEVTPFromReach")


    end subroutine ComputeEVTPFromReach

    !---------------------------------------------------------------------------

    subroutine GenericDecay ()
        
        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------       
        type (T_Property), pointer                  :: Property
        integer                                     :: NodeID
        type (T_node), pointer                      :: CurrNode
        real                                        :: OldMass, MassSink, NewMass
        !Begin------------------------------------------------------------------
        
        
        nullify (Property)
        Property => Me%FirstProperty
        
        do while (associated (Property))
        
            if (Property%ComputeOptions%Generic_Decay) then
        
                do NodeID = 1, Me%TotalNodes

                    if (Me%OpenPointsProcess (NodeID) == OpenPoint) then                    
                        
                        CurrNode => Me%Nodes (NodeID)
                        
                        !Decay occurs as WQ process without volume change
                        !g = g/m3 * m3
                        OldMass = Property%Concentration (NodeID) * CurrNode%VolumeNew
                        
                        !P = P0*exp(-kt)  
                        MassSink = min (OldMass - OldMass * exp(-Property%DecayRate * Me%ExtVar%DT),  OldMass)
                        
                        NewMass = OldMass - MassSink
                        
                        Property%Concentration (NodeID) =  NewMass / CurrNode%VolumeNew
                                                           
                    endif
                                                                      
                end do
                
            endif
            
            Property => Property%Next
            
        enddo

    end subroutine GenericDecay

    !---------------------------------------------------------------------------

    subroutine ColiformDecay ()
        
        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------       
        type (T_Property), pointer                  :: Sal, Temp, Coliforms
        integer                                     :: NodeID, STAT
        real                                        :: T90

        !Get salinity, temperature and FecalColiforms       
        call SearchProperty (Sal, Salinity_     , .true., STAT)
        call SearchProperty (Temp, Temperature_ , .true., STAT)
        call SearchProperty (Coliforms, Fecal_Coliforms_, .true., STAT)

        do NodeID = 1, Me%TotalNodes

            if (Me%OpenPointsProcess (NodeID) == OpenPoint) then                    

                T90 = ComputeT90( Sal%Concentration(NodeID), Temp%Concentration(NodeID),    &  
                                  Me%ShadingFactor * Me%TopRadiation (NodeID),              &
                                  Me%T90Var_Method )
                                        
                Coliforms%Concentration (NodeID) =  Coliforms%Concentration (NodeID) /      &
                                                    (1.0 + Me%ExtVar%DT * (log(10.) / T90))
            endif
                                                              
        end do


    end subroutine ColiformDecay

    !---------------------------------------------------------------------------

    real function ComputeT90 (Sal, Temp, TopRadiation, T90Var_Method)

        !Arguments--------------------------------------------------------------
        real                                        :: Sal, Temp, TopRadiation
        integer                                     :: T90Var_Method  
          
        !Local------------------------------------------------------------------
        real                                        :: Light

        !Begin------------------------------------------------------------------

        select case (T90Var_Method)

        case (Constant)

            ComputeT90 = Me%T90

        case (Canteras)

            ComputeT90 = ComputeT90_Canteras (Temp, Sal, TopRadiation)

        case (Chapra  )

            !Converts W in ly/hr
            Light      = 0.086325 * TopRadiation
            ComputeT90 = ComputeT90_Chapra (Temp, Sal, Light)  

        case default

            write (*,*) 'T90 calculation method unknown'
            stop 'ComputeT90 - ModuleDrainageNetwork - ERR1'

        end select

    end function ComputeT90

    !---------------------------------------------------------------------------

    subroutine ModifyToxicity

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: NodeID
        real                                        :: maxFraction, minEC50

        Me%GlobalToxicity = 0.0
        

        call ComputeToxicityForEachEffluent

        !ComputeGlobalToxicity

do1:    do NodeID = 1, Me%TotalNodes

            maxFraction = null_real
            minEC50     = - null_real

            nullify (Property)
            Property => Me%FirstProperty                                                    

do2:        do while (associated (Property))

if1:            if (Property%ComputeOptions%Toxicity) then
            

if2:                if (Me%GlobalToxicityEvolution == 'MAX') then

                        if (Property%Toxicity%Field (NodeID) .GT. Me%GlobalToxicity (NodeID)) then
                            Me%GlobalToxicity (NodeID) = Property%Toxicity%Field (NodeID)
                        end if

                    else if (Me%GlobalToxicityEvolution == 'SUM') then !if2

                        Me%GlobalToxicity (NodeID) = Me%GlobalToxicity (NodeID)             &
                                                   + Property%Toxicity%Field (NodeID)

                    else if (Me%GlobalToxicityEvolution == 'RISKRATIO') then !if2
                               
                        if (maxFraction .LT. Property%Concentration (NodeID)) then
                            maxFraction = Property%Concentration (NodeID)
                        end if
                    
                        if (minEC50 .GT. Property%Toxicity%EC50) then
                            minEC50 = Property%Toxicity%EC50
                        end if

                    end if if2
                                                     
                end if if1

                Property => Property%Next

            end do do2

            if (Me%GlobalToxicityEvolution == 'RISKRATIO') then
                Me%GlobalToxicity(NodeID) = maxFraction / minEC50
            end if
            
        end do do1

    end subroutine ModifyToxicity
    
    !---------------------------------------------------------------------------
    
    subroutine ComputeToxicityForEachEffluent
    
        !Local--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: NodeID

        nullify (Property)
        Property => Me%FirstProperty                                                    

        do while (associated (Property))

if1:        if (Property%ComputeOptions%Toxicity) then
        
if2:            if (Property%Toxicity%Evolution == Saturation) then

                    do NodeID = 1, Me%TotalNodes
            
                        Property%Toxicity%Field (NodeID) = Property%Concentration (NodeID)      &
                                                         / (Property%Toxicity%EC50              &
                                                         +  Property%Concentration (NodeID))
                    end do

                else if (Property%Toxicity%Evolution == RiskRatio) then !if2

                    do NodeID = 1, Me%TotalNodes
            
                        Property%Toxicity%Field (NodeID) = Property%Concentration (NodeID)      &
                                                         / Property%Toxicity%EC50                                                 
                    end do

                else !if2

                    do NodeID = 1, Me%TotalNodes
            
                        Property%Toxicity%Field (NodeID) = Property%Toxicity%Slope &
                                                         * Property%Concentration (NodeID)

                    end do

                end if if2

            end if if1

            Property => Property%Next
        end do


    end subroutine ComputeToxicityForEachEffluent

    !---------------------------------------------------------------------------

    subroutine ModifyWaterQuality

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        real, dimension(:), pointer                 :: ShortWaveExtinctionField
        integer                                     :: STAT_CALL, NodePos
        real(8)                                     :: PoolVolume, PoolDepth

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ModifyWaterQuality")

        call GetShortWaveExtinctionField(Me%ObjLightExtinction, ShortWaveExtinctionField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyWaterQuality - ModuleDrainageNetwork - ERR01'

        !Updates Me%NodesDWZ
        do NodePos = 1, Me%TotalNodes

            PoolVolume = Me%Nodes(NodePos)%CrossSection%PoolDepth * Me%Nodes(NodePos)%Length * &
                         Me%Nodes(NodePos)%CrossSection%BottomWidth
            if (Me%Nodes(NodePos)%VolumeNew > PoolVolume) then
                PoolDepth = Me%Nodes(NodePos)%CrossSection%PoolDepth
            else
                PoolDepth = Me%Nodes(NodePos)%VolumeNew / (Me%Nodes(NodePos)%Length * Me%Nodes(NodePos)%CrossSection%BottomWidth)
            endif

            Me%NodesDWZ(NodePos) = Me%Nodes(NodePos)%WaterDepth + PoolDepth
        enddo

        if (Me%CurrentTime .GE. Me%Coupled%WQM%NextCompute) then

            PropertyX => Me%FirstProperty
            do while(associated(PropertyX))

                call Modify_Interface(InterfaceID       = Me%ObjInterface,              &
                                      PropertyID        = PropertyX%ID%IDNumber,        &
                                      Concentration     = PropertyX%Concentration,      &
                                      RiverPoints1D     = Me%RiverPoints,               &
                                      OpenPoints1D      = Me%OpenPointsProcess,         &
                                      ShortWaveTop      = Me%ShortWaveField,            &
                                      LightExtCoefField = ShortWaveExtinctionField,     &
                                      DWZ               = Me%NodesDWZ,                  &
                                      STAT              = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'ModifyWaterQuality - ModuleDrainageNetwork - ERR02'

                PropertyX => PropertyX%Next

            end do

            Me%Coupled%WQM%NextCompute = Me%Coupled%WQM%NextCompute + Me%Coupled%WQM%DT_Compute
            
        end if 

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%ComputeOptions%WaterQuality) then
    
!                if (Me%CurrentTime .GE.PropertyX%Evolution%NextCompute) then

               call Modify_Interface(InterfaceID   = Me%ObjInterface,               &
                                     PropertyID    = PropertyX%ID%IDNumber,         &
                                     Concentration = PropertyX%Concentration,       &
                                     DTProp        = Me%ExtVar%DT,                  &
                                     RiverPoints1D = Me%RiverPoints,                &
                                     OpenPoints1D  = Me%OpenPointsProcess,          &
                                     STAT          = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ModifyWaterQuality - ModuleDrainageNetwork - ERR03'
!                endif

            endif
            
            PropertyX=>PropertyX%Next

        enddo


        call UnGetLightExtinction(Me%ObjLightExtinction, ShortWaveExtinctionField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterQuality_Processes - ModuleWaterProperties - ERR06'
        
        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ModifyWaterQuality")


    end subroutine ModifyWaterQuality

    !---------------------------------------------------------------------------

    subroutine ModifyBenthos

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: NodeID, STAT_CALL
        type (T_Node), pointer                      :: CurrNode        
        real(8),dimension(:),pointer                :: WaterVolume
                

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ModifyBenthos")

        if (Me%CurrentTime .GE. Me%Coupled%Benthos%NextCompute) then
                allocate (WaterVolume (1:Me%TotalNodes))
                do NodeID = 1, Me%TotalNodes
                        CurrNode => Me%Nodes (NodeID)
                        WaterVolume(NodeID) = CurrNode%VolumeNew
                enddo

            PropertyX => Me%FirstProperty
            do while(associated(PropertyX))

                if(PropertyX%ComputeOptions%Benthos)then

                    if(PropertyX%ComputeOptions%BottomFluxes)then
                        do NodeID = 1, Me%TotalNodes
                            CurrNode => Me%Nodes (NodeID)
                            if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
                                !kg = kg/m2 * m * m
                                PropertyX%MassInKg(NodeID) = PropertyX%BottomConc(NodeID)      * &
                                                             CurrNode%CrossSection%BottomWidth * &
                                                             CurrNode%Length
                                                             
                            end if
                        enddo
                    else
                        do NodeID = 1, Me%TotalNodes
                            CurrNode => Me%Nodes (NodeID)
                            if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
                                !kg = g/m3 * m3 * 1e-3 kg/g
                                PropertyX%MassInKg(NodeID) = PropertyX%Concentration(NodeID) * &
                                                             CurrNode%VolumeNew              * &
                                                             PropertyX%IScoefficient        
                            end if

                        enddo
                    end if

                endif

                select case(PropertyX%ID%IDNumber)
                
                    case(Temperature_)
                        
                        call Modify_Interface(InterfaceID       = Me%ObjBenthicInterface,       &
                                              PropertyID        = PropertyX%ID%IDNumber,        &
                                              Concentration     = PropertyX%Concentration,      &
                                              RiverPoints1D     = Me%RiverPoints,               &
                                              OpenPoints1D      = Me%OpenPointsProcess,         &
                                              STAT              = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                            &
                            stop 'ModifyBenthos - ModuleDrainageNetwork - ERR01'

                    case(Oxygen_)
                        
                        call Modify_Interface(InterfaceID       = Me%ObjBenthicInterface,       &
                                              PropertyID        = PropertyX%ID%IDNumber,        &
                                              Concentration     = PropertyX%MassInKg,           &
                                              RiverPoints1D     = Me%RiverPoints,               &
                                              OpenPoints1D      = Me%OpenPointsProcess,         &
                                              Oxygen1D          = PropertyX%Concentration,      &
                                              STAT              = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                            &
                            stop 'ModifyBenthos - ModuleDrainageNetwork - ERR01'

                    case default

                        call Modify_Interface(InterfaceID       = Me%ObjBenthicInterface,       &
                                              PropertyID        = PropertyX%ID%IDNumber,        &
                                              Concentration     = PropertyX%MassInKg,           &
                                              RiverPoints1D     = Me%RiverPoints,               &
                                              OpenPoints1D      = Me%OpenPointsProcess,         &
                                              WaterVolume       = WaterVolume,                  &
                                              STAT              = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                            &
                            stop 'ModifyBenthos - ModuleDrainageNetwork - ERR01'
                end select 
                  
                PropertyX => PropertyX%Next

            end do

            Me%Coupled%Benthos%NextCompute = Me%Coupled%Benthos%NextCompute + Me%Coupled%Benthos%DT_Compute
            
        end if

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%ComputeOptions%Benthos) then

                if(PropertyX%ComputeOptions%BottomFluxes)then
                    do NodeID = 1, Me%TotalNodes
                        CurrNode => Me%Nodes (NodeID)
                        if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
                            !kg = kg/m2 * m * m
                            PropertyX%MassInKg(NodeID) = PropertyX%BottomConc(NodeID)      * &
                                                         CurrNode%CrossSection%BottomWidth * &
                                                         CurrNode%Length
                                                         
                        end if
                    enddo
                else
                    do NodeID = 1, Me%TotalNodes
                        CurrNode => Me%Nodes (NodeID)
                        if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
                            !kg = g/m3 * m3 * 1e-3 kg/g
                            PropertyX%MassInKg(NodeID) = PropertyX%Concentration(NodeID) * &
                                                         CurrNode%VolumeNew              * &
                                                         PropertyX%IScoefficient        
                        end if

                    enddo
                end if

                call Modify_Interface(InterfaceID  = Me%ObjBenthicInterface,        &
                                     PropertyID    = PropertyX%ID%IDNumber,         &
                                     Concentration = PropertyX%MassInKg,            &
                                     DTProp        = Me%ExtVar%DT,                  &
                                     RiverPoints1D = Me%RiverPoints,                &
                                     OpenPoints1D  = Me%OpenPointsProcess,          &
                                     WaterVolume   = WaterVolume,                   &
                                     STAT          = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ModifyBenthos - ModuleDrainageNetwork - ERR02'
                    nullify  (WaterVolume)

                if(PropertyX%ComputeOptions%BottomFluxes)then
                    do NodeID = 1, Me%TotalNodes
                        CurrNode => Me%Nodes (NodeID)
                        if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
                            !kg/m2 = kg / m2
                            PropertyX%BottomConc(NodeID) = PropertyX%MassInKg(NodeID)         / &
                                                           (CurrNode%CrossSection%BottomWidth * &
                                                            CurrNode%Length)  
                        end if
                    enddo
                else
                    do NodeID = 1, Me%TotalNodes
                        CurrNode => Me%Nodes (NodeID)
                        if (Me%OpenPointsProcess (NodeID) == OpenPoint) then
                            !g/m3 = kg / m3 / 1e-3 kg/g
                            PropertyX%Concentration(NodeID) = PropertyX%MassInKg(NodeID) / &
                                                              CurrNode%VolumeNew         / &
                                                              PropertyX%IScoefficient
                        end if

                    enddo
                end if
            endif

            PropertyX=>PropertyX%Next

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ModifyBenthos")

    end subroutine ModifyBenthos

    !---------------------------------------------------------------------------

    subroutine ComputeBottomFluxes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleDrainageNetwork", "ComputeBottomFluxes")

        call ModifyShearStress

        call ComputeErosionFluxes

        call ComputeDepositionFluxes

        if (MonitorPerformance) call StopWatch ("ModuleDrainageNetwork", "ComputeBottomFluxes")

    
    end subroutine ComputeBottomFluxes
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ModifyShearStress

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ReachID
        type (T_Reach), pointer                     :: CurrReach
        real                                        :: Chezy
 
        do ReachID = 1, Me%TotalReaches

            if (Me%ComputeFaces(ReachID) == Compute) then
                
                CurrReach => Me%Reaches (ReachID)

                if (CurrReach%HydraulicRadius > AllmostZero) then
                    Chezy = Gravity * CurrReach%Manning**2.0   &
                          / CurrReach%HydraulicRadius** (1./3.)
                else 
                    Chezy = 0.0
                end if

                Me%ShearStress (ReachID) = SigmaDensityReference * Chezy * CurrReach%Velocity**2.0

            else
            
                Me%ShearStress (ReachID) = 0.0
                            
            end if

        end do

    end subroutine ModifyShearStress

    !---------------------------------------------------------------------------
    
    subroutine ComputeErosionFluxes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: NodeID, ReachID
        type (T_Node    ) , pointer                 :: CurrNode
        type (T_Property) , pointer                 :: Property
        real, dimension(:), pointer                 :: SedimentConc
        real                                        :: ErosionRate
        real                                        :: ErodedConc, ErodedMass        
        real                                        :: BottomArea
        real                                        :: aux
        integer                                     :: STAT_CALL


        call SearchProperty(Property, PropertyXIDNumber = TSS_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
         call SearchProperty(Property, PropertyXIDNumber = Cohesive_Sediment_)
        endif
        
        SedimentConc => Property%BottomConc

        nullify (Property)
        Property => Me%FirstProperty                                                    
        do while (associated (Property))

if1:        if (Property%ComputeOptions%BottomFluxes                                       &
                .AND. (Property%ID%IDNumber /= VSS_ .AND. Property%ID%IDNumber /= TSS_)) then

                Property%ErosionRate = 0.0

                do NodeID = 1, Me%TotalNodes

                    ErodedConc  = 0.0        
                    ErosionRate = 0.0


if2:                if (Me%OpenPointsProcess (NodeID) == OpenPoint) then

                        CurrNode => Me%Nodes (NodeID)
                        ReachID  = CurrNode%DownstreamReaches (1)
                    
if3:                    if (Me%ShearStress (ReachID) > Property%ErosionCriticalShear .and.  &
                            CurrNode%VolumeNew > AllmostZero .and.                          &
                            SedimentConc (NodeID) > AllmostZero) then

                            aux = Me%ShearStress (ReachID) / Property%ErosionCriticalShear - 1.0


                            !kg m-2 s-1 = kg m-2 s-1 * m2
                            ErosionRate = Property%ErosionCoefficient                       &
                                        * Property%BottomConc (NodeID)                      &
                                        / SedimentConc (NodeID)                             &
                                        * aux 
                        
                            ErodedConc = ErosionRate * Me%ExtVar%DT
                                                  
                            if (ErodedConc < Property%BottomConc (NodeID) ) then

                                Property%BottomConc (NodeID) = Property%BottomConc (NodeID) &
                                                             - ErodedConc

                            else

                                ErodedConc = Property%BottomConc (NodeID) - Property%BottomMinConc 
                                Property%BottomConc (NodeID) = Property%BottomMinConc

                                ErosionRate = ErodedConc / Me%ExtVar%DT
                            
                            end if

                        

                            BottomArea = CurrNode%CrossSection%BottomWidth * CurrNode%Length
                            ErodedMass = ErodedConc * BottomArea
                           
                            Property%Concentration (NodeID) = Property%Concentration (NodeID) &
                                                            + ErodedMass / CurrNode%VolumeNew &
                                                            / Property%IScoefficient


                            Property%ErosionRate (NodeID) = ErosionRate
                                                              
                        end if if3
                    end if if2

                    if (Property%ComputeOptions%BottomFluxes) then
                        CurrNode => Me%Nodes (NodeID)
                        Property%TotalConc (NodeID) = Property%ConcentrationOld (NodeID) * 10E-3 *CurrNode%WaterDepth &
                                                    + Property%BottomConc (NodeID)                          
                        nullify (CurrNode)
                    end if
                end do
            
            end if if1

            Property => Property%Next
        enddo

    end subroutine ComputeErosionFluxes


    !---------------------------------------------------------------------------
    
    subroutine ComputeDepositionFluxes
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: NodeID, ReachID
        type (T_Node    ) , pointer                 :: CurrNode
        type (T_Property) , pointer                 :: Property
        real                                        :: DepositedMass, Mass, DepositionRate        
        real                                        :: BottomArea, MinimumMass, aux, SPMConc



        nullify (Property)
        Property => Me%FirstProperty                                                    
        do while (associated (Property))

if1:        if (Property%ComputeOptions%BottomFluxes) then

                Property%DepositionRate = 0.0

                do NodeID = 1, Me%TotalNodes

                    DepositedMass  = 0.0
                    DepositionRate = 0.0

if2:                if (Me%OpenPointsProcess (NodeID) == OpenPoint) then                                

                        CurrNode => Me%Nodes (NodeID)
                        ReachID  = CurrNode%DownstreamReaches (1)
                    
if3:                    if (Me%ShearStress (ReachID) < Property%DepositionCriticalShear .and.               &
                            CurrNode%VolumeNew > AllmostZero) then
                        
                            !ModifySettlingVelocity - [m s-1]
                            if (Property%Ws_Type == SPMFunction) then

                                SPMConc = Property%Concentration (NodeID) * Property%ISCoefficient
                                Property%Ws (NodeID) = SettlingVelocity (SPMConc,                           &
                                                                         Property%CHS, Property%KL,         &
                                                                         Property%KL1, Property%M,          &
                                                                         Property%ML, NodeID)
    
                            else

                                Property%Ws (NodeID) = Property%Ws_Value
                            end if
                    

                            
                            BottomArea = CurrNode%CrossSection%BottomWidth * CurrNode%Length

                            aux = 1.0 - Me%ShearStress (ReachID) / Property%DepositionCriticalShear 

                            !kg m-2 s-1 = kg m-3 * m s-1
                            DepositionRate = Property%Concentration (NodeID) * Property%IScoefficient       &
                                           * Property%Ws (NodeID) * aux
                            
                            DepositedMass = DepositionRate * BottomArea * Me%ExtVar%DT
                            Mass = Property%Concentration (NodeID) * Property%IScoefficient * CurrNode%VolumeNew 

                            if (DepositedMass < Mass ) then

                                Property%Concentration (NodeID) = Property%Concentration (NodeID)           &
                                                                - DepositedMass / CurrNode%VolumeNew        &
                                                                / Property%IScoefficient

                            else

                                MinimumMass = Property%MinValue * Property%IScoefficient * CurrNode%VolumeNew

                                DepositedMass = Mass - MinimumMass
                                Property%Concentration (NodeID) = Property%MinValue
                                
                                DepositionRate = DepositedMass / (BottomArea * Me%ExtVar%DT)

                            
                            end if
        
                        
                            Property%BottomConc (NodeID) = Property%BottomConc (NodeID)                     &
                                                         + DepositedMass / BottomArea

                            Property%DepositionRate (NodeID) = DepositionRate

                        end if if3
                    end if if2

                    if (Property%ComputeOptions%BottomFluxes) then
                        CurrNode => Me%Nodes (NodeID)
                        Property%TotalConc (NodeID) = Property%Concentration (NodeID) * 10E-3 *CurrNode%WaterDepth &
                                                    + Property%BottomConc (NodeID)                          
                        nullify (CurrNode)
                    end if
                    
                end do            

            end if if1

            Property => Property%Next
        end do

    end subroutine ComputeDepositionFluxes

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------   

    !Copied from ModuleFreeVerticalMovement
    !Settling velocity computed as function of a hindered settling concentration
    !Used only for cohesive sediments. Units in this formulation in kg/m3
    real function SettlingVelocity (Concentration, CHS, KL, KL1, M, ML, NodeID)
        
        !Arguments-------------------------------------------------
        real,       intent(IN) :: Concentration  !kg/m3
        real,       intent(IN) :: CHS              !Hindered settling concentration
        real,       intent(IN) :: KL, KL1, M, ML
        integer,    intent(IN) :: NodeID

        !Local-----------------------------------------------------
        real                    :: Aux
        
        !Begin-----------------------------------------------------

        !kg/m3
        
        Aux = KL1 * (Concentration - CHS)

        if (Concentration < CHS .and. Concentration >= 0.) then 

            SettlingVelocity = KL*(Concentration)**M

        elseif(Aux < 1. .and. Aux >= 0.) then

            SettlingVelocity = KL*(CHS)**M*(1.0-Aux)**ML 

        elseif(Aux > 1. .and. Concentration < 100000.) then

            SettlingVelocity = 0. !if concentration is to high settling velocity is null

        else

            write(*,*)'Concentration (g/l)          = ', Concentration
            write(*,*)'KL1 * (Concentration - CHS)  = ', Aux
            write(*,*)'NodeID                       = ', NodeID
            stop 'Error computing the settling velocity - SettlingVelocity - ModuleDrainageNetwork'

        endif

    end function SettlingVelocity

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine StoreInitialValues

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: NodeID, ReachID
        type (T_Property), pointer                  :: Property

        !Initial Volumes
        do NodeID = 1, Me%TotalNodes
            Me%Nodes(NodeID)%InitialVolumeNew       = Me%Nodes(NodeID)%VolumeNew
            Me%Nodes(NodeID)%InitialVolumeOld       = Me%Nodes(NodeID)%VolumeOld
            Me%Nodes(NodeID)%InitialWaterDepth      = Me%Nodes(NodeID)%WaterDepth
        end do

        !Initial Flows
        do ReachID = 1, Me%TotalReaches
             Me%Reaches(ReachID)%InitialFlowNew = Me%Reaches(ReachID)%FlowNew
             Me%Reaches(ReachID)%InitialFlowOld = Me%Reaches(ReachID)%FlowOld
        end do
        
        !Initial Concentration
        nullify (Property)
        Property => Me%FirstProperty                                                    
        do while (associated (Property))
            Property%InitialConcentration    = Property%Concentration
            Property%InitialConcentrationOld = Property%ConcentrationOld
            Property => Property%Next
        enddo

        if (Me%CheckMass) then
            Me%InitialTotalOutputVolume = Me%TotalOutputVolume
            Me%InitialTotalFlowVolume   = Me%TotalFlowVolume
            Me%InitialTotalInputVolume  = Me%TotalInputVolume
        end if

        if (Me%Output%ComputeFlowFrequency) then
            do ReachID = 1, Me%TotalReaches
                 Me%Reaches(ReachID)%InitialFlowAccTime = Me%Reaches(ReachID)%FlowAccTime
            end do        
        endif

    end subroutine StoreInitialValues

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine ResetToInitialValues

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                     :: NodeID, ReachID
        type (T_Property), pointer                  :: Property

        !Initial Volumes
        do NodeID = 1, Me%TotalNodes
            Me%Nodes(NodeID)%VolumeNew      = Me%Nodes(NodeID)%InitialVolumeNew
            Me%Nodes(NodeID)%VolumeOld      = Me%Nodes(NodeID)%InitialVolumeOld
            Me%Nodes(NodeID)%WaterDepth     = Me%Nodes(NodeID)%InitialWaterDepth
        end do

        !Initial Flows
        do ReachID = 1, Me%TotalReaches
            Me%Reaches(ReachID)%FlowNew = Me%Reaches(ReachID)%InitialFlowNew
            Me%Reaches(ReachID)%FlowOld = Me%Reaches(ReachID)%InitialFlowOld
        end do
      
        !Initial Concentration
        nullify (Property)
        Property => Me%FirstProperty                                                    
        do while (associated (Property))
            Property%Concentration    = Property%InitialConcentration
            Property%ConcentrationOld = Property%InitialConcentrationOld
            Property => Property%Next
        enddo

        if (Me%CheckMass) then
            Me%TotalOutputVolume = Me%InitialTotalOutputVolume
            Me%TotalFlowVolume   = Me%InitialTotalFlowVolume
            Me%TotalInputVolume  = Me%InitialTotalInputVolume
        end if

        if (Me%Output%ComputeFlowFrequency) then
            do ReachID = 1, Me%TotalReaches
                 Me%Reaches(ReachID)%FlowAccTime = Me%Reaches(ReachID)%InitialFlowAccTime
            end do        
        endif

    end subroutine ResetToInitialValues

    !---------------------------------------------------------------------------    

    subroutine WriteTimeSeries(LocalDT)
    
        !Arguments--------------------------------------------------------------
        real                                    :: LocalDT
    
        if (Me%TimeSerie%ByNodes) then
            call WriteTimeSeriesByNodes(LocalDT)
        else
            call WriteTimeSeriesByProp(LocalDT)
        end if

    end subroutine WriteTimeSeries
    
    !---------------------------------------------------------------------------    
    !---------------------------------------------------------------------------    

    subroutine WriteTimeSeriesByNodes (LocalDT)

        !Arguments--------------------------------------------------------------
        real                                    :: LocalDT

        !Local------------------------------------------------------------------
        type (T_Node    ), pointer                          :: CurrNode
        type (T_Reach   ), pointer                          :: CurrReach 
        type (T_Property), pointer                          :: Property
        integer                                             :: NodePos, i, j        
        real                                                :: PercentageMaxVolume
        real                                                :: PoolDepth, PoolVolume
        integer                                             :: STAT_CALL
        

        i = 0
do1:    do NodePos = 1, Me%TotalNodes    
            
if1:        if (Me%Nodes(NodePos)%TimeSerie) then
                
                i = i + 1
                CurrNode   => Me%Nodes (NodePos)    
                CurrReach  => Me%Reaches (CurrNode%DownstreamReaches(1))            

                PercentageMaxVolume = CurrNode%VolumeNew / CurrNode%VolumeMax * 100.0

                PoolVolume          = CurrNode%CrossSection%PoolDepth * CurrNode%Length * CurrNode%CrossSection%BottomWidth

                !Volume greater then volume in pools
                if (CurrNode%VolumeNew > PoolVolume) then
                    PoolDepth = CurrNode%CrossSection%PoolDepth
                else
                    PoolDepth = CurrNode%VolumeNew / (CurrNode%Length * CurrNode%CrossSection%BottomWidth)
                endif
            
                
                Me%TimeSerie%DataLine (pWaterDepth         ) = CurrNode%WaterDepth
                Me%TimeSerie%DataLine (pWaterLevel         ) = CurrNode%WaterLevel
                Me%TimeSerie%DataLine (pVolume             ) = CurrNode%VolumeNew
                Me%TimeSerie%DataLine (pPercentageMaxVolume) = PercentageMaxVolume
                Me%TimeSerie%DataLine (pFlow               ) = CurrReach%FlowNew
                Me%TimeSerie%DataLine (pVelocity           ) = CurrReach%Velocity
                Me%TimeSerie%DataLine (pVerticalArea       ) = CurrNode%VerticalArea
                Me%TimeSerie%DataLine (pFlowToChannels     ) = Me%RunOffVector (NodePos)
                Me%TimeSerie%DataLine (pGWFlowToChannels   ) = Me%GroundVector (NodePos)
                Me%TimeSerie%DataLine (pPoolDepth          ) = PoolDepth
                Me%TimeSerie%DataLine (pDT                 ) = Me%ExtVar%DT
                Me%TimeSerie%DataLine (pDTLocal            ) = LocalDT

                j = BaseTimeSeries

                if (Me%OutputHydro) then
                    Me%TimeSerie%DataLine (pHydroTimeGradient   ) = CurrReach%HydroTimeGradient
                    Me%TimeSerie%DataLine (pHydroAdvection      ) = CurrReach%HydroAdvection
                    Me%TimeSerie%DataLine (pHydroPressure       ) = CurrReach%HydroPressure
                    Me%TimeSerie%DataLine (pHydroGravity        ) = CurrReach%HydroGravity
                    Me%TimeSerie%DataLine (pHydroFriction       ) = CurrReach%HydroFriction

                    j = j + 5

                endif

if2:            if (Me%HasProperties) then

                    Property => Me%FirstProperty
                    j = j + 1
                                        
                    do while (associated (Property))                

                        if (Property%ComputeOptions%TimeSerie) then
                            Me%TimeSerie%DataLine (j) = Property%Concentration      (NodePos)                            
                            j = j + 1

                            if (Property%ComputeOptions%BottomFluxes) then
                                Me%TimeSerie%DataLine (j) = Property%BottomConc     (NodePos)                            
                                j = j + 1
                            end if

                            if (Property%ComputeOptions%Toxicity) then
                                Me%TimeSerie%DataLine (j) = Property%Toxicity%Field (NodePos)                            
                                j = j + 1
                            end if


                        end if
                
                        Property => Property%Next

                    end do

                   if (Me%ComputeOptions%BottomFluxes) then                       
                       Me%TimeSerie%DataLine (j) = Me%ShearStress (CurrNode%DownstreamReaches(1))                    
                       j = j + 1
                   end if

                   if (Me%ComputeOptions%Toxicity) then
                       Me%TimeSerie%DataLine (j) = Me%GlobalToxicity (NodePos)                    
                       j = j + 1
                   end if

                end if if2
                                
                call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),           &
                                        DataLine  = Me%TimeSerie%DataLine,      &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByNodes - ERR01'
            
            end if if1        
        
        end do do1
    
    end subroutine WriteTimeSeriesByNodes
    
    !---------------------------------------------------------------------------    
    !---------------------------------------------------------------------------    

    subroutine WriteTimeSeriesByProp (LocalDT)

        !Arguments--------------------------------------------------------------
        real                                    :: LocalDT

        !Local------------------------------------------------------------------
        integer                                 :: STAT_CALL
        type (T_Node    ), pointer              :: CurrNode
        type (T_Reach   ), pointer              :: CurrReach 
        type (T_Property), pointer              :: Property
        integer                                 :: nNodes, NodeID, i,j
        integer                                 :: nReaches, ReachID
        real                                    :: PercentageMaxVolume
        real                                    :: PoolDepth, PoolVolume        


do1:    do i = 1, BaseTimeSeries

            select case (Me%TimeSerie%Name (i))            
                case (Char_WaterDepth)

                    nNodes = 1
                    nullify (CurrNode)

                    do NodeID = 1, Me%TotalNodes                            
                        CurrNode   => Me%Nodes (NodeID)                
                        if (CurrNode%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = CurrNode%WaterDepth            
                            nNodes = nNodes + 1
                        end if
                    enddo

                case (Char_WaterLevel)
           
                    nNodes = 1
                    nullify (CurrNode)

                    do NodeID = 1, Me%TotalNodes                            
                        CurrNode   => Me%Nodes (NodeID)                            
                        if (CurrNode%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = CurrNode%WaterLevel         
                            nNodes = nNodes + 1
                        end if
                    enddo    

                case (Char_PercentageMaxVolume)
        
                    nNodes = 1
                    nullify (CurrNode)

                    do NodeID = 1, Me%TotalNodes           
                       CurrNode   => Me%Nodes (NodeID)                            
                        if (CurrNode%TimeSerie) then
                            PercentageMaxVolume = CurrNode%VolumeNew / CurrNode%VolumeMax * 100.0
                            Me%TimeSerie%DataLine (nNodes) = PercentageMaxVolume
                            nNodes = nNodes + 1                
                        endif
                    enddo
  
                case (Char_VerticalArea)
      
                    nNodes = 1      
                    nullify (CurrNode)

                    do NodeID = 1, Me%TotalNodes                
                        CurrNode   => Me%Nodes (NodeID)               
                        if (CurrNode%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = CurrNode%VerticalArea
                            nNodes = nNodes + 1                
                        endif
                    enddo
  
                case (Char_FlowToChannels)
        
                    nNodes = 1      
                    nullify (CurrNode)
                
                    do NodeID = 1, Me%TotalNodes                
                        if (Me%Nodes(NodeID)%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = Me%RunOffVector (NodeID)
                            nNodes = nNodes + 1                
                        endif
                    enddo
  
                case (Char_Volume)

                    nNodes = 1      
                    nullify (CurrNode)

                    do NodeID = 1, Me%TotalNodes                
                        CurrNode   => Me%Nodes (NodeID)                
                        if (CurrNode%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = CurrNode%VolumeNew
                            nNodes = nNodes + 1                
                        endif
                    enddo


                case (Char_Flow)
        
                    nReaches = 1      
                    nullify (CurrReach)
      
                    do ReachID = 1, Me%TotalReaches
                        CurrReach  => Me%Reaches (ReachID)
                        if (CurrReach%TimeSerie) then                 
                            Me%TimeSerie%DataLine (nReaches) = CurrReach%FlowNew
                            nReaches = nReaches + 1                
                        endif
                    enddo
              
                 case (Char_Velocity) 
        
                    nReaches = 1      
                    nullify (CurrReach)

                    do ReachID = 1, Me%TotalReaches
                        CurrReach  => Me%Reaches (ReachID)        
                        if (CurrReach%TimeSerie) then                 
                            Me%TimeSerie%DataLine (nReaches) = CurrReach%Velocity
                            nReaches = nReaches + 1                            
                        endif
                    enddo

                case (Char_GWFlowToChannels)

                    nNodes = 1      
                    nullify (CurrNode)
                
                    do NodeID = 1, Me%TotalNodes                      
                        if (Me%Nodes(NodeID)%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = Me%GroundVector (NodeID)
                            nNodes = nNodes + 1                
                        endif
                    enddo
  
                case (Char_PoolDepth) 


                    nNodes = 1      
                    nullify (CurrNode)
                
                    do NodeID = 1, Me%TotalNodes                
                        if (Me%Nodes(NodeID)%TimeSerie) then
                           CurrNode   => Me%Nodes (NodeID)                
                            !Volume greater then volume in pools
                            PoolVolume = CurrNode%CrossSection%PoolDepth * CurrNode%Length * CurrNode%CrossSection%BottomWidth
                            if (CurrNode%VolumeNew > PoolVolume) then
                                PoolDepth = CurrNode%CrossSection%PoolDepth
                            else
                                PoolDepth = CurrNode%VolumeNew / (CurrNode%Length * CurrNode%CrossSection%BottomWidth)
                            endif
                            Me%TimeSerie%DataLine  (nNodes) = PoolDepth
                            nNodes = nNodes + 1                
                        endif
                    enddo

                case (Char_DT)

                    nNodes = 1      
                    nullify (CurrNode)
                
                    do NodeID = 1, Me%TotalNodes                
                        if (Me%Nodes(NodeID)%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = Me%ExtVar%DT
                            nNodes = nNodes + 1                
                        endif
                    enddo

                case (Char_DTLocal)

                    nNodes = 1      
                    nullify (CurrNode)                
                    do NodeID = 1, Me%TotalNodes                
               
                        if (Me%Nodes(NodeID)%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = LocalDT
                            nNodes = nNodes + 1                
                        endif
                    enddo

                case default
            
                    stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR01'  
        
            end select
        
            call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),          &
                                    DataLine = Me%TimeSerie%DataLine,      &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR01'       

        end do do1

        i = BaseTimeSeries

        !OutputHydro ----------------------------------------------------------

if2:    if (Me%OutputHydro) then

do2:        do j = i+1, i+5

                select case (Me%TimeSerie%Name (j))

                    case (Char_HydroTimeGradient) 
        
                        nReaches = 1      
                        nullify (CurrReach)

                        do ReachID = 1, Me%TotalReaches            
                            CurrReach  => Me%Reaches (ReachID)              
                            if (CurrReach%TimeSerie) then                 
                                Me%TimeSerie%DataLine (nReaches) = CurrReach%HydroTimeGradient
                                nReaches = nReaches + 1                            
                            endif        
                        enddo
      
                    case (Char_HydroAdvection) 
        
                        nReaches = 1      
                        nullify (CurrReach)

                        do ReachID = 1, Me%TotalReaches
                            CurrReach  => Me%Reaches (ReachID)        
                            if (CurrReach%TimeSerie) then                 
                                Me%TimeSerie%DataLine (nReaches) = CurrReach%HydroAdvection
                                nReaches = nReaches + 1                            
                            endif
                        enddo

                    case (Char_HydroPressure) 
        
                        nReaches = 1      
                        nullify (CurrReach)

                        do ReachID = 1, Me%TotalReaches
                            CurrReach  => Me%Reaches (ReachID)        
                            if (CurrReach%TimeSerie) then                 
                                Me%TimeSerie%DataLine (nReaches) = CurrReach%HydroPressure
                                nReaches = nReaches + 1                            
                            endif
                        enddo

                    case (Char_HydroGravity) 
        
                        nReaches = 1      
                        nullify (CurrReach)

                        do ReachID = 1, Me%TotalReaches
                            CurrReach  => Me%Reaches (ReachID)        
                            if (CurrReach%TimeSerie) then                 
                                Me%TimeSerie%DataLine (nReaches) = CurrReach%HydroGravity
                                nReaches = nReaches + 1                            
                            endif
                        enddo

                    case (Char_HydroFriction) 
        
                        nReaches = 1      
                        nullify (CurrReach)

                        do ReachID = 1, Me%TotalReaches
                            CurrReach  => Me%Reaches (ReachID)        
                            if (CurrReach%TimeSerie) then                 
                                Me%TimeSerie%DataLine (nReaches) = CurrReach%HydroFriction
                                nReaches = nReaches + 1                            
                            endif
                        enddo

                    case default
            
                        write(*,*) trim(Me%TimeSerie%Name(j)),'is not OutputHydro'
                        stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR02'  

                end select                            

                call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(j),          &
                                        DataLine = Me%TimeSerie%DataLine,      &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR02'       

            enddo do2
            
            i=i+5

        endif if2


        !Properties-------------------------------------------------------------
if3:    if (Me%TimeSerie%nProp .GT.BaseTimeSeries) then

            Property => Me%FirstProperty
            i = i  + 1
            do while (associated (Property))
                
                if (Property%ComputeOptions%TimeSerie) then
                    
                    nNodes = 1      
                    nullify (CurrNode)
                    
                    do NodeID = 1, Me%TotalNodes                
               
                        if (Me%Nodes(NodeID)%TimeSerie) then
                            Me%TimeSerie%DataLine  (nNodes) = Property%Concentration (NodeID)
                            nNodes = nNodes + 1                
                        endif
                    enddo
                   
                    call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),                   &
                                            DataLine = Me%TimeSerie%DataLine,              &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR03'  
                    
                    i = i + 1 

ifB:                if (Property%ComputeOptions%BottomFluxes) then
                        
                        nNodes = 1      
                        nullify (CurrNode)
                    
                        do NodeID = 1, Me%TotalNodes                
               
                            if (Me%Nodes(NodeID)%TimeSerie) then
                                Me%TimeSerie%DataLine  (nNodes) = Property%BottomConc (NodeID)
                                nNodes = nNodes + 1                
                            endif
                        enddo                    

                        call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),               &
                                                DataLine = Me%TimeSerie%DataLine,          &
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR04'  

                        i = i + 1

                    end if ifB

ifTo:               if (Property%ComputeOptions%SumTotalConc) then
                        
                        nNodes = 1      
                        nullify (CurrNode)
                    
                        do NodeID = 1, Me%TotalNodes                

                            if (Me%Nodes(NodeID)%TimeSerie) then
                                Me%TimeSerie%DataLine  (nNodes) = Property%TotalConc (NodeID)
                                nNodes = nNodes + 1                
                            endif
                        enddo                    

                        call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),               &
                                                DataLine = Me%TimeSerie%DataLine,          &
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - Wassim_ERR05'  

                        i = i + 1

                    end if ifTo
                    
ifLoad:             if (Property%ComputeOptions%ComputeLoad) then
                        
                        nNodes = 1      
                        nullify (CurrNode)
                    
                        do NodeID = 1, Me%TotalNodes                
               
                            if (Me%Nodes(NodeID)%TimeSerie) then
                                Me%TimeSerie%DataLine  (nNodes) = Property%Load (NodeID)
                                nNodes = nNodes + 1                
                            endif
                        enddo                    

                        call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),           &
                                                DataLine = Me%TimeSerie%DataLine,          &
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - Load_ERR0_Wassim'  

                        i = i + 1

                    end if ifLoad
                    
                                        
ifTox:              if (Property%ComputeOptions%Toxicity) then
                        
                        nNodes = 1      
                        nullify (CurrNode)
                    
                        do NodeID = 1, Me%TotalNodes                
               
                            if (Me%Nodes(NodeID)%TimeSerie) then
                                Me%TimeSerie%DataLine  (nNodes) = Property%Toxicity%Field (NodeID)
                                nNodes = nNodes + 1                
                            endif
                        enddo                    

                        call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i + 1),           &
                                                DataLine = Me%TimeSerie%DataLine,          &
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR06'  

                        i = i + 1

                    end if ifTox

                end if
                                
                Property => Property%Next
            
            end do

ifB2:     if (Me%ComputeOptions%BottomFluxes) then

                nNodes = 1      
                nullify (CurrNode)
            
                do NodeID = 1, Me%TotalNodes                
       
                    if (Me%Nodes(NodeID)%TimeSerie) then
                        CurrNode => Me%Nodes (NodeID)
                        ReachID = CurrNode%DownstreamReaches (1)
                        Me%TimeSerie%DataLine  (nNodes) = Me%ShearStress (ReachID)
                        nNodes = nNodes + 1                
                    endif
                enddo                    

                call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),                       &
                                        DataLine = Me%TimeSerie%DataLine,                 &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR07'  

                i = i + 1
            
            end if ifB2

        
ifTox2:     if (Me%ComputeOptions%Toxicity) then

                nNodes = 1      
                nullify (CurrNode)
            
                do NodeID = 1, Me%TotalNodes                
       
                    if (Me%Nodes(NodeID)%TimeSerie) then
                        Me%TimeSerie%DataLine  (nNodes) = Me%GlobalToxicity (NodeID)
                        nNodes = nNodes + 1                
                    endif
                enddo                    

                call WriteTimeSerieLine(Me%TimeSerie%ObjTimeSerie(i),                      &
                                        DataLine = Me%TimeSerie%DataLine,                  &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteTimeSeriesByProp - ERR08'  

                i = i + 1
            
            end if ifTox2
                
        end if if3

        nullify (Property) 


    end subroutine WriteTimeSeriesByProp
    
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine FillOutPutMatrix  (DrainageNetworkID, OutputMatrix,              &
                                  OutPutResultsType, Deposited, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: DrainageNetworkID          
        real   , dimension(:,:), pointer            :: OutputMatrix
        integer                                     :: OutPutResultsType
        logical, intent(IN ) , optional             :: Deposited
        integer, intent(OUT) , optional             :: STAT

        !Local------------------------------------------------------------------
        type (T_Node    ), pointer                  :: CurrNode, OtherNode        
        type (T_Property), pointer                  :: Property
        integer                                     :: STAT_CALL, ready_
        integer                                     :: NodeID, ReachID, OtherNodeID
        real                                        :: Value
        real                                        :: FlowValue, Velocity
        real                                        :: dx, dy, angle
        

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)
        
if0:    if (ready_ .EQ. IDLE_ERR_) then

            OutputMatrix = null_real

            select case (OutPutResultsType)

            case (WaterDepth_)

                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)                    
                    OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = CurrNode%WaterDepth
                end do

            case (WaterLevel_)
                                            
                    OutputMatrix = Me%ChannelsWaterLevel
            
            case (Volume_)

                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)         
                    OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = CurrNode%VolumeNew                   
                end do
    
            case (PercentageMaxVolume_)

                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)
                    Value = CurrNode%VolumeNew / CurrNode%VolumeMax * 100.0
                    OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Value                   
                end do

            case (WaterFluxX_, WaterFluxY_, VelocityU_, VelocityV_)
                  
                !Flow in cell corresponds to the cell outflow
                !Because flow is computed with upstream node properties               
                do NodeID = 1, Me%TotalNodes
                    
                    CurrNode => Me%Nodes (NodeID)
                    
                    Value = 0.0

if2:                if (CurrNode%nDownstreamReaches .NE. 0) then
                    
                        ReachID     = CurrNode%DownstreamReaches (1)
                        OtherNodeID = Me%Reaches(ReachID)%DownstreamNode
                        FlowValue   = Me%Reaches(ReachID)%FlowNew
                        Velocity    = Me%Reaches(ReachID)%Velocity
                                                                                                                    
                        OtherNode => Me%Nodes(OtherNodeID)

                        dx = OtherNode%X - CurrNode%X
                        dy = OtherNode%Y - CurrNode%Y

                        !if (CurrNode%nDownstreamReaches == 0) then                            
                        !    dx = - dx
                        !    dy = - dy
                        !end if

                        if (dx == 0.0) then
                            angle = Pi / 2.0    
                        else
                            angle = atan(dy / dx)                                                                       
                        end if

                        if (dx <= 0.0 .AND. dy < 0.0) then
                           angle = Pi + abs(angle)
                        else if (dx < 0.0 .AND. dy >= 0.0) then
                           angle = Pi - abs(angle)
                        else if (dx > 0.0 .AND. dy < 0.0) then
                           angle = 2.0 * Pi - abs(angle)
                        end if

                        if (OutPutResultsType == WaterFluxX_) then
                            Value = FlowValue * cos (angle)                            
                        else if (OutPutResultsType == WaterFluxY_) then
                            Value = FlowValue * sin (angle)                            
                        else if (OutPutResultsType == VelocityU_) then
                            Value = Velocity * cos (angle)                           
                        else if (OutPutResultsType == VelocityV_) then
                            Value = Velocity * sin (angle)
                        end if
                                                
                        OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Value

                    end if if2

                end do

            case (FlowModulus_)

                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)

                    if (CurrNode%nDownstreamReaches .NE. 0) then      
                                      
                        ReachID = CurrNode%DownstreamReaches (1)
                        OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Me%Reaches(ReachID)%FlowNew
                    
                    end if
                    
                end do
            
            case (VelocityModulus_)

                do NodeID = 1, Me%TotalNodes
                    
                    CurrNode => Me%Nodes (NodeID)
                    if (CurrNode%nDownstreamReaches .NE. 0) then                                            
                        ReachID = CurrNode%DownstreamReaches (1)
                        OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Me%Reaches(ReachID)%Velocity
                    end if
                    
                end do
            
            case (GenericProperty_)

                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)                    
                    OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Me%GlobalToxicity (NodeID)
                end do

            case (ShearStress_)

                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)
                    if (CurrNode%nDownstreamReaches .NE. 0) then                                            
                        ReachID = CurrNode%DownstreamReaches (1)
                        OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Me%ShearStress (ReachID)
                    end if
                end do


            case default
                
                !Find Property with the given ID                
                nullify (Property)
                Property => Me%FirstProperty
                do while (associated (Property))
                    
                    if (present (Deposited)) then

                        if (Property%ID%IDNumber == OutPutResultsType) then
                            !Fills Matrix
                            do NodeID = 1, Me%TotalNodes                    
                                CurrNode => Me%Nodes (NodeID)
                                OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Property%BottomConc (NodeID)
                            end do
                            exit
                        end if

                    else

                        if (Property%ID%IDNumber == OutPutResultsType) then
                            !Fills Matrix
                            do NodeID = 1, Me%TotalNodes                    
                                CurrNode => Me%Nodes (NodeID)
                                OutputMatrix (CurrNode%GridI, CurrNode%GridJ) = Property%Concentration (NodeID)
                            end do
                            exit
                        end if     
                     
                    end if                

                    Property => Property%Next
                enddo

                
            
            end select


            STAT_CALL = SUCCESS_
        else               
            STAT_CALL = ready_
        end if if0

        if (present(STAT)) STAT = STAT_CALL

    end subroutine FillOutPutMatrix


    !--------------------------------------------------------------------------

    subroutine HDF5Output

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL           
        real, dimension(6), target                      :: AuxTime
        real, dimension(:), pointer                     :: TimePointer
        real, dimension(:), pointer                     :: OutputMatrix
        real, dimension(:), pointer                     :: OutputMatrix_bottom
        real, dimension(:), pointer                     :: OutputMatrix_totalconc
        integer                                         :: iReach, NodeID
        type (T_Node), pointer                          :: CurrNode
        type (T_Property), pointer                      :: Property


        if (Me%CurrentTime >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

            allocate (OutputMatrix          (1:Me%TotalReaches))
            allocate (OutputMatrix_bottom   (1:Me%TotalReaches))
            allocate (OutputMatrix_totalconc(1:Me%TotalReaches))

            !Writes current time
            call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2),         &
                                                AuxTime(3), AuxTime(4),         &
                                                AuxTime(5), AuxTime(6))
            TimePointer => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR01'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                 "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D      = TimePointer,                    &
                                 OutputNumber = Me%OutPut%NextOutPut,           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR02'

            !Everything will be written for reaches (graphical reasons)
            !If variable is given by node, use upstream node to set reach property

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, 1, Me%TotalReaches, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleDrainageNetwork - ERR03'
            
            !Writes Reach ID
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = Me%Reaches(CurrNode%DownstreamReaches(1))%ID
                    iReach = iReach + 1
                end if
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/ID", "ReachID",                           &
                                 "-",                                                    &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR04'


            !Writes Flow Modulus
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = Me%Reaches(CurrNode%DownstreamReaches(1))%FlowNew
                    iReach = iReach + 1
                end if
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/channel flow", "channel flow",    &
                                 "m3/s",                                                 &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR04'
            

            !Writes Velocity Modulus
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = Me%Reaches(CurrNode%DownstreamReaches(1))%Velocity
                    iReach = iReach + 1
                end if
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/Velocity", "Velocity",            &
                                 "m/s",                                                  &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR05'


            !Writes Flow accumulation time in percentage
            if (Me%Output%ComputeFlowFrequency) then
                iReach = 1
                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)
                    if (CurrNode%nDownstreamReaches .NE. 0) then      
                        OutputMatrix (iReach) = Me%Reaches(CurrNode%DownstreamReaches(1))%FlowAccPerc
                        iReach = iReach + 1
                    end if
                end do
                call HDF5WriteData  (Me%ObjHDF5, "/Results/flow frequency", "flow frequency",    &
                                     "-",                                                 &
                                     Array1D      = OutputMatrix,                            &
                                     OutputNumber = Me%OutPut%NextOutPut,                    &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR05.1'
            end if

            !Advection
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = Me%Reaches(CurrNode%DownstreamReaches(1))%HydroAdvection
                    iReach = iReach + 1
                end if
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/Advection", "Advection",           &
                                 "m3/s",                                                 &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR05a'

            !Pressure
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = Me%Reaches(CurrNode%DownstreamReaches(1))%HydroPressure
                    iReach = iReach + 1
                end if
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/Pressure", "Pressure",            &
                                 "m3/s",                                                 &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR05b'
            
            !Friction
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = Me%Reaches(CurrNode%DownstreamReaches(1))%HydroFriction
                    iReach = iReach + 1
                end if
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/Friction", "Friction",            &
                                 "-",                                                    &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR05c'

            !Waterdepth
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then
                    OutputMatrix (iReach) = CurrNode%Waterdepth 
                    iReach = iReach + 1
                endif
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/channel water depth",             &
                                 "channel water depth",                                  &  
                                 "m",                                                    &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR06'

            !WaterLevel
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then
                    OutputMatrix (iReach) = CurrNode%WaterLevel 
                    iReach = iReach + 1
                endif
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/channel water level",             &
                                 "channel water level",                                  &  
                                 "m",                                                    &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR06'

            ! Volume
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = CurrNode%VolumeNew
                    iReach = iReach + 1
                endif
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/Volume", "Volume",                &
                                 "m3",                                                   &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR07'


            !%Percentagem Volume
            iReach = 1
            do NodeID = 1, Me%TotalNodes
                CurrNode => Me%Nodes (NodeID)
                if (CurrNode%nDownstreamReaches .NE. 0) then      
                    OutputMatrix (iReach) = CurrNode%VolumeNew / CurrNode%VolumeMax * 100.0
                    iReach = iReach + 1
                endif
            end do
            call HDF5WriteData  (Me%ObjHDF5, "/Results/percentage max volume",           &
                                 "percentage max volume",                                &
                                 "-",                                                    &
                                 Array1D      = OutputMatrix,                            &
                                 OutputNumber = Me%OutPut%NextOutPut,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR07a'

            !Shear Stress
            if (Me%ComputeOptions%BottomFluxes) then
                iReach = 1
                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)
                    if (CurrNode%nDownstreamReaches .NE. 0) then      
                        OutputMatrix (iReach) = Me%ShearStress(iReach)
                        iReach = iReach + 1
                    end if
                end do
                call HDF5WriteData  (Me%ObjHDF5, "/Results/Shear Stress", "Shear Stress",    &
                                    "Pa",                                                   &
                                    Array1D      = OutputMatrix,                            &
                                    OutputNumber = Me%OutPut%NextOutPut,                    &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR08'
            endif
            
         
            !Properties
            nullify (Property)
            Property => Me%FirstProperty
            do while (associated (Property))
                iReach = 1
                do NodeID = 1, Me%TotalNodes
                    CurrNode => Me%Nodes (NodeID)
                    if (CurrNode%nDownstreamReaches .NE. 0) then      
                        OutputMatrix (iReach) = Property%Concentration (iReach)            !MO: IReach or NodeID???                           

                        if(Property%ComputeOptions%BottomFluxes) then  
                            OutputMatrix_bottom (iReach) = Property%BottomConc (iReach)    !MO: IReach or NodeID???          
                        endif

                        if(Property%ComputeOptions%SumTotalConc) then
                            OutputMatrix_totalconc (iReach) = Property%TotalConc (iReach)  !MO: IReach or NodeID???
                        endif

                        iReach = iReach + 1
                    endif
                end do

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//Property%ID%Name,      &
                                Property%ID%Name    ,                                   &
                                "-",                                                    &
                                Array1D      = OutputMatrix,                            &
                                OutputNumber = Me%OutPut%NextOutPut,                    &
                                STAT = STAT_CALL)

                if(Property%ComputeOptions%BottomFluxes) then  
                    call HDF5WriteData  (Me%ObjHDF5,                                 &
                                "/Results/Bottom_"//Property%ID%Name,                   &
                                "Bottom_"//Property%ID%Name,                            &
                                "-",                                                    &
                                Array1D      = OutputMatrix_bottom,                            &
                                OutputNumber = Me%OutPut%NextOutPut,                    &
                                STAT = STAT_CALL)
                endif

                if(Property%ComputeOptions%SumTotalConc) then  
                    call HDF5WriteData  (Me%ObjHDF5,                                 &
                                "/Results/TotalConc_"//Property%ID%Name,                &
                                "TotalConc_"//Property%ID%Name,                         &
                                "-",                                                    &
                                Array1D      = OutputMatrix_totalconc,                            &
                                OutputNumber = Me%OutPut%NextOutPut,                    &
                                STAT = STAT_CALL)
                endif

                if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR51'                 

                Property => Property%Next
            enddo
        


            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1
            deallocate (OutputMatrix)
            deallocate (OutputMatrix_bottom)
            deallocate (OutputMatrix_totalconc)

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleDrainageNetwork - ERR99'

        endif 
        
    end subroutine HDF5Output


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillDrainageNetwork(DrainageNetworkID, STAT)

        !Arguments--------------------------------------------------------------
        integer                             :: DrainageNetworkID              
        integer, optional, intent(OUT)      :: STAT

        !External---------------------------------------------------------------
        integer                             :: ready_              

        !Local------------------------------------------------------------------
        integer                             :: STAT_CALL, nUsers , i
        type (T_Property), pointer          :: Property

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(DrainageNetworkID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            call WriteFinalFile
                        
            call Write_Errors_Messages            
                                                           
            nUsers = DeassociateInstance(mDrainageNetwork_,  Me%InstanceID)

            
            if (nUsers == 0) then
                               
                if (Me%WriteMaxStationValues) call MaxStationValuesOutput

!                !Writes a last output of the time series
!                if (Me%TimeSerie%nNodes .GT.0) then 
!                    call WriteTimeSeries (0.0)
!                    write(*,*)'Writing Time Series Last Time'
!                endif
                
                deallocate (Me%RunOffVector     )
                deallocate (Me%GroundVector     )
                deallocate (Me%DiffuseVector    )
                deallocate (Me%Nodes            )
                deallocate (Me%Reaches          )
                deallocate (Me%ComputeFaces     )
                deallocate (Me%OpenPointsFlow   )
                deallocate (Me%OpenPointsProcess)

                if (Me%ComputeOptions%Toxicity) deallocate (Me%GlobalToxicity)

                if (Me%ObjTime /= 0) then
                    nUsers = DeassociateInstance(mTIME_, Me%ObjTime)
                    if (nUsers == 0) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR01'
                endif
                
                !Kill LightExtinction
                if (Me%ComputeOptions%TopRadiation) then
                    call KillLightExtinction(Me%ObjLightExtinction, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR10' 
                endif

                !Kill WaterQuality    
                if (Me%ComputeOptions%WaterQuality) then
                    call KillInterface (Me%ObjInterface, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR20' 
                endif

                !Kill Benthos  
                if (Me%ComputeOptions%Benthos) then
                    call KillInterface (Me%ObjBenthicInterface, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR30' 
                endif



                !Kill Properties
                if (Me%HasProperties) then

                    Property => Me%FirstProperty
                    do while (associated (Property))

                        deallocate (Property%Concentration           )
                        deallocate (Property%MassCreated             )
                        deallocate (Property%ConcentrationOld        )
                        deallocate (Property%InitialConcentration    )
                        deallocate (Property%InitialConcentrationOld )
                        deallocate (Property%OverLandConc            )
                        deallocate (Property%GWaterConc              )
                        deallocate (Property%TotalConc               )
                        deallocate (Property%Load                    )
                        deallocate (Property%MassInKg                )
                        
                        if(Property%ComputeOptions%BottomFluxes)then
                            deallocate (Property%BottomConc     )
                            deallocate (Property%ErosionRate    )
                            deallocate (Property%DepositionRate )
                            deallocate (Property%Ws             )
                        endif

                        if (Property%ComputeOptions%Toxicity) deallocate (Property%Toxicity%Field)
                        Property => Property%Next
                    end do
                    nullify (Property)
                end if

                !Kills TimeSerie Objects                
                if (Me%TimeSerie%nNodes/=0) then

                    deallocate (Me%TimeSerie%DataLine)

                    if (Me%TimeSerie%ByNodes) then

                        do i = 1, Me%TimeSerie%nNodes
                            call KillTimeSerie(Me%TimeSerie%ObjTimeSerie(i), STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR40'
                        end do                    

                    else

                        do i = 1, Me%TimeSerie%nProp
                            call KillTimeSerie(Me%TimeSerie%ObjTimeSerie(i), STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR50'
                        end do                    

                    end if
                endif         
                               
                deallocate (Me%TimeSerie%ObjTimeSerie)
                deallocate (Me%TimeSerie%Name        )

                if (Me%Downstream%Boundary == ImposedWaterLevel &
                    .AND. Me%Downstream%Evolution == ReadTimeSerie ) then
                    
                    call KillTimeSerie(Me%Downstream%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR60'

                end if


                !KillDischarges
                if (Me%ComputeOptions%Discharges) then

                    nUsers = GetUsersNumber(mDISCHARGES_, Me%ObjDischarges)
                    if (nUsers == 1) then
            
                        call Kill_Discharges(Me%ObjDischarges, STAT = STAT_CALL)
    
                        if (STAT_CALL /= SUCCESS_)                                               &
                            stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR70'

                    else if (nUsers > 1) then

                        nUsers = DeassociateInstance (mDISCHARGES_, Me%ObjDischarges)
                        if (nUsers == 0) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR80'
               
                    else
            
                        stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR90'

                    endif

                endif
                
                if (Me%Output%Yes) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillDrainageNetwork - ModuleDrainageNetwork - ERR100'
                endif            
                
                                                           
                !Deallocates Instance
                call DeallocateInstance ()

                DrainageNetworkID = 0
                STAT_CALL      = SUCCESS_

            end if


        end if cd1


        if (present(STAT)) STAT = STAT_CALL     

    end subroutine KillDrainageNetwork
        
    !---------------------------------------------------------------------------

    subroutine MaxStationValuesOutput


        !Local-----------------------------------------------------------------
        integer                                 :: NodeID, STAT_CALL
        type(T_Node) , pointer                  :: CurrNode
        character(len=PathLength)               :: File
        

        call ReadFileName("ROOT_SRT", File, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MaxStationValuesOutput - ModuleDrainageNetwork - ERR01'
        File= trim(adjustl(File))//"StationsMaxFlow.dat"

        open(UNIT=UnitMax, FILE=File, ACTION='WRITE', STATUS='REPLACE', IOSTAT=STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MaxStationValuesOutput - ModuleDrainageNetwork - ERR02'

        write(UnitMax,*) 'StationName             X              Y        Depth   Flow    Velocity  Time'

        do NodeID = 1, Me%TotalNodes

            CurrNode => Me%Nodes(NodeID)

            if (CurrNode%StationName /= null_str) then

                write(UnitMax,100) CurrNode%StationName, CurrNode%X, CurrNode%Y,               &
                                   CurrNode%Max%Depth, CurrNode%Max%Flow, CurrNode%Max%Vel,    &
                                   CurrNode%Max%Time
            endif

        enddo

        close(UnitMax)
        nullify(CurrNode)

        100 format(A15,1x, 2f15.4, 3e15.4,3x, A19)   

    end subroutine MaxStationValuesOutput

    !---------------------------------------------------------------------------

    subroutine WriteFinalFile

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: FinalFile
        integer                                     :: STAT_CALL
        type (T_Property), pointer                  :: Property
        character(LEN = PathLength)                 :: FileName

        !-----------------------------------------------------------------------

        if (Me%CurrentTime == Me%EndTime) then
            FileName = Me%Files%FinalFile
        else
            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%CurrentTime))//".fin")
        endif            

        call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteFinalFile - ERR01'

        open(Unit = FinalFile, File = FileName, Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteFinalFile - ERR02'
        
        call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleDrainageNetwork - WriteFinalFile - ERR03'

        !Writes Date
        call ExtractDate(Me%CurrentTime, Year_File, Month_File, Day_File,            &
                         Hour_File, Minute_File, Second_File)

        write(FinalFile) Year_File, Month_File, Day_File, Hour_File, Minute_File,    &
                         Second_File

        write(FinalFile)Me%Nodes%WaterLevel
        write(FinalFile)Me%Reaches%FlowNew

        Property => Me%FirstProperty
        do while (associated(Property))
            write (FinalFile) Property%Concentration
            Property => Property%Next
        end do

        write(FinalFile)Me%LastGoodNiter

        Property => Me%FirstProperty
        do while (associated(Property))
            write (FinalFile) Property%BottomConc
            Property => Property%Next
        end do


        call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ModuleDrainageNetwork - WriteFinalFile - ERR04'


    end subroutine WriteFinalFile
    
    !---------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        type (T_DrainageNetwork), pointer          :: AuxObjDrainageNetwork
        type (T_DrainageNetwork), pointer          :: PreviousObjDrainageNetwork

        !Updates pointers
        if (Me%InstanceID == FirstDrainageNetwork%InstanceID) then
            FirstDrainageNetwork => FirstDrainageNetwork%Next
        else
            PreviousObjDrainageNetwork => FirstDrainageNetwork
            AuxObjDrainageNetwork      => FirstDrainageNetwork%Next
            do while (AuxObjDrainageNetwork%InstanceID /= Me%InstanceID)
                PreviousObjDrainageNetwork => AuxObjDrainageNetwork
                AuxObjDrainageNetwork      => AuxObjDrainageNetwork%Next
            enddo

            !Now update linked list
            PreviousObjDrainageNetwork%Next => AuxObjDrainageNetwork%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !---------------------------------------------------------------------------

    subroutine Write_Errors_Messages

        !Local------------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        real                                        :: Total_Mass_Created, Total_Mass
        character (Len = StringLength)              :: str_mass, string_to_be_written
        integer                                     :: NodeID
        type (T_Node), pointer                      :: CurrNode
        real                                        :: BottomArea

            !TotalVolume entering in DNet (by discharges, by Runoff, by GWFlow)
            !TotalWaterVolume left in DNet
            !TotalVolume out of DNet
            
            Property => Me%FirstProperty
            do while (associated (Property))
                
                !Total mass created due to imposed minimum concentration
                !kg                
                Total_Mass_Created = SUM(Property%MassCreated)
                
                str_mass = ''
                
                write(str_mass, '(f20.8)') Total_Mass_Created
      
                string_to_be_written = 'Total mass created on property '                //   &
                                        trim(adjustl(adjustr(Property%ID%name)))//' = ' //   &
                                        trim(adjustl(adjustr(str_mass))) 

                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

                                
                !Total mass of property in the channels
                !kg
                Total_Mass = 0.0
                
                do NodeID = 1, Me%TotalNodes
                    if (Me%Nodes(NodeID)%nDownStreamReaches /= 0) then
                        Total_Mass = Total_Mass + Property%Concentration (NodeID)           &
                                                * Property%ISCoefficient                    &
                                                * Me%Nodes(NodeID)%VolumeNew
                    endif
                end do
                
                
                str_mass = ''
                
                write(str_mass, '(f20.8)') Total_Mass
      
                string_to_be_written = 'Total mass left of property '                   //   &
                                        trim(adjustl(adjustr(Property%ID%name)))//' = ' //   &
                                        trim(adjustl(adjustr(str_mass))) 

                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

                if(Property%ComputeOptions%BottomFluxes)then

                    !Total mass of deposited property in the channels
                    !kg
                    Total_Mass = 0.0

                    do NodeID = 1, Me%TotalNodes
                        if (Me%Nodes(NodeID)%nDownStreamReaches /= 0) then
                            CurrNode => Me%Nodes (NodeID)
                            BottomArea = CurrNode%CrossSection%BottomWidth * CurrNode%Length
                            Total_Mass = Total_Mass + Property%BottomConc (NodeID) * BottomArea
                        endif
                    end do

                    str_mass = ''
                
                    write(str_mass, '(f20.8)') Total_Mass
      
                    string_to_be_written = 'Total mass left of deposited property '         //   &
                                            trim(adjustl(adjustr(Property%ID%name)))//' = ' //   &
                                            trim(adjustl(adjustr(str_mass))) 

                    call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

                end if

                Property => Property%Next
            
            end do

    end subroutine Write_Errors_Messages
                                            
    !---------------------------------------------------------------------------
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMEN

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !---------------------------------------------------------------------------

    subroutine Ready (ObjDrainageNetwork_ID, ready_) 

        !Arguments--------------------------------------------------------------
        integer                                     :: ObjDrainageNetwork_ID
        integer                                     :: ready_

        !-----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjDrainageNetwork_ID > 0) then
            call LocateObjDrainageNetwork (ObjDrainageNetwork_ID)
            ready_ = VerifyReadLock (mDrainageNetwork_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !-----------------------------------------------------------------------

    end subroutine Ready

    !---------------------------------------------------------------------------

    subroutine LocateObjDrainageNetwork (DrainageNetworkID)

        !Arguments--------------------------------------------------------------
        integer                                     :: DrainageNetworkID

        !Local------------------------------------------------------------------

        Me => FirstDrainageNetwork
        do while (associated (Me))
            if (Me%InstanceID == DrainageNetworkID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleDrainageNetwork - LocateObjDrainageNetwork - ERR01'

    end subroutine LocateObjDrainageNetwork

    !---------------------------------------------------------------------------

#ifdef _OPENMI_

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfNodes
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFNODES"::GetNumberOfNodes
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetNumberOfNodes(DrainageNetworkID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetNumberOfNodes = Me%TotalNodes
        else 
            GetNumberOfNodes = - 99.0
        end if
           
    end function GetNumberOfNodes

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetOutletNodeID
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETOUTLETNODEID"::GetOutletNodeID
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetOutletNodeID(DrainageNetworkID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetOutletNodeID = Me%OutletNodePos
        else 
            GetOutletNodeID = - 99.0
        end if

    end function GetOutletNodeID

    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetXCoordinate
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETXCOORDINATE"::GetXCoordinate
    !DEC$ ENDIF
    real(8) function GetXCoordinate(DrainageNetworkID, NodeID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: NodeID
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetXCoordinate = dble(Me%Nodes(NodeID)%X)
        else
            GetXCoordinate = -99.0
        endif
        

    end function GetXCoordinate

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetYCoordinate
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETYCOORDINATE"::GetYCoordinate
    !DEC$ ENDIF
    real(8) function GetYCoordinate(DrainageNetworkID, NodeID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: NodeID
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetYCoordinate = dble(Me%Nodes(NodeID)%Y)
        else
            GetYCoordinate = -99.0
        endif
        
    end function GetYCoordinate


    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetOutletFlow
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETOUTLETFLOW"::GetOutletFlow
    !DEC$ ENDIF
    real(8) function GetOutletFlow(DrainageNetworkID)

        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetOutletFlow = dble(Me%Reaches(Me%OutletReachPos)%FlowNew)  
        else
            GetOutletFlow = -99.0
        endif
        
    end function

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::SetDownStreamWaterLevel
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_SETDOWNSTREAMWATERLEVEL"::SetDownStreamWaterLevel
    !DEC$ ENDIF
    logical function SetDownStreamWaterLevel(DrainageNetworkID, WaterLevel)

        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        real(8)                                     :: WaterLevel
        
        !Local-----------------------------------------------------------------
        type (T_Node), pointer                      :: CurrNode

        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            Me%Downstream%DefaultValue = WaterLevel

            SetDownStreamWaterLevel = .true.
        
        else
        
            SetDownStreamWaterLevel = .false.
        
        endif        


        return

    end function SetDownStreamWaterLevel
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfOutFlowNodes
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFOUTFLOWNODES"::GetNumberOfOutFlowNodes
    !DEC$ ENDIF
    !Return the number of OutflowNodes
    integer function GetNumberOfOutFlowNodes(DrainageNetworkID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        
        !Local-----------------------------------------------------------------
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetNumberOfOutFlowNodes = Me%StormWaterModelLink%nOutflowNodes
        else 
            GetNumberOfOutFlowNodes = -99
        end if
    
    end function GetNumberOfOutFlowNodes
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfInFlowNodes
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFINFLOWNODES"::GetNumberOfInFlowNodes
    !DEC$ ENDIF
    !Return the number of InflowNodes
    integer function GetNumberOfInFlowNodes(DrainageNetworkID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        
        !Local-----------------------------------------------------------------
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetNumberOfInFlowNodes = Me%StormWaterModelLink%nInflowNodes
        else 
            GetNumberOfInFlowNodes = -99
        end if
    
    end function GetNumberOfInFlowNodes

    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetStormWaterOutFlow
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETSTORMWATEROUTFLOW"::GetStormWaterOutFlow
    !DEC$ ENDIF
    !Return the Storm Water Outflow
    logical function GetStormWaterOutFlow(DrainageNetworkID, nOutflowNodes, StormWaterOutflow)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: nOutflowNodes
        real, dimension(nOutflowNodes)              :: StormWaterOutflow
        
        !Local-----------------------------------------------------------------
        integer                                     :: iNode
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            do iNode = 1, nOutflowNodes 
                StormWaterOutflow(iNode) = Me%StormWaterModelLink%Outflow(iNode)
            enddo
        
            GetStormWaterOutFlow = .true.
        
        else 
        
            GetStormWaterOutFlow = .false.
        
        end if
    
    end function GetStormWaterOutFlow
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetStormWaterOutFlowIDs
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETSTORMWATEROUTFLOWIDS"::GetStormWaterOutFlowIDs
    !DEC$ ENDIF
    !Return the Storm Water Outflow IDs
    logical function GetStormWaterOutFlowIDs(DrainageNetworkID, nOutflowNodes, StormWaterOutflowIDs)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: nOutflowNodes
        integer, dimension(nOutflowNodes)           :: StormWaterOutflowIDs
        
        !Local-----------------------------------------------------------------
        integer                                     :: iNode
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            do iNode = 1, nOutflowNodes 
                StormWaterOutflowIDs(iNode) = Me%StormWaterModelLink%OutflowIDs(iNode)
            enddo
        
            GetStormWaterOutFlowIDs = .true.
        
        else 
        
            GetStormWaterOutFlowIDs = .false.
        
        end if
    
    end function GetStormWaterOutFlowIDs    

    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::SetStormWaterInFlow
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_SETSTORMWATERINFLOW"::SetStormWaterInFlow
    !DEC$ ENDIF
    !Return the Storm Water Outflow
    logical function SetStormWaterInFlow(DrainageNetworkID, nInflowNodes, StormWaterInflow)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: nInflowNodes
        real, dimension(nInflowNodes)               :: StormWaterInflow
        
        !Local-----------------------------------------------------------------
        integer                                     :: iNode
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            do iNode = 1, nInflowNodes 
                Me%StormWaterModelLink%Inflow(iNode) = StormWaterInflow(iNode)
            enddo
            
            SetStormWaterInFlow = .true.
        
        else 
        
            SetStormWaterInFlow = .false.
        
        end if
    
    end function SetStormWaterInFlow
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetStormWaterInFlowIDs
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETSTORMWATERINFLOWIDS"::GetStormWaterInFlowIDs
    !DEC$ ENDIF
    !Return the Storm Water Inflow IDs
    logical function GetStormWaterInFlowIDs(DrainageNetworkID, nInflowNodes, StormWaterInflowIDs)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: nInflowNodes
        integer, dimension(nInflowNodes)            :: StormWaterInflowIDs
        
        !Local-----------------------------------------------------------------
        integer                                     :: iNode
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            do iNode = 1, nInflowNodes 
                StormWaterInflowIDs(iNode) = Me%StormWaterModelLink%InflowIDs(iNode)
            enddo
        
            GetStormWaterInFlowIDs = .true.
        
        else 
        
            GetStormWaterInFlowIDs = .false.
        
        end if
    
    end function GetStormWaterInFlowIDs    


    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfProperties
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFPROPERTIES"::GetNumberOfProperties
    !DEC$ ENDIF
    !Return the number of Properties
    integer function GetNumberOfProperties(DrainageNetworkID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            GetNumberOfProperties = Me%PropertiesNumber
        else 
            GetNumberOfProperties = -99
        end if
           
    end function GetNumberOfProperties
    
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetDrainageNetworkPropertyID
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETDRAINAGENETWORKPROPERTYID"::GetDrainageNetworkPropertyID
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetDrainageNetworkPropertyID(DrainageNetworkID, idx)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: idx
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type (T_Property), pointer                  :: Property
        integer                                     :: iProp

        call Ready(DrainageNetworkID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            Property => Me%FirstProperty
            iProp = 1
            do while (associated (Property))
                 
                 if (iProp == idx) then
                 
                    GetDrainageNetworkPropertyID = Property%ID%IDNumber
                    return
                 
                 endif
                 
                 Property => Property%Next
                 iProp = iProp + 1
            enddo
        
            GetDrainageNetworkPropertyID = -99
        else 
            GetDrainageNetworkPropertyID = -99
        end if
           
    end function GetDrainageNetworkPropertyID
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetOutletFlowConcentration
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETOUTLETFLOWCONCENTRATION"::GetOutletFlowConcentration
    !DEC$ ENDIF
    real(8) function GetOutletFlowConcentration(DrainageNetworkID, PropIDNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: PropIDNumber
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type (T_Property), pointer                  :: Property
        type (T_Reach), pointer                     :: OutletReach

        call Ready(DrainageNetworkID, ready_)    

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            call SearchProperty (Property, PropIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then            
                OutletReach => Me%Reaches(Me%OutletReachPos)
                GetOutletFlowConcentration =  dble(Property%Concentration(OutletReach%UpstreamNode))
            else
                GetOutletFlowConcentration = -99.0
            endif
        else
            GetOutletFlowConcentration = -99.0
        endif
        
    end function    
    

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::SetDownStreamConcentration
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_SETDOWNSTREAMCONCENTRATION"::SetDownStreamConcentration
    !DEC$ ENDIF
    logical function SetDownStreamConcentration(DrainageNetworkID, PropIDNumber, Concentration)

        !Arguments-------------------------------------------------------------
        integer                                     :: DrainageNetworkID
        integer                                     :: PropIDNumber
        real(8)                                     :: Concentration
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         
        type (T_Property), pointer                  :: Property

        call Ready(DrainageNetworkID, ready_)    

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            
            call SearchProperty (Property, PropIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then            
                Property%Concentration(Me%OutletNodePos) = Concentration
                SetDownStreamConcentration = .true.
            else
                SetDownStreamConcentration = .false.
            endif
            
        
        else
        
            SetDownStreamConcentration = .false.
        
        endif        


        return

    end function SetDownStreamConcentration    
    
    

#endif


end module ModuleDrainageNetwork

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------











