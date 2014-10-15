!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : PorousMedia
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - Complete Revision, 
! DESCRIPTION   : Simulates Water Flow in variable saturated soils
!
!------------------------------------------------------------------------------

! Keywords read in the Data File
!
! Keyword                   : Data Type         Default     !Comment
!
! BOTTOM_FILE               : char              -           !Path to Bottom Topography File
! START_WITH_FIELD          : logical           1           !Sets Theta initial Field Capacity
! CONTINUOUS                : logical           0           !Continues from previous run
! STOP_ON_WRONG_DATE        : logical           1           !Stops if previous run end is different from actual
!                                                           !Start
! OUTPUT_TIME               : sec. sec. sec.    -           !Output Time
! SURFACE_OUTPUT_TIME       : sec. sec. sec.    -           !Output Time of surface layer
! TIME_SERIE_LOCATION       : char              -           !Path to File which defines Time Series
! CONTINUOUS_OUTPUT_FILE    : logical           1           !Writes "famous" iter.log
! CONDUTIVITYFACE           : integer           1           !Way to interpolate conducivity face
!                                                           !1 - Average, 2 - Maximum, 3 - Minimum, 4 - Weigthed
! HORIZONTAL_K_FACTOR       : real              1.0         !Factor for Horizontal Conductivity = Kh / Kv
! CUT_OFF_THETA_LOW         : real              1e-6        !Disables calculation when Theta is near ThetaR
! CUT_OFF_THETA_HIGH        : real              1e-15       !Set Theta = ThetaS when Theta > ThetaS - CUT_OFF_THETA_HIGH
! MIN_ITER                  : integer           2           !Number of iterations below which the DT is increased
! MAX_ITER                  : integer           3           !Number of iterations above which the DT is decreased
! LIMIT_ITER                : integer           50          !Number of iterations of a time step (for restart)
! THETA_TOLERANCE           : real              0.001       !Converge Parameter
! INCREASE_DT               : real              1.25        !Increase of DT when iter < MIN_ITER
! DECREASE_DT               : real              0.70        !Decrease of DT when iter > MAX_ITER
!
!
!
!<beginproperty>
! NAME                      : Theta / waterlevel 
!
! see Module FillMatrix for more options
!
!<endproperty>

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


Module ModulePorousMedia

    use ModuleGlobalData
    use ModuleStopWatch
    use ModuleFunctions
    use ModuleTime
    use ModuleHDF5
    use ModuleEnterData
    use ModuleProfile,          only : StartProfile, WriteProfile, KillProfile
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData,    &
                                       KillGridData           
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,    &
                                       GetNumberOfTimeSeries, GetTimeSerieLocation,      &
                                       TryIgnoreTimeSerie, CorrectsCellsTimeSerie,       &
                                       GetTimeSerieName, StartTimeSerieInput,            &
                                       GetTimeSerieInitialData, GetTimeSerieValue
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, GetGridCellArea,               &
                                       WriteHorizontalGrid, UnGetHorizontalGrid,         &
                                       GetXYCellZ, GetCoordTypeList, GetGridCoordType,   &
                                       GetGridLatitudeLongitude, GetCellZInterceptByLine, &
                                       GetCellZInterceptByPolygon
    use ModuleBasinGeometry,    only : GetBasinPoints, GetRiverPoints, GetCellSlope,     &
                                       UnGetBasin
    use ModuleGeometry,         only : ConstructGeometry, GetGeometrySize,               &
                                       GetGeometryDistances, GetGeometryKFloor,          &
                                       UnGetGeometry, ComputeInitialGeometry,            &
                                       ComputeVerticalGeometry, GetGeometryVolumes,      &
                                       GetGeometryAreas, KillGeometry 
    use ModuleMap,              only : ConstructMap, GetWaterPoints3D, GetOpenPoints3D,  &
                                       GetComputeFaces3D, UnGetMap,                      &
                                       UpdateComputeFaces3D, KillMap  
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes, UngetBoxDif, &
                                       BoxDif, KillBoxDif                                                      
    use ModuleFillMatrix,       only : ConstructFillMatrix, ModifyFillMatrix,            &
                                       KillFillMatrix, GetIfMatrixRemainsConstant
    use ModuleDrainageNetwork,  only : GetChannelsWaterLevel, GetChannelsBottomLevel,    &
                                       GetChannelsBottomWidth, GetChannelsOpenProcess,   &
                                       GetChannelsNodeLength, UnGetDrainageNetwork
    use ModuleTriangulation,    only : InterPolation, ConstructTriangulation,            &
                                       SetHeightValues, GetNumberOfTriangles,            &
                                       GetTriangleList, GetNumberOfNodes, GetNodesList,  &
                                       KillTriangulation
    use ModuleDischarges        ,only : Construct_Discharges, GetDischargesNumber,       &
                                        GetDischargesGridLocalization,                   &
                                        GetDischargeWaterFlow, GetDischargesIDName,      &
                                        TryIgnoreDischarge, GetDischargeSpatialEmission, &
                                        CorrectsCellsDischarges, Kill_Discharges,        &
                                        SetLocationCellsZ, SetLayer, GetDischargeON,     &
                                        GetDischargeFlowDistribuiton, UnGetDischarges
    use ModuleDrawing
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  ::  ConstructPorousMedia
    private ::      AllocateInstance
    private ::      ReadDataFile
    private ::      ConstructBottomTopography
    private ::      AllocateVariables
    private ::      ReadSoilTypes
    private ::      InitialFields
    private ::      ReadInitialFile
    private ::      ConstructHDF5Output    
    private ::      ConstructTimeSerie

    !Selector
    public  ::  GetNextPorousMediaDT
    public  ::  GetPotentialInfiltration
    public  ::  GetInfiltration
    public  ::  GetEfectiveEVTP
    public  ::  GetGWFlowToChannels
    public  ::  GetGWFlowToChannelsByLayer
    public  ::  GetGWFlowOption
    public  ::  GetGWToChannelsLayers
    public  ::  GetGWLayer
    public  ::  GetGWLayerOld
    public  ::  GetFlowDischarge
    public  ::  GetPMTotalDischargeFlowVolume
    public  ::  GetPorousMediaTotalStoredVolume
    public  ::  GetFluxU
    public  ::  GetFluxV
    public  ::  GetFluxW
    public  ::  GetUnsatU
    public  ::  GetUnsatV
    public  ::  GetUnsatW
    public  ::  GetUnsatWFinal
!    public  ::  GetWaterColumn
    public  ::  GetWaterContent
    public  ::  GetHead
    public  ::  GetThetaR
    public  ::  GetThetaS
    public  ::  GetThetaF
    public  ::  GetOldWaterContent
    public  ::  GetThetaField
    public  ::  GetComputeSoilField
    public  ::  GetLimitThetaLow
    public  ::  GetUnsatK
    public  ::  GetEvaporation
    public  ::  GetTranspiration
    public  ::  GetEVTPVolumes
    public  ::  GetPMBoundaryFlowVolume
    public  ::  GetPMStoredVolume
    public  ::  GetIgnoreWaterColumnOnEVAP
    public  ::  GetBoundaryImposed
    public  ::  GetBoundaryFluxWalls
    public  ::  GetBoundaryFluxBottom
    public  ::  GetBoundaryCells
    public  ::  UnGetPorousMedia
    
    !Modifier
    public  ::  ModifyPorousMedia      
    private ::      VariableSaturatedFlow
    private ::          Condutivity_Face
    private ::              CondutivityAverage
    private ::              CondutivityMaximum
    private ::              CondutivityMinimum
    private ::              CondutivityGeometricAverage
    private ::          SoilWaterVelocity
    private ::          SoilParameters
    private ::          VerticalContinuity
    private ::          CheckStability
    private ::          ComputeNextDT
!    private ::      ExchangeWithDrainageNetwork
    private ::      PorousMediaOutput
    private ::          OutPutTimeSeries
    private ::      CalculateTotalStoredVolume

    !Destructor
    public  ::  KillPorousMedia                                                     
    private ::      DeAllocateInstance
    private ::      WriteFinalFile

    !Management
    private ::      Ready
    private ::          LocateObjPorousMedia 
    
    !Interfaces----------------------------------------------------------------

    interface  UnGetPorousMedia
        module procedure UnGetPorousMedia_R4
        module procedure UnGetPorousMedia_R8
        module procedure UnGetPorousMedia_R
        module procedure UnGetPorousMedia_R1
        module procedure UnGetPorousMedia_RI
        module procedure UnGetPorousMedia_AI
        module procedure UnGetPorousMedia_AI2D
    end interface UnGetPorousMedia

    !Parameter-------------------------------------------------------------------

    !Conductivity Face
    integer, parameter :: Average           = 1
    integer, parameter :: Maximum           = 2
    integer, parameter :: Minimum           = 3
    integer, parameter :: Weighted          = 4
    integer, parameter :: GeometricAvg      = 5

    !Evapotranspiration Method
    integer, parameter :: SingleEvapotranspiration   = 1
    integer, parameter :: SeparateEvapotranspiration = 2
    !Min Thickness UG Watercolumn
    real, parameter    :: MinUGThickness    = 0.10
    
    !Boundary Conditions interpolation methods
    integer, parameter                          :: Triangulation_   = 1
    integer, parameter                          :: InverseWeight_   = 2
    integer, parameter                          :: NoInterpolation_ = 3
    
    !Boundary Piezometers
    character(len=9 ),           parameter      :: FromTimeSerie        = 'TIMESERIE' 
    character(len=12),           parameter      :: SingleValue          = 'CONSTANT'
    
    
    !Water Contents
    character(LEN = StringLength), parameter :: char_Theta        = trim(adjustl('Theta'            ))

    !Conductivity
    character(LEN = StringLength), parameter :: char_SoilID       = trim(adjustl('SoilID'))                        

    !Waterlevel
    character(LEN = StringLength), parameter :: char_waterlevel   = trim(adjustl('waterlevel'       ))
    
    !Infiltration conductivity - for tests only - by default nothing changes
    integer, parameter                       :: SatCond_   = 1
    integer, parameter                       :: UnSatCond_ = 2
    
    !Bottom Boundary Condition
    integer, parameter                       :: NullGradient_ = 1   !Boundary theta is bottom Theta (Null in Theta)
                                                                     !and no hidrostatic pressure (always moving)
                                                                     !so velocity is conductivity (Final Head gradient is 1)
    
    !Drainage Network link formulations - gone to Global data, because of basin use?
    !integer, parameter :: GWFlowToChanByCell_            = 1
    !integer, parameter :: GWFlowToChanByLayer_           = 2
    
    !GW Flow by Cell - Area method
    integer, parameter :: GWFlowAreaWetPerimeter_            = 1
    integer, parameter :: GWFlowAreaWetPerimeterAndAquifer_  = 2
        
    !Types---------------------------------------------------------------------
    type T_OutPut
        type (T_Time), pointer, dimension(:)    :: OutTime              => null()
        type (T_Time), dimension(:), pointer    :: RestartOutTime       => null()
        type (T_Time), dimension(:), pointer    :: SurfaceOutTime       => null()
        integer                                 :: NextOutPut           = null_int
        logical                                 :: Yes                  = .false.
        logical                                 :: TimeSerieON          = .false.
        logical                                 :: ProfileON            = .false.
        logical                                 :: WriteRestartFile     = .false.
        logical                                 :: SurfaceOutput        = .false.
        logical                                 :: RestartOverwrite     = .false.
        logical                                 :: BoxFluxes            = .false.
        integer                                 :: NextRestartOutput    = 1
        integer                                 :: NextSurfaceOutput    = 1
    end type T_OutPut

    type T_IntegrationByHorizon
        character(len=StringLength)                         :: Name             = null_str
        
        integer                                             :: StartLayer,                      &
                                                               EndLayer
        
        real, dimension(:,:), pointer                       :: Old2D            => null(),      &
                                                               Field2D          => null()
            
        real, dimension(:,:,:), pointer                     :: Old3D            => null(),      &
                                                               Field3D          => null() 
    end type T_IntegrationByHorizon
    
    type T_IntegrationInfo
        logical                                             :: Yes              = .false.,      &
                                                               ByLayer          = .false.,      &
                                                               ByHorizon        = .false.
            
        real, dimension(:,:), pointer                       :: Old2D            => null(),      &
                                                               Field2D          => null()
            
        real, dimension(:,:,:), pointer                     :: Old3D            => null(),      &
                                                               Field3D          => null()  
        
        integer                                             :: HorizonsCount    = 0
        type(T_IntegrationByHorizon), dimension(:), pointer :: Horizons         => null()
    end type T_IntegrationInfo
    
    type T_IntegrationOutput
        type (T_Time), dimension(:), pointer                :: OutTime          => null()
        integer                                             :: NextOutPut       = null_int
        logical                                             :: Yes              = .false.,      &
                                                               Initialize       = .true.
        type (T_IntegrationInfo)                            :: WaterContent,                    &
                                                               RelativeWaterContent,            &
                                                               WaterTable,                      &
                                                               Infiltration,                    &
                                                               BoundaryBottom
        real                                                :: AccTime          = 0.0,          &
                                                               OldAccTime       = 0.0
    end type T_IntegrationOutput    
    
    type T_Files
        character(PathLength)                   :: DataFile             = null_str
        character(PathLength)                   :: InitialFile          = null_str
        character(PathLength)                   :: FinalFile            = null_str
        character(PathLength)                   :: TransientHDF         = null_str
        character(PathLength)                   :: BottomFile           = null_str
        character(PathLength)                   :: ASCFile              = null_str
        character(PathLength)                   :: BoxesFile            = null_str
        character(PathLength)                   :: IntegrationHDFFile    = null_str
        integer                                 :: AsciiUnit            = null_int
    end type T_Files    


    type T_ExtVar
        integer, dimension(:,:), pointer        :: BasinPoints          => null()
        
        !ObjGeometry
        real   , pointer, dimension(:,:  )      :: DUX                  => null()
        real   , pointer, dimension(:,:  )      :: DVY                  => null()
        real   , pointer, dimension(:,:  )      :: DZX                  => null()
        real   , pointer, dimension(:,:  )      :: DZY                  => null()
        real   , pointer, dimension(:,:,:)      :: DZZ                  => null()
        real   , pointer, dimension(:,:,:)      :: DWZ                  => null()
        real   , pointer, dimension(:,:  )      :: YY_IE                => null()
        real   , pointer, dimension(:,:  )      :: XX_IE                => null()
        real   , pointer, dimension(:,:,:)      :: CenterCell           => null()
        real   , pointer, dimension(:,:  )      :: Area                 => null()
        
        real   , pointer, dimension(:,:  )      :: Topography           => null()
        real   , pointer, dimension(:,:  )      :: BottomTopoG          => null()      
        
        real   , pointer, dimension(:,:,:)      :: AreaU                => null()
        real   , pointer, dimension(:,:,:)      :: AreaV                => null()

        integer, dimension(:,:), pointer        :: RiverPoints          => null()
        integer, dimension(:,:), pointer        :: KFloor               => null()

        real(8), pointer, dimension(:,:,:)      :: CellVolume           => null()
        real ,   pointer, dimension(:,:,:)      :: SZZ                  => null()
                
        
        !Map 
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D        => null()
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D         => null()
        integer, pointer, dimension(:,:,:)      :: ComputeFacesU3D      => null()
        integer, pointer, dimension(:,:,:)      :: ComputeFacesV3D      => null()
        integer, pointer, dimension(:,:,:)      :: ComputeFacesW3D      => null() 
       
        real(8), dimension(:,:  ), pointer      :: InfiltrationColumn  => null()
        real, dimension(:,:,:), pointer         :: TranspirationFlux   => null()     
        real, dimension(:,:  ), pointer         :: PotentialEvaporationFlux  => null()       
        logical                                 :: ConstructEvaporation      = .false.
        logical                                 :: ConstructTranspiration    = .false.

        !Time
        type (T_Time)                           :: Now
        real                                    :: DT                   = null_real
    end type T_ExtVar

    !type     T_PointF
    !    real                                    :: X            = null_real
    !    real                                    :: Y            = null_real
    !end type T_PointF

    type     T_FromTimeSerie
        integer                                 :: ObjTimeSerie         = 0
        character(len=StringLength)             :: FileName             = null_str
        integer                                 :: DataColumn           = null_int
    end type T_FromTimeSerie

    type     T_Piezometer
        character(len=StringLength)             :: Name                 = null_str
        type(T_PointF)                          :: Location             
        integer                                 :: ID                   = null_int
        real                                    :: DefaultValue         = null_real
        character(len=StringLength)             :: ValueType            = null_str
        logical                                 :: TimeSerieHasData     = .false.
        type(T_FromTimeSerie)                   :: TimeSerie
        type(T_Piezometer), pointer             :: Next
    end type T_Piezometer

    type T_InverseWeight
        real                                    :: MaxDistance          = null_real
        real                                    :: IWDn                 = null_real
    end type T_InverseWeight

    type     T_Nodes
        real, dimension(:), pointer             :: X
        real, dimension(:), pointer             :: Y
        real, dimension(:), pointer             :: Z
    end type T_Nodes
        
    type T_Triangulation
        logical                                 :: FillOutsidePoints    = .false.
        logical                                 :: OutputTriangles      = .false.
        integer                                 :: Instant              = null_int
        character(len=StringLength)             :: FileName             = null_str
        type (T_Nodes)                          :: Nodes
    end type T_Triangulation
    
    type T_Boundary
        logical                                 :: ImposedLevelConstant = .false.
        logical                                 :: ImposedLevelInTime   = .false. 
        integer                                 :: InterpolationMethod  = null_int  
        integer                                 :: NumberOfPiezometers  = null_int
        integer                                 :: ImposeBoundaryBottomCondition = null_int
        real                                    :: BoundaryValue        = null_real
        real                                    :: MaxDtmForBoundary    = null_real  !Max DTM to apply Walls Boundary
        real                                    :: MinThetaFForBoundary = null_real  !Min ThetaF to apply Bottom Boundary
        type(T_PointF), dimension(:,:), pointer :: GridPoint            => null()
        real, dimension(:,:), pointer           :: ImposedBoundaryLevel => null()
        integer, dimension(:,:), pointer        :: BoundaryCells        => null()
        type (T_Triangulation)                  :: Triangulation
        type (T_InverseWeight)                  :: InverseWeight
        type (T_Piezometer), pointer            :: FirstPiezometer
    end type T_Boundary
    
    !Unsaturated Zone Types
    type T_SoilOptions
        logical :: CalcHorizontal                                       = .false.
        logical :: CalcDrainageNetworkFlux                              = .false.
        integer :: CondutivityFace                                      = null_int
        logical :: Continuous                                           = .false.
        logical :: StopOnWrongDate                                      = .false.
        logical :: CheckGlobalMass                                      = .false.
        logical :: StartWithFieldCapacity                               = .false.
        logical :: ComputeSoilField                                     = .false.
        real    :: HCondFactor                                          = null_real
        real    :: FCHCondFactor                                        = null_real
        logical :: LimitEVAPWaterVelocity                               = .false.
        logical :: LimitEVAPHead                                        = .false.
        logical :: IgnoreWaterColumnOnEvap                              = .false.
        real    :: HeadLimit                                            = null_real
        integer :: DNLink                                               = null_int
        integer :: DNLinkAreaMethod                                     = null_int
        logical :: ComputeHydroPressure                                 = .false.
        integer :: InfiltrationConductivity                             = null_int
        logical :: DryChannelsCompletely                                = .false.                   
        logical :: ImposeBoundary                                       = .false. !Boundary imposed
        logical :: ImposeBoundaryWalls                                  = .false. !Boundary imposed in wall
        logical :: ImposeBoundaryBottom                                 = .false. !Boundary imposed in bottom
        logical :: WriteLog                                             = .false.
        logical :: Discharges                                           = .false.
    end type T_SoilOptions

    type T_SoilType
        real                                :: ThetaR           = null_real
        real                                :: ThetaS           = null_real
        real                                :: nfit             = null_real
        real                                :: mfit             = null_real
        real                                :: alfa             = null_real
        real                                :: lfit             = null_real
        real                                :: SatK             = null_real
        real                                :: OverSatSlope     = null_real
    end type T_SoilType

    type T_Retention !Main parameters in the Mualem-van Genuchten retention and conductivity cuves
        real, dimension(:,:,:), allocatable:: ThetaR          !Minimum water content
        real, dimension(:,:,:), allocatable:: ThetaS          !Saturated water content
        real, dimension(:,:,:), allocatable:: ThetaF          !(Theta-ThetaR)/(ThetaS-ThetaR)
    end type T_Retention

    type T_Converge
        integer                                     :: MinIterations           = 1               
        integer                                     :: MaxIterations           = 1024
        logical                                     :: Stabilize               = .true.
        real                                        :: StabilizeFactor         = 0.01        
        real                                        :: DTFactorUp              = 1.25
        real                                        :: DTFactorDown            = 1.25
        real                                        :: StabilizeHardCutLimit   = 128
        real                                        :: DTSplitFactor           = 2.0               
        real                                        :: CurrentDT               = null_real  
        real                                        :: NextDT                  = null_real
        integer                                     :: LastGoodNiteration      = 1
        integer                                     :: NextNiteration          = 1               
        real                                        :: LimitThetaLo            = 1.0e-15
        real                                        :: LimitThetaHi            = 1.0e-15
        real                                        :: LimitThetaHiGWTable     = 0.9995
        real                                        :: ThetaHydroCoef          = 0.98
        real                                        :: VelHydroCoef            = 1.0
        logical                                     :: CheckDecreaseOnly       = .false.
        real                                        :: MinimumValueToStabilize = 0.1  
        integer                                     :: MinToRestart            = 0
        
        real,    allocatable, dimension(:,:,:)      :: ThetaOld   
        real,    allocatable, dimension(:,:,:)      :: ThetaIni   
        real,    allocatable, dimension(:,:,:)      :: HeadIni    
    end type T_Converge
    
    type       T_PorousMedia        
        !Instaces of other Objects
        integer                                 :: InstanceID
        character(len=StringLength)             :: ModelName
        integer                                 :: ObjBasinGeometry         = 0
        integer                                 :: ObjTime                  = 0
        integer                                 :: ObjGeometry              = 0
        integer                                 :: ObjMap                   = 0
        integer                                 :: ObjHorizontalGrid        = 0
        integer                                 :: ObjHorizontalMap         = 0
        integer                                 :: ObjTopography            = 0
        integer                                 :: ObjTimeSerie             = 0
        integer                                 :: ObjHDF5                  = 0
        integer                                 :: ObjIntegrationHDF5       = 0
        integer                                 :: ObjDrainageNetwork       = 0
        integer                                 :: ObjBottomTopography      = 0
        integer                                 :: ObjEnterData             = 0
        integer                                 :: ObjProfile               = 0
        integer                                 :: ObjTriangulation         = 0
        integer                                 :: ObjDischarges            = 0
        integer                                 :: ObjBoxDif                = 0        
        type (T_PropertyID)                     :: ImpermeableFractionID
        
        real,    allocatable, dimension(:,:,:)  :: ThetaField        !!FieldCapacity [m3/m3]                

        type (T_OutPut)                         :: OutPut
        type (T_IntegrationOutput)              :: IntegrationOutput
        type (T_ExtVar)                         :: ExtVar
        type (T_Files)                          :: Files
        type (T_Time)                           :: BeginTime
        type (T_Time)                           :: EndTime
        type (T_Boundary)                       :: Boundary
        real(8),    pointer, dimension(:,:  )   :: WaterColumn              => null()
        real(8),    pointer, dimension(:,:  )   :: Infiltration             => null()
        real(8),    pointer, dimension(:,:  )   :: EfectiveEVTP             => null()
        real(8),    pointer, dimension(:,:  )   :: EfectiveEVAP             => null()
        
        !For Basin Water Balance
        real(8)                                 :: AccEvapFromSoil       = 0.0 !m3
        real(8)                                 :: AccTranspiration      = 0.0 !m3   
        real(8)                                 :: AccBoundaryFlowVolume = 0.0 !m3 
        real(8)                                 :: TotalDischargeFlowVolume

        !Watertable Properties
        real,    dimension(:,:), pointer        :: OldUGWaterLevel2D        => null()
        real,    dimension(:,:), pointer        :: UGWaterLevel2D           => null()
        real,    dimension(:,:), pointer        :: UGWaterDepth2D           => null()
        integer, dimension(:,:), pointer        :: UGCell                   => null()
        integer, dimension(:,:), pointer        :: UGCell_Old               => null()
        
        real,    dimension(:,:,:), allocatable  :: iFlowBoundaryWalls
        real,    dimension(:,:  ), allocatable  :: iFlowBoundaryBottom         
        !Exchange with channels
        real,    dimension(:,:),   allocatable  :: lFlowToChannels         
        real,    dimension(:,:,:), allocatable  :: lFlowToChannelsLayer     
        real,    dimension(:,:),   allocatable  :: iFlowToChannels          
        real,    dimension(:,:,:), allocatable  :: iFlowToChannelsLayer     
        integer, dimension(:,:),   allocatable  :: FlowToChannelsTopLayer   
        integer, dimension(:,:),   allocatable  :: FlowToChannelsBottomLayer 
        
        !Discharge
        real,    dimension(:,:,:), allocatable  :: lFlowDischarge          !Instantaneous Flow of discharges
        real,    dimension(:,:,:), allocatable  :: iFlowDischarge          !Integrated    Flow of discharges

        !Velocities
        real,    dimension(:,:,:), allocatable  :: UnsatVelU               
        real,    dimension(:,:,:), allocatable  :: UnsatVelV                
        real,    dimension(:,:,:), allocatable  :: UnsatVelW               
        real,    dimension(:,:,:), allocatable  :: UnsatVelWFinal          

        !infiltration 
        real,   pointer, dimension(:,:)         :: InfiltrationVelocity     => null()
        real,   pointer, dimension(:,:)         :: ImpermeableFraction      => null()

        !Fluxes
        real(8), dimension(:,:,:), allocatable  :: FluxU                   
        real(8), dimension(:,:,:), allocatable  :: FluxV                   
        real(8), dimension(:,:,:), allocatable  :: FluxW                   
        real(8), dimension(:,:,:), allocatable  :: FluxWFinal            !Flux Corrected with Vertical continuity
        real,    dimension(:,:  ), pointer      :: EvaporationFlux          => null()
        !Flow Properties
        real,    allocatable, dimension(:,:,:)      :: Theta                !water content on each cell [m3/m3]
        real,    allocatable, dimension(:,:,:)      :: Head                 !Suction Head on each cell 
        real,    allocatable, dimension(:,:,:)      :: HydroPressure        !Hydrostatic pressure
        real,    allocatable, dimension(:,:,:)      :: FinalHead            !Sum of Suction, Hydrostatic and Topography

        !Common Properties
        real,    allocatable, dimension(:,:,:)      :: SatK                
        integer, allocatable, dimension(:,:,:)      :: SoilID              
        real,    allocatable, dimension(:,:,:)      :: UnSatK              
        real,    allocatable, dimension(:,:,:)      :: UnSatK_X            
        real,    allocatable, dimension(:,:,:)      :: UnSatK_Y            
        real,    allocatable, dimension(:,:,:)      :: UnSatK_Z            

        !Auxiliar SpeedUp Matrixes          
        logical, allocatable, dimension(:,:,:)      :: CalculateHead       
    
        logical                                     :: TranspirationExists            = .false.
        logical                                     :: EvaporationExists              = .false.

        type (T_Retention       )                   :: RC           !retention curve
        type (T_Converge        )                   :: CV           !Converge data 

        !Unsaturated Options
        type (T_SoilOptions )                       :: SoilOpt

        real(8)                                     :: TotalStoredVolume    = 0.0
        real(8)                                     :: LossToGround         = 0.0

        !Grid size
        type (T_Size3D)                             :: Size,   WorkSize       
        type (T_Size2D)                             :: Size2D
        
        integer                                     :: DomainCellsNumber
        
        !Soil Types
        type (T_SoilType), dimension(:), pointer    :: SoilTypes      => null()
        
        !Properties
        !type (T_Property), pointer              :: FirstProperty    => null()
                
        type(T_PorousMedia), pointer                :: Next           => null()
        
        !Accumulated fluxes for Average flux computation (for Advection diffusion) 
        real(8), dimension(:,:,:), pointer          :: FluxUAcc       => null()
        real(8), dimension(:,:,:), pointer          :: FluxVAcc       => null()
        real(8), dimension(:,:,:), pointer          :: FluxWAcc       => null()
        real(8), dimension(:,:,:), pointer          :: FluxWAccFinal  => null()

        !integer                                 :: ChunkK = 1,  &
        !                                           ChunkJ = 1,  &
        !                                           ChunkI = 1        
    end type  T_PorousMedia

    !Global Module Variables
    type (T_PorousMedia), pointer               :: FirstObjPorousMedia
    type (T_PorousMedia), pointer               :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructPorousMedia(ModelName,                                  &
                                    ObjPorousMediaID,                           &
                                    ComputeTimeID,                              &
                                    HorizontalGridID,                           &
                                    HorizontalMapID,                            &
                                    TopographyID,                               &
                                    BasinGeometryID,                            &
                                    DrainageNetworkID,                          &
                                    CheckGlobalMass,                            &
                                    ConstructEvaporation,                       &
                                    ConstructTranspiration,                     &
                                    GeometryID,                                 &
                                    MapID,                                      &
                                    DischargesID,                               &
                                    STAT)

        !Arguments---------------------------------------------------------------
        character(len=*)                                :: ModelName
        integer                                         :: ObjPorousMediaID 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: TopographyID
        integer                                         :: BasinGeometryID
        integer                                         :: DrainageNetworkID
        logical                                         :: CheckGlobalMass
        logical                                         :: ConstructEvaporation
        logical                                         :: ConstructTranspiration
        integer, intent (OUT)                           :: GeometryID
        integer, intent (OUT)                           :: DischargesID
        integer, intent (OUT)                           :: MapID
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL
        integer                                         :: DummyI
        real                                            :: DummyR

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_


        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mPorousMedia_)) then
            
            nullify (FirstObjPorousMedia)
            call RegisterModule (mPorousMedia_) 
        
        endif

        call Ready(ObjPorousMediaID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ModelName = ModelName

            !Associate External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           ComputeTimeID   )
            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       TopographyID    ) 
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )

            if (DrainageNetworkID /= 0)                                                     &
                Me%ObjDrainageNetwork  = AssociateInstance (MDRAINAGENETWORK_,  DrainageNetworkID )
                
            Me%SoilOpt%CheckGlobalMass       = CheckGlobalMass
            Me%ExtVar%ConstructEvaporation   = ConstructEvaporation
            Me%ExtVar%ConstructTranspiration = ConstructTranspiration

            !Time
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR010'
                    
            call GetComputeTimeLimits   (Me%ObjTime, BeginTime = Me%BeginTime,              &
                                         EndTime = Me%EndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR020'            

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR030'            

           
            !Read data files
            call ReadDataFile
                        
            call ConstructBottomTopography
            
            !Build 3D domain
            call VerticalDiscretization
            
            call ReadConvergenceParameters
            
            !After 3D build send to ModuleBasin IDs to be associated with ModulePorousMediaProperties and 
            !ModuleVegetation
            GeometryID = Me%ObjGeometry
            MapID      = Me%ObjMap

            !Updates Map (OpenPoints) for first output
            call UpdateComputeFaces3D(  Map_ID         = Me%ObjMap,                       &
                                        DummyR         = DummyR,                          &
                                        DummyI         = DummyI,                          &
                                        STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR060'

            call ReadLockExternalVar 
           
            call AllocateVariables
            
            call ReadSoilTypes
            
            !Build Initial fields
            call InitialFields
            
            !if (Me%SoilOpt%ImposeBoundary) call UpdateBoundaryConditions
            !See if boundary piezometers exist
            if (Me%SoilOpt%ImposeBoundary .and. Me%SoilOpt%ImposeBoundaryWalls) then
                call CheckBoundaryCells
                call ReadBoundaryConditions
                !initial values
                if (.not. Me%Boundary%ImposedLevelConstant) call ModifyBoundaryLevel
            endif
            
            if (Me%SoilOpt%Discharges) then
                call ConstructDischarges
            endif
            
            if (Me%SoilOpt%Continuous)  then
                call ReadInitialFile
            endif
            
            !Calculates initial Theta
            if (.not. Me%SoilOpt%Continuous) then
                call StartWithFieldCapacity
            endif    
            
            !vegetation model growth needs field capacity computation
            if (Me%SoilOpt%ComputeSoilField) then
                call ComputeSoilFieldCapacity
            endif   

            !Set initial time steps
            Me%CV%CurrentDT   = Me%ExtVar%DT

            !Calculate initial heads
            call SoilParameters
            
            !Calculates Initial GW Cell
            call CalculateUGWaterLevel
            
            !Check if river bottom is below soil - warning
            if (DrainageNetworkID /= 0) call CheckRiverBelowSoil
            
            if (Me%OutPut%Yes .or. Me%Output%SurfaceOutput) then
                call ConstructHDF5Output
            endif
            
            if (Me%IntegrationOutput%Yes) then
                call ReadIntegrationConfiguration
                call ConstructIntegrationHDF5Output
            endif

            call ConstructTimeSerie

            call ConstructProfileOutput
            
            if (Me%SoilOpt%WriteLog)  call ConstructASCIIOutput

            call StartOutputBoxFluxes

            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR02'

            if (Me%SoilOpt%CheckGlobalMass) then
                call CalculateTotalStoredVolume
            endif

            !First Output
            if (Me%OutPut%Yes .or. Me%OutPut%SurfaceOutput) call PorousMediaOutput                        
            
            call ReadUnLockExternalVar

            !Returns ID
            ObjPorousMediaID          = Me%InstanceID
            DischargesID              = Me%ObjDischarges
            
            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModulePorousMedia - ConstructPorousMedia - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructPorousMedia
 
    !--------------------------------------------------------------------------

    subroutine StartOutputBoxFluxes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer :: STAT_CALL
        integer :: iflag
        logical :: exist, opened
        character(len=StringLength), dimension(:),  pointer :: FluxesOutputList
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList

        !Begin-----------------------------------------------------------------
       
        ! This keyword have two functions if exist fluxes between boxes are compute 
        ! and the value read is the name file where the boxes are defined
        call GetData(Me%Files%BoxesFile,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'BOXFLUXES',                                      &
                     SearchType     = FromFile,                                         &
                     ClientModule   ='ModulePorousMedia',                               &
                     STAT           = STAT_CALL)                                      

        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Subroutine StartOutputBoxFluxes - ModulePorousMedia. ERR02.'

cd6 :   if (iflag .EQ. 1) then

            Me%Output%BoxFluxes = .true.
            
            inquire(FILE = Me%Files%BoxesFile, EXIST = exist)
cd4 :       if (exist) then
                
                inquire(FILE = Me%Files%BoxesFile, OPENED  = opened)
cd5 :           if (opened) then
                    write(*,*    ) 
                    write(*,'(A)') 'BoxFluxesFileName = ', Me%Files%BoxesFile
                    write(*,*    ) 'Already opened.'
                    stop           'Subroutine StartOutputBoxFluxes; ModulePorousMedia. ERR04'    
                end if cd5

                allocate(FluxesOutputList(1), ScalarOutputList(1)) 

                FluxesOutputList = 'soil_water'
                ScalarOutputList = 'soil_water'

                call StartBoxDif(BoxDifID         = Me%ObjBoxDif,                    &
                                 TimeID           = Me%ObjTime,                      &
                                 HorizontalGridID = Me%ObjHorizontalGrid,            &
                                 BoxesFilePath    = Me%Files%BoxesFile,              &
                                 FluxesOutputList = FluxesOutputList,                &
                                 ScalarOutputList = ScalarOutputList,                &
                                 WaterPoints3D    = Me%ExtVar%WaterPoints3D,         &
                                 Size3D           = Me%Size,                         &
                                 WorkSize3D       = Me%WorkSize,                     &
                                 STAT             = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                         &
                    stop 'Subroutine StartOutputBoxFluxes - ModulePorousMedia. ERR15.'

                deallocate(FluxesOutputList, ScalarOutputList) 
                nullify   (FluxesOutputList, ScalarOutputList)
                
            else
                write(*,*) 
                write(*,*)     'Error dont have the file box.'
                write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
                stop           'Subroutine StartOutputBoxFluxes; ModulePorousMedia. ERR03'    
            end if cd4
        else
            Me%Output%BoxFluxes = .false.        
        end if cd6
        
    end subroutine StartOutputBoxFluxes

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
                        Me%Boundary%BoundaryCells(i,j) = BasinPoint
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

    subroutine CheckRiverBelowSoil

        !Local-----------------------------------------------------------------
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel
        integer                                     :: STAT_CALL, i, j
        character (Len = 5)                         :: str_i, str_j
        character (len = StringLength)              :: StrWarning
        !Begin-----------------------------------------------------------------

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRiverBelowSoil - ModulePorousMedia - ERR02'


do1:    do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then
                
                if (ChannelsBottomLevel(i,j) < Me%ExtVar%BottomTopoG(i, j)) then
                    !write (*,*) 
                    !write (*,*) 'Bottom River section is lower than soil profile in cell', i,j
                    write(str_i, '(i4)') i 
                    write(str_j, '(i4)') j 
                    StrWarning =  'Bottom River section is lower than soil profile in cell (i, j): '// &
                                    str_i//','//str_j
                    call SetError(WARNING_, INTERNAL_, StrWarning, OFF) 
                endif
                
            endif
        
        enddo do2
        enddo do1

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRiverBelowSoil - ModulePorousMedia - ERR06'


    end subroutine CheckRiverBelowSoil

    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance        
                                                    
        !Local-----------------------------------------------------------------
        type (T_PorousMedia), pointer       :: NewObjPorousMedia
        type (T_PorousMedia), pointer       :: PreviousObjPorousMedia


        !Allocates new instance
        allocate (NewObjPorousMedia)
        nullify  (NewObjPorousMedia%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjPorousMedia)) then
           
            FirstObjPorousMedia         => NewObjPorousMedia
            Me                          => NewObjPorousMedia
        
        else
        
            PreviousObjPorousMedia      => FirstObjPorousMedia
            Me                          => FirstObjPorousMedia%Next
            do while (associated(Me))
                PreviousObjPorousMedia  => Me
                Me                      => Me%Next
            enddo
            Me                          => NewObjPorousMedia
            PreviousObjPorousMedia%Next => NewObjPorousMedia
        
        endif

        Me%InstanceID = RegisterNewInstance (mPorousMedia_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    
    subroutine ReadDataFile        

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        
        !----------------------------------------------------------------------
        
        !Reads the name of the data file from nomfich
        call ReadFileName ('POROUS_DATA', Me%Files%DataFile, "PorousMedia Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR010'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('POROUS_HDF', Me%Files%TransientHDF, "PorousMedia HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR020'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('POROUS_FIN', Me%Files%FinalFile, "PorousMedia Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR030'

        !Constructs the DataFile
        call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR050'
        
        call GetData(Me%SoilOpt%WriteLog,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'WRITE_LOG',                                          &
                     Default    = .false.,                                              &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR035'        
        
        if (Me%SoilOpt%WriteLog) then
            !Reads the name of the file where to store ASCII data
            call ReadFileName ('POROUS_ASC', Me%Files%ASCFile, "PorousMedia ITER SOL File", Me%EndTime,  &
                               Extension = 'hyt', STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR040'
        endif

        !Botom file
        call GetData(Me%Files%BottomFile,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BOTTOM_FILE',                                      &
                     ClientModule = 'ModulePorousMedia',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR060'

        !General Options        
        call GetData(Me%SoilOpt%StartWithFieldCapacity,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'START_WITH_FIELD',                                   &
                     Default    = .true.,                                               &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR070'

        call GetData(Me%SoilOpt%ComputeSoilField,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'COMPUTE_SOIL_FIELD',                                 &
                     Default    = .false.,                                              &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR080'

        call GetData(Me%SoilOpt%LimitEVAPHead,                                          &     
                     Me%ObjEnterData, iflag,                                            &
                     SearchType     = FromFile,                                         &
                     keyword        ='LIMIT_EVAP_HEAD',                                 &
                     Default        = .false.,                                          &
                     ClientModule   ='ModulePorousMedia',                               &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR100'

        if (Me%SoilOpt%LimitEVAPHead) then
                
            call GetData(Me%SoilOpt%HeadLimit,                                         &     
                         Me%ObjEnterData, iflag,                                       &
                         SearchType     = FromFile,                                    &
                         keyword        ='HEAD_LIMIT',                                 &
                         Default        = -100.0,                                      &
                         ClientModule   ='ModulePorousMedia',                          &
                         STAT           = STAT_CALL)             
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR110'
        
        endif

        call GetData(Me%SoilOpt%Continuous,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'CONTINUOUS',                                         &
                     Default    = .false.,                                              &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR120'
        
        if (Me%SoilOpt%Continuous) then
            !Reads the name of the file where to read initial data
            call ReadFileName ('POROUS_INI', Me%Files%InitialFile, "PorousMedia Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR130'

            call GetData(Me%SoilOpt%StopOnWrongDate,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType = FromFile,                                         &
                         keyword    = 'STOP_ON_WRONG_DATE',                             &
                         Default    = .true.,                                           &                                           
                         ClientModule ='ModulePorousMedia',                             &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR140'

        endif

        call GetData(Me%SoilOpt%IgnoreWaterColumnOnEvap,                                &     
                     Me%ObjEnterData, iflag,                                            &
                     SearchType     = FromFile,                                         &
                     keyword        ='IGNORE_WATER_COLUMN_ON_EVAP',                     &
                     Default        = .true.,                                           &
                     ClientModule   ='ModulePorousMedia',                               &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR150'

        !Output Options--------------------------------------------------------

        !Gets Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%ExtVar%Now,                                 &
                           EndTime     = Me%EndTime,                                    &
                           keyword     = 'OUTPUT_TIME',                                 &
                           SearchType  = FromFile,                                      &
                           OutPutsTime = Me%OutPut%OutTime,                             &
                           OutPutsOn   = Me%OutPut%Yes,                                 &
                           STAT        = STAT_CALL)
        Me%OutPut%NextOutPut = 1
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR160'

        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%ExtVar%Now,                                &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR170'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModulePorousMedia',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModulePorousMedia - ERR180'

        !Output for surface output
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%ExtVar%Now,                                &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'SURFACE_OUTPUT_TIME',                        &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%SurfaceOutTime,                     &
                           OutPutsOn    = Me%OutPut%SurfaceOutput,                      &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR190'

        !IN PROGRESS
		!Sets Integrated Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime = Me%ExtVar%Now,                                 &
                           EndTime     = Me%EndTime,                                    &
                           keyword     = 'INTEGRATION_TIME',						    &
                           SearchType  = FromFile,                                      &
                           OutPutsTime = Me%IntegrationOutput%OutTime,                  &
                           OutPutsOn   = Me%IntegrationOutput%Yes,					    &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR191' 
        Me%IntegrationOutput%NextOutput = 1
        
        if (Me%IntegrationOutput%Yes) then
            call ReadFileName('POROUS_INT_HDF', Me%Files%IntegrationHDFFile,	        &
                              Message = "Porous Media Integration HDF File",		    &
                              STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR192'
            
            
            
        endif        
        
        
        !Checks consistency
        if (Me%OutPut%Yes .and. Me%OutPut%SurfaceOutput) then
            write(*,*)'Only normal output or 2D output can be active'
            write(*,*)'OUTPUT_TIME or SURFACE_OUTPUT_TIME'
            stop 'ReadDataFile - ModulePorousMedia - ERR200'
        endif

        !Directional Options---------------------------------------------------        

        call GetData(Me%SoilOpt%CalcHorizontal,                                 &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CALC_HORIZONTAL',                         &
                     Default        =.TRUE.,                                    &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR210'


        call GetData(Me%SoilOpt%CalcDrainageNetworkFlux,                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CALC_DRAINAGE_FLUX',                      &
                     Default        =.TRUE.,                                    &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR220'
       
        ! 1 - AVERAGE of the conductivity in the cells
        ! 2 - MAXIMUM of the conductivity in the cells
        ! 3 - MINIMUM of the conductivity in the cells
        ! 4 - WEIGTHED of the conductivity in the cells
        ! 5 - GEOMETRIC AVERAGE of the conductivity in the cells

        call GetData(Me%SoilOpt%CondutivityFace,                                &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CONDUTIVITYFACE',                         &
                     Default        = 1    ,                                    &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR230'


        call GetData(Me%SoilOpt%HCondFactor,                                    &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='HORIZONTAL_K_FACTOR',                     &
                     Default        = 1.0,                                      &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR240'

        call GetData(Me%SoilOpt%FCHCondFactor,                                  &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='FC_K_FACTOR',                             &
                     Default        = Me%SoilOpt%HCondFactor,                   &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR250'                
        
        call GetData(Me%SoilOpt%LimitEVAPWaterVelocity,                         &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'LIMIT_EVAP_WATER_VEL',                       &
                     Default    = .false.,                                      &                                           
                     ClientModule ='ModulePorousMedia',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR260'

        call GetData(Me%SoilOpt%DNLink,                                         &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'DN_LINK',                                    &
                     Default    = GWFlowToChanByCell_,                          &                                           
                     ClientModule ='ModulePorousMedia',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR270'
        
        if ((Me%SoilOpt%DNLink /= GWFlowToChanByLayer_) .and. (Me%SoilOpt%DNLink /= GWFlowToChanByCell_)) then
            write(*,*)' DN_LINK uncorrectly defined - 2 for GW flow to drainage network by layers '
            write(*,*)' and 1 for GW flow for each cell'
            stop 'ReadDataFile - ModulePorousMedia - ERR271'
        endif
        
        !New method
        if (Me%SoilOpt%DNLink == GWFlowToChanByCell_) then
            call GetData(Me%SoilOpt%DNLinkAreaMethod,                               &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType = FromFile,                                     &
                         keyword    = 'DN_LINK_AREA_METHOD',                        &
                         Default    = GWFlowAreaWetPerimeter_,                      &                                           
                         ClientModule ='ModulePorousMedia',                         &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR275'            
            if ((Me%SoilOpt%DNLinkAreaMethod /= GWFlowAreaWetPerimeter_) .and.    &
                (Me%SoilOpt%DNLinkAreaMethod /= GWFlowAreaWetPerimeterAndAquifer_)) then
                write(*,*)' DN_LINK_AREA_METHOD uncorrectly defined '
                write(*,*)' 1 - area for Flux is around wet perimeter '
                write(*,*)' 2 - area for Flux is around wet perimeter + aquifer if aquifer higher then river '
                stop 'ReadDataFile - ModulePorousMedia - ERR276'
            endif
        
        endif
        
        call GetData(Me%SoilOpt%ComputeHydroPressure,                           &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'COMPUTE_HYDRO_PRESSURE',                     &
                     Default    = .true.,                                       &                                           
                     ClientModule ='ModulePorousMedia',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR280'
        
        !Impose boundary saturated levels in walls
        call GetData(Me%SoilOpt%ImposeBoundaryWalls,                            &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'IMPOSE_BOUNDARY_VALUE',                      &
                     Default    = .false.,                                      & 
                     ClientModule ='ModulePorousMedia',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR290'
        
        if (Me%SoilOpt%ImposeBoundaryWalls) then
            
            Me%SoilOpt%ImposeBoundary = .true.
            
            call GetData(Me%Boundary%MaxDtmForBoundary,                         &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'MAX_DTM_FOR_BOUNDARY',                 &
                         ClientModule = 'ModulePorousMedia',                    &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR291'        

            if (iflag == 0) then
                write(*,*)'MAX_DTM_FOR_BOUNDARY must be defined in module PorousMedia'
                stop 'ReadDataFile - ModulePorousMedia - ERR0292'
            endif
            
            !Moved to after checking if piezometers exist
!            call GetData(Me%Boundary%BoundaryValue,                                 &
!                         Me%ObjEnterData, iflag,                                    &
!                         SearchType = FromFile,                                     &
!                         keyword    = 'BOUNDARY_VALUE',                             &
!                         ClientModule ='ModulePorousMedia',                         &
!                         STAT       = STAT_CALL)            
!            if (STAT_CALL /= SUCCESS_) stop 'GetUnSaturatedOptions - ModulePorousMedia - ERR271' 
!
!            if (iflag == 0) then
!                write(*,*)'BOUNDARY_VALUE must be defined in module PorousMedia'
!                stop 'ReadDataFile - ModulePorousMedia - ERR0230'
!            endif            
            
        endif

        !Impose boundary in bottom
        call GetData(Me%SoilOpt%ImposeBoundaryBottom,                           &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'IMPOSE_BOUNDARY_BOTTOM',                     &
                     Default    = .false.,                                      & 
                     ClientModule ='ModulePorousMedia',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR293'
        
        if (Me%SoilOpt%ImposeBoundaryBottom) then
            
            Me%SoilOpt%ImposeBoundary = .true.
            
            call GetData(Me%Boundary%MinThetaFForBoundary,                      &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'MIN_THETAF_FOR_BOUNDARY',              &
                         Default      = 0.0,                                    &
                         ClientModule = 'ModulePorousMedia',                    &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR294'        

            call GetData(Me%Boundary%ImposeBoundaryBottomCondition,             &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'IMPOSE_BOUNDARY_BOTTOM_CONDITION',     &
                         Default      = NullGradient_,                          &
                         ClientModule = 'ModulePorousMedia',                    &
                         SearchType   = FromFile,                               &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR295'
            
            if (Me%Boundary%ImposeBoundaryBottomCondition /= NullGradient_) then
                write(*,*)'IMPOSE_BOUNDARY_BOTTOM_CONDITION for now can only be 1 (NullGradient)'
                stop 'ReadDataFile - ModulePorousMedia - ERR0296'
            endif
            
        endif

        !Discharges
        call GetData(Me%SoilOpt%Discharges,                                 &
                     Me%ObjEnterData, iflag,                                &  
                     keyword      = 'DISCHARGES',                           &
                     ClientModule = 'ModulePorousMedia',                    &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR300'

        !Conductivity used for infiltration - for tests only - by default nothing changes
        call GetData(Me%SoilOpt%InfiltrationConductivity,                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType = FromFile,                                     &
                     keyword    = 'INFIL_CONDUCTIVITY',                         &
                     Default    = SatCond_,                                     &                                           
                     ClientModule ='ModulePorousMedia',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR310'
        
        if ((Me%SoilOpt%InfiltrationConductivity .ne. SatCond_) .and.            &
            (Me%SoilOpt%InfiltrationConductivity .ne. UnSatCond_)) then

            write(*,*)
            write(*,*)'Method not known for infiltration conductivity'
            write(*,*)'Please check INFIL_CONDUCTIVITY keyword'
            stop 'ReadDataFile - ModulePorousMedia - ERR320'
        
        endif
    
    end subroutine ReadDataFile
    
    !--------------------------------------------------------------------------
    
    subroutine ReadConvergenceParameters
    
        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL,               &
                                                       iflag,                   &
                                                       MIN_ITER_flag,           &
                                                       MAX_ITER_flag,           &
                                                       LIMIT_ITER_flag,         &
                                                       THETA_TOLERANCE_flag,    &
                                                       INCREASE_DT_flag,        &
                                                       DECREASE_DT_flag,        &
                                                       CUT_OFF_THETA_GW_flag    
                                                            
        real                                        :: dummy_real
        integer                                     :: dummy_int
        !integer                                     :: dummy_int,       &
        !                                               ChunkKFactor,    &
        !                                               ChunkJFactor,    &
        !                                               ChunkIFactor
        
        !----------------------------------------------------------------------    
        
        !----------------------------------------------------------------------
        !Find deprecated keywords in data file
        !----------------------------------------------------------------------
        call GetData(dummy_int,                                                 &
                     Me%ObjEnterData, MIN_ITER_flag,                            &
                     SearchType     = FromFile,                                 &
                     keyword        ='MIN_ITER',                                &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR010")
        call GetData(dummy_int,                                                 &
                     Me%ObjEnterData, MAX_ITER_flag,                            &
                     SearchType     = FromFile,                                 &
                     keyword        ='MAX_ITER',                                &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR020")
        call GetData(dummy_int,                                                 &
                     Me%ObjEnterData, LIMIT_ITER_flag,                          &
                     SearchType     = FromFile,                                 &
                     keyword        ='LIMIT_ITER',                              &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR030")
        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, THETA_TOLERANCE_flag,                     &
                     SearchType     = FromFile,                                 &
                     keyword        ='THETA_TOLERANCE',                         &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR040")
        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, INCREASE_DT_flag,                         &
                     SearchType     = FromFile,                                 &
                     keyword        ='INCREASE_DT',                             &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR050")
        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, DECREASE_DT_flag,                         &
                     SearchType     = FromFile,                                 &
                     keyword        ='DECREASE_DT',                             &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR060")
        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, CUT_OFF_THETA_GW_flag,                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CUT_OFF_THETA_HIGH_GW_TABLE',             &                     
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR061")             

        if (MIN_ITER_flag > 0 .or. MAX_ITER_flag > 0 .or. LIMIT_ITER_flag > 0 .or. &
            THETA_TOLERANCE_flag > 0 .or. INCREASE_DT_flag > 0 .or. DECREASE_DT_flag > 0 .or. &
            CUT_OFF_THETA_GW_flag > 0) then
            
            write (*,*) '======================================================================='
            write (*,*) 'The following deprecated keywords were found in Porous Media data file:'
            write (*,*) ''
            
            if (MIN_ITER_flag > 0) &
                write(*,*) 'MIN_ITER                    : Use MIN_ITERATIONS instead.'
            if (MAX_ITER_flag > 0) &
                write(*,*) 'MAX_ITER'
            if (LIMIT_ITER_flag > 0) &
                write(*,*) 'LIMIT_ITER                  : Use MAX_ITERATIONS instead.'            
            if (THETA_TOLERANCE_flag > 0) &
                write(*,*) 'THETA_TOLERANCE             : Use STABILIZE_FACTOR instead.'
            if (INCREASE_DT_flag > 0) &
                write(*,*) 'INCREASE_DT                 : Use DT_FACTOR_UP instead.'
            if (DECREASE_DT_flag > 0) &
                write(*,*) 'DECREASE_DT                 : Use DT_FACTOR_DOWN instead.'
            if (CUT_OFF_THETA_GW_flag > 0) &
                write(*,*) 'CUT_OFF_THETA_HIGH_GW_TABLE : Use GW_SAT_FACTOR instead.'
                
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR070")                              
        endif

        !----------------------------------------------------------------------
        !Read convergence options
        !---------------------------------------------------------------------- 
        
        !call GetData(ChunkKFactor,                                              &
        !             Me%ObjEnterData, iflag,                                    &  
        !             keyword      = 'CHUNK_K_FACTOR',                           &
        !             ClientModule = 'ModulePorousMedia',                        &
        !             Default      = 3,                                          &
        !             SearchType   = FromFile,                                   &
        !             STAT         = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) & 
        !    call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR071")
        !
        !call GetData(ChunkJFactor,                                              &
        !             Me%ObjEnterData, iflag,                                    &  
        !             keyword      = 'CHUNK_I_FACTOR',                           &
        !             ClientModule = 'ModulePorousMedia',                        &
        !             Default      = 10,                                         &
        !             SearchType   = FromFile,                                   &
        !             STAT         = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) & 
        !    call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR071")
        !    
        !call GetData(ChunkIFactor,                                              &
        !             Me%ObjEnterData, iflag,                                    &  
        !             keyword      = 'CHUNK_J_FACTOR',                           &
        !             ClientModule = 'ModulePorousMedia',                        &
        !             Default      = 10,                                         &
        !             SearchType   = FromFile,                                   &
        !             STAT         = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) & 
        !    call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR071")            

        Me%DomainCellsNumber = CountDomainPoints(ChunkKFactor, ChunkJFactor, ChunkIFactor)
        
        call GetData(Me%CV%Stabilize,                                           &
                     Me%ObjEnterData, iflag,                                    &  
                     keyword      = 'STABILIZE',                                &
                     ClientModule = 'ModulePorousMedia',                        &
                     SearchType   = FromFile,                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) & 
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR080")        
        if (iflag <= 0) then
            write(*,*) 'WARNING: Missing STABILIZE keyword in Porous Media input data file.'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR081")
        endif         
        if (Me%CV%Stabilize) then                
            !Maximun change of water content (in %) allowed in one time step.
            call GetData(Me%CV%StabilizeFactor,                                     &
                         Me%ObjEnterData, iflag,                                    &  
                         keyword      = 'STABILIZE_FACTOR',                         &
                         ClientModule = 'ModulePorousMedia',                        &
                         SearchType   = FromFile,                                   &
                         Default      = 0.1,                                        &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) & 
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR082")

            if (Me%CV%StabilizeFactor < 0.0 .or. Me%CV%StabilizeFactor > 1.0) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR083")
                
            call GetData(Me%CV%MinimumValueToStabilize,                     &
                         Me%ObjEnterData, iflag,                            &
                         SearchType   = FromFile,                           &
                         keyword      = 'STABILIZE_MIN_FACTOR',             &
                         default      = 0.0,                                &
                         ClientModule = 'ModulePorousMedia',                &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR084")
            if (Me%CV%MinimumValueToStabilize < 0.0) then
                write (*,*)'Invalid Minimun to Stabilize value [STABILIZE_MIN_FACTOR]'
                write (*,*)'Value must be >= 0'            
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR085")
            endif 
            
            call GetData(dummy_real,                                            &
                         Me%ObjEnterData, iflag,                                &  
                         keyword      = 'STABILIZE_RESTART_FACTOR',             &
                         ClientModule = 'ModuleRunOff',                         &
                         SearchType   = FromFile,                               &
                         Default      = 0.,                                     &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR086")
            Me%CV%MinToRestart = max(int(Me%DomainCellsNumber * dummy_real), 0)
            
            call GetData(Me%CV%CheckDecreaseOnly,                                   &
                         Me%ObjEnterData, iflag,                                    &  
                         keyword      = 'CHECK_DEC_ONLY',                           &
                         ClientModule = 'ModulePorousMedia',                        &
                         SearchType   = FromFile,                                   &
                         Default      = .false.,                                    &
                         STAT         = STAT_CALL)                                  
            if (STAT_CALL /= SUCCESS_) &
                call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR087")
        endif        

       !Number of iterations threshold for starting to ask for a lower DT 
        call GetData(Me%CV%MinIterations,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MIN_ITERATIONS',                          &
                     Default        = 1,                                        &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR090")
        if (Me%CV%MinIterations < 1) then
            write (*,*)'Invalid Minimun Iterations value [MIN_ITERATIONS]'
            write (*,*)'Value must be greater than 0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR091")
        endif                                 

        !Number of iterations threshold that causes the model to stop
        call GetData(Me%CV%MaxIterations,                                       &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='MAX_ITERATIONS',                          &
                     Default        = 1024,                                     &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR100")
        if (Me%CV%MaxIterations < Me%CV%MinIterations) then
            write (*,*)'Invalid Maximun Iterations value [MAX_ITERATIONS]'
            write (*,*)'Value must be greater than the value of MIN_ITERATIONS'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR101")              
        endif
                            
        !% of the maximun iterations that causes the DT to be cut to the value of one internal time step
        call GetData(dummy_real,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'DT_CUT_FACTOR',                    &
                     default      = 0.1,                                &
                     ClientModule = 'ModulePorousMedia',                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR110") 
        if (dummy_real <= 0.0 .or. dummy_real > 1.0) then
            write (*,*)'Invalid DT Cut Factor [DT_CUT_FACTOR]'
            write (*,*)'Value must be >= 0.0 and < 1.0'        
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR111") 
        endif
        Me%CV%StabilizeHardCutLimit = dummy_real * Me%CV%MaxIterations
        
       !Internal Time Step Split
        call GetData(Me%CV%DTSplitFactor,                                   &
                     Me%ObjEnterData, iflag,                                &
                     keyword      = 'DT_SPLIT_FACTOR',                      &
                     ClientModule = 'ModulePorousMedia',                    &
                     SearchType   = FromFile,                               &
                     Default      = 2.0,                                    &
                     STAT         = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'ReadConvergenceParameters - ModulePorousMedia - ERR120'        
        if (Me%CV%DTSplitFactor <= 1.0) then
            write (*,*)'Invalid DT Split Factor [DT_SPLIT_FACTOR]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR121")              
        endif            

        call GetData(dummy_real,                                                &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR',                               &
                     Default        = 1.25,                                     &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR130")             
        if (dummy_real <= 1.0) then
            write (*,*)'Invalid DT Factor [DT_FACTOR]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR131")              
        endif            
        
        call GetData(Me%CV%DTFactorUp,                                          &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR_UP',                            &
                     Default        = dummy_real,                               &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR140")  
        if (Me%CV%DTFactorUp <= 1.0) then
            write (*,*)'Invalid DT Factor Up [DT_FACTOR_UP]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR141")              
        endif                  
                
        call GetData(Me%CV%DTFactorDown,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DT_FACTOR_DOWN',                          &
                     Default        = dummy_real,                               &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR150")  
        if (Me%CV%DTFactorDown <= 1.0) then
            write (*,*)'Invalid DT Factor Down [DT_FACTOR_DOWN]'
            write (*,*)'Value must be greater then 1.0'
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR151")
        endif                     
                
        ! This value says how near the ThetaR the calculation is disconected. 
        ! Disables calculation when Theta is near ThetaR
        call GetData(Me%CV%LimitThetaLo,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CUT_OFF_THETA_LOW',                       &
                     Default        = 1.0e-15,                                  &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR160") 
      
        ! This value says when Theta is converted to ThetaS
        ! Set Theta = ThetaS
        call GetData(Me%CV%LimitThetaHi,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CUT_OFF_THETA_HIGH',                      &
                     Default        = 1.0e-15,                                  &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR170") 
        
        !avoid instabilities in searching for saturated water table. using low values as 1e-15 
        !usually causes variations in water table depth of order of meters from iteration to another
        !because of small theta variations (speially problematic in river cells)
        call GetData(Me%CV%LimitThetaHiGWTable,                                 &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='GW_SAT_FACTOR',                           &
                     Default        = 0.99,                                     &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR180") 
            
        ! This value says from which thetaS hydrostatic pressure is to be consider
        ! Set Theta = ThetaS
        call GetData(Me%CV%ThetaHydroCoef,                                      &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='THETA_HYDRO_COEF',                        &
                     Default        = 0.98,                                     &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR190") 

        call GetData(Me%CV%VelHydroCoef,                                        &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='VEL_HYDRO_COEF',                          &
                     Default        = 1.00,                                     &
                     ClientModule   ='ModulePorousMedia',                       &
                     STAT           = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, KEYWORD_, "ReadConvergenceParameters - ModulePorousMedia - ERR200")
        
        !----------------------------------------------------------------------
    
    end subroutine ReadConvergenceParameters
    
    !--------------------------------------------------------------------------
    
    integer function CountDomainPoints (ChunkKFactor, ChunkJFactor, ChunkIFactor)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: ChunkKFactor,    &
                                                       ChunkJFactor,    &
                                                       ChunkIFactor        
        
        !Local----------------------------------------------------------------- 
        integer                                     :: i, j, k
        integer                                     :: count_i, count_j, count_k, STAT_CALL
        logical                                     :: j_has, i_has
        
        !Begin-----------------------------------------------------------------       
                
        CountDomainPoints = 0
        count_i = 0
        count_j = 0
        count_k = 0      
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CountDomainPoints - ModulePorousMedia - ERR010'    
        
        call GetGeometryKFloor(Me%ObjGeometry,                                          &
                               Z    = Me%ExtVar%KFloor,                                 &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "CountDomainPoints; ModulePorousMedia. ERR011")              
        
        !Initializes Water Column
        do j = Me%Size%JLB, Me%Size%JUB
        
            j_has = .false.
            
            do i = Me%Size%ILB, Me%Size%IUB

                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                    j_has = .true.
                    i_has = .false.                    
                
                    do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB - 1
                        i_has = .true.            
                        count_k = count_k + 1
                    enddo 
                    
                    if (i_has) count_i = count_i + 1
                    
                endif

            enddo
        
            if (j_has) count_j = count_j + 1
        
        enddo
        
        ChunkK = max((Me%Size%KUB - Me%Size%KLB) / ChunkKFactor, 1)
        !Me%ChunkJ = max((Me%Size%JUB - Me%Size%JLB) / ChunkJFactor, 1)
        !Me%ChunkI = max((Me%Size%IUB - Me%Size%ILB) / ChunkIFactor, 1)
        
        !write (*,*) 'CHUNK factors:', ChunkKFactor, ChunkJFactor, ChunkIFactor, Me%ChunkK, Me%ChunkJ, Me%ChunkI
        !stop
        
        CountDomainPoints = count_k
        
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CountDomainPoints - ModulePorousMedia - ERR020'       
        
        call UnGetGeometry(Me%ObjGeometry, Me%ExtVar%KFloor, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "CountDomainPoints; ModulePorousMedia. ERR021")         
           
    end function CountDomainPoints
    
    !-------------------------------------------------------------------------
    
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
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModulePorousMedia - ERR01'

        open(Unit = InitialFile, File = Me%Files%InitialFile, Form = 'UNFORMATTED', status = 'OLD', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModulePorousMedia - ERR02'

        !Reads Date
        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModulePorousMedia - ERR03'
        
        DT_error = EndTimeFile - BeginTime

        !Avoid rounding erros
        if (abs(DT_error) >= 0.01) then
            
            write(*,*) 'The end time of the previous run is different from the start time of this run'
            write(*,*) 'Date in the file'
            write(*,*) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
            write(*,*) 'DT_error', DT_error
            if (Me%SoilOpt%StopOnWrongDate) stop 'ReadInitialFileOld - ModulePorousMedia - ERR04'   

        endif

        read(InitialFile)Me%Theta

        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFileOld - ModulePorousMedia - ERR05'
        

    end subroutine ReadInitialFile

    !--------------------------------------------------------------------------

    subroutine ConstructBottomTopography        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Constructs GridData
        call ConstructGridData      (Me%ObjBottomTopography, Me%ObjHorizontalGrid,  &
                                     FileName = Me%Files%BottomFile,                &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR05'


    end subroutine ConstructBottomTopography

    !--------------------------------------------------------------------------
    
    subroutine VerticalDiscretization
                
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !get the topography
        call GetGridData      (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR00'

         !Constructs GridData
        call ConstructGeometry (GeometryID       = Me%ObjGeometry,                 &
                                GridDataID       = Me%ObjBottomTopography,         &
                                HorizontalGridID = Me%ObjHorizontalGrid,           &
                                HorizontalMapID  = Me%ObjHorizontalMap,            &
                                ActualTime       = Me%ExtVar%Now,                  &
                                SurfaceElevation = Me%ExtVar%Topography,           &
                                BathymTopoFactor = -1.0,                           &
                                STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - VerticalDiscretization - ERR01'
        
        !Map - Soil Column            
        call ConstructMap       (Map_ID           = Me%ObjMap,                      &
                                 GeometryID       = Me%ObjGeometry,                 &
                                 HorizontalMapID  = Me%ObjHorizontalMap,            &
                                 TimeID           = Me%ObjTime,                     &
                                 STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - VerticalDiscretization - ERR02'
        
        !Get water points
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - VerticalDiscretization - ERR03'


        !Geometry Size
        call GetGeometrySize    (Me%ObjGeometry,             &    
                                 Size     = Me%Size,         &
                                 WorkSize = Me%WorkSize,     &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR05'

        Me%Size2D%ILB = Me%Size%ILB
        Me%Size2D%IUB = Me%Size%IUB
        Me%Size2D%JLB = Me%Size%JLB
        Me%Size2D%JUB = Me%Size%JUB      

        !Initial geometry
        call ComputeInitialGeometry(Me%ObjGeometry,                                 &
                                    WaterPoints3D       = Me%ExtVar%WaterPoints3D,  &
                                    SurfaceElevation    = Me%ExtVar%Topography,     &
                                    ContinuesCompute    = .false.,                  &
                                    ActualTime          = Me%ExtVar%Now,            &
                                    STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR06'
        
       
        !Checks Vertical Discretization
        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  DWZ         = Me%ExtVar%DWZ,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia. ERR07'

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%DWZ, STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia. ERR09'

        call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "VerticalDiscretization; ModulePorousMedia. ERR10") 

        call UngetGridData (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR11'

    end subroutine VerticalDiscretization

    !--------------------------------------------------------------------------

    subroutine ConstructAsciiOutPut            

        !Local-----------------------------------------------------------------
        integer               :: status
        integer               :: STAT_CALL                
        integer               :: Counter
        character(LEN=4)      :: Number

        call UnitsManager(Me%Files%AsciiUnit, OPEN_FILE, STAT = status) 
        if (status /= SUCCESS_) stop "ConstructAsciiOutPut - ModulePorousMedia - ERR01"

        Counter  = 1
do1:     do
            Number = '    '
            write(Number, fmt='(i4)')Counter
            open(UNIT   = Me%Files%AsciiUnit,                                      &
!                 FILE   = '..\res\iter.soi_'//trim(adjustl(Number))//'.log', &
                 FILE   = Me%Files%ASCFile, &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do1
            else
                Counter = Counter + 1
            end if
        enddo do1

        write (Me%Files%AsciiUnit, FMT=*) 'YY     MM   DD   HH   MM     SS       Iter  Time_Step '
        write (Me%Files%AsciiUnit, FMT=*) '                                                 s    '
    
    end subroutine ConstructAsciiOutPut
    
    !--------------------------------------------------------------------------

    subroutine AllocateVariables        
        
        !Local-----------------------------------------------------------------        
        integer                                         :: ILB, IUB, JLB,  JUB 
        integer                                         :: KLB, KUB 

        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KUB = Me%Size%KUB
        KLB = Me%Size%KLB
        
               
        !Water Content---------------------------------------------------------
        allocate (Me%Theta          (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%Head           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%HydroPressure  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%FinalHead      (ILB:IUB,JLB:JUB,KLB:KUB))
        
        allocate(Me%UGWaterLevel2D  (ILB:IUB,JLB:JUB)) 
        allocate(Me%UGWaterDepth2D  (ILB:IUB,JLB:JUB))
        allocate(Me%UGCell          (ILB:IUB,JLB:JUB)) 
        allocate(Me%UGCell_Old      (ILB:IUB,JLB:JUB)) 
        allocate(Me%WaterColumn     (ILB:IUB,JLB:JUB))
        allocate(Me%Infiltration    (ILB:IUB,JLB:JUB))
        allocate(Me%EfectiveEVTP    (ILB:IUB,JLB:JUB))
        allocate(Me%ImpermeableFraction (ILB:IUB,JLB:JUB))

        if (Me%SoilOpt%ComputeSoilField) then
            allocate (Me%ThetaField          (ILB:IUB,JLB:JUB,KLB:KUB))
        endif
          
        Me%Theta                = null_real
        Me%Head                 = null_real
        Me%HydroPressure        = null_real
        Me%FinalHead            = null_real

        Me%UGWaterLevel2D       = null_real
        Me%UGWaterDepth2D       = null_real
        Me%UGCell               = null_int
        Me%UGCell_Old           = null_int
        Me%WaterColumn          = null_real
        Me%Infiltration         = null_real
        Me%EfectiveEVTP         = null_real
        Me%ImpermeableFraction  = null_real

        !Conductivities--------------------------------------------------------
        allocate (Me%SatK               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK             (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK_X           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK_Y           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%UnsatK_Z           (ILB:IUB,JLB:JUB,KLB:KUB))

        !SoilID
        allocate (Me%SoilID             (ILB:IUB,JLB:JUB,KLB:KUB))        
        
        !Speed Up Matrtix
        allocate (Me%CalculateHead      (ILB:IUB,JLB:JUB,KLB:KUB))
        
        
        Me%SatK            = null_real
        Me%UnsatK          = null_real
        Me%UnsatK_X        = null_real
        Me%UnsatK_Y        = null_real
        Me%UnsatK_Z        = null_real
        
        Me%SoilID          = null_int
        
        Me%CalculateHead   = .false.

        !Retention Curve-------------------------------------------------------
        allocate (Me%RC%ThetaR (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%RC%ThetaS (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%RC%ThetaF (ILB:IUB,JLB:JUB,KLB:KUB))
        
        Me%RC%ThetaR        = null_real
        Me%RC%ThetaS        = null_real
        Me%RC%ThetaF        = null_real
        
        !Converge method arrays------------------------------------------------
        allocate(Me%CV%HeadIni          (ILB:IUB,JLB:JUB,KLB:KUB))

        allocate(Me%CV%ThetaIni         (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%CV%ThetaOld         (ILB:IUB,JLB:JUB,KLB:KUB))
                
        !Velocities------------------------------------------------------------
        allocate(Me%UnsatVelU           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%UnsatVelV           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%UnsatVelW           (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%UnsatVelWFinal      (ILB:IUB,JLB:JUB,KLB:KUB))
        
        allocate(Me%FluxU               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxV               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxW               (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxWFinal          (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxUAcc            (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxVAcc            (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxWAcc            (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FluxWAccFinal       (ILB:IUB,JLB:JUB,KLB:KUB))
        
        if (Me%ExtVar%ConstructEvaporation) then
            allocate(Me%EvaporationFlux (ILB:IUB,JLB:JUB))
            allocate(Me%EfectiveEVAP    (ILB:IUB,JLB:JUB))
        endif
        
        allocate(Me%InfiltrationVelocity(ILB:IUB,JLB:JUB        ))

        Me%UnsatVelU            = 0.0
        Me%UnsatVelV            = 0.0
        Me%UnsatVelW            = 0.0
        Me%UnsatVelWFinal       = 0.0
        
        Me%FluxU                = 0.0
        Me%FluxV                = 0.0
        Me%FluxW                = 0.0
        Me%FluxWFinal           = 0.0
        if (Me%ExtVar%ConstructEvaporation) then
            Me%EvaporationFlux      = 0.0
            Me%EfectiveEVAP         = null_real
        endif

        Me%InfiltrationVelocity = null_real

        !Flow to Channel
        allocate(Me%iFlowToChannels                 (ILB:IUB,JLB:JUB))
        allocate(Me%iFlowToChannelsLayer    (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%lFlowToChannels                 (ILB:IUB,JLB:JUB))
        allocate(Me%lFlowToChannelsLayer    (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%FlowToChannelsTopLayer          (ILB:IUB,JLB:JUB))
        allocate(Me%FlowToChannelsBottomLayer       (ILB:IUB,JLB:JUB))
        
        allocate(Me%lFlowDischarge    (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%iFlowDischarge    (ILB:IUB,JLB:JUB,KLB:KUB))
        Me%lFlowDischarge        = 0.0   !Sets values initially to zero, so 
        Me%iFlowDischarge        = 0.0   !model can run without Dis
                
        if (Me%SoilOpt%ImposeBoundary) then
            if (Me%SoilOpt%ImposeBoundaryWalls) then
                allocate(Me%Boundary%ImposedBoundaryLevel (ILB:IUB,JLB:JUB))
                Me%Boundary%ImposedBoundaryLevel = null_real
                allocate(Me%Boundary%BoundaryCells        (ILB:IUB,JLB:JUB))
                Me%Boundary%BoundaryCells        = 0
                allocate(Me%iFlowBoundaryWalls              (ILB:IUB,JLB:JUB,KLB:KUB))
                Me%iFlowBoundaryWalls            = 0.0
            endif
            if (Me%SoilOpt%ImposeBoundaryBottom) then
                allocate(Me%iFlowBoundaryBottom                (ILB:IUB,JLB:JUB))
                Me%iFlowBoundaryBottom          = 0.0
            endif
        endif
        Me%iFlowToChannels      = 0.0
        Me%iFlowToChannelsLayer = 0.0
        Me%lFlowToChannels      = 0.0
        Me%lFlowToChannelsLayer = 0.0
        Me%FlowToChannelsBottomLayer = Me%WorkSize%KLB
        Me%FlowToChannelsTopLayer    = Me%WorkSize%KUB       

        Me%FluxUAcc      = 0.0
        Me%FluxVAcc      = 0.0
        Me%FluxWAcc      = 0.0
        Me%FluxWAccFinal = 0.0
    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ConstructDischarges

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        character(len=StringLength)                 :: DischargeName
        real                                        :: CoordinateX, CoordinateY
        logical                                     :: CoordinatesON, IgnoreOK
        integer                                     :: Id, Jd, Kd, dn, DischargesNumber, nc
        integer                                     :: STAT_CALL
        type (T_Lines),   pointer                   :: LineX
        type (T_Polygon), pointer                   :: PolygonX
        integer, dimension(:),   pointer            :: VectorI, VectorJ, VectorK
        integer                                     :: SpatialEmission, nCells, DischVertical

        call Construct_Discharges(Me%ObjDischarges, Me%ObjTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR01' 
                
        call GetDischargesNumber(Me%ObjDischarges, DischargesNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR02' 

        do dn = 1, DischargesNumber

            call GetDischargesGridLocalization(Me%ObjDischarges, dn,            &
                                               DischVertical = DischVertical,   &
                                               KGrid         = Kd,              &            
                                               CoordinateX   = CoordinateX,     &
                                               CoordinateY   = CoordinateY,     & 
                                               CoordinatesON = CoordinatesON,   &
                                               STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR03' 
                    
            call GetDischargesIDName (Me%ObjDischarges, dn, DischargeName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR03' 

            if (CoordinatesON) then
                
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordinateX, CoordinateY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR04' 

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreDischarge(Me%ObjDischarges, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR05' 

                    if (IgnoreOK) then
                        write(*,*) 'Discharge outside the domain - ',trim(DischargeName),' - ',trim(Me%ModelName)
                        cycle
                    else
                        stop 'ModulePorousMedia - ConstructDischarges - ERR06' 
                    endif

                endif

                call CorrectsCellsDischarges(Me%ObjDischarges, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR07' 
                    
            endif

            call GetDischargeSpatialEmission(Me%ObjDischarges, dn, LineX, PolygonX, &
                                             SpatialEmission, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR08' 
            
            
            
            !ATTENTION - NEED TO VERIFY IF DISCHARGES ARE COLLINEAR.
            !Do not allow with properties since the flow used in PMP is not distributed by discharges
            !and will be accounted with flow duplicating
            if (SpatialEmission == DischPoint_) then
 
                call GetDischargesGridLocalization(Me%ObjDischarges, dn,            &
                                                   DischVertical = DischVertical,   &
                                                   Igrid         = Id,              &
                                                   JGrid         = Jd,              &
                                                   STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR09' 

                if (Me%ExtVar%BasinPoints(Id,Jd) /= WaterPoint) then
                    call TryIgnoreDischarge(Me%ObjDischarges, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR10' 

                    write(*,*) 'Discharge outside the domain I=',Id,' J=',Jd,'Model name=',trim(Me%ModelName)

                    if (IgnoreOK) then
                        write(*,*) 'Discharge in a land cell - ',trim(DischargeName),' - ',trim(Me%ModelName)
                        cycle
                    else
                        stop 'ModulePorousMedia - ConstructDischarges - ERR11' 
                    endif
                endif

                nCells    = 1
                allocate(VectorI(nCells), VectorJ(nCells), VectorK(nCells))
                VectorJ(nCells) = Jd
                VectorI(nCells) = Id

            else

                if (SpatialEmission == DischLine_) then
                    call GetCellZInterceptByLine(Me%ObjHorizontalGrid, LineX,               &
                                                 Me%ExtVar%BasinPoints, VectorI, VectorJ,   &
                                                 VectorK, nCells, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR12' 

                    if (nCells < 1) then
                        write(*,*) 'Discharge line intercept 0 cells'       
                        stop 'ModulePorousMedia - ConstructDischarges - ERR13' 
                    endif

                endif 


                if (SpatialEmission == DischPolygon_) then
                    call GetCellZInterceptByPolygon(Me%ObjHorizontalGrid, PolygonX,         &
                                                 Me%ExtVar%BasinPoints, VectorI, VectorJ,   &
                                                 VectorK, nCells, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR14' 

                    if (nCells < 1) then
                        write(*,*) 'Discharge contains 0 center cells'       
                        write(*,*) 'Or the polygon is to small and is best to a discharge in a point or'
                        write(*,*) 'the polygon not define properly'
                        stop 'ModulePorousMedia - ConstructDischarges - ERR15' 
                    endif

                endif


            endif

c1:         select case (DischVertical)
                        
                case (DischLayer_)

                    VectorK(:) = Kd

                case (DischDepth_)
        
                    write(*,*) "VERTICAL DISCHARGE option not active - Depth  =",DischDepth_
                    write(*,*) 'This option is not active'
                    stop 'ModulePorousMedia - ConstructDischarges - ERR170'

                case (DischBottom_)

                    call GetGeometryKFloor(Me%ObjGeometry,                  &
                                           Z = Me%ExtVar%KFloor,            &
                                           STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                              &
                        stop 'ModulePorousMedia - ConstructDischarges - ERR180.'

n1:                 do nC =1, nCells

                        VectorK(nC) = Me%ExtVar%Kfloor(VectorI(nC), VectorJ(nC))

                    enddo n1
    
                    call UnGetGeometry(Me%ObjGeometry,                      &
                                           Me%ExtVar%KFloor,                &
                                           STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                              &
                        stop 'ModulePorousMedia - ConstructDischarges - ERR190.'


                case (DischSurf_)
    
                    VectorK(:) = Me%WorkSize%KUB
    
                case (DischUniform_)
                    !do not do nothing                     
                case default
                    write(*,*) "VERTICAL DISCHARGE option not known ", DischVertical

                    write(*,*) "The known options are : "," Bottom=",DischBottom_," Surface=",DischSurf_,  &
                                                          " Layer =",DischLayer_, " Depth  =",DischDepth_, &
                                                          " Uniform=",DischUniform_
                    stop 'ModulePorousMedia - ConstructDischarges - ERR200'

            end select c1
                        
            if (SpatialEmission /= DischPoint_) then


                call SetLocationCellsZ (Me%ObjDischarges, dn, nCells, VectorI, VectorJ, VectorK, STAT= STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR210' 

            else
                if (DischVertical == DischBottom_ .or. DischVertical == DischSurf_) then
                    call SetLayer (Me%ObjDischarges, dn, VectorK(nCells), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ConstructDischarges - ERR220' 
                endif
                deallocate(VectorI, VectorJ, VectorK)
            endif

        enddo

   
    end subroutine ConstructDischarges
    
    !--------------------------------------------------------------------------
    
    subroutine ReadSoilTypes
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, nSoils
        logical                                     :: SoilTypeFound
        integer                                     :: ClientNumber, iflag, SoilID
        real                                        :: thf
        
        !Counts the number of SoilTypes
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        nSoils = 0
doS:    do
            !Gets soils type
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginsoiltype>',                    &
                                        block_end       = '<endsoiltype>',                      &
                                        BlockFound      = SoilTypeFound,                        &   
                                        STAT            = STAT_CALL)
SF:         if (STAT_CALL == SUCCESS_ .and. SoilTypeFound) then
                nSoils = nSoils + 1
           else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR02'
               
                exit doS

            end if SF
        enddo doS            
        
        allocate(Me%SoilTypes(nSoils))
        
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)

doH:    do 
            !Gets soils type
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginsoiltype>',                    &
                                        block_end       = '<endsoiltype>',                      &
                                        BlockFound      = SoilTypeFound,                        &   
                                        STAT            = STAT_CALL)
HF:         if (STAT_CALL == SUCCESS_ .and. SoilTypeFound) then

                !Reads ID
                call GetData(SoilID, Me%ObjEnterData,  iflag,                                   &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'ID',                                             &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR10'
                if (iflag /= 1) then
                    write(*,*)'Missing ID in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR15'
                endif

                !Reads Theta R
                call GetData(Me%SoilTypes(SoilID)%ThetaR, Me%ObjEnterData,  iflag,              &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'THETA_R',                                        &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR20'
                if (iflag /= 1) then
                    write(*,*)'Missing THETA_R in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR30'
                endif
    
                !Reads Theta S
                call GetData(Me%SoilTypes(SoilID)%ThetaS, Me%ObjEnterData,  iflag,              &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'THETA_S',                                        &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR30'
                if (iflag /= 1) then
                    write(*,*)'Missing THETA_S in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR40'
                endif

                !Reads NFit
                call GetData(Me%SoilTypes(SoilID)%NFit, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'N_FIT',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR50'
                if (iflag /= 1) then
                    write(*,*)'Missing N_FIT in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR60'
                endif


                !Reads LFit
                call GetData(Me%SoilTypes(SoilID)%LFit, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'L_FIT',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR70'
                if (iflag /= 1) then
                    write(*,*)'Missing L_FIT in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR80'
                endif
                
                !Reads Alfa
                call GetData(Me%SoilTypes(SoilID)%Alfa, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'ALPHA',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR90'
                if (iflag /= 1) then
                    write(*,*)'Missing ALPHA in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR100'
                endif

                !Reads SatK
                call GetData(Me%SoilTypes(SoilID)%SatK, Me%ObjEnterData,  iflag,                &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'SAT_K',                                          &
                             ClientModule   = 'ModulePorousMedia',                              &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR110'
                if (iflag /= 1) then
                    write(*,*)'Missing SAT_K in Soil Type definition'
                    stop 'ReadSoilTypes - ModulePorousMedia - ERR120'
                endif

                !Calculates MFIT
                Me%SoilTypes(SoilID)%MFit = 1.0 - (1.0 / Me%SoilTypes(SoilID)%NFit)
                
                !Calculates Oversat Head Slope
                thf = ThetaF_ (Me%SoilTypes(SoilID)%ThetaS - Me%CV%LimitThetaHi, SoilID)
                Me%SoilTypes(SoilID)%OverSatSlope = - Head_(thf, SoilID) / Me%CV%LimitThetaHi

           else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadSoilTypes - ModulePorousMedia - ERR16'
               
                exit doH

            end if HF
        enddo doH
                
    
    end subroutine ReadSoilTypes

    !--------------------------------------------------------------------------   
  
    subroutine ReadIntegrationPropertyConfig (Info, Is3D, BlockBegin, BlockEnd, AllocateArrays, ClientNumber)
        
        !Arguments-------------------------------------------------------------
        type (T_IntegrationInfo), pointer           :: Info
        logical, intent(IN)                         :: Is3D
        character(*), intent(IN)                    :: BlockBegin
        character(*), intent(IN)                    :: BlockEnd
        logical, intent(IN), optional               :: AllocateArrays
        integer, intent(IN )                        :: ClientNumber
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, block_i
        logical                                     :: BlockFound
        integer                                     :: iflag
        integer                                     :: ILB, IUB, JLB,  JUB 
        integer                                     :: KLB, KUB
        logical                                     :: AllocateArrays_
    
        !----------------------------------------------------------------------
        if (present(AllocateArrays)) then
            AllocateArrays_ = AllocateArrays
        else
            AllocateArrays_ = .true.
        endif
        
        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KUB = Me%Size%KUB
        KLB = Me%Size%KLB
        
        call ExtractBlockFromBlock (Me%ObjEnterData,                                            &
                                    ClientNumber        = ClientNumber,                         &
                                    block_begin         = BlockBegin,                           &
                                    block_end           = BlockEnd,                             &
                                    BlockInBlockFound   = BlockFound,                           &
                                    STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR010'
        if (BlockFound) then
            Info%Yes = .true.
            
            call GetData (Info%ByLayer, Me%ObjEnterData,  iflag,                                &                
                          SearchType     = FromBlockInBlock,                                    &
                          keyword        = 'BY_LAYER',                                          &
                          default        = .false.,                                             &
                          ClientModule   = 'ModulePorousMedia',                                 &
                          STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR020'
            if (Info%ByLayer .and. AllocateArrays_) then
                if (Is3D) then
                    allocate (Info%Old3d (ILB:IUB,JLB:JUB,KLB:KUB))
                    Info%Old3d = 0.0
                    
                    allocate (Info%Field3d (ILB:IUB,JLB:JUB,KLB:KUB))
                    Info%Field3d = 0.0
                else
                    allocate (Info%Old2d (ILB:IUB,JLB:JUB))
                    Info%Old2d = 0.0
                    
                    allocate (Info%Field2d (ILB:IUB,JLB:JUB))
                    Info%Field2d = 0.0
                endif
            endif
            
            call GetData (Info%ByHorizon, Me%ObjEnterData,  iflag,                              &
                          SearchType     = FromBlockInBlock,                                    &
                          keyword        = 'BY_HORIZON',                                        &
                          default        = .false.,                                             &
                          ClientModule   = 'ModulePorousMedia',                                 &
                          STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR030'
            if (Info%ByHorizon) then
                
                call GetNumberOfBlocks (Me%ObjEnterData, '<begininthorizon>', '<endinthorizon>',&
                                        FromBlockInBlock_,                                      &
                                        Info%HorizonsCount,                                     &
                                        ClientNumber, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR040'        
                if (Info%HorizonsCount <= 0) stop ''                                
                
                allocate (Info%Horizons (Info%HorizonsCount))
                                
                do block_i=1, Info%HorizonsCount
                
                    call ExtractBlockFromBlockFromBlock (Me%ObjEnterData,                           &
                                                         ClientNumber        = ClientNumber,        &
                                                         block_begin         = '<begininthorizon>', &
                                                         block_end           = '<endinthorizon>',   &
                                                         BlockInBlockInBlockFound = BlockFound,     &
                                                         STAT                = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR050'
                    if (BlockFound) then
                    
                        call GetData (Info%Horizons(block_i)%Name,                              &
                                      Me%ObjEnterData,  iflag,                                  &
                                      SearchType     = FromBlockInBlockInBlock,                 &
                                      keyword        = 'NAME',                                  &
                                      ClientModule   = 'ModulePorousMedia',                     &
                                      STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR060' 
                    
                        call GetData (Info%Horizons(block_i)%StartLayer,                        &
                                      Me%ObjEnterData,  iflag,                                  &
                                      SearchType     = FromBlockInBlockInBlock,                 &
                                      keyword        = 'START',                                 &
                                      ClientModule   = 'ModulePorousMedia',                     &
                                      STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR070' 
                    
                        call GetData (Info%Horizons(block_i)%EndLayer,                          &
                                      Me%ObjEnterData,  iflag,                                  &
                                      SearchType     = FromBlockInBlockInBlock,                 &
                                      keyword        = 'END',                                   &
                                      ClientModule   = 'ModulePorousMedia',                     &
                                      STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR080'
                    
                        if (AllocateArrays_) then
                            if (Is3D) then
                                allocate (Info%Horizons(block_i)%Old3D (ILB:IUB,JLB:JUB,KLB:KUB))
                                Info%Horizons(block_i)%Old3D = 0.0
                                
                                allocate (Info%Horizons(block_i)%Field3D (ILB:IUB,JLB:JUB,KLB:KUB))
                                Info%Horizons(block_i)%Field3D = 0.0
                            else
                                allocate (Info%Horizons(block_i)%Old2D (ILB:IUB,JLB:JUB))
                                Info%Horizons(block_i)%Old2D = 0.0
                                
                                allocate (Info%Horizons(block_i)%Field2D (ILB:IUB,JLB:JUB))
                                Info%Horizons(block_i)%Field2D = 0.0                                
                            endif
                        endif
                    !    
                    !else
                    !    stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR090'
                    endif
                    
                enddo
                               
            endif
                       
        else
            Info%Yes = .false.
        endif
        
        call RewindBlock(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadIntegrationPropertyConfig - ModulePorousMedia - ERR100'
    
    end subroutine ReadIntegrationPropertyConfig
    
    !--------------------------------------------------------------------------    
    
    subroutine ReadIntegrationConfiguration
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: BlockFound
        integer                                     :: ClientNumber
        type (T_IntegrationInfo), pointer           :: aux
        integer                                     :: ILB, IUB, JLB,  JUB 
        integer                                     :: KLB, KUB

        !----------------------------------------------------------------------
             
        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KUB = Me%Size%KUB
        KLB = Me%Size%KLB        
        
        Me%IntegrationOutput%AccTime = 0.0
        
        call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
        
        call ExtractBlockFromBuffer (Me%ObjEnterData,                                           &
                                     ClientNumber    = ClientNumber,                            &
                                     block_begin     = '<beginintegration>',                    &
                                     block_end       = '<endintegration>',                      &
                                     BlockFound      = BlockFound,                              &
                                     STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop ''
        if (.not. BlockFound) &
            stop ''
                   
        aux => Me%IntegrationOutput%WaterContent
        call ReadIntegrationPropertyConfig (aux,                                                &
                                            .true.,                                             &
                                            '<beginwatercontent>',                              &
                                            '<endwatercontent>', ClientNumber = ClientNumber)
               
        aux => Me%IntegrationOutput%RelativeWaterContent
        call ReadIntegrationPropertyConfig (aux,                                                &
                                            .true.,                                             &
                                            '<beginrelativewatercontent>',                      &
                                            '<endrelativeatercontent>',                         &
                                            AllocateArrays = .false.,                           &
                                            ClientNumber = ClientNumber)
        
        aux => Me%IntegrationOutput%WaterTable
        call ReadIntegrationPropertyConfig (aux,                                                &
                                            .false.,                                            &
                                            '<beginwatertable>',                                &
                                            '<endwatertable>', ClientNumber = ClientNumber) 
        if (aux%yes) &
            allocate (Me%OldUGWaterLevel2D (ILB:IUB,JLB:JUB))
        
        !Not implemented yet
        aux => Me%IntegrationOutput%Infiltration
        call ReadIntegrationPropertyConfig (aux,                                                &
                                            .false.,                                            &
                                            '<begininfiltration>',                              &
                                            '<endinfiltration>', ClientNumber = ClientNumber)         

        aux => Me%IntegrationOutput%BoundaryBottom
        call ReadIntegrationPropertyConfig (aux,                                                &
                                            .false.,                                            &
                                            '<beginboundarybottom>',                            &
                                            '<endboundarybottom>', ClientNumber = ClientNumber)
 
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadIntegrationConfiguration - ModulePorousMedia - ERR030'                
    
    end subroutine ReadIntegrationConfiguration

    !--------------------------------------------------------------------------       
    
    subroutine InitialFields

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, WLNumber
        integer                                     :: KLB, KUB, ClientNumber
        logical                                     :: HorizonFound, BlockFound
        logical                                     :: AllOK
        integer                                     :: nProps, iflag
        type (T_PropertyID)                         :: WaterLevelID
        !type (T_PropertyID)                         :: ImpermeableFractionID
        integer, allocatable, dimension(:)          :: LayerControl
        integer, dimension(:,:,:), pointer          :: AuxPointsToFill
        real, dimension(:,:,:), pointer             :: AuxSoilID
        integer                                     :: i, j, k
        character(LEN = StringLength)               :: string
        type (T_PropertyID)                         :: ID
        
        allocate(LayerControl(Me%WorkSize%KLB: Me%WorkSize%KUB))
        LayerControl = 0

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)

doH:    do 
            !Gets soils horizons
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginhorizon>',                     &
                                        block_end       = '<endhorizon>',                       &
                                        BlockFound      = HorizonFound,                         &   
                                        STAT            = STAT_CALL)
HF:         if (STAT_CALL == SUCCESS_ .and. HorizonFound) then

                !Reads lower layer of horizon
                call GetData(KLB, Me%ObjEnterData,  iflag,                                      &
                            SearchType     = FromBlock,                                         &
                            keyword        = 'KLB',                                             &
                            ClientModule   = 'ModuleFillMatrix',                                &
                            STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR08'
                if (iflag /= 1) then
                    write(*,*)'Missing KLB in Horizon definition'
                    stop 'InitialFields - ModulePorousMedia - ERR09'
                endif

                !Reads upper layer of horizon       
                call GetData(KUB, Me%ObjEnterData,  iflag,                                      &
                            SearchType     = FromBlock,                                         &
                            keyword        = 'KUB',                                             &
                            ClientModule   = 'ModuleFillMatrix',                                &
                            STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR10'
                if (iflag /= 1) then
                    write(*,*)'Missing KUB in Horizon definition'
                    stop 'InitialFields - ModulePorousMedia - ERR11'
                endif

                do k = KLB, KUB
                    if (LayerControl(k) /= 0) then
                        write(*,*)'Inconsistent horizon definition. Layer:',k
                    else
                        LayerControl(k) = 1
                    endif
                enddo
                
                !Allocates Aux Matrix
                allocate(AuxPointsToFill (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))

                do i = Me%Size%ILB, Me%Size%IUB
                do j = Me%Size%JLB, Me%Size%JUB
                do k = Me%Size%KLB, Me%Size%KUB
                    if (k >= KLB .and. k <= KUB) then
                        AuxPointsToFill(i, j, k) = Me%ExtVar%WaterPoints3D(i, j, k)
                    else
                        AuxPointsToFill(i, j, k) = 0
                    endif
                enddo
                enddo
                enddo                

                nProps         = 0  !Control Variable
doSP:           do
                    !Gets properties of horizon
                    call ExtractBlockFromBlock(Me%ObjEnterData,                             &
                                                ClientNumber        = ClientNumber,         &
                                                block_begin         = '<beginproperty>',    &
                                                block_end           = '<endproperty>',      &
                                                BlockInBlockFound   = BlockFound,           &   
                                                STAT                = STAT_CALL)
                    if (STAT_CALL .EQ. SUCCESS_ .and. BlockFound) then                                                  

                        !Get Name of Property
                        call GetData(String, Me%ObjEnterData, iflag,        &
                                    SearchType     = FromBlockInBlock,     &
                                    keyword        ='NAME',                &
                                    ClientModule   ='ModulePOrousMedia',   &
                                    STAT           = STAT_CALL)  

                        select case (trim(adjustl(string)))
            
                        !Water contents
                        case (char_Theta        )                 
        
                            call ConstructFillMatrix(PropertyID           = ID,                               &
                                                     EnterDataID          = Me%ObjEnterData,                  &
                                                     TimeID               = Me%ObjTime,                       &
                                                     HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                                     GeometryID           = Me%ObjGeometry,                   &
                                                     ExtractType          = FromBlockInBlock,                 &
                                                     PointsToFill3D       = AuxPointsToFill,                  &
                                                     Matrix3D             = GetPointer(Me%Theta),             &
                                                     TypeZUV              = TypeZ_,                           &
                                                     STAT                 = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR97'
        
                            call KillFillMatrix (ID%ObjFillMatrix, STAT = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR97a'

                        !SoilID
                        case (char_SoilID)

                            !Construct Fill Matrix cant read integer. So allocate aux matrix here.
                            allocate(AuxSoilID (Me%Size%ILB:Me%Size%IUB,                                      &
                                                Me%Size%JLB:Me%Size%JUB,                                      &
                                                Me%Size%KLB:Me%Size%KUB))
                                                
                            AuxSoilID = null_real

                            call ConstructFillMatrix(PropertyID           = ID,                               &
                                                     EnterDataID          = Me%ObjEnterData,                  &
                                                     TimeID               = Me%ObjTime,                       &
                                                     HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                                     GeometryID           = Me%ObjGeometry,                   &
                                                     ExtractType          = FromBlockInBlock,                 &
                                                     PointsToFill3D       = AuxPointsToFill,                  &
                                                     Matrix3D             = AuxSoilID,                        &
                                                     TypeZUV              = TypeZ_,                           &
                                                     STAT                 = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR98'
        
                            call KillFillMatrix (ID%ObjFillMatrix, STAT = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR98a'

                            do k = Me%Size%KLB, Me%Size%KUB
                            do j = Me%Size%JLB, Me%Size%JUB
                            do i = Me%Size%ILB, Me%Size%IUB
                                if (AuxPointsToFill(i, j, k) == 1) then
                                    Me%SoilID(i, j, k) = NINT(AuxSoilID(i, j, k))
                                endif
                            enddo
                            enddo
                            enddo
                            
                            deallocate (AuxSoilID)

                        !Invalid
                        case default

                            write(*,*)'Invalid Property', trim(adjustl(string))
                            stop 'InitialFields - ModulePorousMedia - ERR99'

                        end select                    
                    
                        nProps = nProps + 1
                    else
                        if (nProps /= 2) then
                            write(*,*)'Soil hydraulic properties incorrected defined'
                            stop 'InitialFields - ModulePorousMedia - ERR15'
                        endif
                        exit doSP
                    endif
                enddo doSP
                
                deallocate (AuxPointsToFill) 
                !Fills Constant Soil Variables
                

            else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR16'
               
                exit doH

            end if HF
        enddo doH

        AllOK = .true.
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            if (LayerControl(k) /= 1) then
                AllOK = .false.
            endif
        enddo
        
        deallocate(LayerControl)

        if (.not. AllOK) then
            write(*,*)'Inconsistent horizon definition.'
            stop 'InitialFields - ModulePorousMedia - ERR61'
        endif
        
        do k = Me%Size%KLB, Me%Size%KUB
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExtVar%Waterpoints3D(i, j, k) == 1) then
                if (Me%SoilID(i, j, k) < 0) then
                    write(*,*)'Soils not defined for [i,j,k]', i, j, k
                    AllOK = .false.
                endif
            endif
        enddo
        enddo
        enddo
        
        if (.not. AllOK) then
            write(*,*)'Inconsistent soil definition.'
            stop 'InitialFields - ModulePorousMedia - ERR61a'
        endif
        

        !Sets Matrixes of ThetaR and ThetaS        
        do k = Me%Size%KLB, Me%Size%KUB
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExtVar%Waterpoints3D(i, j, k) == 1) then
                Me%RC%ThetaS(i, j, k) = Me%SoilTypes(Me%SoilID(i, j, k))%ThetaS
                Me%RC%ThetaR(i, j, k) = Me%SoilTypes(Me%SoilID(i, j, k))%ThetaR
                Me%SatK     (i, j, k) = Me%SoilTypes(Me%SoilID(i, j, k))%SatK
            endif
        enddo
        enddo
        enddo
        
      
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR70'

        !Constructs Water Level
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                    ClientNumber    = WLNumber,                             &
                                    block_begin     = '<beginwaterlevel>',                  &
                                    block_end       = '<endwaterlevel>',                    &
                                    BlockFound      = BlockFound,                           &   
                                    STAT            = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if (.not. BlockFound) then
                write(*,*)'Missing Block <beginwaterlevel> / <endwaterlevel>'
                stop 'InitialFields - ModulePorousMedia - ERR80'
            endif
            
            call ConstructFillMatrix  ( PropertyID           = WaterLevelID,                &
                                        EnterDataID          = Me%ObjEnterData,             &
                                        TimeID               = Me%ObjTime,                  &
                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                        ExtractType          = FromBlock,                   &
                                        PointsToFill2D       = Me%ExtVar%BasinPoints,       &
                                        Matrix2D             = Me%UGWaterLevel2D,           &
                                        TypeZUV              = TypeZ_,                      &
                                        STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR81'
            
            call KillFillMatrix       (WaterLevelID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR82'

        else
            stop 'InitialFields - ModulePorousMedia - ERR90'
        endif
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR70'

        !Constructs Impermeable Fraction
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                    ClientNumber    = WLNumber,                             &
                                    block_begin     = '<beginimpermeablefraction>',         &
                                    block_end       = '<endimpermeablefraction>',           &
                                    BlockFound      = BlockFound,                           &   
                                    STAT            = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if (.not. BlockFound) then
                write(*,*)'Missing Block <beginimpermeablefraction> / <endimpermeablefraction>'
                stop 'InitialFields - ModulePorousMedia - ERR100'
            endif
            
            call ConstructFillMatrix  ( PropertyID           = Me%ImpermeableFractionID,    &
                                        EnterDataID          = Me%ObjEnterData,             &
                                        TimeID               = Me%ObjTime,                  &
                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                        ExtractType          = FromBlock,                   &
                                        PointsToFill2D       = Me%ExtVar%BasinPoints,       &
                                        Matrix2D             = Me%ImpermeableFraction,      &
                                        TypeZUV              = TypeZ_,                      &
                                        STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR110'


            if (.not. Me%ImpermeableFractionID%SolutionFromFile) then
                call KillFillMatrix (Me%ImpermeableFractionID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR0130'
            endif

        else
            stop 'InitialFields - ModulePorousMedia - ERR140'
        endif

        call Block_Unlock(Me%ObjEnterData, WLNumber, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'InitalFields - ModulePorousMedia - ERR150'
        
    end subroutine InitialFields

    !--------------------------------------------------------------------------

    subroutine UpdateBoundaryConditions
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                    :: i,j,k
        !Begin-----------------------------------------------------------------
        
!        allocate (Mapping (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
!        if (STAT_CALL /= SUCCESS_) stop UpdateBoundaryConditions - ModulePorousMedia - ERR01
!    
!        SetMatrixValueAllocatable(Mapping, Me%Size, dble(0.0), Me%ExtVar%WaterPoints3D)
        
!        if (Me%SoilOpt%ImposeBoundaryLevelX) then
!            if (Me%SoilOpt%ImposedLevelLowerX .gt. -95.) then
!                j = Me%Size%JLB
!                do i= Me%Size%ILB, Me%Size%IUB
!                do k= Me%Size%KLB, Me%Size%KUB
!                    !in saturated cells assuming hydrostatic pressure
!                    !in unsaturated cells gravity balances suction (field capcity) and
!                    !final head is uniform and equal to level
!                    Me%FinalHead(i,j,k) = Me%SoilOpt%ImposedLevelLowerX 
!                    
!                    !only mapping and compute faces in cells between level gradient 
!                    if ((-Me%ExtVar%SZZ(i,j,k) .le. max (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelLowerX)) &
!                       (-Me%ExtVar%SZZ(i,j,k) .ge. min (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelLowerX)))
!                        Mapping (i,j,k)     = 1
!                    endif
!                enddo
!                enddo
!            endif          
!            if (Me%SoilOpt%ImposedLevelUpperX .gt. -95.) then
!                j = Me%Size%JUB
!                do i= Me%Size%ILB, Me%Size%IUB
!                do k= Me%Size%KLB, Me%Size%KUB
!                    Me%FinalHead(i,j,k) = Me%SoilOpt%ImposedLevelUpperX
!                    
!                    if ((-Me%ExtVar%SZZ(i,j,k) .le. max (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelUpperX)) &
!                       (-Me%ExtVar%SZZ(i,j,k) .ge. min (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelUpperX)))
!                        Mapping (i,j,k)     = 1
!                    endif               
!                enddo
!                enddo
!            endif                     
!        endif
!        if (Me%SoilOpt%ImposeBoundaryLevelY) then
!            if (Me%SoilOpt%ImposedLevelLowerY .gt. -95.) then
!                i = Me%Size%ILB
!                do j= Me%Size%JLB, Me%Size%JUB
!                do k= Me%Size%KLB, Me%Size%KUB
!
!                    Me%FinalHead(i,j,k) = Me%SoilOpt%ImposedLevelLowerY 
!                    
!                    if ((-Me%ExtVar%SZZ(i,j,k) .le. max (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelLowerY)) &
!                       (-Me%ExtVar%SZZ(i,j,k) .ge. min (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelLowerY)))
!                        Mapping (i,j,k)     = 1
!                    endif               
!                enddo
!                enddo
!            endif          
!            if (Me%SoilOpt%ImposedLevelUpperY .gt. -95.) then
!                i = Me%Size%IUB
!                do j= Me%Size%JLB, Me%Size%JUB
!                do k= Me%Size%KLB, Me%Size%KUB
!
!                    Me%FinalHead(i,j,k) = Me%SoilOpt%ImposedLevelUpperY 
!                    
!                    if ((-Me%ExtVar%SZZ(i,j,k) .le. max (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelUpperY)) &
!                       (-Me%ExtVar%SZZ(i,j,k) .ge. min (Me%UGWaterLevel2D(i,j), Me%SoilOpt%ImposedLevelUpperY)))
!                        Mapping (i,j,k)     = 1
!                    endif  
!                enddo
!                enddo
!            endif                     
!        endif 

        do j= Me%Size%JLB, Me%Size%JUB
        do i= Me%Size%ILB, Me%Size%IUB
            if (Me%ExtVar%BasinPoints(i, j) == 0) then           
                do k= Me%Size%KLB, Me%Size%KUB

                    Me%FinalHead(i,j,k) = Me%Boundary%ImposedBoundaryLevel(i,j)
                    
                enddo
            endif
        enddo
        enddo
        
        !update mapping with boundary cells so that compute faces in boundary are open
        !(only the ones between saturated levels)
!        call UpDateWaterPoints3D(Me%ObjMap, Mapping, STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop UpdateBoundaryConditions - ModulePorousMedia - ERR01
        
                  
    end subroutine UpdateBoundaryConditions

    !--------------------------------------------------------------------------

    subroutine ReadBoundaryConditions
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        logical                                     :: BlockFound
        integer                                     :: ClientNumber
        !Begin-----------------------------------------------------------------

        !it will be changed to true if level is constant(in space and time - no piezometer)
        Me%Boundary%ImposedLevelConstant = .false.
        
        !it will be changed to true if at least a piezometer exists and has timeserie
        Me%Boundary%ImposedLevelInTime = .false.
        
        !Search for boundary block
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR10'

        !Constructs Impermeable Fraction
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                    ClientNumber    = ClientNumber,                         &
                                    block_begin     = '<begin_boundary>',                   &
                                    block_end       = '<end_boundary>',                     &
                                    BlockFound      = BlockFound,                           &   
                                    STAT            = STAT_CALL)
        if (STAT_CALL == SUCCESS_ .and. BlockFound) then
        
            call GetData (Me%Boundary%InterpolationMethod, Me%ObjEnterData, iflag,      &
                          SearchType   = FromBlock_,                                     &
                          keyword      ='INTERPOLATION_METHOD',                         &
                          ClientModule ='ModulePorousMedia',                            &
                          STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR20'
            if (iflag == 0) then
                write(*,*)'Interpolation Method not given'
                write(*,*)'Use Keyword INTERPOLATION_METHOD '   
                write(*,*)' [1 - Triangulation / 2 - IWD]'
                stop 'ReadBoundaryConditions - ModulePorousMedia - ERR30'
            endif
            
            !info needed for interpolation
            call DefineGridPoints
            
            select case (Me%Boundary%InterpolationMethod)

            case (Triangulation_)

                call GetData(Me%Boundary%Triangulation%FillOutsidePoints,                   &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock_,                                     &
                             keyword      = 'FILL_OUTSIDE_POINTS',                          &
                             ClientModule = 'ModulePorousMedia',                            &
                             default      = .false.,                                        &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR40'


                call GetData (Me%Boundary%Triangulation%OutputTriangles, Me%ObjEnterData, iflag,  &
                              SearchType   = FromBlock_,                                    &
                              keyword      ='WRITE_TRIANGLES',                              &
                              default      = .false.,                                       &
                              ClientModule ='ModulePorousMedia',                            &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR50'

                call GetData (Me%Boundary%Triangulation%FileName, Me%ObjEnterData, iflag,   &
                              SearchType   = FromBlock_,                                    &
                              keyword      ='TRIANGLES_FILE',                               &
                              ClientModule ='ModulePorousMedia',                            &
                              STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR60'
                
                if(Me%Boundary%Triangulation%OutputTriangles) then
                    !first instant for triangles output
                    Me%Boundary%Triangulation%Instant = 1
                endif
                
            case (InverseWeight_)

                call GetData(Me%Boundary%InverseWeight%MaxDistance,                         &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock_,                                     &
                             keyword      = 'MAX_DISTANCE',                                 &
                             ClientModule = 'ModulePorousMedia',                            &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR70'
                if (iflag == 0) then
                    write(*,*)'Max Distance not given'
                    write(*,*)'Use Keyword MAX_DISTANCE'
                    stop 'ReadBoundaryConditions - ModulePorousMedia - ERR70'
                endif
                
                call GetData(Me%Boundary%InverseWeight%IWDn,                                &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock_,                                     &
                             keyword      = 'IWD_N',                                        &
                             default      = 2.0,                                            &
                             ClientModule = 'ModulePorousMedia',                            &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR80'
            
            case (NoInterpolation_)  !one piezometer for all domain
            
                !do nothing
            
            case default
                write(*,*)'Invalid Interpolation Method for boundary condition'
                write(*,*)'Use Keyword INTERPOLATION_METHOD  '
                write(*,*)'[1 - Triangulation / 2 - IWD / 3 - No interpolation]'
                stop 'ReadBoundaryConditions - ModulePorousMedia - ERR90'
            end select

            call ConstructPiezometers(ClientNumber)
            
            !verify after piezometers build
            if ((Me%Boundary%InterpolationMethod == NoInterpolation_) .and. (Me%Boundary%NumberOfPiezometers > 1)) then
                write(*,*)'If using no interpolation for boundary method need'
                write(*,*)'to define only one piezometer'
                write(*,*)'Remove piezometers or change interpolation method for boundary condition'
                stop 'ReadBoundaryConditions - ModulePorousMedia - ERR100'            
            endif

        else

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadBoundaryConditions - ModulePorousMedia - ERR110'
            
        endif
        
        !separate from the last if block because of SearchType different
        if (.not. BlockFound) then
        
            !boundary values are given by the constant value everywhere
            Me%Boundary%ImposedLevelConstant = .true.
            
            call GetData(Me%Boundary%BoundaryValue,                                 &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType = FromFile,                                     &
                         keyword    = 'BOUNDARY_VALUE',                             &
                         ClientModule ='ModulePorousMedia',                         &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'GetUnSaturatedOptions - ModulePorousMedia - ERR120' 

            if (iflag == 0) then
                write(*,*)'if using boundary, BOUNDARY_VALUE must be defined in module PorousMedia'
                stop 'ReadDataFile - ModulePorousMedia - ERR0110'
            endif
        
            !boundary values are given by the constant value everywhere
            call SetMatrixValue (Me%Boundary%ImposedBoundaryLevel, Me%Size2D, Me%Boundary%BoundaryValue, Me%Boundary%BoundaryCells)

        endif
        
    end subroutine ReadBoundaryConditions

    !--------------------------------------------------------------------------


    subroutine ConstructPiezometers(ClientNumber)
        
        !Arguments-----------------------------------------------------
        integer                                     :: ClientNumber
        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        logical                                     :: BlockFound
        type(T_Piezometer), pointer                 :: NewPiezometer

        !Begin---------------------------------------------------------

        nullify(Me%Boundary%FirstPiezometer, NewPiezometer)
        Me%Boundary%NumberOfPiezometers = 0
        
do1 :   do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,              &
                                        '<<begin_piezometer>>', '<<end_piezometer>>',  &
                                        BlockFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPiezometers - ModulePorousMedia - ERR10'

if1 :       if(STAT_CALL .EQ. SUCCESS_) then 
   
if2 :           if (BlockFound) then

                    call AddPiezometer(Me%Boundary%FirstPiezometer, NewPiezometer)

                    call GetData(NewPiezometer%Name,                                &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlockInBlock,                   &
                                 keyword      = 'NAME',                             &
                                 ClientModule = 'ModulePorousMedia',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPiezometers - ModulePorousMedia - ERR20'

                    call GetData(NewPiezometer%Location%X,                      &
                                 Me%ObjEnterData, iflag,                        &
                                 SearchType   = FromBlockInBlock,               &
                                 keyword      = 'COORD_X',                      &
                                 ClientModule = 'ModulePorousMedia',            &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPiezometers - ModulePorousMedia - ERR30'

                    call GetData(NewPiezometer%Location%Y,                      &
                                 Me%ObjEnterData, iflag,                        &
                                 SearchType   = FromBlockInBlock,               &
                                 keyword      = 'COORD_Y',                      &
                                 ClientModule = 'ModulePorousMedia',            &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPiezometers - ModulePorousMedia - ERR40'
                    

                    call GetData(NewPiezometer%ValueType,                           &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlockInBlock,                   &
                                 keyword      = 'VALUE_TYPE',                       &
                                 Default      = SingleValue,                        &
                                 ClientModule = 'ModulePorousMedia',                &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPiezometers - ModulePorousMedia - ERR50'

                    select case(trim(NewPiezometer%ValueType))
                        
                        case(SingleValue)
                            
                            call GetData(NewPiezometer%DefaultValue,                &
                                         Me%ObjEnterData, iflag,                    &
                                         SearchType   = FromBlockInBlock,           &
                                         keyword      = 'DEFAULTVALUE',             &
                                         ClientModule = 'ModulePorousMedia',        &
                                         STAT         = STAT_CALL)        
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructPiezometers - ModulePorousMedia - ERR60'

                        case(FromTimeSerie)
                        
                            Me%Boundary%ImposedLevelInTime = .true.

                            call ReadPiezometerTimeSerie(NewPiezometer)

                    end select


                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop 'ConstructPiezometers - ModulePorousMedia - ERR70'
                        
                    exit do1    !No more blocks

                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'ConstructPiezometers - ModulePorousMedia - ERR80'
                    
            end if if1
        end do do1

        
    end subroutine ConstructPiezometers

    !------------------------------------------------------------------    

    subroutine AddPiezometer (FirstPiezometer, ObjPiezometer)

        !Arguments-------------------------------------------------------------
        type (T_Piezometer), pointer                   :: FirstPiezometer
        type (T_Piezometer), pointer                   :: ObjPiezometer

        !Local-----------------------------------------------------------------
        type (T_Piezometer), pointer                   :: NewPiezometer
        type (T_Piezometer), pointer                   :: PreviousPiezometer
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewPiezometer)
        nullify  (NewPiezometer%Next)

        Me%Boundary%NumberOfPiezometers = Me%Boundary%NumberOfPiezometers + 1

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstPiezometer)) then
            FirstPiezometer         => NewPiezometer
            ObjPiezometer           => NewPiezometer
        else
            PreviousPiezometer      => FirstPiezometer
            ObjPiezometer           => FirstPiezometer%Next
            do while (associated(ObjPiezometer))
                PreviousPiezometer  => ObjPiezometer
                ObjPiezometer       => ObjPiezometer%Next
            enddo
            ObjPiezometer           => NewPiezometer
            PreviousPiezometer%Next => NewPiezometer
        endif


    end subroutine AddPiezometer
    
    !--------------------------------------------------------------------------

    subroutine ReadPiezometerTimeSerie(NewPiezometer)   
    
        !Arguments-------------------------------------------------------------
        type (T_Piezometer), pointer                   :: NewPiezometer

        !Local-----------------------------------------------------------------
        integer                                        :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------


        call GetData(NewPiezometer%TimeSerie%FileName,                  &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'FillMatrix',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadPiezometerTimeSerie - ModulePorousMedia - ERR01'


        call GetData(NewPiezometer%TimeSerie%DataColumn,                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'DATA_COLUMN',                      &
                     ClientModule = 'FillMatrix',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadPiezometerTimeSerie - ModulePorousMedia - ERR02'


        call StartTimeSerieInput(NewPiezometer%TimeSerie%ObjTimeSerie,  &
                                 NewPiezometer%TimeSerie%FileName,      &
                                 Me%ObjTime,                            &
                                 CheckDates = .false.,                  &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadPiezometerTimeSerie - ModulePorousMedia - ERR03'



    end subroutine ReadPiezometerTimeSerie
    
    !------------------------------------------------------------------

    subroutine ModifyBoundaryLevel

        !Local-----------------------------------------------------------------
        type (T_Piezometer), pointer                :: Piezometer
        integer                                     :: i,j, STAT_CALL, CHUNK
        real                                        :: Nominator, Denominator, dist
!        character(5)                                :: AuxChar
        !Begin-----------------------------------------------------------------

        call SetMatrixValue (Me%Boundary%ImposedBoundaryLevel, Me%Size2D, null_real, Me%ExtVar%BasinPoints)

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        
        select case (Me%Boundary%InterpolationMethod)

        case (Triangulation_)

            if (Me%Boundary%NumberOfPiezometers < 3) then
                write (*,*)'Need at least 3 piezometers'
                write (*,*)'for triangulation'
                stop 'ModifyBoundaryLevel - ModulePorousMedia - ERR01'
            endif

            allocate(Me%Boundary%Triangulation%Nodes%X(1:Me%Boundary%NumberOfPiezometers))
            allocate(Me%Boundary%Triangulation%Nodes%Y(1:Me%Boundary%NumberOfPiezometers))
            allocate(Me%Boundary%Triangulation%Nodes%Z(1:Me%Boundary%NumberOfPiezometers))

            Piezometer => Me%Boundary%FirstPiezometer
        
            i = 0

            do while(associated(Piezometer))

                i = i + 1
                
                call UpDatePiezometerValue(Piezometer, Me%ExtVar%Now)
                Me%Boundary%Triangulation%Nodes%X(i) = Piezometer%Location%X
                Me%Boundary%Triangulation%Nodes%Y(i) = Piezometer%Location%Y
                Me%Boundary%Triangulation%Nodes%Z(i) = Piezometer%DefaultValue
            
                Piezometer => Piezometer%Next

            end do

            !Constructs Triangulation
            call ConstructTriangulation (Me%ObjTriangulation,                           &
                                         Me%Boundary%NumberOfPiezometers,               &
                                         Me%Boundary%Triangulation%Nodes%X,             &
                                         Me%Boundary%Triangulation%Nodes%Y,             &
                                         Me%Boundary%Triangulation%OutputTriangles,     &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyBoundaryLevel - ModulePorousMedia - ERR10'


            call SetHeightValues(Me%ObjTriangulation, Me%Boundary%Triangulation%Nodes%Z, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyBoundaryLevel - ModulePorousMedia - ERR20'

            !This cycle can not be parallelized because of the function Interpolation
!            !$OMP PARALLEL PRIVATE(I,J)
!            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if(Me%Boundary%BoundaryCells(i, j) == BasinPoint) then

                    Me%Boundary%ImposedBoundaryLevel(i, j) = InterPolation(Me%ObjTriangulation,            &
                                                              Me%Boundary%GridPoint(i,j)%X,                &
                                                              Me%Boundary%GridPoint(i,j)%Y,                &
                                                              Me%Boundary%Triangulation%FillOutsidePoints, &
                                                              Default = null_real,                         &
                                                              STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyBoundaryLevel - ModulePorousMedia - ERR30'
                
                end if

            enddo
            enddo
!            !$OMP END DO NOWAIT
!            !$OMP END PARALLEL

            !Only write at the end since triangles are always the same (time span was removed)
            if(Me%Boundary%Triangulation%OutputTriangles .and. Me%ExtVar%Now == Me%EndTime) then
                call WriteTriangles(trim(adjustl(Me%Boundary%Triangulation%FileName)))
            endif

            !Kills Triangulation
            call KillTriangulation (Me%ObjTriangulation, STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyBoundaryLevel - ModulePorousMedia - ERR40'

            !Deallocates Matrixes 
            deallocate(Me%Boundary%Triangulation%Nodes%X)
            deallocate(Me%Boundary%Triangulation%Nodes%Y)
            deallocate(Me%Boundary%Triangulation%Nodes%Z)

        case (InverseWeight_)

            !Interpolates Values
            
            !$OMP PARALLEL PRIVATE(I,J, Nominator, Denominator, Piezometer, dist)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if(Me%Boundary%BoundaryCells(i, j) == BasinPoint) then

                    Nominator   = 0.0
                    Denominator = 0.0
                    Piezometer  => Me%Boundary%FirstPiezometer
DoPiezometers:      do while(associated(Piezometer))
                        !if (Piezometer%TimeSerieHasData) then
                            dist = sqrt((Me%Boundary%GridPoint(i,j)%X - Piezometer%Location%X)**2.0 + &
                                        (Me%Boundary%GridPoint(i,j)%Y - Piezometer%Location%Y)**2.0)

                            if (dist > 0.0) then    
                                if (dist < Me%Boundary%InverseWeight%MaxDistance) then 
                                    Nominator   = Nominator   + Piezometer%DefaultValue / (dist ** Me%Boundary%InverseWeight%IWDn)
                                    Denominator = Denominator + 1.0 / (dist ** Me%Boundary%InverseWeight%IWDn)
                                endif
                            else
                                Nominator   = Piezometer%DefaultValue
                                Denominator = 1.0
                                exit DoPiezometers
                            endif
                        !endif
                        Piezometer => Piezometer%Next
                    enddo DoPiezometers

                    if (Denominator .lt. AllMostZero) then
                        write (*,*)'Insufficient data avaliable'
                        write (*,*)'Increase MAX_DISTANCE or get more piezometers '
                        write (*,*)'Point [i, j]', i, j
                        stop 'ModifyBoundaryLevel - ModulePorousMedia - ERR01'
                    endif
                    Me%Boundary%ImposedBoundaryLevel(i, j) = Nominator / Denominator

                end if

            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

        case (NoInterpolation_)  !one piezometer (timeserie everywehre)
        
        
            Piezometer => Me%Boundary%FirstPiezometer
        
            call UpDatePiezometerValue(Piezometer, Me%ExtVar%Now)
            
            !boundary values are given by the timeserie value everywhere
            call SetMatrixValue (Me%Boundary%ImposedBoundaryLevel, Me%Size2D, Piezometer%DefaultValue, Me%Boundary%BoundaryCells)
                    
        end select

    end subroutine ModifyBoundaryLevel

    !--------------------------------------------------------------------------

    subroutine UpDatePiezometerValue(Piezometer, CurrentTime)
        
        !Arguments-------------------------------------------------------------
        type(T_Piezometer), pointer                 :: Piezometer
        type(T_Time)                                :: CurrentTime

        !Local-----------------------------------------------------------------
        logical                                     :: TimeCycle
        type (T_Time)                               :: Time1, Time2, InitialDate
        real                                        :: Value1, Value2
!        real                                        :: dt1, dt2
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        
        select case(trim(Piezometer%ValueType))
            
            case(SingleValue)

                !Remains Constant
                !Piezometer%TimeSerieHasData = .true.

            case(FromTimeSerie)

                call GetTimeSerieInitialData(Piezometer%TimeSerie%ObjTimeSerie, InitialDate, &
                                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'UpDatePiezometerValue - ModulePorousMedia - ERR01'

                
                if (CurrentTime >= InitialDate) then

                    !Gets Value for current Time
                    call GetTimeSerieValue (Piezometer%TimeSerie%ObjTimeSerie,                &
                                            CurrentTime,                                      &
                                            Piezometer%TimeSerie%DataColumn,                  &
                                            Time1, Value1, Time2, Value2, TimeCycle,          &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'UpDatePiezometerValue - ModulePorousMedia - ERR10'

                    if (TimeCycle) then
                        Piezometer%DefaultValue     = Value1
                        !Piezometer%TimeSerieHasData = .true.
                    else
                        
                        !Do not use TimeSpan. Use all values to interpolate
!                        dt1 = CurrentTime - Time1
!                        dt2 = Time2 - CurrentTime
!
!                        if (dt1 <= Me%Time%MaxTimeSpan .and. dt2 <= Me%Time%MaxTimeSpan) then

                            !Piezometer%TimeSerieHasData = .true.

                            !Interpolates Value for current instant
                            call InterpolateValueInTime(CurrentTime, Time1, Value1, Time2, Value2, Piezometer%DefaultValue)
                        
!                        else
!
!                            Piezometer%TimeSerieHasData = .false.
!
!                        endif
                    
                                                
                    endif

                else
                    write(*,*) 'Piezometer time serie does not have data' 
                    write(*,*) 'for the beggining of the simulation'
                    write(*,*) 'Piezometer name: ', Piezometer%Name
                    stop 'UpDatePiezometerValue - ModulePorousMedia - ERR20'
                    !Piezometer%TimeSerieHasData = .false.

                endif

            case default

        end select

    end subroutine UpDatePiezometerValue

    !------------------------------------------------------------------

    subroutine WriteTriangles (TrianglesFileName)

        !Arguments-------------------------------------------------------------
        character(len=*)                :: TrianglesFileName
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_CALL, nNodes
        integer                         :: UnitNumber, iT, nTriangles
        real,    dimension(:), pointer  :: XT, YT, ZT
        integer, dimension(:), pointer  :: V1, V2, V3

        !Begin-----------------------------------------------------------------

        !Get the number of triangles
        call GetNumberOfTriangles   (Me%ObjTriangulation, nTriangles)

        !Allocates space for the Triangle vertices and gets them
        allocate(V1(nTriangles))
        allocate(V2(nTriangles))
        allocate(V3(nTriangles))

        call GetTriangleList (Me%ObjTriangulation, v1, v2, v3, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - ModulePorousMedia - ERR10'


        !Gets nodes effictive used and the reordered nodes 
        call GetNumberOfNodes (Me%ObjTriangulation, nNodes, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - ModulePorousMedia - ERR20'

        allocate(XT(nNodes))
        allocate(YT(nNodes))
        allocate(ZT(nNodes))

        call GetNodesList   (Me%ObjTriangulation, XT, YT, ZT, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteTriangles - ModulePorousMedia - ERR30'


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

    !--------------------------------------------------------------------------

    subroutine DefineGridPoints

        !Local-----------------------------------------------------------------
        integer                             :: i, j, CHUNK
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW
        integer                             :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG, NLRD
        integer                             :: GRID_COORD, CoordType, STAT_CALL

        !Begin-----------------------------------------------------------------

        allocate(Me%Boundary%GridPoint(Me%WorkSize%ILB-1:Me%WorkSize%IUB+1, &
                                       Me%WorkSize%JLB-1:Me%WorkSize%JUB+1))

        !Gets Coordinates in use
        call GetGridCoordType(Me%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGridPoints - ModulePorousMedia - ERR10'
        
        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        if    (CoordType == SIMPLE_GEOG)then
            
            call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                 &
                                          GridLatitudeConn  = Me%ExtVar%YY_IE,  &
                                          GridLongitudeConn = Me%ExtVar%XX_IE,  &
                                          STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGridPoints - ModulePorousMedia - ERR20'

        elseif(CoordType == UTM        .or. CoordType == MIL_PORT .or. &
               CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(Me%ObjHorizontalGrid,                        &
                                   XX_IE = Me%ExtVar%XX_IE,                     &
                                   YY_IE = Me%ExtVar%YY_IE,                     &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGridPoints - ModulePorousMedia - ERR30'

        else

            write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
            stop 'DefineGridPoints - ModulePorousMedia - ERR40'

        end if

        CHUNK = ChunkI !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do i = Me%WorkSize%ILB,  Me%WorkSize%IUB
        do j = Me%WorkSize%JLB , Me%WorkSize%JUB
                       
            XSW = Me%ExtVar%XX_IE(i, j)
            YSW = Me%ExtVar%YY_IE(i, j)
            XSE = Me%ExtVar%XX_IE(i, j + 1)
            YSE = Me%ExtVar%YY_IE(i, j + 1)
            XNE = Me%ExtVar%XX_IE(i + 1, j + 1)
            YNE = Me%ExtVar%YY_IE(i + 1, j + 1)
            XNW = Me%ExtVar%XX_IE(i + 1, j)
            YNW = Me%ExtVar%YY_IE(i + 1, j)

            Me%Boundary%GridPoint(i,j)%X = (XSW+XNW+XNE+XSE) / 4.
            Me%Boundary%GridPoint(i,j)%Y = (YSW+YNW+YNE+YSE) / 4.


        end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGridPoints - ModulePorousMedia - ERR50'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGridPoints - ModulePorousMedia - ERR60'

    end subroutine DefineGridPoints
    
    !------------------------------------------------------------------
    
    subroutine StartWithFieldCapacity

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i,j,k
        real                                        :: Hinf
        integer                                     :: WTCell
        real                                        :: inf_border 
        
        do j= Me%WorkSize%JLB, Me%WorkSize%JUB
        do i= Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == 1) then
            
                !Searches cell were initial WT is locateboundary
                WTCell  = 0
                do k = Me%WorkSize%KUB, Me%ExtVar%KFloor(i, j), -1
                    inf_border = - Me%ExtVar%SZZ(i,j,k-1)
                    if (Me%UGWaterLevel2D(i,j) .gt. inf_border ) then
                        WTCell  = k
                        exit
                    endif                
                enddo        
                
                if (WTCell < Me%ExtVar%KFloor(i, j)) WTCell = Me%ExtVar%KFloor(i, j)
            
                !Saturates cells below WTable
                do k = Me%ExtVar%KFloor(i, j), WTCell-1
                    Me%Theta(i, j, k) = Me%RC%ThetaS(i, j, k)
                enddo


                if (Me%SoilOpt%StartWithFieldCapacity) then
                
                    !Sets cell above WT to Field Capacity
                    do k = WTCell, Me%WorkSize%KUB

                        if (k == WTCell) then
                            !Half of the cell distance
                            Hinf = - Me%ExtVar%DWZ(i,j,k) * 0.5
                        else
                            Hinf = - (Me%ExtVar%DZZ(i,j,k-1) - Hinf)
                        endif
                        
                        Me%Theta(i, j, k) = Theta_(Hinf, Me%SoilID(i, j, k))
                        
                    enddo
                
                endif
                    
            endif
        enddo
        enddo        
                                

    end subroutine StartWithFieldCapacity

    !--------------------------------------------------------------------------

    subroutine ComputeSoilFieldCapacity

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: i,j,k
        real                                                :: Hinf
        !Begin-----------------------------------------------------------------

        do j= Me%WorkSize%JLB, Me%WorkSize%JUB
        do i= Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == 1) then
            
                do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB

                    if (k == Me%ExtVar%KFloor(i, j)) then
                        !Half of the cell distance
                        Hinf = - Me%ExtVar%DWZ(i,j,k) * 0.5
                    else
                        Hinf = - (Me%ExtVar%DZZ(i,j,k-1) - Hinf)
                    endif
                    
                    Me%ThetaField(i, j, k) = Theta_(Hinf, Me%SoilID(i, j, k))
                    
                enddo
                
                    
            endif
        enddo
        enddo        

    end subroutine ComputeSoilFieldCapacity

    !--------------------------------------------------------------------------

    subroutine ConstructIntegrationHDF5Output
    
        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB,KLB,KUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        real, dimension(:, :), pointer                      :: BottomData

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        call GetHDF5FileAccess (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5 (Me%ObjIntegrationHDF5, trim(Me%Files%IntegrationHDFFile)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR010'

        !Write the Horizontal Grid
        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjIntegrationHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR020'

        !Sets limits for next write operations
        call HDF5SetLimits (Me%ObjIntegrationHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR030'

        !Writes the Grid
        call HDF5WriteData (Me%ObjIntegrationHDF5, "/Grid", "Topography", "m",                      &
                            Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR040'

        call GetGridData(Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR050'

        call HDF5WriteData (Me%ObjIntegrationHDF5, "/Grid", "Bathymetry", "m",                      &
                            Array2D = BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR060'

        call HDF5WriteData (Me%ObjIntegrationHDF5, "/Grid", "BasinPoints", "-",                     &
                            Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR070'

        call HDF5WriteData (Me%ObjIntegrationHDF5, "/Grid", "WaterPoints3D", "-",                  &
                            Array3D = Me%ExtVar%WaterPoints3d, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR080'
                
        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjIntegrationHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR090'

        call UnGetGridData (Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructIntegrationHDF5Output - ModulePorousMedia - ERR100'    
    
    end subroutine ConstructIntegrationHDF5Output
    
    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Output        

        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB,KLB,KUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        real, dimension(:, :), pointer                      :: BottomData

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR03'

        !Writes the Grid
        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "Topography", "m",                    &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR04'

        call GetGridData(Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR05'

        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR06'

        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "BasinPoints", "-",                   &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR07'

        !Water Points
        if (Me%OutPut%Yes) then
            call HDF5WriteData   ( Me%ObjHDF5,  "/Grid", "WaterPoints3D", "-",                  &
                                   Array3D = Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR08'
        endif              
                
        !Flushes All pending HDF5 commands
        call HDF5FlushMemory    (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR09'

        call UnGetGridData(Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR10'

    end subroutine ConstructHDF5Output

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        integer                                             :: nProperties
        integer                                             :: STAT_CALL
        integer                                             :: iflag, i
        character(len=StringLength)                         :: TimeSerieLocationFile
        integer                                             :: TimeSerieNumber, dn, Id, Jd
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        character(len=StringLength)                         :: TimeSerieName
        character(len= 20)                                  :: comment
        
        
        !Begin-----------------------------------------------------------------
        
        nProperties = 11
        if (Me%ObjDrainageNetwork /= 0) then
            nProperties = nProperties + 1
            if (Me%SoilOpt%DNLink == GWFlowToChanByLayer_) then
                nProperties = nProperties + 1
            endif
        endif       
        if (Me%ExtVar%ConstructEvaporation) then
            nProperties = nProperties + 1
        endif
        if (Me%ExtVar%ConstructTranspiration) then
            nProperties = nProperties + 1
        endif


        !Allocates PropertyList
        allocate(PropertyList(nProperties))

        !Fills up PropertyList
        PropertyList(1)  = 'Theta'
        PropertyList(2)  = 'relative water content'
        PropertyList(3)  = 'VelW [m/s]'
        PropertyList(4)  = 'VelW_Corr [m/s]'
        PropertyList(5)  = 'InF_Vel [m/s]'
        PropertyList(6)  = 'Head [m]'
        PropertyList(7)  = 'Conductivity [m/s]'
        PropertyList(8)  = 'level water table [m]'
        PropertyList(9)  = 'water table depth [m]'
        PropertyList(10)  = 'Hydro Pressure [m]'
        PropertyList(11) = 'Final Head [m]'
!        PropertyList(12) = 'Infiltration Column [m]'
        i = 11
        if (Me%ObjDrainageNetwork /= 0) then
            i = i + 1
            PropertyList(i)  = 'GW flow to river total [m3/s]'        
            if( Me%SoilOpt%DNLink == GWFlowToChanByLayer_) then
                i = i + 1
                PropertyList(i)  = 'GW flow to river layer [m3/s]'
            endif
        endif        
        if (Me%ExtVar%ConstructEvaporation) then
            i = i + 1
            PropertyList(i) = 'Surface Evaporation Flux [m3/s]'
        endif
        if (Me%ExtVar%ConstructTranspiration) then
            i = i + 1
            PropertyList(i) = 'Transpiration Flux [m3/s]'
        endif

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMedia',                                &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR01' 

        if (iflag == 1) then
            Me%OutPut%TimeSerieON = .true.
        else
            Me%OutPut%TimeSerieON = .false.
        endif

        write(comment, *) Me%DomainCellsNumber

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srp",                                        &
                            WaterPoints3D = Me%ExtVar%WaterPoints3D,                    &
                            ModelName = Me%ModelName,                                   &
                            Comment = 'Number of cells in 3D domain: '//trim(comment),  &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR02' 

        !Deallocates PropertyList
        deallocate(PropertyList)
       
        
        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR03'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR04'
            
            call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR04'
            
i1:         if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR05'

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR06'

                    if (IgnoreOK) then
                        write(*,*) 'Time Serie outside the domain - ',trim(TimeSerieName),' - ',trim(Me%ModelName)
                        cycle
                    else
                        stop 'ConstructTimeSerie - PorousMedia - ERR07'
                    endif

                endif


                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR08'

            endif i1

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      LocalizationI   = Id,                             &
                                      LocalizationJ   = Jd,                             & 
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR09'

            if (Me%ExtVar%WaterPoints3D(Id, Jd, Me%WorkSize%KUB) /= WaterPoint) then
                 write(*,*) 'Time Serie in a land cell - ',trim(TimeSerieName),' - ',trim(Me%ModelName)
            endif


        enddo
        
       
       
    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------

    subroutine ConstructProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL, iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        character(len=StringLength), dimension(:,:), pointer:: PropertyList
        integer                                             :: nProperties

        nProperties = 3 !Theta, ThetaF, Head

        !Allocates PropertyList
        allocate(PropertyList(nProperties, 2))

        !Fills up PropertyList
        PropertyList(1, 1) = 'PropI'
        PropertyList(2, 1) = 'PropII'

        PropertyList(1, 2) = "m3/m3"
        PropertyList(2, 2) = "m3/m3"


        !----------------------------------------------------------------------

        call GetData(Me%OutPut%ProfileON,                                               &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OUTPUT_PROFILE',                                   &
                     ClientModule = 'ModulePorousMedia',                                &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR01' 


        if (Me%OutPut%ProfileON) then

            call GetData(TimeSerieLocationFile,                                             &
                         Me%ObjEnterData,iflag,                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModulePorousMedia',                                &
                         Default      = Me%Files%DataFile,                                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR02' 
            
            !Starts Profile for Theta / ThetaF
            call StartProfile  (ProfileID       = Me%ObjProfile,                            &
                                ObjTime         = Me%ObjTime,                               &
                                ProfileDataFile = trim(TimeSerieLocationFile),              &
                                WaterPoints2D   = Me%ExtVar%BasinPoints,                    &
                                nProperties     = 3,                                        &
                                PropertyList    = PropertyList,                             &
                                KUB             = Me%WorkSize%KUB,                          &
                                ClientName      = "PorousMedia",                            &
                                STAT            = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR03' 

        endif

        deallocate (PropertyList)
        
    end subroutine ConstructProfileOutput


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !---------------------------------------------------------------------------

    subroutine GetNextPorousMediaDT (ObjPorousMediaID, DT, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real, intent(OUT)                               :: DT
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            DT        = Me%CV%NextDT

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNextPorousMediaDT

    !--------------------------------------------------------------------------
    
    subroutine GetPotentialInfiltration (ObjPorousMediaID, PotInfiltration, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: PotInfiltration
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            PotInfiltration => Me%ExtVar%InfiltrationColumn

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetPotentialInfiltration

    !--------------------------------------------------------------------------    

    subroutine GetInfiltration (ObjPorousMediaID, Infiltration, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: Infiltration
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            Infiltration => Me%Infiltration

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetInfiltration

    !--------------------------------------------------------------------------

    subroutine GetEfectiveEVTP (ObjPorousMediaID, EfectiveEVTP,STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: EfectiveEVTP
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            EfectiveEVTP => Me%EfectiveEVTP

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetEfectiveEVTP

    !--------------------------------------------------------------------------

    subroutine GetGWFlowToChannels (ObjPorousMediaID, FlowToChannels, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real, dimension(:, :), pointer                  :: FlowToChannels
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            FlowToChannels => Me%iFlowToChannels

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWFlowToChannels

    !--------------------------------------------------------------------------

    subroutine GetGWFlowToChannelsByLayer (ObjPorousMediaID, FlowToChannelsLayer, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real, dimension(:, :, :), pointer               :: FlowToChannelsLayer
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            FlowToChannelsLayer => Me%iFlowToChannelsLayer

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWFlowToChannelsByLayer

    !--------------------------------------------------------------------------

    subroutine GetGWToChannelsLayers (ObjPorousMediaID, GWToChannelsBottomLayer, GWToChannelsTopLayer, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, dimension(:, :), pointer               :: GWToChannelsBottomLayer
        integer, dimension(:, :), pointer               :: GWToChannelsTopLayer
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            GWToChannelsBottomLayer => Me%FlowToChannelsBottomLayer

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            GWToChannelsTopLayer => Me%FlowToChannelsTopLayer

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWToChannelsLayers

    !--------------------------------------------------------------------------

    subroutine GetGWFlowOption (ObjPorousMediaID, DrainageNetworkLink, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer                                         :: DrainageNetworkLink
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DrainageNetworkLink = Me%SoilOpt%DNLink

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWFlowOption

    !--------------------------------------------------------------------------

    subroutine GetGWLayer (ObjPorousMediaID, GWLayer, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, dimension(:,:), pointer                :: GWLayer
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            GWLayer => Me%UGCell

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWLayer

    !--------------------------------------------------------------------------

    subroutine GetGWLayerOld (ObjPorousMediaID, GWLayerOld, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, dimension(:,:), pointer                :: GWLayerOld
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        

        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            GWLayerOld => Me%UGCell_Old

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGWLayerOld

    !--------------------------------------------------------------------------

    subroutine GetFlowDischarge (ObjPorousMediaID, FlowDischarge, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real, dimension(:, :, :), pointer               :: FlowDischarge
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjPorousMediaID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            FlowDischarge => Me%iFlowDischarge

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine GetFlowDischarge    
    
    !--------------------------------------------------------------------------

    subroutine GetPMTotalDischargeFlowVolume (ObjPorousMediaID, Volume, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8)                                         :: Volume
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        call Ready(ObjPorousMediaID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !call Read_Lock(mPorousMedia_, Me%InstanceID)
            Volume = Me%TotalDischargeFlowVolume

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
    
    end subroutine GetPMTotalDischargeFlowVolume
    
    !--------------------------------------------------------------------------
    
    subroutine GetPorousMediaTotalStoredVolume (ObjPorousMediaID, TotalStoredVolume, LossToGround, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8)                                         :: TotalStoredVolume
        real(8), optional                               :: LossToGround
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            TotalStoredVolume = Me%TotalStoredVolume

            if (present(LossToGround)) LossToGround      = Me%LossToGround

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetPorousMediaTotalStoredVolume

    !-------------------------------------------------------------------------
    subroutine GetGeometryInstance (ObjPorousMediaID, GeometryID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, intent(OUT), optional                  :: STAT
        integer                                         :: GeometryID

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            !call Read_Lock(mPorousMedia_, Me%InstanceID)  !Changed 

            GeometryID = Me%ObjGeometry

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetGeometryInstance

    !--------------------------------------------------------------------------


    subroutine GetFluxU (ObjPorousMediaID, FlowU, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:,:,:), pointer              :: FlowU
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            FlowU => Me%FluxUAcc

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFluxU

    !--------------------------------------------------------------------------

    subroutine GetFluxV (ObjPorousMediaID, FlowV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:,:,:), pointer              :: FlowV
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            FlowV => Me%FluxVAcc

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFluxV

    !--------------------------------------------------------------------------

    subroutine GetFluxW (ObjPorousMediaID, FluxWFinal, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:,:,:), pointer              :: FluxWFinal
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ! Lucia
            !FlowWOld => Me%FluxWold
            
            !changed by Eduardo Jauch
            FluxWFinal => Me%FluxWAccFinal

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetFluxW

    !--------------------------------------------------------------------------

    subroutine GetUnsatV (ObjPorousMediaID, UnsatV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,dimension(:,:,:), pointer                  :: UnsatV                        
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            !UnsatW => Me%UnsatVelWOld
            UnsatV => Me%UnsatVelV

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetUnsatV
    
    !--------------------------------------------------------------------------

    subroutine GetUnsatU (ObjPorousMediaID, UnsatU, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,dimension(:,:,:), pointer                  :: UnsatU                        
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            !UnsatW => Me%UnsatVelWOld
            UnsatU => Me%UnsatVelU

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetUnsatU
    
    !--------------------------------------------------------------------------

    subroutine GetUnsatW (ObjPorousMediaID, UnsatW, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,dimension(:,:,:), pointer                  :: UnsatW                        
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            UnsatW => Me%UnsatVelW

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetUnsatW

    !---------------------------------------------------------------------------

    subroutine GetUnsatWFinal (ObjPorousMediaID, UnsatW, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,dimension(:,:,:), pointer                  :: UnsatW                        
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            UnsatW => Me%UnsatVelWFinal

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetUnsatWFinal

    !---------------------------------------------------------------------------

    subroutine GetWaterColumn (ObjPorousMediaID, WaterColumn, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                          :: ObjPorousMediaID
        real(8), pointer, dimension(:,:) :: WaterColumn
        integer, intent(OUT), optional   :: STAT

        !Local-----------------------------------------------------------------
        integer :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            WaterColumn => Me%WaterColumn

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_    
    
    end subroutine GetWaterColumn

    !---------------------------------------------------------------------------

    subroutine GetWaterContent (ObjPorousMediaID, WC, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: WC
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            WC => Me%Theta

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetWaterContent

    !--------------------------------------------------------------------------

    subroutine GetOldWaterContent (ObjPorousMediaID, WCold, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: WCold
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            WCold => Me%CV%ThetaIni !Me%CV%ThetaOld

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetOldWaterContent

    !--------------------------------------------------------------------------

    subroutine GetHead (ObjPorousMediaID, Head, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: Head
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            Head => Me%Head

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetHead

    !--------------------------------------------------------------------------

    subroutine GetThetaR (ObjPorousMediaID, ThetaR, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: ThetaR
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ThetaR => Me%RC%ThetaR

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetThetaR

    !--------------------------------------------------------------------------

    subroutine GetThetaS (ObjPorousMediaID, ThetaS, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: ThetaS
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ThetaS => Me%RC%ThetaS

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetThetaS

    !--------------------------------------------------------------------------

    subroutine GetThetaF (ObjPorousMediaID, ThetaF, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: ThetaF
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ThetaF => Me%RC%ThetaF

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetThetaF

    !--------------------------------------------------------------------------

    subroutine GetUnsatK (ObjPorousMediaID, UnsatK, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: UnsatK
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            UnsatK => Me%UnsatK

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetUnsatK

    !--------------------------------------------------------------------------

    subroutine GetEvaporation (ObjPorousMediaID, Evaporation, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), pointer, dimension(:,:)                :: Evaporation
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            Evaporation => Me%EfectiveEVAP

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetEvaporation

    !--------------------------------------------------------------------------
    
    subroutine GetTranspiration (ObjPorousMediaID, Transpiration, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real   ,    pointer, dimension(:,:,:)           :: Transpiration
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            Transpiration => Me%ExtVar%TranspirationFlux

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetTranspiration

    !--------------------------------------------------------------------------        

    subroutine GetEVTPVolumes(ID, Evaporation, Transpiration, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ID
        real(8), intent(OUT), optional                  :: Evaporation
        real(8), intent(OUT), optional                  :: Transpiration
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !Begin------------------------------------------------------------------
        STAT_CALL = UNKNOWN_

        call Ready(ID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            if (present (Evaporation  )) Evaporation   = Me%AccEvapFromSoil
            if (present (Transpiration)) Transpiration = Me%AccTranspiration
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL
        !-----------------------------------------------------------------------
    
    end subroutine GetEVTPVolumes

    !---------------------------------------------------------------------------    
    
    subroutine GetPMBoundaryFlowVolume(ID, AccBoundaryFlowVolume, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ID
        real(8), intent(OUT)                            :: AccBoundaryFlowVolume
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !Begin------------------------------------------------------------------
        STAT_CALL = UNKNOWN_

        call Ready(ID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            AccBoundaryFlowVolume = Me%AccBoundaryFlowVolume
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL
        !-----------------------------------------------------------------------
    
    end subroutine GetPMBoundaryFlowVolume

    !---------------------------------------------------------------------------    

    subroutine GetPMStoredVolume (ID, StoredVolume, STAT)    

        !Arguments--------------------------------------------------------------
        integer                                         :: ID
        real(8), intent(OUT)                            :: StoredVolume        
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        integer                                         :: i, j, k

        !Begin------------------------------------------------------------------
        STAT_CALL = UNKNOWN_

        call Ready(ID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            StoredVolume = 0.0
            
            call GetWaterPoints3D(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPMStoredVolume - ModulePorousMedia - ERR010'  
            
            call GetGeometryVolumes(Me%ObjGeometry, VolumeZ = Me%ExtVar%CellVolume, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPMStoredVolume - ModulePorousMedia - ERR020'
            
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                
                    StoredVolume = StoredVolume +                               &
                                   Me%Theta(i,j,k) * Me%ExtVar%CellVolume(i,j,k)
                    
                endif

            enddo
            enddo
            enddo      
            
            call UnGetGeometry(Me%ObjGeometry, Me%ExtVar%CellVolume, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPMStoredVolume - ModulePorousMedia - ERR030'
            
            call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetPMStoredVolume - ModulePorousMedia - ERR040'            
            
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL
        !-----------------------------------------------------------------------
        
    end subroutine GetPMStoredVolume
    
    !---------------------------------------------------------------------------
    
    subroutine GetIgnoreWaterColumnOnEVAP (ObjPorousMediaID, IgnoreWaterColumnOnEvap, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        logical, intent(OUT)                            :: IgnoreWaterColumnOnEvap
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            IgnoreWaterColumnOnEvap = Me%SoilOpt%IgnoreWaterColumnOnEvap

            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetIgnoreWaterColumnOnEVAP

    !--------------------------------------------------------------------------    

    subroutine GetBoundaryImposed (ObjPorousMediaID, BoundaryOpen, BoundaryWallsOpen, BoundaryBottomOpen, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        logical, intent(OUT)                            :: BoundaryOpen
        logical, intent(OUT)                            :: BoundaryWallsOpen
        logical, intent(OUT)                            :: BoundaryBottomOpen
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryOpen       = Me%SoilOpt%ImposeBoundary
            BoundaryWallsOpen  = Me%SoilOpt%ImposeBoundaryWalls
            BoundaryBottomOpen = Me%SoilOpt%ImposeBoundaryBottom
            
            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetBoundaryImposed

    !--------------------------------------------------------------------------   

    subroutine GetBoundaryFluxWalls (ObjPorousMediaID, BoundaryFlux, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: BoundaryFlux
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            BoundaryFlux => Me%iFlowBoundaryWalls

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryFluxWalls

    !--------------------------------------------------------------------------   

    subroutine GetBoundaryFluxBottom (ObjPorousMediaID, BoundaryFlux, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:)                :: BoundaryFlux
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            BoundaryFlux => Me%iFlowBoundaryBottom

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryFluxBottom

    !--------------------------------------------------------------------------   

    subroutine GetBoundaryCells (ObjPorousMediaID, BoundaryCells, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer,   pointer, dimension(:,:)              :: BoundaryCells
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BoundaryCells => Me%Boundary%BoundaryCells

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoundaryCells

    !--------------------------------------------------------------------------
    

    subroutine GetThetaField (ObjPorousMediaID, ThetaField, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real,    pointer, dimension(:,:,:)              :: ThetaField
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMedia_, Me%InstanceID)
            
            ThetaField => Me%ThetaField

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetThetaField

    !--------------------------------------------------------------------------

    subroutine GetComputeSoilField (ObjPorousMediaID, ComputeSoilField, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        logical                                         :: ComputeSoilField
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            ComputeSoilField = Me%SoilOpt%ComputeSoilField

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetComputeSoilField

    !--------------------------------------------------------------------------

    subroutine GetLimitThetaLow (ObjPorousMediaID, LimitThetaLow, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real                                            :: LimitThetaLow
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            LimitThetaLow = Me%CV%LimitThetaLo

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetLimitThetaLow

    !--------------------------------------------------------------------------

    subroutine UnGetPorousMedia_R4(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R4")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R4

    !--------------------------------------------------------------------------

    subroutine UnGetPorousMedia_R8(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R8
       
    !----------------------------------------------------------------------------
    subroutine UnGetPorousMedia_R(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(8), pointer, dimension(:,:,:)               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R
    !-----------------------------------------------------------------------------

    subroutine UnGetPorousMedia_R1(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real(4), pointer, dimension(:,:,:)              :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMedia_R1

!-------------------------------------------------------------------------------

    subroutine UnGetPorousMedia_RI(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer                                         :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            array = 0
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_RI")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        end subroutine UnGetPorousMedia_RI

    ! ---------------------------------------------------------------------!

     subroutine UnGetPorousMedia_AI(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, pointer, dimension(:,:,:)              :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_AI")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        end subroutine UnGetPorousMedia_AI

    ! ---------------------------------------------------------------------!

     subroutine UnGetPorousMedia_AI2D(ObjPorousMediaID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        integer, pointer, dimension(:,:)                :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMedia_, Me%InstanceID, "UnGetPorousMedia_AI2D")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        end subroutine UnGetPorousMedia_AI2D

    ! ---------------------------------------------------------------------!
      
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyPorousMedia(ObjPorousMediaID,                         &
                                 InfiltrationColumn,                       &
                                 PotentialEvaporation,                     &
                                 ActualTranspiration,                      &
                                 STAT)

        !Arguments---------------------------------------------------------------
        integer,                            intent(IN)           :: ObjPorousMediaID       !IN
        real(8), dimension(:,:  ), pointer                       :: InfiltrationColumn     !IN
        real,    dimension(:,:  ), pointer, optional             :: PotentialEvaporation   !IN
        real,    dimension(:,:,:), pointer, optional             :: ActualTranspiration    !IN  
        integer, optional,                  intent(OUT)          :: STAT                   !OUT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        integer                                     :: DummyI
        real                                        :: DummyR
        !------------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "ModifyPorousMedia")

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then                        

            !Updates time
            call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMedia - ModulePorousMedia - ERR01'
            
            !Gets Time Step
            call GetComputeTimeStep (Me%ObjTime, DT = Me%ExtVar%DT)                


            !Update unsaturated computefaces
            call UpdateComputeFaces3D(  Map_ID         = Me%ObjMap,                       &
                                        DummyR         = DummyR,                          &
                                        DummyI         = DummyI,                          &
                                        STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMedia - ModulePorousMedia - ERR02'


            !Sets External Variables
            call ReadLockExternalVar
            
            !Check impermeability solution from file
            if (Me%ImpermeableFractionID%ObjFillMatrix /= 0 .and. Me%ImpermeableFractionID%SolutionFromFile) then
                call ModifyFillMatrix (FillMatrixID   = Me%ImpermeableFractionID%ObjFillMatrix,  &
                                       Matrix2D       = Me%ImpermeableFraction,                  &
                                       PointsToFill2D = Me%ExtVar%BasinPoints,                   &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMedia - ModulePorousMedia - ERR010'            
            endif
            
            !Points to Arguments
            Me%TranspirationExists = .false.
            Me%EvaporationExists   = .false.
            if (present(ActualTranspiration)) then
                Me%TranspirationExists               = .true.
                !m3/s
                Me%ExtVar%TranspirationFlux          => ActualTranspiration
            endif
            if (present(PotentialEvaporation)) then
                Me%EvaporationExists                 = .true.
                !m/s
                Me%ExtVar%PotentialEvaporationFlux   => PotentialEvaporation
            endif
            
            Me%ExtVar%InfiltrationColumn => InfiltrationColumn
       
            !Calculate flow in unsaturated part of soil
            call VariableSaturatedFlow (InfiltrationColumn)
            
            if (Me%SoilOpt%ImposeBoundary) then
                if (Me%SoilOpt%ImposeBoundaryWalls .and. Me%Boundary%ImposedLevelInTime) call ModifyBoundaryLevel
                call ModifyBoundaryFlux
            endif
            
            call CalculateUGWaterLevel
            
            call ComputeIntegration (Me%ExtVar%DT)

            !Output
            if (Me%OutPut%Yes .or. Me%OutPut%SurfaceOutput) call PorousMediaOutput                        
            if (Me%OutPut%TimeSerieON)                      call OutPutTimeSeries
            if (Me%OutPut%ProfileON  )                      call ProfileOutput
            if (Me%Output%BoxFluxes  )                      call ComputeBoxesWaterFluxes
            if (Me%IntegrationOutput%Yes)                   call IntegrationOutput
            
            !Restart Output
            if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%Now == Me%EndTime)) then
                if(Me%ExtVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    call WriteFinalFile
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                endif
            endif

            
            if (Me%SoilOpt%CheckGlobalMass) then
                call CalculateTotalStoredVolume
            endif

            call ReadUnLockExternalVar
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "ModifyPorousMedia")

    end subroutine ModifyPorousMedia

    !--------------------------------------------------------------------------

    subroutine ComputeIntegration(LocalDT)
        
        !Argument--------------------------------------------------------------
        real                                    :: LocalDT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        type (T_IntegrationInfo), pointer           :: Info
        type (T_IntegrationByHorizon), pointer      :: Horizon
        integer                                     :: hor_i
        real                                        :: aux, dwz
        
        !----------------------------------------------------------------------        

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
        
        Me%IntegrationOutput%AccTime = Me%IntegrationOutput%AccTime + LocalDT
        
        if (Me%IntegrationOutput%BoundaryBottom%yes) then
            Info => Me%IntegrationOutput%BoundaryBottom
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%ExtVar%BasinPoints (i, j) == 1) then
                    Info%Field2D (i,j) = Info%Field2D (i,j) + &
                                         (Me%iFlowBoundaryBottom(i,j) * LocalDT)
                else
                    Info%Field2D (i,j) = null_real
                endif
            enddo
            enddo 
        endif
        
        if (Me%IntegrationOutput%WaterTable%Yes) then    
            Info => Me%IntegrationOutput%WaterTable
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%ExtVar%BasinPoints (i, j) == 1) then
                    Info%Field2D (i,j) = Info%Field2D (i,j) + &
                                         ((Me%OldUGWaterLevel2D(i,j) + Me%UGWaterLevel2D(i,j)) / 2) * LocalDT
                else
                    Info%Field2D (i,j) = null_real
                endif
            enddo
            enddo 
        endif
        
        if (Me%IntegrationOutput%WaterContent%Yes) then    
            Info => Me%IntegrationOutput%WaterContent
            
            if (Info%ByLayer) then
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%ExtVar%BasinPoints (i, j) == 1) then
                        Info%Field3D (i,j,k) = Info%Field3D (i,j,k) + &
                                                ((Me%CV%ThetaIni(i,j,k) + Me%Theta(i,j,k)) / 2) * LocalDT
                    else
                        Info%Field3D (i,j,k) = null_real
                    endif
                enddo
                enddo 
                enddo
            endif
            
            !if (Info%ByHorizon) then
            !    do hor_i = 1 to Info%HorizonsCount
            !        Horizon => Info%Horizons(hor_i)
            !        
            !        do j = JLB, JUB
            !        do i = ILB, IUB
            !            
            !            if (Me%ExtVar%BasinPoints(i,j) == 1) then
            !                
            !                aux = 0.0
            !                dwz = 0.0
            !                do k = Horizon%StartLayer, Horizon%EndLayer                                
            !                    if (Me%ExtVar%AterPoints3D(i,j,k) == 1) then
            !                        aux = aux + ((Me%CV%ThetaIni(i,j,k) + Me%Theta(i,j,k)) / 2) * Me%ExtVar%DWZ(i, j, k)
            !                        dwz = dwz + Me%ExtVar%DWZ(i, j, k)
            !                    endif                                
            !                enddo
            !                
            !                Horizon%Field2D (i,j) = (aux / dwz) * LocalDT
            !            else
            !                Horizon%Field2D (i,j) = null_real
            !            endif
            !            
            !        enddo
            !        enddo
            !    enddo
            !endif
        endif        
    
    end subroutine ComputeIntegration
    
    !--------------------------------------------------------------------------

    subroutine VariableSaturatedFlow(InfiltrationColumn)
    
        !Arguments-------------------------------------------------------------
        real(8), dimension(:,:), pointer            :: InfiltrationColumn
        
        !Local-----------------------------------------------------------------
        logical                                     :: Restart
        integer                                     :: Niteration, iteration
        real                                        :: SumDT
        real                                        :: Zero = 0.0
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "VariableSaturatedFlow")

        !Stores initial values
        call SetMatrixValueAllocatable (Me%CV%ThetaIni,      Me%Size,   Me%Theta,          Me%ExtVar%WaterPoints3D)
        call SetMatrixValueAllocatable (Me%CV%HeadIni,       Me%Size,   Me%Head,           Me%ExtVar%WaterPoints3D)
        call SetMatrixValue (Me%WaterColumn,      Me%Size2D, InfiltrationColumn,Me%ExtVar%BasinPoints)

        !Time Integrated Values
        call SetMatrixValue (Me%Infiltration, Me%Size2D, dble(0.0), Me%ExtVar%BasinPoints)
        call SetMatrixValue (Me%EfectiveEVTP, Me%Size2D, dble(0.0), Me%ExtVar%BasinPoints)
        call SetMatrixValueAllocatable (Me%iFlowToChannels,  Me%Size2D, Zero, Me%ExtVar%BasinPoints)
        call SetMatrixValueAllocatable (Me%iFlowToChannelsLayer,  Me%Size, Zero, Me%ExtVar%WaterPoints3D)
        call SetMatrixValueAllocatable (Me%iFlowDischarge,Me%Size, Zero, Me%ExtVar%WaterPoints3D)
        if (Me%EvaporationExists) then
            call SetMatrixValue (Me%EfectiveEVAP, Me%Size2D, dble(0.0), Me%ExtVar%BasinPoints)
        endif  
       
        if (Me%CV%NextNIteration > 1 .and. Me%ExtVar%DT < (Me%CV%CurrentDT * Me%CV%NextNIteration)) then
            Niteration = max(aint(Me%ExtVar%DT / Me%CV%CurrentDT), 1.0)
        else          
            Niteration = Me%CV%NextNIteration    !DB
        endif
        
        iteration         = 1     
        Me%CV%CurrentDT   = Me%ExtVar%DT / Niteration
        SumDT             = 0.0        

        Me%FluxWAcc = 0.
        Me%FluxVAcc = 0.
        Me%FluxUAcc = 0.
        Me%FluxWAccFinal = 0.
        
        Me%AccEvapFromSoil  = 0.0
        Me%AccTranspiration = 0.0
        
        Me%TotalDischargeFlowVolume = 0.0
        
dConv:  do while (iteration <= Niteration)
        
            !Convergence Test
            call SetMatrixValueAllocatable (Me%CV%ThetaOld, Me%Size, Me%Theta, Me%ExtVar%WaterPoints3D)

            !Calculates Face Conductivities
            call Condutivity_Face
            
            call ComputeFinalHead

            !Calculates Water velocity
            call SoilWaterVelocity
           
            !Calculates Water Flux
            call SoilWaterFlux
            
            !Calculates Flux to channels
            if (Me%ObjDrainageNetwork /= 0 .and. Me%SoilOpt%CalcDrainageNetworkFlux) then
                !remove/add fluxes on top of aquifer
                if (Me%SoilOpt%DNLink == GWFlowToChanByCell_) then
                    !Areas for flux are around river wet perimeter
                    if (Me%SoilOpt%DNLinkAreaMethod == GWFlowAreaWetPerimeter_) then
                        call ExchangeWithDrainageNet_1
                    !Areas are around river wet perimeter and aquifer if aquifer higher then river
                    elseif (Me%SoilOpt%DNLinkAreaMethod == GWFlowAreaWetPerimeterAndAquifer_) then
                        call ExchangeWithDrainageNet_2
                    endif
                !remove/add fluxes on all layers that are influenced by level difference
                elseif (Me%SoilOpt%DNLink == GWFlowToChanByLayer_) then                    
                    call ExchangeWithDrainageNet_3
                endif
            endif
            
            if (Me%EvaporationExists) then
                call EvaporationFlux    ()
            endif

            !Inputs Water from discharges - moved to CalculateNewTheta
!            if (Me%SoilOpt%Discharges) then
!                call ModifyWaterDischarges               
!            endif
            
            !Calculates New Theta
            call CalculateNewTheta  ()

            !Checks for variation of theta values
            call CheckStability (Restart)

            !Vertical Continuty            
            if (Restart) then
                            
                Niteration        = Me%CV%NextNIteration
                Me%CV%CurrentDT   = Me%ExtVar%DT / Niteration
                call WriteDTLog_ML ('ModulePorousMedia', Niteration, Me%CV%CurrentDT)                    
                iteration         = 1
                SumDT             = 0.0
                
                !Restores Initial Values
                call SetMatrixValueAllocatable (Me%Theta,          Me%Size,   Me%CV%ThetaIni,      Me%ExtVar%WaterPoints3D)
                call SetMatrixValueAllocatable (Me%Head,           Me%Size,   Me%CV%HeadIni,       Me%ExtVar%WaterPoints3D)
                call SetMatrixValue (Me%WaterColumn,    Me%Size2D, InfiltrationColumn,  Me%ExtVar%BasinPoints)

                !Resets Time Integrated Values
                call SetMatrixValue (Me%Infiltration,   Me%Size2D, dble(0.0),           Me%ExtVar%BasinPoints)
                call SetMatrixValue (Me%EfectiveEVTP,   Me%Size2D, dble(0.0),           Me%ExtVar%BasinPoints)
                call SetMatrixValueAllocatable (Me%iFlowToChannels,Me%Size2D, Zero,                Me%ExtVar%BasinPoints)
                call SetMatrixValueAllocatable (Me%iFlowToChannelsLayer,Me%Size, Zero,           Me%ExtVar%WaterPoints3D)
                call SetMatrixValueAllocatable (Me%iFlowDischarge,Me%Size, Zero,           Me%ExtVar%WaterPoints3D)
                if (Me%EvaporationExists) then
                    call SetMatrixValue (Me%EfectiveEVAP,     Me%Size2D, dble(0.0),         Me%ExtVar%BasinPoints)
                endif                 !Resets Accumulated flows

                Me%FluxWAcc = 0.
                Me%FluxVAcc = 0.
                Me%FluxUAcc = 0.
                Me%FluxWAccFinal = 0.
                
                Me%AccEvapFromSoil  = 0.0
                Me%AccTranspiration = 0.0
                
                Me%TotalDischargeFlowVolume = 0.0
                                              
            else
                
                !Moves water in oversaturated cells up and down
                call VerticalContinuity

                !Calulates Heads / Conductivities from new Theta values
                call SoilParameters

                !Removes water due to infiltration
                call IntegrateValuesInTime(SumDT)
                
                !Accumulate flows
                call AccumulateFlows
                
                
                SumDT       = SumDT + Me%CV%CurrentDT
                iteration   = iteration + 1
                
            endif            
            
        enddo dConv
        
        call CalculateMeanFlows (Niteration)
        call InsertInfiltrationOnFluxMatrix
        
        if (Me%SoilOpt%WriteLog) call LogDT (Niteration)
        
        call ComputeNextDT(Niteration)

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "VariableSaturatedFlow")   

    end subroutine VariableSaturatedFlow
    
    !--------------------------------------------------------------------------

    subroutine ModifyWaterDischarges ()

        !Arguments--------------------------------------------------------------
        !real                                    :: LocalDT

        !Local------------------------------------------------------------------
        integer                                 :: iDis, nDischarges
        integer                                 :: i, j, k, kd, kmin, kmax
        real                                    :: SurfaceElevation    
        real                                    :: DischargeFlow, MaxFlow
        integer                                 :: STAT_CALL
        logical                                 :: IgnoreOK
        integer                                 :: DischVertical
        integer                                 :: nCells, n
        integer                                 :: FlowDistribution 
        real,    dimension(:    ), pointer      :: DistributionCoef
        integer, dimension(:    ), pointer      :: VectorI, VectorJ, VectorK
        real                                    :: AuxFlowIJ, SectionHeight, AuxFlowK
                
        !Sets to 0
        call SetMatrixValueAllocatable(Me%lFlowDischarge, Me%Size, 0.0, Me%ExtVar%WaterPoints3D)

        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMedia - ModifyWaterDischarges - ERR01'

        do iDis = 1, nDischarges

            call GetDischargeON(Me%ObjDischarges,iDis, IgnoreOK, STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'ModulePorousMedia - ModifyWaterDischarges - ERR010'

            if (IgnoreOK) cycle

            call GetDischargesGridLocalization(Me%ObjDischarges,                        &
                                               DischargeIDNumber = iDis,                &
                                               Igrid = i,                               &
                                               JGrid = j,                               &
                                               KGrid = kd,                              &
                                               DischVertical = DischVertical,           &  
                                               STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMedia - ModifyWaterDischarges - ERR020'
            
            !do not process runoff discharges (K=0). if uniform, K_layer is not used but instead k_min and K_max 
            !and user may forget K_layer zero
            if ((DischVertical == DischUniform_) .or. (kd /= 0)) then
                
                !real(8) to real as expected in GetDischargeWaterFlow
                SurfaceElevation = Me%ExtVar%Topography(i,j)
                call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                        Me%ExtVar%Now, iDis,                            &
                                        SurfaceElevation,                               &
                                        DischargeFlow, STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMedia - ModifyWaterDischarges - ERR030'


                call GetDischargeFlowDistribuiton(Me%ObjDischarges, iDis, nCells, FlowDistribution, &
                                                  VectorI, VectorJ, VectorK, kmin, kmax, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMedia - ModifyWaterDischarges - ERR040'

                !Horizontal distribution
i1:             if (nCells > 1) then
                    allocate(DistributionCoef(1:nCells))
i2:                 if      (FlowDistribution == DischByCell_ ) then
                    
                        DistributionCoef(1:nCells) = 1./float(nCells)

                    else i2
                    
                        stop 'ModulePorousMedia - ModifyWaterDischarges - ERR050'

                    endif i2
                endif i1
                
                AuxFlowIJ = DischargeFlow
                
 dn:            do n=1, nCells
 
                    if (nCells > 1) then
                        i         = VectorI(n)
                        j         = VectorJ(n)
                        kd        = VectorK(n)
                        
                        !For every cell get the total flow and multiply it by distribution coef
                        call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                                Me%ExtVar%Now, iDis,                            &
                                                SurfaceElevation,                               &
                                                AuxFlowIJ,                                      &
                                                FlowDistribution  = DistributionCoef(n),        &
                                                STAT = STAT_CALL)
                        if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMedia - ModifyWaterDischarges - ERR070'


                    endif

                    if (DischVertical == DischUniform_) then

                        if (kmin == FillValueInt) kmin = Me%ExtVar%KFloor(i, j)
                        if (kmax == FillValueInt) kmax = Me%WorkSize%KUB
                        SectionHeight = 0                                                
                        
                        do k=kmin, kmax                            
                            SectionHeight = SectionHeight + Me%ExtVar%DWZ(i, j, k)                        
                        enddo
                    else
            
                        kmin = kd; kmax = kd

                    endif

dk:                 do k=kmin, kmax

                        if (Me%ExtVar%WaterPoints3D(i, j, k) /= WaterPoint)  Cycle

                        if (DischVertical == DischUniform_) then
                        
                            AuxFlowK = Me%ExtVar%DWZ(i, j, k) / SectionHeight * AuxFlowIJ

                        else

                            AuxFlowK = AuxFlowIJ

                        endif
 
                        !each additional flow can remove all water left
                        if (AuxFlowK .lt. 0.0) then
                            !m3/s = m3 /s
                            MaxFlow = - ((Me%Theta(i,j,k) - Me%RC%ThetaR(i, j,k)) *            &
                                          Me%ExtVar%CellVolume(i,j,k)) / Me%CV%CurrentDT
                          
                            if (abs(AuxFlowK) .gt. abs(MaxFlow)) then
                                AuxFlowK = MaxFlow 
                            endif
                        endif

                        Me%lFlowDischarge(i, j, k)  = Me%lFlowDischarge(i, j, k) + AuxFlowK

                        Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +          &
                                             Me%lFlowDischarge(i, j, k) * Me%CV%CurrentDT) /              &
                                             Me%ExtVar%CellVolume(i, j, k)


 
                    enddo dk

                enddo dn

                if (nCells > 1) deallocate(DistributionCoef)

                call UnGetDischarges(Me%ObjDischarges, VectorI, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModulePorousMedia - ModifyWaterDischarges - ERR070'

                call UnGetDischarges(Me%ObjDischarges, VectorJ, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModulePorousMedia - ModifyWaterDischarges - ERR080'

                call UnGetDischarges(Me%ObjDischarges, VectorK, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModulePorousMedia - ModifyWaterDischarges - ERR090'


                !if (Me%CheckMass) Me%TotalInputVolume = Me%TotalInputVolume + Me%DischargesFlow(iDis) * LocalDT
 
            endif
           
        enddo

    end subroutine ModifyWaterDischarges  
    
    !--------------------------------------------------------------------------
   
    subroutine InsertInfiltrationOnFluxMatrix
    
        !Local-----------------------------------------------------------------
        integer :: i, j, k
        real    :: infiltration_flux

        !----------------------------------------------------------------------
        
        k = Me%WorkSize%KUB + 1
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB 
        
            if (Me%ExtVar%BasinPoints (i, j) == 1) then
        
                infiltration_flux = Me%Infiltration(i, j) * Me%ExtVar%Area(i, j) / Me%ExtVar%DT
                Me%FluxWAccFinal(i, j, k) = Me%FluxWAccFinal(i, j, k) - Infiltration_flux
                
            endif
        
        enddo
        enddo           
        
        !----------------------------------------------------------------------
    
    end subroutine InsertInfiltrationOnFluxMatrix
    
    !--------------------------------------------------------------------------
    
    subroutine CalculateMeanFlows (iterations)
    
        !Arguments-------------------------------------------------------------
        integer :: iterations

        !Local-----------------------------------------------------------------
        integer :: i, j, k

        !----------------------------------------------------------------------
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB            
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        
            if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                
               Me%FluxUAcc(i, j, k) = Me%FluxUAcc(i, j, k) / iterations                 
            
            endif
            
            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then
            
                Me%FluxVAcc(i, j, k) = Me%FluxVAcc(i, j, k) / iterations   
            
            endif

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
            
                Me%FluxWAcc(i, j, k) = Me%FluxWAcc(i, j, k) / iterations
                Me%FluxWAccFinal(i, j, k) = Me%FluxWAccFinal(i, j, k) / iterations
            
            endif
            
        enddo
        enddo
        enddo
        
!        Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) = Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) / iterations
             
    
    end subroutine CalculateMeanFlows
 
    !--------------------------------------------------------------------------
    
    subroutine AccumulateFlows

        !Local-----------------------------------------------------------------
        integer :: i, j, k

        !----------------------------------------------------------------------
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB            
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        
            if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                
               Me%FluxUAcc(i, j, k) = Me%FluxUAcc(i, j, k) + Me%FluxU(i, j, k)                 
            
            endif
            
            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then
            
                Me%FluxVAcc(i, j, k) = Me%FluxVAcc(i, j, k) + Me%FluxV(i, j, k)    
            
            endif

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
            
                Me%FluxWAcc(i, j, k) = Me%FluxWAcc(i, j, k) + Me%FluxW(i, j, k)
                Me%FluxWAccFinal(i, j, k) = Me%FluxWAccFinal(i, j, k) + Me%FluxWFinal(i, j, k)
            
            endif
            
        enddo
        enddo
        enddo
        
        Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) = Me%FluxWAcc(i, j, Me%WorkSize%KUB+1) + Me%FluxW(i, j, Me%WorkSize%KUB+1)
             
        !----------------------------------------------------------------------
        
    end subroutine AccumulateFlows
    
    !--------------------------------------------------------------------------

    subroutine EffectiveVelocity

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
    
    
    end subroutine EffectiveVelocity

    !--------------------------------------------------------------------------

    subroutine ComputeNextDT(Niter)

        !Arguments-------------------------------------------------------------
        integer                                     :: Niter

        !Local-----------------------------------------------------------------
        real                                        :: CurrentDT, MaxDT
        integer                                     :: STAT_CALL
        logical                                     :: VariableDT

!        Me%NextDT = Me%ExtVar%DT
!        
!        if (Niter <= Me%CV%MinIterations) then
!            Me%NextDT = Me%NextDT * Me%CV%DTFactorUp
!        else if (Niter > Me%CV%MaxIter) then
!            if (Me%CV%Stabilize .and. (Me%NextNIteration >= Me%CV%StabilizeHardCutLimit)) then
!                Me%NextDT = Me%CV%CurrentDT
!                Me%NextNIteration = 1
!            elseif (Me%NextNIteration > Me%LastGoodNiteration) then            
!                Me%NextDT = Me%NextDT * Me%CV%DTFactorDown
!            else
!                Me%NextDT = Me%NextDT * (Me%CV%DTFactorUp * 0.5)
!            endif
!        endif
!        
!        if (Me%CV%CurrentDT > Me%NextDT) then
!            Me%NextDT = Me%CV%CurrentDT
!            Me%NextNIteration = Niter
!        endif
!        
!        Me%LastGoodNiteration = Niter

        call GetVariableDT(Me%ObjTime, VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModulePorousMedia -  ERR010'

        call GetMaxComputeTimeStep(Me%ObjTime, MaxDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNextDT - ModulePorousMedia -  ERR020'
        
        if (VariableDT) then

            if (Niter == 1) then
            
                Me%CV%NextDT         = Me%ExtVar%DT * Me%CV%DTFactorUp
                Me%CV%NextNIteration = Niter
                
            elseif (Niter <= Me%CV%MinIterations) then                            
            
                if (Niter > Me%CV%LastGoodNiteration) then

                    Me%CV%NextDT         = Me%ExtVar%DT
                    Me%CV%NextNIteration = Niter

                else
                
                    Me%CV%NextDT         = Me%ExtVar%DT * Me%CV%DTFactorUp
                    Me%CV%NextNIteration = Niter

                endif
                
            else
            
                if (Niter >= Me%CV%StabilizeHardCutLimit) then
                
                    Me%CV%NextDT         = (Me%ExtVar%DT / Niter) * Me%CV%MinIterations
                    Me%CV%NextNIteration = Me%CV%MinIterations
                    
                elseif (Niter > Me%CV%LastGoodNiteration) then
                
                    Me%CV%NextDT = Me%ExtVar%DT / Me%CV%DTFactorDown
                    Me%CV%NextNIteration = max(int(Me%CV%NextDT / Me%CV%CurrentDT), 1)
                    
                else
                
                    Me%CV%NextDT = Me%ExtVar%DT
                    Me%CV%NextNIteration = max(min(int(Niter / Me%CV%DTSplitFactor), Niter - 1), 1)
                    
                endif 
                               
            endif
            
            CurrentDT = Me%CV%NextDT / Me%CV%NextNIteration                                                                       
            
            if (MaxDT < Me%CV%NextDT) then 
                Me%CV%NextDT         = MaxDT
                Me%CV%NextNIteration = max(int(Me%CV%NextDT/CurrentDT), 1)
            endif
                       
        else
        
            Me%CV%NextDT         = Me%ExtVar%DT
            Me%CV%NextNIteration = Niter            
        
        endif
        
        Me%CV%LastGoodNiteration = Niter
        Me%CV%CurrentDT          = Me%CV%NextDT / Me%CV%NextNIteration
   
    end subroutine ComputeNextDT 

    !----------------------------------------------------------------------------

    subroutine SoilWaterVelocity

        !Arguments---------------------------------------------------------------
                  
        !Local-------------------------------------------------------------------
        integer                                     :: I, J, K                
        integer                                     :: Chunk
        !------------------------------------------------------------------------          

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "SoilWaterVelocity")

!        call FinalHead

        !$OMP PARALLEL PRIVATE(I,J,K)

        !Horizontal Velocities
        if (Me%SoilOpt%CalcHorizontal) then
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB            

                if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                    Me%UnsatVelV(I,J,K) = BuckinghamDarcyEquation               &
                                       (con      = Me%UnsatK_Y(I,  J,K),        &
                                        hinf     = Me%FinalHead (I-1,J,K),      &
                                        hsup     = Me%FinalHead (I,  J,K),      &
                                        delta    = Me%ExtVar%DZY (I-1,J  ))
                end if

                if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                    
                    Me%UnsatVelU(I,J,K) = BuckinghamDarcyEquation               &
                                       (con      = Me%UnsatK_X(I,J,K  ),        &
                                        hinf     = Me%FinalHead (I,J-1,K),      &
                                        hsup     = Me%FinalHead (I,J,  K),      &
                                        delta    = Me%ExtVar%DZX (I,J-1  ))
                end if
            end do
            end do
            end do
            !$OMP END DO NOWAIT
        endif
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB            
        
            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
            
                Me%UnsatVelW(I,J,K) = BuckinghamDarcyEquation           &
                                   (con      = Me%UnsatK_Z(I,J,K  ),    &
                                    hinf     = Me%FinalHead (I,J,K-1),  &
                                    hsup     = Me%FinalHead (I,J,K  ),  &
                                    delta    = Me%ExtVar%DZZ (I,J,K-1))
                Me%UnsatVelWFinal(I,J,K) = Me%UnsatVelW(I,J,K)
            else

                Me%UnsatVelW(I,J,K) = 0.0
                Me%UnsatVelWFinal(I,J,K) = Me%UnsatVelW(I,J,K)
            end if
        end do
        end do
        end do
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        call InfiltrationVelocity


        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "SoilWaterVelocity")

    end subroutine SoilWaterVelocity

    !--------------------------------------------------------------------------
    
    subroutine ComputeFinalHead

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: CHUNK
        real                                        :: AccumPressure, Coef
        real                                        :: CenterVelocityW, ThetaInterpol
        integer                                     :: KUB
        !Begin-----------------------------------------------------------------
        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        KUB   = Me%WorkSize%KUB

        !$OMP PARALLEL PRIVATE(I,J,K,AccumPressure,Coef,ThetaInterpol,CenterVelocityW)
        !$OMP DO SCHEDULE    (DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints (i, j) == 1) then
            
                !Initial pressure = watercolumn - move downwards
                AccumPressure = Me%WaterColumn(i, j)
                
                do k = Me%WorkSize%KUB, Me%ExtVar%KFloor(i, j), -1
                    
                    !Evaluate if cell is almost saturated (for stability reasons a range of theta has to be given)
                    if (Me%Theta(i, j, k) .gt. Me%RC%ThetaS(i, j, k) * Me%CV%ThetaHydroCoef) then
                        
                        !Vertical velocity at the center of the cell and from previous step 
                        !without sub VerticalContinuity or sub ExchangeWithDrainageNetwork (layered flow) changes (UnsatVelWFinal)
                        CenterVelocityW  = (Me%UnsatVelW(i, j, k) + Me%UnsatVelW(i, j, k+1)) / 2.0
                        
                        !CenterVelocityU         = (Me%UnsatVelU(i, j, k) + Me%UnsatVelU(i, j+1, k)) / 2.0
                        !CenterVelocityV         = (Me%UnsatVelV(i, j, k) + Me%UnsatVelV(i+1, j, k)) / 2.0
                        !CenterVelocityW         = (Me%UnsatVelW(i, j, k) + Me%UnsatVelW(i, j, k+1)) / 2.0
                        !Modulus                 = sqrt(CenterVelocityU**2.0 + CenterVelocityV**2.0 + CenterVelocityW**2.0)
                        
                        !Compute hydro pressure if the velocity is upwards or if is downwards and lower than saturated conductivity
                        if ((CenterVelocityW .gt. 0.0)                                                               &
                             .or. (CenterVelocityW .le. 0.0 .and. CenterVelocityW .ge. -1.0 * Me%SatK(i, j, k))) then
                            
                            !Interpolation between 0 and 1 for theta (because hydro pressure evaluation starts
                            !at thetaS * HydroCoef and to give a proportional factor). Exponential for stability reasons
                            ThetaInterpol = (LinearInterpolation(Me%RC%ThetaS(i, j, k) * Me%CV%ThetaHydroCoef, 0.0,  &
                                                                Me%RC%ThetaS(i, j, k), 1.0, Me%Theta(i, j, k)))**3

                            !Hydro pressure coef dependent on velocity (remind, upward or downward and lower eq than SatK)
                            Coef =  (1.0 + CenterVelocityW / Me%SatK(i, j, k)) * ThetaInterpol
                            
                        else
                            !If velocity downwards and bigger than SatK there will be no hydro pressure 
                            Coef = 0.0
                        endif
                        
                        AccumPressure = AccumPressure + 0.5 * Me%ExtVar%DWZ(i, j, k) * Coef
                        Me%HydroPressure(i,j,k) = AccumPressure
                        AccumPressure = AccumPressure + 0.5 * Me%ExtVar%DWZ(i, j, k) * Coef
                        
                    else

                        Me%HydroPressure(i,j,k) = 0.0
                        AccumPressure = 0.0
                    
                    endif
                
                enddo
                
                !Final Head = Topography + Hydropressure + Soil Suction
                do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB
                    Me%FinalHead(i, j, k) = Me%ExtVar%CenterCell(i, j, k) + Me%HydroPressure(i,j,k) + Me%Head(i, j, k)
                enddo

            endif
            
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine ComputeFinalHead

    !----------------------------------------------------------------------------
       
    subroutine InfiltrationVelocity
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, KUB
        real                                        :: hsup_aux, DeltaSup_aux, hinf_aux
        real                                        :: Conductivity
        !Begin-----------------------------------------------------------------
        
        KUB = Me%WorkSize%KUB
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints (i, j) == 1) then
                
                !dz in surface interface
                DeltaSup_aux = Me%ExtVar%DWZ(i, j, KUB)/2.0 + Me%WaterColumn(i, j)
                
                !head in surface
                hsup_aux = Me%ExtVar%Topography(i, j) + Me%WaterColumn(i, j)
                
                !head in top cell - hydropressure had to be removed in order to accurate 
                !estimate infiltration velocity because hydropressure is computed for stability reasons 
                !even when cell is not saturated (satK * hydroCoef) - hydropressure and suction exist at 
                !the same time increasing Final Head in top cell and decreasing infil. vel.
                hinf_aux = Me%FinalHead(i, j, KUB) - Me%HydroPressure(i,j,KUB)
                
                !for tests only - by default nothing changes - saturation conduct
                if (Me%SoilOpt%InfiltrationConductivity == UnSatCond_) then
                    Conductivity = Me%UnSatK (i, j, KUB)
                elseif (Me%SoilOpt%InfiltrationConductivity == SatCond_) then
                    Conductivity = Me%SatK (i, j, KUB)
                endif
                
                Me%UnsatVelW(i, j, KUB+1) = BuckinghamDarcyEquation                             &
!                                               (con      = Me%SatK (i, j, KUB),                 &
                                                (con      = Conductivity,                       &
!                                                hinf     = Me%FinalHead(i, j, KUB),            &
                                                hinf     = hinf_aux,                            &
                                                hsup     = hsup_aux,                            &
!                                                delta    = Me%ExtVar%DWZ(i, j, KUB)/2.0)
                                                delta    = DeltaSup_aux)
 
                Me%UnsatVelWFinal(I,J,KUB+1) = Me%UnsatVelW(I,J,KUB+1)

                   
                if (-1.0 * Me%UnsatVelW(i, j, KUB+1) > Me%WaterColumn(i, j) / Me%CV%CurrentDT) then  
                    
                    Me%UnsatVelW(i, j, KUB+1)  = -1.0 * Me%WaterColumn(i, j) / Me%CV%CurrentDT

                    Me%UnsatVelWFinal(I,J,KUB+1) = Me%UnsatVelW(I,J,KUB+1)
               
                endif
                

            endif
            

        enddo
        enddo
    
    end subroutine InfiltrationVelocity

    !----------------------------------------------------------------------------

    subroutine SoilWaterFlux

        !Arguments---------------------------------------------------------------
                  
        !Local-------------------------------------------------------------------
        integer                                     :: I, J, K                
        integer                                     :: Chunk

        !------------------------------------------------------------------------          

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        )
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "SoilWaterFlux")

        !$OMP PARALLEL PRIVATE(I,J,K)

        !Horizontal Velocities
        if (Me%SoilOpt%CalcHorizontal) then
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB            

                if (Me%ExtVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then
                    Me%FluxU    (I,J,K) = Me%UnsatVelU(I,J,K) * Me%ExtVar%AreaU(I,J,K)
                end if

                if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then
                    Me%FluxV    (I,J,K) = Me%UnsatVelV(I,J,K) * Me%ExtVar%AreaV(I,J,K)
                end if

            end do
            end do
            end do
            !$OMP END DO NOWAIT
        endif
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then                                       
                    Me%FluxW(i,j,k)      = Me%UnsatVelW(i,j,k) * Me%ExtVar%Area(i,j)
                    Me%FluxWFinal(i,j,k) = Me%FluxW(i,j,k)
            end if
            
        end do
        end do
        end do
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "SoilWaterFlux")

    end subroutine SoilWaterFlux
    
    !--------------------------------------------------------------------------
    
    subroutine Condutivity_Face                          

        !Local-----------------------------------------------------------------                        

        !----------------------------------------------------------------------        

        select case (Me%SoilOpt%CondutivityFace)

        case (Average )

            call CondutivityAverage

        case (Maximum )

            call CondutivityMaximum            

        case (Minimum )

            call CondutivityMinimum

        case (Weighted)

            write (*,*)'Not Implemented'
            stop 'Condutivity_Face - ModulePorousMedia - ERR01'

        case (GeometricAvg)

            call CondutivityGeometricAverage

        end select

    end subroutine Condutivity_Face

    !--------------------------------------------------------------------------

    subroutine CondutivityAverage

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        )
        
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = (Me%UnSatK(i,j-1,k) * Me%ExtVar%DUX(i,j  )       +  &
                                         Me%UnSatK(i,j  ,k) * Me%ExtVar%DUX(i,j-1))      /  &
                                        (Me%ExtVar%DUX (i,j-1  ) + Me%ExtVar%DUX(i,j  ))      *  &
                                         Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = (Me%UnSatK(i-1,j,k) * Me%ExtVar%DVY(i  ,j)       +  &
                                         Me%UnSatK(i,j  ,k) * Me%ExtVar%DVY(i-1,j))      /  &
                                        (Me%ExtVar%DVY (i-1,j  ) + Me%ExtVar%DVY(i,j  ))      *  &
                                         Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then

                Me%UnsatK_Z(i, j, k) = (Me%UnSatK(i,j,k-1) * Me%ExtVar%DWZ(i,j,k  )      +  &
                                        Me%UnSatK(i,j,k)   * Me%ExtVar%DWZ(i,j,k-1))     /  &
                                       (Me%ExtVar%DWZ  (i,j,k-1) + Me%ExtVar%DWZ(i,j  ,k))

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityAverage

    !--------------------------------------------------------------------------

    subroutine CondutivityMaximum

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        )
        
        !X
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = max (Me%UnSatK(i,j-1,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT
        
        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = max (Me%UnSatK(i-1,j,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then


                Me%UnsatK_Z(i, j, k) = max (Me%UnSatK(i,j,k-1), Me%UnSatK(i,j,k))

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityMaximum

    !--------------------------------------------------------------------------

    subroutine CondutivityMinimum

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        )
        
        !X
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = min (Me%UnSatK(i,j-1,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = min (Me%UnSatK(i-1,j,k), Me%UnSatK(i,j,k)) *  &
                                                Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then
                Me%UnsatK_Z(i, j, k) = min (Me%UnSatK(i,j,k-1), Me%UnSatK(i,j,k))
            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityMinimum

    !--------------------------------------------------------------------------

    subroutine CondutivityGeometricAverage

        !Local-----------------------------------------------------------------                        
        integer                                     :: i, j, k
        integer                                     :: Chunk        

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        )
        

        !X
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) .EQ. Compute) then
                
                 Me%UnsatK_X(i, j, k) = (Me%UnSatK(i,j-1,k) ** Me%ExtVar%DUX(i,j  )     *  &
                                            Me%UnSatK(i,j  ,k) ** Me%ExtVar%DUX(i,j-1))    **  &
                                            (1.0/(Me%ExtVar%DUX (i,j-1) + Me%ExtVar%DUX(i,j)))   *  &
                                            Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                 Me%UnsatK_Y(i, j, k) = (Me%UnSatK(i-1,j,k) ** Me%ExtVar%DVY(i  ,j)       *   &
                                            Me%UnSatK(i,j  ,k) ** Me%ExtVar%DVY(i-1,j))      **  &
                                            (1.0/(Me%ExtVar%DVY (i-1,j) + Me%ExtVar%DVY(i,j)))     *   &
                                            Me%SoilOpt%HCondFactor

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !Z
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then

                Me%UnsatK_Z(i, j, k) = (Me%UnSatK(i,j,k-1) ** Me%ExtVar%DWZ(i,j,k  )      *  &
                                           Me%UnSatK(i,j,k)   ** Me%ExtVar%DWZ(i,j,k-1))     ** &
                                          (1.0/(Me%ExtVar%DWZ  (i,j,k-1) + Me%ExtVar%DWZ(i,j  ,k)))

            endif

        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

    end subroutine CondutivityGeometricAverage

    !--------------------------------------------------------------------------

    subroutine EvaporationFlux()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------                 
        integer                                     :: i, j, k 
        real                                        :: WaterVolume
        real                                        :: VelocityVolume
        real(8)                                     :: EvapoVolume, SoilVolume
        real                                        ::  NewTheta 
        real                                        :: TotalCol, HeadLimit
        !Begin-----------------------------------------------------------------

        !Set EvapFlux to zero
        Me%EvaporationFlux(:,:) = 0.0

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == 1) then
        
                if ((Me%SoilOpt%IgnoreWaterColumnOnEvap) .or. (Me%WaterColumn(i, j) < AllmostZero)) then             
                   
                    k = Me%WorkSize%KUB   ! evaporation just at the surface
                    
                    !m = m/s * s
                    TotalCol = Me%ExtVar%PotentialEvaporationFlux(i, j) * Me%CV%CurrentDT    ! available water for evaporation
                    ! m3               
                    WaterVolume = TotalCol * Me%ExtVar%Area (i,j)

                    !Velocity Volume
                    if (Me%SoilOpt%LimitEVAPWaterVelocity) then
                        VelocityVolume = Me%UnSatK (i,j,k) * Me%CV%CurrentDT * Me%ExtVar%Area (i,j)
                        EvapoVolume    = min(WaterVolume, VelocityVolume)
                    else
                        EvapoVolume    = WaterVolume
                    endif
                            
                    !Avaliable Soil Water volume in layer
                    SoilVolume  = (Me%Theta(i,j,k) - Me%RC%ThetaR(i,j,k)) * Me%ExtVar%CellVolume(i,j,k) 
                    EvapoVolume = min(EvapoVolume, SoilVolume)
                            
                  
                    !Estimates new Theta
                    NewTheta = Me%Theta(i,j,k) - EvapoVolume/ Me%ExtVar%CellVolume(i,j,k)
                            
                    if (Me%SoilOpt%LimitEVAPHead) then     
                        HeadLimit= Me%Soilopt%HeadLimit
                       
                        if (NewTheta > Theta_(HeadLimit, Me%SoilID(i,j,k)))   then        
                                
                            !Evaporation Flux
                            Me%EvaporationFlux(i,j) = EvapoVolume / Me%CV%CurrentDT
                        else 
                    
                            Me%EvaporationFlux(i,j) = 0.0

                        endif
                    else
                        
                        !Just uses new theta if not all dry... stability reasons
                        if (NewTheta > Me%RC%ThetaR(i,j,k) + Me%CV%LimitThetaLo) then
                        
                            !Evaporation Flux
                           Me%EvaporationFlux(i,j) = EvapoVolume / Me%CV%CurrentDT

                        else 
                    
                            Me%EvaporationFlux(i,j) = 0.0

                        endif 
                    endif           
                
                endif
            
            endif
            
        enddo
        enddo            


    end subroutine EvaporationFlux

    !--------------------------------------------------------------------------

    subroutine CalculateNewTheta()

        !Arguments-------------------------------------------------------------  
!        logical                             :: OverSaturation
        !Local-----------------------------------------------------------------        
        integer                             :: CHUNK, I, J, K
        
        
        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        )
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "CalculateNewTheta")

        !$OMP PARALLEL PRIVATE(I,J,K)

        !Horizontal Fluxes
        if (Me%SoilOpt%CalcHorizontal) then

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                    
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +        &
                                        (Me%FluxU(i,j,k)  * Me%ExtVar%ComputeFacesU3D(i,j,k) -      &
                                         Me%FluxU(i,j+1,k)* Me%ExtVar%ComputeFacesU3D(i,j+1,k)) *   &
                                         Me%CV%CurrentDT) / Me%ExtVar%CellVolume(i, j, k)
                    
                endif
                
            enddo
            enddo
            enddo
            !$OMP END DO
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +        &
                                        (Me%FluxV(i,j,k)  * Me%ExtVar%ComputeFacesV3D(i,j,k) -      &
                                         Me%FluxV(i+1,j,k)* Me%ExtVar%ComputeFacesV3D(i+1,j,k)) *   &
                                         Me%CV%CurrentDT) / Me%ExtVar%CellVolume(i, j, k)
                    
                endif
                
            enddo
            enddo
            enddo
            !$OMP END DO

        endif
               
        !Flux W may be corrected in river points on the saturated zone to compensate exchange with river
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
            
                Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) +        &
                                     (Me%FluxWFinal(i,j,k)  * Me%ExtVar%ComputeFacesW3D(i,j,k) -     &
                                      Me%FluxWFinal(i,j,k+1)* Me%ExtVar%ComputeFacesW3D(i,j,k+1)) *  &
                                      Me%CV%CurrentDT) / Me%ExtVar%CellVolume(i, j, k)
                              
            endif 
            
        enddo
        enddo
        enddo
        !$OMP END DO

        
        !Transpiration        
        if (Me%TranspirationExists) then

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -        & 
                                         Me%ExtVar%TranspirationFlux(i, j, k) * Me%CV%CurrentDT)  / &
                                         Me%ExtVar%CellVolume(i,j,k)
                  
                endif
                
            enddo
            enddo
            enddo
            !$OMP END DO
            
        endif

        if (Me%EvaporationExists) then

            k = Me%WorkSize%KUB

            !Evaporation
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(I,J) == WaterPoint) then
                    Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -       & 
                                         Me%EvaporationFlux(i, j) * Me%CV%CurrentDT)  /            &
                                         Me%ExtVar%CellVolume(i,j,k)
                endif
            enddo
            enddo
            !$OMP END DO
        endif
        
        !Infiltration
        k = Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%BasinPoints(i,j) == WaterPoint) then
            
                Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -                &
                                     Me%UnsatVelWFinal(i, j, k+1) * Me%ExtVar%Area(i, j) *              &
                                     (1.0 - Me%ImpermeableFraction(i, j)) * Me%CV%CurrentDT) /          &
                                     Me%ExtVar%CellVolume(i, j, k)
                
            endif
            
        enddo
        enddo
        !$OMP END DO
        
        if (Me%ObjDrainageNetwork /= 0) then
        
            if (Me%SoilOpt%DNLink == GWFlowToChanByCell_) then
            
                !Exchange with River
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    if (Me%ExtVar%RiverPoints(i,j) == WaterPoint) then
                    
                        k = Me%UGCell(i,j)
                        Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -                &
                                             Me%lFlowToChannels(i, j) * Me%CV%CurrentDT) /                      &
                                             Me%ExtVar%CellVolume(i, j, k)
                       
                    endif
                    
                enddo
                enddo
                !$OMP END DO
            
            elseif (Me%SoilOpt%DNLink == GWFlowToChanByLayer_) then

                !Exchange with River
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    if (Me%ExtVar%RiverPoints(i,j) == WaterPoint) then
                    
                        do k = Me%FlowToChannelsBottomLayer(i,j), Me%FlowToChannelsTopLayer(i,j)
                            Me%Theta(i, j, k) = (Me%Theta(i, j, k) * Me%ExtVar%CellVolume(i, j, k) -                &
                                                 Me%lFlowToChannelsLayer(i, j, k) * Me%CV%CurrentDT) /              &
                                                 Me%ExtVar%CellVolume(i, j, k)
                                                    
                        enddo
                        
                    endif
                    
                enddo
                enddo
                !$OMP END DO
                
            endif
            
        endif
        
        !$OMP END PARALLEL        
        
        !Discharges - modify specific cells - more eficient than do's cicles fo all cells 
        if (Me%SoilOpt%Discharges) then

            call ModifyWaterDischarges
            
        endif        
        
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "CalculateNewTheta")
        

    end subroutine CalculateNewTheta

    !--------------------------------------------------------------------------

    subroutine IntegrateValuesInTime(SumDT)
    
        !Arguments-------------------------------------------------------------  
        real                                :: SumDT

        !Local-----------------------------------------------------------------        

        !Updates water column and infiltration
        call UpdateWaterColumnInfiltration
        
        !Updates Efective EVTP
        call UpdateEfectiveEVTP
        
        !Integrates Flow
        if (Me%ObjDrainageNetwork /= 0) then
            call IntegrateFlow (SumDT)
        endif

        if (Me%SoilOpt%Discharges) then
            call IntegrateDischargeFlow (SumDT)
        endif

    end subroutine IntegrateValuesInTime

    !--------------------------------------------------------------------------

    subroutine IntegrateDischargeFlow (SumDT)

        !Arguments-------------------------------------------------------------
        real                                        :: SumDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: CHUNK
        real(8)                                     :: sum

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        sum = Me%TotalDischargeFlowVolume
        
        !Integrates Flow Discharges        
        !$OMP PARALLEL PRIVATE(I,J) REDUCTION(+:sum)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%iFlowDischarge(i, j, k) = (Me%iFlowDischarge(i, j, k) * SumDT + Me%lFlowDischarge(i, j, k) * Me%CV%CurrentDT) / &
                                         (SumDT + Me%CV%CurrentDT)
            sum = sum + Me%lFlowDischarge(i, j, k) * Me%CV%CurrentDT
        enddo
        enddo
        !$OMP END DO NOWAIT      
        enddo
        !$OMP END PARALLEL  
        
        Me%TotalDischargeFlowVolume = sum

    end subroutine IntegrateDischargeFlow

    !--------------------------------------------------------------------------

    subroutine UpdateWaterColumnInfiltration 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: chunk
        real                                        :: dh
     
        !Begin-----------------------------------------------------------------
        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J, dh)

        !Update Water Column
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB,     Me%WorkSize%JUB
        do I = Me%WorkSize%ILB,     Me%WorkSize%IUB
                                    
            if (Me%ExtVar%BasinPoints (i,j) == 1) then
                
                !Variation in height
!                dh                  = Me%UnsatVelW(i, j, Me%WorkSize%KUB+1) * (1.0 - Me%ImpermeableFraction(i, j)) * &
!                                      Me%CV%CurrentDT
                
                !m/s * m2AreaPerm * s = m3 H20 in flux
                !m3 H20 / m2AreaTotal = m uniformly distributed in area
                dh                  = -1.0 * Me%UnsatVelWFinal(i, j, Me%WorkSize%KUB+1) * (1.0 - Me%ImpermeableFraction(i, j)) * &
                                      Me%CV%CurrentDT
                
                !Just reduce water column due to Infiltration (PermeableFraction only)
                Me%WaterColumn(i,j) = Me%WaterColumn(i,j)  - dh
                
                Me%Infiltration(i,j)= Me%Infiltration(i,j) + dh
                
                if (abs(Me%WaterColumn(i, j)) < AllmostZero) Me%WaterColumn(i, j) = 0.0
                
                if (Me%WaterColumn(i, j) < -1.0e-10) then
                    write(*,*)'Bug Infiltration', i, j, Me%WaterColumn(i, j)
                endif
                
            endif
            
        enddo
        enddo            
        !$OMP END DO
        !$OMP END PARALLEL


    end subroutine UpdateWaterColumnInfiltration

    !--------------------------------------------------------------------------

    subroutine UpdateEfectiveEVTP

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: chunk
        real(8)                                     :: sum
     
        !Begin-----------------------------------------------------------------
        
        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)        )
        sum = Me%AccTranspiration
        
        if (Me%TranspirationExists) then
            !$OMP PARALLEL PRIVATE(I,J,K) REDUCTION(+:sum)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB,     Me%WorkSize%KUB
            do J = Me%WorkSize%JLB,     Me%WorkSize%JUB
            do I = Me%WorkSize%ILB,     Me%WorkSize%IUB
                                
                if (Me%ExtVar%Waterpoints3D (i,j,k) == 1) then
          
                    ! m = m + m3/s * s / m2
                    Me%EfectiveEVTP(i,j) = Me%EfectiveEVTP(i,j) + Me%ExtVar%TranspirationFlux(i, j, k) * Me%CV%CurrentDT/ &
                                           Me%ExtVar%Area(i, j)
                    
                    ! m3 = m3 + m3/s * s
                    sum = sum + Me%ExtVar%TranspirationFlux(i, j, k) * Me%CV%CurrentDT
                    
                endif

            enddo
            enddo            
            enddo  
            !$OMP END DO
            !$OMP END PARALLEL     
        
            Me%AccTranspiration = sum
            
        endif
        
        sum = Me%AccEvapFromSoil
        
        if (Me%EvaporationExists) then
            !$OMP PARALLEL PRIVATE(I,J,K) REDUCTION(+:sum)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB,     Me%WorkSize%JUB
            do I = Me%WorkSize%ILB,     Me%WorkSize%IUB
                                
                if (Me%ExtVar%BasinPoints (i,j) == 1) then
          
                    ! m = m + m3/s * s / m2
                    Me%EfectiveEVTP(i,j) = Me%EfectiveEVTP(i,j) + Me%EvaporationFlux(i, j) * Me%CV%CurrentDT / &
                                           Me%ExtVar%Area(i, j)
                    ! m = m + m3/s * s / m2
                    Me%EfectiveEVAP(i,j) = Me%EfectiveEVAP(i,j) + Me%EvaporationFlux(i, j) * Me%CV%CurrentDT / &
                                           Me%ExtVar%Area(i, j)
                    
                    ! m3 = m3 + m3/s * s
                    sum = sum + Me%EvaporationFlux(i, j) * Me%CV%CurrentDT
                endif
        
            enddo
            enddo            
              
            !$OMP END DO
            !$OMP END PARALLEL     
        
            Me%AccEvapFromSoil = sum
            
        endif
        

    end subroutine UpdateEfectiveEVTP

    !--------------------------------------------------------------------------
    
    subroutine IntegrateFlow(SumDT)
    
        !Arguments-------------------------------------------------------------
        real                                        :: SumDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: chunk
     
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%KLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,K)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%RiverPoints (i,j) == 1) then
                if (Me%SoilOpt%DNLink == GWFlowToChanByCell_) then

                    Me%iFlowToChannels(i, j) = (Me%iFlowToChannels(i, j) * SumDT +                                  &
                                                        Me%lFlowToChannels(i, j) * Me%CV%CurrentDT) /               &
                                                       (SumDT + Me%CV%CurrentDT)
                elseif (Me%SoilOpt%DNLink == GWFlowToChanByLayer_) then    
                    do k = Me%FlowToChannelsBottomLayer(i,j), Me%FlowToChannelsTopLayer(i,j)
                        Me%iFlowToChannelsLayer(i, j, k) = (Me%iFlowToChannelsLayer(i, j, k) * SumDT +                  &
                                                            Me%lFlowToChannelsLayer(i, j, k) * Me%CV%CurrentDT) /       &
                                                           (SumDT + Me%CV%CurrentDT)
                        Me%iFlowToChannels(i, j) = (Me%iFlowToChannels(i, j) * SumDT +                                  &
                                                            Me%lFlowToChannels(i, j) * Me%CV%CurrentDT) /               &
                                                           (SumDT + Me%CV%CurrentDT)
                    enddo
                endif 
            endif
            
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine IntegrateFlow

    !--------------------------------------------------------------------------

    subroutine SoilParameters

        !Arguments-------------------------------------------------------------        

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: Chunk

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "SoilParameters")
                
        !$OMP PARALLEL PRIVATE(I,J,K)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

cd1 :   if (Me%ExtVar%BasinPoints(i, j) == 1) then

            do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB

            !Verifies if cell is saturated (or very close). If so set it to
            !saturation
            !Otherwise calculate ThetaF Value                           !0.00001
            if (abs(Me%Theta(i,j,k) - Me%RC%ThetaS(i,j,k)) < Me%CV%LimitThetaHi) then

                Me%Theta      (i, j, k)         = Me%RC%ThetaS(i,j,k)
                Me%Head       (i, j, k)         = 0.0 
                Me%RC%ThetaF  (i, j, k)         = 1.0
                Me%UnSatK     (i, j, k)         = Me%SatK(i, j, k)
                Me%CalculateHead(i, j, k)       = .false.

            !Close or below residual value
            else if (Me%Theta(i,j,k) < Me%RC%ThetaR(i,j,k) + Me%CV%LimitThetaLo) then
            
                !This creates mass...
                Me%RC%ThetaF    (i, j, k) = Me%CV%LimitThetaLo / (Me%SoilTypes(Me%SoilID(I,J,K))%ThetaS   &
                                                                 - Me%SoilTypes(Me%SoilID(I,J,K))%ThetaR)
                Me%CalculateHead(i, j, k) = .true.
                
                call SetError(WARNING_, INTERNAL_, "Mass Created, SoilParameters", OFF)

            !Normal Case
            else
                Me%RC%ThetaF  (i, j, k)   = ThetaF_ (Me%Theta (i, j, k), Me%SoilID(i, j, k))
                Me%CalculateHead(i, j, k) = .true.
            endif
            
            end do

        end if cd1
        
        end do
        end do
        !$OMP END DO


        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

cd2 :   if (Me%ExtVar%BasinPoints(i, j) == 1) then

            do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB

            if (Me%CalculateHead(i, j, k)) then

                !Over saturation
                if (Me%RC%ThetaF(i, j, k) > 1.0) then

                    Me%Head   (i, j, k) = 0.0
                    Me%UnSatK (i, j, k) = Me%SatK(i, j, k)

                !0 < Theta < 1
                else if (Me%RC%ThetaF(i, j, k) > 0.0 .and. Me%RC%ThetaF(i, j, k) < 1.0) then

                    Me%Head  (i, j, k) = Head_   (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))

                    Me%UnSatK(i, j, k) = UnsatK_ (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))



                !Theta <= 0
                else

                    Me%Head(I,J,K) = null_real

                endif
                
            endif

            end do

        end if cd2
        
        end do
        end do
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
            

        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "SoilParameters")

    end subroutine SoilParameters

    !--------------------------------------------------------------------------    
    
    subroutine CheckCompleteStupidSolution(IsStupid, CorrectTheta)

  
        !Arguments-------------------------------------------------------------
        logical, intent(OUT)                        :: IsStupid        
        logical, intent(IN)                         :: CorrectTheta

        !Local-----------------------------------------------------------------        
        integer                                     :: I, J, K
        integer                                     :: nStupid

        !----------------------------------------------------------------------               

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "CheckCompleteStupidSolution")

        nStupid = 0

        !Tests variations
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == 1) then
                if (Me%Theta(i, j, k) < Me%RC%ThetaR(i, j, k)) then
                    nStupid = nStupid + 1
                endif
            endif
        enddo
        enddo
        enddo    
        
        if (nStupid > 0) then
            IsStupid = .true.
        else
            IsStupid = .false.
        endif
        
        if (CorrectTheta .and. IsStupid) then
            IsStupid = .false.
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == 1) then
                    if (Me%Theta(i, j, k) < Me%RC%ThetaR(i, j, k)) then
                        Me%Theta(i, j, k) = Me%RC%ThetaR(i, j, k)
                        call SetError(WARNING_, INTERNAL_, "Mass Created, STUPID SOLUTION", OFF)
                    endif
                endif
            enddo
            enddo
            enddo    
        endif
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "CheckCompleteStupidSolution")
        
    end subroutine CheckCompleteStupidSolution        

    !--------------------------------------------------------------------------    

    subroutine CheckStability(Restart)
        
        !Arguments-------------------------------------------------------------
        logical                                     :: Restart        

        !Local-----------------------------------------------------------------        
        integer                                     :: I, J, K
        real                                        :: variation
        integer                                     :: n_restart

        !----------------------------------------------------------------------               

        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "CheckStability")
        
        Restart   = .false.
        n_restart = 0

        !Verifies negative volumes
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                do K = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB - 1
            
                    if (Me%Theta(I,J,K) < -1.0 * AllmostZero) then
                        Restart = .true.
                        exit do1
                    elseif (Me%Theta(I,J,K) < 0.0) then
                        Me%Theta(I,J,K) = 0.0
                    endif
                
                enddo
                
            endif
            
        enddo
        enddo do1

        if ((.not. Restart) .and. Me%CV%Stabilize) then

            !Tests variations
do2:        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
                if (Me%ExtVar%BasinPoints(I,J) == 1) then
                
                    !In the upper layer the module can't check the convergence criteria because it is imposing
                    !a boundary value... 
                    !If the water column isn't suficient, DT for infiltration velocity and DT to converge are cutted
                    !in the same way.... 
                    do K = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB - 1
                        if ((.not. Me%CV%CheckDecreaseOnly) .or. (Me%CV%ThetaOld(I,J,K) > Me%Theta(I,J,K))) then
                            if (Me%CV%ThetaOld(I,J,K) >= Me%CV%MinimumValueToStabilize * Me%RC%ThetaS(I,J,K)) then

                                variation = abs(Me%Theta(I,J,K) - Me%CV%ThetaOld(I,J,K)) / Me%CV%ThetaOld(I,J,K)                 
                    
                                if (variation > Me%CV%StabilizeFactor) then  
                                    !Debug routine - may be usefull for using in debug situation
                                    !call DebugStability (i,j,k,variation)
                                                      
                                    n_restart = n_restart + 1
                                endif
                            endif
                        endif
                    enddo
                endif
                
            enddo
            enddo do2    
        
            if (n_restart > Me%CV%MinToRestart) then
                Restart = .true.
            endif        
        
        endif
        
        if (Restart) then        
            Me%CV%NextNiteration = max(int(Me%CV%NextNiteration * Me%CV%DTSplitFactor), Me%CV%NextNiteration + 1)
                 
            if (Me%CV%NextNiteration >= Me%CV%MaxIterations) then
                 write(*,*)'Number of iterations above maximum: ', Me%CV%NextNiteration
                 stop 'CheckStability - ModulePosourMedia - ERR010'
            endif                          
        endif              
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "CheckStability")

    end subroutine CheckStability
    
    !--------------------------------------------------------------------------

    subroutine DebugStability(i,j,k, variation)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: I, J, K
        real                                        :: variation
        !Local-----------------------------------------------------------------
        character (Len = 5)                         :: str_i, str_j, str_k
        character (Len = 15)                        :: str_1, str_2, str_3
        character (len = StringLength)              :: string_to_be_written 
        
        write(str_i, '(i3)') i 
        write(str_j, '(i3)') j
        write(str_k, '(i3)') k
        write(str_1, '(ES10.3)') Me%CV%ThetaOld(I,J,K)  
        write(str_2, '(ES10.3)') Me%Theta(I,J,K)   
        write(str_3, '(ES10.3)') variation                            
        
        string_to_be_written = ' '//str_i//','//str_j//' '//str_k//' '//str_1//' '//str_2//' '//str_3
        
        call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)           
    
    
    end subroutine DebugStability
    
    !--------------------------------------------------------------------------

    subroutine VerticalContinuity
    
        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        integer                                     :: chunk
        real                                        :: ExcessVolume !, dh

        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,K, ExcessVolume)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == 1) then
            
                do k = Me%WorkSize%KUB, (Me%ExtVar%KFloor(i, j)+1), -1               
                    if (Me%Theta(i,j,k) .gt. Me%RC%ThetaS(i,j,k)) then
                        !If cell is oversaturated, set to saturation and put water downwards
                        ExcessVolume         = (Me%Theta(i,j,k) - Me%RC%ThetaS (i,j,k)) * Me%ExtVar%CellVolume(i,j,k)
                        Me%Theta(i,j,k-1)    = ((Me%Theta(i,j,k-1) * Me%ExtVar%CellVolume(i,j,k-1)) + ExcessVolume)   &
                                               / Me%ExtVar%CellVolume(i,j,k-1)
                        Me%FluxWFinal(i,j,k) = Me%FluxWFinal(i,j,k) + (-ExcessVolume / Me%CV%CurrentDT)
                        !m/s
                        Me%UnsatVelWFinal(i,j,k)  = Me%FluxWFinal(i,j,k) / Me%ExtVar%Area(i,j)
                        Me%Theta(i,j,k)      = Me%RC%ThetaS (i,j,k)
                    endif
                enddo
                
                !Invert process and put the rest upwards
                do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB-1
                    if (Me%Theta(i,j,k) .gt. Me%RC%ThetaS(i,j,k)) then
                        !If cell is oversaturated, set to saturation and put water upwards
                        ExcessVolume           = (Me%Theta(i,j,k) - Me%RC%ThetaS (i,j,k)) * Me%ExtVar%CellVolume(i,j,k)
                        Me%Theta(i,j,k+1)      = ((Me%Theta(i,j,k+1) * Me%ExtVar%CellVolume(i,j,k+1)) + ExcessVolume)   &
                                                 / Me%ExtVar%CellVolume(i,j,k+1)
                        Me%FluxWFinal(i,j,k+1) = Me%FluxWFinal(i,j,k+1) + (ExcessVolume / Me%CV%CurrentDT)  
                        !m/s 
                        Me%UnsatVelWFinal(i,j,k+1)  = Me%FluxWFinal(i,j,k+1) / Me%ExtVar%Area(i,j)
                        Me%Theta(i,j,k)        = Me%RC%ThetaS (i,j,k)
                    endif
                enddo
                
                !Put remaing volume to watercolumn
                k = Me%WorkSize%KUB
                if (Me%Theta(i,j,k) .gt. Me%RC%ThetaS(i,j,k)) then
                    
                    if (Me%ImpermeableFraction(i, j) == 1.0) then
                        write (*,*)
                        write (*,*) 'Error. Trying to expel water from soil - mandatory to avoid'
                        write (*,*) 'over saturation - but the cell is impermeable. In cell', i, j 
                    endif
                    
                    !If cell is oversaturated, set to saturation and put water on water column
                    ExcessVolume           = ((Me%Theta(i,j,k) - Me%RC%ThetaS (i,j,k)) * Me%ExtVar%CellVolume(i,j,k))
!                    Me%FluxWFinal(i,j,k+1) = Me%FluxWFinal(i,j,k+1) + (ExcessVolume / Me%CV%CurrentDT)
                    Me%FluxWFinal(i,j,k+1) = (Me%UnsatVelWFinal(i,j,k+1) * Me%ExtVar%Area(i,j)                        &
                                              * (1.0 - Me%ImpermeableFraction(i, j)))                                 &
                                              + (ExcessVolume / Me%CV%CurrentDT)
                    !m/s
                    Me%UnsatVelWFinal(i,j,k+1)  = Me%FluxWFinal(i,j,k+1) /                                            &
                                                  (Me%ExtVar%Area(i,j) * (1.0 - Me%ImpermeableFraction(i, j)))

                                                  
!                    dh                     = ExcessVolume  / Me%ExtVar%Area(i, j)
!                    Me%WaterColumn  (i,j)  = Me%WaterColumn(i,j) + dh
!                    Me%Infiltration (i,j)  = Me%Infiltration (i,j) - dh
                    Me%Theta(i,j,k)        = Me%RC%ThetaS (i,j,k)
                endif
            endif
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine VerticalContinuity

    !--------------------------------------------------------------------------

    subroutine CorrectFinalHead
    
        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        !integer                                     :: chunk
        real                                        :: ExcessFlow, TotalAreaForFlow
        real                                        :: Distance, ExcessPressure

        
!        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

!        !$OMP PARALLEL PRIVATE(I,J,K, ExcessVolume, dh)
!        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == 1) then
            
                do k = Me%WorkSize%KUB, (Me%ExtVar%KFloor(i, j)), -1               
                
                    !dQ = dv * Atotal = -k dP/dxi
                    !dP/dxi = 
                    
                    if (Me%Theta(i,j,k) .gt. Me%RC%ThetaS(i,j,k)) then
                        !m3/s
                        ExcessFlow         = (Me%Theta(i,j,k) - Me%RC%ThetaS (i,j,k)) * Me%ExtVar%CellVolume(i,j,k)   &
                                              / Me%CV%CurrentDT
                        !m2                                              
                        TotalAreaForFlow   = Me%ExtVar%AreaU(i,j,k  )   * Me%ExtVar%ComputeFacesU3D(i,j,k)              &
                                             + Me%ExtVar%AreaU(i,j+1,k) * Me%ExtVar%ComputeFacesU3D(i,j+1,k)            &
                                             + Me%ExtVar%AreaV(i,j,k  ) * Me%ExtVar%ComputeFacesV3D(i,j,k)              &
                                             + Me%ExtVar%AreaV(i+1,j,k) * Me%ExtVar%ComputeFacesV3D(i+1,j,k)            &
                                             + Me%ExtVar%Area(i,j     ) * Me%ExtVar%ComputeFacesW3D(i,j,k)              &
                                             + Me%ExtVar%Area(i,j     ) * Me%ExtVar%ComputeFacesW3D(i,j,k+1)                     
                        
                        !Distance for Pressure dissipation
                        Distance           = Me%ExtVar%CellVolume(i,j,k) / TotalAreaForFlow
                        
                        !m = (m3/s * m ) / (m/s * m2)
                        ExcessPressure     = (ExcessFlow * Distance) / (Me%SatK(i,j,k) * TotalAreaForFlow)
                        
                        Me%FinalHead(i,j,k) = Me%FinalHead(i,j,k) + ExcessPressure
                        
                    endif
                enddo
            endif
        enddo
        enddo
!        !$OMP END DO NOWAIT
!        !$OMP END PARALLEL

    end subroutine CorrectFinalHead

    !--------------------------------------------------------------------------
    
    subroutine LogDT (iteration)

        !Arguments-------------------------------------------------------------
        integer                                     :: iteration

        !Local-----------------------------------------------------------------
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second

        call ExtractDate(Me%ExtVar%Now, Year, Month, Day, Hour, Minute, Second)
        
        write (Me%Files%AsciiUnit, fmt=1000) Year, Month, Day, Hour, Minute, Second, &
                                             iteration, Me%CV%CurrentDT

        1000 format(f5.0, f5.0, f5.0, f5.0, f5.0, f12.5, i3, f12.5)

    end subroutine LogDT

    !--------------------------------------------------------------------------

    subroutine ModifyBoundaryFlux
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, di, dj, CHUNK !, Sum, SumFaces
        !logical                                     :: NearBoundary
        !real                                        :: AverageArea, AverageDist
        real                                        :: OldVolume, AreaZX, AreaZY, BoundaryFinalHead
        real                                        :: NewTheta, ConductivityFace, sum
        !Begin-----------------------------------------------------------------

        
        Me%AccBoundaryFlowVolume = 0.0
        sum = 0.0
        
        !Impose lateral aquifer levels
        if (Me%SoilOpt%ImposeBoundaryWalls) then

            call SetMatrixValueAllocatable (Me%iFlowBoundaryWalls, Me%Size, 0.0, Me%ExtVar%WaterPoints3D)


            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

            !$OMP PARALLEL PRIVATE(I,J,K,di,dj,BoundaryFinalHead,ConductivityFace,AreaZX,AreaZY,OldVolume,NewTheta) REDUCTION(+:sum)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                !Process cells collumns that are boundary and that have topography lower than maximum
                if (Me%Boundary%BoundaryCells(i,j) == BasinPoint     &
                     .and. Me%ExtVar%Topography(i,j) < Me%Boundary%MaxDtmForBoundary) then
                 
                    !Final Head is in all boundary cells the same as BoundaryValue (field capacity)
                    !Limit it to bottom or it may create high fluxes
                    BoundaryFinalHead = max(Me%Boundary%ImposedBoundaryLevel(i,j), Me%ExtVar%BottomTopoG(i, j))

    do1:            do k = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB
                        
                        !do not process in saturated cells below BoundaryValue and any above maximum level 
                        !so flux occurs only between areas where saturation occurs in boundary an not inside
                        !or not inside but in boundary
                        if (     (- Me%ExtVar%SZZ(i,j,k-1)                                                     &
                                   .ge. max(Me%Boundary%ImposedBoundaryLevel(i,j), Me%UGWaterLevel2D(i, j)))   &
                            .or. (- Me%ExtVar%SZZ(i,j,k)                                                       &
                                   .lt. min (Me%Boundary%ImposedBoundaryLevel(i,j), Me%UGWaterLevel2D(i, j)))) then
                            cycle do1
                        endif

                        !Conductivity is the maximum (saturation front entering or exiting domain)
                        ConductivityFace  = Me%SatK(i,j,k) * Me%SoilOpt%HCondFactor
                        
                        !U direction - use middle area because in closed faces does not exist AreaU
                        !if gradient positive, than flow negative (exiting soil)
                        AreaZX = Me%ExtVar%DVY(i,j) * Me%ExtVar%DWZ(i,j,k)
                        do dj = 0, 1
                            if ((Me%ExtVar%ComputeFacesU3D(i,j+dj,k) == 0)) then
                                Me%iFlowBoundaryWalls(i,j,k) = Me%iFlowBoundaryWalls(i,j,k) + AreaZX &
                                                          *     BuckinghamDarcyEquation   &
                                                  (con      = ConductivityFace,           &
                                                   hinf     = BoundaryFinalHead,          &
                                                   hsup     = Me%FinalHead (i  ,  j,k),   &
                                                   delta    = Me%ExtVar%DUX(i,j)       )
                            endif
                        enddo

                        !V direction - use middle area because in closed faces does not exist AreaV
                        AreaZY = Me%ExtVar%DUX(i,j) * Me%ExtVar%DWZ(i,j,k)                      
                        do di = 0, 1
                            if ((Me%ExtVar%ComputeFacesV3D(i+di,j,k) == 0)) then
                                Me%iFlowBoundaryWalls(i,j,k) = Me%iFlowBoundaryWalls(i,j,k) + AreaZY &
                                                          *     BuckinghamDarcyEquation   &
                                                  (con      = ConductivityFace,           &
                                                   hinf     = BoundaryFinalHead,          &
                                                   hsup     = Me%FinalHead (i  ,  j,k),   &
                                                   delta    = Me%ExtVar%DVY(i,j)       )
                            endif
                        enddo
                                               
                        !m3H2O = m3H20/m3cell * m3cell
                        OldVolume = Me%Theta(i,j,k) * Me%ExtVar%CellVolume(i,j,k)
                        
                        !m3H20/m3cell = (m3H20 + m3/s * s) / m3cell
                        NewTheta = (OldVolume + Me%iFlowBoundaryWalls(i,j,k) * Me%ExtVar%DT) / &
                                    Me%ExtVar%CellVolume(i,j,k)
                        
                        if (NewTheta .gt. (Me%RC%ThetaS(i,j,k) - Me%CV%LimitThetaHi)) then
                            NewTheta = Me%RC%ThetaS(i,j,k)
                        elseif (Me%Theta(i,j,k) .lt. (Me%RC%ThetaR(i,j,k) + Me%CV%LimitThetaLo)) then
                            NewTheta = Me%RC%ThetaR(i,j,k)
                        endif
                        
                        Me%Theta(i,j,k) = NewTheta
                        
                        !m3/s = m3 / s - if new > old positive flow
                        Me%iFlowBoundaryWalls(i,j,k)  = ((Me%Theta(i,j,k) * Me%ExtVar%CellVolume(i,j,k))  &
                                                    - OldVolume) / Me%ExtVar%DT
                        sum = sum + (Me%iFlowBoundaryWalls(i,j,k) * Me%ExtVar%DT)
                        
                        Me%RC%ThetaF  (i, j, k)   = ThetaF_ (Me%Theta (i, j, k), Me%SoilID(i, j, k))
                        Me%Head       (i, j, k)   = Head_   (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))
                        Me%UnSatK     (i, j, k)   = UnsatK_ (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))                        
                        
                    enddo do1

                endif
               
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        
            Me%AccBoundaryFlowVolume = sum
        
        endif
        
        !Impose bottom flux - exterior theta = interior theta -> flux = conductivity * Area
        if (Me%SoilOpt%ImposeBoundaryBottom .and. Me%Boundary%ImposeBoundaryBottomCondition == NullGradient_) then

            call SetMatrixValueAllocatable (Me%iFlowBoundaryBottom, Me%Size2D, 0.0, Me%ExtVar%BasinPoints)

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

            sum = Me%AccBoundaryFlowVolume

            !$OMP PARALLEL PRIVATE(I,J,K,di,dj,BoundaryFinalHead,ConductivityFace,AreaZX,AreaZY,OldVolume,NewTheta) REDUCTION(+:sum)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExtVar%BasinPoints (i, j) == 1) then
                
                    k = Me%ExtVar%KFloor(i,j)
                    
                    !Process cells that are bottom and that have ThetaF higher than minimum
                    if (Me%RC%ThetaF(i, j, k) > Me%Boundary%MinThetaFForBoundary) then
                     
                        !Conductivity is the maximum (saturation front entering or exiting domain)
                        ConductivityFace  = Me%UnSatK(i,j,k)
                        
                        !flow negative (exiting soil)
                        !Assuming that exterior and interior theta are the same (null_gradient)
                        !than Final Head gradient is dz/dz = 1 and BuckingamDarcy velocity is 
                        !interior conductivity and water exists soil (gravity). 
                        !Assuming also no hydrostatic pressure (since water is moving through bottom)
                        !m3/s = m2 * m/s
                        Me%iFlowBoundaryBottom(i,j) =  - Me%ExtVar%Area(i,j) * Me%UnSatK(i,j,k)

                        !m3H2O = m3H20/m3cell * m3cell
                        OldVolume = Me%Theta(i,j,k) * Me%ExtVar%CellVolume(i,j,k)
                        
                        !m3H20/m3cell = (m3H20 + m3/s * s) / m3cell
                        NewTheta = (OldVolume + Me%iFlowBoundaryBottom(i,j) * Me%ExtVar%DT) / &
                                    Me%ExtVar%CellVolume(i,j,k)
                        
                        if (Me%Theta(i,j,k) .lt. (Me%RC%ThetaR(i,j,k) + Me%CV%LimitThetaLo)) then
                            NewTheta = Me%RC%ThetaR(i,j,k)
                            
                            !m3/s = m3 / s - if new > old positive flow
                            Me%iFlowBoundaryBottom(i,j)  = ((NewTheta * Me%ExtVar%CellVolume(i,j,k))  &
                                                           - OldVolume) / Me%ExtVar%DT                            
                        endif
                        
                        Me%Theta(i,j,k) = NewTheta
                        
                        sum = sum + (Me%iFlowBoundaryBottom(i,j) * Me%ExtVar%DT)
                        
                        Me%RC%ThetaF  (i, j, k)   = ThetaF_ (Me%Theta (i, j, k), Me%SoilID(i, j, k))
                        Me%Head       (i, j, k)   = Head_   (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))
                        Me%UnSatK     (i, j, k)   = UnsatK_ (Me%RC%ThetaF (i, j, k), Me%SoilID(i, j, k))                            

                    endif
                
                endif
               
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        
            Me%AccBoundaryFlowVolume = sum
        
        endif
    
    end subroutine ModifyBoundaryFlux
    
    !--------------------------------------------------------------------------

    subroutine CalculateUGWaterLevel
    
        !Arguments-------------------------------------------------------------
            
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, CHUNK
        real                                        :: CellBottomLevel, DZInCell
        real                                        :: FieldHead, FieldTheta

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,K, CellBottomLevel, DZInCell, FieldHead, FieldTheta)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == 1) then

                !Find 1. non saturated cell            
doK:            do K = Me%ExtVar%KFloor(i, j), Me%WorkSize%KUB
                    !Do not use LimitThetaHigh for stability reasons 
                    !(usually too small producing GWLevels going up and down from iteration to iteration)
                    !if (Me%Theta(i,j,k) < Me%RC%ThetaS(i, j, k) - Me%CV%LimitThetaHi) then
                    !if (Me%Theta(i,j,k) < Me%RC%ThetaS(i, j, k) - Me%CV%LimitThetaHiGWTable) then
                    if (Me%Theta(i,j,k) < (Me%CV%LimitThetaHiGWTable * Me%RC%ThetaS(i, j, k))) then
                        exit doK
                    endif
                enddo doK
                
                k = min(k, Me%WorkSize%KUB)
                
                CellBottomLevel = - Me%ExtVar%SZZ(i,j,k-1)
                FieldHead       = - 0.5 * Me%ExtVar%DWZ(i, j, k)
                FieldTheta      = Theta_(FieldHead, Me%SoilID(i,j,k))
                
                if (Me%Theta(i,j,k) < FieldTheta) then
                    DZInCell        = 0.0
                else
                    DZInCell        = LinearInterpolation(FieldTheta, 0.0, Me%RC%ThetaS(i,j,k), &
                                                          Me%ExtVar%DWZ(i,j,k), Me%Theta(i,j,k))
                endif
                
                if (Me%IntegrationOutput%WaterTable%Yes) &
                    Me%OldUGWaterLevel2D(i, j) = Me%UGWaterLevel2D(i, j)
                
                Me%UGWaterLevel2D(i, j) = CellBottomLevel + DZInCell
                Me%UGWaterDepth2D(i, j) = Me%ExtVar%Topography(i, j) -Me%UGWaterLevel2D(i, j)
                !PMP need the same UGCell that PM used for computing DN flux
                Me%UGCell_Old    (i,j)  = Me%UGCell (i, j)
                Me%UGCell        (i, j) = k
                
            endif
            
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
                
    end subroutine CalculateUGWaterLevel

    !--------------------------------------------------------------------------

    real function Head_(thf, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: thf
        integer, intent(IN)                         :: SoilID

        !Local-------------------------------------------------------------------

        Head_= - abs( ( ( thf**( -1.0 / Me%SoilTypes(SoilID)%MFit)) -1.0 )**( 1.0 / Me%SoilTypes(SoilID)%NFit ) ) / &
                Me%SoilTypes(SoilID)%Alfa

    end function Head_   

    !--------------------------------------------------------------------------

    real function Head_slope_(th, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: th
        integer, intent(IN)                         :: SoilID
        
        !Local-------------------------------------------------------------------
        real                                        :: thf

        !------------------------------------------------------------------------

        if (th < Me%SoilTypes(SoilID)%ThetaS - Me%CV%LimitThetaHi) then

            if (th > Me%SoilTypes(SoilID)%ThetaR) then
                thf = ThetaF_(th, SoilID)
            else 
                thf = Me%CV%LimitThetaLo / (Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR)
            endif

            Head_slope_= (1./Me%SoilTypes(SoilID)%Alfa/Me%SoilTypes(SoilID)%NFit) * &
                        (thf**(-1./Me%SoilTypes(SoilID)%MFit)-1.) ** ((1.-Me%SoilTypes(SoilID)%NFit) / &
                        Me%SoilTypes(SoilID)%NFit) * (1./Me%SoilTypes(SoilID)%MFit) * &
                        thf ** (-(1.+ Me%SoilTypes(SoilID)%MFit) / Me%SoilTypes(SoilID)%MFit) * 1. /   &
                        ( Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR)
            return
        else
            
            Head_slope_ = Me%SoilTypes(SoilID)%OverSatSlope
            return

        end if
            

    end function Head_slope_
    
    !--------------------------------------------------------------------------

    real function Theta_ (head, SoilID)

        !Arguments-------------------------------------------------------------------
        real, intent(IN)                            :: head
        integer, intent(IN)                         :: SoilID

        !----------------------------------------------------------------------

        Theta_ = Me%SoilTypes(SoilID)%ThetaR + ( Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR ) / &
                ((1.0 + (Me%SoilTypes(SoilID)%alfa * (- head) ) ** Me%SoilTypes(SoilID)%Nfit ) **              &
                Me%SoilTypes(SoilID)%MFit )

    end function Theta_

    !----------------------------------------------------------------------------
    
    real function ThetaF_ (th, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: th
        integer, intent(IN)                         :: SoilID
        
        ThetaF_ = (th - Me%SoilTypes(SoilID)%ThetaR) / (Me%SoilTypes(SoilID)%ThetaS - Me%SoilTypes(SoilID)%ThetaR)     

    end function ThetaF_ 

    !----------------------------------------------------------------------------

    real function UnsatK_(thf, SoilID)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: thf
        integer, intent(IN)                         :: SoilID

        !Local-----------------------------------------------------------------

        UnsatK_ = Me%SoilTypes(SoilID)%SatK*((thf)**Me%SoilTypes(SoilID)%lfit) * &
                (1.0d0-(1.0d0-thf**(1.0d0/Me%SoilTypes(SoilID)%MFit))**Me%SoilTypes(SoilID)%MFit)**2.0d0
                    
    end function UnsatK_

    !--------------------------------------------------------------------------

!    subroutine ExchangeWithDrainageNetwork2
!
!        !Arguments-------------------------------------------------------------
!        
!        !Local-----------------------------------------------------------------
!        integer                                     :: i, j, STAT_CALL
!        real                                        :: Area
!        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
!        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
!        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
!        real,   dimension(:, :), pointer            :: ChannelsNodeLength
!        integer,dimension(:, :), pointer            :: ChannelsOpenProcess
!
!        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'
!
!        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'
!
!        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'
!
!        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'
!
!        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'
!
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!
!            if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then
!            
!                !Channel will loose water
!                if (ChannelsWaterLevel(i, j) > Me%UGWaterLevel2D(i, j)) then
!
!                    
!                    if (ChannelsOpenProcess(i, j) == 1) then
!
!                        !The Area is always two times the lateral area and
!                        !one time the bottom area
!                        Area  = (2. * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))  & 
!                                     + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
!                                     
!                        !Positive Flow to channel when Channel "gains" water
!                        ![m3/s]                    [m2]   [m/s]                     
!                        Me%lFlowToChannels(i, j) = -1.0 * Area * Me%UnSatK(i, j, Me%UGCell(i,j))
!
!                        if (-1. * Me%lFlowToChannels(i, j) * Me%ExtVar%DT > (ChannelsWaterLevel(i, j)      &
!                            - ChannelsBottomLevel(i, j)) * Area) then
!                            Me%lFlowToChannels(i, j) = -0.5 * (ChannelsWaterLevel(i, j)                    &
!                                                       - ChannelsBottomLevel(i, j)) * Area / Me%ExtVar%DT
!                            call SetError(WARNING_, INTERNAL_, "Flow to channel corrected. FLOW TO CHANNEL", OFF)
!                        endif
!                        
!                    else
!                    
!                        Me%lFlowToChannels(i, j) = 0.0
!                    
!                    endif
!                        
!                    
!                else
!            
!                    Area = (2. * (Me%UGWaterLevel2D(i, j) - ChannelsBottomLevel(i, j)) + &
!                            ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
!
!                    !Positive Flow to channel when Channel "gains" water
!                    ![m3/s]                    [m2]   [m/s]                     
!                    Me%lFlowToChannels(i, j) = Area * Me%UnSatK(i, j, Me%UGCell(i,j))
!                    
!                endif
!
!            else
!
!                Me%lFlowToChannels(i, j) = 0.0
!
!            endif
!
!        enddo
!        enddo
!        
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR10'
!
!    end subroutine ExchangeWithDrainageNetwork2
    
    !--------------------------------------------------------------------------
    
    !This subroutine assumes that near the river the aquifer level and river level converge to river level
    !so the flux area is always the river section. 
    !Because of this when there is no water in river or river level is lower than soil bottom there is no flux
    !so this was replaced with ExchangeWithDrainageNet_2
    subroutine ExchangeWithDrainageNet_1

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL
        real                                        :: dH, dX, TotalArea !, LayerArea
        real                                        :: MinHead, FieldTheta !, TopHeightInCell
        real                                        :: InfiltrationVolume, ChannelsVolume
        real                                        :: MaxFlow !, VerticalFluxVariation
        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
        real,   dimension(:, :), pointer            :: ChannelsNodeLength
        integer,dimension(:, :), pointer            :: ChannelsOpenProcess
        !logical                                     :: FoundUpperFlowCell, FoundLowerFlowCell
        !Begin----------------------------------------------------------------

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'

        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'

        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'
        
       
        !!Computation of TotalFlow - always
do1:    do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then
                if ((ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) > 0.0) then            
                
    !                if (ChannelsBottomLevel(i,j) < Me%ExtVar%BottomTopoG(i, j)) then
    !                    write (*,*) 
    !                    write (*,*) 'Bottom River section is lower than soil profile in cell', i,j
    !                    write (*,*) 'Increase bottom soil depth or decrease river depth'
    !                    stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'
    !                endif
                
                    !Not Tested yet if bottom channel lower than soil
                    if (ChannelsBottomLevel(i,j) < Me%ExtVar%BottomTopoG(i, j)) then
                    
                        !if the channel has water above soil bottom
                        if (ChannelsWaterLevel(i, j) > Me%ExtVar%BottomTopoG(i, j)) then
                        
                            !Area is only lateral from bottom
                            TotalArea = 2. * (ChannelsWaterLevel(i, j) - Me%ExtVar%BottomTopoG(i, j)) * ChannelsNodeLength(i, j)
                        else
                            !There is no channel water above soil bottom and there will be no flow.
                            TotalArea = 0.0
                        endif    
                    else                
                        !Computing Area for flux in normal case - is always two times the lateral area and one time the bottom area
                        TotalArea  = (2. * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))  &
                                      + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
                    endif
                
                    !Computing height gradient dH [m]
                    !Negative dh -> Flux from channels to porous media
                    !Positive dh -> Flux from porous media to channels 
                    dH    = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j)                
                
                    !spatial step for Height Gradient (dH) [m]
                    dX    = max(ChannelsBottomWidth(i, j) / 2.0, 1.0)
                    
                    !if flow from river put the fux on first unsaturated cell on top of saturated aquifer (UGCell)
                    !if flow from soil remove from upper cell of saturated aquifer (UGCell - 1). 
                    !This avoids removing water from low water content layer above aquifer 
                    !(can cause instability if using minimum conductivity)
                    !Do not test this now. Need to update Routine CalculateNewTheta and adapt 
                    !PorousMediaProperties use of flux also for the right flow at the right cell                    
!                    if (dH < 0) then
                        k = Me%UGCell(i,j)
!                    else
!                        k = min (Me%UGCell(i,j) - 1, Me%ExtVar%KFloor(i,j))
!                    endif
                    
                    ![m3/s]                   = [m/m] * [m2] * [m/s] * [] 
                    !Me%lFlowToChannels(i, j) = (dH / dX ) * TotalArea * Me%UnSatK(i, j, Me%UGCell(i,j)) * Me%SoilOpt%HCondFactor
                    !Me%lFlowToChannels(i, j) = (dH / dX ) * TotalArea * Me%SatK(i, j, Me%UGCell(i,j))             &
                    !                             * Me%SoilOpt%FCHCondFactor !Me%SoilOpt%HCondFactor
                    Me%lFlowToChannels(i, j) = (dH / dX ) * TotalArea * Me%SatK(i, j, k)                  &
                                                 * Me%SoilOpt%FCHCondFactor 
                
                    !If the channel looses water (infiltration), then set max flux so that volume in channel does not get 
                    !negative.                
                    if (dH < 0) then


                        ![m3]
                        ChannelsVolume     = (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))               &
                                              * ChannelsBottomWidth(i, j) * ChannelsNodeLength(i, j)
                        InfiltrationVolume = -1. * Me%lFlowToChannels(i, j) * Me%CV%CurrentDT
                        !InfiltrationVolume = -1. * Me%lFlowToChannels(i, j) * Me%ExtVar%DT                

                        if (InfiltrationVolume > 0.5 * ChannelsVolume) then
                            Me%lFlowToChannels(i, j) = -0.50 * ChannelsVolume / Me%CV%CurrentDT
                            !Me%lFlowToChannels(i, j) = -0.50 * ChannelsVolume / Me%ExtVar%DT
                        endif

                            
                        !This will only infiltrate water when there is more then enough... Increase numerical stability
                        if (.not. Me%SoilOpt%DryChannelsCompletely) then
                    
                            if (ChannelsOpenProcess(i, j) == 0) then 
                                Me%lFlowToChannels(i, j) = 0.0
                            endif
                    
                        endif
                    

                
                    !If soil looses water set flow so that cell stays at least with field theta
                    !in the top cell of aquifer head at field capacity in cell center 
                    !will be half the height
                    elseif (dH > 0) then
                        !k           = Me%UGCell(i,j)
                        MinHead     = -0.5 * Me%ExtVar%DWZ(i, j, k)
                        FieldTheta  = Theta_(MinHead, Me%SoilID(i,j,k))
                        
                        !Maxflow was erroneous. It has to use actual Theta and not ThetaS because content may not be
                        !completely saturated       
                        !m3/s   = (-) * m3 / s
                        if (Me%Theta (i,j,k) <= FieldTheta) then
                            MaxFlow = 0.
                        else
                            MaxFlow   = (Me%Theta (i,j,k) - FieldTheta) * Me%ExtVar%CellVolume(i,j,k) / Me%CV%CurrentDT
                            !MaxFlow   = (Me%Theta (i,j,k) - FieldTheta) * Me%ExtVar%CellVolume(i,j,k) / Me%ExtVar%DT
                            !MaxFlow   = (Me%RC%ThetaS (i,j,k) - FieldTheta) * Me%ExtVar%CellVolume(i,j,k) / Me%ExtVar%DT
                        endif
                        
                        if (Me%lFlowToChannels(i,j) > MaxFlow) then
                            Me%lFlowToChannels(i,j) = MaxFlow
                        endif
                    endif
                
                else

                    Me%lFlowToChannels(i, j) = 0.0
                endif
            else
                Me%lFlowToChannels(i, j) = 0.0
            endif
            
        enddo do2
        enddo do1

        
!        !!Computation of Flow by layers - from total flow computed above.
!        if (Me%SoilOpt%DNLink == GWFlowToChanByLayer_) then
!
!            call SetMatrixValueAllocatable (Me%lFlowToChannelsLayer, Me%Size, 0.0, Me%ExtVar%WaterPoints3D)
!           
!do3:        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!do4:        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!                
!                !Flow is zero if water column in river is lower than soil bottom
!                if ((Me%ExtVar%RiverPoints(i, j) == OpenPoint) .and.                                  &
!                    (ChannelsWaterLevel(i,j) .gt. Me%ExtVar%BottomTopoG(i, j))) then
!                    
!                    if (ChannelsBottomLevel(i,j) < Me%ExtVar%BottomTopoG(i, j)) then
!                        FoundLowerFlowCell = .true.
!                        Me%FlowToChannelsBottomLayer(i,j) = Me%ExtVar%KFloor(i,j)
!                    else
!                        FoundLowerFlowCell    = .false.
!                    endif
!                    FoundUpperFlowCell    = .false.
!                    VerticalFluxVariation = 0.0
!                    !Search for lower and upper cell computation - water column in river
!                    !And compute area for flow inside layer
!do5:                do K = Me%ExtVar%KFloor(i,j), Me%WorkSize%KUB
!                        
!                        if (.not. FoundLowerFlowCell) then
!                        
!                            !SZZ(i,j,k) is -altitude at the top face. But has also the bottom face for the kfloor
!                            if((-Me%ExtVar%SZZ(i,j,k-1) .le. ChannelsBottomLevel(i, j)) .and.          &
!                               (-Me%ExtVar%SZZ(i,j,k  ) .gt. ChannelsBottomLevel(i, j))) then
!                            
!                                FoundLowerFlowCell                = .true.
!                                Me%FlowToChannelsBottomLayer(i,j) = k
!                                
!                                !Area for bottom river inside cell (2*height in cell + bottom width) * length
!                                TopHeightInCell = min(-Me%ExtVar%SZZ(i,j,k), ChannelsWaterLevel(i, j))
!                                LayerArea       = (2 * (TopHeightInCell - ChannelsBottomLevel(i, j))     &
!                                                  + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
!                                
!                                !Check if Top channel height is also inside the channel bottom cell
!                                if (ChannelsWaterLevel(i, j) .lt. -Me%ExtVar%SZZ(i,j,k)) then
!                                    FoundUpperFlowCell             = .true.
!                                    Me%FlowToChannelsTopLayer(i,j) = k
!                                endif
!                                
!                            else
!                                !Cell Lower than channel bottom
!                                LayerArea = 0.0
!                            endif
!                        else
!                            if (.not. FoundUpperFlowCell) then
!                                if((-Me%ExtVar%SZZ(i,j,k-1) .lt. ChannelsWaterLevel(i, j)) .and.               &
!                                   (-Me%ExtVar%SZZ(i,j,k  ) .ge. ChannelsWaterLevel(i, j))) then
!                                
!                                    FoundUpperFlowCell             = .true.
!                                    Me%FlowToChannelsTopLayer(i,j) = k
!                                    
!                                    !Area for river inside cell (2*height in cell) * length
!                                    LayerArea = (2 * (ChannelsWaterLevel(i, j) - (-Me%ExtVar%SZZ(i,j,k-1))))   &
!                                                 * ChannelsNodeLength(i, j)
!                                    
!                                else
!                                    !cells between river bottom and water level
!                                    !Area for river inside cell (2*cellheight) * length
!                                    LayerArea = 2 * (Me%ExtVar%DWZ(i,j,k)) * ChannelsNodeLength(i, j)
!                                
!                                endif
!                            else
!                                if (K .ge. Me%UGCell(i,j)) then
!                                    !if found top river cell and higher than GW cell do not need to continue 
!                                    !(all computations may cease).
!                                    exit do5
!                                else
!                                    !cells higher than river level - no flow 
!                                    LayerArea = 0.0
!                                endif
!                            endif
!                        endif
!                        
!                        !Water level greater than soil top - top k undefined
!                        if (ChannelsWaterLevel(i, j) .gt. -Me%ExtVar%SZZ(i,j,Me%WorkSize%KUB)) then
!                            Me%FlowToChannelsTopLayer(i,j) = Me%WorkSize%KUB
!                        endif
!                        
!                        !Total area computed above - area related to total flow: Me%lFlowToChannels(i,j)
!                        TotalArea  = (2. * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))  &
!                                      + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)  
!                                                      
!                        !The layer flow is proportional to the [layer area] / [total flow area]                       
!                        ![m3/s] = [m3/s] * [m2l]/[m2total]
!                        Me%lFlowToChannelsLayer(i, j, k) = Me%lFlowToChannels(i,j) * (LayerArea / TotalArea)
!                        
!                        !if flow is computed in saturated area then vertical flows have to be compensated so
!                        !that water content is maintained
!                        !This will result in a transport in vertical up to the groundwater top cell
!                        if ((K .ge. Me%FlowToChannelsBottomLayer(i,j)) .and. (K .lt. Me%UGCell(i,j))) then
!                            
!                            !Vertical flux (upper face) has to compensate the lower face flux and the flux with river. 
!                            !if flow to channels positive (to river) vertical flow is negative (downwards to compensate exit)
!                            !if flow to channels negative (to soil) vertical flow is positive (upwards to compensate entrance)
!                            ![m3/s] in iteration (dt is the same for both fluxes)
!                            VerticalFluxVariation = VerticalFluxVariation - Me%lFlowToChannelsLayer(i, j, k)
!                            Me%FluxWFinal(i,j,k+1)  = Me%FluxWFinal(i,j,k+1) + VerticalFluxVariation
!                            
!                            !m/s = m/s + (m3/s / m2)
!                            Me%UnsatVelWFinal(i,j,k+1)  = Me%UnsatVelWFinal(i,j,k+1)                         &
!                                                          + (Me%FluxWFinal(i,j,k+1) / Me%ExtVar%Area(i,j))
!                        endif
!                    
!                    enddo do5
!                
!                endif
!
!            enddo do4
!            enddo do3
!        
!        endif
        
        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR10'

    end subroutine ExchangeWithDrainageNet_1

    !--------------------------------------------------------------------------

    !This is the same as ExchangeWithDrainageNet_1 (flux in soil occurs at one cell top of aquifer) 
    !but with areas for flux changed to account the river wet area + aquifer if
    !aquifer is higher than river (lateral water entering trough the margins)
    subroutine ExchangeWithDrainageNet_2

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL, CHUNK
        real                                        :: dH, dX, TotalArea
        real                                        :: TotalHeight, Toplevel, BottomLevel, ChannelColumn
        real                                        :: MinHead, FieldTheta
        real                                        :: InfiltrationVolume, ChannelsVolume
        real                                        :: MaxFlow
        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
        real,   dimension(:, :), pointer            :: ChannelsNodeLength
        integer,dimension(:, :), pointer            :: ChannelsOpenProcess
        !Begin----------------------------------------------------------------

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'

        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'

        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'
       
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
       
        !$OMP PARALLEL PRIVATE(I,J,K,dH,ChannelColumn,Toplevel,BottomLevel,TotalHeight,TotalArea,dX, &
        !$OMP& ChannelsVolume,InfiltrationVolume,MinHead,FieldTheta,MaxFlow)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)       
do1:    do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do I = Me%WorkSize%ILB, Me%WorkSize%IUB

if1:        if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then

                !Computing height gradient dH [m]
                !Negative dh -> Flux from channels to porous media
                !Positive dh -> Flux from porous media to channels 
                dH            = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j) 
                ChannelColumn = ChannelsWaterLevel(i,j) - ChannelsBottomLevel(i, j)
                
                !Flux only occurs if there is gradient or water
if2:            if ((dH > 0.0) .or. (dH < 0 .and. ChannelColumn > 0.)) then

                    !Maximum between the aquifer and river, limited by soil top
                    Toplevel    = min(max(Me%UGWaterLevel2D(i, j), ChannelsWaterLevel(i, j)), Me%ExtVar%Topography(i, j))
                    !Maximum between the river bottom and porous media bottom
                    BottomLevel = max(ChannelsBottomLevel(i, j), Me%ExtVar%BottomTopoG(i, j)) 
                    !Hegight for flux
                    TotalHeight = TopLevel - BottomLevel
                    
                    !Compute total area for flux
                    !if soil bottom above channel bottom do not consider bottom area
                    if (Me%ExtVar%BottomTopoG(i, j) > ChannelsBottomLevel(i, j)) then 

                        !Computing Area for flux - is two times the lateral area
                        TotalArea  = (2. * TotalHeight) * ChannelsNodeLength(i, j)
                                                                
                    else 
                        
                        !Computing Area for flux - is two times the lateral area and one time the bottom area
                        TotalArea  = (2. * TotalHeight + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
                   
                    endif
                
                    !spatial step for Height Gradient (dH) [m]
                    !review this!
                    dX    = max(ChannelsBottomWidth(i, j) / 2.0, 1.0)
                
                    !Layer to compute flux - first non saturated layer
                    k = Me%UGCell(i,j)
                    
                    ![m3/s]                   = [m/m] * [m2] * [m/s] * [] 
                    Me%lFlowToChannels(i, j) = (dH / dX ) * TotalArea * Me%SatK(i, j, k) * Me%SoilOpt%FCHCondFactor 
                
                    !If the channel looses water (infiltration), then set max flux so that volume in channel does not get 
                    !negative.   
                    if (dH < 0) then             

                        ![m3]
                        ChannelsVolume     = (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))               &
                                              * ChannelsBottomWidth(i, j) * ChannelsNodeLength(i, j)
                        InfiltrationVolume = -1. * Me%lFlowToChannels(i, j) * Me%CV%CurrentDT

                        if (InfiltrationVolume > 0.5 * ChannelsVolume) then
                            Me%lFlowToChannels(i, j) = -0.50 * ChannelsVolume / Me%CV%CurrentDT
                        endif
                            
                        !This will only infiltrate water when there is more then enough... Increase numerical stability
                        if (.not. Me%SoilOpt%DryChannelsCompletely) then
                    
                            if (ChannelsOpenProcess(i, j) == 0) then 
                                Me%lFlowToChannels(i, j) = 0.0
                            endif
                    
                        endif
                
                    !If soil looses water set flow so that cell stays at least with field theta
                    !in the top cell of aquifer head at field capacity in cell center 
                    !will be half the height
                    elseif (dH > 0) then
                        !k           = Me%UGCell(i,j)
                        MinHead     = -0.5 * Me%ExtVar%DWZ(i, j, k)
                        FieldTheta  = Theta_(MinHead, Me%SoilID(i,j,k))
                        
                        !m3/s   = (-) * m3 / s
                        if (Me%Theta (i,j,k) <= FieldTheta) then
                            MaxFlow = 0.
                        else
                            MaxFlow   = (Me%Theta (i,j,k) - FieldTheta) * Me%ExtVar%CellVolume(i,j,k) / Me%CV%CurrentDT
                        endif
                        
                        if (Me%lFlowToChannels(i,j) > MaxFlow) then
                            Me%lFlowToChannels(i,j) = MaxFlow
                        endif
                    endif
                
                else
                    Me%lFlowToChannels(i, j) = 0.0
                endif if2
            
            endif if1
            
        enddo do2
        enddo do1
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR10'

    end subroutine ExchangeWithDrainageNet_2

    !--------------------------------------------------------------------------
    
    !In this subroutine the area for flux accounts the river wet area + aquifer if
    !aquifer is higher than river (lateral water entering trough the margins)  
    !The flux is computed per layer and not from/to top of aquifer (as in 
    ! ExchangeWithDrainageNet_2 and ExchangeWithDrainageNet_1) 
    !The latter routines have the disvantage
    !that top of aquifer may be much lower than river (flow from river) and flux needs to pass
    !by all other layers and not go directly to aquifer
    !or flux may be limited if low water content exists just above aquifer (flow from soil).
    !In saturated areas the vertical flux is updated to mantain saturation.    
    subroutine ExchangeWithDrainageNet_3

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL, CHUNK
        real                                        :: dH, dX, TotalArea
        real                                        :: TotalHeight, Toplevel, BottomLevel, ChannelColumn
        real                                        :: LayerArea, LayerHeight, sumLayerArea
        real                                        :: MaxFlow, sumFlow, VerticalFluxVariation
        logical                                     :: AccountBottomSurface
        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
        real,   dimension(:, :), pointer            :: ChannelsNodeLength
        integer,dimension(:, :), pointer            :: ChannelsOpenProcess
        !Begin----------------------------------------------------------------

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'

        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'

        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'
       
        !compute where are layer limts to compute fluxes. optimize the code
        call GetLayersLimits_3
       
        !Nuliffy layer fluxes
        call SetMatrixValueAllocatable (Me%lFlowToChannelsLayer, Me%Size, 0.0, Me%ExtVar%WaterPoints3D)
        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K,dH,ChannelColumn,Toplevel,BottomLevel,TotalHeight,TotalArea,dX, &
        !$OMP& sumLayerArea,sumFlow,VerticalFluxVariation,LayerHeight,AccountBottomSurface,LayerArea,MaxFlow)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)       
do1:    do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then

                !Computing height gradient dH [m]
                !Negative dh -> Flux from channels to porous media
                !Positive dh -> Flux from porous media to channels 
                dH            = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j) 
                ChannelColumn = ChannelsWaterLevel(i,j) - ChannelsBottomLevel(i, j)
                
                !Flux only occurs if there is gradient or water                
                if ((dH > 0.0) .or. (dH < 0 .and. ChannelColumn > 0.)) then
                
                    !Maximum between the aquifer and river, limited by soil top
                    Toplevel    = min(max(Me%UGWaterLevel2D(i, j), ChannelsWaterLevel(i, j)), Me%ExtVar%Topography(i, j))
                    !Maximum between the river bottom and porous media bottom
                    BottomLevel = max(ChannelsBottomLevel(i, j), Me%ExtVar%BottomTopoG(i, j)) 
                    !Hegight for flux
                    TotalHeight = TopLevel - BottomLevel

                    !Compute total area for flux
                    !if soil bottom above channel bottom do not consider bottom area
                    if (Me%ExtVar%BottomTopoG(i, j) > ChannelsBottomLevel(i, j)) then 

                        !Computing Area for flux - is two times the lateral area
                        TotalArea  = (2. * TotalHeight) * ChannelsNodeLength(i, j)
                                                                
                    else 
                        
                        !Computing Area for flux - is two times the lateral area and one time the bottom area
                        TotalArea  = (2. * TotalHeight + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)

                    endif

                    !spatial step for Height Gradient (dH) [m]
                    !review this!
                    dX    = max(ChannelsBottomWidth(i, j) / 2.0, 1.0)
                
                    sumLayerArea          = 0.0
                    sumFlow               = 0.0
                    VerticalFluxVariation = 0.0
                   
                   !go trough all active flux layers (computed in the routine GetLayersLimits
do3:                do K = Me%FlowToChannelsBottomLayer(i,j), Me%FlowToChannelsTopLayer(i,j)
                        
                        !if layers equal the height is level difference
                        if (Me%FlowToChannelsTopLayer(i,j) == Me%FlowToChannelsBottomLayer(i,j)) then
                            
                            LayerHeight = TotalHeight
                            
                            if (Me%ExtVar%BottomTopoG(i, j) > ChannelsBottomLevel(i, j)) then
                                AccountBottomSurface = .false.
                            else
                                AccountBottomSurface = .true.
                            endif                            
                        
                        !in top layer (the layer where TopLevel was found) 
                        !the layer height for flux is depenedent on top level
                        elseif (K == Me%FlowToChannelsTopLayer(i,j)) then
                            
                            LayerHeight = TopLevel - (-Me%ExtVar%SZZ(i,j,k-1))
                            AccountBottomSurface = .false. 
                        
                        !in bottom layer (the layer where TopLevel was found)
                        ! the layer height for flux is depenedent on bottom level and 
                        !need to check if account bottom section width or not (flux trough the bottom)
                        elseif (K == Me%FlowToChannelsBottomLayer(i,j)) then
                            
                            LayerHeight = -Me%ExtVar%SZZ(i,j,k) - BottomLevel                                            
                            
                            if (Me%ExtVar%BottomTopoG(i, j) > ChannelsBottomLevel(i, j)) then
                                AccountBottomSurface = .false.
                            else
                                AccountBottomSurface = .true.
                            endif
                        
                        !otherwise the height is the cell height
                        else
                            
                            LayerHeight = Me%ExtVar%DWZ(i,j,k)
                            AccountBottomSurface = .false.
                        
                        endif

                        if (AccountBottomSurface) then
                            LayerArea   = (2. * LayerHeight + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
                        else
                            LayerArea   = 2. * LayerHeight * ChannelsNodeLength(i, j)
                        endif 
                        
                        sumLayerArea = sumLayerArea + LayerArea
                        
                        !flux in layer is total flux * fraction LayerArea/TotalArea
                        !m3/s = - * m2total * m2layer/m2total * m/s * -
                        Me%lFlowToChannelsLayer(i, j, k)  = (dH / dX ) * TotalArea * (LayerArea/TotalArea) *     &
                                                             Me%SatK(i, j, k) * Me%SoilOpt%FCHCondFactor 
                                                             
                        !limit flux
                        if (dH < 0) then
                            !limit to water existing in virtual river layer (the water in between layer thickness)
                            !m3/s = m * m * m / s
                            MaxFlow = - LayerHeight * ChannelsBottomWidth(i, j) * ChannelsNodeLength(i, j) / Me%CV%CurrentDT
                        else
                            !limit to residual
                            !m3/s = m3H20/m3cell * m3cell / s
                            MaxFlow = (Me%Theta(i,j,k) - Me%RC%ThetaR(i,j,k)) * Me%ExtVar%CellVolume(i,j,k) / Me%CV%CurrentDT
                        endif
                        
                        if (abs(Me%lFlowToChannelsLayer(i, j, k)) > abs(MaxFlow)) then
                            Me%lFlowToChannelsLayer(i, j, k) = MaxFlow
                        endif
                        
                        sumFlow = sumFlow + Me%lFlowToChannelsLayer(i, j, k)
                        
                        Me%lFlowToChannels(i, j) = sumFlow

                        !to maintain saturation vertical flux has to compensate
                        !and go up to top of aquifer
                        if (K .lt. Me%UGCell(i,j)) then
                            
                            !Vertical flux (upper face) has to compensate the lower face flux and the flux with river. 
                            !if flow to channels positive (to river) vertical flow is negative (downwards to compensate exit)
                            !if flow to channels negative (to soil) vertical flow is positive (upwards to compensate entrance)
                            ![m3/s] in iteration (dt is the same for both fluxes)
                            VerticalFluxVariation = VerticalFluxVariation - Me%lFlowToChannelsLayer(i, j, k)
                            Me%FluxWFinal(i,j,k+1)  = Me%FluxWFinal(i,j,k+1) + VerticalFluxVariation
                            
                            !m/s = m/s + (m3/s / m2)
                            Me%UnsatVelWFinal(i,j,k+1)  = Me%UnsatVelWFinal(i,j,k+1)                         &
                                                          + (Me%FluxWFinal(i,j,k+1) / Me%ExtVar%Area(i,j))
                        endif
                                        
                    enddo do3                    
                    
                    !Test if method was OK if sum of areas give the total area
                    if (abs(sumLayerArea - TotalArea) .ge. 1E-10) then
                        write(*,*)
                        write(*,*) 'DN_LINK : 2 in Porous Media got inconsistent results'
                        write(*,*)
                        stop 'ExchangeWithDrainageNet_3 - ModulePorousMedia - ERR01'
                    endif
                        
                endif
            
            endif
            
        enddo do2
        enddo do1
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR10'

    end subroutine ExchangeWithDrainageNet_3

    !--------------------------------------------------------------------------

    subroutine GetLayersLimits_2

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL, CHUNK
        real                                        :: dH, ChannelColumn
        real                                        :: TopLevel, Bottomlevel
        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
        real,   dimension(:, :), pointer            :: ChannelsNodeLength
        integer,dimension(:, :), pointer            :: ChannelsOpenProcess
        logical                                     :: FoundUpperFlowCell, FoundLowerFlowCell
        !Begin----------------------------------------------------------------

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR01'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR02'

        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR03'

        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR04'

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR05'

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,K,FoundUpperFlowCell,FoundLowerFlowCell,dH,ChannelColumn,Toplevel,BottomLevel)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3:    do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do4:    do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if ((Me%ExtVar%RiverPoints(i, j) == OpenPoint) ) then
                
                FoundUpperFlowCell = .false.  
                FoundLowerFlowCell = .false.  

                dH            = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j) 
                ChannelColumn = ChannelsWaterLevel(i,j) - ChannelsBottomLevel(i, j)
                
                !if flow exists
                !if ((dH > 0.0) .or. (dH < 0 .and. ChannelsOpenProcess)) then
                if ((dH > 0.0) .or. (dH < 0 .and. ChannelColumn > 0.)) then
                
                    !Maximum between the aquifer and river, limited by soil top
                    TopLevel    = min(max(ChannelsWaterLevel(i, j), Me%UGWaterLevel2D(i, j)), Me%ExtVar%Topography(i, j))
                    !Minimum between the aquifer and river, limited by channel bottom
                    BottomLevel = max(min(ChannelsWaterLevel(i, j), Me%UGWaterLevel2D(i, j)), ChannelsBottomLevel(i, j), &
                                  Me%ExtVar%BottomTopoG(i, j))
                    
                    !check first for limits
                    !out of soil boundaries
                    if (TopLevel == Me%ExtVar%Topography(i, j)) then
                        Me%FlowToChannelsTopLayer(i,j) = Me%WorkSize%KUB
                        FoundUpperFlowCell = .true.    
                    endif                        
                    
                    !go trough layes and find TopLevel and BottomLevel
do1:                do k = Me%WorkSize%KUB, Me%ExtVar%KFloor(i,j), -1
                        if (.not. FoundUpperFlowCell) then 
                            if((-Me%ExtVar%SZZ(i,j,k-1) .le. TopLevel) .and.          &
                               (-Me%ExtVar%SZZ(i,j,k  ) .gt. TopLevel)) then
                               Me%FlowToChannelsTopLayer(i,j) = k
                               FoundUpperFlowCell = .true.
                            endif
                        endif
                        if (.not. FoundLowerFlowCell) then
                            if((-Me%ExtVar%SZZ(i,j,k-1) .le. BottomLevel) .and.       &
                               (-Me%ExtVar%SZZ(i,j,k  ) .gt. BottomLevel)) then
                               Me%FlowToChannelsBottomLayer(i,j) = k
                               FoundLowerFlowCell = .true.
                            endif                            
                        endif
                        if (FoundUpperFlowCell .and. FoundLowerFlowCell) exit do1
                    enddo do1
                    
                    if (.not. FoundUpperFlowCell) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR01'
                    if (.not. FoundlowerFlowCell) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR010'
                
                endif
            endif
            
        enddo do4
        enddo do3
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL


        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_2 - ModulePorousMedia - ERR10'

    end subroutine GetLayersLimits_2

    !---------------------------------------------------------------------------

    subroutine GetLayersLimits_3

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL, CHUNK
        real                                        :: dH, ChannelColumn
        real                                        :: TopLevel, Bottomlevel
        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
        real,   dimension(:, :), pointer            :: ChannelsNodeLength
        integer,dimension(:, :), pointer            :: ChannelsOpenProcess
        logical                                     :: FoundUpperFlowCell, FoundLowerFlowCell
        !Begin----------------------------------------------------------------

        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR01'

        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR02'

        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR03'

        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR04'

        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR05'

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(I,J,K,FoundUpperFlowCell,FoundLowerFlowCell,dH,ChannelColumn,Toplevel,BottomLevel)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3:    do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do4:    do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if ((Me%ExtVar%RiverPoints(i, j) == OpenPoint) ) then
                
                FoundUpperFlowCell = .false.  
                FoundLowerFlowCell = .false.  

                dH            = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j) 
                ChannelColumn = ChannelsWaterLevel(i,j) - ChannelsBottomLevel(i, j)
                
                !if flow exists
                !if ((dH > 0.0) .or. (dH < 0 .and. ChannelsOpenProcess)) then
                if ((dH > 0.0) .or. (dH < 0 .and. ChannelColumn > 0.)) then
                
                    !Maximum between the aquifer and river, limited by soil top
                    TopLevel    = min(max(ChannelsWaterLevel(i, j), Me%UGWaterLevel2D(i, j)), Me%ExtVar%Topography(i, j))
                    !Maximum between the channel bottom and soil bottom
                    BottomLevel = max(ChannelsBottomLevel(i, j), Me%ExtVar%BottomTopoG(i, j)) 
                    
                    !check first for limits
                    !out of soil boundaries
                    if (TopLevel == Me%ExtVar%Topography(i, j)) then
                        Me%FlowToChannelsTopLayer(i,j) = Me%WorkSize%KUB
                        FoundUpperFlowCell = .true.    
                    endif                        
                    
                    !go trough layes and find TopLevel and BottomLevel
do1:                do k = Me%WorkSize%KUB, Me%ExtVar%KFloor(i,j), -1
                        if (.not. FoundUpperFlowCell) then 
                            if((-Me%ExtVar%SZZ(i,j,k-1) .le. TopLevel) .and.          &
                               (-Me%ExtVar%SZZ(i,j,k  ) .gt. TopLevel)) then
                               Me%FlowToChannelsTopLayer(i,j) = k
                               FoundUpperFlowCell = .true.
                            endif
                        endif
                        if (.not. FoundLowerFlowCell) then
                            if((-Me%ExtVar%SZZ(i,j,k-1) .le. BottomLevel) .and.       &
                               (-Me%ExtVar%SZZ(i,j,k  ) .gt. BottomLevel)) then
                               Me%FlowToChannelsBottomLayer(i,j) = k
                               FoundLowerFlowCell = .true.
                            endif                            
                        endif
                        if (FoundUpperFlowCell .and. FoundLowerFlowCell) exit do1
                    enddo do1
                    
                    if (.not. FoundUpperFlowCell) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR01'
                    if (.not. FoundlowerFlowCell) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR010'
                
                endif
            endif
            
        enddo do4
        enddo do3
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR06'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR07'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR08'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR09'

        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GetLayersLimits_3 - ModulePorousMedia - ERR10'

    end subroutine GetLayersLimits_3
    
    !---------------------------------------------------------------------------
    
!    subroutine ExchangeWithDrainageNet_Corr
!
!        !Arguments-------------------------------------------------------------
!        
!        !Local-----------------------------------------------------------------
!        integer                                     :: i, j, k, STAT_CALL
!        real                                        :: dH, dX, Area
!        real                                        :: MinHead, FieldTheta, MaxFlow
!        real,   dimension(:, :), pointer            :: ChannelsWaterLevel
!        real,   dimension(:, :), pointer            :: ChannelsBottomLevel 
!        real,   dimension(:, :), pointer            :: ChannelsBottomWidth
!        real,   dimension(:, :), pointer            :: ChannelsNodeLength
!        integer,dimension(:, :), pointer            :: ChannelsOpenProcess
!        real                                        :: ChannelsVolume, InfiltrationVolume
!        !Begin-----------------------------------------------------------------
!
!        call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'
!
!        call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'
!
!        call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'
!
!        call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'
!
!        call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'
!
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!
!            if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then
!                
!             !   if (ChannelsBottomLevel(i,j) < Me%ExtVar%BottomTopoG(i, j)) then
!             !       write (*,*) 
!             !       write (*,*) 'Bottom River section is lower than soil profile in cell', i,j
!             !       write (*,*) 'Increase bottom soil depth or decrease river depth'
!             !       stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'
!             !   endif
!                    
!                !Computing Area for flux - is always two times the lateral area and one time the bottom area
!                Area  = (2. * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))  &
!                            + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
!                
!                !Computing height gradient dH [m]
!                !Negative dh -> Flux from channels to porous media
!                !Positive dh -> Flux from porous media to channels 
!                dH    = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j)                
!                
!                !spatial step for Height Gradient (dH) [m]
!                dX = max(ChannelsBottomWidth(i, j) / 2.0, 1.0)
!                
!                ![m3/s]                   = [m/m] * [m2] * [m/s] * [] 
!                Me%lFlowToChannels(i, j) = (dH / dX ) * Area * Me%UnSatK(i, j, Me%UGCell(i,j))            &
!                                            * Me%SoilOpt%FCHCondFactor !Me%SoilOpt%HCondFactor
!
!                !If the channel looses water (infiltration), then set max flux so that volume in channel does not get 
!                !negative                
!                if (dH < 0) then
!                    ![m3]
!                    ChannelsVolume     = (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) * Area
!                    InfiltrationVolume = -1. * Me%lFlowToChannels(i, j) * Me%ExtVar%DT                
!                    
!                    if (InfiltrationVolume > ChannelsVolume) then
!                        Me%lFlowToChannels(i, j) = -0.5 * ChannelsVolume / Me%ExtVar%DT
!                        !write(*,*)'FlowToChannels corrected - ModulePorousMedia'
!                    endif
!                
!                !If soil looses water set flow so that cell stays at least with field theta
!                elseif (dH > 0) then
!                    k           = Me%UGCell(i,j)
!                    MinHead     = -0.5 * Me%ExtVar%DWZ(i, j, k)
!                    FieldTheta  = Theta_(MinHead, Me%SoilID(i,j,k))
!                    !m3/s   = (-) * m3 / s
!                    MaxFlow   = (Me%RC%ThetaS (i,j,k) - FieldTheta) * Me%ExtVar%CellVolume(i,j,k) / Me%ExtVar%DT
!
!                    if (Me%lFlowToChannels(i,j) > MaxFlow) then
!                        Me%lFlowToChannels(i,j) = MaxFlow
!                    endif
!                endif
!                
!
!            else
!
!                Me%lFlowToChannels(i, j) = 0.0
!
!            endif
!
!        enddo
!        enddo
!        
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR010'
!
!        call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR11'
!
!    end subroutine ExchangeWithDrainageNet_Corr
!
!    !--------------------------------------------------------------------------
!
!   subroutine ExchangeWithDrainageNet_Orig
!
!       !Arguments-------------------------------------------------------------
!       
!       !Local-----------------------------------------------------------------
!       integer                                     :: i, j, k, STAT_CALL
!       real                                        :: dH, Area
!       real                                        :: MinHead, FieldTheta, MaxVolume
!       real,   dimension(:, :), pointer            :: ChannelsWaterLevel
!       real,   dimension(:, :), pointer            :: ChannelsBottomLevel
!       real,   dimension(:, :), pointer            :: ChannelsBottomWidth
!       real,   dimension(:, :), pointer            :: ChannelsNodeLength
!       integer,dimension(:, :), pointer            :: ChannelsOpenProcess
!
!       call GetChannelsWaterLevel  (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR01'
!
!       call GetChannelsBottomLevel (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR02'
!
!       call GetChannelsBottomWidth (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR03'
!
!       call GetChannelsOpenProcess (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR04'
!
!       call GetChannelsNodeLength  (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR05'
!
!       do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!       do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!
!           if (Me%ExtVar%RiverPoints(i, j) == OpenPoint) then
!           
!           
!               !Channel will loose water
!               if (ChannelsWaterLevel(i, j) > Me%UGWaterLevel2D(i, j)) then
!
!                   !The Area is always two times the lateral area and
!                   !one time the bottom area
!                   Area  = (2. * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j))  &
!                                + ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
!
!                   if (Me%UGWaterLevel2D(i, j) > ChannelsBottomLevel(i, j)) then
!
!                       !Negative dh -> Flux from channels to porous media
!                       dH    = Me%UGWaterLevel2D(i, j) - ChannelsWaterLevel(i, j)
!
!                   else
!
!                       !Negative dh -> Flux from channels to porous media
!                       dH    = ChannelsBottomLevel(i, j) - ChannelsWaterLevel(i, j)
!
!                   endif
!
!               else
!           
!                   !REVER AQUI
!                   if (Me%UGWaterLevel2D(i, j) - Me%ExtVar%BottomTopoG(i, j) > MinUGThickness) then
!
!                       !Postive dh -> Flux from porous media to channels
!                       !This logical has to be rechecked - What happens if channels are deeper then soil layer??????
!                       dH   = Me%UGWaterLevel2D(i, j) - max(ChannelsWaterLevel(i, j), Me%ExtVar%BottomTopoG(i, j))
!                       !2* Lateral Area + Bottom Area
!                       Area = (2. * (Me%UGWaterLevel2D(i, j) - ChannelsBottomLevel(i, j)) + &
!                               ChannelsBottomWidth(i, j)) * ChannelsNodeLength(i, j)
!
!                   else
!
!                       dH   = 0.0
!                       Area = 0.0
!
!                   endif
!
!               endif
!
!               !Positive Flow to channel when Channel "gains" water
!               ![m3/s]                  [m]  [m2]   [m/s]                     [m]
!               !Me%lFlowToChannels(i, j) = dH * Area * Me%SoilOpt%HCondFactor * Me%UnSatK(i, j, Me%UGCell(i,j)) / &
!               !                           max(ChannelsBottomWidth(i, j) / 2.0, 1.0)               
!               Me%lFlowToChannels(i, j) = dH * Area * Me%SoilOpt%FCHCondFactor * Me%UnSatK(i, j, Me%UGCell(i,j)) / &
!                                          max(ChannelsBottomWidth(i, j) / 2.0, 1.0)
!               
!               
!               !Me%lFlowToChannels(i, j) = Area * Me%SatK(i, j, Me%UGCell(i,j))
!
!               !If the channel looses water (infiltration), then set max flux so that volume in channel does not get
!               !negative
!               if (dH < 0) then
!                   if (-1. * Me%lFlowToChannels(i, j) * Me%ExtVar%DT >            &
!                   (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) * Area) then
!                       Me%lFlowToChannels(i, j) = -0.5 * (ChannelsWaterLevel(i, j) - ChannelsBottomLevel(i, j)) &
!                       * Area / Me%ExtVar%DT
!                       write(*,*)'FlowToChannels corrected - ModulePorousMedia'
!                   endif
!               endif
!               
!               !If soil looses water set flow so that cell stays at least with field theta
!               !THIS DOESN'T WORK
!               if (dH > 0) then
!                   k           = Me%UGCell(i,j)
!                   MinHead     = -0.5 * Me%ExtVar%DWZ(i, j, k)
!                   FieldTheta  = Theta_(MinHead, Me%SoilID(i,j,k))
!                   MaxVolume   = (Me%RC%ThetaS (i,j,k) - Me%Theta(i,j,k)) * Me%ExtVar%CellVolume(i,j,k)
!                   Me%lFlowToChannels(i,j) = min(Me%lFlowToChannels(i,j), MaxVolume / Me%ExtVar%DT)
!               endif
!               
!
!           else
!
!               Me%lFlowToChannels(i, j) = 0.0
!
!           endif
!
!       enddo
!       enddo
!       
!       call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsWaterLevel, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR06'
!
!       call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomLevel, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR07'
!
!       call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsBottomWidth, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR08'
!
!       call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsOpenProcess, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR09'
!
!       call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsNodeLength, STAT = STAT_CALL)
!       if (STAT_CALL /= SUCCESS_) stop 'ExchangeWithDrainageNetwork - ModulePorousMedia - ERR10'
!
!   end subroutine ExchangeWithDrainageNet_Orig

   !--------------------------------------------------------------------------

    function BuckinghamDarcyEquation(con, hinf, hsup, delta)
    real ::  BuckinghamDarcyEquation

        !Local-------------------------------------------------------------------
        real, intent(IN) :: con
        real, intent(IN) :: hinf,  hsup
        real, intent(IN) :: delta
        !------------------------------------------------------------------------

        BuckinghamDarcyEquation = -1.0 * con * ((hsup - hinf) / delta )
                       
        !------------------------------------------------------------------------

    end function BuckinghamDarcyEquation

    !----------------------------------------------------------------------------

    subroutine IntegrationOutput
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, hor_i
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePointer       
        real, dimension(:,:,:), pointer             :: Modulus3D
        real, dimension(:,:,:), pointer             :: CenterU3D, CenterV3D, CenterW3D
        real, dimension(:,:), pointer               :: SurfaceSlice
        real, dimension(:,:), pointer               :: Output2D, ThetaFOutput
        real, dimension(:,:,:), pointer             :: Output3D
        type (T_IntegrationInfo), pointer           :: Info
        type (T_IntegrationByHorizon), pointer      :: Horizon   
        real                                        :: aux, dwz

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
           
        if (Me%ExtVar%Now >= Me%IntegrationOutput%OutTime(Me%IntegrationOutput%NextOutPut)) then

            !----------------------------------------------------------------------------
            !Common Output---------------------------------------------------------------
            !----------------------------------------------------------------------------

            !Time
            call ExtractDate (Me%ExtVar%Now , AuxTime(1), AuxTime(2),           &
                                              AuxTime(3), AuxTime(4),           &
                                              AuxTime(5), AuxTime(6))
            TimePointer => AuxTime

            call HDF5SetLimits  (Me%ObjIntegrationHDF5, 1, 6, STAT = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR010'


            call HDF5WriteData (Me%ObjIntegrationHDF5, "/Time", "Time",         &
                                "YYYY/MM/DD HH:MM:SS",                          &
                                Array1D      = TimePointer,                     &
                                OutputNumber = Me%IntegrationOutput%NextOutPut, &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR020'


            !Limits 
            call HDF5SetLimits (Me%ObjIntegrationHDF5, ILB, IUB, JLB, JUB, KLB-1, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR030'


            !Vertical 
            call HDF5WriteData (Me%ObjIntegrationHDF5,  "/Grid/VerticalZ",      &
                                "Vertical",   "m"              ,                &
                                Array3D      = Me%ExtVar%SZZ  ,                 &
                                OutputNumber = Me%IntegrationOutput%NextOutPut, &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR040'


            !Limits 
            call HDF5SetLimits (Me%ObjIntegrationHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR050'

            !Open Points
            call HDF5WriteData (Me%ObjIntegrationHDF5,  "/Grid/OpenPoints" ,    &
                                "OpenPoints", "-"                 ,             &
                                Array3D      = Me%ExtVar%OpenPoints3D ,         &
                                OutputNumber = Me%IntegrationOutput%NextOutPut, &
                                STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR060'

            allocate (Output2D (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate (Output3D (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
            
            !----------------------------------------------------------------------------
            !Bottom boundary output------------------------------------------------------
            !----------------------------------------------------------------------------            
            if (Me%IntegrationOutput%BoundaryBottom%Yes) then
                Info => Me%IntegrationOutput%BoundaryBottom
                             
                call HDF5WriteData (Me%ObjIntegrationHDF5, "/Results/bottom boundary volume",   &
                                    "bottom boundary volume", "m3",                             &
                                    Array2D      = Output2D,                                    &
                                    OutputNumber = Me%IntegrationOutput%NextOutPut,             &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR070'

                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%ExtVar%BasinPoints (i, j) == 1) then
                        Output2D (i,j) = Info%Field2D (i,j) / Me%IntegrationOutput%AccTime
                    else
                        Output2D (i,j) = null_real
                    endif
                enddo
                enddo                  
                !UGWaterDepth
                call HDF5WriteData (Me%ObjIntegrationHDF5, "/Results/bottom boundary Flow",     &
                                    "bottom boundary Flow", "m3/s",                             &
                                    Array2D      = Output2D,                                    &
                                    OutputNumber = Me%IntegrationOutput%NextOutPut,             &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR080'
                
                Info%Field2D = 0.0                                
                
            endif
            
            !----------------------------------------------------------------------------
            !Saturated output------------------------------------------------------------
            !----------------------------------------------------------------------------
            if (Me%IntegrationOutput%WaterTable%Yes) then    
                Info => Me%IntegrationOutput%WaterTable
                
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%ExtVar%BasinPoints (i, j) == 1) then
                        Output2D (i,j) = min (Info%Field2D (i,j) / Me%IntegrationOutput%AccTime, Me%ExtVar%Topography (i,j))
                    else
                        Output2D (i,j) = null_real
                    endif
                enddo
                enddo                
                !UGWaterLevel                
                call HDF5WriteData (Me%ObjIntegrationHDF5, "/Results/water level",          &
                                    "water level", "m",                                     &
                                    Array2D      = Output2D,                                &
                                    OutputNumber = Me%IntegrationOutput%NextOutPut,         &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR090'

                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%ExtVar%BasinPoints (i, j) == 1) then
                        Output2D (i,j) = max (Me%ExtVar%Topography (i,j) - Output2D (i,j), 0.0)
                    else
                        Output2D (i,j) = null_real
                    endif
                enddo
                enddo                  
                !UGWaterDepth
                call HDF5WriteData (Me%ObjIntegrationHDF5, "/Results/water table depth",    &
                                    "water table depth", "m",                               &
                                    Array2D      = Output2D,                                &
                                    OutputNumber = Me%IntegrationOutput%NextOutPut,         &
                                    STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR100'

                Info%Field2D = 0.0                               
                
            endif

            !----------------------------------------------------------------------------
            !Unsaturated output----------------------------------------------------------
            !----------------------------------------------------------------------------
                       
            Info => Me%IntegrationOutput%WaterContent
            if (Info%Yes) then
                if (Info%ByLayer) then

                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB
                        if (Me%ExtVar%WaterPoints3d (i,j,k) == 1) then
                            Output3D (i,j,k) = Info%Field3D (i,j,k) / Me%IntegrationOutput%AccTime
                        else
                            Output3D (i,j,k) = null_real
                        endif
                    enddo
                    enddo
                    enddo
                    
                    call HDF5WriteData (Me%ObjIntegrationHDF5, "/Results/water content",    &
                                        "water content", "m3water/m3soil",                  &
                                        Array3D      = Output3D,                            &
                                        OutputNumber = Me%IntegrationOutput%NextOutPut,     &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR110'
                
                    if (Info%ByHorizon) then
                    
                        do hor_i = 1, Info%HorizonsCount
                            Horizon => Info%Horizons(hor_i)
                        
                            do j = JLB, JUB
                            do i = ILB, IUB
                                aux = 0.0
                                dwz = 0.0
                                
                                if (Me%ExtVar%BasinPoints(i,j) == 1) then
                                    
                                    do k = Horizon%StartLayer, Horizon%EndLayer
                                        if (Me%ExtVar%WaterPoints3D (i,j,k) == 1) then
                                            aux = aux + Output3D (i,j,k) * Me%ExtVar%DWZ(i,j,k)
                                            dwz = dwz + Me%ExtVar%DWZ(i, j, k)
                                        endif                                
                                    enddo
                            
                                    if (dwz > 0.0) then
                                        Output2D (i,j) = aux / dwz
                                    else
                                        Output2D (i,j) = null_real
                                    endif
                                else
                                    Output2D (i,j) = null_real
                                endif                            
                            enddo
                            enddo
                    
                            call HDF5WriteData (Me%ObjIntegrationHDF5,                                  &
                                                "/Results/water content ["//trim(Horizon%Name)//"]",    &
                                                "water content", "m3water/m3soil",                      &
                                                Array2D      = Output2D,                                &
                                                OutputNumber = Me%IntegrationOutput%NextOutPut,         &
                                                STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR120'                        
                        enddo
                    endif                     
                    
                    Info%Field3D = 0.0
                    
                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB
                        if (Me%ExtVar%WaterPoints3d (i,j,k) == 1) then
                            Output3D (i,j,k) = ThetaF_ (Output3D (i,j,k), Me%SoilID(i,j,k))
                        else
                            Output3D (i,j,k) = null_real
                        endif
                    enddo
                    enddo
                    enddo
                        
                    call HDF5WriteData (Me%ObjIntegrationHDF5, "/Results/relative water content",   &
                                        "relative water content", "m3water/m3soil",                 &
                                        Array3D      = Output3D,                                    &
                                        OutputNumber = Me%IntegrationOutput%NextOutPut,             &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR130'     
                    
                    if (Info%ByHorizon) then
                    
                        do hor_i = 1, Info%HorizonsCount
                            Horizon => Info%Horizons(hor_i)
                        
                            do j = JLB, JUB
                            do i = ILB, IUB
                                aux = 0.0
                                dwz = 0.0
                                
                                if (Me%ExtVar%BasinPoints(i,j) == 1) then
                                    
                                    do k = Horizon%StartLayer, Horizon%EndLayer
                                        if (Me%ExtVar%WaterPoints3D (i,j,k) == 1) then
                                            aux = aux + Output3D (i,j,k) * Me%ExtVar%DWZ(i,j,k)
                                            dwz = dwz + Me%ExtVar%DWZ(i, j, k)
                                        endif                                
                                    enddo
                            
                                    if (dwz > 0.0) then
                                        Output2D (i,j) = aux / dwz
                                    else
                                        Output2D (i,j) = null_real
                                    endif
                                else
                                    Output2D (i,j) = null_real
                                endif                            
                            enddo
                            enddo
                    
                            call HDF5WriteData (Me%ObjIntegrationHDF5,                                          &
                                                "/Results/relative water content ["//trim(Horizon%Name)//"]",   &
                                                "relative water content", "m3water/m3soil",                     &
                                                Array2D      = Output2D,                                        &
                                                OutputNumber = Me%IntegrationOutput%NextOutPut,                 &
                                                STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR140'                        
                        enddo
                    endif                     
                endif                                
            endif            
            
            Me%IntegrationOutput%NextOutPut = Me%IntegrationOutput%NextOutPut + 1
            Me%IntegrationOutput%AccTime = 0.0

            !----------------------------------------------------------------------------
            !Write everything to disk----------------------------------------------------
            !----------------------------------------------------------------------------

            call HDF5FlushMemory (Me%ObjIntegrationHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'IntegrationOutput - ModulePorousMedia - ERR0150'
            
            deallocate (Output2D)
            deallocate (Output3D)

        endif
    
    end subroutine IntegrationOutput

    !----------------------------------------------------------------------------
    
    subroutine PorousMediaOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePointer       
        real, dimension(:,:,:), pointer             :: Modulus3D
        real, dimension(:,:,:), pointer             :: CenterU3D, CenterV3D, CenterW3D
        real, dimension(:,:), pointer               :: SurfaceSlice

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !NOTE: This test below is only valid if either normal output or 2D is active. 
        !That way during startup it is checked if only one of the output types is activated
        !Both can't be active at the same time beause of the TIME output

        if (Me%OutPut%Yes) then
            
            !3D "Normal" Output
            if (Me%ExtVar%Now >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

                !----------------------------------------------------------------------------
                !Common Output---------------------------------------------------------------
                !----------------------------------------------------------------------------

                !Time
                call ExtractDate   (Me%ExtVar%Now , AuxTime(1), AuxTime(2),         &
                                                    AuxTime(3), AuxTime(4),         &  
                                                    AuxTime(5), AuxTime(6))
                TimePointer => AuxTime

                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR01'


                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                     "YYYY/MM/DD HH:MM:SS",                         &
                                     Array1D      = TimePointer,                    &
                                     OutputNumber = Me%OutPut%NextOutPut,           &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR02'


                !Limits 
                call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB-1, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR03'


                !Vertical 
                call HDF5WriteData  ( Me%ObjHDF5,  "/Grid/VerticalZ",               & 
                                     "Vertical",   "m"              ,               & 
                                      Array3D      = Me%ExtVar%SZZ  ,               &
                                      OutputNumber = Me%OutPut%NextOutPut,          &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR010'


                !Limits 
                call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR020'

                !Open Points
                call HDF5WriteData   ( Me%ObjHDF5,  "/Grid/OpenPoints" ,            &
                                      "OpenPoints", "-"                 ,           &
                                       Array3D      = Me%ExtVar%OpenPoints3D ,      &
                                       OutputNumber = Me%OutPut%NextOutPut,         &
                                       STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR030'


                !----------------------------------------------------------------------------
                !Saturated Output------------------------------------------------------------
                !----------------------------------------------------------------------------
                
                !UGWaterLevel
                call HDF5WriteData   (Me%ObjHDF5, "/Results/water level",           &
                                      "water level", "m",                           &
                                      Array2D      = Me%UGWaterLevel2D,             &
                                      OutputNumber = Me%OutPut%NextOutPut,          &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR040'

                !UGWaterDepth
                call HDF5WriteData   (Me%ObjHDF5, "/Results/water table depth",     &
                                      "water table depth", "m",                     &
                                      Array2D      = Me%UGWaterDepth2D,             &
                                      OutputNumber = Me%OutPut%NextOutPut,          &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR050'

                if (Me%SoilOpt%ImposeBoundary .and. Me%SoilOpt%ImposeBoundaryWalls) then
                    !Imposedlevel on boundary cells
                    call HDF5WriteData   (Me%ObjHDF5, "/Results/imposed boundary level",   &
                                          "imposed boundary level", "m",                   &
                                          Array2D      = Me%Boundary%ImposedBoundaryLevel, &
                                          OutputNumber = Me%OutPut%NextOutPut,             &
                                          STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR060'                
                endif

                !----------------------------------------------------------------------------
                !Unsaturated Output----------------------------------------------------------
                !----------------------------------------------------------------------------
                call HDF5WriteData ( Me%ObjHDF5,    "/Results/water content",       &
                                    "water content", "m3water/m3soil"            ,  &
                                     Array3D      =  GetPointer(Me%Theta) ,         &
                                     OutputNumber =  Me%OutPut%NextOutPut    ,      &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR070'
                
                call HDF5WriteData ( Me%ObjHDF5,    "/Results/relative water content",   &
                                    "relative water content", "m3water/m3water"      ,   &
                                     Array3D      =  GetPointer(Me%RC%ThetaF)  ,         &
                                     OutputNumber =  Me%OutPut%NextOutPut    ,      &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR080'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/Head"            ,       &
                                    "Head"     , 'm'                            ,   &
                                     Array3D      = GetPointer(Me%Head)   ,         &
                                     OutputNumber = Me%OutPut%NextOutPut        ,   & 
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR090'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/FinalHead"    ,          &
                                    "FinalHead"     , 'm'                       ,   &
                                     Array3D      = GetPointer(Me%FinalHead),       &
                                     OutputNumber = Me%OutPut%NextOutPut,           & 
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR0100'

                call HDF5WriteData ( Me%ObjHDF5, "/Results/HydroPressure"    ,   &
                                    "HydroPressure"     , 'm'                            ,   &
                                     Array3D      = GetPointer(Me%HydroPressure),   &
                                     OutputNumber = Me%OutPut%NextOutPut        ,   & 
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR0110'

                !Write unsat velocities
                allocate(CenterU3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB), &
                         CenterV3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB), &
                         CenterW3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB), &
                         Modulus3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                CenterU3D  = 0.0
                CenterV3D  = 0.0
                CenterW3D  = 0.0
                Modulus3D  = 0.0
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        CenterU3D (i, j, k) = (Me%UnsatVelU(i, j, k) + Me%UnsatVelU(i, j+1, k)) / 2.0
                        CenterV3D (i, j, k) = (Me%UnsatVelV(i, j, k) + Me%UnsatVelV(i+1, j, k)) / 2.0
!                        CenterW3D (i, j, k) = (Me%UnsatVelW(i, j, k) + Me%UnsatVelW(i, j, k+1)) / 2.0
                        CenterW3D (i, j, k) = (Me%UnsatVelWFinal(i, j, k) + Me%UnsatVelWFinal(i, j, k+1)) / 2.0
                        Modulus3D (i, j, k) = sqrt(CenterU3D (i, j, k)**2.0 + CenterV3D (i, j, k)**2.0 + CenterW3D(i,j,k)**2.0)
                    endif
                enddo
                enddo            
                enddo

                call HDF5WriteData (Me%ObjHDF5,                                          &
                                    "/Results/"//trim(GetPropertyName (VelocityU_)),     &
                                    trim(GetPropertyName (VelocityU_)),                  &
                                    "m/s",                                               &
                                     Array3D      = CenterU3D,                           &
                                     OutputNumber = Me%OutPut%NextOutPut,                &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR0120'
            
                call HDF5WriteData (Me%ObjHDF5,                                          &
                                    "/Results/"//trim(GetPropertyName (VelocityV_)),     &
                                    trim(GetPropertyName (VelocityV_)),                  &
                                    "m/s",                                               &
                                     Array3D      = CenterV3D,                      &
                                     OutputNumber = Me%OutPut%NextOutPut,           &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR0130'

                call HDF5WriteData (Me%ObjHDF5,                                          &
                                    "/Results/"//trim(GetPropertyName (VelocityW_)),     &
                                    trim(GetPropertyName (VelocityW_)),                  &
                                    "m/s",                                               &
                                     Array3D      = CenterW3D,                      &
                                     OutputNumber = Me%OutPut%NextOutPut,           &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR0140'

                call HDF5WriteData (Me%ObjHDF5,                                                &
                                    "/Results/"//trim(GetPropertyName (VelocityModulus_)),     &
                                    trim(GetPropertyName (VelocityModulus_)),                  &
                                    "m/s",                                               &
                                     Array3D      = Modulus3D,                      &
                                     OutputNumber = Me%OutPut%NextOutPut,           &  
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR0150'


                deallocate (CenterU3D, CenterV3D, CenterW3D, Modulus3D)

                Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

                !----------------------------------------------------------------------------
                !Write everything to disk----------------------------------------------------
                !----------------------------------------------------------------------------

                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR0160'


            endif
            
        endif
        
        !2D Surface output
        if (Me%OutPut%SurfaceOutput) then
            
            if (Me%ExtVar%Now >= Me%OutPut%SurfaceOutTime(Me%OutPut%NextSurfaceOutput)) then


                !----------------------------------------------------------------------------
                !Common Output---------------------------------------------------------------
                !----------------------------------------------------------------------------

                !Time
                call ExtractDate   (Me%ExtVar%Now , AuxTime(1), AuxTime(2),         &
                                                    AuxTime(3), AuxTime(4),         &  
                                                    AuxTime(5), AuxTime(6))
                TimePointer => AuxTime

                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR0170'


                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                     "YYYY/MM/DD HH:MM:SS",                         &
                                     Array1D      = TimePointer,                    &
                                     OutputNumber = Me%OutPut%NextSurfaceOutput,    &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR0180'


                !Limits 
                call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR0190'


                !Write unsat velocities
                allocate(SurfaceSlice(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                
                !Theta
                do j = JLB, JUB
                do i = ILB, IUB
                    SurfaceSlice(i, j) =  Me%RC%ThetaF(i, j, KUB)
                enddo            
                enddo                
                
                !Please don't change the name of the data sets, since they are used by 
                !MOHID Land Operational
                call HDF5WriteData (Me%ObjHDF5, "/Results/relative water content",                  &
                                    "relative water content",                                       &
                                    "m3water/m3water",                                              &
                                    Array2D      =  SurfaceSlice,                                   &
                                    OutputNumber =  Me%OutPut%NextSurfaceOutput,                    &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_Format - ModulePorousMedia - ERR0210'
                
                deallocate(SurfaceSlice)
                
                !Please don't change the name of the data sets, since they are used by 
                !MOHID Land Operational
                call HDF5WriteData   (Me%ObjHDF5, "/Results/water table depth",     &
                                      "water table depth", "m",                     &
                                      Array2D      = Me%UGWaterDepth2D,             &
                                      OutputNumber = Me%OutPut%NextSurfaceOutput,   &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR0220'


                Me%OutPut%NextSurfaceOutput = Me%OutPut%NextSurfaceOutput + 1
                

                !----------------------------------------------------------------------------
                !Write everything to disk----------------------------------------------------
                !----------------------------------------------------------------------------

                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaOutput - ModulePorousMedia - ERR0230'
    
            
            endif


        endif

       
    end subroutine PorousMediaOutput

    !--------------------------------------------------------------------------

    subroutine OutPutTimeSeries

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, STAT_CALL
        !Begin-----------------------------------------------------------------

        !Calculates real Infiltration Velocity
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == 1) then
                Me%InfiltrationVelocity (i, j) = - Me%Infiltration(i, j) / Me%ExtVar%DT
            endif
        enddo
        enddo

        !Theta
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%Theta),                              &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR010'

        !ThetaF
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%RC%ThetaF),                                      &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR020'

        !velocity for now......
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%UnsatVelW),                          &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR030'

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%UnsatVelWFinal),                     &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR040'

        !infiltration velocity for now
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%InfiltrationVelocity,                           &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR050'
        

        !Head
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%Head),                               &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR060'

        !Conductivity
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%UnSatK),                             &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR070'

        !Level Water Table
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%UGWaterLevel2D,                                 &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR080'

        !Depth Water Table
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data2D = Me%UGWaterDepth2D,                                 &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR090'

        !Hydro Pressure
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%HydroPressure),                      &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR100'

        !Final Head
        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D = GetPointer(Me%FinalHead),                          &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR110'

        !Porous Media Water column - to compare to Basin water column
!        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
!                            Data2D = PMWaterColumn,                                     &
!                            STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR091'

        if (Me%ObjDrainageNetwork /= 0) then
        
            !Flow to channels total
            call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                                Data2D = GetPointer(Me%iFlowToChannels),                    &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR120' 
                           
            if (Me%SoilOpt%DNLink == GWFlowToChanByLayer_) then
                !Flow to channels by layer
                call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                                    Data3D = GetPointer(Me%iFlowToChannelsLayer),               &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR130'
   
            endif
        endif
        
        if (Me%ExtVar%ConstructEvaporation) then
            !Evaporation
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D = Me%EvaporationFlux,                            &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR140'
        endif

        if (Me%ExtVar%ConstructTranspiration) then
            !Transpiration
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data3D = Me%ExtVar%TranspirationFlux,                   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModulePorousMedia - ERR150'        
        endif


    end subroutine OutPutTimeSeries

    !--------------------------------------------------------------------------

    subroutine ProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Theta
        call WriteProfile(Me%ObjProfile,                                        &
                          Data3D = GetPointer(Me%Theta),                        &
                          SZZ    = Me%ExtVar%SZZ,                               &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR010'

        !ThetaF
        call WriteProfile(Me%ObjProfile,                                        &
                          Data3D = GetPointer(Me%RC%ThetaF),                                &
                          SZZ    = Me%ExtVar%SZZ,                               &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR020'

        !Head
        call WriteProfile(Me%ObjProfile,                                        &
                          Data3D = GetPointer(Me%Head),                         &
                          SZZ    = Me%ExtVar%SZZ,                               &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR030'

    end subroutine ProfileOutput

    !--------------------------------------------------------------------------

    subroutine ComputeBoxesWaterFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, CHUNK, i, j, k
        real, dimension(:,:,:), pointer         :: WaterVolume
        real(8), dimension (:,:,:), pointer     :: FluxU, FluxV, FluxW    
        !----------------------------------------------------------------------
       
        if (MonitorPerformance) call StartWatch ("ModulePorousMedia", "ComputeBoxesWaterFluxes")
        
        !from allocatable to pointer to fit to BoxDif 
        FluxU => Me%FluxU
        FluxV => Me%FluxV
        FluxW => Me%FluxWFinal
        
        call BoxDif(Me%ObjBoxDif,                                                    &
                    FluxU,                                                           &
                    FluxV,                                                           &
                    FluxW,                                                           & 
                    'soil_water',                                                    &
                    Me%ExtVar%WaterPoints3D,                                         &
                    STAT = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                 &
           stop 'Subroutine ComputeBoxesWaterFluxes - ModulePorousMedia. ERR01'
        
        nullify (FluxU)
        nullify (FluxV)
        nullify (FluxW)
        
        allocate(WaterVolume(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB, Me%WorkSize%KLB:Me%WorkSize%KUB))
        WaterVolume = null_real
        
        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                
                ![m3] = [m3H20/m3soil] * [m3soil]
                WaterVolume(i,j,k) = Me%Theta(i,j,k) * Me%ExtVar%CellVolume(i,j,k)
                
            endif

        enddo
        enddo
        enddo      
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        
        call BoxDif(Me%ObjBoxDif,                                                    &
                    WaterVolume,                                                     &
                    'soil_water',                                                    &
                    Me%ExtVar%WaterPoints3D,                                         &
                    STAT = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                 &
           stop 'Subroutine ComputeBoxesWaterFluxes - ModulePorousMedia. ERR02'

        deallocate (WaterVolume)
        
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMedia", "ComputeBoxesWaterFluxes")

    end subroutine ComputeBoxesWaterFluxes

    !--------------------------------------------------------------------------

    subroutine CalculateTotalStoredVolume(WriteOut)

        !Arguments-------------------------------------------------------------
        logical, optional                           :: WriteOut !Debug only

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        
        Me%TotalStoredVolume = 0.0

        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
            if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then

                Me%TotalStoredVolume = Me%TotalStoredVolume + Me%Theta(i,j,k)             * &
                                       Me%ExtVar%CellVolume(i,j,k)
            endif

        enddo
        enddo
        enddo

        if (present(WriteOut)) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%BasinPoints(i, j) == 1) then
                    Me%TotalStoredVolume = Me%TotalStoredVolume + Me%WaterColumn(i, j)* Me%ExtVar%Area(i, j)
                endif
            enddo
            enddo

            write(*,*)Me%TotalStoredVolume
        endif
      
    end subroutine CalculateTotalStoredVolume

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillPorousMedia(ObjPorousMediaID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjPorousMediaID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, STAT_CALL           
        type(T_Piezometer), pointer         :: Piezometer

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mPorousMedia_,  Me%InstanceID)

            if (nUsers == 0) then
                
                !Write Output for continuous computation
                call WriteFinalFile

                !Kills the TimeSerie
                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR01'
                endif

                !Kills Profile Output
                if (Me%OutPut%ProfileON) then
                    call KillProfile(Me%ObjProfile, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR010'
                endif
                                    
                !Deassociates External Instances
                if (Me%ObjDrainageNetwork /= 0) then
                    nUsers = DeassociateInstance (mDRAINAGENETWORK_, Me%ObjDrainageNetwork)
                    if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR020'
                endif                
                
                if (Me%ImpermeableFractionID%ObjFillMatrix /= 0) then
                    call KillFillMatrix (Me%ImpermeableFractionID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillPorousMedia - PorousMedia - ERR030'                
                endif


                if (Me%Output%BoxFluxes) then
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillPorousMedia - PorousMedia - ERR35'
                endif
                
                if (Me%SoilOpt%ImposeBoundary) then
                    
                    Piezometer => Me%Boundary%FirstPiezometer

                    do while(associated(Piezometer))  

                        if (Piezometer%TimeSerie%ObjTimeSerie /= 0) then
                            call KillTimeSerie(Piezometer%TimeSerie%ObjTimeSerie, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillPorousMedia - ModulePorousMedia - ERR40' 
                        endif

                        Piezometer => Piezometer%Next
                                               
                    enddo
                endif
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMedia - ERR560'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMedia - ERR70'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjTopography)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMedia - ERR80'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMedia - ERR90'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMedia - ERR100'
                
                call KillMap (Me%ObjMap, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - ModulePorousMedia - ERR110'

                call KillGeometry (Me%Objgeometry, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - ModulePorousMedia - ERR120'

                call KillGridData (Me%ObjBottomTopography, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - ModulePorousMedia - ERR130'                
                
                if (Me%OutPut%Yes .or. Me%Output%SurfaceOutput) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - ModulePorousMedia - ERR140'
                endif

                if (Me%IntegrationOutput%Yes) then
                    call KillHDF5 (Me%ObjIntegrationHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - ModulePorousMedia - ERR141'
                endif

                !Deallocates Instance
                call DeallocateInstance ()

                ObjPorousMediaID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
                   
    end subroutine KillPorousMedia        

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

        if (Me%ExtVar%Now == Me%EndTime) then
            FileName = Me%Files%FinalFile
        else
            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%ExtVar%Now))//".fin")
        endif            


        call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFileOld - ModulePorousMedia - ERR01'

        open(Unit = FinalFile, File = FileName, Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFileOld - ModulePorousMedia - ERR02'

        !Writes Date
        call ExtractDate(Me%ExtVar%Now, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)
        write(FinalFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File

        write(FinalFile)Me%Theta

        call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFileOld - ModulePorousMedia - ERR03'
        

    end subroutine WriteFinalFile

    !------------------------------------------------------------------------    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_PorousMedia), pointer          :: AuxObjPorousMedia
        type (T_PorousMedia), pointer          :: PreviousObjPorousMedia

        !Updates pointers
        if (Me%InstanceID == FirstObjPorousMedia%InstanceID) then
            
            FirstObjPorousMedia => FirstObjPorousMedia%Next
        
        else

            PreviousObjPorousMedia => FirstObjPorousMedia
            AuxObjPorousMedia      => FirstObjPorousMedia%Next
            
            do while (AuxObjPorousMedia%InstanceID /= Me%InstanceID)                
                PreviousObjPorousMedia => AuxObjPorousMedia
                AuxObjPorousMedia      => AuxObjPorousMedia%Next            
            enddo

            !Now update linked list
            PreviousObjPorousMedia%Next => AuxObjPorousMedia%Next

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

    subroutine Ready (ObjPorousMedia_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMedia_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

        if (ObjPorousMedia_ID > 0) then
            call LocateObjPorousMedia (ObjPorousMedia_ID)
            ready_ = VerifyReadLock (mPorousMedia_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjPorousMedia (ObjPorousMediaID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaID        

        Me => FirstObjPorousMedia
        do while (associated (Me))
            if (Me%InstanceID == ObjPorousMediaID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModulePorousMedia - LocateObjPorousMedia - ERR01'

    end subroutine LocateObjPorousMedia

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar                

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL


        !Time------------------------------------------------------------------
        
        call GetComputeCurrentTime  (Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR01'

        call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR02'
        


        !Topography------------------------------------------------------------

        call GetGridData  (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR03'


        !BottomTopography------------------------------------------------------
        
        call GetGridData  (Me%ObjBottomTopography, Me%ExtVar%BottomTopoG, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR04'



        !Basin Geometry--------------------------------------------------------                        
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR05'  

        call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR06'


        !Drainage network
        
                      

        !Horizontal Grid-------------------------------------------------------
        
        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                    &
                               DUX = Me%ExtVar%DUX, DVY = Me%ExtVar%DVY,                &
                               DZX = Me%ExtVar%DZX, DZY = Me%ExtVar%DZY,                &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR07'
        

        call GetGridCellArea (Me%ObjHorizontalGrid,                                     & 
                              GridCellArea = Me%ExtVar%Area,                            & 
                              STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR08'



        !Geometry--------------------------------------------------------------

        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  SZZ         = Me%ExtVar%SZZ,                          &
                                  DZZ         = Me%ExtVar%DZZ,                          &
                                  DWZ         = Me%ExtVar%DWZ,                          &
                                  ZCellCenter = Me%ExtVar%CenterCell,                   &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR09")


        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ    = Me%ExtVar%CellVolume,                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR010")


        call GetGeometryKFloor(Me%ObjGeometry,                                          &
                               Z    = Me%ExtVar%KFloor,                                 &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR11")

        call GetGeometryAreas(Me%ObjGeometry,                                           &
                              Me%ExtVar%AreaU, Me%ExtVar%AreaV,                         &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR12")

!        call GetGeometryAreas(Me%ObjGeometry,                                           &
!                              Me%ExtVar%AreaV, Me%ExtVar%AreaV,                         &
!                              STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                      &
!            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR13")

        


        !Map-------------------------------------------------------------------
        
        call GetWaterPoints3D(Me%ObjMap,                                                &
                              Me%ExtVar%WaterPoints3D,                                  &
                              STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR14")


        call GetOpenPoints3D(Me%ObjMap,                                                 &
                             Me%ExtVar%OpenPoints3D,                                    &
                             STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR15")


        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesU3D = Me%ExtVar%ComputeFacesU3D,             &
                               ComputeFacesV3D = Me%ExtVar%ComputeFacesV3D,             &
                               ComputeFacesW3D = Me%ExtVar%ComputeFacesW3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModulePorousMedia. ERR16")


        !Uderground HorizontalMap----------------------------------------------


    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar                

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        
        !Topography------------------------------------------------------------

        call UngetGridData (Me%ObjTopography, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR01'


        !BottomTopography------------------------------------------------------
        
        call UngetGridData (Me%ObjBottomTopography, Me%ExtVar%BottomTopoG, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR02'

        
        !Basin Geometry--------------------------------------------------------                        
          
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR03'

        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR04'               


        !Horizontal Grid-------------------------------------------------------
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR06'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR09'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR10'
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%Area, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMedia - ERR11'               
        


        !Geometry--------------------------------------------------------------

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%SZZ,          STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR12")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%DZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR14") 

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%DWZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR15") 

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%CenterCell,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR15a") 


        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%CellVolume,   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR16")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%KFloor,       STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR18")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%AreaU,        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR19")

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%AreaV,       STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR20")
        

        
        !Map-------------------------------------------------------------------
        
        call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR21") 
        
        call UnGetMap(Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR22") 


        call UnGetMap(Me%ObjMap, Me%ExtVar%ComputeFacesU3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR23") 


        call UnGetMap(Me%ObjMap, Me%ExtVar%ComputeFacesV3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR24") 


        call UnGetMap(Me%ObjMap, Me%ExtVar%ComputeFacesW3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModulePorousMedia. ERR25")         


    end subroutine ReadUnLockExternalVar
    
end module ModulePorousMedia

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 

