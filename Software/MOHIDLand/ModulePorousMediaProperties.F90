!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land
! MODULE        : PorousMediaProperties
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Fev 2010 David Brito
! REVISION      : 
! DESCRIPTION   : Module to handle Properties in PorousMedia
!
!------------------------------------------------------------------------------
!
!
!Units in porous media properties
!   Transported properties (soluble)  : g/m3 (or mg/L)                     (needs to convert concentrations to SedimentQuality and 
!                                                                           PREEQC at entrance and exit)
!   Adsorbed properties (non soluble) : mg/kgsoil                          (needs to convert concentrations to PREEQC )
!   Solid Phases (non soluble)        : mg/kgsoil                          (used in PhreeqC to make equilibrium with soil solution)
!   Gas Phases (non soluble)          : mol (mass); atm (partial pressure) (used in PhreeqC to make equilibrium with soil solution)
!   Gas Phases (assumed soluble)      : mg/L                               (used in SedimentQuality - N2, CH4 and CO2)
!
!   Soil Dry Density                  : kg/m3                              
!   H+ and Ionic Strenght             : mol/L                              (used in SedimentQuality)
!   Microorganisms Population         : #org/kgsoil                        (used in SedimentQuality)
!
! ADVDIFF_EXPLICIT              : 0/1                [1]        !REMARK: Horizontal diffusion is always explicit
!                                                               !(1 - adv and diff are explicit in all directions; 0 - adv and diff 
!                                                               !are implicit in vertical, horizontal adv may be explicit or impl) 
!   ADVDIFF_ADVECTION_H_IMP_EXP : 0/1                [1]        !(read if ADVDIFF_EXPLICIT : 0; 0 - horiz adv implicit; 
!                                                                1 - horiz adv explicit)
! NEW_FORMULATION               : 0/1                [0]        !if 1 then spatial methods will be the same for all properties
!     ADVDIFF_METHOD_H          : integer      [UpwindOrder1]   !Spatial methods for horizontal advection
!                                                               !UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4,
!                                                                CentralDif = 5, LeapFrog = 6    
!     ADVDIFF_METHOD_V          : integer      [UpwindOrder1]   !Spatial methods for vertical advection
!                                                               !UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4,
!                                                                CentralDif = 5, LeapFrog = 6!
! <beginproperty>
!   ADVECTION_DIFFUSION         : 0/1               [0]         !Property advection - diffusion
!       ADVDIFF_METHOD_H        : integer      [UpwindOrder1]   !Spatial methods for horizontal advection
!                                                               !UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4,
!                                                                CentralDif = 5, LeapFrog = 6    
!       ADVDIFF_METHOD_V        : integer      [UpwindOrder1]   !Spatial methods for vertical advection
!                                                               !UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4,
!                                                                CentralDif = 5, LeapFrog = 6
!       ADVDIFF_TVD_LIMIT_H     : integer        [Superbee]     !Horizontal advection non-linear stability conditions
!                                                                MinMod = 1, VanLeer = 2, Muscl = 3, Superbee = 4, PDM = 5
!       ADVDIFF_TVD_LIMIT_V     : integer        [Superbee]     !Vertical advection non-linear stability conditions
!                                                               !MinMod = 1, VanLeer = 2, Muscl = 3, Superbee = 4, PDM = 5
!       ADVDIFF_VOLUME_RELATION_MAX : real          5.          !The relation between adjacent volumes above which 
!                                                               !the advection is upwind

!   SOIL_CHEMISTRY              : 0/1               [0]         !Use PREEQC model to change property (source/sink model)
!   SOIL_QUALITY                : 0/1               [0]         !Use SedimentQuality model to change property (source/sink model)
!   PARTITION                   : 0/1               [0]         !Compute partition between dissolved-particulate phases
!       PARTITION_COUPLE        : char               +          !Name of the property (oposite phase) to compute partition
!       PARTITION_FRACTION      : real               -          !Percentage of mass of a property in a determined phase 
!       PARTITION_RATE          : real            [1 s-1]       !Kinetic rate of partition to reach equilibrium
! <endproperty>
!
!------------------------------------------------------------------------------

Module ModulePorousMediaProperties

    use ModuleGlobalData
    use ModuleStopWatch
    use ModuleFunctions
    use ModuleTime
    use ModuleHDF5
    use ModuleEnterData
    use ModuleProfile,            only : StartProfile, WriteProfile, KillProfile
    use ModuleGridData,           only : ConstructGridData, GetGridData, UngetGridData,    &
                                         KillGridData, WriteGridData
    use ModuleTimeSerie,          only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,    &
                                         GetNumberOfTimeSeries, GetTimeSerieLocation,      &
                                         TryIgnoreTimeSerie, CorrectsCellsTimeSerie,       &
                                         GetTimeSerieName
    use ModuleHorizontalGrid,     only : GetHorizontalGrid, GetGridCellArea,               &
                                         WriteHorizontalGrid, UnGetHorizontalGrid, GetXYCellZ
    use ModuleBasinGeometry,      only : GetBasinPoints, GetRiverPoints,  UnGetBasin 
                                       
    use ModuleFillMatrix,         only : ConstructFillMatrix, GetDefaultValue,             &
                                         KillFillMatrix, ModifyFillMatrix
    use ModuleGeometry
    use ModuleMap
    use ModuleBoxDif,             only : StartBoxDif, GetBoxes, GetNumberOfBoxes, UngetBoxDif,     &
                                          BoxDif, KillBoxDif        
    use ModulePorousMedia,        only : GetOldWaterContent, GetWaterContent, GetFluxU,    &
                                         GetFluxV, GetFluxW, GetUnsatW, GetUnsatV,         &
                                         GetUnsatU, UnGetPorousMedia, GetUnsatWFinal,      &
                                         GetThetaS, GetGWFlowToChannels, GetThetaF,        &
                                         GetGWLayer, GetGWLayerOld, GetPotentialInfiltration, &
                                         GetGWFlowToChannelsByLayer, GetGWToChannelsLayers,&
                                         GetGWFlowOption, GetTranspiration, GetEvaporation,&
                                         GetBoundaryImposed, GetBoundaryCells,             &
                                         GetBoundaryFluxWalls, GetBoundaryFluxBottom,      &
                                         GetFlowDischarge
                                         
    use ModuleChainReactions,     only : StartChainReactions, SetCRPropertyConcentration,  &
                                         GetCRPropertiesList, UnGetChainReactions,         &
                                         ModifyChainReactions, KillChainReactions,         &
                                         InitCRSoilPhase
                                         
#ifdef _PHREEQC_                                          
    use ModuleInterface,          only : ConstructInterface, Modify_Interface, GetPhreeqCID, &
                                         GetRateFlux 
#else
    use ModuleInterface,          only : ConstructInterface, Modify_Interface, GetRateFlux     
#endif
                                         
    use ModuleAdvectionDiffusion, only : StartAdvectionDiffusion, AdvectionDiffusion,      &
                                         GetAdvFlux, GetDifFlux, GetBoundaryConditionList, &
                                         UngetAdvectionDiffusion, KillAdvectionDiffusion

#ifdef _PHREEQC_     
    use ModulePhreeqC
#endif

#ifdef _ENABLE_CUDA
    use ModuleCuda
#endif _ENABLE_CUDA

    use ModuleDischarges        ,only : Construct_Discharges, GetDischargesNumber,       &
                                        GetDischargesGridLocalization,                   &
                                        GetDischargeWaterFlow, GetDischargesIDName,      &
                                        TryIgnoreDischarge, GetDischargeSpatialEmission, &
                                        CorrectsCellsDischarges, Kill_Discharges,        &
                                        GetDischargeConcentration,                       & 
                                        GetDischargeFlowDistribuiton, UngetDischarges

    !griflet
    !$ use omp_lib

   implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !PorousMediaProperties
    public  :: ConstructPorousMediaProperties
    private ::      AllocateInstance
    private ::      ReadFileNames
    private ::      ReadGlobalOptions
#ifdef _PHREEQC_    
    private ::      ConstructPhreeqCModels
    private ::          LoadPhreeqCModelInfo
#endif
    private ::      Construct_PropertyList
    private ::          Construct_Property
    private ::              Construct_PropertyState
    private ::              Construct_PropertyValues
    private ::              Construct_PropertyEvolution
    private ::                  ReadAdvectionDiffusionParam
    private ::                  ConstructPropertyDiffusivity 
    private ::              Construct_PropertyOutPut   
    private ::      ConstructPartition 
    private ::      ConstructHDF
    private ::      ConstructTimeSerie
    private ::      CoupleSoilQuality
#ifdef _PHREEQC_       
    private ::      CoupleSoilChemistry
    private ::          SetupSoilChemistry
#endif
    private ::      CoupleChainReactions
    private ::      CheckFlowDirections
    
    !Selector
    public  :: GetPMPMassBalance
    public  :: GetPMPnProperties
    public  :: GetPMPPropertiesIDByIdx    
    public  :: GetPMPConcentration
    public  :: GetPMPConcentrationOld
    public  :: GetPMPCoupled
    public  :: GetECw
    public  :: CheckPMPProperty
    public  :: SetDNConcPMP              !PMP gets conc from Drainage network
    public  :: SetInfColConcPMP          !PMP gets infcol conc from basin
    public  :: SetVegetationPMProperties !PMP gets conc from vegetation
    public  :: SetWindVelocity                 
    public  :: UnGetPorousMediaProperties
    
    !Modifier
    public  :: ModifyPorousMediaProperties
    private ::      ActualizePropertiesFromFile
    private ::      InterfaceFluxes
    private ::      AdvectionDiffusionProcesses 
    private ::          ComputeVolumes
    private ::          ComputeThetaAtFaces
    private ::              ComputeThetaAtFacesByAvg
    private ::              ComputeThetaAtFacesByMin
    private ::      SoilQualityProcesses
#ifdef _PHREEQC_    
    private ::      SoilChemistryProcesses
#endif        
    private ::      ChainReactionsProcesses
    private ::      Partition_Processes
    private ::      ComputeECw
    private ::      OutPut_TimeSeries
    private ::      Output_HDF
    
    !Destructor
    public  :: KillPorousMediaProperties                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjPorousMediaProperties 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetPorousMediaProperties3D_I
    private :: UnGetPorousMediaProperties3D_R8
    private :: UnGetPorousMediaProperties3D_R4
    interface  UnGetPorousMediaProperties
        module procedure UnGetPorousMediaProperties3D_I
        module procedure UnGetPorousMediaProperties3D_R8
        module procedure UnGetPorousMediaProperties3D_R4
    end interface  UnGetPorousMediaProperties


    !Parameters-----------------------------------------------------------------

    real,    parameter :: WaterReferenceDensity = 998. ![kg/m3]
    
    !For integration with options from GWflow to river (if flow is by layer or not)
    integer, parameter :: Layer_ = 2
    
    integer, parameter :: DirectionX = 1
    integer, parameter :: DirectionY = 2

    character(LEN = StringLength), parameter :: prop_block_begin     = '<beginproperty>'
    character(LEN = StringLength), parameter :: prop_block_end       = '<endproperty>'
    character(LEN = *),            parameter :: ModuleName           = 'PorousMediaProperties'    
    
    integer, parameter                :: AdvDif_ModulePMP_  = 1
    integer, parameter                :: AdvDif_ModuleAD_   = 2
    
    !Boundary condition
    integer, parameter                      :: ImposedValue_  = 1
    integer, parameter                      :: NullGradient_  = 2
    
    !Decay Equations
    integer, parameter                      :: Peyrard_       = 1
 
    !Types---------------------------------------------------------------------
    
    private :: T_PorousMediaProperties

    type T_RelatedID
        integer                       :: IDNumber = -1
        character(LEN = StringLength) :: name     = ''   
    end type T_RelatedID

    type T_ID
        integer                       :: IDNumber    = null_int
        character(LEN = StringLength) :: name        = null_str
        character(LEN = StringLength) :: description = null_str
        character(LEN = StringLength) :: units       = null_str
    end type T_ID

    type T_Property_3D
        type(T_PropertyID)               :: ID
        real, pointer, dimension (:,:,:) :: Field    => null()
        real                             :: Scalar   = null_real
    end type T_Property_3D


    type T_ExtVar
        !Map
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D   => null()
        integer, pointer, dimension(:,:,:)          :: OpenPoints3D    => null()
        integer, pointer, dimension(:,:,:)          :: LandPoints3D    => null()
        integer, dimension(:,:), pointer            :: BasinPoints     => null()
        integer, dimension(:,:), pointer            :: RiverPoints     => null()
        integer, pointer, dimension(:,:,:)          :: ComputeFacesU3D => null()
        integer, pointer, dimension(:,:,:)          :: ComputeFacesV3D => null()
        integer, pointer, dimension(:,:,:)          :: ComputeFacesW3D => null() !from basin
        integer, pointer, dimension(:,:)            :: KFloor          => null()
        real                                        :: PorousMediapropDT = null_real
        type(T_Time)                                :: Now              
        type(T_Time)                                :: BeginTime
        type(T_Time)                                :: EndTime
   
        ! from porousMedia
        real,    dimension(:,:,:), pointer          :: UnSatW              => null()
        real,    dimension(:,:,:), pointer          :: UnSatV              => null()
        real,    dimension(:,:,:), pointer          :: UnSatU              => null()
        real,    pointer, dimension(:,:,:)          :: WaterContent        => null() 
        real,    pointer, dimension(:,:,:)          :: WaterContentOld     => null()
        real(8), pointer, dimension(:,:)            :: WaterColumn         => null()
        real(8), pointer, dimension(:,:)            :: WaterColumnOld      => null()
        real(8), pointer, dimension(:,:)            :: InfiltrationColumn  => null()
        real(8), pointer, dimension(:,:,:)          :: CellVolume          => null()
        real(8), dimension(:,:,:), pointer          :: FluxU               => null()
        real(8), dimension(:,:,:), pointer          :: FluxV               => null()
        real(8), dimension(:,:,:), pointer          :: FluxW               => null()
        real,    pointer, dimension(:,:,:)          :: ThetaS              => null()
        real,    pointer, dimension(:,:,:)          :: ThetaF              => null()
        real   , pointer, dimension(:,:  )          :: Area                => null()
        real   , pointer, dimension(:,:,:)          :: AreaU               => null()
        real   , pointer, dimension(:,:,:)          :: AreaV               => null()
        real                                        :: DT                  = null_real
        real   , pointer, dimension(:,:,:)          :: DWZ                 => null()
        real   , pointer, dimension(:,:,:)          :: DZZ                 => null()
        real   , pointer, dimension(:,:  )          :: DVY                 => null()
        real   , pointer, dimension(:,:  )          :: DZX                 => null()
        real   , pointer, dimension(:,:  )          :: DUX                 => null()
        real   , pointer, dimension(:,:  )          :: DZY                 => null()
        real   , pointer, dimension(:,:  )          :: Topography          => null()
        real ,   pointer, dimension(:,:,:)          :: SZZ                 => null()
        real,    pointer, dimension(:,:)            :: FlowToChannels      => null()
        real,    pointer, dimension(:,:,:)          :: FlowToChannelsLayer => null()
        integer, pointer, dimension(:,:)            :: GWFlowBottomLayer   => null()
        integer, pointer, dimension(:,:)            :: GWFlowTopLayer      => null()
        integer, pointer, dimension(:,:)            :: GWLayerOld          => null()
        logical                                     :: BoundaryImposed     = .false.
        logical                                     :: BoundaryImposedWalls = .false.
        logical                                     :: BoundaryImposedBottom = .false.
        integer, pointer, dimension(:,:)            :: BoundaryCells       => null()
        real   , pointer, dimension(:,:,:)          :: BoundaryFluxWalls   => null()
        real   , pointer, dimension(:,:  )          :: BoundaryFluxBottom  => null()
        
        !from vegetation
        logical                                     :: ComputeVegInterfaceFluxes = .false.
        logical, dimension(:,:  ), pointer          :: SoilFluxesActive          => null()
        real,    dimension(:,:  ), pointer          :: GrazingBiomass            => null()
        real,    dimension(:,:  ), pointer          :: GrazingNitrogen           => null()
        real,    dimension(:,:  ), pointer          :: GrazingPhosphorus         => null()
        real,    dimension(:,:  ), pointer          :: HarvestKillAerialBiomass  => null()
        real,    dimension(:,:  ), pointer          :: HarvestKillNitrogen       => null()
        real,    dimension(:,:  ), pointer          :: HarvestKillPhosphorus     => null()
        real,    dimension(:,:  ), pointer          :: HarvestKillRootBiomass    => null()
        real,    dimension(:,:  ), pointer          :: DormancyBiomass           => null()
        real,    dimension(:,:  ), pointer          :: DormancyNitrogen          => null()
        real,    dimension(:,:  ), pointer          :: DormancyPhosphorus        => null()
        real,    dimension(:,:  ), pointer          :: FertilNitrateSurface      => null()
        real,    dimension(:,:  ), pointer          :: FertilNitrateSubSurface   => null()
        real,    dimension(:,:  ), pointer          :: FertilAmmoniaSurface      => null()
        real,    dimension(:,:  ), pointer          :: FertilAmmoniaSubSurface   => null()
        real,    dimension(:,:  ), pointer          :: FertilOrganicNSurface     => null()
        real,    dimension(:,:  ), pointer          :: FertilOrganicNSubSurface  => null()
        real,    dimension(:,:  ), pointer          :: FertilOrganicPSurface     => null()
        real,    dimension(:,:  ), pointer          :: FertilOrganicPSubSurface  => null()
        real,    dimension(:,:  ), pointer          :: FertilMineralPSurface     => null()
        real,    dimension(:,:  ), pointer          :: FertilMineralPSubSurface  => null()
        real,    dimension(:,:  ), pointer          :: NitrogenFraction          => null()
        real,    dimension(:,:  ), pointer          :: PhosphorusFraction        => null()
        real,    dimension(:,:,:), pointer          :: NitrogenUptake            => null()
        real,    dimension(:,:,:), pointer          :: PhosphorusUptake          => null()
        real,    dimension(:,:  ), pointer          :: RootDepth                 => null()
        integer, dimension(:,:  ), pointer          :: TranspirationBottomLayer  => null()
        logical                                     :: Grazing                   = .false.
        logical                                     :: HarvestKill               = .false.
        logical                                     :: Dormancy                  = .false.
        logical                                     :: Fertilization             = .false.
        logical                                     :: ModelWater                = .false.
        logical                                     :: ModelNitrogen             = .false.
        logical                                     :: ModelPhosphorus           = .false.
        logical                                     :: GrowthModel               = .false.
        logical                                     :: CoupledVegetation         = .false.
        real                                        :: VegetationDT              = null_real 
        integer                                     :: NutrientUptakeMethod      = null_int
        
        logical                                     :: CoupledDN                 = .false.
        !from basin
        real,    dimension(:,:  ), pointer          :: WindVelocity2D            => null()!m/s
        real,    dimension(:,:,:), pointer          :: WindVelocity3D            => null()!km/day
     end type T_ExtVar

    type T_OutPut
        type (T_Time), pointer, dimension(:)    :: OutTime              => null()
        type (T_Time), dimension(:), pointer    :: RestartOutTime       => null()
        integer                                 :: NextOutPut           = null_int
        integer                                 :: Number               = null_int
        logical                                 :: Yes                  = .false.
        logical                                 :: TimeSerie_ON         = .false.
        logical                                 :: HDF_ON               = .false.
        logical                                 :: Boxes_ON             = .false.
        logical                                 :: RateFluxes           = .false.
        logical                                 :: Profile_ON           = .false.
        logical                                 :: WriteRestartFile     = .false.       
        logical                                 :: RestartOverwrite     = .false.
        integer                                 :: NextRestartOutput    = 1  
        logical                                 :: AverageConc_ON       = .false.
        logical                                 :: AverageDecay_ON      = .false. 
        logical                                 :: IntegratedDecay_ON   = .false.     
    end type T_OutPut

    type T_AdvectionDiffusion   
        !--For AdvectionDiffusion module use
        integer                                :: BoundaryCondition     = null_int
        real                                   :: SchmidtNumberH        = null_real
        real                                   :: SchmidtCoefV          = null_real
        real                                   :: SchmidtBackgroundV    = null_real
        real                                   :: DiffusionH_imp_exp    = null_real
        real                                   :: ImplicitH_direction   = null_real
        logical                                :: Nulldif               = .false.
        logical                                :: NumericStability      = .false.
        real                                   :: VolumeRelMax          = null_real
        integer                                :: AdvMethodH            = null_int
        integer                                :: TVDLimitationH        = null_int
        integer                                :: AdvMethodV            = null_int
        integer                                :: TVDLimitationV        = null_int
        logical                                :: Upwind2H              = .false.
        logical                                :: Upwind2V              = .false.
        logical                                :: Adv_Dif_Explicit      = .false.                
        !--For both models use
        real                                   :: Molecular_Diff_Coef   = null_real
                           
    end type T_AdvectionDiffusion

    type       T_Partition
        character (LEN = StringLength)          :: Couple               = null_str
        integer                                 :: Couple_ID            = FillValueInt
        real                                    :: Fraction             = FillValueReal
        real                                    :: Rate                 = FillValueReal
        real                                    :: EmpiricCoef          = FillValueReal
        real                                    :: SedimentRefConc      = FillValueReal
        logical                                 :: UseSedimentRefConc   = .false.
        logical                                 :: SalinityEffect       = .false.
        
        !Isothermic Adsortion
!        logical                                 :: NonLinear
!        character(LEN = StringLength)           :: NonLinear_ks_Units
!        type(T_Property_2D)                     :: Nu            
!        type(T_Property_2D)                     :: Be          
!        type(T_Property_2D)                     :: ks
!        type(T_Property_2D)                     :: PartitionRate
!        type(T_Property_2D)                     :: Fraction2D 
!        character (LEN = StringLength)          :: Partition_Couple        
    end type T_Partition

    type T_Boundary
        integer                                 :: BoundaryCondition   = null_int
        real                                    :: DefaultBoundary     = null_real
    end type T_Boundary

    type T_Evolution
        logical                                 :: Variable             = .false.
        real                                    :: DTInterval           = null_real
        logical                                 :: DTIntervalAssociated = .false.        
        type(T_Time)                            :: LastCompute
        type(T_Time)                            :: NextCompute
        logical                                 :: SoilQuality          = .false.
        logical                                 :: SoilChemistry        = .false.
        logical                                 :: Partitioning         = .false.
        logical                                 :: CationExchangeProcess = .false.
        logical                                 :: ChemEquilibriumProcess = .false.
        logical                                 :: TransportedInEVTP    = .false.    
        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: SoilWaterFluxes      = .false.
        logical                                 :: Macropores           = .false.
        logical                                 :: MinConcentration     = .false.
        logical                                 :: WarnOnNegativeValues = .false.                
        logical                                 :: UseMaxForUptakeConc  = .false.      
        real                                    :: MaxForUptakeConc     = 0.0   
        logical                                 :: Decay                = .false.
        logical                                 :: DecayEquationCoupled = .false.
        integer                                 :: DecayEquation        = null_int
        real                                    :: DecayHalfSaturation  = null_real
        real                                    :: DecayRate            = null_real !day-1
        logical                                 :: DecayMass            = .false.
        logical                                 :: Discharges           = .false.        
        type (T_AdvectionDiffusion)             :: AdvDiff
        type (T_Partition)                      :: Partition
        type (T_Boundary)                       :: Boundary        
    end type T_Evolution
    
    type T_ThetaAtFaces
        real, dimension(:,:,:), pointer         :: ThetaW              => null()
        real, dimension(:,:,:), pointer         :: ThetaV              => null()
        real, dimension(:,:,:), pointer         :: ThetaU              => null()
    end type T_ThetaAtFaces
    
    type T_MassBalance
        real(8)                                 :: TotalStoredMass     = null_real
        real(8)                                 :: TranspiredMass      = null_real
        real(8)                                 :: DNExchangeMass      = null_real
        real(8)                                 :: RPExchangeMass      = null_real
        real(8)                                 :: TotalDischargeMass  = null_real        
    end type T_MassBalance
    
    type T_Files
        character(PathLength)                   :: InitialFile         = null_str
        character(PathLength)                   :: DataFile            = null_str
        character(PathLength)                   :: FinalFile           = null_str
        character(PathLength)                   :: TransientHDF        = null_str
        character(PathLength)                   :: DataSedimentQualityFile = null_str
        character(PathLength)                   :: ChainReactionsDataFile = null_str
        character(len=StringLength)             :: BoxesFile           = null_str
        integer                                 :: AsciiUnit           = null_int
    end type T_Files    

    type T_Property
        type (T_PropertyID)                     :: ID
        real, dimension(:,:,:), pointer         :: Concentration            => null()
#ifdef _USE_PAGELOCKED
        type(C_PTR)                             :: ConcentrationPtr         => null()
#endif _USE_PAGELOCKED
        real, dimension(:,:,:), pointer         :: ConcentrationOld         => null()
        real, dimension(:,:), pointer           :: ConcentrationOnInfColumn => null()
        real, dimension(:,:), pointer           :: ConcentrationDN          => null()
        real, dimension(:,:,:), pointer         :: ConcInInterfaceDN        => null()
        real, dimension(:,:,:), pointer         :: ConcInBoundary           => null()
        real, dimension(:,:  ), pointer         :: AverageAquiferConc       => null()
        real, dimension(:,:  ), pointer         :: AverageVadozeConc        => null()
        real, dimension(:,:  ), pointer         :: AverageAquiferDecay      => null()
        real, dimension(:,:  ), pointer         :: AverageVadozeDecay       => null()  
        real, dimension(:,:  ), pointer         :: IntegratedDecay          => null() 
        character(PathLength)                   :: IntegratedDecayFile      = null_str     
        real, dimension(:,:),   pointer         :: PesticideFlux            => null()
        logical, pointer, dimension(:,:,:)      :: UptakeActive             => null()
        real, pointer, dimension(:,:,:)         :: Mass_Created             => null()
        real(8)                                 :: TotalStoredMass          = null_real
        real, pointer, dimension(:,:,:)         :: ViscosityU               => null()
        real, pointer, dimension(:,:,:)         :: ViscosityV               => null()
        type (T_Property), pointer              :: Next, Prev               => null()
        logical                                 :: Particulate              = .false.
        logical                                 :: Pesticide                = .false.
        type (T_Evolution)                      :: Evolution       
        type (T_MassBalance)                    :: MB
        real, pointer, dimension(:,:,:)         :: Diffusivity              => null()
        real, pointer, dimension(:,:,:)         :: Diff_Turbulence_H        => null()
        real, pointer, dimension(:,:,:)         :: Diff_Turbulence_V        => null()
        real, pointer, dimension(:,:,:)         :: Viscosity                => null()
        logical                                 :: Old                      = .false.
        real                                    :: MinValue                 = FillValueReal
        logical                                 :: TimeSerie                = .false.
        logical                                 :: BoxTimeSerie             = .false.
        logical                                 :: BoxTimeSerie2D           = .false.
        logical                                 :: OutputHDF                = .false. 
        logical                                 :: OutputAverageConc        = .false. 
        logical                                 :: OutputAverageDecay       = .false.
        logical                                 :: OutputIntegratedDecay    = .false. 
        logical                                 :: UseToCalcECw             = .false.
        real                                    :: ECwFactor                = 0.0
        real, pointer, dimension(:,:,:)         :: PropertyDecay            => null()
    end type T_Property

    type       T_SedimentRate
        type (T_ID)                             :: ID
        type (T_ID)                             :: FirstProp
        type (T_ID)                             :: SecondProp
        real, pointer, dimension(:,:,:)         :: Field                    => null()
        type(T_SedimentRate), pointer           :: Next                     => null()
        type(T_SedimentRate), pointer           :: Prev                     => null()
    end type T_SedimentRate


    type T_Coupled
        logical                                 :: SoilQuality          = .false. !Sediment source/sink model (Sediment Quality)
        real                                    :: SoilQuality_DT       = null_real
        type (T_Time)                           :: SoilQuality_NextCompute

#ifdef _PHREEQC_        
        logical                                 :: SoilChemistry        = .false.  !Chemical reactions model (PhreeqC)
        real                                    :: SoilChemistry_DT     = null_real
        type (T_Time)                           :: SoilChemistry_NextCompute
#endif        

        logical                                 :: ChainReactions       = .false. !Simple generic reactions ZERO & FIRST order rate

        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: Partition            = .false.
        logical                                 :: MinConcentration     = .false.
        logical                                 :: WarnOnNegativeValues = .false.
        logical                                 :: Decay                = .false.
        logical                                 :: SedQualityOxygenForcing = .false.
        logical                                 :: Discharges           = .false.
        
    end type T_Coupled
    
    !Implicit coef for thomas matrix
    type       T_DEF
        real   , pointer, dimension(: , : , :)  :: D   => null()
        real(8), pointer, dimension(: , : , :)  :: E   => null()
        real   , pointer, dimension(: , : , :)  :: F   => null()
    end type T_DEF
    
    !Explicit coefs
    type       T_A_B_C_Explicit
!        real(8), pointer, dimension(: , : , :)  :: CoefA_W
!        real(8), pointer, dimension(: , : , :)  :: CoefB_W
!        real(8), pointer, dimension(: , : , :)  :: CoefC_W
!        real(8), pointer, dimension(: , : , :)  :: CoefA_U
!        real(8), pointer, dimension(: , : , :)  :: CoefB_U
!        real(8), pointer, dimension(: , : , :)  :: CoefC_U
!        real(8), pointer, dimension(: , : , :)  :: CoefA_V
!        real(8), pointer, dimension(: , : , :)  :: CoefB_V
!        real(8), pointer, dimension(: , : , :)  :: CoefC_V
        real, pointer, dimension(: , : , :)  :: CoefInterfRunoff  => null() !transport from and to Runoff
        real, pointer, dimension(: , : , :)  :: CoefInterfDN      => null() !transport from and to River
        real, pointer, dimension(: , : , :)  :: CoefInterfTransp  => null() !transport to plants
        real, pointer, dimension(: , : , :)  :: CoefInterfEvap    => null() !transport to evaporation
        real, pointer, dimension(: , : , :)  :: CoefInterfBoundaryWalls => null() !transport trough boundary walls
        real, pointer, dimension(: , : , :)  :: CoefInterfBoundaryBottom => null() !transport trough boundary bottom
    end type T_A_B_C_Explicit

    type       T_FluxCoef
        real   , pointer, dimension(: , : , :)  :: C_flux   => null()  !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: D_flux   => null() !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: E_flux   => null() !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: F_flux   => null() !Coeficient to calculate AdvFlux and DifFlux
    end type T_FluxCoef

    type  T_Fluxes
        real, pointer, dimension(:,:,:)         :: AdvFluxX   => null()
        real, pointer, dimension(:,:,:)         :: AdvFluxY   => null()
        real, pointer, dimension(:,:,:)         :: AdvFluxZ   => null()
        
        real, pointer, dimension(:,:,:)         :: DifFluxX   => null()
        real, pointer, dimension(:,:,:)         :: DifFluxY   => null()
        real, pointer, dimension(:,:,:)         :: DifFluxZ   => null()
        
        real, pointer, dimension(:,:,:)         :: MassFluxesX   => null()
        real, pointer, dimension(:,:,:)         :: MassFluxesY   => null()
        real, pointer, dimension(:,:,:)         :: MassFluxesZ   => null()
    end type T_Fluxes

#ifdef _PHREEQC_  
    type T_PhreeqCModel
        integer               :: ID           = null_int
        integer               :: KUB          = null_int
        integer               :: KLB          = null_int
        character(PathLength) :: Database     = ''
        character(PathLength) :: DatabaseAux  = ''
        integer               :: ObjID        =  0
        integer               :: ObjInterface =  0 
    end type T_PhreeqCModel  
    
    type T_PhreeqC
        real, dimension(:,:,:), pointer              :: CellSoilMass              => null()
        real, dimension(:,:,:), pointer              :: CellWaterMass             => null()
        type(T_Property), pointer                    :: SoilDryDensity            => null()    
        type(T_PhreeqCModel), dimension(:), pointer  :: Models                    => null()
        integer                                      :: NumberOfModels            =  0
        integer, dimension(:,:,:), pointer           :: Filter                    => null()            
    end type T_PhreeqC
#endif   
    
    type T_PorousMediaProperties
        integer                                      :: ObjTime               = 0
        integer                                      :: ObjHorizontalGrid     = 0
        integer                                      :: ObjHorizontalMap      = 0
        integer                                      :: ObjAdvectionDiffusion = 0
        integer                                      :: ObjBasinGeometry      = 0
        integer                                      :: ObjPorousMedia        = 0
        integer                                      :: ObjGeometry           = 0
        integer                                      :: ObjMap                = 0
        integer                                      :: ObjGridData           = 0
        integer                                      :: ObjEnterData          = 0
        integer                                      :: ObjtimeSerie          = 0
        integer                                      :: ObjSedimentQuality    = 0
        integer                                      :: ObjHDF5               = 0
        integer                                      :: ObjBottomTopography   = 0
        integer                                      :: ObjProfile            = 0
        integer                                      :: ObjInterface          = 0
        integer                                      :: ObjChainReactions     = 0 !ChainReactionsID
        integer                                      :: ObjBoxDif             = 0
#ifdef _ENABLE_CUDA
        integer                                      :: ObjCuda               = 0
#endif _ENABLE_CUDA        
        integer                                      :: ObjDischarges         = 0

        real,    pointer, dimension(:,:,:)           :: CellWaterVolume       => null() !Used SoilChemistry & ChainReactions        
        integer, pointer, dimension(:)               :: PropertiesList        => null() !List with ID of all properties used 
 
        type (T_ExtVar)                              :: ExtVar
        type (T_Files)                               :: Files
        type (T_OutPut)                              :: OutPut
        type (T_Property), pointer                   :: FirstProperty    => null() 
        type (T_Property), pointer                   :: LastProperty     => null()   
        type(T_SedimentRate), pointer                :: FirstSedimentRate => null()
        type(T_SedimentRate), pointer                :: LastSedimentRate  => null()
        type (T_PorousMediaProperties), pointer      :: Next             => null()
        type (T_Coupled)                             :: Coupled
        type (T_Time)                                :: LastOutputHDF5
        type(T_D_E_F)                                :: COEF3
        type(T_A_B_C_Explicit)                       :: COEFExpl 
        type(T_FluxCoef)                             :: COEF3_VertAdv            !Vertical advection coeficients
        type(T_FluxCoef)                             :: COEF3_HorAdvXX           !Horizont advection coeficients
        type(T_FluxCoef)                             :: COEF3_HorAdvYY           !Horizont advection coeficients
        type(T_Fluxes)                               :: Fluxes
        real, pointer, dimension(: , : , :)          :: TICOEF3          => null()      
#ifdef _USE_PAGELOCKED
        type(C_PTR)                                  :: TICOEF3Ptr
#endif _USE_PAGELOCKED      
!        real(8), pointer, dimension(:)               :: VECG                     !Auxiliar thomas arrays 
!        real(8), pointer, dimension(:)               :: VECW                     !Auxiliar thomas arrays     

        logical                                      :: PorousMediaProperties  = .false.
        real,    pointer, dimension(:,:,:)           :: Volume                 => null()
        integer                                      :: PropertiesNumber       = 0
        integer                                      :: SedimentRatesNumber    = 0
        integer                                      :: NumberPropForBoxes     = 0
        real   , pointer, dimension(:,:,:)           :: DissolvedToParticulate3D => null()
        real                                         :: ResidualTime           = null_real
        
        integer                                      :: InstanceID             = null_int
        type (T_Size3D)                              :: Size, WorkSize
        type (T_Size2D)                              :: Size2D

        type(T_Property_3D)                          :: Disper_Trans
        type(T_Property_3D)                          :: Disper_Longi
               
        logical                                      :: AdvDiff_Explicit        = .false.    ! 1 All explicit; 0 Vertical Implicit
        logical                                      :: AdvDiff_AdvectionH_ImpExp = .false.  ! 1 Horiz Adv Expl; 0 Horiz Adv Impl
        logical                                      :: AdvDiff_CheckCoefs      = .false.
        logical                                      :: Vertical1D              = .false.
        logical                                      :: XZFlow                  = .false.
        
        logical                                      :: NewFormulation          = .false.   !New formulation for advection 
        integer                                      :: AdvDiff_AdvMethodV      = null_int  !
        integer                                      :: AdvDiff_AdvMethodH      = null_int  !methods are general for all properties
       

        type(T_ThetaAtFaces)                         :: ThetaAtFaces
        integer                                      :: ThetaAtFacesMethod      = null_int
         
#ifdef _PHREEQC_
        type(T_PhreeqC)                              :: PhreeqC
#endif        

        logical                                      :: CalculateECw =  .false.
        real, dimension(:,:,:), pointer              :: ECw          => null()
        real, dimension(:,:,:), pointer              :: CellMass     => null()

        !--For PorousMediaProperties Advection-Diffusion Method
        real,    pointer, dimension(:,:,:)           :: DifusionNumber  => null()
        real,    pointer, dimension(:,:,:)           :: ReynoldsMNumber => null()  

        
        logical                                      :: CheckGlobalMass = .false.     
        
        real(8), pointer, dimension(:,:,:)           :: WaterVolume     => null()
        real(8), pointer, dimension(:,:,:)           :: FluxWCorr       => null()
        real,    pointer, dimension(:,:,:)           :: WaterContentBT  => null()  !Water Content Before Transport
        
        logical                                      :: DTIntervalAssociated     = .false.

        integer                                      :: nPropWithDischarge   = 0
              
        !griflet, openmp
        type(T_THOMAS), pointer                      :: THOMAS          => null()
        integer                                      :: MaxThreads      = null_int

    end type  T_PorousMediaProperties

    !Global Module Variables
    type (T_PorousMediaProperties), pointer          :: FirstObjPorousMediaProperties => null()
    type (T_PorousMediaProperties), pointer          :: Me                            => null()

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructPorousMediaProperties(ObjPorousMediaPropertiesID,                 &
                                              ComputeTimeID,                              &
                                              GridDataID,                                 &
                                              HorizontalGridID,                           &
                                              HorizontalMapID,                            &
                                              BasinGeometryID,                            &
                                              PorousMediaID,                              &
                                              GeometryID,                                 &
                                              MapID,                                      &
                                              DischargesID,                               &
                                              CoupledDN,                                  &
                                              CheckGlobalMass,                            &
#ifdef _ENABLE_CUDA
                                              CudaID,                                   &
#endif _ENABLE_CUDA    
                                              STAT)
     
        !Arguments---------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID 
        integer                                         :: ComputeTimeID
        integer                                         :: GridDataID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: BasinGeometryID
        integer                                         :: PorousMediaID
        integer                                         :: GeometryID
        integer                                         :: MapID
        integer                                         :: DischargesID
        logical, optional                               :: CoupledDN
        logical                                         :: CheckGlobalMass
#ifdef _ENABLE_CUDA
        integer                                         :: CudaID
#endif _ENABLE_CUDA
        integer, optional, intent(OUT)                  :: STAT 
        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_,STAT_CALL
        !------------------------------------------------------------------------
                                    

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mPorousMediaProperties_)) then
            nullify (FirstObjPorousMediaProperties)
            call RegisterModule (mPorousMediaProperties_) 
        endif

        call Ready(ObjPorousMediaPropertiesID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then


            call AllocateInstance
            
            !Associate External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           ComputeTimeID   )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID      )            
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjPorousMedia    = AssociateInstance (mPOROUSMEDIA_,    PorousMediaID   )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMap_,            MapID           )
#ifdef _ENABLE_CUDA
            Me%ObjCuda           = AssociateInstance (mCUDA_,           CudaID          )
#endif
        
            if (present(CoupledDN)) then
                Me%ExtVar%CoupledDN = CoupledDN
            endif
            Me%CheckGlobalMass       = CheckGlobalMass
            
            call CheckBoundary
            
            call ReadFileNames


            !Constructs the DataFile
            call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR01'                
           
            call ReadGlobalOptions

            call Construct_PropertyList
            
            !Rates output for sediment quality model
            call Construct_SedimentRateList

#ifdef _PHREEQC_ 
            if (Me%Coupled%SoilChemistry) then           
                call ConstructPhreeqCModels
            endif
#endif    

            if (Me%Coupled%Partition) then
                call ConstructPartition
            end if     

            if (Me%Coupled%Discharges) then
                call ConstructDischarges(DischargesID)
            endif

            
            call AllocateVariables
       
            call ConstructHDF    
    
            call ConstructTimeSerie
            
            if (Me%Output%RateFluxes .or. Me%Output%Boxes_ON) then
                call StartOutputBoxFluxes
            endif
            
            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR010'


            !Couple nutrient, carbon and oxygen sources and sinks model
            if (Me%Coupled%SoilQuality) then
                call CoupleSoilQuality
            endif
            
#ifdef _PHREEQC_            
            !Couple soil chemical model
            if (Me%Coupled%SoilChemistry) then
                call CoupleSoilChemistry
            endif
#endif

            if (Me%Coupled%AdvectionDiffusion) then
                !See if simulation is Vertical1D or XZ direction - to spare resources on computation
                call CheckFlowDirections
            endif
              
            if (Me%Coupled%ChainReactions) then
                call CoupleChainReactions
            endif
                       
            !Message to user
            call ConstructLog
            

            if (Me%CheckGlobalMass) then
                call CalculateTotalStoredMass
            endif


            !First Output
            if (Me%Output%HDF_ON) then
                call ReadLockExternalVar
                call OutPut_HDF
                call ReadUnLockExternalVar
            endif



            !Returns ID
            ObjPorousMediaPropertiesID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR060' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructPorousMediaProperties

    !--------------------------------------------------------------------------

    subroutine ConstructDischarges(DischargesID)

        !Arguments--------------------------------------------------------------
        integer                                     :: DischargesID 
        !Local------------------------------------------------------------------
        !integer                                     :: STAT_CALL
        !integer                                     :: nDischarges, iDis
        !character(len=StringLength)                 :: DischargeName
        !real                                        :: CoordinateX, CoordinateY
        !logical                                     :: CoordinatesON, IgnoreOK
        !integer                                     :: Jd
         

        !ObjDischarges comes from ModueRunoff
!        call Construct_Discharges(Me%ObjDischarges,                              &
!                                  Me%ObjTime,                                    &
!                                  STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR01'  

        !ATTENTION - DIFFERENT DISCHARGES CAN NOT OCCUR AT THE SAME CELL
        !Test this here in PMP because is only problem if properties are used
        !or create duplications of integrated discharge in porous media

        if (DischargesID == 0)  then                                                
            write(*,*)'You need to define DISCHARGES : 1 in the POROUS MEDIA input file' 
            stop      'ModulePorousMediaProperties - ConstructDischarges - ERR01'
        else            
            Me%ObjDischarges     = AssociateInstance (mDISCHARGES_,     DischargesID    )
        endif
           
    
    endsubroutine ConstructDischarges            
    
    !--------------------------------------------------------------------------
     
    subroutine CheckBoundary

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                        :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        !check if boundary is imposed
        call GetBoundaryImposed (Me%ObjPorousMedia,               &
                                 Me%ExtVar%BoundaryImposed,       &
                                 Me%ExtVar%BoundaryImposedWalls,  &
                                 Me%ExtVar%BoundaryImposedBottom, &
                                 STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckBoundary - ModulePorousMedia - ERR01'
        
        if (Me%ExtVar%BoundaryImposedWalls) then
            call GetBoundaryCells (Me%ObjPorousMedia, Me%ExtVar%BoundaryCells, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckBoundary - ModulePorousMedia - ERR010'
            
        endif
        
    end subroutine CheckBoundary

    !--------------------------------------------------------------------------
    
    subroutine CheckFlowDirections

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                                     :: YFlow, XFlow, ZFlow  
        integer                                     :: ZComputeCells, XComputeCells
        integer                                     :: YComputeCells
        !Begin-----------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB


        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckDirections - ModulePorousMediaProperties - ERR0140'

       
        ZFlow = .false.
        XFlow = .false.
        YFlow = .false.
        
        ZComputeCells = 0
        YComputeCells = 0
        XComputeCells = 0
        
        !Z evaluation
doi:    do I = ILB, IUB
        do J = JLB, JUB
            
            ZComputeCells = 0
            
            do K = KLB, KUB
                
                if (Me%ExtVar%WaterPoints3D(i, j, k) == WaterPoint) then     
                    
                    ZComputeCells = ZComputeCells + 1
                    
                    if (ZComputeCells .gt. 1) then
                        ZFlow = .true.
                        exit doi
                    endif
                    
                endif
            
            enddo
        enddo
        enddo   doi


        !X evaluation
doi2:   do I = ILB, IUB
            
            XComputeCells = 0
            
            do J = JLB, JUB
                
                if (Me%ExtVar%WaterPoints3D(i, j, k) == WaterPoint) then     
                    
                    XComputeCells = XComputeCells + 1
                    
                    if (XComputeCells .gt. 1) then
                        XFlow = .true.
                        exit doi2
                    endif
                    
                endif
                
            enddo
        enddo   doi2


        !Y evaluation
doi3:   do J = JLB, JUB
            
            YComputeCells = 0
            
            do I = ILB, IUB
                
                if (Me%ExtVar%WaterPoints3D(i, j, k) == WaterPoint) then     
                    
                    YComputeCells = YComputeCells + 1
                    
                    if (YComputeCells .gt. 1) then
                        YFlow = .true.
                        exit doi3
                    endif
                    
                endif
                
            enddo
        enddo   doi3

        !1D flow in Z direction
        Me%Vertical1D = .false.
        !2D flow in XZ direction
        Me%XZFlow     = .false.
        
        !If in the end Xflow exists but never found Yflow then is a 2D XZflow
        if (ZFlow .and. XFlow .and. (.not. YFlow)) then
            Me%XZFlow = .true.
        endif
        !If in the end only Zflow exists is a 1D vertical flow
        if (ZFlow .and. (.not. XFlow) .and. (.not. YFlow)) then
            Me%Vertical1D = .true.
        endif


        call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR160'


        
    end subroutine CheckFlowDirections

    !--------------------------------------------------------------------------
    
    subroutine ConstructLog


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty


        write(*, *)"---------------- POROUS MEDIA PROPERTIES -----------------"
        write(*, *)
        write(*, *)"Num of Properties : ", Me%PropertiesNumber
        write(*, *)

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))

            write(*, *)"Property            : ", trim(CurrentProperty%ID%Name)
            write(*, *)"---Adv. Diff.       : ", CurrentProperty%Evolution%AdvectionDiffusion
            write(*, *)"---Sediment Quality : ", CurrentProperty%Evolution%SoilQuality
            write(*, *)"---PREEQC           : ", CurrentProperty%Evolution%SoilChemistry
            write(*, *)"---Partitioning     : ", CurrentProperty%Evolution%Partitioning
            write(*, *)"---Particulate      : ", CurrentProperty%Particulate
            write(*, *)"---Discharges       : ", CurrentProperty%Evolution%Discharges
            write(*, *)

            CurrentProperty=>CurrentProperty%Next
        enddo

    end subroutine ConstructLog
    
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_PorousMediaProperties), pointer                         :: NewObjPorousMediaProperties
        type (T_PorousMediaProperties), pointer                         :: PreviousObjPorousMediaProp


        !Allocates new instance
        allocate (NewObjPorousMediaProperties)
        nullify  (NewObjPorousMediaProperties%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjPorousMediaProperties)) then
            FirstObjPorousMediaProperties         => NewObjPorousMediaProperties
            Me                    => NewObjPorousMediaProperties
        else
            PreviousObjPorousMediaProp      => FirstObjPorousMediaProperties
            Me                    => FirstObjPorousMediaProperties%Next
            do while (associated(Me))
                PreviousObjPorousMediaProp  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjPorousMediaProperties
            PreviousObjPorousMediaProp%Next => NewObjPorousMediaProperties
        endif

        Me%InstanceID = RegisterNewInstance (mPorousMediaProperties_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    subroutine ReadFileNames

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
!        integer                                     :: iflag

        !Reads the name of the data file from nomfich
        call ReadFileName ('POROUS_PROP_DATA', Me%Files%DataFile, "PorousMedia Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR010'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('POROUS_PROP_HDF', Me%Files%TransientHDF, "PorousMedia HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR020'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('POROUS_PROP_FIN', Me%Files%FinalFile, "PorousMedia Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR030'
   
    end subroutine ReadFileNames
    
    !--------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        type(T_Property_3D), pointer            :: Scalar3D
        !Begin-----------------------------------------------------------------

        !Geometry Size
        call GetGeometrySize    (Me%ObjGeometry,             &    
                                 Size     = Me%Size,         &
                                 WorkSize = Me%WorkSize,     &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR010'

        Me%Size2D%ILB = Me%Size%ILB
        Me%Size2D%IUB = Me%Size%IUB
        Me%Size2D%JLB = Me%Size%JLB
        Me%Size2D%JUB = Me%Size%JUB

        call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR020'

        call GetComputeTimeLimits(Me%ObjTime,                      &
                                  EndTime   = Me%ExtVar%EndTime,   &
                                  BeginTime = Me%ExtVar%BeginTime, &
                                  STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR030'

        ! Sets the last output equal to zero 
        call SetDate(Me%LastOutPutHDF5, 0, 0, 0, 0, 0, 0)

        call GetData(Me%Coupled%ChainReactions,                    &   
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'CHAINREACTIONS_MODULE',       &
                     Default      = .false.,                       &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR040'

        if (Me%Coupled%ChainReactions) then
                call GetData(Me%Files%ChainReactionsDataFile,      &   
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'CHAINREACTIONS_FILE',         &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR050'
        if (iflag .EQ. 0) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR060'
        endif

        !Want explicit or implict model? By default is explicit
        call GetData(Me%AdvDiff_Explicit,                          &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'ADVDIFF_EXPLICIT',            &
                     Default      = .true.,                        &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)  &
            stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR070'
        
        if (.not. Me%AdvDiff_Explicit) then
            !If implicit model, only vertical is forced to be implicit (adv + diff) and horizontal is explicit (adv + diff)
            !Also want advection horiz to be implicit? By default is explicit
            call GetData(Me%AdvDiff_AdvectionH_ImpExp,                 &
                         Me%ObjEnterData, iflag,                       &
                         SearchType   = FromFile,                      &
                         keyword      = 'ADVDIFF_ADVECTION_H_IMP_EXP',   &
                         Default      = .true.,                        &
                         ClientModule = 'ModulePorousMediaProperties', &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)  &
                stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR080'
            
            !verify why using horizontal implicit, concentrations are not maintained
            !when it should for instance with rain concentrations equal to amount 
            !existing at the beggining in soil    
            if (.not. Me%AdvDiff_AdvectionH_ImpExp) then
                write(*,*)
                write(*,*)'For now ADVDIFF_ADVECTION_H_IMP_EXP in PMP'
                write(*,*)'can only have value 1 (explicit in horizontal)'
                stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR085'
            endif
        endif 
                
        ! 1 - MINIMUN of the Theta in the cells
        ! 2 - AVERAGE of the Theta in the cells
                
        !Want Theta on faces using min or average? Default is min
        call GetData(Me%ThetaAtFacesMethod,                        &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'THETAFACE',                   &
                     Default      = 1,                             &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)  &
            stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR090'                                
                
        !Want new formulation with advection methods not depending on property
        call GetData(Me%NewFormulation,                            &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'NEW_FORMULATION',             &
                     Default      = .false.,                       &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)  &
            stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR100'                

        if (Me%NewFormulation) then
            !horizontal advection method that will be general for all props
            call GetData(Me%AdvDiff_AdvMethodH,                        &
                         Me%ObjEnterData, iflag,                       &
                         SearchType   = FromFile,                      &
                         keyword      = 'ADVDIFF_METHOD_H',            &
                         Default      = UpwindOrder1,                  &
                         ClientModule = 'ModulePorousMediaProperties', &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)  &
                stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR110'
            
            !vertical advection method that will be general for all props
            call GetData(Me%AdvDiff_AdvMethodV,                        &
                         Me%ObjEnterData, iflag,                       &
                         SearchType   = FromFile,                      &
                         keyword      = 'ADVDIFF_METHOD_V',            &
                         Default      = UpwindOrder1,                  &
                         ClientModule = 'ModulePorousMediaProperties', &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)  &
                stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR120'
        
        endif
        
        Scalar3D => Me%Disper_Longi
        call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                              block_begin = '<begin_dispersion_long>',          &
                              block_end   = '<end_dispersion_long>')
        
        Scalar3D => Me%Disper_Trans
        call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                              block_begin = '<begin_dispersion_trans>',         &
                              block_end   = '<end_dispersion_trans>')            

    end subroutine ReadGlobalOptions        

    !--------------------------------------------------------------------------

#ifdef _PHREEQC_
    subroutine ConstructPhreeqCModels
    
        !Local-----------------------------------------------------------------        
        integer :: status
        integer :: client
        logical :: found
        integer :: index

        !----------------------------------------------------------------------
        call GetNumberOfBlocks (Me%ObjEnterData,           &
                                '<beginphreeqcmodel>',     &
                                '<endphreeqcmodel>',       &
                                FromFile_,                 &
                                Me%PhreeqC%NumberOfModels, &
                                STAT = status)
        if (status /= SUCCESS_) &
            stop 'ConstructPhreeqCModels - ModulePorousMediaProperties - ERR010'

        allocate (Me%PhreeqC%Models (1:Me%PhreeqC%NumberOfModels))
        
        do index = 1, Me%PhreeqC%NumberOfModels

            call ExtractBlockFromBuffer(Me%ObjEnterData,                      &
                                        ClientNumber = client,                &
                                        block_begin  = '<beginphreeqcmodel>', &
                                        block_end    = '<endphreeqcmodel>',   &
                                        BlockFound   = found,                 &
                                        STAT         = status)
            if (status /= SUCCESS_) &
                stop 'ConstructPhreeqCModels - ModulePorousMediaProperties - ERR020'

            if (found) then

                call LoadPhreeqCModelInfo (index)                
                
            else

                stop 'ConstructPhreeqCModels - ModulePorousMediaProperties - ERR030'

            endif

        enddo   

        call Block_Unlock(Me%ObjEnterData, client, status)
        if (status .NE. SUCCESS_) &
            stop 'ConstructPhreeqCModels - ModulePorousMediaProperties - ERR040'
        !----------------------------------------------------------------------

    endsubroutine ConstructPhreeqCModels

    !--------------------------------------------------------------------------

    subroutine LoadPhreeqCModelInfo (Index)
    
        !Arguments-------------------------------------------------------------
        integer, intent(IN) :: Index
        
        !Local----------------------------------------------------------------
        integer               :: status
        integer               :: iflag
        
        !----------------------------------------------------------------------
        call GetData(Me%PhreeqC%Models(Index)%Database,            &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'DATABASE',                    &
                     default      = '',                            &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = status)
        if (status .NE. SUCCESS_) &
            stop 'LoadPhreeqCModelInfo - ModulePorousMediaProperties - ERR010'
            
        call GetData(Me%PhreeqC%Models(Index)%DatabaseAux,         &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'DATABASE_AUX',                &
                     default      = '',                            &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = status)
        if (status .NE. SUCCESS_) &
            stop 'LoadPhreeqCModelInfo - ModulePorousMediaProperties - ERR020'
                        
        call GetData(Me%PhreeqC%Models(Index)%KUB,                 &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'KUB',                         &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = status)
        if (status .NE. SUCCESS_) &
            stop 'LoadPhreeqCModelInfo - ModulePorousMediaProperties - ERR030'
        if (iflag .NE. 1) &
            stop 'LoadPhreeqCModelInfo - ModulePorousMediaProperties - ERR040'
        
        call GetData(Me%PhreeqC%Models(Index)%KLB,                 &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'KLB',                         &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = status)
        if (status .NE. SUCCESS_) &
            stop 'LoadPhreeqCModelInfo - ModulePorousMediaProperties - ERR050' 
        if (iflag .NE. 1) &
            stop 'LoadPhreeqCModelInfo - ModulePorousMediaProperties - ERR060'                                          
        !----------------------------------------------------------------------
        
    endsubroutine LoadPhreeqCModelInfo 
#endif
    !--------------------------------------------------------------------------

    subroutine AllocateVariables        
        
        !Local-----------------------------------------------------------------        
        integer                                         :: ILB, IUB, JLB,  JUB 
        integer                                         :: KLB, KUB, IJKLB, IJKUB 

        !griflet: openmp
        integer                                 :: m
        type(T_VECGW), pointer                  :: VECGW

        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KUB = Me%Size%KUB
        KLB = Me%Size%KLB
        
               
        !Water Content---------------------------------------------------------

        allocate (Me%ExtVar%WindVelocity3D   (ILB:IUB,JLB:JUB,KLB:KUB))
        
        allocate (Me%CellWaterVolume (ILB:IUB,JLB:JUB,KLB:KUB))
        Me%CellWaterVolume = 0.
                
#ifdef _PHREEQC_
        if (Me%Coupled%SoilChemistry) then
            allocate (Me%PhreeqC%CellWaterMass (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%PhreeqC%CellSoilMass  (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%PhreeqC%Filter        (ILB:IUB,JLB:JUB,KLB:KUB))
        
            Me%PhreeqC%CellWaterMass = 0.
            Me%PhreeqC%CellSoilMass  = 0.
            Me%PhreeqC%Filter        = 0                
        endif
#endif

        if (Me%Coupled%AdvectionDiffusion) then
            allocate (Me%WaterVolume          (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%FluxWCorr                (ILB:IUB,JLB:JUB,KLB:KUB))

            allocate (Me%WaterContentBT       (ILB:IUB,JLB:JUB,KLB:KUB))
            Me%WaterContentBT = 0.
            
            allocate (Me%ThetaAtFaces%ThetaW (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%ThetaAtFaces%ThetaU (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%ThetaAtFaces%ThetaV (ILB:IUB,JLB:JUB,KLB:KUB))

            Me%WaterVolume          = 0.
            Me%FluxWCorr            = 0.
            
    !        allocate (Me%COEFExpl%CoefA_W          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefB_W          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefC_W          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefA_U          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefB_U          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefC_U          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefA_V          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefB_V          (ILB:IUB,JLB:JUB,KLB:KUB))
    !        allocate (Me%COEFExpl%CoefC_V          (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%COEFExpl%CoefInterfRunoff (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%COEFExpl%CoefInterfDN     (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%COEFExpl%CoefInterfTransp (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%COEFExpl%CoefInterfEvap   (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%COEFExpl%CoefInterfBoundaryWalls (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%COEFExpl%CoefInterfBoundaryBottom (ILB:IUB,JLB:JUB,KLB:KUB))
            
#ifdef _USE_PAGELOCKED
            ! Allocate pagelocked memory to optimize CUDA transfers
            call Alloc3DPageLocked(Me%ObjCuda, Me%TICOEF3Ptr, Me%TICOEF3, IUB + 1, JUB + 1, KUB + 1)        
#else
            allocate(Me%TICOEF3                 (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB))
#endif _USE_PAGELOCKED            

            allocate(Me%COEF3_VertAdv%C_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%COEF3_VertAdv%D_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%COEF3_VertAdv%E_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%COEF3_VertAdv%F_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))

            Me%COEF3_VertAdv%C_Flux     = Null_real
            Me%COEF3_VertAdv%D_Flux     = Null_real
            Me%COEF3_VertAdv%E_Flux     = Null_real
            Me%COEF3_VertAdv%F_Flux     = Null_real
            
            if (.not. Me%Vertical1D) then
                allocate(Me%COEF3_HorAdvXX%C_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3_HorAdvXX%D_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3_HorAdvXX%E_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3_HorAdvXX%F_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
                
                allocate(Me%COEF3_HorAdvYY%C_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3_HorAdvYY%D_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3_HorAdvYY%E_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3_HorAdvYY%F_Flux    (ILB:IUB, JLB:JUB, KLB:KUB)) 
                
                Me%COEF3_HorAdvXX%C_Flux     = Null_real
                Me%COEF3_HorAdvXX%D_Flux     = Null_real
                Me%COEF3_HorAdvXX%E_Flux     = Null_real
                Me%COEF3_HorAdvXX%F_Flux     = Null_real 

                Me%COEF3_HorAdvYY%C_Flux     = Null_real
                Me%COEF3_HorAdvYY%D_Flux     = Null_real
                Me%COEF3_HorAdvYY%E_Flux     = Null_real
                Me%COEF3_HorAdvYY%F_Flux     = Null_real                                
            
            endif  
            
            if (.not. Me%AdvDiff_Explicit) then
                
#ifdef _USE_PAGELOCKED
                ! Allocate pagelocked memory to optimize CUDA transfers
                call Alloc3DPageLocked(Me%ObjCuda, Me%COEF3%DPtr, Me%COEF3%D, IUB + 1, JUB + 1, KUB + 1)        
                call Alloc3DPageLocked(Me%ObjCuda, Me%COEF3%EPtr, Me%COEF3%E, IUB + 1, JUB + 1, KUB + 1)        
                call Alloc3DPageLocked(Me%ObjCuda, Me%COEF3%FPtr, Me%COEF3%F, IUB + 1, JUB + 1, KUB + 1)        
#else
                allocate(Me%COEF3%D                 (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB))
                allocate(Me%COEF3%E                 (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB))
                allocate(Me%COEF3%F                 (ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB))
#endif _USE_PAGELOCKED

!                allocate(Me%VECG                    (IJKLB:IJKUB))
!                allocate(Me%VECW                    (IJKLB:IJKUB))

                Me%COEF3%D                  = 0.0
                Me%COEF3%E                  = 1.0
                Me%COEF3%F                  = 0.0
                Me%TICOEF3                  = 0.0
!                Me%VECG                     = Null_real 
!                Me%VECW                     = Null_real 
           
            endif
            
            if (Me%Output%Boxes_ON) then
                allocate(Me%Fluxes%AdvFluxX    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Fluxes%AdvFluxY    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Fluxes%AdvFluxZ    (ILB:IUB, JLB:JUB, KLB:KUB))
                
                allocate(Me%Fluxes%DifFluxX    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Fluxes%DifFluxY    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Fluxes%DifFluxZ    (ILB:IUB, JLB:JUB, KLB:KUB))
                
                allocate(Me%Fluxes%MassFluxesX (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Fluxes%MassFluxesY (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%Fluxes%MassFluxesZ (ILB:IUB, JLB:JUB, KLB:KUB))
            
                Me%Fluxes%AdvFluxX     = Null_real
                Me%Fluxes%AdvFluxY     = Null_real
                Me%Fluxes%AdvFluxZ     = Null_real
                
                Me%Fluxes%DifFluxX     = Null_real
                Me%Fluxes%DifFluxY     = Null_real
                Me%Fluxes%DifFluxZ     = Null_real
                
                Me%Fluxes%MassFluxesX  = Null_real
                Me%Fluxes%MassFluxesY  = Null_real
                Me%Fluxes%MassFluxesZ  = Null_real
            
            endif
            
            !griflet: BEGIN this is the alternate version that allows parallel openmp
            IJKLB = min (ILB, JLB, KLB)
            IJKUB = max (IUB, JUB, KUB)

            Me%MaxThreads = 1
            !$ Me%MaxThreads = omp_get_max_threads()

            allocate(Me%THOMAS)
            allocate(Me%THOMAS%COEF3)
            allocate(Me%THOMAS%VEC(1:Me%MaxThreads))

            do m = 1, Me%MaxThreads
                
                VECGW => Me%THOMAS%VEC(m)

                allocate(VECGW%G(IJKLB:IJKUB))
                allocate(VECGW%W(IJKLB:IJKUB))

            enddo

            Me%THOMAS%COEF3%D => Me%COEF3%D
            Me%THOMAS%COEF3%E => Me%COEF3%E
            Me%THOMAS%COEF3%F => Me%COEF3%F
            Me%THOMAS%TI      => Me%TICOEF3

            !griflet: END

        endif
        
        if (Me%CalculateECw) then
            allocate (Me%ECw(ILB:IUB,JLB:JUB,KLB:KUB))
        endif
        
        
    end subroutine AllocateVariables

    !--------------------------------------------------------------------------
  
    subroutine ConstructScalar3D(Scalar3D, ExtractType, ClientNumber, block_begin, block_end)

        !Arguments-------------------------------------------------------------
        type(T_Property_3D), pointer        :: Scalar3D
        integer, intent(in)                 :: ExtractType
        integer, intent(in), optional       :: ClientNumber
        character(len=*)                    :: block_begin, block_end

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        logical                             :: BlockFound
        integer                             :: BlockClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        
        !----------------------------------------------------------------------
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB


        select case(ExtractType)

            case(FromBlock)

                call ExtractBlockFromBuffer(Me%ObjEnterData, BlockClientNumber, block_begin, block_end,  &
                                            BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR01'


            case(FromBlockInBlock)

                if(.not. present(ClientNumber))then
                    stop 'ConstructScalar3D - ModuleSoilProperties - ERR02'
                end if
                
                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber, block_begin, block_end,  &
                                           BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR03'

        end select

        if(BlockFound)then

            allocate(Scalar3D%Field(ILB:IUB, JLB:JUB, KLB:KUB))

            call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR03.1'

            call ConstructFillMatrix  (PropertyID           = Scalar3D%ID,                      &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = ExtractType,                      &
                                       PointsToFill3D       = Me%ExtVar%WaterPoints3D,          &
                                       Matrix3D             = Scalar3D%Field,                   &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR04'


            call GetDefaultValue(Scalar3D%ID%ObjFillMatrix, Scalar3D%Scalar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR05'

            call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR5.1'

            call KillFillMatrix(Scalar3D%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR06'

            if(ExtractType == FromBlockInBlock)then
                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR07'
            end if


            if(ExtractType == FromBlock)then
                call Block_Unlock(Me%ObjEnterData, BlockClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR08'

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR09'
            end if
        
        else
            write(*,*) 'Block not present:', block_begin, block_end
            stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR10'
        end if

   
    end subroutine ConstructScalar3D

    !--------------------------------------------------------------------------
  
    subroutine Construct_PropertyList

        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_Property), pointer          :: NewProperty
        type (T_Property), pointer          :: PropertyX 
        integer                             :: Index       

        !------------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                    &
                                        ClientNumber    = ClientNumber,     &
                                        block_begin     = prop_block_begin, &
                                        block_end       = prop_block_end,   &
                                        BlockFound      = BlockFound,       &
                                        STAT            = STAT_CALL)
cd1 :       if (STAT_CALL .EQ. SUCCESS_) then    

cd2 :           if (BlockFound) then                                                  
                    
                    !Construct a New Property 
                    Call Construct_Property(NewProperty)

                    !Add new Property to the SoilProperties List 
                    Call Add_Property(NewProperty)

                else cd2

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_PropertyList - ModulePorousMediaProeprties - ERR01'
                    exit do1    !No more blocks
                
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_PropertyList - ModulePorousMediaProeprties - ERR02'
            
            else cd1
                
                stop 'Construct_PropertyList - ModulePorousMediaProeprties - ERR03'
            
            end if cd1
        
        end do do1

        !Allocate space for PropertiesList
        nullify (Me%PropertiesList)
        allocate (Me%PropertiesList(Me%PropertiesNumber))
        
        !Fill PropertiesList with the ID of all properties currently in use
        Index = 1
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
        
            Me%PropertiesList(Index) = PropertyX%ID%IDNumber
            
            Index = Index + 1
            PropertyX => PropertyX%Next
        
        enddo

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------------    
    
    subroutine Construct_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer           :: NewProperty

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_)stop 'Construct_Property - ModulePorousMediaProeprties - ERR00'
        
        nullify(NewProperty%Prev, NewProperty%Next)
        nullify(NewProperty%Concentration)
        nullify(NewProperty%ConcentrationOnInfColumn)

        call ConstructPropertyID            (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call Construct_PropertyState        (NewProperty)

        call Construct_PropertyValues       (NewProperty)

        call Construct_PropertyEvolution    (NewProperty)

        call Construct_PropertyOutPut       (NewProperty)

    end subroutine Construct_Property
    
    !-------------------------------------------------------------------------------    
    
    subroutine Add_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),              pointer     :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstProperty)) then
            Me%PropertiesNumber     = 1
            Me%FirstProperty        => NewProperty
            Me%LastProperty         => NewProperty
        else
            NewProperty%Prev        => Me%LastProperty
            Me%LastProperty%Next    => NewProperty
            Me%LastProperty         => NewProperty
            Me%PropertiesNumber     = Me%PropertiesNumber + 1
        end if 


    end subroutine Add_Property 

    !-------------------------------------------------------------------------- 

    subroutine Construct_PropertyState(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL, iflag

        !----------------------------------------------------------------------
        

        !<BeginKeyword>
            !Keyword          : PARTICULATE
            !<BeginDescription>
            !<EndDescription>
            !Type             : logical   
            !Default          : Dissolved
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : From Block
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Particulate,                                            &
                     Me%ObjEnterData,  iflag,                                            &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'PARTICULATE',                                       &
                     ClientModule = 'ModulePorousMediaProeprties',                       &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR010'

        if (NewProperty%Particulate)then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR020'
            end if
        endif
        
        call GetData(NewProperty%UseToCalcECw,                                           &
                     Me%ObjEnterData,  iflag,                                            &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'USE_TO_CALC_ECW',                                   &
                     default      = .false.,                                             &
                     ClientModule = 'ModulePorousMediaProeprties',                       &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR030'

        if (NewProperty%UseToCalcECw) then
            call GetData(NewProperty%ECwFactor,                                              &
                         Me%ObjEnterData,  iflag,                                            &
                         SearchType   = FromBlock,                                           &
                         keyword      = 'ECW_FACTOR',                                        &
                         ClientModule = 'ModulePorousMediaProeprties',                       &
                         STAT         = STAT_CALL)
            if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR040'
            Me%CalculateECw = .true.
        endif
        
    end subroutine Construct_PropertyState

    !--------------------------------------------------------------------------

    subroutine Construct_PropertyEvolution(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: iflag, Decayflag
        real                                        :: ErrorAux, AuxFactor, DTAux
        real                                        :: ModelDT
        !----------------------------------------------------------------------

        call GetData(NewProperty%Evolution%MaxForUptakeConc,                             &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'MAX_CONC_FOR_UPTAKE',                               &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = 0.0,                                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR010'
        if (iflag .NE. 1) then
            NewProperty%Evolution%UseMaxForUptakeConc = .false.
        else
            NewProperty%Evolution%UseMaxForUptakeConc = .true.
        end if

        call GetData(NewProperty%Evolution%AdvectionDiffusion,                           &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'ADVECTION_DIFFUSION',                               &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR020'
        
        !particulate properties (defined by user or by MOHID) can not have movement in soil 
        !(particulates do not enter or leave soil) ModuleBasin will not give the concentration 
        !of this species (it will remain FillValueReal) and advection in soil has problems
        if (NewProperty%Evolution%AdvectionDiffusion) then 
            if (NewProperty%Particulate) then
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' has PARTICULATE option ON'
                write(*,*) 'and can not have ADVECTION_DIFFUSION ON in Porus Media Properties'
                stop 'Construct_PropertyEvolution - ModulePorousMediaProeprties - ERR030'      
            elseif (Check_Particulate_Property(NewProperty%ID%IDNumber)) then
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' has not PARTICULATE option ON'
                write(*,*) 'but is recognized by the model as being particulate tupe'
                write(*,*) 'and can not have ADVECTION_DIFFUSION ON in Porus Media Properties'
                stop 'Construct_PropertyEvolution - ModulePorousMediaProeprties - ERR035' 
            endif     
        endif  

        if (NewProperty%Evolution%AdvectionDiffusion) then
            Me%Coupled%AdvectionDiffusion  = .true.
            NewProperty%Evolution%Variable = .true.
            
            call ReadAdvectionDiffusionParam (NewProperty)
            call ConstructPropertyDiffusivity (NewProperty)
        endif
        
!        if (NewProperty%Evolution%AdvectionDiffusion) then
!
!            call ReadAdvectionDiffusionParam (NewProperty)
!        
!        end if
!
!        call ConstructPropertyDiffusivity (NewProperty)
        
        !properties that are transported in transpiration and evaporation 
        !except nitrate and inorganic phosphorus that are accounted already (e.g. temperature and oxygen)
        !if the flux is not accounted than the properties increase when transpiration and evaporation occurs
        if (NewProperty%Evolution%AdvectionDiffusion) then
            call GetData(NewProperty%Evolution%TransportedInEVTP,                            &
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromBlock,                                           &
                         keyword      = 'TRANSPORTED_IN_EVTP',                               &
                         ClientModule = 'ModulePorousMediaProperties',                       &
                         Default      = .false.,                                             &
                         STAT         = STAT_CALL)
            if(STAT_CALL .NE. SUCCESS_)                                                      &
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR037'
            
            !do not account for nitrate and inorganic phosphorus that are already accounted 
            !these properties are special because the uptake may not be Q*C (as SWAT generates)
            !and as reality roots can create gradients to assimilate more nutrients    
            if ((NewProperty%Evolution%TransportedInEVTP)       &
                 .and. ((NewProperty%ID%IDNumber == Nitrate_)    &
                        .or. (NewProperty%ID%IDNumber == Inorganic_Phosphorus_))) then
                write(*,*) 'Using TRANSPORTED_IN_EVTP option in nitrate and/or inorganic phosphorus property'
                write(*,*) 'These properties cant have this option active'
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR038'                
            endif                
        endif
        
        !By default the boundary value is imposed (this is only used for "wall" level imposed
        !in bottom flux it is free flux trough bottom and does not come up so does not need
        !exterior prop
        call GetData(NewProperty%Evolution%Boundary%BoundaryCondition,                   &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'BOUNDARY_CONDITION',                                &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = ImposedValue_,                                       &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR039'

        if (NewProperty%Evolution%Boundary%BoundaryCondition /= ImposedValue_   .and.    &
            NewProperty%Evolution%Boundary%BoundaryCondition /= NullGradient_        ) then 
            write(*,*) ' Boundary Condition can only be ImposedValue = 1 or'
            write(*,*) ' NullGradient = 2. Check BOUNDARY_CONDITION keyword'
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR40'
        endif
        
        !By default if not given the property enters with zero concentration.
        call GetData(NewProperty%Evolution%Boundary%DefaultBoundary,                     &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DEFAULTBOUNDARY',                                   &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = 0.0,                                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR041'
        
        !<BeginKeyword>
            !Keyword          : SOIL_QUALITY
            !<BeginDescription>
               ! Property has the Soil quality model (sediment quality) as sink and source
            !<EndDescription>
            !Type             : Boolean
            !Default          : .false.
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%SoilQuality,                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SOIL_QUALITY',                                      &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR045'

        if (NewProperty%Evolution%SoilQuality) then
            Me%Coupled%SoilQuality     = .true.
            NewProperty%Evolution%Variable = .true.
        endif


#ifdef _PHREEQC_
        !<BeginKeyword>
            !Keyword          : SOIL_CHEMISTRY
            !<BeginDescription>
               ! Property has the Soil chemistry model (PHREEQC)
            !<EndDescription>
            !Type             : Boolean
            !Default          : .false.
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%SoilChemistry,                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SOIL_CHEMISTRY',                                    &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR050'

        if (NewProperty%Evolution%SoilChemistry) then
            Me%Coupled%SoilChemistry       = .true.
            NewProperty%Evolution%Variable = .true.
            
!            call ReadChemistryParameters (NewProperty)
        endif
#endif

        !<BeginKeyword>
            !Keyword          : PARTITION
            !<BeginDescription>       
               ! This property has partition as a sink and source
            !<EndDescription>
            !Type             : Logical 
            !Default          : FALSE
            !File keyword     : DISPQUAL
            !Multiple Options : NO, WQM
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>       
        call GetData(NewProperty%Evolution%Partitioning,                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'PARTITION',                                         &
                     Default      = .false.,                                             &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Subroutine Construct_PropertyEvolution - ModuleRunoffProperties - ERR060'

        if(NewProperty%Evolution%Partitioning) then
            NewProperty%Evolution%Variable = .true.
            Me%Coupled%Partition           = .true.
        endif
        
        if(NewProperty%Evolution%Partitioning)                                           &
            call Read_Partition_Parameters(NewProperty)


        call GetData(NewProperty%Evolution%Decay,                                        &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DECAY',                                             &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR60'

        if (NewProperty%Evolution%Decay) then
            Me%Coupled%Decay       = .true.
        endif

        if (NewProperty%Evolution%Decay) then

            !to output decay rate computed (g/m3.day-1 or mg/kg.day-1)
            allocate(NewProperty%PropertyDecay(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,              &
                                               Me%Size%KLB:Me%Size%KUB), STAT = STAT_CALL)
            if(STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModulePorousMediaProeprties - ERR68'

           
            NewProperty%Evolution%DecayEquationCoupled = .false.             

            !Decay equation. For now exists Peyrard model for denitrification but structure is created
            !for introducing other equations
            call GetData(NewProperty%Evolution%DecayEquation,                                &
                         Me%ObjEnterData,Decayflag,                                          &
                         SearchType   = FromBlock,                                           &
                         keyword      = 'DECAY_EQUATION',                                    &
                         ClientModule = 'ModulePorousMediaProperties',                       &
                         default      = Peyrard_,                                            &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                     &
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR65'            
            
            !if flag different from zero keyword is present
            if (Decayflag /= 0) then
                if (NewProperty%Evolution%DecayEquation /= Peyrard_) then
                    write(*,*)
                    write(*,*) 'In Porous Media Properties the DECAY_EQUATION in property block'
                    write(*,*) 'allowed values are: 1 - Peyrard'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR67'
                endif            
                
                NewProperty%Evolution%DecayEquationCoupled = .true.
                                
                !search for property half sauration
                if (NewProperty%Evolution%DecayEquation == Peyrard_ .and. NewProperty%ID%IDNumber == Nitrate_) then
                
                    !Property half saturation value for decay (same units as property)
                    call GetData(NewProperty%Evolution%DecayHalfSaturation,                          &
                                 Me%ObjEnterData,iflag,                                              &
                                 SearchType   = FromBlock,                                           &
                                 keyword      = 'DECAY_HALF_SATURATION',                             &
                                 ClientModule = 'ModulePorousMediaProperties',                       &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                                     &
                        stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR69' 
                    
                    if (NewProperty%Evolution%DecayHalfSaturation .le. 0.) then
                        write(*,*)
                        write(*,*) 'In Porous Media Properties the DECAY_HALF_SATURATION in property block'
                        write(*,*) 'for property', trim(NewProperty%ID%Name)
                        write(*,*) 'needs to be defined and > 0'
                        stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR70'                        
                    endif
                endif
                
            else !flag = 0 -> keyword not present
                !if decay equation keyword is not present for the property than search for the decay rate value
                !and first order decay will be implemented
                                      
                !Decay rate k (day-1) in first order decay P = Po*exp(-kt)
                call GetData(NewProperty%Evolution%DecayRate,                                    &
                             Me%ObjEnterData,iflag,                                              &
                             SearchType   = FromBlock,                                           &
                             keyword      = 'DECAY_RATE',                                        &
                             ClientModule = 'ModulePorousMediaProperties',                       &
                             default      = 0.,                                                  &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                     &
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR71'            

                !Apply first order decay to mass or concentration?
                call GetData(NewProperty%Evolution%DecayMass,                                    &
                             Me%ObjEnterData,iflag,                                              &
                             SearchType   = FromBlock,                                           &
                             keyword      = 'DECAY_MASS',                                        &
                             ClientModule = 'ModulePorousMediaProperties',                       &
                             default      = .false.,                                             &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                     &
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR72'   
            
            endif
            
        endif
        
        call GetData(NewProperty%Evolution%Discharges,                              &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'DISCHARGES',                                 &
                     ClientModule   = 'ModulePorousMediaProperties',                &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR80' 
        
        if (NewProperty%Evolution%Discharges) then
            Me%Coupled%Discharges = .true.
            Me%nPropWithDischarge = Me%nPropWithDischarge + 1
        endif        

        !Property time step
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR080'

            ModelDT = Me%ExtVar%DT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DTINTERVAL',                                    &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModulePorousMediaProperties',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR090'
                                       

            if (iflag == 1) then
                NewProperty%Evolution%DTIntervalAssociated = .true.
                Me%DTIntervalAssociated                    = .true.                           

                if (NewProperty%Evolution%DTInterval < ModelDT) then
                    write(*,*) 
                    write(*,*) 'Property time step is smaller then model time step'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR100'

                elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                    !Property time step must be a multiple of the model time step
                    auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                    Erroraux = auxFactor - int(auxFactor)
                    if (Erroraux /= 0) then
                        write(*,*) 
                        write(*,*) 'Property time step must be a multiple of model time step.'
                        write(*,*) 'Please review your input data.'
                        stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR110'
                    endif

                    !Run period in seconds
                    DTaux = Me%ExtVar%EndTime - Me%ExtVar%Now

                    !The run period   must be a multiple of the Property DT
                    auxFactor = DTaux / NewProperty%Evolution%DTInterval

                    ErrorAux = auxFactor - int(auxFactor)
                    if (ErrorAux /= 0) then

                        write(*,*) 
                        write(*,*) 'Property time step is not a multiple of model time step.'
                        stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR120'
                    end if
                endif

                NewProperty%Evolution%NextCompute = Me%ExtVar%Now + NewProperty%Evolution%DTInterval
            
            endif
            
        else

            call null_time(NewProperty%Evolution%NextCompute)

            NewProperty%Evolution%DTInterval = FillValueReal

        endif

    end subroutine Construct_PropertyEvolution     

    !--------------------------------------------------------------------------

    subroutine Construct_SedimentRateList

        !External----------------------------------------------------------------
        integer                                :: ClientNumber
        integer                                :: STAT_CALL
        logical                                :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_SedimentRate),    pointer      :: NewSedimentRate

        !------------------------------------------------------------------------
        
        Me%Output%RateFluxes = .false.
        
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                    &
                                        ClientNumber    = ClientNumber,     &
                                        block_begin     = '<beginSQrate>',  &
                                        block_end       = '<endSQrate>',    &
                                        BlockFound      = BlockFound,       &
                                        STAT            = STAT_CALL)
            if(STAT_CALL .EQ. SUCCESS_)then    
                if (BlockFound) then                                                  
                    
                    Me%Output%RateFluxes = .true.
                    
                    !Construct a New Sediment Rate
                    Call Construct_SedimentRate(NewSedimentRate)

                    !Add new Rate to the Sediment Rates List 
                    Call Add_SedimentRate(NewSedimentRate)

                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Construct_SedimentRateList - ModulePorousMediaProperties - ERR01'

                    exit do1    !No more blocks
                end if


            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_SedimentRateList - ModulePorousMediaProperties - ERR02'
            else
                stop       'Construct_SedimentRateList - ModulePorousMediaProperties - ERR03'
            end if
        end do do1
         
    end subroutine Construct_SedimentRateList

    !----------------------------------------------------------------------------

    subroutine Construct_SedimentRate(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_SedimentRate), pointer       :: NewSedimentRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewSedimentRate, STAT = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                          &
            stop 'Construct_SedimentRate - ModulePorousMediaProperties - ERR01' 

        nullify(NewSedimentRate%Field, NewSedimentRate%Prev, NewSedimentRate%Next)

        call Construct_SedimentRateID       (NewSedimentRate)

        call Construct_SedimentRateValues   (NewSedimentRate)


    end subroutine Construct_SedimentRate

     !--------------------------------------------------------------------------
    
    !This subroutine reads all the information needed to construct the property ID          
    subroutine Construct_SedimentRateID(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_SedimentRate), pointer       :: NewSedimentRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: iflag, PropNumber
        logical                             :: CheckName
        type (T_Property), pointer          :: PropertyX
      
        !----------------------------------------------------------------------
        
        !First Property defined in a rate relation
        call GetData(NewSedimentRate%FirstProp%name,                                    &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'FIRSTPROP',                                        &
                     ClientModule = 'ModulePorousMediaProperties',                      &
                     SearchType   = FromBlock,                                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR01' 
        if (iflag==0)                                                                   &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR02' 

        !Check if the property name is valid
        CheckName = CheckPropertyName(NewSedimentRate%FirstProp%Name, Number = PropNumber)
        if (CheckName) then
            NewSedimentRate%FirstProp%IDnumber = PropNumber
        else
            write(*,*)
            write(*,*) 'The first property name is not recognised by the model.'
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR03' 
        end if 

        call Search_Property(PropertyX, PropertyXID = PropNumber, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                     &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR04' 
        
        !second Property defined in a rate relation
        call GetData(NewSedimentRate%SecondProp%name,                                    &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'SECONDPROP',                                        &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     SearchType   = FromBlock,                                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR05' 
        if (iflag==0)                                                                    &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR06' 
      
        ! Check if the property name is valid OR not
        CheckName = CheckPropertyName(NewSedimentRate%SecondProp%name, Number = PropNumber)
        if (CheckName) then
            NewSedimentRate%SecondProp%IDnumber = PropNumber
        else
            write(*,*)
            write(*,*) 'The Second property name is not recognised by the model.'
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR07' 
        end if
        
        call Search_Property(PropertyX,PropertyXID = PropNumber, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR08' 
  
  
        !Rate description ex: zooplankton grazing over phytoplankton
        call GetData(NewSedimentRate%ID%Description,                                     &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'DESCRIPTION',                                       &
                     Default      = 'No description was given.',                         &
                     ClientModule = 'ModulePorousMediaProperties',                          &
                     SearchType   = FromBlock,                                           &
                     STAT         = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR09' 

        call GetData(NewSedimentRate%ID%Name,                                              &
                     Me%ObjEnterData, iflag,                                               &
                     keyword      = 'NAME',                                                &
                     ClientModule = 'ModulePorousMediaProperties',                            &
                     SearchType   = FromBlock,                                             &
                     Default      = 'No name was given to sediment rate.',                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                         &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR10' 
        if (iflag==0)                                                                      &
            stop 'Construct_SedimentRateID - ModulePorousMediaProperties - ERR11' 

    end subroutine Construct_SedimentRateID

    !--------------------------------------------------------------------------

    subroutine Construct_SedimentRateValues(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_Sedimentrate), pointer       :: NewSedimentRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        allocate(NewSedimentRate%Field(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_SedimentRateValues - ModulePorousMediaProperties - ERR01' 
        NewSedimentRate%Field(:,:,:) = FillValueReal

    end subroutine Construct_SedimentRateValues
    
    !--------------------------------------------------------------------------

    subroutine Add_SedimentRate(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_SedimentRate), pointer       :: NewSedimentRate

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstSedimentRate)) then
            Me%SedimentRatesNumber      = 1
            Me%FirstSedimentRate        => NewSedimentRate
            Me%LastSedimentRate         => NewSedimentRate
        else
            NewSedimentRate%Prev        => Me%LastSedimentRate
            Me%LastSedimentRate%Next    => NewSedimentRate
            Me%LastSedimentRate         => NewSedimentRate
            Me%SedimentRatesNumber      = Me%SedimentRatesNumber + 1
        end if 

    end subroutine Add_SedimentRate 

    
    !--------------------------------------------------------------------------

    subroutine Read_Partition_Parameters(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: iflag
        real                                    :: DefaultFraction
        
        !Begin-----------------------------------------------------------------

        
        !<BeginKeyword>
            !Keyword          : PARTITION_FRACTION
            !<BeginDescription>       
               ! 
               ! Partition fraction
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 0.0
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        
        if(NewProperty%Particulate)then
            DefaultFraction = 0.9
        else
            DefaultFraction = 0.1
        end if

        call GetData(NewProperty%Evolution%Partition%Fraction,                          &
                     Me%ObjEnterData, iflag,                                            &
                     keyword    = 'PARTITION_FRACTION',                                 & 
                     default    = DefaultFraction,                                      &
                     SearchType = FromBlock,                                            &
                     ClientModule = 'ModuleRunoffProperties',                           &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR10'

        !<BeginKeyword>
            !Keyword          : PARTITION_RATE
            !<BeginDescription>       
               ! 
               ! Partition transfer rate between the particulate and the dissolved phase
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 1.0
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%Partition%Rate,                              &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'PARTITION_RATE',                                   & 
                     default      = 1.,                                                 &
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR20'


        !<BeginKeyword>
            !Keyword          : USE_SED_REF_CONC
            !<BeginDescription>       
               ! 
               ! Use Reference cohesive sediment concentration method
               ! 
            !<EndDescription>
            !Type             : Real 
            !Default          : 1.0
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%Partition%UseSedimentRefConc,                 &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'USE_SED_REF_CONC',                                  & 
                     default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR30'

        if(NewProperty%evolution%Partition%UseSedimentRefConc)then

            !<BeginKeyword>
                !Keyword          : SED_REF_CONC
                !<BeginDescription>       
                   ! 
                   ! Reference cohesive sediment concentration 
                   ! 
                !<EndDescription>
                !Type             : Real 
                !Default          : 1.0
                !File keyword     : DISPQUAL
                !Multiple Options : Do not have
                !Search Type      : FromBlock
                !Begin Block      : <beginproperty>
                !End Block        : <endproperty>
            !<EndKeyword>

            call GetData(NewProperty%Evolution%Partition%SedimentRefConc,               &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'SED_REF_CONC',                                 & 
                         default      = 1.,                                             &
                         SearchType   = FromBlock,                                      &
                         ClientModule = 'ModuleRunoffProperties',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR30'

        end if
        
        !<BeginKeyword>
            !Keyword          : PARTITION_COUPLE
            !<BeginDescription>       
               ! 
               ! Name of property (dissolved/particulated) to couple  
               ! 
            !<EndDescription>
            !Type             : Character 
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%Partition%Couple,                             &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'PARTITION_COUPLE',                                  & 
                     ClientModule ='ModuleRunoffProperties',                              &
                     SearchType   = FromBlock,                                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR40'
        if (iflag .NE. 1)                                                                &
            stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR50'

!        !<BeginKeyword>
!            !Keyword          : SALINITY_EFFECT
!            !<BeginDescription>       
!               ! Verifies if the user wants to compute partition coefficient between the 
!               ! particulate and the dissolved phase function of salinity
!               ! 
!            !<EndDescription>
!            !Type             : Boolean 
!            !Default          : .false.
!            !File keyword     : DISPQUAL
!            !Multiple Options : Do not have
!            !Search Type      : FromBlock
!            !Begin Block      : <beginproperty>
!            !End Block        : <endproperty>
!        !<EndKeyword>
!        call GetData(NewProperty%Evolution%Partition%SalinityEffect,                     &
!                     Me%ObjEnterData, iflag,                                             &
!                     keyword      = 'SALINITY_EFFECT',                                   & 
!                     default      = OFF,                                                 &
!                     SearchType   = FromBlock,                                           &
!                     ClientModule = 'ModuleRunoffProperties',                             &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR60'
!
!
!cd1:    if (NewProperty%Evolution%Partition%SalinityEffect) then
!
!            !<BeginKeyword>
!                !Keyword          : EMPIRIC_COEF
!                !<BeginDescription>       
!                   ! 
!                   ! Do not have
!                   ! 
!                !<EndDescription>
!                !Type             : Real 
!                !Default          : -14.505
!                !File keyword     : DISPQUAL
!                !Multiple Options : Do not have
!                !Search Type      : FromBlock
!                !Begin Block      : <beginproperty>
!                !End Block        : <endproperty>
!            !<EndKeyword> 
!            call GetData(NewProperty%Evolution%Partition%EmpiricCoef,                    &
!                         Me%ObjEnterData, iflag,                                         &
!                         keyword      = 'EMPIRIC_COEF',                                  & 
!                         default      = -14.505 ,                                        &
!                         SearchType   = FromBlock,                                       &
!                         ClientModule = 'ModuleRunoffProperties',                         &
!                         STAT         = STAT_CALL)
!            if (STAT_CALL .NE. SUCCESS_)                                                 &
!                stop 'Read_Partition_Parameters - ModuleRunoffProperties - ERR70'
!
!        endif cd1


    end subroutine Read_Partition_Parameters 
    
    !-------------------------------------------------------------------------- 
    
    subroutine ReadAdvectionDiffusionParam (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !External--------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: MassConservation
        integer                   :: ImposedValue
        integer                   :: NullGradient, CyclicBoundary
        integer                   :: Orlanski, MassConservNullGrad

        !Local-----------------------------------------------------------------
        integer                   :: iflag !, BoundaryCondition

        !----------------------------------------------------------------------

        call GetData(NewProperty%Evolution%AdvDiff%Molecular_Diff_Coef, &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'ADVDIFF_MOLECULAR_DIFF_COEF',      &
                     Default      = 0.0,                                &
                     ClientModule = 'ModulePorousMediaProperties',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR10'

            
        !This keyword has no effect because is not used in ModuleAdvectionDiffusion
        call GetData(NewProperty%Evolution%AdvDiff%NumericStability, &
                     Me%ObjEnterData, iflag,                         &
                     SearchType   = FromBlock,                       &
                     keyword      = 'ADVDIFF_NUM_STABILITY',         &
                     Default      = .FALSE.,                         &
                     ClientModule = 'ModulePorousMediaProperties',   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR20'
        
        !Parameters to compute diffusivity - not needed because we are computing it according to Jury
        call GetData(NewProperty%Evolution%AdvDiff%SchmidtNumberH, &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_SCHMIDT_NUMBER_H',    &
                     Default      = 1.0,                           &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR30'

        !Parameters to compute diffusivity - not needed because we are computing it according to Jury
        call GetData(NewProperty%Evolution%AdvDiff%SchmidtCoefV,   &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_SCHMIDT_COEF_V',      &
                     Default      = 1.0,                           &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR40'

        !Parameters to compute diffusivity - not needed because we are computing it according to Jury
        call GetData(NewProperty%Evolution%AdvDiff%SchmidtBackgroundV, &
                     Me%ObjEnterData, iflag,                           &
                     SearchType   = FromBlock,                         &
                     keyword      = 'ADVDIFF_SCHMIDT_BACKGROUND_V',    &
                     Default      = 0.,                                &
                     ClientModule = 'ModulePorousMediaProperties',     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR50'

        !Zero diffusivity - not needed because we are computing it according to Jury
        call GetData(NewProperty%Evolution%AdvDiff%NullDif,        &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_NULLDIF',             &
                     Default      = .false.,                       &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR60'
        
        !Open Boundary condition - in runoff is known
        call GetBoundaryConditionList(MassConservation    = MassConservation,    &
                                      ImposedValue        = ImposedValue,        &
                                      NullGradient        = NullGradient,        &
                                      Orlanski            = Orlanski,            &
                                      MassConservNullGrad = MassConservNullGrad, &
                                      CyclicBoundary      = CyclicBoundary)

!        call GetData(BoundaryCondition,                            &
!                     Me%ObjEnterData,  iflag,                      &
!                     SearchType   = FromBlock,                     &
!                     keyword      = 'ADVDIFF_BOUNDARY_CONDITION',  &
!                     Default      = MassConservation,              &
!                     ClientModule = 'ModulePorousMediaProperties', &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR70'
!
!        ! By default it's imposed a value dependent only from the exterior
!        ! value and of the decay time. However this method doesn't conserve mass
!        ! when the water fluxes near the frontier are dominant
!
!        if (BoundaryCondition /= MassConservation     .and. &
!            BoundaryCondition /= ImposedValue         .and. &
!            BoundaryCondition /= NullGradient         .and. &
!            BoundaryCondition /= CyclicBoundary       .and. &
!            BoundaryCondition /= Orlanski             .and. &
!            BoundaryCondition /= MassConservNullGrad) &
!            stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR80'
!
!        NewProperty%Evolution%AdvDiff%BoundaryCondition = BoundaryCondition

        !By default the horizontal Diffusion discretization is explicit
        NewProperty%Evolution%AdvDiff%DiffusionH_imp_exp  = ExplicitScheme

        NewProperty%Evolution%AdvDiff%ImplicitH_Direction = DirectionX
        
        !Spatial Method for horizontal advection
        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodH,     &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_METHOD_H',            &
                     Default      = UpwindOrder1,                  &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR90'

        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationH, &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_TVD_LIMIT_H',         &
                     Default      = Superbee,                      &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR100'

        !Spatial Method for horizontal advection
        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodV,     &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_METHOD_V',            &
                     Default      = UpwindOrder1,                  &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR110'

        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationV, &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_TVD_LIMIT_V',         &
                     Default      = Superbee,                      &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR120'


        call GetData(NewProperty%Evolution%AdvDiff%VolumeRelMax,   &
                     Me%ObjEnterData, iflag,                       &
                     Keyword      = 'ADVDIFF_VOLUME_RELATION_MAX', &
                     Default      = 5.,                            &
                     SearchType   = FromBlock,                     &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR130'


        if (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodH == P2_TVD) then
            NewProperty%Evolution%AdvDiff%Upwind2H = .true.
        else
            NewProperty%Evolution%AdvDiff%Upwind2H = .false.
        endif

        if (NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder2 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder3 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodV == P2_TVD) then
            NewProperty%Evolution%AdvDiff%Upwind2V = .true.
        else
            NewProperty%Evolution%AdvDiff%Upwind2V = .false.
        endif


        if (.not. Me%AdvDiff_Explicit .and.&
           (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3)) then

            write(*,*) 'If the advection of mass in the horizontal is implicit'
            write(*,*) 'the advection method can not be a second or third order upwind'
            stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR140.'

        endif

        if (.not. Me%AdvDiff_Explicit .and.&
           (NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder2 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder3)) then

            write(*,*) 'If the advection of mass in the vertical is implicit'
            write(*,*) 'the advection method can not be a second or third order upwind'
            stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR150.'

        endif


    end subroutine ReadAdvectionDiffusionParam

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyDiffusivity (NewProperty)

        !Arguments---------------------------------------------------------
        type(T_Property),    pointer                :: NewProperty

        !Local-------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: STAT_CALL
        
        !Begin-------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
        allocate (NewProperty%Diffusivity (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR010'
        NewProperty%Diffusivity       = 0. 
        
           
        allocate (NewProperty%Viscosity (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR020'
        NewProperty%Viscosity         = 0.
        
!        allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR30'
!        NewProperty%Diff_Turbulence_H = 0.
!
!        allocate (NewProperty%Diff_Turbulence_V (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR40'
!        NewProperty%Diff_Turbulence_V = 0.                   
        
        allocate (NewProperty%ViscosityU (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR050'
        
        NewProperty%ViscosityU  = 0.0
        
        allocate (NewProperty%ViscosityV (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR060'
    
        NewProperty%ViscosityV  = 0.0    
            
        allocate (NewProperty%UptakeActive (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR070'
        
        NewProperty%UptakeActive  = .false.
        
    end subroutine ConstructPropertyDiffusivity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine Construct_PropertyValues(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),              pointer      :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        integer                                     :: ILB,IUB
        integer                                     :: JLB,JUB
        integer                                     :: KLB,KUB
        integer                                     :: WorkSizeILB, WorkSizeIUB
        integer                                     :: WorkSizeJLB, WorkSizeJUB
        integer                                     :: WorkSizeKLB, WorkSizeKUB
        
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        WorkSizeILB = Me%WorkSize%ILB
        WorkSizeIUB = Me%WorkSize%IUB
        WorkSizeJLB = Me%WorkSize%JLB
        WorkSizeJUB = Me%WorkSize%JUB
        WorkSizeKLB = Me%WorkSize%KLB
        WorkSizeKUB = Me%WorkSize%KUB
        
#ifdef _USE_PAGELOCKED
        ! Allocate pagelocked memory to optimize CUDA transfers
        call Alloc3DPageLocked(Me%ObjCuda, NewProperty%ConcentrationPtr, NewProperty%Concentration, IUB + 1, JUB + 1, KUB + 1)
#else
        allocate(NewProperty%Concentration(ILB:Pad(ILB, IUB), JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR10'
#endif _USE_PAGELOCKED        
        NewProperty%Concentration(:,:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOld(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR20'
        NewProperty%ConcentrationOld(:,:,:) = FillValueReal

        !by default zero because not all properties need to be written from basin (e.g. property not existing in runoff)
        allocate(NewProperty%ConcentrationOnInfColumn(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR40'
        NewProperty%ConcentrationOnInfColumn(:,:) = 0.
        
        if (Me%ExtVar%CoupledDN) then
            allocate(NewProperty%ConcentrationDN(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR50'
            NewProperty%ConcentrationDN(:,:) = FillValueReal
            

        endif
        
        call GetData(NewProperty%Pesticide,                                                 &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'PESTICIDE',                                            &
                     default      = .false.,                                                &
                     ClientModule = 'ModulePorousMediaProperties',                          &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR60'

        !by default zero because user may copy property blocks and leave pesticide ON in a property that is not pesticide.
        if (NewProperty%Pesticide) then
            allocate(NewProperty%PesticideFlux(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR70'
            NewProperty%PesticideFlux(:,:) = 0.
        endif

        
        !it has to be always allocated
        allocate(NewProperty%ConcInInterfaceDN(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR80'
        NewProperty%ConcInInterfaceDN(:,:,:) = FillValueReal            

        allocate(NewProperty%ConcInBoundary(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR90'
        NewProperty%ConcInBoundary(:,:,:) = FillValueReal            
        
        call GetData(NewProperty%MinValue,                                                  &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'MIN_VALUE',                                            &
                     ClientModule = 'ModulePorousMediaProperties',                          &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR100'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
            Me%Coupled%MinConcentration = .true.
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        call GetData(NewProperty%Evolution%WarnOnNegativeValues,                         &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType     = FromBlock,                                         &
                     keyword        = 'WARN_ON_NEGATIVE_VALUES',                         &
                     Default        = .false.,                                           &                     
                     ClientModule   = 'ModulePorousMediaProperties',                     &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ModulePorousMediaProperties - ConstructPropertyValues - ERR105' 
        
        if (NewProperty%Evolution%WarnOnNegativeValues) Me%Coupled%WarnOnNegativeValues = .true.

!        if(NewProperty%Evolution%MinConcentration)then
        !Mass_Created is also used for when concentration is negative
        allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)&
            stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR110'
        NewProperty%Mass_Created(:,:,:) = 0.
!        endif
        
        !if vegetation nitrogen uptake than nitrate need mass created
        if ((NewProperty%ID%IDNumber == Nitrate_) .and. (Me%ExtVar%ModelNitrogen)       &
             .and. (.not. NewProperty%Evolution%MinConcentration)) then
            write(*,*) 'Using vegetation nitrogen uptake so MIN_VALUE should be active'
            write(*,*) 'in nitrate property in porous media properties'
            stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR111'
        endif

        !if vegetation phosphorus uptake than inorganic phosphorus need mass created
        if ((NewProperty%ID%IDNumber == Inorganic_Phosphorus_) .and. (Me%ExtVar%ModelPhosphorus)       &
             .and. (.not. NewProperty%Evolution%MinConcentration)) then
            write(*,*) 'Using vegetation phosphorus uptake so MIN_VALUE should be active'
            write(*,*) 'in inorganic phosphorus property in porous media properties'
            stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR112'
        endif
        
        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModulePorousMediaProperties',                      &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR120'
          
        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not. NewProperty%Old) then

            !Get water points
            call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR130'

            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill3D       = Me%ExtVar%WaterPoints3D,          &
                                       Matrix3D             = NewProperty%Concentration,        &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR140'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR0150'
            end if

            call SetMatrixValue(NewProperty%ConcentrationOld, Me%Size, NewProperty%Concentration,Me%ExtVar%WaterPoints3D)

            call CheckFieldConsistence (NewProperty)

            call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR160'

        else

            ! If the property is old then the program is going to try to find a property
            ! with the same name in the Water properties initial file written in HDF format  
            call ReadOldConcBoundariesHDF(NewProperty)

        end if   

    end subroutine Construct_PropertyValues

      !--------------------------------------------------------------------------

    subroutine ConstructPartition

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: Property
        type(T_Property), pointer                   :: DissolvedProperty
!        type(T_Property), pointer                   :: CohesiveSediment
        type(T_Property), pointer                   :: ParticulateProperty
        character(len=StringLength)                 :: PartPropName
        real                                        :: TotalPartition
        integer                                     :: Couple_ID
        integer                                     :: Error

        !Begin-----------------------------------------------------------------

        Property => Me%FirstProperty
            
do1:    do while(associated(Property))

            if (Property%Evolution%Partitioning .and. .not. Property%Particulate) then

                DissolvedProperty => Property 

                PartPropName = trim(DissolvedProperty%Evolution%Partition%Couple)
                
                !Get the ID of the associated particulate
                if (.not. CheckPropertyName(PartPropName, Couple_ID)) then
                    write(*,*)
                    write(*,*) 'Particulate property name is not recognised by the model'
                    stop       'ConstructPartition - ModulePorousMediaProperties - ERR010'
                else
                    
                    DissolvedProperty%Evolution%Partition%Couple_ID = Couple_ID

                end if
                
                !Search for the particulate
                call Search_Property(ParticulateProperty, PropertyXID = Couple_ID, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    write(*,*)
                    write(*,*) 'Particulate property was not found'
                    stop 'ConstructPartition - ModulePorousMediaProperties - ERR10'
                else
                    if(.not. Check_Particulate_Property (Couple_ID)) then
                        write(*,*)
                        write(*,*) 'Couple property', trim(GetPropertyName(Couple_ID))
                        write(*,*) 'is not recognized by the model as being particulate'
                        stop 'ConstructPartition - ModulePorousMediaProperties - ERR020'                        
                    endif
                endif
                
                !Check if coupled names of each other match
                if(trim(ParticulateProperty%Evolution%Partition%Couple) .ne.                &
                   trim(DissolvedProperty%ID%Name)) then
                    write(*,*)
                    write(*,*) 'Particulate property pair and dissolved property names do not match'
                    stop 'ConstructPartition - ModulePorousMediaProperties - ERR030'
                endif
                
                if(ParticulateProperty%Evolution%DTInterval .ne.                            &
                   DissolvedProperty%Evolution%DTInterval) then
                    write(*,*)
                    write(*,*) 'Particulate property and dissolved property DTs do not match'
                    stop 'ConstructPartition - ModulePorousMediaProperties - ERR040'
                endif
                
                if(DissolvedProperty%Evolution%Partition%Rate /=                            &
                   ParticulateProperty%Evolution%Partition%Rate   )then
                    write(*,*)'Particulate and dissolved phases must have equal partition rates'
                    stop 'ConstructPartition - ModulePorousMediaProperties - ERR050'
                end if

                TotalPartition = DissolvedProperty%Evolution%Partition%Fraction  +          &
                                 ParticulateProperty%Evolution%Partition%Fraction

                Error = abs(1. - TotalPartition)
                     
                if(Error > 0.001)then
                    write(*,*)'Particulate and dissolved phases fractions must sum iqual to 1.'
                    stop 'ConstructPartition - ModulePorousMediaProperties - ERR060'
                end if

                ! .EQV. = Logical equivalence: the expression is true if both A and B 
                !  are true, or both are false.  
!                if (.NOT. (DissolvedProperty%Evolution%Partition%SalinityEffect .EQV.       &
!                           ParticulateProperty%Evolution%Partition%SalinityEffect))         &
!                    stop 'ConstructPartition - ModuleWaterProperties - ERR60'

!                if (DissolvedProperty%Evolution%Partition%EmpiricCoef .ne.                  &
!                    ParticulateProperty%Evolution%Partition%EmpiricCoef )                   &
!                    stop 'ConstructPartition - ModuleWaterProperties - ERR70'

!                if (.NOT. (DissolvedProperty%Evolution%Partition%UseSedimentRefConc .EQV.   &
!                           ParticulateProperty%Evolution%Partition%UseSedimentRefConc)) then
!                    write(*,*)
!                    write(*,*) 'Both particulate property and dissolved property must use or not sediment ref conc'
!                    stop 'ConstructPartition - ModulePorousMediaProperties - ERR070'
!                endif
!
!                if(DissolvedProperty%Evolution%Partition%UseSedimentRefConc)then
!                
!                    call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) then                                      
!                        write(*,*)
!                        write(*,*) 'If using sediment ref conc need cohesive sediment property'
!                        stop 'ConstructPartition - ModulePorousMediaProperties - ERR080'
!                    endif
!                endif
!                
!                if(DissolvedProperty%Evolution%Partition%SedimentRefConc /=                 &
!                   ParticulateProperty%Evolution%Partition%SedimentRefConc   )then
!                    write(*,*)'Particulate and dissolved phases must have equal cohesive sediment'
!                    write(*,*)'reference concentration'
!                    stop 'ConstructPartition - ModulePorousMediaProperties - ERR090'
!                end if


            end if

            nullify(DissolvedProperty, ParticulateProperty)

            Property => Property%Next


        end do do1
          
        nullify(Property)


    end subroutine ConstructPartition

    !--------------------------------------------------------------------------

    subroutine ConstructProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL, iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        character(len=StringLength), dimension(:,:), pointer:: PropertyList
        integer                                             :: nProperties
        integer                                             :: n
        type (T_Property), pointer                          :: PropertyX

        nProperties = Me%PropertiesNumber 

        !Allocates PropertyList
        allocate(PropertyList(nProperties, 2))
       
        n=1
        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))

            !Fills up PropertyList
            PropertyList(n, 1) = trim(PropertyX%ID%Name)
            PropertyList(n, 2) = "m3/m3"

            n=n+1

            PropertyX=>PropertyX%Next

        enddo

        !----------------------------------------------------------------------

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMedia',                                &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMediaProperties - ERR02' 
        
        !Starts Profile for Theta / ThetaF
        call StartProfile  (ProfileID       = Me%ObjProfile,                            &
                            ObjTime         = Me%ObjTime,                               &
                            ProfileDataFile = trim(TimeSerieLocationFile),              &
                            WaterPoints2D   = Me%ExtVar%BasinPoints,                    &
                            nProperties     = Me%PropertiesNumber ,                     &
                            PropertyList    = PropertyList,                             &
                            KUB             = Me%WorkSize%KUB,                          &
                            ClientName      = "PorousMedia",                            &
                            STAT            = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMediaProperties - ERR03' 


        deallocate (PropertyList)
        
    end subroutine ConstructProfileOutput

    !---------------------------------------------------------------------------

    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),    pointer        :: NewProperty

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------


        call GetData(NewProperty%TimeSerie,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'TIME_SERIE',                                        &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR01'
        

        call GetData(NewProperty%BoxTimeSerie,                                           &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'BOX_TIME_SERIE',                                    &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR02'
        
        if (NewProperty%BoxTimeSerie) then
            Me%Output%Boxes_ON = .true.
            Me%NumberPropForBoxes = Me%NumberPropForBoxes + 1
        endif

        call GetData(NewProperty%BoxTimeSerie2D,                                           &
                     Me%ObjEnterData, iflag,                                               &
                     Keyword      = 'BOX_TIME_SERIE2D',                                    &
                     Default      = .false.,                                               &
                     SearchType   = FromBlock,                                             &
                     ClientModule = 'ModulePorousMediaProperties',                         &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR03'

        call GetData(NewProperty%OutputHDF,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'OUTPUT_HDF',                                        &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR04'
        
        !output the average conc in aquifer and in vadoze zone for each soil collumn
        !if (NewProperty%OutputHDF) then
            call GetData(NewProperty%OutputAverageConc,                                      &
                         Me%ObjEnterData, iflag,                                             &
                         Keyword      = 'OUTPUT_AVERAGE_CONC',                               &
                         ClientModule = 'ModulePorousMediaProperties',                       &
                         Default      = .false.,                                             &
                         SearchType   = FromBlock,                                           &
                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR05'
            
            if (NewProperty%OutputAverageConc) then
            
                Me%Output%AverageConc_ON = .true.
                
                allocate(NewProperty%AverageAquiferConc(Me%WorkSize%ILB:Me%WorkSize%IUB,     &
                                                        Me%WorkSize%JLB:Me%WorkSize%JUB),    &
                                                        STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR10'
                NewProperty%AverageAquiferConc(:,:) = FillValueReal            

                allocate(NewProperty%AverageVadozeConc(Me%WorkSize%ILB:Me%WorkSize%IUB,      &
                                                       Me%WorkSize%JLB:Me%WorkSize%JUB),     &
                                                       STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR20'
                NewProperty%AverageAquiferConc(:,:) = FillValueReal  
            endif
            
            !compute average decay
            if (NewProperty%Evolution%Decay) then
            
                call GetData(NewProperty%OutputAverageDecay,                                     &
                             Me%ObjEnterData, iflag,                                             &
                             Keyword      = 'OUTPUT_AVERAGE_DECAY',                              &
                             ClientModule = 'ModulePorousMediaProperties',                       &
                             Default      = .false.,                                             &
                             SearchType   = FromBlock,                                           &
                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR030'            
                
                if (NewProperty%OutputAverageDecay) then
                
                    Me%Output%AverageDecay_ON = .true.
                
                    allocate(NewProperty%AverageAquiferDecay(Me%WorkSize%ILB:Me%WorkSize%IUB,     &
                                                            Me%WorkSize%JLB:Me%WorkSize%JUB),    &
                                                            STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR40'
                    NewProperty%AverageAquiferDecay(:,:) = FillValueReal  
                    
                    allocate(NewProperty%AverageVadozeDecay(Me%WorkSize%ILB:Me%WorkSize%IUB,      &
                                                           Me%WorkSize%JLB:Me%WorkSize%JUB),     &
                                                           STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR50'
                    NewProperty%AverageAquiferDecay(:,:) = FillValueReal 
                
                endif
                
                !integrate decay (all column) over the simulation period
                call GetData(NewProperty%OutputIntegratedDecay,                                  &
                             Me%ObjEnterData, iflag,                                             &
                             Keyword      = 'OUTPUT_INTEGRATED_DECAY',                           &
                             ClientModule = 'ModulePorousMediaProperties',                       &
                             Default      = .false.,                                             &
                             SearchType   = FromBlock,                                           &
                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR060'            
                
                if (NewProperty%OutputIntegratedDecay) then
                
                    Me%Output%IntegratedDecay_ON = .true.
                
                    allocate(NewProperty%IntegratedDecay(Me%WorkSize%ILB:Me%WorkSize%IUB,        &
                                                            Me%WorkSize%JLB:Me%WorkSize%JUB),    &
                                                            STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR40'
                    NewProperty%IntegratedDecay(:,:) = 0.0  
                   
                    !Gets the root path from the file nomfich.dat
                    call ReadFileName("ROOT_SRT", NewProperty%IntegratedDecayFile, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunOff - ERR0750'
                    NewProperty%IntegratedDecayFile = trim(adjustl(NewProperty%IntegratedDecayFile))   &
                                                     //trim(NewProperty%ID%Name)//"_IntegratedDecay.dat"              
                
                endif                
                
            endif
            
        !endif
        
    end subroutine Construct_PropertyOutPut
   
   !---------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        integer                                             :: nProperties
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: n
        integer                                             :: TimeSerieNumber, dn, Id, Jd
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        character(len=StringLength)                         :: TimeSerieName
        
        !Begin------------------------------------------------------------------
        
        !Counts the number of Properties which has timeserie option set to true (2x for inf col concentration)
        PropertyX => Me%FirstProperty
        nProperties = 0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                if (PropertyX%Evolution%AdvectionDiffusion) then
                    nProperties = nProperties + 2
                else
                    nProperties = nProperties + 1
                endif
                if (PropertyX%Evolution%Decay) then
                    nProperties = nProperties + 1
                    if (PropertyX%OutputAverageDecay) then
                        nProperties = nProperties + 2
                    endif
                    if (PropertyX%OutputIntegratedDecay) then
                        nProperties = nProperties + 1
                    endif
                endif
                if (PropertyX%OutputAverageConc) then
                    nProperties = nProperties + 2
                endif
            endif
            PropertyX => PropertyX%Next
        enddo

        if (Me%CalculateECw) then
            allocate(PropertyList(nProperties + 1))
        else
            allocate(PropertyList(nProperties))
        endif

!        !Allocates PropertyList
!        allocate(PropertyList(nProperties))
        
        !Property names
        n=1
        PropertyX  => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                PropertyList(n)  = trim(PropertyX%ID%Name)//'['//trim(PropertyX%ID%Units)//']'
                n=n+1
                
                if (PropertyX%Evolution%AdvectionDiffusion) then
                    PropertyList(n)  = trim(PropertyX%ID%Name)//'_in_InfilColumn['//trim(PropertyX%ID%Units)//']'
                    n=n+1
!                    PropertyList(n)  = trim(PropertyX%ID%Name)//'_Diffusivity[m2/s]'
!                    n=n+1
!                    PropertyList(n)  = trim(PropertyX%ID%Name)//'_DiffusivitySurface[m2/s]'
!                    n=n+1                
                endif
                if (PropertyX%Evolution%Decay) then
                    PropertyList(n)  = trim(PropertyX%ID%Name)//'_DecayRate['//trim(PropertyX%ID%Units)//'.day-1]'
                    n=n+1
                    if (PropertyX%OutputAverageDecay) then
                        PropertyList(n)  = trim(PropertyX%ID%Name)//'_AvrgAquifDecay['//trim(PropertyX%ID%Units)//'.day-1]'
                        n=n+1  
                        PropertyList(n)  = trim(PropertyX%ID%Name)//'_AvrgVadozeDecay['//trim(PropertyX%ID%Units)//'.day-1]'
                        n=n+1                                          
                    endif
                    if (PropertyX%OutputIntegratedDecay) then
                        PropertyList(n)  = trim(PropertyX%ID%Name)//'_IntegratedDecay[kg/ha]'
                        n=n+1                                        
                    endif                    
                endif
                if (PropertyX%OutputAverageConc) then
                    PropertyList(n)  = trim(PropertyX%ID%Name)//'_AvrgAquiferConc['//trim(PropertyX%ID%Units)//']'
                    n=n+1 
                    PropertyList(n)  = trim(PropertyX%ID%Name)//'_AvrgVadozeConc['//trim(PropertyX%ID%Units)//']'
                    n=n+1                                         
                endif                
            endif
            PropertyX=>PropertyX%Next
        enddo

        if (Me%CalculateECw) then
            PropertyList (n) = 'ECw[dS/L]'
        endif

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMediaProperties',                      &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR010' 

        if (iflag == 1) then
            Me%OutPut%TimeSerie_ON = .true.
        else
            Me%OutPut%TimeSerie_ON = .false.
        endif
        
        !Get water points
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ContructTimeSerie - ModulePorousMediaProperties - ERR020'


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srpp",                                        &
                            WaterPoints3D    = Me%ExtVar%WaterPoints3D,                 &
                            UseTabulatedData = .false.,                                 &
                            STAT             = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR030' 

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR040'

        !Deallocates PropertyList
        deallocate(PropertyList)


        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR03'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR04'
            
            call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR04'
            
i1:         if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR05'

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR06'

                    if (IgnoreOK) then
                        write(*,*) 'Time Serie outside the domain - ',trim(TimeSerieName)
                        cycle
                    else
                        stop 'ConstructTimeSerie - PorousMedia - ERR07'
                    endif

                endif


                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR08'

            endif i1

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      LocalizationI   = Id,                             &
                                      LocalizationJ   = Jd,                             & 
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR09'

            call GetWaterPoints3D (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR090'

            if (Me%ExtVar%WaterPoints3D(Id, Jd, Me%WorkSize%KUB) /= WaterPoint) then
                 write(*,*) 'Time Serie in a land cell - ',trim(TimeSerieName)
            endif
            
            call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModulePorousMediaProperties - ERR0140'


        enddo

       
    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------
    
    subroutine ConstructHDF
        
        !External-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        logical                                     :: OutputON
        integer :: STAT_CALL

        !Begin-----------------------------------------------------------------

        nullify(Me%OutPut%OutTime)

        OutputON        = OFF

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))
            
            if(CurrentProperty%OutputHDF) OutputON = ON
            CurrentProperty => CurrentProperty%Next

        enddo

        if(OutputON)then
  
            call GetOutPutTime(Me%ObjEnterData,                              &
                               CurrentTime = Me%ExtVar%BeginTime,            &
                               EndTime     = Me%ExtVar%EndTime,              &
                               keyword     = 'OUTPUT_TIME',                  &
                               SearchType  = FromFile,                       &
                               OutPutsTime = Me%OutPut%OutTime,              &
                               OutPutsOn   = Me%OutPut%HDF_ON,               &
                               OutPutsNumber = Me%OutPut%Number,             &
                               STAT        = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                       &
                stop 'ConstructHDF - ModulePorousMediaProperties - ERR01' 

            if (Me%OutPut%HDF_ON) then

                Me%OutPut%NextOutPut = 1

                call Open_HDF5_OutPut_File

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructHDF - ModulePorousMediaProperties - ERR02'
            endif 

            !Output for restart
            call GetOutPutTime(Me%ObjEnterData,                                             &
                               CurrentTime  = Me%ExtVar%Now,                                &
                               EndTime      = Me%ExtVar%EndTime,                            &
                               keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                               SearchType   = FromFile,                                     &
                               OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                               OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR11a'

        endif

    end subroutine ConstructHDF

    !--------------------------------------------------------------------------

    subroutine StartOutputBoxFluxes

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
        integer                                             :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                                             :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        type(T_Property    ),                       pointer :: PropertyX
        type(T_SedimentRate),                       pointer :: SedimentRateX
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList
        character(len=StringLength), dimension(:),  pointer :: FluxesOutputList
        integer                                             :: nScalars, n, nFluxes

        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        ! This keyword have two functions if exist fluxes between boxes are compute 
        ! and the value read is the name file where the boxes are defined
        call GetData(Me%Files%BoxesFile,                                            &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'BOXFLUXES',                                    &
                     ClientModule = 'ModulePorousMediaProperties',                  &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR01'
        if (iflag .EQ. 0)                                                           &
            stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR02'    
        
        inquire(File = Me%Files%BoxesFile, Exist = exist)
        if (exist) then
            inquire(File = Me%Files%BoxesFile, Opened  = Opened)
            if (opened) then
                write(*,*    ) 
                write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Files%BoxesFile))
                write(*,*    ) 'Already opened.'
                stop           'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR03'    
            end if
        else
            write(*,*) 
            write(*,*)     'Could not find the boxes file.'
            write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
            stop           'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR04'    
        end if
        
        !Output Rates and Properties inside box (.bxm)
        nScalars = Me%NumberPropForBoxes + Me%SedimentRatesNumber
        !Output Properties fluxes between boxes
        nFluxes  = Me%NumberPropForBoxes
            
        allocate(ScalarOutputList(nScalars), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR05'

        allocate(FluxesOutputList(nFluxes), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR06'

        n = 0
        if (Me%Output%Boxes_ON) then
            PropertyX => Me%FirstProperty
            do while(associated(PropertyX))
                if (PropertyX%BoxTimeSerie) then
                    n = n + 1
                    ScalarOutputList(n) = "soil_"//trim(PropertyX%ID%Name)
                    FluxesOutputList(n) = "soil_"//trim(PropertyX%ID%Name)
                endif
                PropertyX => PropertyX%Next
            end do        
        endif        
        if (Me%Output%RateFluxes) then
            SedimentRateX => Me%FirstSedimentRate
            do while(associated(SedimentRateX))
                n = n + 1
                ScalarOutputList(n) ="soil_"//trim(SedimentRateX%ID%Name)
                SedimentRateX => SedimentRateX%Next
            end do
        endif

        
        call GetWaterPoints3D (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR090'


        call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                 &
                         TimeID             = Me%ObjTime,                   &
                         HorizontalGridID   = Me%ObjHorizontalGrid,         &
                         BoxesFilePath      = Me%Files%BoxesFile,           &
                         FluxesOutputList   = FluxesOutputList,             &
                         ScalarOutputList   = ScalarOutputList,             &
                         WaterPoints3D      = Me%ExtVar%WaterPoints3D,      &
                         Size3D             = Me%Size,                      &
                         WorkSize3D         = Me%WorkSize,                  &
                         STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR07'

        deallocate(ScalarOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR09'

        deallocate(FluxesOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR10'


        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR0140'
    
        if (Me%Output%Boxes_ON) then

            allocate(Me%CellMass(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'StartOutputBoxFluxes - ModulePorousMediaProperties - ERR170'
            Me%CellMass(:,:,:) = 0.            
        
        endif
        
        
    end subroutine StartOutputBoxFluxes

    !--------------------------------------------------------------------------

     subroutine Open_HDF5_OutPut_File        

        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB,KLB,KUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        !Begin-----------------------------------------------------------------

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
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR10'
     
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR020'

        !Sets limits for next write operations
        call HDF5SetLimits (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR030'
        
        call GetGridData (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR040'
        
        call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR050'  
        
        !Writes the Grid
        call HDF5WriteData (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                            Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR060'

        !WriteBasinPoints
        call HDF5WriteData (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                            Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR070'

        call HDF5SetLimits (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR080'

        call GetWaterPoints3D (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR090'
        
        !WaterPoints3D
        call HDF5WriteData (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",       &
                            Array3D = Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR100'

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR110'       

        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR120'  

        call UnGetGridData              (Me%ObjGridData, Me%ExtVar%Topography,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR130'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR0140'

    end subroutine Open_HDF5_OutPut_File   
   
    !--------------------------------------------------------------------------

   
    subroutine ReadOldConcBoundariesHDF(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        character (Len=StringLength)                :: PropertyName
        logical                                     :: EXIST
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ
        logical, save                               :: FileNameRead = .false.                                     
        !----------------------------------------------------------------------

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

        !----------------------------------------------------------------------

        !Reads name of the file from Nomfich
        if (.not. FileNameRead) then
        
            !Reads the name of the file where to store final data
            call ReadFileName ('POROUS_PROP_INI', Me%Files%InitialFile, "PorousMedia Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ReadOldConcBoundariesHDF - ERR00'

            FileNameRead = .true.

        endif


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
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR01'


            PropertyName = trim(adjustl(NewProperty%ID%name))

            NewProperty%Concentration(:,:,:) = FillValueReal


            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB,                     &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR02'

            call HDF5ReadData   (ObjHDF5, "/Results/"//trim(adjustl(NewProperty%ID%Name)),        &
                                 trim(adjustl(NewProperty%ID%Name)),                              &
                                 Array3D = NewProperty%Concentration,                    &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR03'


            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR06'

        else
            
            write(*,*)
            stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR07'

        end if cd0

    end subroutine ReadOldConcBoundariesHDF


    !--------------------------------------------------------------------------


    subroutine CheckFieldConsistence(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer               :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: Counter
        integer                                 :: i,j,k
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                                 :: StopSubroutine = .false.
        integer                                 :: UnitAux, kaux
        real                                    :: Aux, Sum
        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

            
        !Verification if the values read are lower than zero in water points
        do I = ILB, IUB
        do J = JLB, JUB
        do K = KLB, KUB
            
            if (Me%ExtVar%WaterPoints3D(i, j, k) == WaterPoint) then
                               
                if (NewProperty%Concentration(i, j, k) < 0.) then

                    StopSubroutine = .true.
                    Aux            = -1
                    kaux           = k + 1

                    do while (Aux < 0)

                        Aux = NewProperty%Concentration(i, j, kaux) 

                        kaux = kaux + 1

                        if (kaux > KUB)  then 

                            Counter = 0
                            Sum     = 0
                            
                            if (NewProperty%Concentration(i-1, j, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i-1, j, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i-1, j, k)
                            
                            endif

                        
                            if (NewProperty%Concentration(i+1, j, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i+1, j, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i+1, j, k)
                            
                            endif

                            if (NewProperty%Concentration(i, j-1, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i, j-1, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i, j-1, k)
                            
                            endif

                            if (NewProperty%Concentration(i, j+1, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i, j+1, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i, j+1, k)
                            
                            endif
                                  

                            if (Counter > 0) then                                        

                                Aux = Sum / real(Counter)
                                exit

                            else

                                stop 'Subroutine CheckFieldConsistence; PorousMediaProperties. ERR01.'

                            endif

                        endif

                    enddo

                    NewProperty%Concentration(i, j, k) = Aux

                endif

            else

                NewProperty%Concentration(i, j, k) = FillValueReal

            endif

        enddo
        enddo
        enddo

        if (StopSubroutine) then                                                   
            
            call UnitsManager(UnitAux, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - PorousMediaProperties - ERR02' 

            open(UnitAux, FILE = trim(NewProperty%ID%name)//'.new',                 &
                 FORM = 'FORMATTED', STATUS = 'UNKNOWN', IOSTAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - PorousMediaProperties - ERR03' 

            write(UnitAux,*) '<ConcentrationBegin>'
           
            do I = ILB, IUB
            do J = JLB, JUB
            do K = KLB, KUB
            
                write(UnitAux,*) NewProperty%Concentration(i, j, k)

            enddo
            enddo
            enddo

            write(UnitAux,*) '<ConcentrationEnd>'

            call UnitsManager(UnitAux, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - PorousMediaProperties - ERR04' 

            write(*,*) 'A new concentration file was created for property: ', trim(NewProperty%ID%Name)
            write(*,*) 'Run again with this new file ', trim(NewProperty%ID%name)//'.new'
            stop 'CheckFieldConsistence - PorousMediaProperties - ERR05'  

        endif

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------


    subroutine CoupleSoilQuality        

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX, Oxygen
        integer, pointer, dimension(:)                      :: SoilQualityPropertyList
        integer                                             :: STAT_CALL
        real                                                :: SoilQualityDT
        integer                                             :: nProp = 0 

        !Begin------------------------------------------------------------------

        !Counts the number of Properties which has WaterQuality option set to true
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilQuality) then
                nProp = nProp + 1
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates Array to hold IDs
        allocate (SoilQualityPropertyList(1:nProp))

        !Fills Array
        PropertyX => Me%FirstProperty
        nProp = 0
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilQuality) then
                nProp = nProp + 1
                SoilQualityPropertyList(nProp) = PropertyX%ID%IDNumber
            endif
            PropertyX => PropertyX%Next
        enddo

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModulePorousMediaProperties - ERR01'

        !Start Interface
        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = SedimentQualityModel,          &
                                DT                  = SoilQualityDT,                 &
                                PropertiesList      = SoilQualityPropertyList,       &
                                WaterPoints3D       = Me%ExtVar%WaterPoints3D,       &
                                Size3D              = Me%WorkSize,                   &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleSoilQuality - ModulePorousMediaProperties - ERR02'


        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModulePorousMediaProperties - ERR03'


        deallocate (SoilQualityPropertyList)

        Me%Coupled%SoilQuality_DT          = SoilQualityDT 
        Me%Coupled%SoilQuality_NextCompute = Me%ExtVar%Now    

        nullify (Me%DissolvedToParticulate3D)
        allocate(Me%DissolvedToParticulate3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, &
                 Me%Size%KLB:Me%Size%KUB))
        Me%DissolvedToParticulate3D(:,:,:) = null_real

        Me%ResidualTime = 0.
        
        call SearchProperty(Oxygen, Oxygen_ , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR45'                    
        if (Oxygen%ID%SolutionFromFile)   then
            Me%Coupled%SedQualityOxygenForcing = .true.
        endif     
        
        
    end subroutine CoupleSoilQuality

    !--------------------------------------------------------------------------


#ifdef _PHREEQC_
    !--------------------------------------------------------------------------
    subroutine CoupleSoilChemistry        

        !Local-----------------------------------------------------------------
        type(T_Property), pointer       :: PropertyX
        integer, pointer, dimension(:)  :: SoilChemistryPropertyList
        integer                         :: STAT_CALL
        real                            :: SoilChemistryDT
        integer                         :: nProp = 0 
        type(T_PhreeqCModel), pointer   :: phreeqc
        integer                         :: index
        integer                         :: i, j, k

        !Begin------------------------------------------------------------------
        !Counts the number of Properties which has SoilChemistry option set to true
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilChemistry) then
                nProp = nProp + 1
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates Array to hold IDs
        allocate (SoilChemistryPropertyList(1:nProp))

        !Fills Array
        PropertyX => Me%FirstProperty
        nProp = 0
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilChemistry) then
                nProp = nProp + 1
                SoilChemistryPropertyList(nProp) = PropertyX%ID%IDNumber
            endif
            PropertyX => PropertyX%Next
        enddo       
        
        do index = 1, Me%PhreeqC%NumberOfModels       
            
            phreeqc => Me%PhreeqC%Models(index) 
               
            call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR01'

            !Start Interface
            call ConstructInterface(InterfaceID         = phreeqc%ObjInterface,          &
                                    TimeID              = Me%ObjTime,                    &
                                    SinksSourcesModel   = PhreeqCModel,                  &
                                    DT                  = SoilChemistryDT,               &
                                    PropertiesList      = SoilChemistryPropertyList,     &
                                    WaterPoints3D       = Me%ExtVar%WaterPoints3D,       &
                                    PhreeqCDatabase     = phreeqc%Database,              &
                                    PhreeqCDatabaseAux  = phreeqc%DatabaseAux,           &
                                    PhreeqCModelID      = phreeqc%ID,                    &
                                    Size3D              = Me%WorkSize,                   &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR02'

            call SetupSoilChemistry  

            do k = phreeqc%KLB, phreeqc%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%PhreeqC%Filter(i, j, k) = phreeqc%ID
            enddo
            enddo
            enddo   
            
            call UnGetMap (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR03'
            
        enddo
        
        deallocate (SoilChemistryPropertyList)

        Me%Coupled%SoilChemistry_DT          = SoilChemistryDT 
        Me%Coupled%SoilChemistry_NextCompute = Me%ExtVar%Now    
            
    end subroutine CoupleSoilChemistry
    !--------------------------------------------------------------------------
    
    subroutine SetupSoilChemistry   

        !Local-----------------------------------------------------------------
        type(T_Property), pointer :: PropertyX
        logical                   :: found
    
        !Begin-----------------------------------------------------------------           
        found = .false.
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX) .AND. (.NOT. found))
            if (PropertyX%ID%IDNumber .EQ. SoilDryDensity_) then
                Me%PhreeqC%SoilDryDensity => PropertyX  
                
                !Done to ensure that this property will not be sended to PhreeqC
                Me%PhreeqC%SoilDryDensity%Evolution%SoilChemistry = .false.
                
                found = .true.
                    
            endif

            PropertyX => PropertyX%Next
        enddo

        if (.not. found) then !stop 'PhreeqC needs property soil dry density - PorousMediaProperties - SetupSoilChemistry - ERR04'
        
            Me%PhreeqC%SoilDryDensity => null()            
        
        end if
    
    end subroutine SetupSoilChemistry
   
    !--------------------------------------------------------------------------
    
#endif

    !--------------------------------------------------------------------------
    subroutine CoupleChainReactions
        
        !Local-----------------------------------------------------------------
        integer                   :: STAT
        integer                   :: PropertiesCount
        type(T_Property), pointer :: SoilDensity

        !Begin-----------------------------------------------------------------
        call StartChainReactions (Me%ObjChainReactions,             &   
                                  Me%PropertiesList,                &
                                  Me%ObjTime,                       &
                                  Me%ObjHorizontalGrid,             &
                                  Me%ObjHorizontalMap,              &
                                  Me%ObjBasinGeometry,              &
                                  Me%ObjGeometry,                   &
                                  Me%ObjMap,                        &
                                  2,                                & 
                                  "PorousMediaProperties",          &
                                  Me%Files%ChainReactionsDataFile,  &
                                  PropertiesCount,                  &
                                  STAT = STAT)
        if (STAT /= SUCCESS_) &
            stop 'PorousMediaProperties - CoupleChainReactions - ERR010'
            
        if (PropertiesCount .EQ. 0) then
            
            call KillChainReactions(Me%ObjChainReactions, STAT = STAT)
            if (STAT /= SUCCESS_) &
                stop 'PorousMediaProperties - CoupleChainReactions - ERR020'
                
            Me%Coupled%ChainReactions = .false.
            
        else
        
            call SetChainReactionsProperties
            
            call Search_Property(SoilDensity, SoilDryDensity_, STAT) 
            if (STAT == SUCCESS_) then            
                call InitCRSoilPhase (Me%ObjChainReactions, SoilDensity3D = SoilDensity%Concentration, STAT = STAT)
                if (STAT .NE. SUCCESS_) &
                    stop 'PorousMediaProperties - CoupleChainReactions - ERR030' 
            endif     
                               
        endif        
        !----------------------------------------------------------------------
        
    end subroutine CoupleChainReactions
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine Search_Property(PropertyX, PropertyXID, STAT)

        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer             :: PropertyX
        integer,                    intent (IN)         :: PropertyXID
        integer         , optional, intent (OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber==PropertyXID) then
                exit        
            else
                PropertyX => PropertyX%Next                 
            end if    
        end do    

       if (associated(PropertyX)) then

            STAT_ = SUCCESS_  

        else
            STAT_ = NOT_FOUND_ERR_  
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine Search_Property

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
!
    !---------------------------------------------------------------------------

    subroutine GetPMPCoupled(PorousMediaPropertiesID, &
                             SoilQuality,             &
#ifdef _PHREEQC_
                             SoilChemistry,           &
#endif                             
                             STAT) 

        !Arguments-------------------------------------------------------------
        integer                        :: PorousMediaPropertiesID
        integer, optional, intent(OUT) :: STAT
        logical, optional, intent(OUT) :: SoilQuality        
#ifdef _PHREEQC_
        logical, optional, intent(OUT) :: SoilChemistry
#endif

        !External--------------------------------------------------------------

        integer :: ready_              

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(SoilQuality   )) SoilQuality    = Me%Coupled%SoilQuality
            
#ifdef _PHREEQC_            
            if (present(SoilChemistry )) SoilChemistry  = Me%Coupled%SoilChemistry
#endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetPMPCoupled
    
    !--------------------------------------------------------------------------
    
    subroutine GetECw (PorousMediaPropertiesID, ECw, STAT) 

        !Arguments-------------------------------------------------------------
        integer                         :: PorousMediaPropertiesID
        real, dimension(:,:,:), pointer :: ECw
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer :: STAT_              !Auxiliar local variable
        
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPOROUSMEDIAPROPERTIES_, Me%InstanceID)

            if (associated(Me%ECw)) then
                ECw   => Me%ECw
                STAT_ =  SUCCESS_
            else
                STAT_ =  NOT_ASSOCIATE_
            endif
            
        else
         
            STAT_ = ready_
            
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetECw

    !--------------------------------------------------------------------------

    subroutine GetPMPConcentration(PorousMediaPropertiesID, ConcentrationX, PropertyXIDNumber, &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        real, pointer, dimension(:,:,:)             :: ConcentrationX
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

        call Ready(PorousMediaPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mPOROUSMEDIAPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                ConcentrationX => PropertyX%concentration

                if (present(PropertyXUnits)) then 
                   UnitsSize      = LEN (PropertyXUnits)
                   PropertyXUnits = PropertyX%ID%Units(1:UnitsSize)
                end if

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetPMPConcentration

    !--------------------------------------------------------------------------------


    subroutine GetPMPConcentrationOld(PorousMediaPropertiesID, ConcentrationXOld, PropertyXIDNumber, &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        real, pointer, dimension(:,:,:)             :: ConcentrationXOld
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

        call Ready(PorousMediaPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mPOROUSMEDIAPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                ConcentrationXOld => PropertyX%ConcentrationOld

                if (present(PropertyXUnits)) then 
                   UnitsSize      = LEN (PropertyXUnits)
                   PropertyXUnits = PropertyX%ID%Units(1:UnitsSize)
                end if

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetPMPConcentrationOld

    !--------------------------------------------------------------------------------

    subroutine GetPMPTotalStoredMass(PorousMediaPropertiesID, TotalStoredMass, PropertyXIDNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        real(8)                                     :: TotalStoredMass
        integer,                      intent(IN )   :: PropertyXIDNumber
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_CALL              
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: STAT_    

        !------------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mPOROUSMEDIAPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                
                TotalStoredMass = PropertyX%TotalStoredMass

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetPMPTotalStoredMass

    !--------------------------------------------------------------------------------

    subroutine SetWindVelocity (PorousMediaPropertiesID, WindModulus, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        real, dimension(:,:), pointer               :: WindModulus

        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            

            Me%ExtVar%WindVelocity2D   => WindModulus


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_

    end subroutine SetWindVelocity 

    !---------------------------------------------------------------------------

 
    !--------------------------------------------------------------------------    
    !non matrixes values (modelling options) should pass to PMP in construction phase
    subroutine SetVegetationPMProperties(PorousMediaPropertiesID,                  &
                                         SoilFluxesActive,                         &
                                         GrazingBiomass,                           &
                                         GrazingNitrogen,                          &
                                         GrazingPhosphorus,                        &
                                         HarvestKillAerialBiomass,                 &
                                         HarvestKillNitrogen,                      &
                                         HarvestKillPhosphorus,                    &
                                         HarvestKillRootBiomass,                   &
                                         DormancyBiomass,                          &
                                         DormancyNitrogen,                         &
                                         DormancyPhosphorus,                       &
                                         FertilNitrateSurface,                     &
                                         FertilNitrateSubSurface,                  &
                                         FertilAmmoniaSurface,                     &
                                         FertilAmmoniaSubSurface,                  &
                                         FertilOrganicNSurface,                    &
                                         FertilOrganicNSubSurface,                 &
                                         FertilOrganicPSurface,                    &
                                         FertilOrganicPSubSurface,                 &
                                         FertilMineralPSurface,                    &
                                         FertilMineralPSubSurface,                 &
                                         NitrogenUptake,                           &
                                         PhosphorusUptake,                         &
                                         Grazing,                                  &
                                         HarvestKill,                              &
                                         Dormancy,                                 &
                                         Fertilization,                            &
                                         NutrientFluxesWithSoil,                   &
                                         RootDepth,                                &
                                         ModelWater,                               &
                                         ModelNitrogen,                            &
                                         ModelPhosphorus,                          &
                                         GrowthModel,                              &
                                         CoupledVegetation,                        &
                                         NitrogenFraction,                         &
                                         PhosphorusFraction,                       &
                                         VegetationDT,                             &
                                         TranspirationBottomLayer,                 &
                                         NutrientUptakeMethod,                     &
                                         PesticideIDNumber,                        &
                                         PesticideSoilFlux,                        &
                                         STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        logical, dimension(:,:), pointer, optional  :: SoilFluxesActive
        real, dimension(:,:), pointer, optional     :: GrazingBiomass
        real, dimension(:,:), pointer, optional     :: GrazingNitrogen
        real, dimension(:,:), pointer, optional     :: GrazingPhosphorus
        real, dimension(:,:), pointer, optional     :: HarvestKillAerialBiomass
        real, dimension(:,:), pointer, optional     :: HarvestKillNitrogen
        real, dimension(:,:), pointer, optional     :: HarvestKillPhosphorus
        real, dimension(:,:), pointer, optional     :: HarvestKillRootBiomass
        real, dimension(:,:), pointer, optional     :: DormancyBiomass
        real, dimension(:,:), pointer, optional     :: DormancyNitrogen
        real, dimension(:,:), pointer, optional     :: DormancyPhosphorus
        real, dimension(:,:), pointer, optional     :: FertilNitrateSurface
        real, dimension(:,:), pointer, optional     :: FertilNitrateSubSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSubSurface
        real, dimension(:,:), pointer, optional     :: FertilMineralPSurface
        real, dimension(:,:), pointer, optional     :: FertilMineralPSubSurface
        real, dimension(:,:), pointer, optional     :: NitrogenFraction
        real, dimension(:,:), pointer, optional     :: PhosphorusFraction
        real, dimension(:,:), pointer, optional     :: RootDepth
        real, dimension(:,:,:),pointer,optional     :: NitrogenUptake
        real, dimension(:,:,:),pointer,optional     :: PhosphorusUptake
        integer, dimension(:,:), pointer, optional  :: TranspirationBottomLayer
        logical,  optional                          :: Grazing
        logical,  optional                          :: HarvestKill
        logical,  optional                          :: Dormancy
        logical,  optional                          :: Fertilization
        logical,  optional                          :: NutrientFluxesWithSoil
        logical,  optional                          :: ModelWater
        logical,  optional                          :: ModelNitrogen
        logical,  optional                          :: ModelPhosphorus
        logical,  optional                          :: GrowthModel
        logical,  optional                          :: CoupledVegetation
        real,     optional                          :: VegetationDT
        integer,  optional                          :: NutrientUptakeMethod
        integer, optional                           :: PesticideIDNumber
        real, dimension(:,:), pointer, optional     :: PesticideSoilFlux
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        type (T_Property), pointer                  :: PropertyX
     
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
            if (present(TranspirationBottomLayer)) then
                Me%ExtVar%TranspirationBottomLayer => TranspirationBottomLayer
            endif
            
            if (present(NutrientFluxesWithSoil  )) then
                Me%ExtVar%ComputeVegInterfaceFluxes   = NutrientFluxesWithSoil
            endif

            if (present(SoilFluxesActive        )) then
                Me%ExtVar%SoilFluxesActive         => SoilFluxesActive
            endif
            if (present(Grazing)) then
                Me%ExtVar%Grazing                  = Grazing
                if (present(GrazingBiomass          )) Me%ExtVar%GrazingBiomass           => GrazingBiomass
                if (present(GrazingNitrogen         )) Me%ExtVar%GrazingNitrogen          => GrazingNitrogen
                if (present(GrazingPhosphorus       )) Me%ExtVar%GrazingPhosphorus        => GrazingPhosphorus
            endif
            if (present(HarvestKill)) then
                Me%ExtVar%HarvestKill               = HarvestKill
                if (present(HarvestKillAerialBiomass )) Me%ExtVar%HarvestKillAerialBiomass  => HarvestKillAerialBiomass
                if (present(HarvestKillNitrogen      )) Me%ExtVar%HarvestKillNitrogen       => HarvestKillNitrogen
                if (present(HarvestKillPhosphorus    )) Me%ExtVar%HarvestKillPhosphorus     => HarvestKillPhosphorus
                if (present(HarvestKillRootBiomass   )) Me%ExtVar%HarvestKillRootBiomass    => HarvestKillRootBiomass
            endif
            if (present(Dormancy)) then
                Me%ExtVar%Dormancy                 = Dormancy
                if (present(DormancyBiomass         )) Me%ExtVar%DormancyBiomass          => DormancyBiomass
                if (present(DormancyNitrogen        )) Me%ExtVar%DormancyNitrogen         => DormancyNitrogen
                if (present(DormancyPhosphorus      )) Me%ExtVar%DormancyPhosphorus       => DormancyPhosphorus
            endif
            if (present(Fertilization)) then
                Me%ExtVar%Fertilization            = Fertilization
                if (present(FertilNitrateSurface    )) Me%ExtVar%FertilNitrateSurface     => FertilNitrateSurface
                if (present(FertilNitrateSubSurface )) Me%ExtVar%FertilNitrateSubSurface  => FertilNitrateSubSurface
                if (present(FertilAmmoniaSurface    )) Me%ExtVar%FertilAmmoniaSurface     => FertilAmmoniaSurface
                if (present(FertilAmmoniaSubSurface )) Me%ExtVar%FertilAmmoniaSubSurface  => FertilAmmoniaSubSurface
                if (present(FertilOrganicNSurface   )) Me%ExtVar%FertilOrganicNSurface    => FertilOrganicNSurface
                if (present(FertilOrganicNSubSurface)) Me%ExtVar%FertilOrganicNSubSurface => FertilOrganicNSubSurface
                if (present(FertilOrganicPSurface   )) Me%ExtVar%FertilOrganicPSurface    => FertilOrganicPSurface
                if (present(FertilOrganicPSubSurface)) Me%ExtVar%FertilOrganicPSubSurface => FertilOrganicPSubSurface
                if (present(FertilMineralPSurface   )) Me%ExtVar%FertilMineralPSurface    => FertilMineralPSurface
                if (present(FertilMineralPSubSurface)) Me%ExtVar%FertilMineralPSubSurface => FertilMineralPSubSurface
            endif

            if (present(NitrogenUptake          )) Me%ExtVar%NitrogenUptake           => NitrogenUptake
            if (present(PhosphorusUptake        )) Me%ExtVar%PhosphorusUptake         => PhosphorusUptake
            if (present(RootDepth               )) Me%ExtVar%RootDepth                => RootDepth
            
            if (present(ModelWater              )) Me%ExtVar%ModelWater               =  ModelWater
            if (present(ModelNitrogen           )) Me%ExtVar%ModelNitrogen            =  ModelNitrogen
            if (present(ModelPhosphorus         )) Me%ExtVar%ModelPhosphorus          =  ModelPhosphorus

            if (present(GrowthModel             )) Me%ExtVar%GrowthModel              =  GrowthModel
            if (present(CoupledVegetation       )) Me%ExtVar%CoupledVegetation        =  CoupledVegetation

            if (present(VegetationDT            )) Me%ExtVar%VegetationDT             =  VegetationDT

            if (present(NitrogenFraction        )) Me%ExtVar%NitrogenFraction         => NitrogenFraction
            if (present(PhosphorusFraction      )) Me%ExtVar%PhosphorusFraction       => PhosphorusFraction
            
            if (present(NutrientUptakeMethod    )) Me%ExtVar%NutrientUptakeMethod     =  NutrientUptakeMethod

            if (present(PesticideIDNumber       )) then
                call SearchProperty(PropertyX, PropertyXIDNumber = PesticideIDNumber, STAT = STAT_)
                if (STAT_ == SUCCESS_) then
                    if (.not. present(PesticideSoilFlux)) stop 'SetVegetationPMProperties - ModulePorousMediaProperties - ERR10'
                    if (.not. associated(PropertyX%PesticideFlux)) then
                        write(*,*) 'Property', GetPropertyName(PesticideIDNumber) 
                        write(*,*) 'is a pesticide and needs PESTICIDE keyword in Porous Media Properties'
                        stop 'SetVegetationPMProperties - PorousMediaProperties - ERR01'
                    else    
                        PropertyX%PesticideFlux => PesticideSoilFlux
                    endif
                else
                    write(*,*)
                    write(*,*) 'Not found pesticide', GetPropertyName(PesticideIDNumber) 
                    write(*,*) 'in Porous Media Properties property list'
                    stop ' SetVegetationPMProperties - PorousMediaProperties - ERR10'
                endif    
            endif
            
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_

    end subroutine SetVegetationPMProperties 

    !---------------------------------------------------------------------------
    
    subroutine SetDNConcPMP (PorousMediaPropertiesID, PropertyID, DNConcentration, ChannelsID, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: PropertyID   
        real, dimension (:), pointer                    :: DNConcentration
        integer, dimension(:, :), pointer               :: ChannelsID
        integer                                         :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j
        type(T_Property), pointer                       :: PropertyX
   
        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
        
            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetDNConcPMP - ModulePorousMediaProperties - ERR01'
            
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then
            
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%RiverPoints(i,j) == WaterPoint) then
                        PropertyX%ConcentrationDN (i,j) = DNConcentration (ChannelsID(i,j))
                    endif
                enddo
                enddo                

            else
                write(*,*) 'Looking for Drainage Network Property in Porous Media Properties', GetPropertyName(PropertyID)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetDNConcPMP - ModulePorousMediaProperties - ERR010'
            end if

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetDNConcPMP - ModulePorousMediaProperties - ERR040'               


        else
            STAT_ = ready_
        end if

        STAT = STAT_        
                     
    end subroutine SetDNConcPMP

    !---------------------------------------------------------------------------
    

    subroutine SetInfColConcPMP (PorousMediaPropertiesID, PropertyXIDNumber, ConcentrationX, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: PropertyXIDNumber   
        real, dimension (:,:), pointer                  :: ConcentrationX
        integer                                         :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetInfColConcPMP - ModulePorousMediaProperties - ERR01'
        
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_)
            if (STAT_ == SUCCESS_) then
            
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
                        PropertyX%ConcentrationOnInfColumn (i,j) = ConcentrationX (i,j)
                    endif
                enddo
                enddo                

            else
                write(*,*) 'Looking for Atmosphere Property in Porous Media Properties', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetInfColConcPMP - ModulePorousMediaProperties - ERR010'
            end if
           
            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetInfColConcPMP - ModulePorousMediaProperties - ERR040'               


        else
            STAT_ = ready_
        end if

        STAT = STAT_        
                     
    end subroutine SetInfColConcPMP
    
    !---------------------------------------------------------------------------

    subroutine GetPMPMassBalance (PorousMediaPropertiesID,                  &
                                  PropertyID,                               &
                                  TotalStoredMass,                          &
                                  TranspiredMass,                           &
                                  DNExchangeMass,                           &
                                  RPExchangeMass,                           &
                                  STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: PropertyID
        real(8), optional                               :: TotalStoredMass, TranspiredMass
        real(8), optional                               :: DNExchangeMass, RPExchangeMass
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type (T_Property), pointer                      :: PropertyX
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
           
            nullify(PropertyX) 
           
            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then 
                if (present (TotalStoredMass)) TotalStoredMass = PropertyX%MB%TotalStoredMass
                if (present (TranspiredMass))  TranspiredMass  = PropertyX%MB%TranspiredMass
                if (present (DNExchangeMass))  DNExchangeMass  = PropertyX%MB%DNExchangeMass
                if (present (RPExchangeMass))  RPExchangeMass  = PropertyX%MB%RPExchangeMass
                
                STAT_CALL = SUCCESS_
            else
                stop 'GetPMPMassBalance - ModulePorousMediaProperties - ERR01'
            endif
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetPMPMassBalance

    !---------------------------------------------------------------------------

    subroutine CheckPMPProperty (PorousMediaPropertiesID,                    &
                                  PropertyID,                               &
                                  STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: PropertyID
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type (T_Property), pointer                      :: PropertyX
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            nullify(PropertyX) 

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then 
                if (.not. PropertyX%Evolution%AdvectionDiffusion) then
                    write(*,*)
                    write(*,*)'Property', GetPropertyName(PropertyID)
                    write(*,*)'has advection diffusion inactive in PorousMediaProperties Module'
                    write(*,*)'and it is unconsistent with activation in other Modules'
                    stop 'CheckPMPProperty - ModulePorousMediaProperties - ERR01' 
                else               
                    STAT_CALL = SUCCESS_
                endif
            else
                write(*,*)
                write(*,*)'Could not find property', GetPropertyName(PropertyID)
                write(*,*)'in PorousMediaProperties Module'
                stop 'CheckPMPProperty - ModulePorousMediaProperties - ERR010'
            endif
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine CheckPMPProperty

    !---------------------------------------------------------------------------

    subroutine GetPMPnProperties (PorousMediaPropertiesID, nProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: nProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            nProperties       = Me%PropertiesNumber
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetPMPnProperties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetPMPPropertiesIDByIdx (PorousMediaPropertiesID, Idx, ID, PropAdvDiff, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer, intent(IN)                             :: Idx
        integer, intent(OUT)                            :: ID
        logical, intent(OUT)                            :: PropAdvDiff
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, i
        type (T_Property), pointer                      :: CurrProp

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrProp => Me%FirstProperty
            do i = 1, idx - 1
                CurrProp => CurrProp%Next
            enddo

            ID          = CurrProp%ID%IDNumber
            PropAdvDiff = CurrProp%Evolution%AdvectionDiffusion
            
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetPMPPropertiesIDByIdx

    !---------------------------------------------------------------------------

    subroutine UnGetPorousMediaProperties3D_I(ObjPorousMediaPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID, "UnGetPorousMediaProperties3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMediaProperties3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetPorousMediaProperties3D_R4(ObjPorousMediaPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(4), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID,  "UnGetPorousMediaProperties3D_R4")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMediaProperties3D_R4


    !--------------------------------------------------------------------------

    subroutine UnGetPorousMediaProperties3D_R8(ObjPorousMediaPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID,  "UnGetPorousMediaProperties3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMediaProperties3D_R8


    !--------------------------------------------------------------------------

!    subroutine UnGetPorousMediaProperties3D_R8i(ObjPorousMediaPropertiesID, Array, STAT)
!
!        !Arguments-------------------------------------------------------------
!        integer                                         :: ObjPorousMediaPropertiesID
!        real(8), dimension(:, :), pointer               :: Array
!        integer, intent(OUT), optional                  :: STAT
!
!        !Local-----------------------------------------------------------------
!        integer                                         :: STAT_, ready_
!
!        !----------------------------------------------------------------------
!
!        STAT_ = UNKNOWN_
!
!        call Ready(ObjPorousMediaPropertiesID, ready_)
!
!        if (ready_ .EQ. READ_LOCK_ERR_) then
!
!            nullify(Array)
!            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID,  "UnGetPorousMediaProperties3D_R8")
!
!
!            STAT_ = SUCCESS_
!        else               
!            STAT_ = ready_
!        end if
!
!        if (present(STAT)) STAT = STAT_
!
!    end subroutine UnGetPorousMediaProperties3D_R8i


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyPorousMediaProperties(ObjPorousMediaPropertiesID,          &
                                           WaterColumn, WaterColumnOld,         &
                                           STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaPropertiesID
        real(8), dimension(:,:), pointer            :: WaterColumn
        real(8), dimension(:,:), pointer            :: WaterColumnOld
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_,STAT_CALL
        type(T_Property), pointer                   :: PropertyX
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyPorousMediaProperties")


        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then


            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMediaProperties - ModulePorousMediaProperties - ERR02'
            
            !Actualize the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMediaProperties - ModulePorousMediaProperties - ERR03'
                      
            call ReadLockExternalVar
            
            !New Water column to see if there will be diffusion with runoff
            Me%ExtVar%WaterColumn => WaterColumn
            !Old Water column - only used for pesticide add before transport
            Me%ExtVar%WaterColumnOLd => WaterColumnOLd

            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))
                
                call SetMatrixValue (PropertyX%ConcentrationOld, Me%Size, PropertyX%Concentration, Me%ExtVar%WaterPoints3D)
                
                if (Me%CheckGlobalMass) then
                    PropertyX%MB%TotalStoredMass     = 0.0
                    PropertyX%MB%TranspiredMass      = 0.0
                    PropertyX%MB%DNExchangeMass      = 0.0
                    PropertyX%MB%RPExchangeMass      = 0.0
                    PropertyX%MB%TotalDischargeMass  = 0.0
                endif
                
                PropertyX => PropertyX%Next
                
            end do                        
            
            !Actualize properties if evolution from file
            call ActualizePropertiesFromFile
            
            !Nutrient sources from vegetation. Sinks are explicit or implict cared in transport
            call InterfaceFluxes

            !Water discharges fluxes - explicit because may exist more than one discharge in the same cell
            if (Me%Coupled%Discharges) call ModifyDischarges
            
            if (Me%Coupled%AdvectionDiffusion) then
                call AdvectionDiffusionProcesses
            endif

            if (Me%Coupled%ChainReactions) then
                call ChainReactionsProcesses
            endif

            if (Me%Coupled%SoilQuality) then
                call SoilQualityProcesses
            endif
            
#ifdef _PHREEQC_
            if (Me%Coupled%SoilChemistry) then
                call SoilChemistryProcesses
            endif
#endif            

            if (Me%Coupled%Decay) then
                call DecayProcesses
            endif

            if (Me%Coupled%Partition)then
                call Partition_Processes 
            endif           

            if (Me%CalculateECw) then
                call ComputeECw
            endif
            
            !Compute average conc for each soil column
            if (Me%Output%AverageConc_ON) then
                call AverageConcentration
            endif
            !Compute average decay for each soil column
            if (Me%Coupled%Decay .and. (Me%Output%AverageDecay_ON .or. Me%Output%IntegratedDecay_ON)) then
                call DecayProcessing
            endif            

            if (Me%Output%Timeserie_ON) then
                call OutPut_TimeSeries
            endif

            if (Me%Output%HDF_ON) then
                call OutPut_HDF
            endif
            
            if (Me%Output%Boxes_ON) then
                call Output_Boxes_Mass
            endif

    !       call ProfileOutput    em teste no construct e no kill

            !Restart Output
            if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%Now == Me%ExtVar%EndTime)) then
                if(Me%ExtVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    call WriteFinalFile
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                endif
            endif
            
            if (Me%DTIntervalAssociated) then
                call Actualize_Time_Evolution
            endif
            
            if (Me%CheckGlobalMass) then
                call CalculateTotalStoredMass
            endif
        
            call ReadUnlockExternalVar

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyPorousMediaProperties")

    end subroutine ModifyPorousMediaProperties

    !-----------------------------------------------------------------------------

    subroutine AverageConcentration

        !Arguments----------------------------------------------------------------
        !Local--------------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX
        !Begin--------------------------------------------------------------------
        
        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))
            
            if (PropertyX%OutputAverageConc) then
            
                call ComputeAverageConc (PropertyX)
            
            endif
            
            PropertyX => PropertyX%Next
            
        end do                        

    
    end subroutine AverageConcentration

    !-----------------------------------------------------------------------------

    subroutine DecayProcessing

        !Arguments----------------------------------------------------------------
        !Local--------------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX
        !Begin--------------------------------------------------------------------
        
        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))
            
            if (PropertyX%OutputAverageDecay) then
            
                call ComputeAverageDecay (PropertyX)
            
            endif

            if (PropertyX%OutputIntegratedDecay) then
            
                call ComputeIntegratedDecay (PropertyX)
            
            endif
            
            PropertyX => PropertyX%Next
            
        end do                        

    
    end subroutine DecayProcessing

    !-----------------------------------------------------------------------------

    subroutine ModifyDischarges ()

        !Arguments--------------------------------------------------------------
        !Local------------------------------------------------------------------
        integer                                 :: iDis, nDischarges
        integer                                 :: i, j, k, n, kd
        integer                                 :: STAT_CALL
        real, dimension(:, :, :), pointer       :: FlowDischarge
        real, dimension(:, :   ), pointer       :: DischargesConc
        type (T_Property), pointer              :: Property
        integer                                 :: iProp
        real(8)                                 :: VolumeOld
        integer, dimension(:),   pointer        :: VectorI, VectorJ, VectorK
        integer                                 :: nCells, DischVertical
        integer                                 :: kmin, kmax, FlowDistribution

        !Get integrated flow from porous media to be sure using same values
        !ATTENTION - DIFFERENT DISCHARGES CAN NOT OCCUR AT THE SAME CELL
        !Test this here in PMP because is only problem if properties are used
        !or create duplications of integrated discharge (per discharge) in porous media
        call GetFlowDischarge (Me%ObjPorousMedia, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ModifyDischarges - ERR01'
        
        !Get concentration
        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ModifyDischarges - ERR010'
        
        allocate(DischargesConc(nDischarges, Me%nPropWithDischarge))
        DischargesConc = null_real
        
        do iDis = 1, nDischarges

            call GetDischargesGridLocalization(Me%ObjDischarges,                        &
                                               DischargeIDNumber = iDis,                &
                                               Igrid = i,                               &
                                               JGrid = j,                               &
                                               KGrid = kd,                              &
                                               DischVertical = DischVertical,           &
                                               STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ModifyDischarges - ERR020'
            
            !do not process runoff discharges (K=0). if uniform, K_layer is not used but instead k_min and K_max 
            !and user may forget K_layer zero            
            if ((DischVertical == DischUniform_) .or. (kd /= 0)) then

                call GetDischargeFlowDistribuiton(Me%ObjDischarges, iDis, nCells, FlowDistribution, &
                                                  VectorI, VectorJ, VectorK, kmin, kmax, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ModifyDischarges - ERR040'             
               
 dn:            do n=1, nCells
 
                    if (nCells > 1) then
                        i         = VectorI(n)
                        j         = VectorJ(n)
                        kd        = VectorK(n)                       
                    endif

                    if (DischVertical == DischUniform_) then

                        if (kmin == FillValueInt) kmin = Me%ExtVar%KFloor(i, j)
                        if (kmax == FillValueInt) kmax = Me%WorkSize%KUB

                    else
            
                        kmin = kd; kmax = kd

                    endif

dk:                 do k=kmin, kmax

                        if (Me%ExtVar%WaterPoints3D(i, j, k) /= WaterPoint)  Cycle

                        nullify (Property)
                        Property => Me%FirstProperty
                        iProp = 0
                        do while (associated (Property))
                            
                            if (Property%Evolution%Discharges) then 
                                
                                iProp = iProp + 1

                                !Gets Discharge Concentration for this cycle of iter
                                call GetDischargeConcentration (Me%ObjDischarges,                           &
                                                                Me%ExtVar%Now,                              &
                                                                iDis, DischargesConc(iDis, iProp),          &
                                                                Property%ID%IDNumber,                       &
                                                                STAT = STAT_CALL)
                                if (STAT_CALL/=SUCCESS_) then
                                    if (STAT_CALL == NOT_FOUND_ERR_) then 
                                        !When a property is not found associated to a discharge
                                        !by default is consider that the concentration is zero
                                        DischargesConc(iDis, iProp) = 0.
                                    else
                                        stop 'ModulePorousMediaProperties - ModifyDischarges - ERR030'
                                    endif
                                endif

                                !In case of negative discharge flux for mass balance is done using old concentration
                                !and before concentration is updated in routine DischargeProperty
                                !Do not move this computation to after DischargeProperty
                                !In case of positive use dicharge concentration
                                if (Me%CheckGlobalMass) then
                                    if (FlowDischarge(i,j,k) .lt. 0.0) then                        
                                        !kg = kg + m3/s * s * g/m3 * 1e-3kg/g
                                        Property%MB%TotalDischargeMass = Property%MB%TotalDischargeMass + (FlowDischarge(i,j,k)   &
                                                                         * Me%ExtVar%DT * Property%ConcentrationOld(i,j,k))
                                    else
                                        !kg = kg + m3/s * s * g/m3 * 1e-3kg/g
                                        Property%MB%TotalDischargeMass = Property%MB%TotalDischargeMass + (FlowDischarge(i,j,k)  &
                                                                         * Me%ExtVar%DT * DischargesConc(iDis, iProp))
                                    
                                    endif
                                endif
                                
                                !initial volume in runoff - no discharge
                                !m3 = m3H20/m3cell * m3cell
                                VolumeOld = Me%ExtVar%WaterContentOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k)
                                
                                !Update old concentration (same as doing new concentrarion and then old = new before other processes)
                                call DischargeProperty (FlowDischarge(i,j,k), DischargesConc(iDis, iProp),        &
                                                        i, j, k, VolumeOld,   Property, Me%ExtVar%DT) !, .false.)
                                
                            end if
                                            
                            Property => Property%Next

                        enddo
                    
                    enddo dk
                    
                enddo dn

                call UnGetDischarges(Me%ObjDischarges, VectorI, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModulePorousMedia - ModifyWaterDischarges - ERR070'

                call UnGetDischarges(Me%ObjDischarges, VectorJ, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModulePorousMedia - ModifyWaterDischarges - ERR080'

                call UnGetDischarges(Me%ObjDischarges, VectorK, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModulePorousMedia - ModifyWaterDischarges - ERR090'
          
           endif
           
        enddo
        
        deallocate (DischargesConc)
        
        call UngetPorousMedia (Me%ObjPorousMedia, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ModifyDischarges - ERR040'

    end subroutine ModifyDischarges  
    
    !--------------------------------------------------------------------------

    subroutine DischargeProperty (DischargeFlow, DischargeConc, Igrid, Jgrid, Kgrid, VolumeOld,         &
                                   Property, LocalDT) !, Accumulate)
        !Arguments--------------------------------------------------------------
        real                                        :: DischargeFlow, DischargeConc
        real(8)                                     :: VolumeOld
        integer                                     :: Igrid, Jgrid, Kgrid
        type (T_Property), pointer                  :: Property
        real                                        :: LocalDT
        !logical                                     :: Accumulate

        !Local------------------------------------------------------------------
        real(8)                                     :: DischargeVolume
        real(8)                                     :: OldMass, NewMass
        real                                        :: Concentration
        
        Concentration = Property%ConcentrationOld(Igrid,Jgrid,Kgrid) 

        if (abs(DischargeFlow) > AllmostZero) then
            
            ![m3] = [s] * [m3/s]
            DischargeVolume  = dble(LocalDT)*dble(DischargeFlow)
            
            ![g] = [g/m3] * [m3]
            OldMass          = dble(Concentration) * VolumeOld            
        
            if      (DischargeFlow > 0.0) then

                !Explicit discharges input 
                ![g] = [g] + [m3] * [g/m3]
                NewMass          = OldMass + DischargeVolume * dble(DischargeConc)                                       
                
                ![g/m3] = [g] / (m3 + m3/s * s)
                Concentration = NewMass / (VolumeOld + DischargeFlow * LocalDT)

            elseif (DischargeFlow < 0.0 .and. VolumeOld > 0.0) then
                    
                !If the discharge flow is negative (Output) then the concentration
                !to consider is the concentration of the cell where the discharge
                !is located

                !There is no accumulation since there is no particulate properties dissolved
                
!                if (Accumulate) then
!                    NewMass          = OldMass
!                else
                    NewMass          = OldMass * (1.0 + DischargeVolume / VolumeOld)
!                endif
                
                !if water remains
                if (abs(DischargeVolume) < VolumeOld) then
                   
                   ![g/m3] = [g] / [m3]
                    Concentration    = NewMass / (VolumeOld + DischargeVolume)
                
                else   !if all water exits node than accumulated mass needs to be accounted in bottom!
                    
                    Concentration  = 0.0
                    
!                    if (Accumulate) then
!                        ![kg/m2] = [kg/m2] + [g] * 1e-3 [kg/g] / m2
!                        Property%BottomConcentration(Igrid,Jgrid) = Property%BottomConcentration(Igrid,Jgrid) +    &
!                                                                    (NewMass * 1e-3 / Me%ExtVar%Area(Igrid,Jgrid))
!                    endif
                endif
            endif

        else
            
            !Do Nothing            

        endif
                    
        !Update concOld with discharge
        Property%ConcentrationOld(Igrid,Jgrid,KGrid) = Concentration


    end subroutine DischargeProperty

    !---------------------------------------------------------------------------

    subroutine ComputeECw 
    
        !Local--------------------------------------------------------------------
        integer                   :: i, j, k
        type(T_Property), pointer :: PropertyX

        !-------------------------------------------------------------------------        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,k) == WaterPoint) then             
                Me%ECw(i, j, k) = 0.0
                PropertyX => Me%FirstProperty
                do while (associated (PropertyX))
                    Me%ECw(i, j, k) = Me%ECw(i, j, k) + &
                          (PropertyX%Concentration(i, j, k) * PropertyX%ECwFactor)                
                    PropertyX => PropertyX%Next                    
                end do
            endif
        enddo
        enddo
        enddo
        !-------------------------------------------------------------------------

    endsubroutine ComputeECw

    !-----------------------------------------------------------------------------

    subroutine ComputeVolumes

        !Local-----------------------------------------------------------------
        integer                                   :: i, j, k, CHUNK, STAT_CALL        
        real, dimension(:,:,:), pointer           :: FlowDischarge

        !----------------------------------------------------------------------

        !The discharge flows have to be added because they are accounted separately 
        !see ModifyDischarges. Advection/Diffusion will onlye act between WaterContentBT 
        !(Before Transport) and WaterColumnFinal (After Transport and other fluxes)

        !Get integrated flow from runoff to be sure using same values
        call GetFlowDischarge (Me%ObjPorousMedia, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ComputeVolumes - ERR01'
        
        
        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        !$OMP PARALLEL PRIVATE(I,J,K)
        
        !Compute volumes and correct top flux taking FluxW(KUB+1) because it would be interpreted by module advection diffusion
        !as an additional water flux with the conc of C(i,j,k)
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,k) == WaterPoint) then             
                 
                Me%WaterVolume(i,j,k)        = Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)
                
                !m3H20/m3cell = m3H20/m3cell + m3H20/s * s /m3cell
                Me%WaterContentBT(i,j,k)     = Me%ExtVar%WaterContentOld(i,j,k) +   &
                                            (FlowDischarge(i,j,k) * Me%ExtVar%DT / Me%ExtVar%Cellvolume(i,j,k))
                
            endif
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL 

        call UngetPorousMedia (Me%ObjPorousMedia, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ComputeVolumes - ERR040'

               
        !Correct fluxw - take FluxW(KUB+1) because it would be interpreted by module advection diffusion
        !as an additional water flux with the conc of C(i,j,k)
        
        k = Me%WorkSize%KUB
        call SetMatrixValue (Me%FluxWCorr, Me%Size, Me%ExtVar%FluxW)
         
         CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(I,J)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%OpenPoints3D(i,j,k) == OpenPoint) then  
                Me%FluxWCorr(i,j,k+1) = 0.0
            endif
        enddo
        enddo        
        !$OMP END DO
        !$OMP END PARALLEL 

   
    end subroutine ComputeVolumes
    
    !----------------------------------------------------------------------

    subroutine ComputeThetaAtFaces
    
        select case (Me%ThetaAtFacesMethod)

        case (1) !Minimun

            call ComputeThetaAtFacesByMin

        case (2) !Average

            call ComputeThetaAtFacesByAvg            

        end select
    
    end subroutine ComputeThetaAtFaces

    !----------------------------------------------------------------------

    subroutine ComputeThetaAtFacesByAvg
    
        !Local-----------------------------------------------------------------
        integer                         :: I,J,K, CHUNK
        real, pointer, dimension(:,:,:) :: ThetaOld, Theta
        !Begin-----------------------------------------------------------------
        
        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K)

        !ThetaOld  => Me%ExtVar%WaterContentOld
        ThetaOld  => Me%WaterContentBT
        Theta     => Me%ExtVar%WaterContent

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                !Z direction
                
                !condition so that does not produce NAN, in WaterContetFace in the boundary
                !in boundary, diffusivity is zero (computefaces) but in release version model when evaluates disp/watercontent_face
                !produces NAN. In debug version this does not appen and the result is the total evaluation (zero).
                if (Me%ExtVar%WaterPoints3D(I,J,K-1) /= WaterPoint) then
                    if (Me%AdvDiff_Explicit) then
                        Me%ThetaAtFaces%ThetaW(i, j, k) = ThetaOld(i,j,k)
                    else
                        Me%ThetaAtFaces%ThetaW(i, j, k) = Theta(i,j,k)
                    endif
                else
                    if (Me%AdvDiff_Explicit) then
                        Me%ThetaAtFaces%ThetaW(i, j, k) = ((ThetaOld(i,j,k)      * Me%ExtVar%DWZ(i,j,k))    + &
                                                           (ThetaOld(i,j,k-1)    * Me%ExtVar%DWZ(i,j,k-1))) / &
                                                          (Me%ExtVar%DWZ(i,j,k) + Me%ExtVar%DWZ(i,j,k-1))
                    else
                        Me%ThetaAtFaces%ThetaW(i, j, k) = ((Theta(i,j,k)         * Me%ExtVar%DWZ(i,j,k))    + &
                                                           (Theta(i,j,k-1)       * Me%ExtVar%DWZ(i,j,k-1))) / &
                                                          (Me%ExtVar%DWZ(i,j,k) + Me%ExtVar%DWZ(i,j,k-1))
                    endif
                endif                               
                
                !U direction
                
                if (Me%ExtVar%WaterPoints3D(I,J-1,K) /= WaterPoint) then
                    Me%ThetaAtFaces%ThetaU(i, j, k) = ThetaOld(i,j,k)
                else                                                                                       
                    Me%ThetaAtFaces%ThetaU(i, j, k) = ((ThetaOld(i,j,k)   * Me%ExtVar%DUX(i,j))    + &
                                                       (ThetaOld(i,j-1,k) * Me%ExtVar%DUX(i,j-1))) / &
                                                      (Me%ExtVar%DUX(i,j) + Me%ExtVar%DUX(i,j-1))
                endif
                
                !V direction
                
                if (Me%ExtVar%WaterPoints3D(I-1,J,K) /= WaterPoint) then
                    Me%ThetaAtFaces%ThetaV(i, j, k) = ThetaOld(i,j,k)
                else                                                           
                    Me%ThetaAtFaces%ThetaV(i, j, k) = ((ThetaOld(i,j,k)   * Me%ExtVar%DVY(i,j))    + &
                                                       (ThetaOld(i-1,j,k) * Me%ExtVar%DVY(i-1,j))) / &
                                                      (Me%ExtVar%DVY(i,j) + Me%ExtVar%DVY(i-1,j))
                endif

            endif     
                                                                               
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            K = Me%WorkSize%KUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
            
                !UpperFace computation
                Me%ThetaAtFaces%ThetaW(i, j, k) = ThetaOld(i,j,k)
                                
            endif                                                                        
            
        enddo
        enddo
            
    end subroutine ComputeThetaAtFacesByAvg

    !----------------------------------------------------------------------
    
    subroutine ComputeThetaAtFacesByMin
    
        !Local-----------------------------------------------------------------
        integer                         :: I,J,K, CHUNK
        real, pointer, dimension(:,:,:) :: ThetaOld, Theta
        !Begin-----------------------------------------------------------------
        
        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K)

        !ThetaOld  => Me%ExtVar%WaterContentOld
        ThetaOld  => Me%WaterContentBT
        Theta     => Me%ExtVar%WaterContent

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                !Z direction
                
                !condition so that does not produce NAN, in WaterContetFace in the boundary
                !in boundary, diffusivity is zero (computefaces) but in release version model when evaluates disp/watercontent_face
                !produces NAN. In debug version this does not appen and the result is the total evaluation (zero).
                if (Me%ExtVar%WaterPoints3D(I,J,K-1) /= WaterPoint) then
                    if (Me%AdvDiff_Explicit) then
                        Me%ThetaAtFaces%ThetaW(i, j, k) = ThetaOld(i,j,k)
                    else
                        Me%ThetaAtFaces%ThetaW(i, j, k) = Theta(i,j,k)
                    endif
                else
                    if (Me%AdvDiff_Explicit) then
                        Me%ThetaAtFaces%ThetaW(i, j, k) = min(ThetaOld(i,j,k), ThetaOld(i,j,k-1))
                    else
                        Me%ThetaAtFaces%ThetaW(i, j, k) = min(Theta(i,j,k), Theta(i,j,k-1))
                    endif
                endif                               
                
                !U direction
                
                if (Me%ExtVar%WaterPoints3D(I,J-1,K) /= WaterPoint) then
                    Me%ThetaAtFaces%ThetaU(i, j, k) = ThetaOld(i,j,k)
                else                                                                                       
                    Me%ThetaAtFaces%ThetaU(i, j, k) = min(ThetaOld(i,j,k),ThetaOld(i,j-1,k))
                endif
                
                !V direction
                
                if (Me%ExtVar%WaterPoints3D(I-1,J,K) /= WaterPoint) then
                    Me%ThetaAtFaces%ThetaV(i, j, k) = ThetaOld(i,j,k)
                else                                                           
                    Me%ThetaAtFaces%ThetaV(i, j, k) = min(ThetaOld(i,j,k),ThetaOld(i-1,j,k))
                endif

            endif     
                                                                               
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            K = Me%WorkSize%KUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
            
                !UpperFace computation
                Me%ThetaAtFaces%ThetaW(i, j, k) = ThetaOld(i,j,k)
                                
            endif                                                                        
            
        enddo
        enddo
            
    end subroutine ComputeThetaAtFacesByMin

    !----------------------------------------------------------------------
    
    subroutine ActualizePropertiesFromFile
    
        !Local--------------------------------------------------------------------        
        type (T_Property), pointer :: PropertyX 
        integer                    :: STAT_CALL   
        
        !-------------------------------------------------------------------------    

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%ID%SolutionFromFile) then
            
                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
                                       Matrix3D       = PropertyX%Concentration,    &
                                       PointsToFill3D = Me%ExtVar%WaterPoints3D,    &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModulePorousMediaProperties - ERR01'
            
            endif
            
            PropertyX => PropertyX%Next
            
        enddo
    
        !-------------------------------------------------------------------------    
    
    end subroutine ActualizePropertiesFromFile
    
    !-------------------------------------------------------------------------    

    subroutine InterfaceFluxes !routine for mass sources - from vegetation (OM, fertilization..)
        !Local--------------------------------------------------------------------
        !Begin--------------------------------------------------------------------

        if (Me%ExtVar%CoupledVegetation) then
!            if (Me%ExtVar%ComputeVegInterfaceFluxes) then
                call VegetationInterfaceFluxes
 !           endif
        endif

!        if (Me%ExtVar%CoupledDN) then
!            call DrainageNetworkInterfaceFluxes
!        endif

!        if (Me%ExtVar%CoupledRunoff) then
!            if (Me%ExtVar%ComputeRunoffInterfaceFluxes) then
!                call RunoffInterfaceFluxes
!            endif
!        endif


    end subroutine InterfaceFluxes
 
    !-----------------------------------------------------------------------------
    ! This routine solves mass sources due to vegetation. - needs BIG revision it is highly inefficient. David B.
    
    subroutine VegetationInterfaceFluxes 

        !Local--------------------------------------------------------------------
        integer                                     :: i, j, k!, CHUNK
        real                                        :: Area, RootDepth
        logical                                     :: FoundEnd
        real                                        :: BottomDepth, TopDepth
        real                                        :: GrazingNotCarbon, GrazingNitrogen, GrazingPhosphorus
        real                                        :: GrazingBiomass, GrazingCarbon
        real                                        :: DormancyNotCarbon, DormancyNitrogen, DormancyPhosphorus
        real                                        :: DormancyBiomass, DormancyCarbon
        real                                        :: HarvestKillNotCarbon, HarvestKillNitrogen, HarvestKillPhosphorus
        real                                        :: HarvestKillAerialBiomass, HarvestKillCarbon
        real                                        :: HarvestKillRootNotCarbon, HarvestKillRootNitrogen, HarvestKillRootPhosphorus
        real                                        :: HarvestKillRootBiomass, HarvestKillRootCarbon
        real                                        :: NitrogenFraction, PhosphorusFraction, RootDistribution
        real                                        :: FertilizationAmmonia, FertilizationNitrate
        real                                        :: FertilizationOrganicN, FertilizationOrganicP
        real                                        :: FertilizationMineralP
        real                                        :: NitrogenUptake, PhosphorusUptake
        real                                        :: ModelDT, VegDT, CellWaterVolume, CellSoilMass
        type (T_Property), pointer                  :: Property
        integer                                     :: STAT_CALL
        character(len=5)                            :: char_i, char_j, char_k
        character(len=15)                           :: char_conc        
        character (len = StringLength)              :: StrWarning 

        !Begin--------------------------------------------------------------------

                
        !!CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR01'  

        call GetGeometryKFloor(Me%ObjGeometry,                                          &
                               Z    = Me%ExtVar%KFloor,                                 &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "InterfaceFluxes - ModulePorousMediaProperties. ERR10")

        !s
        ModelDT         = Me%ExtVar%DT
        VegDT           = Me%ExtVar%VegetationDT

        if (Me%ExtVar%ComputeVegInterfaceFluxes) then

            !!! $OMP PARALLEL PRIVATE(I,J,K)
            !!! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
    !        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            !if (Me%ExtVar%BasinPoints(i,j) == BasinPoint .and. Me%ExtVar%SoilFluxesActive(i,j)) then
            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then    
                Area         = Me%ExtVar%Area(i,j)
                RootDepth    = Me%ExtVar%RootDepth(i,j)
                FoundEnd     = .false.
                BottomDepth  = 0.
          
    do3:        do K = Me%WorkSize%KUB, Me%ExtVar%KFloor(i,j), -1                
                
                    
                    if (FoundEnd) then
                        exit do3
                    endif
                    
                    TopDepth    = BottomDepth
                    BottomDepth = BottomDepth + Me%ExtVar%DWZ(i,j,k)
                    !If found root, let compute, will exit next iteration
                    if (BottomDepth .ge. RootDepth) then
                        FoundEnd = .true.
                        BottomDepth = RootDepth
                    endif

                    GrazingCarbon               = 0.0
                    GrazingBiomass              = 0.0
                    GrazingNotCarbon            = 0.0
                    GrazingNitrogen             = 0.0
                    GrazingPhosphorus           = 0.0
                    DormancyBiomass             = 0.0
                    DormancyCarbon              = 0.0
                    DormancyNotCarbon           = 0.0
                    DormancyNitrogen            = 0.0
                    DormancyPhosphorus          = 0.0
                    HarvestKillNotCarbon        = 0.0
                    HarvestKillNitrogen         = 0.0
                    HarvestKillPhosphorus       = 0.0                
                    FertilizationAmmonia        = 0.0
                    FertilizationNitrate        = 0.0
                    FertilizationOrganicN       = 0.0
                    FertilizationOrganicP       = 0.0
                    FertilizationMineralP       = 0.0
                    HarvestKillRootCarbon       = 0.0
                    HarvestKillRootNotCarbon    = 0.0
                    HarvestKillRootNitrogen     = 0.0
                    HarvestKillRootPhosphorus   = 0.0 
                    HarvestKillAerialBiomass    = 0.0 
                    HarvestKillCarbon           = 0.0   
                    NitrogenUptake              = 0.0
                    PhosphorusUptake            = 0.0
                                  
                    
                    if (Me%ExtVar%GrowthModel) then

                        !Fluxes only occuring in surface (aerial biomass residue from grazing, dormancy and HarvestKill; 
                        !surface fertilization)
                        if (k == Me%WorkSize%KUB) then
                    
                            !Grazing
                            GrazingNotCarbon  = 0.0
                            GrazingNitrogen   = 0.0
                            GrazingPhosphorus = 0.0
                            if (Me%ExtVar%Grazing) then
                    
                                if (Me%ExtVar%ModelNitrogen) then
                                
                                    !     mgN        = KgN/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                    GrazingNitrogen  = Me%ExtVar%GrazingNitrogen(i,j) * 1e6 * Area / 10000.
                                
                                    GrazingNotCarbon = GrazingNotCarbon + GrazingNitrogen
                    
                                endif                    
                                if (Me%ExtVar%ModelPhosphorus) then

                                    !      mgP       = KgP/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                    GrazingPhosphorus = Me%ExtVar%GrazingPhosphorus(i,j) * 1e6 * Area / 10000.

                                    GrazingNotCarbon  = GrazingNotCarbon + GrazingPhosphorus
                    
                                endif                          
                    
                                !      mg       = Kg/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                GrazingBiomass = Me%ExtVar%GrazingBiomass(i,j) * 1e6 * Area / 10000.
                    
                                GrazingCarbon  = GrazingBiomass - GrazingNotCarbon

                            endif

            !                !Dormancy
                            DormancyNotCarbon  = 0.0
                            DormancyNitrogen   = 0.0
                            DormancyPhosphorus = 0.0
                            if (Me%ExtVar%Dormancy) then
                    
                                if (Me%ExtVar%ModelNitrogen) then

                                    !      mgN       = KgN/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                    DormancyNitrogen  = Me%ExtVar%DormancyNitrogen(i,j) * 1e6 * Area / 10000.

                                    DormancyNotCarbon = DormancyNotCarbon + DormancyNitrogen
                    
                                endif                    
                                if (Me%ExtVar%ModelPhosphorus) then

                                    !      mgP       = KgP/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                    DormancyPhosphorus = Me%ExtVar%DormancyPhosphorus(i,j) * 1e6 * Area / 10000.

                                    DormancyNotCarbon  = DormancyNotCarbon + DormancyPhosphorus
                    
                                endif                          
                    
                                !      mg       = Kg/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                DormancyBiomass = Me%ExtVar%DormancyBiomass(i,j) * 1e6 * Area / 10000.
                    
                                DormancyCarbon  = DormancyBiomass - DormancyNotCarbon

                            endif


            !                !HarvestKill
                            HarvestKillNotCarbon  = 0.0
                            HarvestKillNitrogen   = 0.0
                            HarvestKillPhosphorus = 0.0
                            if (Me%ExtVar%HarvestKill) then
                    
                                if (Me%ExtVar%ModelNitrogen) then

                                    !      mgN       = KgN/ha * 1E6ug/kg * (m2) * 1ha/10000m2                     
                                    HarvestKillNitrogen  = Me%ExtVar%HarvestKillNitrogen(i,j) * 1e6 * Area / 10000.

                                    HarvestKillNotCarbon = HarvestKillNotCarbon + HarvestKillNitrogen
                    
                                endif                    
                                if (Me%ExtVar%ModelPhosphorus) then

                                    !      mgP       = KgP/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                    HarvestKillPhosphorus = Me%ExtVar%HarvestKillPhosphorus(i,j) * 1e6 * Area / 10000.

                                    HarvestKillNotCarbon  = HarvestKillNotCarbon + HarvestKillPhosphorus
                    
                                endif                          
                    
                                HarvestKillAerialBiomass = Me%ExtVar%HarvestKillAerialBiomass(i,j) * 1e6 * Area / 10000.
                    
                                !      mg       = Kg/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                HarvestKillCarbon  = HarvestKillAerialBiomass - HarvestKillNotCarbon

                            endif

            !                !Fertilization in Surface
                            FertilizationAmmonia    = 0.0
                            FertilizationNitrate    = 0.0
                            FertilizationOrganicN   = 0.0
                            FertilizationOrganicP   = 0.0
                            FertilizationMineralP   = 0.0
                            if (Me%ExtVar%Fertilization) then
                    

                                if (Me%ExtVar%ModelNitrogen) then
                        
                                    !       gN       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2                     
                                    FertilizationNitrate  = Me%ExtVar%FertilNitrateSurface(i,j) * 1e3 * Area / 10000.
                                    FertilizationAmmonia  = Me%ExtVar%FertilAmmoniaSurface(i,j) * 1e3 * Area / 10000.
                                    !      mgN
                                    FertilizationOrganicN = Me%ExtVar%FertilOrganicNSurface(i,j) * 1e6 * Area / 10000.
                        
                                endif                    
                                if (Me%ExtVar%ModelPhosphorus) then

                                    !      mgP       = KgP/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                    FertilizationOrganicP = Me%ExtVar%FertilOrganicPSurface(i,j) * 1e6 * Area / 10000.
                                    !      gP
                                    FertilizationMineralP = Me%ExtVar%FertilMineralPSurface(i,j) * 1e3 * Area / 10000.
                        
                                endif                          
                    
                            endif
                            
                   
                        !Fluxes only occuring in subsurface (fertilization in sub surface)
                        elseif (k == Me%WorkSize%KUB - 1) then

            !                !Fertilization in SubSurface
                            FertilizationAmmonia    = 0.0
                            FertilizationNitrate    = 0.0
                            FertilizationOrganicN   = 0.0
                            FertilizationOrganicP   = 0.0
                            FertilizationMineralP   = 0.0
                            if (Me%ExtVar%Fertilization) then
                    

                                if (Me%ExtVar%ModelNitrogen) then
                        
                                    !      gN       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2                     
                                    FertilizationNitrate  = Me%ExtVar%FertilNitrateSubSurface(i,j) * 1e3 * Area / 10000.
                                    FertilizationAmmonia  = Me%ExtVar%FertilAmmoniaSubSurface(i,j) * 1e3 * Area / 10000.
                                    !      mgN       = KgN/ha * 1E6mg/kg * (m2) * 1ha/10000m2 
                                    FertilizationOrganicN = Me%ExtVar%FertilOrganicNSubSurface(i,j) * 1e6 * Area / 10000.
                        
                                endif                    
                                if (Me%ExtVar%ModelPhosphorus) then

                                    !      mgP       = KgP/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                    FertilizationOrganicP = Me%ExtVar%FertilOrganicPSubSurface(i,j) * 1e6 * Area / 10000.
                                    !      gP
                                    FertilizationMineralP = Me%ExtVar%FertilMineralPSubSurface(i,j) * 1e3 * Area / 10000.
                        
                                endif                          
                    
                            endif

                        endif

                        !Root death to soil (occurrs in all layers until root end )
                        if (RootDepth .gt. 0.0) then
                            
                            if (Me%ExtVar%HarvestKill) then
                                      !      mg       = Kg/ha * 1E6mg/kg * (m2) * 1ha/10000m2                     
                                HarvestKillRootBiomass = Me%ExtVar%HarvestKillRootBiomass(i,j) * 1e6 * Area / 10000.
                                !Root distribution (based on SWAT formulation)
                                RootDistribution = (1.0 - exp(-10.0 * BottomDepth/RootDepth))/(1.0 - exp(-10.0))                   &
                                                   - (1.0 - exp(-10.0 * TopDepth/RootDepth))/(1.0 - exp(-10.0))

                                HarvestKillRootNotCarbon  = 0.0
                                HarvestKillRootNitrogen   = 0.0
                                HarvestKillRootPhosphorus = 0.0
                                if (Me%ExtVar%ModelNitrogen) then
                                    NitrogenFraction        = Me%ExtVar%NitrogenFraction (i,j)
                                    HarvestKillRootNitrogen  = HarvestKillRootBiomass * NitrogenFraction * RootDistribution
                                    HarvestKillRootNotCarbon = HarvestKillRootNotCarbon + HarvestKillRootNitrogen
                                endif

                                if (Me%ExtVar%ModelPhosphorus) then
                                    PhosphorusFraction        = Me%ExtVar%PhosphorusFraction (i,j)
                                    HarvestKillRootPhosphorus  = HarvestKillRootBiomass * PhosphorusFraction * RootDistribution
                                    HarvestKillRootNotCarbon   = HarvestKillRootNotCarbon + HarvestKillRootPhosphorus
                                endif
                        
                                !     mg
                                HarvestKillRootCarbon = HarvestKillRootBiomass - HarvestKillRootNotCarbon
                            endif
                        endif
                    endif

                    !Transpiration sink has to be done by mass balance here and not in advection diffusion
                    !as flux dependent because this flux may be disconnected from Q*C
                    !It will occur only for the specified species (nitrate and dissolved phosphorus)
                    !Plant Uptake (occurrs in all layers until root end )
                    if (Me%ExtVar%ModelNitrogen) then
                        
                        !      gN       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2                     
                        NitrogenUptake = Me%ExtVar%NitrogenUptake(i,j,k) * 1e3 * Area / 10000.
                        
                    endif
    
                    if (Me%ExtVar%ModelPhosphorus) then
                        
                        !      gP       = KgP/ha * 1E3g/kg * (m2) * 1ha/10000m2                     
                        PhosphorusUptake = Me%ExtVar%PhosphorusUptake(i,j,k) * 1e3 * Area / 10000.
    
                    endif

                    
                    !!NOW UPDATE CONCENTRATIONS WITH THE MASS FLUXES COMPUTED. ConcentrationOld used because done 
                    !previous to transport that changes to new concentration
                    

                    !m3             = m3H20/m3cell * m3cell
                    CellWaterVolume = Me%ExtVar%WaterContentOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 

                    !Soil mass to compute organic and microorganisms pools
                    call SearchProperty(Property, SoilDryDensity_        , .false., STAT = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) then
                        write(*,*)
                        write(*,*) 'need property soil dry density in porous media properties'
                        stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
                    endif
                    !kgsoil         = kg/m3  * m3cell
                    CellSoilMass    = Property%ConcentrationOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
                    
                    if (Me%ExtVar%GrowthModel) then

                        ! Property Calculation
                        !!Carbon
                        call SearchProperty (Property, RefreactaryOrganicC_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) then
                            write(*,*)
                            write(*,*) 'need property particulated refractory organic carbon in porous media properties'
                            stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR100'
                        endif                
                        !         mg/kgsoil            = mg/kgsoil + mg / kgsoil
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingCarbon + DormancyCarbon  &
                                                            + HarvestKillCarbon + HarvestKillRootCarbon) * ModelDT / VegDT)        &
                                                            / CellSoilMass)

        !                call SearchProperty (Property, LabileOrganicC_        , .false., STAT = STAT_CALL)    
        !                if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
        !
        !                Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicC)          &
        !                                                * ModelDT / VegDT) / CellSoilMass)
       
                    endif

                    !!Nitrogen
                    if (Me%ExtVar%ModelNitrogen) then
                        
                        if (Me%ExtVar%GrowthModel) then
                            call SearchProperty (Property, RefreactaryOrganicN_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property particulated refractory organic nitrogen in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR110'
                            endif                          
                            !         mg/kgsoil 
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingNitrogen           &
                                                               + DormancyNitrogen + HarvestKillNitrogen                          &
                                                               + HarvestKillRootNitrogen) * ModelDT / VegDT)                     &
                                                               / CellSoilMass)


                            call SearchProperty (Property, PON_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property particulate organic nitrogen in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR120'
                            endif                        
                            !         mg/kgsoil 
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicN)   &
                                                               * ModelDT / VegDT) / CellSoilMass)


                            call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property nitrate in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR130'
                            endif    
                            !         g/m3                = g/m3 + g / m3H20
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationNitrate     & 
                                                                - NitrogenUptake) * ModelDT / VegDT) / CellWaterVolume)
    !                                                             * ModelDT / VegDT) / CellWaterVolume)
                            
                            !avoid negative because vegetation check to mass available is done daily
                            !meanwhile the mass may exit cell by transport or transformation
                            if (Property%ConcentrationOld (i,j,k) .lt. 0.0) then

                                write(char_i, '(i4)')i
                                write(char_j, '(i4)')j
                                write(char_k, '(i4)')k
                                write(char_conc, '(ES10.3)') Property%ConcentrationOld(i,j,k)                                 
                                StrWarning = 'Soil nitrate negative concentration corrected'// &
                                             ' because vegetation uptake in cell(i,j,k)'//     &
                                               char_i//','//char_j//','//char_k//' '//char_conc

                                call SetError(WARNING_, INTERNAL_, StrWarning, OFF)                                

                                !g = g + g/m3 * m3
                                Property%Mass_created(i, j, k) = Property%Mass_Created(i, j, k)   +   &
                                                               (- Property%ConcentrationOld(i, j, k)) &
                                                               *  (Me%ExtVar%WaterContentOld(i,j,k)   &
                                                               * Me%ExtVar%CellVolume (i, j, k))                                
                                
                                !uptake is waht exists. However this fix will not be updated in plant uptake!
                                !g  = g/m3 * m3
                                NitrogenUptake = NitrogenUptake - (- Property%ConcentrationOld(i, j, k))  &
                                                 *  (Me%ExtVar%WaterContentOld(i,j,k)                     &
                                                 * Me%ExtVar%CellVolume (i, j, k))     
                                
                                Property%ConcentrationOld (i,j,k) = 0.0

  !                              write(*,*) 'WARNING: '
  !                              write(*,*) 'Soil nitrate concentration negative corrected '
  !                              write(*,*) 'because vegetation uptake in cell', i, j, k 
  !                              stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR135'
  
                            endif
                            

                            if (Me%CheckGlobalMass) then
                                !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
                                Property%MB%TranspiredMass = Property%MB%TranspiredMass + (NitrogenUptake                     &
                                                                     * 1E-3 * ModelDT / VegDT)
                            endif

                            call SearchProperty (Property, Ammonia_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property ammonia in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR140'
                            endif    
                            !         g/m3                = g/m3 + g / m3H20 
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationAmmonia)   &
                                                               * ModelDT / VegDT) / CellWaterVolume)
                        else

                            call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR142'
    
                            !         g/m3                = g/m3 + g / m3H20 
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) - (((NitrogenUptake)        &
                                                               * ModelDT / VegDT) / CellWaterVolume)

                            !avoid negative because vegetation check to mass available is done daily
                            !meanwhile the mass may exit cell by transport or transformation
                            if (Property%ConcentrationOld (i,j,k) .lt. 0.0) then

                                write(char_i, '(i4)')i
                                write(char_j, '(i4)')j
                                write(char_k, '(i4)')k
                                write(char_conc, '(ES10.3)') Property%ConcentrationOld(i,j,k)                                 
                                StrWarning = 'Soil nitrate negative concentration corrected'// &
                                             ' because vegetation uptake in cell(i,j,k)'//     &
                                               char_i//','//char_j//','//char_k//' '//char_conc
                                               
                                call SetError(WARNING_, INTERNAL_, StrWarning, OFF)                 
                                
                                !g = g + g/m3 * m3
                                Property%Mass_created(i, j, k) = Property%Mass_Created(i, j, k)   +   &
                                                               (- Property%ConcentrationOld(i, j, k)) &
                                                               *  (Me%ExtVar%WaterContentOld(i,j,k)   &
                                                               * Me%ExtVar%CellVolume (i, j, k))                                
                                !g  = g/m3 * m3 - this will not be updated in plant!
                                NitrogenUptake = NitrogenUptake - (- Property%ConcentrationOld(i, j, k))  &
                                                 *  (Me%ExtVar%WaterContentOld(i,j,k)                     &
                                                 * Me%ExtVar%CellVolume (i, j, k))        
                                
                                Property%ConcentrationOld (i,j,k) = 0.0

  !                              write(*,*) 'WARNING: '
  !                              write(*,*) 'Soil nitrate concentration negative corrected '
  !                              write(*,*) 'because vegetation uptake in cell', i, j, k
  !                              stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR135'
                            endif

                            if (Me%CheckGlobalMass) then
                                !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
                                Property%MB%TranspiredMass = Property%MB%TranspiredMass + (NitrogenUptake                      &
                                                                    *  1E-3 * ModelDT / VegDT)
                            endif
                                                               
                        endif

                    endif
                    
                    !Phosphorus
                    if (Me%ExtVar%ModelPhosphorus) then

                        if (Me%ExtVar%GrowthModel) then

                            call SearchProperty (Property, RefreactaryOrganicP_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property particulated refractory organic phosphorus in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR150'
                            endif    
                            !         mg/kgsoil  
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingPhosphorus         &
                                                               + DormancyPhosphorus + HarvestKillPhosphorus                      &
                                                               + HarvestKillRootNitrogen) * ModelDT / VegDT)                     &
                                                               / CellSoilMass)


                            call SearchProperty (Property, POP_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property particulate organic phosphorus in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR160'
                            endif 
                            !         mg/kgsoil  
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicP)    &
                                                               * ModelDT / VegDT) / CellSoilMass)


                            call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property inorganic phosphorus in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR170'
                            endif 
                            !         g/m3                    = g/m3 + g / m3H20  
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationMineralP   &
                                                                - PhosphorusUptake) * ModelDT / VegDT) / CellWaterVolume)   
    !                                                             * ModelDT / VegDT) / CellWaterVolume)   
                            
                            !avoid negative because vegetation check to mass available is done daily
                            !meanwhile the mass may exit cell by transport or transformation
                            if (Property%ConcentrationOld (i,j,k) .lt. 0.0) then

                                write(char_i, '(i4)')i
                                write(char_j, '(i4)')j
                                write(char_k, '(i4)')k
                                write(char_conc, '(ES10.3)') Property%ConcentrationOld(i,j,k)                                 
                                StrWarning = 'Soil inorganic phosphorus negative concentration corrected'// &
                                             ' because vegetation uptake in cell(i,j,k)'//                  &
                                               char_i//','//char_j//','//char_k//' '//char_conc
                                               
                                call SetError(WARNING_, INTERNAL_, StrWarning, OFF)                 
                            
                                !g = g + g/m3 * m3
                                Property%Mass_created(i, j, k) = Property%Mass_Created(i, j, k)   +   &
                                                               (- Property%ConcentrationOld(i, j, k)) &
                                                               *  (Me%ExtVar%WaterContentOld(i,j,k)   &
                                                               * Me%ExtVar%CellVolume (i, j, k))                                
                                !g  = g/m3 * m3 - this will not be updated in plant!
                                PhosphorusUptake = PhosphorusUptake - (- Property%ConcentrationOld(i, j, k))  &
                                                 *  (Me%ExtVar%WaterContentOld(i,j,k)                         &
                                                 * Me%ExtVar%CellVolume (i, j, k))      
                                
                                Property%ConcentrationOld (i,j,k) = 0.0

  !                              write(*,*) 'WARNING: '
  !                              write(*,*) 'Soil dissolved phosphorus concentration negative corrected '
  !                              write(*,*) 'because vegetation uptake in cell', i, j, k
  !                              stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR135'
                            endif

                            if (Me%CheckGlobalMass) then
                                !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
                                Property%MB%TranspiredMass = Property%MB%TranspiredMass + (PhosphorusUptake                    &
                                                                    *  1E-3 * ModelDT / VegDT)
                            endif                                                            
                                                                
                        else
                            
                            call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                            if (STAT_CALL /= SUCCESS_) then
                                write(*,*)
                                write(*,*) 'need property inorganic phosphorus in porous media properties'
                                stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR180'
                            endif 
                            !         g/m3                    = g/m3 + g / m3H20  
                            Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) - (((PhosphorusUptake)      &
                                                               * ModelDT / VegDT) / CellWaterVolume) 


                            !avoid negative because vegetation check to mass available is done daily
                            !meanwhile the mass may exit cell by transport or transformation
                            if (Property%ConcentrationOld (i,j,k) .lt. 0.0) then
                            
                                write(char_i, '(i4)')i
                                write(char_j, '(i4)')j
                                write(char_k, '(i4)')k
                                write(char_conc, '(ES10.3)') Property%ConcentrationOld(i,j,k)                                 
                                StrWarning = 'Soil inorganic phosphorus negative concentration corrected'// &
                                             ' because vegetation uptake in cell(i,j,k)'//                  &
                                               char_i//','//char_j//','//char_k//' '//char_conc
                                               
                                call SetError(WARNING_, INTERNAL_, StrWarning, OFF)
                                                            
                                !g = g + g/m3 * m3
                                Property%Mass_created(i, j, k) = Property%Mass_Created(i, j, k)   +   &
                                                               (- Property%ConcentrationOld(i, j, k)) &
                                                               *  (Me%ExtVar%WaterContentOld(i,j,k)   &
                                                               * Me%ExtVar%CellVolume (i, j, k))                                
                                !g  = g/m3 * m3 - this will not be updated in plant!
                                PhosphorusUptake = PhosphorusUptake - (- Property%ConcentrationOld(i, j, k))  &
                                                 *  (Me%ExtVar%WaterContentOld(i,j,k)                         &
                                                 * Me%ExtVar%CellVolume (i, j, k))      
                                
                                Property%ConcentrationOld (i,j,k) = 0.0

  !                              write(*,*) 'WARNING: '
  !                              write(*,*) 'Soil dissolved phosphorus concentration negative corrected '
  !                              write(*,*) 'because vegetation uptake in cell', i, j, k
  !                              stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR135'
                            endif                                                             
                                                               
                            if (Me%CheckGlobalMass) then
                                !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
                                Property%MB%TranspiredMass = Property%MB%TranspiredMass + (PhosphorusUptake                   &
                                                                    * 1E-3 * ModelDT / VegDT)
                            endif                                                                                      
                         
                        endif             

                    endif

                enddo do3

            endif
            
            enddo
            enddo
        endif

        !Pesticide only in surface and only if no water column
        if (Me%ExtVar%GrowthModel) then

            Property => Me%FirstProperty

            do while (associated(Property)) 
                
                if (Property%Pesticide) then
                
                    !if (Property%Evolution%AdvectionDiffusion) then

                        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                            
                            !Here it is needed the ne Water
                            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint .and. Me%ExtVar%WaterColumnOld(i,j) .le. AlmostZero) then
                            
                                k = Me%WorkSize%KUB
                                !m3             = m3H20/m3cell * m3cell
                                CellWaterVolume = Me%ExtVar%WaterContentOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
                                !g/m3                          = g/m3  + (kg/ha * dt/vegdt * 1E3g/kg * (m2) * 1ha/10000m2) / m3 
                                Property%ConcentrationOld(i,j,k) = Property%ConcentrationOld(i,j,k) + (Property%PesticideFlux(i,j) &
                                                                    * ModelDT / VegDT * 1e3 * Me%ExtVar%Area(i,j) / 10000.)        &
                                                                     / CellWaterVolume
                                                                
                            endif
                        enddo
                        enddo
!                    else
!                        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!                        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!                            
!                            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
!                            
!                                k = Me%WorkSize%KUB
!                                !m3             = m3H20/m3cell * m3cell
!                                CellWaterVolume = Me%ExtVar%WaterContentOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
!                                !g/m3                          = g/m3  + (kg/ha * dt/vegdt * 1E3g/kg * (m2) * 1ha/10000m2) / m3 
!                                Property%Concentration(i,j,k) = Property%ConcentrationOld(i,j,k) + (Property%PesticideFlux(i,j) &
!                                                                    * ModelDT / VegDT * 1e3 * Me%ExtVar%Area(i,j) / 10000.)     &
!                                                                     / CellWaterVolume
!                                                                
!                            endif
!                        enddo
!                        enddo                    
!                    endif         
                endif
            
                Property => Property%Next
                
            enddo                    
            
        endif

        !Properties that do not have advection difusion will not be updated (advection diffusion updates conc old
        !to conc new). So this properties are here updated to conc new
        if (Me%ExtVar%GrowthModel) then

            Property => Me%FirstProperty

            do while (associated(Property)) 
                
                if (.not. Property%Evolution%AdvectionDiffusion) then

                    call SetMatrixValue (Property%Concentration, Me%Size, Property%ConcentrationOld)
                    
                endif         
            
                Property => Property%Next
                
            enddo                    
            
        endif


        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR1070'  

        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%KFloor,       STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_, "InterfaceFluxes - ModulePorousMediaProperties. ERR1080")


    end subroutine VegetationInterfaceFluxes

    !-----------------------------------------------------------------------------
    
    subroutine AdvectionDiffusionProcesses

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        
        !begin--------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "AdvectionDiffusionProcesses")
        

       !Compute water volume and remove infiltration from fluxZ 
       !(it would create error in AdvectionDiffusion routines)
       !Also remove discharges because they were explicity computed
        call ComputeVolumes
        
        call ComputeThetaAtFaces
        
        if (Me%NewFormulation) then
            call ComputeAdvectionTerms
        endif
                 
             
        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%Evolution%AdvectionDiffusion) then
            
                !Restart matrixes for computation 
                call RestartVariables(PropertyX)
                
                !Compute the coefs needed for transport. Coefs assciated to all PorousMedia fluxes
                call ModifyCoefs(PropertyX)
            
                !Update property values based on the new coefs computed
                call ModifyPropertyValues (PropertyX)
               
                !Update property mass fluxes between modules
                if (Me%CheckGlobalMass) then
                    call ModifyInterfaceMassFluxes (PropertyX)
                endif
            
            endif


            PropertyX => PropertyX%Next

        enddo

        if (Me%Coupled%MinConcentration)     call SetLimitsConcentration  
        if (Me%Coupled%WarnOnNegativeValues) call WarnOnNegativeValues    ('After Advection Diffusion')
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "AdvectionDiffusionProcesses")


    end subroutine AdvectionDiffusionProcesses

    !-----------------------------------------------------------------------------

    subroutine RestartVariables(PropertyX)
        
        !Argument-----------------------------------------------------------------
         type (T_Property), pointer                  :: PropertyX

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        !begin--------------------------------------------------------------------
        
        CurrProperty => PropertyX
        
!        call SetMatrixValue (Me%COEFExpl%CoefA_W, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefB_W, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefC_W, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefA_U, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefB_U, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefC_U, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefA_V, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefB_V, Me%Size, 0.0)
!        call SetMatrixValue (Me%COEFExpl%CoefC_V, Me%Size, 0.0)
        call SetMatrixValue (Me%COEFExpl%CoefInterfRunoff, Me%Size, 0.0)
        call SetMatrixValue (Me%COEFExpl%CoefInterfDN, Me%Size, 0.0)
        call SetMatrixValue (Me%COEFExpl%CoefInterfTransp, Me%Size, 0.0)
        call SetMatrixValue (Me%COEFExpl%CoefInterfEvap, Me%Size, 0.0)
        call SetMatrixValue (Me%COEFExpl%CoefInterfBoundaryWalls, Me%Size, 0.0)
        call SetMatrixValue (Me%COEFExpl%CoefInterfBoundaryBottom, Me%Size, 0.0)
        call SetMatrixValue (CurrProperty%ConcInInterfaceDN, Me%Size, 0.0)
        call SetMatrixValue (CurrProperty%ConcInBoundary, Me%Size, 0.0)
        
        call SetMatrixValue (Me%TICOEF3, Me%Size, 0.0)
        
        if (.not. Me%NewFormulation) then
            call SetMatrixValue (Me%COEF3_VertAdv%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%F_flux, Me%Size, 0.0)
            
            if (.not. Me%Vertical1D) then
                call SetMatrixValue (Me%COEF3_HorAdvXX%C_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_HorAdvXX%D_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_HorAdvXX%E_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_HorAdvXX%F_flux, Me%Size, 0.0)
                
                call SetMatrixValue (Me%COEF3_HorAdvYY%C_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_HorAdvYY%D_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_HorAdvYY%E_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_HorAdvYY%F_flux, Me%Size, 0.0)                
                
            endif
        endif
        
        if (.not. Me%AdvDiff_Explicit) then
            call SetMatrixValue (Me%COEF3%D, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3%E, Me%Size, dble(1.0))
            call SetMatrixValue (Me%COEF3%F, Me%Size, 0.0)

        endif

      if (Me%Output%Boxes_ON) then
            
            call SetMatrixValue (Me%Fluxes%AdvFluxX, Me%Size, dble(0.0))
            call SetMatrixValue (Me%Fluxes%AdvFluxY, Me%Size, dble(0.0))
            call SetMatrixValue (Me%Fluxes%AdvFluxZ, Me%Size, dble(0.0))

            call SetMatrixValue (Me%Fluxes%DifFluxX, Me%Size, dble(0.0))
            call SetMatrixValue (Me%Fluxes%DifFluxY, Me%Size, dble(0.0))
            call SetMatrixValue (Me%Fluxes%DifFluxZ, Me%Size, dble(0.0))
            
            call SetMatrixValue (Me%Fluxes%MassFluxesX, Me%Size, dble(0.0))
            call SetMatrixValue (Me%Fluxes%MassFluxesY, Me%Size, dble(0.0))
            call SetMatrixValue (Me%Fluxes%MassFluxesZ, Me%Size, dble(0.0))
        end if 
    
    end subroutine RestartVariables

    !-----------------------------------------------------------------------------

    subroutine ComputeAdvectionTerms         

        !Local--------------------------------------------------------------------
        character (Len = StringLength)           :: Direction
        !Begin--------------------------------------------------------------------


            
         if (.not. Me%Vertical1D) then

            Direction = 'Horizontal'
            
            select case (Me%AdvDiff_AdvMethodH)
            
                case (UpwindOrder1)
                    
                    call AdvectionUpwindOrder1 (Direction)
                
                case (CentralDif)
                
                    call AdvectionCentralDifferences(Direction)
                    
                case default
                    
                    write(*,*) 'Undefined method for horizontal advection check ADVDIFF_METHOD_H'
                    stop 'AdvectionDiffusionProcesses - Porous Media Properties - ERR01'
                
             end select
             
        endif

        Direction = 'Vertical'

        select case (Me%AdvDiff_AdvMethodV)
        
            case (UpwindOrder1)
                
                call AdvectionUpwindOrder1 (Direction)
            
            case (CentralDif)
            
                call AdvectionCentralDifferences(Direction)
                
            case default
                
                write(*,*) 'Undefined method for vertical advection check ADVDIFF_METHOD_V'
                stop 'AdvectionDiffusionProcesses - Porous Media Properties - ERR010'
            
         end select        
    
    end subroutine ComputeAdvectionTerms

    !-----------------------------------------------------------------------------

    subroutine AdvectionUpwindOrder1(Direction)
        
        !Argument-----------------------------------------------------------------
        character(len=StringLength)         :: Direction
        !Local--------------------------------------------------------------------
        integer                             :: i,j,k, CHUNK                             
       
        !begin-------------------------------------------------------------------- 

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "AdvectionUpwindOrder1")

        
        !By default at this point all adv fluxes are zero (routine RestartVariables) so only values different from zero are defined
        
        if (Direction == 'Horizontal') then

            call SetMatrixValue (Me%COEF3_HorAdvXX%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%F_flux, Me%Size, 0.0)
            
            call SetMatrixValue (Me%COEF3_HorAdvYY%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%F_flux, Me%Size, 0.0)                
            
            
            CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB,Me%WorkSize%KUB)
            !$OMP PARALLEL PRIVATE(i,j,k) 
     
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :       do k = Me%WorkSize%KLB, Me%WorkSize%KUB
do2 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                    
                    !x direction
                    if (Me%ExtVar%ComputeFacesU3D(i,j,k) == 1) then
                        if (Me%ExtVar%FluxU(i,j,k) .gt. 0.0) then
                            Me%COEF3_HorAdvXX%D_flux(i,j,k) = Me%ExtVar%FluxU(i,j,k)
                        else
                            Me%COEF3_HorAdvXX%E_flux(i,j,k) = Me%ExtVar%FluxU(i,j,k)
                        endif    
                    endif
                    
                    !y direction
                    if (Me%ExtVar%ComputeFacesV3D(i,j,k) == 1) then                
                        if (Me%ExtVar%FluxV(i,j,k) .gt. 0.0) then
                            Me%COEF3_HorAdvYY%D_flux(i,j,k) = Me%ExtVar%FluxV(i,j,k)
                        else
                            Me%COEF3_HorAdvYY%E_flux(i,j,k) = Me%ExtVar%FluxV(i,j,k)
                        endif   
                    endif                 
                    
                endif

            end do do1
            end do do2
            end do do3
            !$OMP END DO
                
            !$OMP END PARALLEL
            
        elseif (Direction == 'Vertical') then ! z direction

            call SetMatrixValue (Me%COEF3_VertAdv%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%F_flux, Me%Size, 0.0)
            
            
            CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB,Me%WorkSize%KUB)
            !$OMP PARALLEL PRIVATE(i,j,k) 
     
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do6 :       do k = Me%WorkSize%KLB, Me%WorkSize%KUB
do5 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do4 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                    
                    !z direction
                    if (Me%ExtVar%ComputeFacesW3D(i,j,k) == 1) then                
                        if (Me%ExtVar%FluxW(i,j,k) .gt. 0.0) then
                            Me%COEF3_VertAdv%D_flux(i,j,k) = Me%ExtVar%FluxW(i,j,k)
                        else
                            Me%COEF3_VertAdv%E_flux(i,j,k) = Me%ExtVar%FluxW(i,j,k)
                        endif   
                    endif                                 
                    
                endif

            end do do4
            end do do5
            end do do6
            !$OMP END DO
                
            !$OMP END PARALLEL
        
        else
            stop 'AdvectionUpwindOrder1 - Porous Media Properties - ERR010'
        endif

        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "AdvectionUpwindOrder1")
          
           
    end subroutine AdvectionUpwindOrder1
    
    !-----------------------------------------------------------------------------

    subroutine AdvectionCentralDifferences(Direction)

        !Argument-----------------------------------------------------------------
        character(len=StringLength)         :: Direction
        
        !Local--------------------------------------------------------------------
        integer                             :: i,j,k, CHUNK   
        real                                :: Aux1, Aux2
        real, dimension(:,:  ), pointer     :: DUX, DVY
        real, dimension(:,:,:), pointer     :: DWZ                      
       
        !begin-------------------------------------------------------------------- 

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "AdvectionCentralDifferences")

        
        !By default at this point all adv fluxes are zero (routine RestartVariables) so only values different from zero are defined
        
        !xx and yy direction
        if (Direction == 'Horizontal') then

            call SetMatrixValue (Me%COEF3_HorAdvXX%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%F_flux, Me%Size, 0.0)
            
            call SetMatrixValue (Me%COEF3_HorAdvYY%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%F_flux, Me%Size, 0.0)                
            
        
            DUX => Me%ExtVar%DUX
            DVY => Me%ExtVar%DVY
            
            CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB,Me%WorkSize%KUB)
            !$OMP PARALLEL PRIVATE(i,j,k, Aux1, Aux2) 
     
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :       do k = Me%WorkSize%KLB, Me%WorkSize%KUB
do2 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                    
                    !x direction
                    if (Me%ExtVar%ComputeFacesU3D(i,j,k) == 1) then
                        Aux1 = DUX(i,j  )/(DUX(i,j) + DUX(i,j-1))
                        Aux2 = DUX(i,j-1)/(DUX(i,j) + DUX(i,j-1))
                        Me%COEF3_HorAdvXX%D_flux(i,j,k) = Me%ExtVar%FluxU(i,j,k) * Aux1
                        Me%COEF3_HorAdvXX%E_flux(i,j,k) = Me%ExtVar%FluxU(i,j,k) * Aux2
                    endif
                    
                    !y direction
                    if (Me%ExtVar%ComputeFacesV3D(i,j,k) == 1) then                
                        Aux1 = DVY(i  ,j)/(DVY(i,j) + DVY(i-1,j))
                        Aux2 = DVY(i-1,j)/(DVY(i,j) + DVY(i-1,j))
                        Me%COEF3_HorAdvYY%D_flux(i,j,k) = Me%ExtVar%FluxV(i,j,k) * Aux1
                        Me%COEF3_HorAdvYY%E_flux(i,j,k) = Me%ExtVar%FluxV(i,j,k) * Aux2
                    endif                 
                    
                    
                endif

            end do do1
            end do do2
            end do do3
            !$OMP END DO
                
            !$OMP END PARALLEL
        
        elseif (Direction == 'Vertical') then 

            call SetMatrixValue (Me%COEF3_VertAdv%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_VertAdv%F_flux, Me%Size, 0.0)
       
            !Z direction can't be parallelized in the same cycle as XX and YY (K)
            DWZ => Me%ExtVar%DWZ

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(i,j,k,Aux1,Aux2) 
     
do6 :       do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do5 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do4 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                    
                    !z direction
                    if (Me%ExtVar%ComputeFacesW3D(i,j,k) == 1) then 
                        Aux1 = DWZ(i,j,k  )/(DWZ(i,j,k) + DWZ(i,j,k-1))
                        Aux2 = DWZ(i,j,k-1)/(DWZ(i,j,k) + DWZ(i,j,k-1))
                        Me%COEF3_VertAdv%D_flux(i,j,k) = Me%ExtVar%FluxW(i,j,k) * Aux1
                        Me%COEF3_VertAdv%E_flux(i,j,k) = Me%ExtVar%FluxW(i,j,k) * Aux2
                    endif                              
                    
                endif

            end do do4
            end do do5
            !$OMP END DO
            end do do6
                
            !$OMP END PARALLEL
        else
            stop 'AdvectionCentralDifferences - Porous Media Properties - ERR010'        
        endif
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "AdvectionCentralDifferences")
           
           
    end subroutine AdvectionCentralDifferences
    
    !-----------------------------------------------------------------------------

    subroutine ModifyCoefs(PropertyX)

        !Argument-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
       
        !begin--------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyCoefs")
        
        !In this subroutine are computed all the sources/sinks of water and mass
        !that exist in the Porous Media Module (for new theta computation)
               
        !Diffusivity will be used in Infiltration coef and in AdvectionDiffusion
        call ModifyDiffusivity_New(PropertyX)

        !Fluxes in X, Y and Z direction explicit or implicit. 
        call ModifyAdvectionDiffusionCoefs(PropertyX)

        !Evaporation
        !Evaporation fluxes do not take mass and do not need to be accounted (e.g. salinity)
        !But some properties as temperature and oxygen are transported - see below
        
        !Transpiration fluxes - in cells along roots - take all dissolved properties
        !Vegetation may remove mass not associated to Q*C and be selective for species
        !so this option was abandoned and nitrate and dissolved phosphorus (only) removed in
        !routine VegetationInterfaceFluxes as mass sink
        !if (Me%ExtVar%CoupledVegetation .and. Me%ExtVar%ModelWater ) then
        !    call ModifyVegetationCoefs(PropertyX)
        !endif
        
        !Some special dissolved properties are transported with evaporation and transpiration and 
        !if not accounted will produce increase concentration in soil (e.g. temperature and oxygen)
        if (PropertyX%Evolution%TransportedInEVTP) then
            call ModifyEVTPCoefs()
        endif
        
        !Infiltration - in surface cells
        call ModifyInfiltrationCoefs(PropertyX)
        
        !Fluxes with Drainage network - in the cells that link with river
        if (Me%ExtVar%CoupledDN) then
            call ModifyDrainageNetworkCoefs(PropertyX)
        endif
        
        !BoundaryFluxes
        if (Me%ExtVar%BoundaryImposed) then
            call ModifyBoundaryCoefs(PropertyX)
        endif
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyCoefs")

        
    end subroutine ModifyCoefs
    
    !-----------------------------------------------------------------------------
    
    subroutine ModifyAdvectionDiffusionCoefs(PropertyX)

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                       :: PropertyX
        real                                             :: ImpExp_AdvXX, ImpExp_AdvYY           
        real                                             :: AdvectionV_Imp_Exp  
        real                                             :: DiffusionV_Imp_Exp  
        real                                             :: AdvectionH_Imp_Exp  
        integer                                          :: di,    dj    
        integer                                          :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                          :: ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS
                        
        !begin--------------------------------------------------------------------
 
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyAdvectionDiffusionCoefs")
        
        !Explicit coefs in X and Y direction. Z direction computed if vertical is explicit
        !call ModifyExplicitCoefs (PropertyX)
        
       !check spatial and temporal options
        call CheckTransportOptions(PropertyX, AdvectionV_Imp_Exp, DiffusionV_Imp_Exp,      &
                                     AdvectionH_Imp_Exp, ImpExp_AdvXX, ImpExp_AdvYY)
        
        !Horizontal fluxes computed if not vertical 1D solution
        if (.not. Me%Vertical1D) then
            
            !Routines from ModuleAdvectionDiffusion
            call HorizontalDiffusion(PropertyX)                              !always explicit
            call HorizontalAdvection(PropertyX, ImpExp_AdvXX, ImpExp_AdvYY)  !explicit or implicit

        endif
        
        !Vertical fluxes in vertical grid
        if (Me%WorkSize%KUB .gt. 1) then    
            
            !Routines from ModuleAdvectionDiffusion
            call VerticalDiffusion(PropertyX, DiffusionV_Imp_Exp)            !implicit or explicit                       
            call VerticalAdvection(PropertyX, AdvectionV_Imp_Exp)            !implicit or explicit                           
         
         endif   

!        THOMASZ WAS MOVED TO ROUTINE THAT COMPUTES
!        CONCENTRATIONS BASED ON ALL PROCESSES (ModifyPropertyValues).
!        MAYBE THE PROBLEM WITH VERTICAL ADVECTION IMPLICIT IS DUE TO LACK OF
!        FIRST THOMASZ HERE AND RESET COEFS and THEN USE AGAIN (THOMASZ OR THOMAS3D??)

        ILBWS = Me%WorkSize%ILB
        IUBWS = Me%WorkSize%IUB
        JLBWS = Me%WorkSize%JLB
        JUBWS = Me%WorkSize%JUB
        KLBWS = Me%WorkSize%KLB
        KUBWS = Me%WorkSize%KUB
        ILB   = Me%Size%ILB
        IUB   = Me%Size%IUB
        JLB   = Me%Size%JLB 
        JUB   = Me%Size%JUB
        KLB   = Me%Size%KLB
        KUB   = Me%Size%KUB


cd3:    if (KUBWS == 1 .and. ImpExp_AdvXX == ImplicitScheme) then !ImplicitScheme = 0

            di = 0
            dj = 1

            !griflet: old call
            !call THOMAS_3D(ILBWS, IUBWS,                                                &
            !               JLBWS, JUBWS,                                                &
            !               KLBWS, KUBWS,                                                &
            !               di, dj,                                                      &
            !               Me%COEF3%D,                                                  &
            !               Me%COEF3%E,                                                  &
            !               Me%COEF3%F,                                                  &
            !               Me%TICOEF3,                                                  &
            !               PropertyX%Concentration,                                     &
            !               Me%VECG,                                                     &
            !               Me%VECW)      
            !griflet: new  call
            call THOMAS_3D(ILBWS, IUBWS,                                                &
                           JLBWS, JUBWS,                                                &
                           KLBWS, KUBWS,                                                &
                           di, dj,                                                      &
                           Me%THOMAS,                                                   &
                           PropertyX%Concentration                                      &     
#ifdef _ENABLE_CUDA
                           , Me%ObjCuda,                                                &
                           .FALSE.                                                      &
#endif _ENABLE_CUDA
                           )
                            
        else if (KUBWS == 1 .and. ImpExp_AdvYY == ImplicitScheme) then cd3 !ImplicitScheme = 0

            di = 1
            dj = 0
            
            !griflet: old call
            !call THOMAS_3D(JLBWS, JUBWS,                                                &
            !               ILBWS, IUBWS,                                                &
            !               KLBWS, KUBWS,                                                &
            !               di, dj,                                                      &
            !               Me%COEF3%D,                                                  &
            !               Me%COEF3%E,                                                  &
            !               Me%COEF3%F,                                                  &
            !               Me%TICOEF3,                                                  &
            !               PropertyX%Concentration,                                     &
            !               Me%VECG,                                                     &
            !               Me%VECW)      
            !griflet: new call                           
            call THOMAS_3D(JLBWS, JUBWS,                                                &
                           ILBWS, IUBWS,                                                &
                           KLBWS, KUBWS,                                                &
                           di, dj,                                                      &
                           Me%THOMAS,                                                   &
                           PropertyX%Concentration                                      &      
#ifdef _ENABLE_CUDA
                           , Me%ObjCuda,                                                &
                           .FALSE.                                                      &
#endif _ENABLE_CUDA
                           )
        endif cd3
!        else cd3
! 
!            ! If the model is 3D the vertical diffusion must be implicit so is necessary to 
!            ! compute the vertical diffusion  implicitly
!            
!            !griflet: old call   
!            !CALL THOMASZ(ILBWS, IUBWS,                                                  &
!            !             JLBWS, JUBWS,                                                  &
!            !             KLBWS, KUBWS,                                                  &
!            !             Me%COEF3%D,                                                    &
!            !             Me%COEF3%E,                                                    &
!            !             Me%COEF3%F,                                                    &
!            !             Me%TICOEF3,                                                    &
!            !             PropertyX%Concentration,                                       &
!            !             Me%VECG,                                                       &
!            !             Me%VECW)      
!        
!            !griflet: new call
!            call THOMASZ(ILBWS, IUBWS,                                                  &
!                         JLBWS, JUBWS,                                                  &
!                         KLBWS, KUBWS,                                                  &
!                         Me%THOMAS,                                                     &
!                         PropertyX%Concentration                                        &
!#ifdef _ENABLE_CUDA
!                         , Me%ObjCuda,                                                  &
!                         .FALSE.                                                        &
!#endif _ENABLE_CUDA
!                        )
!
!        endif cd3


cd5 :   if (Me%Output%Boxes_ON) then
            if (AdvectionV_Imp_Exp > 0.0 .and. KUBWS > 1)                        &
                call CalcVerticalAdvFlux(PropertyX, AdvectionV_Imp_Exp)

            if (DiffusionV_Imp_Exp > 0.0 .and. KUBWS > 1)                        &
                call CalcVerticalDifFlux (PropertyX, DiffusionV_Imp_Exp)

            if ((.not. Me%Vertical1D) .and. (ImpExp_AdvXX == ImplicitScheme))    &
                call CalcHorizontalAdvFluxXX(PropertyX, ImpExp_AdvXX)

            if ((.not. Me%Vertical1D) .and. (.not. Me%XZFlow) .and. (ImpExp_AdvYY == ImplicitScheme))  &
                call CalcHorizontalAdvFluxYY(PropertyX, ImpExp_AdvYY)

            call Output_Boxes_Fluxes (PropertyX)
            
        end if cd5


        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyAdvectionDiffusionCoefs")
        
        
    end subroutine ModifyAdvectionDiffusionCoefs
    
    !-----------------------------------------------------------------------------
    
    subroutine CheckTransportOptions(PropertyX, AdvectionV_Imp_Exp, DiffusionV_Imp_Exp,      &
                                     AdvectionH_Imp_Exp, ImpExp_AdvXX, ImpExp_AdvYY)
        !Local--------------------------------------------------------------------
        type (T_Property), pointer                       :: PropertyX
        real, intent(OUT)                                :: ImpExp_AdvXX, ImpExp_AdvYY           
        real, intent(OUT)                                :: AdvectionV_Imp_Exp  
        real, intent(OUT)                                :: DiffusionV_Imp_Exp  
        real, intent(OUT)                                :: AdvectionH_Imp_Exp          
        !begin--------------------------------------------------------------------

        if (Me%AdvDiff_Explicit) then
            
            !if Explicit - all explicit 
            AdvectionV_Imp_Exp = ExplicitScheme
            DiffusionV_Imp_Exp = ExplicitScheme
            AdvectionH_Imp_Exp = ExplicitScheme
        
        else
            !if Implicit - vertical advection and diffusion implicit
            AdvectionV_Imp_Exp = ImplicitScheme
            DiffusionV_Imp_Exp = ImplicitScheme
            
            !check if also horizontal advection is implicit
            if (.not. Me%AdvDiff_AdvectionH_ImpExp) then
                AdvectionH_Imp_Exp = ImplicitScheme
            else
                AdvectionH_Imp_Exp = ExplicitScheme
            endif
            
            !Horizontal diffusion is always explicit
            
        endif
            
        if(AdvectionH_Imp_Exp == ImplicitScheme) then

            if(PropertyX%Evolution%AdvDiff%ImplicitH_Direction == DirectionX)then
                                           
                !Direction X implicit
                ImpExp_AdvXX = ImplicitScheme 
                ImpExp_AdvYY = ExplicitScheme 

                PropertyX%Evolution%AdvDiff%ImplicitH_Direction = DirectionY

            else 
            
                !Direction Y implicit
                ImpExp_AdvXX = ExplicitScheme 
                ImpExp_AdvYY = ImplicitScheme 

                PropertyX%Evolution%AdvDiff%ImplicitH_Direction = DirectionX

            endif 
    
        else ! Horizontal Advection Explicit

            ImpExp_AdvXX = ExplicitScheme 
            ImpExp_AdvYY = ExplicitScheme 

        endif

    
    end subroutine CheckTransportOptions

    !-----------------------------------------------------------------------------

    subroutine HorizontalDiffusion(CurrProp) !Routine from ModuleAdvectionDiffusion

        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp
        !----------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "HorizontalDiffusion")

        
        call HorizontalDiffusionXX(CurrProp)
        
        if (.not. Me%XZFlow) call HorizontalDiffusionYY(CurrProp)
                
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "HorizontalDiffusion")


    end subroutine HorizontalDiffusion

    !--------------------------------------------------------------------------

    subroutine HorizontalDiffusionXX(CurrProp) !Routine from ModuleAdvectionDiffusion

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp
        !Local-----------------------------------------------------------------
        real(8)                                     :: DTPropDouble 
        real(8)                                     :: AuxJ
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "HorizontalDiffusionXX")


        DTPropDouble = dble(Me%ExtVar%DT) 

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB,Me%WorkSize%KUB)
        !$OMP PARALLEL PRIVATE(i,j,k,AuxJ) 
 
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
do2 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j, k) == 1) then

                AuxJ = CurrProp%ViscosityU      (i,j  ,k) &
                       * Me%ExtVar%AreaU        (i,j  ,k) &
                       * Me%ThetaAtFaces%ThetaU (i,j  ,k) &
                       / Me%ExtVar%DZX          (i,j-1  )                    

                Me%TICOEF3(i,j-1,k) = Me%TICOEF3(i,j-1,k) + AuxJ * DTPropDouble /      &
                                      Me%WaterVolume(i, j-1, k) *                      &
                                     (CurrProp%Concentration(i,j,k) - CurrProp%Concentration(i,j-1,k))


                Me%TICOEF3(i,j  ,k) = Me%TICOEF3(i,j  ,k) - AuxJ * DTPropDouble /       &
                                      Me%WaterVolume(i, j  , k) *                       &
                                     (CurrProp%Concentration(i,j,k) - CurrProp%Concentration(i,j-1,k))


            endif

        end do do1
        end do do2
        end do do3
        !$OMP END DO
            
        !$OMP END PARALLEL

        if (Me%Output%Boxes_ON) call CalcHorizontalDifFluxXX(CurrProp)
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "HorizontalDiffusionXX")
    
    end subroutine HorizontalDiffusionXX

    !--------------------------------------------------------------------------

    subroutine HorizontalDiffusionYY(CurrProp)  !Routine from ModuleAdvectionDiffusion

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp
        !Local-----------------------------------------------------------------
        real(8)                                     :: DTPropDouble 
        real(8)                                     :: AuxI
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "HorizontalDiffusionYY")


        DTPropDouble = dble(Me%ExtVar%DT) 

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(i,j,k,AuxI) 

do3 :   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesV3D(i  , j, k) == 1) then

                       
                AuxI = CurrProp%ViscosityV      (i  ,j,k) &
                       * Me%ExtVar%AreaV        (i  ,j,k) &
                       * Me%ThetaAtFaces%ThetaV (i  ,j,k) &
                       / Me%ExtVar%DZY          (i-1,j  )

                Me%TICOEF3(i-1,j,k) = Me%TICOEF3(i-1,j,k) + AuxI * DTPropDouble /       &
                                      Me%WaterVolume(i-1, j, k) *                      &
                                     (CurrProp%Concentration(i,j,k) - CurrProp%Concentration(i-1,j,k))


                Me%TICOEF3(i,j  ,k) = Me%TICOEF3(i,j  ,k) - AuxI * DTPropDouble /       &
                                      Me%WaterVolume(i  , j, k) *                       &
                                     (CurrProp%Concentration(i,j,k) - CurrProp%Concentration(i-1,j,k))
            endif
        end do do1
        end do do2
        !$OMP END DO
        end do do3

        !$OMP END PARALLEL

        if (Me%Output%Boxes_ON) call CalcHorizontalDifFluxYY(CurrProp)
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "HorizontalDiffusionYY")


    end subroutine HorizontalDiffusionYY

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvection(CurrProp, ImpExp_AdvXX, ImpExp_AdvYY)  !Routine from ModuleAdvectionDiffusion

        !Arguments-------------------------------------------------------------
        real                                :: ImpExp_AdvXX, ImpExp_AdvYY
        type (T_Property), pointer          :: CurrProp

        !Local-----------------------------------------------------------------
        integer                             :: di,    dj    
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "HorizontalAdvection")

        ILBWS = Me%WorkSize%ILB
        IUBWS = Me%WorkSize%IUB

        JLBWS = Me%WorkSize%JLB
        JUBWS = Me%WorkSize%JUB

        KLBWS = Me%WorkSize%KLB
        KUBWS = Me%WorkSize%KUB

        ILB   = Me%Size%ILB
        IUB   = Me%Size%IUB

        JLB   = Me%Size%JLB 
        JUB   = Me%Size%JUB

        KLB   = Me%Size%KLB
        KUB   = Me%Size%KUB

        
        call HorizontalAdvectionXX(CurrProp, ImpExp_AdvXX)

        if (.not. Me%XZFlow) call HorizontalAdvectionYY(CurrProp, ImpExp_AdvYY)

cd1:    if (ImpExp_AdvYY == ImplicitScheme .or. ImpExp_AdvXX == ImplicitScheme) then 

cd3D:       if (KUBWS > 1) then

                !if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorAdvection-THOMAS3D")


cd2:            if (ImpExp_AdvXX == ImplicitScheme) then 

                    di = 0
                    dj = 1

                    !griflet: old call
                    !call THOMAS_3D(ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS, di, dj,    &
                    !     Me%COEF3%D,                                                    &
                    !     Me%COEF3%E,                                                    &
                    !     Me%COEF3%F,                                                    &
                    !     Me%TICOEF3,                                                    &
                    !     CurrProp%Concentration,                                        &
                    !     Me%VECG,                                                       &
                    !     Me%VECW)     
                    !griflet: new call 
                    call THOMAS_3D(ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS, di, dj,    &
                         Me%THOMAS,                                                     &
                         CurrProp%Concentration                                         &     
#ifdef _ENABLE_CUDA
                         , Me%ObjCuda,                                                  &
                         .FALSE.                                                        &
#endif _ENABLE_CUDA
                         )
                         
                else if (ImpExp_AdvYY == ImplicitScheme) then cd2

                    di = 1
                    dj = 0

                    !griflet: old call
!                    call THOMAS_3D(JLBWS, JUBWS, ILBWS, IUBWS, KLBWS, KUBWS, di, dj,    &
!                         Me%COEF3%D,                                                    &
!                         Me%COEF3%E,                                                    &
!                         Me%COEF3%F,                                                    &
!                         Me%TICOEF3,                                                    &
!                         CurrProp%Concentration,                                        &
!                         Me%VECG,                                                       &
!                         Me%VECW)  
                    !griflet: new call    
                    call THOMAS_3D(JLBWS, JUBWS, ILBWS, IUBWS, KLBWS, KUBWS, di, dj,    &
                         Me%THOMAS,                                                     &
                         CurrProp%Concentration                                         &     
#ifdef _ENABLE_CUDA
                         , Me%ObjCuda,                                                  &
                         .FALSE.                                                        &
#endif _ENABLE_CUDA
                         )
                         
                endif cd2

                
                call SetMatrixValue (Me%COEF3%D, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3%E, Me%Size, dble(1.0))
                call SetMatrixValue (Me%COEF3%F, Me%Size, 0.0)
                call SetMatrixValue (Me%TICOEF3, Me%Size, CurrProp%Concentration)

                !if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorAdvection-THOMAS3D")

            endif cd3D
     

        endif cd1

        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "HorizontalAdvection")


    end subroutine HorizontalAdvection

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvectionXX(CurrProp, ImpExp_AdvXX)    !Routine from ModuleAdvectionDiffusion

        !Arguments--------------------------------------------------------------
        type (T_Property), pointer          :: CurrProp
        real                                :: ImpExp_AdvXX

        !Local-----------------------------------------------------------------               

        real(8) :: AdvFluxX, DT1, DT2

        integer :: i,     j,     k                             
        integer :: ILB, IUB, JLB, JUB, KLB, KUB
        integer :: CHUNK
        !----------------------------------------------------------------------


        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "HorizontalAdvectionXX")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
        
        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxX,DT2,DT1)
        
        if (.not. Me%NewFormulation) then
        
            CHUNK = ChunkI !CHUNK_I(ILB, IUB)
                
k1:         do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
i1:         do i = ILB, IUB


                call ComputeAdvection1D_V2(JLB+1, JUB+1, Me%ExtVar%DT,                  &
                                        Me%ExtVar%DUX               (i,:),              &
                                        CurrProp%Concentration      (i,:,k),            &
                                        Me%ExtVar%FluxU             (i,:,k),            &
                                        Me%WaterVolume              (i,:,k),            & 
                                        Me%ExtVar%OpenPoints3D      (i,:,k),            &
                                        Me%COEF3_HorAdvXX%C_flux    (i,:,k),            &
                                        Me%COEF3_HorAdvXX%D_flux    (i,:,k),            &
                                        Me%COEF3_HorAdvXX%E_flux    (i,:,k),            &
                                        Me%COEF3_HorAdvXX%F_flux    (i,:,k),            &
                                        CurrProp%Evolution%AdvDiff%AdvMethodH,          &
                                        CurrProp%Evolution%AdvDiff%TVDLimitationH,      &
                                        CurrProp%Evolution%AdvDiff%VolumeRelMax,        &
                                        CurrProp%Evolution%AdvDiff%Upwind2H)

            end do i1
            !$OMP END DO
            end do k1
            
            !! $OMP END PARALLEL
        
        endif

        CHUNK = ChunkK !CHUNK_K(KLB, KUB)

cd6:    if (ImpExp_AdvXX == ExplicitScheme)  then !ExplicitScheme = 0

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
dok3 :      do k = KLB, KUB
doj3 :      do j = JLB, JUB
doi3 :      do i = ILB, IUB

            if (Me%ExtVar%ComputeFacesU3D(i, j  , k) == 1) then

                AdvFluxX =    (Me%COEF3_HorAdvXX%C_flux(i,   j, k)                          &
                            *  CurrProp%Concentration  (i, j-2, k)                          &
                            +  Me%COEF3_HorAdvXX%D_flux(i,   j, k)                          &
                            *  CurrProp%Concentration  (i, j-1, k)                          &
                            +  Me%COEF3_HorAdvXX%E_flux(i,   j, k)                          &
                            *  CurrProp%Concentration  (i,   j, k)                          &
                            +  Me%COEF3_HorAdvXX%F_flux(i,   j, k)                          &
                            *  CurrProp%Concentration  (i, j+1, k))

                Me%TICOEF3(i,j  ,k) = Me%TICOEF3(i,j  ,k) + AdvFluxX * Me%ExtVar%DT / Me%WaterVolume(i,j  ,k)
                Me%TICOEF3(i,j-1,k) = Me%TICOEF3(i,j-1,k) - AdvFluxX * Me%ExtVar%DT / Me%WaterVolume(i,j-1,k)

            endif

            end do doi3
            end do doj3
            end do dok3
            !$OMP END DO NOWAIT 

        else if (ImpExp_AdvXX == ImplicitScheme) then cd6 !ImplicitScheme = 1

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
dok4 :      do k = KLB, KUB
doj4 :      do j = JLB, JUB
doi4 :      do i = ILB, IUB


            if (Me%ExtVar%ComputeFacesU3D(i, j  , k) == 1) then

                DT2 = Me%ExtVar%DT / Me%WaterVolume(i,j  ,k)
                DT1 = Me%ExtVar%DT / Me%WaterVolume(i,j-1,k)

                Me%COEF3%D(i,j  ,k) = Me%COEF3%D(i,j  ,k) - Me%COEF3_HorAdvXX%D_flux(i,   j, k) * DT2
                Me%COEF3%E(i,j  ,k) = Me%COEF3%E(i,j  ,k) - Me%COEF3_HorAdvXX%E_flux(i,   j, k) * DT2

                Me%COEF3%E(i,j-1,k) = Me%COEF3%E(i,j-1,k) + Me%COEF3_HorAdvXX%D_flux(i,   j, k) * DT1
                Me%COEF3%F(i,j-1,k) = Me%COEF3%F(i,j-1,k) + Me%COEF3_HorAdvXX%E_flux(i,   j, k) * DT1


            endif


            end do doi4
            end do doj4
            end do dok4
            !$OMP END DO NOWAIT 

        else cd6

            stop 'sub. ModulePorousMediaProperties - HorizontalAdvectionXX - ERR01'
        
        endif cd6

        !$OMP END PARALLEL
        
        if (Me%Output%Boxes_ON .and. ImpExp_AdvXX == ExplicitScheme) call CalcHorizontalAdvFluxXX(CurrProp, ImpExp_AdvXX)

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "HorizontalAdvectionXX")


    end subroutine HorizontalAdvectionXX

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvectionYY(CurrProp, ImpExp_AdvYY)    !Routine from ModuleAdvectionDiffusion

        !Arguments--------------------------------------------------------------
        real                                :: ImpExp_AdvYY
        type (T_Property), pointer          :: CurrProp

        !Local-----------------------------------------------------------------               
        real(8)                             :: AdvFluxY, DT1, DT2
        integer                             :: i,     j,     k                             
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "HorizontalAdvectionYY")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
        
        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxY,DT2,DT1)
        
        if (.not. Me%NewFormulation) then
        
            CHUNK = ChunkJ !CHUNK_J(JLB, JUB)

k1:         do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
j1:         do j = JLB, JUB


                call ComputeAdvection1D_V2(ILB+1, IUB+1, Me%ExtVar%DT,                              &
                                        Me%ExtVar%DVY               (:,j),                          &
                                        CurrProp%Concentration      (:,j,k),                        &
                                        Me%ExtVar%FluxV             (:,j,k),                        &
                                        Me%WaterVolume              (:,j,k),                        & 
                                        Me%ExtVar%OpenPoints3D      (:,j,k),                        &
                                        Me%COEF3_HorAdvYY%C_flux    (:,j,k),                        &
                                        Me%COEF3_HorAdvYY%D_flux    (:,j,k),                        &
                                        Me%COEF3_HorAdvYY%E_flux    (:,j,k),                        &
                                        Me%COEF3_HorAdvYY%F_flux    (:,j,k),                        &
                                        CurrProp%Evolution%AdvDiff%AdvMethodH,                      &
                                        CurrProp%Evolution%AdvDiff%TVDLimitationH,                  &
                                        CurrProp%Evolution%AdvDiff%VolumeRelMax,                    &
                                        CurrProp%Evolution%AdvDiff%Upwind2H)

            end do j1
            !$OMP END DO
            end do k1
            
            !! $OMP END DO NOWAIT
            !! $OMP END PARALLEL
        
        endif
        
        CHUNK = ChunkK !CHUNK_K(KLB, KUB)
        
cd6:    if (ImpExp_AdvYY == ExplicitScheme)  then !ExplicitScheme = 0
           
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
dok3 :      do k = KLB, KUB
doj3 :      do j = JLB, JUB
doi3 :      do i = ILB, IUB

            if (Me%ExtVar%ComputeFacesV3D(i, j  , k) == 1) then

                AdvFluxY =    (Me%COEF3_HorAdvYY%C_flux(  i, j, k)                          &
                            *  CurrProp%Concentration  (i-2, j, k)                          &
                            +  Me%COEF3_HorAdvYY%D_flux(  i, j, k)                          &
                            *  CurrProp%Concentration  (i-1, j, k)                          &
                            +  Me%COEF3_HorAdvYY%E_flux(  i, j, k)                          &
                            *  CurrProp%Concentration  (  i, j, k)                          &
                            +  Me%COEF3_HorAdvYY%F_flux(  i, j, k)                          &
                            *  CurrProp%Concentration  (i+1, j, k))

                Me%TICOEF3(i  ,j,k) = Me%TICOEF3(i  ,j,k) + AdvFluxY * Me%ExtVar%DT / Me%WaterVolume(i  ,j,k)
                Me%TICOEF3(i-1,j,k) = Me%TICOEF3(i-1,j,k) - AdvFluxY * Me%ExtVar%DT / Me%WaterVolume(i-1,j,k)

            endif

            end do doi3
            end do doj3
            end do dok3
            !$OMP END DO NOWAIT

        else if (ImpExp_AdvYY == ImplicitScheme) then cd6 !ImplicitScheme = 1

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
dok4 :      do k = KLB, KUB
doj4 :      do j = JLB, JUB
doi4 :      do i = ILB, IUB


            if (Me%ExtVar%ComputeFacesV3D(i, j  , k) == 1) then

                DT2 = Me%ExtVar%DT / Me%WaterVolume(i  ,j  ,k)
                DT1 = Me%ExtVar%DT / Me%WaterVolume(i-1,j  ,k)

                Me%COEF3%D(i,j  ,k) = Me%COEF3%D(i,j  ,k) - Me%COEF3_HorAdvYY%D_flux(i,   j, k) * DT2
                Me%COEF3%E(i,j  ,k) = Me%COEF3%E(i,j  ,k) - Me%COEF3_HorAdvYY%E_flux(i,   j, k) * DT2

                Me%COEF3%E(i-1,j,k) = Me%COEF3%E(i-1,j,k) + Me%COEF3_HorAdvYY%D_flux(i,   j, k) * DT1
                Me%COEF3%F(i-1,j,k) = Me%COEF3%F(i-1,j,k) + Me%COEF3_HorAdvYY%E_flux(i,   j, k) * DT1

            endif


            end do doi4
            end do doj4
            end do dok4
            !$OMP END DO NOWAIT

        else cd6

            stop 'sub. HorizontalAdvectionYY - ModuleAdvectionDiffusion - ERR01'
        
        endif cd6

        !$OMP END PARALLEL
        
        if (Me%Output%Boxes_ON .and. ImpExp_AdvYY == ExplicitScheme) call CalcHorizontalAdvFluxXX(CurrProp, ImpExp_AdvYY)
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "HorizontalAdvectionYY")

    end subroutine HorizontalAdvectionYY
    
    !--------------------------------------------------------------------------

    subroutine VerticalDiffusion(CurrProperty, DiffusionV_Imp_Exp) !Routine from ModuleAdvectionDiffusion

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer           :: CurrProperty
        real                                :: DiffusionV_Imp_Exp
        !Local-----------------------------------------------------------------               
        
        real(8) :: AuxK, Aux1, Aux2

        integer :: i, j, k

        integer :: CHUNK

        !----------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "VerticalDiffusion")

        
        if (DiffusionV_Imp_Exp == ExplicitScheme) then
        
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,AuxK,Aux1,Aux2)

do2 :       do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%ComputeFacesW3D(i, j, k  ) == 1)  then
                
                    ! [m^2/s * m * m / m]
                    AuxK  =  CurrProperty%Diffusivity(i,j,k  ) &
                           * Me%ExtVar%Area          (i,j    ) &
                           * Me%ThetaAtFaces%ThetaW  (i,j,k  ) &
                           / Me%ExtVar%DZZ           (i,j,k-1)

                    ![m^3/s * s / m^3]
                    Aux1 = AuxK * dble(Me%ExtVar%DT) / Me%WaterVolume(i, j, k-1)
                    Aux2 = AuxK * dble(Me%ExtVar%DT) / Me%WaterVolume(i, j, k  ) 

                    Me%TICOEF3(i,j,k-1) = Me%TICOEF3(i,j,k-1) + Aux1 *                                  &
                                         (CurrProperty%Concentration(i,j,k)-CurrProperty%Concentration(i,j,k-1))


                    Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k  ) - Aux2 *                                  &
                                         (CurrProperty%Concentration(i,j,k)-CurrProperty%Concentration(i,j,k-1))

                endif

            end do do1
            end do do3
            !$OMP END DO
            end do do2

            !$OMP END PARALLEL
            
            
        elseif (DiffusionV_Imp_Exp == ImplicitScheme) then

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,AuxK,Aux1,Aux2)

do5 :       do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do6 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do4 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%ComputeFacesW3D(i, j, k  ) == 1)  then
                
                    ! [m^2/s * m * m / m]
                    AuxK  =  CurrProperty%Diffusivity(i,j,k  ) &
                           * Me%ExtVar%Area          (i,j    ) &
                           * Me%ThetaAtFaces%ThetaW  (i,j,k  ) &
                           / Me%ExtVar%DZZ           (i,j,k-1)

                    ![m^3/s * s / m^3]
                    Aux1 = AuxK * dble(Me%ExtVar%DT) / Me%WaterVolume(i, j, k-1) 
                    Aux2 = AuxK * dble(Me%ExtVar%DT) / Me%WaterVolume(i, j, k  ) 

                    Me%COEF3%E(i,j,k-1) = Me%COEF3%E(i,j,k-1) + Aux1 

                    Me%COEF3%F(i,j,k-1) = Me%COEF3%F(i,j,k-1) - Aux1 



                    Me%COEF3%D(i,j,k  ) = Me%COEF3%D(i,j,k  ) - Aux2 

                    Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k  ) + Aux2 


                endif

            end do do4
            end do do6
            !$OMP END DO
            end do do5

            !$OMP END PARALLEL

        endif

        if (Me%Output%Boxes_ON .and. DiffusionV_Imp_Exp < 1.)                              &
            call CalcVerticalDifFlux(CurrProperty, 1. - DiffusionV_Imp_Exp)        
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "VerticalDiffusion")


    end subroutine VerticalDiffusion                                                          

    !--------------------------------------------------------------------------

    subroutine VerticalAdvection(CurrProp, AdvectionV_Imp_Exp) !Routine from ModuleAdvectionDiffusion

        !Argiments-------------------------------------------------------------
        type(T_Property), pointer           :: CurrProp
        real                                :: AdvectionV_Imp_Exp
        !Local-----------------------------------------------------------------               
        real(8)             :: AdvFluxZ, DT1, DT2
        integer             :: i, j, k  
        integer             :: CHUNK                           
        !----------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "VerticalAdvection")

       
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxZ,DT1,DT2)

        !CHUNK = ChunkJ !CHUNK_J(Me%Size%JLB, Me%Size%JUB)
    
        !! $OMP PARALLEL SHARED(CHUNK) PRIVATE(I,J)

        if (.not. Me%NewFormulation) then 

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
j1:         do j = Me%WorkSize%JLB, Me%WorkSize%JUB
i1:         do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%OpenPoints3D (i,j,Me%WorkSize%KUB) == 1) then
                    call ComputeAdvection1D_V2( Me%WorkSize%KLB+1,                          &
                                                Me%WorkSize%KUB+1,                          &
                                                Me%ExtVar%DT,                               &
                                                Me%ExtVar%DWZ               (i,j,:),        &
                                                CurrProp%Concentration      (i,j,:),        &
                                                Me%FluxWCorr                (i,j,:),        &
                                                Me%WaterVolume              (i,j,:),        & 
                                                Me%ExtVar%OpenPoints3D      (i,j,:),        &
                                                Me%COEF3_VertAdv%C_flux     (i,j,:),        &
                                                Me%COEF3_VertAdv%D_flux     (i,j,:),        &
                                                Me%COEF3_VertAdv%E_flux     (i,j,:),        &
                                                Me%COEF3_VertAdv%F_flux     (i,j,:),        &
                                                CurrProp%Evolution%AdvDiff%AdvMethodV,      &
                                                CurrProp%Evolution%AdvDiff%TVDLimitationV,  &
                                                CurrProp%Evolution%AdvDiff%VolumeRelMax,    &   
                                                CurrProp%Evolution%AdvDiff%Upwind2V)
                endif

            end do i1
            end do j1
            
            !$OMP END DO
            !! $OMP END PARALLEL
        
        endif

cd6:    if (AdvectionV_Imp_Exp == ExplicitScheme)  then !ExplicitScheme = 0

dok3 :      do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj3 :      do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi3 :      do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(i, j  , k) == 1) then
                
                AdvFluxZ =    (Me%COEF3_VertAdv%C_flux(i,   j,   k)                     &
                            *  CurrProp%Concentration (i,   j, k-2)                     &
                            +  Me%COEF3_VertAdv%D_flux(i,   j,   k)                     &
                            *  CurrProp%Concentration (i,   j, k-1)                     &
                            +  Me%COEF3_VertAdv%E_flux(i,   j,   k)                     &
                            *  CurrProp%Concentration (i,   j,   k)                     &
                            +  Me%COEF3_VertAdv%F_flux(i,   j,   k)                     &
                            *  CurrProp%Concentration (i,   j, k+1))

                Me%TICOEF3(i,j,k-1) = Me%TICOEF3(i,j,k-1) - AdvFluxZ *                  &
                                      Me%ExtVar%DT / Me%WaterVolume(i,j,k-1)
                Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k  ) + AdvFluxZ *                  &
                                      Me%ExtVar%DT / Me%WaterVolume(i,j,k  )


            endif

            end do doi3
            end do doj3
            !$OMP END DO
            end do dok3

        else if (AdvectionV_Imp_Exp == ImplicitScheme) then cd6 !ImplicitScheme = 1

dok4 :      do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj4 :      do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi4 :      do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%ComputeFacesW3D(i, j  , k) == 1) then

                DT1 = Me%ExtVar%DT / Me%WaterVolume(i,j,k-1)
                DT2 = Me%ExtVar%DT / Me%WaterVolume(i,j,k  )

                Me%COEF3%D(i,j,k  ) = Me%COEF3%D(i,j,k  ) - Me%COEF3_VertAdv%D_flux(i,   j, k) * DT2
                Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k  ) - Me%COEF3_VertAdv%E_flux(i,   j, k) * DT2

                Me%COEF3%E(i,j,k-1) = Me%COEF3%E(i,j,k-1) + Me%COEF3_VertAdv%D_flux(i,   j, k) * DT1
                Me%COEF3%F(i,j,k-1) = Me%COEF3%F(i,j,k-1) + Me%COEF3_VertAdv%E_flux(i,   j, k) * DT1

            endif

            end do doi4
            end do doj4
            !$OMP END DO           
            end do dok4

        else cd6

            stop 'sub. VerticalAdvection - ModuleAdvectionDiffusion - ERR01'
        
        endif cd6

        !$OMP END PARALLEL

        if (Me%Output%Boxes_ON .and. AdvectionV_Imp_Exp < 1.)    &
            call CalcVerticalAdvFlux(CurrProp, 1. - AdvectionV_Imp_Exp)
            
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "VerticalAdvection")


    end subroutine VerticalAdvection

    !--------------------------------------------------------------------------

    subroutine CalcVerticalAdvFlux(CurrProp, Weigth)


        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp
        real, intent(IN)                            :: Weigth !Refers to the wigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "CalcVerticalAdvFlux")

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(k,j,i)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

        if (Me%ExtVar%ComputeFacesW3D(i, j, k) == 1) then

            Me%Fluxes%AdvFluxZ(i, j, k) =                                               &
                          Me%Fluxes%AdvFluxZ(i, j, k)                                   &
                        + Weigth                                                        &
                        *(Me%COEF3_VertAdv%C_flux(i,j,k)                                &
                        * CurrProp%Concentration(i, j, k-2)                             &
                        + Me%COEF3_VertAdv%D_flux(i,j,k)                                &
                        * CurrProp%Concentration(i, j, k-1)                             &
                        + Me%COEF3_VertAdv%E_flux(i,j,k)                                &
                        * CurrProp%Concentration(i, j, k )                              &
                        + Me%COEF3_VertAdv%F_flux(i,j,k)                                &
                        * CurrProp%Concentration(i, j, k+1))

        endif

        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "CalcVerticalAdvFlux")

        !----------------------------------------------------------------------

    end subroutine CalcVerticalAdvFlux

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalAdvFluxXX(CurrProp, Weigth)

        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp    
        real, intent(IN) :: Weigth !Refers to the wigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------

        integer :: i,     j,     k  
        integer :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "CalcHorizontalAdvFluxXX")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(i,j,k)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExtVar%ComputeFacesU3D(i, j, k) == 1) then

            Me%Fluxes%AdvFluxX(i, j, k) =                                               &
                          Me%Fluxes%AdvFluxX(i, j, k)                                   &
                        + Weigth                                                        &
                        * (Me%COEF3_HorAdvXX%C_flux(i,   j, k)                          &
                        *  CurrProp%Concentration  (i, j-2, k)                          &
                        +  Me%COEF3_HorAdvXX%D_flux(i,   j, k)                          &
                        *  CurrProp%Concentration  (i, j-1, k)                          &
                        +  Me%COEF3_HorAdvXX%E_flux(i,   j, k)                          &
                        *  CurrProp%Concentration  (i,   j, k)                          &
                        +  Me%COEF3_HorAdvXX%F_flux(i,   j, k)                          &
                        *  CurrProp%Concentration  (i, j+1, k))

        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "CalcHorizontalAdvFluxXX")

        !----------------------------------------------------------------------

    end subroutine CalcHorizontalAdvFluxXX

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalAdvFluxYY(CurrProp, Weigth)

        !External--------------------------------------------------------------
    
        real, intent(IN) :: Weigth !Refers to the wigth of Implicit-Explicit calculations
        type (T_Property), pointer                  :: CurrProp    

        !Local-----------------------------------------------------------------

        integer :: i,     j,     k  
        integer :: CHUNK
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "CalcHorizontalAdvFluxYY")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,k)
        
dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExtVar%ComputeFacesV3D(i  , j, k) == 1) then
            Me%Fluxes%AdvFluxY(i, j, k) =                                               &
                          Me%Fluxes%AdvFluxY(i, j, k)                                   &
                        + Weigth                                                        &
                        * (Me%COEF3_HorAdvYY%C_flux(i  , j, k)                          &
                        *  CurrProp%Concentration  (i-2, j, k)                          &
                        +  Me%COEF3_HorAdvYY%D_flux(i,   j, k)                          &
                        *  CurrProp%Concentration  (i-1, j, k)                          &
                        +  Me%COEF3_HorAdvYY%E_flux(i,   j, k)                          &
                        *  CurrProp%Concentration  (i,   j, k)                          &
                        +  Me%COEF3_HorAdvYY%F_flux(i,   j, k)                          &
                        *  CurrProp%Concentration  (i+1, j, k))
        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "CalcHorizontalAdvFluxYY")

        !----------------------------------------------------------------------

    end subroutine CalcHorizontalAdvFluxYY

    !--------------------------------------------------------------------------

    subroutine CalcVerticalDifFlux(CurrProp, Weigth)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp           
        real, intent(IN)                            :: Weigth !Refers to the wigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "CalcVerticalDifFlux")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i, j, k)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExtVar%ComputeFacesW3D(i, j, k  ) == 1) then
            Me%Fluxes%DifFluxZ(i, j, k) =                                               &
                          Me%Fluxes%DifFluxZ(i, j, k)                                   &
                        - Weigth                                                        &
                        * CurrProp%Diffusivity(i,j,k)                                   &
                        * Me%ExtVar%DUX(i,j)                                            &
                        * Me%ExtVar%DVY(i,j)                                            &
                        / Me%ExtVar%DZZ(i,j,k-1)                                        &
                        *(CurrProp%Concentration(i, j, k  )                             &
                        - CurrProp%Concentration(i, j, k-1))
        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "CalcVerticalDifFlux")

        !----------------------------------------------------------------------

    end subroutine CalcVerticalDifFlux

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalDifFluxXX(CurrProp)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp    

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "CalcHorizontalDifFluxXX")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(i,j,k)
        
dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%ComputeFacesU3D(i, j  , k) == 1) then                           
                Me%Fluxes%DifFluxX(i, j, k) =                                           &
                              Me%Fluxes%DifFluxX      (i,  j, k)                        &
                            - CurrProp%ViscosityU     (i,  j, k)                        &
                            * Me%ExtVar%AreaU         (i,  j, k)                        &
                            / Me%ExtVar%DZX           (i,j-1   )                        &
                            *(CurrProp%Concentration  (i,  j, k)                        &
                            - CurrProp%Concentration  (i,j-1, k))

            endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "CalcHorizontalDifFluxXX")


    end subroutine CalcHorizontalDifFluxXX

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalDifFluxYY(CurrProp)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp    

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "CalcHorizontalDifFluxYY")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,k)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExtVar%ComputeFacesV3D(i  , j, k) == 1) then
            Me%Fluxes%DifFluxY(i, j, k) =                                               &
                          Me%Fluxes%DifFluxY      (i  , j, k)                           &
                        - CurrProp%ViscosityV     (i  , j, k)                           &
                        * Me%ExtVar%AreaV         (i  , j, k)                           &
                        / Me%ExtVar%DZY           (i-1, j   )                           &
                        *(CurrProp%Concentration  (i  , j, k)                           &
                        - CurrProp%Concentration  (i-1, j, k))
        endif
        end do doi1
        end do doj1
        !$OMP END DO
        end do dok1
        
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "CalcHorizontalDifFluxYY")
        

    end subroutine CalcHorizontalDifFluxYY

    !--------------------------------------------------------------------------
!    subroutine ModifyExplicitCoefs(PropertyX) 
!    
!        !Arguments-------------------------------------------------------------
!        type (T_Property), pointer                  :: PropertyX
!
!
!        !Local-----------------------------------------------------------------
!        type (T_Property), pointer                  :: CurrProperty
!        integer                                     :: i, j, k !, CHUNK
!        real(8)                                     :: Area_Top_W, Area_Bottom_W, Area_Top_U, Area_Bottom_U
!        real(8)                                     :: Area_Top_V, Area_Bottom_V
!        real(8)                                     :: AdvTermA_W, AdvTermA_U, AdvTermA_V
!        real(8)                                     :: AdvTermB_Top_W, AdvTermB_Top_U, AdvTermB_Top_V
!        real(8)                                     :: AdvTermB_Bottom_W, AdvTermB_Bottom_U, AdvTermB_Bottom_V
!        real(8)                                     :: AdvTermC_W, AdvTermC_U, AdvTermC_V        
!        real(8)                                     :: AdvTermD_W, AdvTermD_U, AdvTermD_V
!        real(8)                                     :: DifTerm_Top_W, DifTerm_Bottom_W, DifTerm_Top_U
!        real(8)                                     :: DifTerm_Bottom_U, DifTerm_Top_V, DifTerm_Bottom_V
!        real(8)                                     :: aux 
!        real(8), pointer, dimension(:,:,:)          :: FluxW, FluxU, FluxV
!        real   , pointer, dimension(:,:,:)          :: DWZ, DZZ, Theta, ThetaOld 
!        real   , pointer, dimension(:,:  )          :: DZX, DZY, DXX, DYY 
!        logical                                     :: ComputeCoefC_W, ComputeCoefInterfRunoff
!        !Begin-----------------------------------------------------------------
!   
!        !!CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
!        
!        CurrProperty => PropertyX
!        
!        !Change diffusivity
!        call ModifyDiffusivity_New(CurrProperty)
!        
!        FluxW     => Me%ExtVar%FluxW
!        FluxU     => Me%ExtVar%FluxU
!        FluxV     => Me%ExtVar%FluxV        
!        DWZ       => Me%ExtVar%DWZ
!        DZZ       => Me%ExtVar%DZZ
!        DZX       => Me%ExtVar%DZX
!        DZY       => Me%ExtVar%DZY
!        DXX       => Me%ExtVar%DXX
!        DYY       => Me%ExtVar%DYY
!        Theta     => Me%ExtVar%WaterContent
!        ThetaOld  => Me%ExtVar%WaterContentOld
!        
!        !!! $OMP PARALLEL PRIVATE(I,J,K)
!        !!! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
!                
!                Area_Top_W         = Me%ExtVar%Area(i, j)
!                Area_Bottom_W      = Me%ExtVar%Area(i, j)
!                Area_Top_U         = DYY(i,j+1) * DWZ(i,j,k)
!                Area_Bottom_U      = DYY(i,j  ) * DWZ(i,j,k)
!                Area_Top_V         = DXX(i+1,j) * DWZ(i,j,k)
!                Area_Bottom_V      = DXX(i,j  ) * DWZ(i,j,k)                
!              
!                !Auxuliar value for transport - units of flow^-1
!                !s/m3
!                aux             = (Me%ExtVar%DT/(Theta(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
!
!                !!FLUXES IN X, Y AND Z DIRECTION - including the infiltration flux
!                if (Me%AdvDiff_SpatialMethod==AdvDif_CentralDif_) then ! diferenças centrais
!                    
!                    !in top cell check if flux entering or exiting to define coefB_Top_W or coefD_W
!                    if (K == Me%WorkSize%KUB) then
!                        
!                        ComputeCoefInterfRunoff     = .true.
!                        ComputeCoefC_W              = .false.
!                        Me%COEFExpl%coefC_W(i,j,k)  = 0.0
!                        DZZ(i,j,k)                  = DWZ(i,j,k)/2. + Me%ExtVar%InfiltrationColumn(i,j) / 2.
!                        
!                        if (FluxW(i,j,k+1) .lt. 0.0) then !Negative - downwards, entering the soil 
!                            
!                            AdvTermB_Top_W = 0.0
!                            AdvTermD_W     = (aux * FluxW(i,j,k+1))
!                            
!                        else !Positive (upwards, exiting the soil) or zero
!                            
!                            AdvTermB_Top_W = (aux * FluxW(i,j,k+1))
!                            AdvTermD_W     = 0.0
!
!                        endif                      
!                    else 
!                        
!                        ComputeCoefInterfRunoff                 = .false.
!                        ComputeCoefC_W                          = .true.
!                        Me%COEFExpl%CoefInterfRunoff(i,j,k)     = 0.0
!                        AdvTermB_Top_W                          = (aux * FluxW(i,j,k+1) / 2.)
!                        AdvTermC_W                              = (aux * FluxW(i,j,k+1) / 2.)
!                        
!                    endif
!                    
!                    AdvTermA_W        = (aux * FluxW(i,j,k  ) / 2.)
!                    AdvTermA_U        = (aux * FluxU(i,j,k  ) / 2.)
!                    AdvTermA_V        = (aux * FluxV(i,j,k  ) / 2.)
!                   !AdvTermB_Top_W    = already defined
!                    AdvTermB_Bottom_W = (aux * FluxW(i,j,k  ) / 2.)
!                    AdvTermB_Top_U    = (aux * FluxU(i,j+1,k) / 2.)
!                    AdvTermB_Bottom_U = (aux * FluxU(i,j  ,k) / 2.)
!                    AdvTermB_Top_V    = (aux * FluxV(i+1,j,k) / 2.)
!                    AdvTermB_Bottom_V = (aux * FluxV(i  ,j,k) / 2.)
!                   !AdvTermC_W        = already defined
!                    AdvTermC_U        = (aux * FluxU(i,j+1,k) / 2.)
!                    AdvTermC_V        = (aux * FluxV(i+1,j,k) / 2.)
!                    
!                        
!                elseif (Me%AdvDiff_SpatialMethod==AdvDif_Upwind_) then ! upwind
!
!                    if (K == Me%WorkSize%KUB) then
!                        
!                        ComputeCoefInterfRunoff    = .true.
!                        ComputeCoefC_W             = .false.
!                        Me%COEFExpl%coefC_W(i,j,k) = 0.0
!                        DZZ(i,j,k)                 = DWZ(i,j,k)/2. + Me%ExtVar%InfiltrationColumn(i,j) / 2.
!                    
!                    else 
!                        
!                        ComputeCoefInterfRunoff                 = .false.
!                        ComputeCoefC_W                          = .true.
!                        Me%COEFExpl%CoefInterfRunoff(i,j,k)     = 0.0
!                    
!                    endif
!                    
!                    !DirecW face k
!                    if (FluxW(i,j,k) .lt. 0.0) then !Bottom face, Negative - exiting
!                        AdvTermA_W        = 0.0
!                        AdvTermB_Bottom_W = aux * FluxW(i,j,k)
!                    else !Positive - entering or zero.
!                        AdvTermA_W        = aux * FluxW(i,j,k)
!                        AdvTermB_Bottom_W = 0.0
!                    endif
!
!                    !DirecU face j
!                    if (FluxU(i,j,k) .lt. 0.0) then !Left face, Negative - exiting
!                        AdvTermA_U        = 0.0
!                        AdvTermB_Bottom_U = aux * FluxU(i,j,k)
!                    else !Positive - entering or zero.
!                        AdvTermA_U        = aux * FluxU(i,j,k)
!                        AdvTermB_Bottom_U = 0.0
!                    endif
!
!                    !DirecV face i
!                    if (FluxV(i,j,k) .lt. 0.0) then !Left face, Negative - exiting
!                        AdvTermA_V        = 0.0
!                        AdvTermB_Bottom_V = aux * FluxV(i,j,k)
!                    else !Positive - entering or zero.
!                        AdvTermA_V        = aux * FluxV(i,j,k)
!                        AdvTermB_Bottom_V = 0.0
!                    endif
!                    
!                    !DirecW face k+1
!                    if (FluxW(i,j,k+1) .lt. 0.0) then !Top face, Negative - entering
!                        AdvTermC_W        = aux * FluxW(i,j,k+1)
!                        AdvTermD_W        = aux * FluxW(i,j,k+1)
!                        AdvTermB_Top_W    = 0.0
!                    else !Positive - exiting or zero.
!                        AdvTermC_W        = 0.0
!                        AdvTermD_W        = 0.0
!                        AdvTermB_Top_W    = aux * FluxW(i,j,k+1)
!                    endif
!
!                    !DirecU face j+1
!                    if (FluxU(i,j+1,k) .lt. 0.0) then !Right face, Negative - entering
!                        AdvTermC_U        = aux * FluxU(i,j+1,k)
!                        AdvTermD_U        = aux * FluxU(i,j+1,k)
!                        AdvTermB_Top_U    = 0.0
!                    else !Positive - exiting or zero.
!                        AdvTermC_U        = 0.0
!                        AdvTermD_U        = 0.0
!                        AdvTermB_Top_U    = aux * FluxU(i,j+1,k)
!                    endif                    
!
!                    !DirecV face i+1
!                    if (FluxV(i+1,j,k) .lt. 0.0) then !Right face, Negative - entering
!                        AdvTermC_V        = aux * FluxV(i+1,j,k)
!                        AdvTermD_V        = aux * FluxV(i+1,j,k)
!                        AdvTermB_Top_V    = 0.0
!                    else !Positive - exiting or zero.
!                        AdvTermC_V        = 0.0
!                        AdvTermD_V        = 0.0
!                        AdvTermB_Top_V    = aux * FluxV(i+1,j,k)
!                    endif                    
!                        
!
!                endif
!                
!                DifTerm_Top_W    = CurrProperty%Diffusivity(i,j,k+1) * Area_Top_W    * aux / DZZ(i,j,k  )
!                DifTerm_Bottom_W = CurrProperty%Diffusivity(i,j,k  ) * Area_Bottom_W * aux / DZZ(i,j,k-1)
!                DifTerm_Top_U    = CurrProperty%ViscosityU(i,j+1,k)  * Area_Top_U    * aux / DZX(i,j  )
!                DifTerm_Bottom_U = CurrProperty%ViscosityU(i,j,k  )  * Area_Bottom_U * aux / DZX(i,j-1)
!                DifTerm_Top_V    = CurrProperty%ViscosityV(i+1,j,k)  * Area_Top_V    * aux / DZY(i  ,j)
!                DifTerm_Bottom_V = CurrProperty%ViscosityV(i,j,k  )  * Area_Bottom_V * aux / DZY(i-1,j) 
!                
!                Me%COEFExpl%coefA_W(i,j,k) = AdvTermA_W                                                           &
!                                             + DifTerm_Bottom_W 
!
!                Me%COEFExpl%coefA_U(i,j,k) = AdvTermA_U                                                           &
!                                             + DifTerm_Bottom_U 
!                
!                Me%COEFExpl%coefA_V(i,j,k) = AdvTermA_V                                                           &
!                                             + DifTerm_Bottom_V 
!            
!                Me%COEFExpl%coefB_W(i,j,k) = - AdvTermB_Top_W                                                     &
!                                             + AdvTermB_Bottom_W                                                  &
!                                             - DifTerm_Bottom_W                                                   &
!                                             - DifTerm_Top_W
!
!                Me%COEFExpl%coefB_U(i,j,k) = - AdvTermB_Top_U                                                     &
!                                             + AdvTermB_Bottom_U                                                  &
!                                             - DifTerm_Bottom_U                                                   &
!                                             - DifTerm_Top_U
!                
!                Me%COEFExpl%coefB_V(i,j,k) = - AdvTermB_Top_V                                                     &
!                                             + AdvTermB_Bottom_V                                                  &
!                                             - DifTerm_Bottom_V                                                   &
!                                             - DifTerm_Top_V         
!                
!                if (ComputeCoefC_W) then    
!                    Me%COEFExpl%coefC_W(i,j,k) = - AdvTermC_W                                                     &
!                                                 + DifTerm_Top_W
!                endif
!                 
!                Me%COEFExpl%coefC_U(i,j,k) = - AdvTermC_U                                                         &
!                                             + DifTerm_Top_U                   
!
!                Me%COEFExpl%coefC_V(i,j,k) = - AdvTermC_V                                                         &
!                                             + DifTerm_Top_V    
!                if (ComputeCoefInterfRunoff) then    
!                    Me%COEFExpl%CoefInterfRunoff(i,j,k) = - AdvTermD_W                                            &
!                                                          + DifTerm_Top_W 
!                endif
!                
!            
!            endif
!        enddo
!        enddo
!        enddo
!        !!! $OMP END DO
!        !!! $OMP END PARALLEL
!
!
!    end subroutine ModifyExplicitCoefs
    
    !---------------------------------------------------------------------------

    subroutine ModifyEVTPCoefs()

        !Arguments----------------------------------------------------------------
        !type (T_Property), pointer                  :: PropertyX
        
        !Local--------------------------------------------------------------------
        real, pointer, dimension(:,:,:)             :: TranspirationLayer
        real(8), pointer, dimension(:,:  )          :: Evaporation
        logical                                     :: TranspirationExists, EvaporationExists
        integer                                     :: i,j,k, CHUNK, STAT_CALL
        real(8)                                     :: aux
        
        !Begin--------------------------------------------------------------------
        
        TranspirationExists = .false.
        EvaporationExists   = .false.
        
        !Effective transpiration used in PM (computed in ModuleVegetation) in [m3/s]
        call GetTranspiration(Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop ' ModifyEVTPCoefs - ModulePorousMediaProperties - ERR01'
        
        if (associated(TranspirationLayer)) TranspirationExists = .true.
        
        !Effective evaporation computed in PM in [m]
        call GetEvaporation(Me%ObjPorousMedia, Evaporation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop ' ModifyEVTPCoefs - ModulePorousMediaProperties - ERR010'
        
        if (associated(Evaporation)) EvaporationExists = .true.
        
        !Computation
        if (TranspirationExists) then
        
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB) 
          
            !$OMP PARALLEL PRIVATE(i,j,k,aux)
            
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if ((Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint)) then
                    
                    !Auxuliar value for transport - units of flow^-1
                    !s/m3
                    aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                    
                    ! Positive flow -> looses mass
                    Me%COEFExpl%CoefInterfTransp(i,j,k) = - aux * TranspirationLayer(i,j,k)
                
                endif
                
            enddo
            enddo
            !$OMP END DO
            enddo
            
            !$OMP END PARALLEL

        endif

        if (EvaporationExists) then
        
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB) 
          
            !$OMP PARALLEL PRIVATE(i,j,k,aux)
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                K = Me%WorkSize%KUB
                if ((Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint)) then
                    
                    !Auxuliar value for transport - units of flow^-1
                    !s/m3
                    aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                    
                    ! Positive flow -> looses mass
                    !evaporation flux [m3/s] = [m*m2/s]
                    Me%COEFExpl%CoefInterfEvap(i,j,k) =  - aux * (Evaporation(i,j) * Me%ExtVar%Area(i,j) / Me%ExtVar%DT)
                
                endif
                
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

        endif
        
        call UngetPorousMedia (Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyEVTPCoefs - ModulePorousMediaProperties - ERR050'

        call UngetPorousMedia (Me%ObjPorousMedia, Evaporation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyEVTPCoefs - ModulePorousMediaProperties - ERR060'

    
    end subroutine ModifyEVTPCoefs

    !-----------------------------------------------------------------------------

    subroutine ModifyDrainageNetworkCoefs (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k, STAT_CALL, CHUNK
        real(8)                                     :: aux 
        integer                                     :: GWFlowLink
        character (Len = 5)                         :: str_i, str_j, str_k
        character (Len = 15)                        :: str_conc
        character (len = StringLength)              :: string_to_be_written
        
        !Begin-----------------------------------------------------------------
  
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyDrainageNetworkCoefs")

   
        call GetGWFlowOption (Me%ObjPorousMedia, GWFlowLink, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop ' ModifyDrainageNetworkCoefs - ModulePorousMediaProperties - ERR01'
        
        CurrProperty => PropertyX
        
       !Flux between river and runoff in layers
       
        if (GWFlowLink == Layer_) then

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,str_i,str_j,str_k,str_conc,string_to_be_written,aux)
            
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if ((Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) .and. (Me%ExtVar%RiverPoints(I,J) == WaterPoint)) then   
                           
                    if((K .ge. Me%ExtVar%GWFlowBottomLayer(i,j)) .and. (K .le. Me%ExtVar%GWFlowTopLayer(i,j))) then
                        
                        !Auxuliar value for transport - units of flow^-1
                        !s/m3
                        aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                        
                        ! Positive flow -> looses mass
                        Me%COEFExpl%CoefInterfDN(i,j,k) = - aux * Me%ExtVar%FlowToChannelsLayer(i,j,k)

                      
                        ! mass going to channel -> conc from soil
                        if (Me%ExtVar%FlowToChannelsLayer(i,j,k) .gt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) =  CurrProperty%ConcentrationOld(i,j,k)                         
                            
                            if (CurrProperty%ConcentrationOld(i,j,k) < 0.0) then
                                write(str_i, '(i3)') i 
                                write(str_j, '(i3)') j
                                write(str_k, '(i3)') k
                                write(str_conc, '(ES10.3)') CurrProperty%ConcentrationOld(i,j,k)                               
                                string_to_be_written = trim(CurrProperty%ID%Name) //' in PMP is < 0 going from soil to river' &
                                                       //' in cell (i,j,k) with conc: '   &
                                                      //str_i//','//str_j//','//str_k//' '//str_conc
                                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)
                            endif   
                                                  
                        !mass coming from channel -> conc from DN
                        elseif (Me%ExtVar%FlowToChannelsLayer(i,j,k) .lt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) = CurrProperty%ConcentrationDN(i,j)

                            if (CurrProperty%ConcentrationDN(i,j) < 0.0) then
                                write(str_i, '(i3)') i 
                                write(str_j, '(i3)') j
                                write(str_k, '(i3)') k
                                write(str_conc, '(ES10.3)') CurrProperty%ConcentrationDN(i,j)                               
                                string_to_be_written = trim(CurrProperty%ID%Name)// ' in PMP is < 0 going from river to soil' &
                                                        //' in cell (i,j,k) with conc: '   &
                                                      //str_i//','//str_j//','//str_k//' '//str_conc
                                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

                            endif
                            
                        endif
                    
                    endif
               endif
            
            enddo
            enddo
            !$OMP END DO
            enddo
            
            !$OMP END PARALLEL
                           
        else
            
            !Flux between river and runoff one value for each soil column
            !water removed or added in top cell of saturated zone

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyDrainageNetworkCoefs - ModulePorousMediaProperties - ERR010'
            
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,str_i,str_j,str_k,str_conc,string_to_be_written,aux)
            
            !do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%BasinPoints(I,J) == WaterPoint .and. Me%ExtVar%RiverPoints(I,J) == WaterPoint) then   

                    K = Me%ExtVar%GWlayerOld(i,j)                
                    
                    if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                    
                        !Auxuliar value for transport - units of flow^-1
                        !s/m3
                        aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))

                        ! Positive flow -> looses mass
                        Me%COEFExpl%CoefInterfDN(i,j,k) = - aux * Me%ExtVar%FlowToChannels(i,j)

                        
                        ! mass going to channel -> conc from soil
                        if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) =  CurrProperty%ConcentrationOld(i,j,k)

                            if (CurrProperty%ConcentrationOld(i,j,k) < 0.0) then
                                write(str_i, '(i3)') i 
                                write(str_j, '(i3)') j
                                write(str_k, '(i3)') k
                                write(str_conc, '(ES10.3)') CurrProperty%ConcentrationOld(i,j,k)                               
                                string_to_be_written = trim(CurrProperty%ID%Name) //' in PMP is < 0 going from soil to river' &
                                                       //' in cell (i,j,k) with conc: '   &
                                                      //str_i//','//str_j//','//str_k//' '//str_conc
                                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)
                            endif  

                        !mass coming from channel -> conc from DN
                        elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) = CurrProperty%ConcentrationDN(i,j)

                            if (CurrProperty%ConcentrationDN(i,j) < 0.0) then
                                write(str_i, '(i3)') i 
                                write(str_j, '(i3)') j
                                write(str_k, '(i3)') k
                                write(str_conc, '(ES10.3)') CurrProperty%ConcentrationDN(i,j)                               
                                string_to_be_written = trim(CurrProperty%ID%Name)// ' in PMP is < 0 going from river to soil' &
                                                        //' in cell (i,j,k) with conc: '   &
                                                      //str_i//','//str_j//','//str_k//' '//str_conc
                                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

                            endif

                        endif
                    endif
                endif
            enddo
            enddo
            !$OMP END DO
            !enddo  
            !$OMP END PARALLEL
       
            call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyDrainageNetworkCoefs - ModulePorousMediaProperties - ERR020'  
       
       
        endif
   
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyDrainageNetworkCoefs")

   
   end subroutine ModifyDrainageNetworkCoefs

    !---------------------------------------------------------------------------

    subroutine ModifyInfiltrationCoefs (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        real(8), pointer, dimension(:,:,:)          :: FluxW
        integer                                     :: i, j, k, CHUNK
        real                                        :: aux, TopDiffusion, DZZ

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyInfiltrationCoefs")
        
       
        CurrProperty => PropertyX
        
        FluxW => Me%ExtVar%FluxW
       
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
      
        !$OMP PARALLEL PRIVATE(i,j,k,aux,TopDiffusion,DZZ)
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            K = Me%WorkSize%KUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then   
                
                !Auxuliar value for transport - units of flow-1
                !s/m3
                aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                
                DZZ = 0.5 * (Me%ExtVar%DWZ(i,j,k) + Me%ExtVar%WaterColumn(i,j))
                
                ! - = m2/s * m2 /m * s/m3
                TopDiffusion = CurrProperty%Diffusivity(i,j,k+1) &
                               * Me%ExtVar%Area(i,j)             &
                               * Me%ThetaAtFaces%ThetaW(i,j,k  ) &
                               * aux / DZZ
                
                !positive flow  -> exiting soil
                Me%COEFExpl%CoefInterfRunoff(i,j,k) = - aux * FluxW(i,j,k+1) + TopDiffusion
                
             
           endif
        
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyInfiltrationCoefs")

   
   end subroutine ModifyInfiltrationCoefs

    !---------------------------------------------------------------------------

!    subroutine ModifyVegetationCoefs (PropertyX)
!    
!        !Arguments-------------------------------------------------------------
!        type (T_Property), pointer                  :: PropertyX
!
!
!        !Local-----------------------------------------------------------------
!        type (T_Property), pointer                  :: CurrProperty
!        real, pointer, dimension(:,:,:)             :: TranspirationLayer
!        integer                                     :: i, j, k, STAT_CALL, CHUNK
!!        real                                        :: NitrogenUptake, PhosphorusUptake
!!        real                                        :: CorrectedTranspirationFlux
!        real                                        :: aux
!
!        !Begin-----------------------------------------------------------------
!
!        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyVegetationCoefs")
!
!        
!        !all dissolved properties with advection diffusion will be transported with the transpiration
!        !if in vegetation module is used swat uptake, then nitrate and phosphorus uptake has to be changed.
!        !see this in the future
!        
!        CurrProperty => PropertyX
!        
!        call GetTranspiration(Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop ' ModifyVegetationCoefs - ModulePorousMediaProperties - ERR01'
!
!        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
!      
!        !$OMP PARALLEL PRIVATE(i,j,k,aux)
!        
!        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
!        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)        
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!            
!            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then   
!                
!                if(K .ge. Me%ExtVar%TranspirationBottomLayer(i,j)) then
!                    
!                    !Auxuliar value for transport - units of flow-1
!                    !s/m3
!                    aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
!                    
!                    !Positive flow - looses mass
!                    Me%COEFExpl%CoefInterfTransp(i,j,k) = - aux * TranspirationLayer(i,j,k)
!                    
!                    
!!                    !Change nitrate and diss phosphorus flux because it may be disconnected from flow (if SWAT model original
!!                    ! formulation is used)
!!                    !in this case the transpiration flux * conc does not match the mass flux so here is changed
!!                    !Concentration old was the concentration used in module vegetation to the uptake so it is maintained
!!                    ! for consistency
!!                    if (Me%ExtVar%ModelNitrogen .and. CurrProperty%ID%IDNumber == Nitrate_) then
!!                        !      gN/s       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2 / s                     
!!                        NitrogenUptake = Me%ExtVar%NitrogenUptake(i,j,k) * 1e3 * Me%ExtVar%Area(i,j) * 1e-4      &
!!                                         / Me%ExtVar%VegetationDT
!!                        
!!                        !AssociatedFlux - g/s / g/m3 = g/s * m3/g = m3/s
!!                        CorrectedTranspirationFlux        = NitrogenUptake * CurrProperty%ConcentrationOld(i,j,k)
!!                        Me%COEFExpl%CoefInterfTransp(i,j,k) = - aux * CorrectedTranspirationFlux
!! 
!!                    elseif (Me%ExtVar%ModelPhosphorus .and. CurrProperty%ID%IDNumber == Inorganic_Phosphorus_) then
!!                        !      gN/s       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2 / s                     
!!                        PhosphorusUptake = Me%ExtVar%PhosphorusUptake(i,j,k) * 1e3 * Me%ExtVar%Area(i,j) * 1e-4    &
!!                                           / Me%ExtVar%VegetationDT
!!                        
!!                        !AssociatedFlux - g/s / g/m3 = g/s * m3/g = m3/s
!!                        CorrectedTranspirationFlux        = PhosphorusUptake * CurrProperty%ConcentrationOld(i,j,k)
!!                        Me%COEFExpl%CoefInterfTransp(i,j,k) = - aux * CorrectedTranspirationFlux
!!                            
!!                    endif
!
!                endif
!           endif
!        
!        enddo
!        enddo
!        !$OMP END DO
!        enddo
!        
!        !$OMP END PARALLEL        
!        
!        call UngetPorousMedia (Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ModifyVegetationCoefs - ModulePorousMediaProperties - ERR050'
!                           
!   
!        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyVegetationCoefs")
!
!   
!   end subroutine ModifyVegetationCoefs

    !---------------------------------------------------------------------------

    subroutine ModifyBoundaryCoefs(PropertyX)

        !Arguments----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        
        !Local--------------------------------------------------------------------
        integer                                     :: i,j,k, CHUNK, STAT_CALL
        real(8)                                     :: aux
        
        !Begin--------------------------------------------------------------------
        
        
        if (Me%ExtVar%BoundaryImposedWalls) then
        
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB) 
          
            !$OMP PARALLEL PRIVATE(i,j,k,aux)
            
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if ((Me%ExtVar%BoundaryCells(I,J) == WaterPoint)) then
                    
                    !Auxuliar value for transport - units of flow^-1
                    !s/m3
                    aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                    
                    ! Positive flow -> gains mass
                    Me%COEFExpl%CoefInterfBoundaryWalls(i,j,k) = aux * Me%ExtVar%BoundaryFluxWalls(i,j,k)
                    
                    !negative flow - exiting - Prop is from soil if explicit
                    if (Me%ExtVar%BoundaryFluxWalls(i,j,k) .le. 0.0) then
                        
                        PropertyX%ConcInBoundary(i,j,k) = PropertyX%ConcentrationOld(i,j,k)
                    
                    else !entering - Prop is imposed
                        
                        if (PropertyX%Evolution%Boundary%BoundaryCondition == ImposedValue_) then
                            
                            PropertyX%ConcInBoundary(i,j,k) = PropertyX%Evolution%Boundary%DefaultBoundary
                        
                        elseif (PropertyX%Evolution%Boundary%BoundaryCondition == NullGradient_) then
                            
                            PropertyX%ConcInBoundary(i,j,k) = PropertyX%ConcentrationOld(i,j,k)
                        
                        endif
                    endif
                    
                endif
                
            enddo
            enddo
            !$OMP END DO
            enddo
            
            !$OMP END PARALLEL
        
        endif


        if (Me%ExtVar%BoundaryImposedBottom) then

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyBoundaryCoefs - ModulePorousMediaProperties - ERR005'

            call GetGeometryKFloor(Me%ObjGeometry,                                          &
                                   Z    = Me%ExtVar%KFloor,                                 &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, INTERNAL_, "ModifyBoundaryCoefs - ModulePorousMediaProperties. ERR10")
        
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB) 
          
            !$OMP PARALLEL PRIVATE(i,j,k,aux)
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExtVar%BasinPoints (i,j) == 1) then
                
                    K = Me%ExtVar%KFloor(i,j)
                    
                    if ((Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint)) then
                        
                        !Auxuliar value for transport - units of flow^-1
                        !s/m3
                        aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                        
                        ! Is always negative flow -> losses mass
                        Me%COEFExpl%CoefInterfBoundaryBottom(i,j,k) = aux * Me%ExtVar%BoundaryFluxBottom(i,j)
                        
                        !it does not need concentration definiton. it is always negative so removes what is in soil
                                          
                    endif
                
                endif
                
            enddo
            enddo
            !$OMP END DO            
            !$OMP END PARALLEL

            call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%KFloor,       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                call SetError(FATAL_, INTERNAL_, "ModifyBoundaryCoefs - ModulePorousMediaProperties. ERR20")
        
            call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyBoundaryCoefs - ModulePorousMediaProperties - ERR030'         
        
        endif

        
    
    end subroutine ModifyBoundaryCoefs

    !-----------------------------------------------------------------------------

    subroutine ModifyPropertyValues(PropertyX)
        
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k, CHUNK
        real(8)                                     :: coefB, CoefInterfDN, CoefInterfTransp
        real(8)                                     :: CoefInterfRunoff, CoefInterfEvap
        real(8)                                     :: CoefInterfBoundaryW, CoefInterfBoundaryB
        real(8), pointer, dimension(:,:,:)          :: FluxW
        !Begin-----------------------------------------------------------------    

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyPropertyValues")

        CurrProperty => PropertyX
        
        if (Me%AdvDiff_Explicit) then

            CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,CoefInterfRunoff,CoefInterfDN,CoefInterfTransp, &
            !$OMP& CoefInterfEvap,CoefInterfBoundaryW,CoefInterfBoundaryB,coefB)
               
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)              
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then   
                    
                    CoefInterfRunoff    = Me%COEFExpl%CoefInterfRunoff(i,j,k)
                    CoefInterfDN        = Me%COEFExpl%CoefInterfDN(i,j,k)
                    CoefInterfTransp    = Me%COEFExpl%CoefInterfTransp(i,j,k)
                    CoefInterfEvap      = Me%COEFExpl%CoefInterfEvap(i,j,k)
                    CoefInterfBoundaryW = Me%COEFExpl%CoefInterfBoundaryWalls(i,j,k)
                    CoefInterfBoundaryB = Me%COEFExpl%CoefInterfBoundaryBottom(i,j,k)
                    
                    !CoefB = Me%ExtVar%WaterContentOld(i,j,k)/Me%ExtVar%WaterContent(i,j,k)
                    CoefB = Me%WaterContentBT(i,j,k)/Me%ExtVar%WaterContent(i,j,k)
                    
                    Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k) + CoefB * CurrProperty%ConcentrationOld(i,j,k  )                  &
                                                      + CoefInterfRunoff * CurrProperty%ConcentrationOnInfColumn(i,j)         &
                                                      + CoefInterfDN     * CurrProperty%ConcInInterfaceDN(i,j,k)              &
                                                      + CoefInterfTransp * CurrProperty%ConcentrationOld(i,j,k  )             &
                                                      + CoefInterfEvap   * CurrProperty%ConcentrationOld(i,j,k  )             &
                                                      + CoefInterfBoundaryW * CurrProperty%ConcInBoundary(i,j,k)              &
                                                      + CoefInterfBoundaryB * CurrProperty%ConcentrationOld(i,j,k)
                                                                     
                    CurrProperty%Concentration(i,j,k) = Me%TICOEF3(i,j,k  )

                endif
           enddo
           enddo
           enddo
           !$OMP END DO
           !$OMP END PARALLEL  
           
        else !implicit
            
            FluxW => Me%ExtVar%FluxW

            CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,CoefInterfRunoff,CoefInterfDN,CoefInterfTransp, &
            !$OMP& CoefInterfEvap,CoefInterfBoundaryW,CoefInterfBoundaryB,coefB)
               
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)                          
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then           
                    
                    CoefInterfRunoff    = Me%COEFExpl%CoefInterfRunoff(i,j,k)
                    CoefInterfDN        = Me%COEFExpl%CoefInterfDN(i,j,k)
                    CoefInterfTransp    = Me%COEFExpl%CoefInterfTransp(i,j,k)
                    CoefInterfEvap      = Me%COEFExpl%CoefInterfEvap(i,j,k)
                    CoefInterfBoundaryW = Me%COEFExpl%CoefInterfBoundaryWalls(i,j,k)
                    CoefInterfBoundaryB = Me%COEFExpl%CoefInterfBoundaryBottom(i,j,k)
                    
                    !Water exiting soil to runoff - Compute Conc Implicitly
                    if ((CoefInterfRunoff /= 0.0) .and. (FluxW(i,j,k+1) .gt. 0.0)) then
                        !Put flux in implicit
                        Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k) - CoefInterfRunoff
                        !And nullify in independent term
                        CoefInterfRunoff    = 0.0
                    endif
                    
                    !Water exiting soil to DN - Comput Conc Implicitly
                    if (CoefInterfDN .lt. 0.0) then
                        Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k) - CoefInterfDN
                        CoefInterfDN        = 0.0
                    endif
                    
                    !Water exiting soil to plants - Comput Conc Implicitly
                    if (CoefInterfTransp .lt. 0.0) then
                        Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k) - CoefInterfTransp
                        CoefInterfTransp    = 0.0
                    endif

                    !Water exiting soil to atmosphere by evaporation - Comput Conc Implicitly
                    if (CoefInterfEvap .lt. 0.0) then
                        Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k) - CoefInterfEvap
                        CoefInterfEvap      = 0.0
                    endif

                    !Water exiting soil to trough wall boundary - Comput Conc Implicitly
                    if (CoefInterfBoundaryW .lt. 0.0) then
                        Me%COEF3%E(i,j,k  )  = Me%COEF3%E(i,j,k) - CoefInterfBoundaryW
                        CoefInterfBoundaryW  = 0.0
                    endif
                    
                    !Water exiting soil to trough soil boundary (always) - Comput Conc Implicitly
                    if (CoefInterfBoundaryB .lt. 0.0) then
                        Me%COEF3%E(i,j,k  )  = Me%COEF3%E(i,j,k) - CoefInterfBoundaryB
                        CoefInterfBoundaryB  = 0.0
                    endif                    

                    !CoefB = Me%ExtVar%WaterContentOld(i,j,k)/Me%ExtVar%WaterContent(i,j,k)
                    CoefB = Me%WaterContentBT(i,j,k)/Me%ExtVar%WaterContent(i,j,k)
                    
                    Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k)                                                      &
                                        + coefB * CurrProperty%ConcentrationOld(i,j,k)                           &
                                        + CoefInterfRunoff * CurrProperty%ConcentrationOnInfColumn(i,j)          &
                                        + CoefInterfDN     * CurrProperty%ConcInInterfaceDN(i,j,k)               & 
                                        + CoefInterfTransp * CurrProperty%ConcentrationOld(i,j,k  )              &
                                        + CoefInterfEvap   * CurrProperty%ConcentrationOld(i,j,k  )              &
                                        + CoefInterfBoundaryW * CurrProperty%ConcInBoundary(i,j,k)               &
                                        + CoefInterfBoundaryB * CurrProperty%ConcentrationOld(i,j,k)
                                        

                endif
            enddo
            enddo
            enddo

           !$OMP END DO
           !$OMP END PARALLEL  
           
           !3D model or 1D vertical
           !griflet: old call
!            call THOMASZ(Me%WorkSize%ILB, Me%WorkSize%IUB,                              &
!                         Me%WorkSize%JLB, Me%WorkSize%JUB,                              &
!                         Me%WorkSize%KLB, Me%WorkSize%KUB,                              &
!                         Me%COEF3%D,                                                    &
!                         Me%COEF3%E,                                                    &
!                         Me%COEF3%F,                                                    &
!                         Me%TICOEF3,                                                    &
!                         CurrProperty%Concentration,                                    &
!                         Me%VECG,                                                       &
!                         Me%VECW)      
            !griflet: new call
            ! Use CUDA to solve Thomas
            call THOMASZ(Me%WorkSize%ILB, Me%WorkSize%IUB,                              &
                         Me%WorkSize%JLB, Me%WorkSize%JUB,                              &
                         Me%WorkSize%KLB, Me%WorkSize%KUB,                              &
                         Me%THOMAS,                                                     &
                         CurrProperty%Concentration                                     &
#ifdef _ENABLE_CUDA
                         , Me%ObjCuda,                                                  &
                         .FALSE.                                                        &
#endif _ENABLE_CUDA
                         )
                 
        endif

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyPropertyValues")

    
    end subroutine ModifyPropertyValues
    
    !---------------------------------------------------------------------------

    subroutine ModifyInterfaceMassFluxes(PropertyX) 

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k, STAT_CALL !, CHUNK
        integer                                     :: GWFlowLink
        real(8), pointer, dimension(:,:,:)          :: FluxW
!        real, pointer, dimension(:,:,:)             :: TranspirationLayer
        !Begin-----------------------------------------------------------------    
        
        CurrProperty => PropertyX
        
        !!Drainage network interface mass balance 
        if (Me%ExtVar%CoupledDN) then
            call GetGWFlowOption (Me%ObjPorousMedia, GWFlowLink, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop ' ModifyAdvectionDiffusion - ModulePorousMediaProperties - ERR01'
            
            CurrProperty => PropertyX
            
           !Flux between river and runoff in layers
           
            if (GWFlowLink == Layer_) then

                !!! $OMP PARALLEL PRIVATE(I,J,K)
                !!! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                    if ((Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) .and. (Me%ExtVar%RiverPoints(I,J) == BasinPoint)) then   
                               
                        if((K .ge. Me%ExtVar%GWFlowBottomLayer(i,j)) .and. (K .le. Me%ExtVar%GWFlowTopLayer(i,j))) then
                            
                            ! mass going to channel -> conc from soil
                            if (Me%ExtVar%FlowToChannelsLayer(i,j,k) .gt. 0.0) then
                                
                                !Global Mass Exchange
                                ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
                                CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
!                                  + (Me%ExtVar%FlowToChannelsLayer(i,j,k) * CurrProperty%ConcentrationOld(i,j,k)   &
                                  + (Me%ExtVar%FlowToChannelsLayer(i,j,k) * CurrProperty%Concentration(i,j,k)      &
                                     * 1e-3 * Me%ExtVar%DT)
                            
                            !mass coming from channel -> conc from DN
                            elseif (Me%ExtVar%FlowToChannelsLayer(i,j,k) .lt. 0.0) then
                                
                                !Global Mass Exchange
                                ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
                                CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
                                  + (Me%ExtVar%FlowToChannelsLayer(i,j,k) * CurrProperty%ConcentrationDN(i,j)      &
                                     * 1e-3 * Me%ExtVar%DT)  
                                
                            endif
                        
                        endif
                   endif
                
                enddo
                enddo
                enddo
                               
            else
                
                !Flux between river and runoff one value for each soil column
                !water removed or added in top cell of saturated zone

                call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyInterfaceMassFluxes - ModulePorousMediaProperties - ERR010'

                !!! $OMP PARALLEL PRIVATE(I,J,K)
                !!! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
                !do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%BasinPoints(I,J) == WaterPoint .and. Me%ExtVar%RiverPoints(I,J) == WaterPoint) then   
                    
                        K = Me%ExtVar%GWlayerOld(i,j)                    
                        
                        if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                            
                            ! mass going to channel -> conc from soil
                            if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
                                
                                !Global Mass Exchange
                                ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
                                CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
    !                                  + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%ConcentrationOld(i,j,k)          &
                                  + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%Concentration(i,j,k)             &
                                     * 1e-3 * Me%ExtVar%DT)
                            
                            !mass coming from channel -> conc from DN
                            elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
                                
                                !Global Mass Exchange
                                ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
                                CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
                                  + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%ConcentrationDN(i,j)             &
                                     * 1e-3 * Me%ExtVar%DT)  
                                
                            endif
                        endif
                    endif
                enddo
                enddo
                !enddo  

                call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyInterfaceMassFluxes - ModulePorousMediaProperties - ERR020' 

            endif
        endif
        
        FluxW     => Me%ExtVar%FluxW
        
        !! Infiltration mass balance
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB 
            K = Me%WorkSize%KUB                       
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then                      
                if (FluxW(i,j,k+1) .lt. 0.0) then !Negative - downwards, entering the soil 
                    !Global Mass Exchange
                    ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
                    CurrProperty%MB%RPExchangeMass =  CurrProperty%MB%RPExchangeMass                                &
                                                     - (FluxW(i,j,k+1) * CurrProperty%ConcentrationOnInfColumn(i,j) &
                                                     * 1e-3 * Me%ExtVar%DT)  
                
                else !Positive (upwards, exiting the soil) or zero
                    !Global Mass Exchange
                    ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
                    CurrProperty%MB%RPExchangeMass =  CurrProperty%MB%RPExchangeMass                          &
!                                                     - (FluxW(i,j,k+1) * CurrProperty%ConcentrationOld(i,j,k) &
                                                     - (FluxW(i,j,k+1) * CurrProperty%Concentration(i,j,k)    &
                                                     * 1e-3 * Me%ExtVar%DT)  
                endif
            endif
        enddo
        enddo
        
!        !!Vegetation mass balance
!        if (Me%ExtVar%CoupledVegetation .and. Me%ExtVar%ModelWater) then
!    
!            call GetTranspiration(Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop ' ModifyInterfaceMassFluxes - ModulePorousMediaProperties - ERR01'
!    
!    
!            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
!            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!                
!                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then   
!                           
!                    if(K .ge. Me%ExtVar%TranspirationBottomLayer(i,j) ) then
!                        
!                        !kg = [kg] + [m3/s] * [g/m3] * [1e-3kg/g] * [s] 
!                        CurrProperty%MB%TranspiredMass = CurrProperty%MB%TranspiredMass +                                   &
!                                                          (TranspirationLayer(i,j,k) * CurrProperty%Concentration(i,j,k)    &
!                                                            *  1E-3 * Me%ExtVar%DT)
!                        
!    
!                    endif
!               endif
!            
!            enddo
!            enddo
!            enddo
!
!            call UngetPorousMedia (Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ModifyInterfaceMassFluxes - ModulePorousMediaProperties - ERR050'
!
!        
!        endif

    
    end subroutine ModifyInterfaceMassFluxes
    
    !---------------------------------------------------------------------------

    subroutine ModifyDiffusivity_New(CurrProperty)
        
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        !Local-----------------------------------------------------------------
        integer                                     :: I,J,K, CHUNK
        real(8), pointer, dimension(:,:  )          :: WaterCol
        real   , pointer, dimension(:,:,:)          :: Porosity, UnsatW, UnsatU, UnsatV
        real                                        :: WaterContent_Face, Porosity_Face
        real                                        :: DiffCoef
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ModifyDiffusivity_New")

        
        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K,DiffCoef,WaterContent_Face,Porosity_Face)

        Porosity  => Me%ExtVar%ThetaS
        
        !diffusion on top boundary only occurs if there will be water at the end of the time step in runoff
        WaterCol  => Me%ExtVar%WaterColumn         
        UnsatW    => Me%ExtVar%UnsatW
        DiffCoef  = CurrProperty%Evolution%AdvDiff%Molecular_Diff_Coef
        
        UnsatU    => Me%ExtVar%UnsatU
        UnsatV    => Me%ExtVar%UnsatV
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                !Z direction
                
                WaterContent_Face = Me%ThetaAtFaces%ThetaW(i,j,k)                
                Porosity_Face     = min(Porosity(i,j,k),Porosity(i,j,k-1))
                
                !Only compute diffusivity in compute faces W (faces tath are not boundaries), else is zero
                !m2/s = m2/s + m/s * m
                CurrProperty%Diffusivity(i,j,k)  = Me%ExtVar%ComputeFacesW3D(i,j,k) *                                       &
                                                   (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face))                &
                                                   + (abs(UnsatW(i,j,k)) * Me%Disper_Longi%Field(i,j,k) / WaterContent_Face)
                
                
                !U direction
                
                WaterContent_Face = Me%ThetaAtFaces%ThetaU(i,j,k)
                Porosity_Face     = min(Porosity(i,j,k),Porosity(i,j-1,k))
                
                !Only compute diffusivity in compute faces U (faces tath are not boundaries), else is zero
                !m2/s = m2/s + m/s * m
                CurrProperty%ViscosityU(i,j,k)  = Me%ExtVar%ComputeFacesU3D(i,j,k) *                                         &
                                                  (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face))                  &
                                                  + (abs(UnsatU(i,j,k)) * Me%Disper_Trans%Field(i,j,k) / WaterContent_Face)
                
                
                !V direction
                
                WaterContent_Face = Me%ThetaAtFaces%ThetaV(i,j,k)
                Porosity_Face     = min(Porosity(i,j,k),Porosity(i-1,j,k))
                
                !Only compute diffusivity in compute faces V (faces tath are not boundaries), else is zero
                CurrProperty%ViscosityV(i,j,k)  = Me%ExtVar%ComputeFacesV3D(i,j,k) *                                         &
                                                  (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face))                  &
                                                  + (abs(UnsatV(i,j,k)) * Me%Disper_Trans%Field(i,j,k) / WaterContent_Face)
                                                     

            endif                                                                        
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            K = Me%WorkSize%KUB
            
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
            
                !UpperFace computation
                WaterContent_Face = Me%ThetaAtFaces%ThetaW(i,j,k) 
                Porosity_Face     = Porosity(i,j,k)
                
                !Only compute diffusivity in top face if there is water in the runoff at the end of the timestep
                if (WaterCol(i,j) .gt. 0.) then
                    !m2/s = m2/s + ((m3/s) /m2) * m
                    CurrProperty%Diffusivity(i,j,k+1)  = (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face)) 
                else
                    CurrProperty%Diffusivity(i,j,k+1)  = 0.0
                endif
            endif                                                                        
        enddo
        enddo

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ModifyDiffusivity_New")

    
    end subroutine ModifyDiffusivity_New

    !---------------------------------------------------------------------------

    subroutine SoilQualityProcesses

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property),          pointer     :: PropertyX
        type (T_Property),          pointer     :: SoilDryDensity, Salinity, pH
        type (T_Property),          pointer     :: IonicStrength, PhosphorusAdsortionIndex
        type (T_Property),          pointer     :: Oxygen
        character (Len = StringLength)          :: WarningString
        type(T_SedimentRate),       pointer     :: SedimentRateX
        real                                    :: DT
!        type (T_SoilRate),      pointer     :: SoilRateX
!        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB 
!        integer                                 :: i, j, k
        
        !Begin-----------------------------------------------------------------
        
!        WIUB = Me%WorkSize%IUB
!        WJUB = Me%WorkSize%JUB
!        WILB = Me%WorkSize%ILB
!        WJLB = Me%WorkSize%JLB
!        WKUB = Me%WorkSize%KUB
!        WKLB = Me%WorkSize%KLB
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "SoilQualityProcesses")

        call ComputeDissolvedToParticulate3D
        
        !Properties not modified by sediment quality (not state variables) but needed in argument
        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property soil dry density not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR00'
        endif

        call SearchProperty(Salinity, Salinity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property salinity not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR10'
        endif

        call SearchProperty(pH, pH_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property pH not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR20'
        endif

        call SearchProperty(IonicStrength, IonicStrength_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property ionic strength not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR30'

        endif

        call SearchProperty(PhosphorusAdsortionIndex, PhosphorusAdsortionIndex_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property phosphorus adsortion index not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR40'
        endif

        call SearchProperty(Oxygen, Oxygen_  , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property oxygen not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR45'
        endif


        call ComputeWindVelocity
 

        
        if (Me%ExtVar%Now .GE. Me%Coupled%SoilQuality_NextCompute) then
            
            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))
                
                !Change property units for Module Sediment Quality
                !Gases O2, CO2, CH4
                if ((PropertyX%ID%IDNumber == Oxygen_)        .OR.   &
                    (PropertyX%ID%IDNumber == CarbonDioxide_) .OR.   &
                    (PropertyX%ID%IDNumber == Methane_)            ) then
                    
                    WarningString = "Entering"
                    call ChangePropertyUnits (PropertyX, WarningString)
                endif
                
                !if forcing with oxygen send it to sediment quality
                if (.not. Me%Coupled%SedQualityOxygenForcing) then

                    call Modify_Interface(InterfaceID               = Me%ObjInterface,                         &
                                          PropertyID                = PropertyX%ID%IDNumber,                   &
                                          Concentration             = PropertyX%Concentration,                 &
                                          WaterPoints3D             = Me%ExtVar%WaterPoints3D,                 &
                                          OpenPoints3D              = Me%ExtVar%OpenPoints3D,                  &
                                          WaterPercentage           = Me%ExtVar%ThetaF,                        &
                                          DissolvedToParticulate3D  = Me%DissolvedToParticulate3D,             &
                                          SoilDryDensity            = SoilDryDensity%Concentration,            &
                                          Salinity                  = Salinity%Concentration,                  &
                                          pH                        = pH%Concentration,                        &
                                          IonicStrength             = IonicStrength%Concentration,             &
                                          PhosphorusAdsortionIndex  = PhosphorusAdsortionIndex%Concentration,  &
                                          WindVelocity              = Me%ExtVar%WindVelocity3D,                &
                                          STAT                      = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                                               &
                        stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR050'
                
                else
                    call Modify_Interface(InterfaceID               = Me%ObjInterface,                         &
                                          PropertyID                = PropertyX%ID%IDNumber,                   &
                                          Concentration             = PropertyX%Concentration,                 &
                                          WaterPoints3D             = Me%ExtVar%WaterPoints3D,                 &
                                          OpenPoints3D              = Me%ExtVar%OpenPoints3D,                  &
                                          WaterPercentage           = Me%ExtVar%ThetaF,                        &
                                          DissolvedToParticulate3D  = Me%DissolvedToParticulate3D,             &
                                          SoilDryDensity            = SoilDryDensity%Concentration,            &
                                          Salinity                  = Salinity%Concentration,                  &
                                          pH                        = pH%Concentration,                        &
                                          IonicStrength             = IonicStrength%Concentration,             &
                                          PhosphorusAdsortionIndex  = PhosphorusAdsortionIndex%Concentration,  &
                                          WindVelocity              = Me%ExtVar%WindVelocity3D,                &
                                          Oxygen                    = Oxygen%Concentration,                    &
                                          STAT                      = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                                               &
                        stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR055'                
                
                endif

                !Change property units from Module Sediment Quality
                !Gases O2, CO2, CH4
                if ((PropertyX%ID%IDNumber == Oxygen_)        .OR.   &
                    (PropertyX%ID%IDNumber == CarbonDioxide_) .OR.   &
                    (PropertyX%ID%IDNumber == Methane_)            ) then
                    
                    WarningString = "Exiting"
                    call ChangePropertyUnits (PropertyX, WarningString)
                endif

                PropertyX => PropertyX%Next
                

            end do
            
            Me%Coupled%SoilQuality_NextCompute = Me%Coupled%SoilQuality_NextCompute +       &
                                                     Me%Coupled%SoilQuality_DT

        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if (PropertyX%Evolution%SoilQuality) then

                !if DTInterval, only update at given time
                if (PropertyX%Evolution%DTIntervalAssociated) then
                    DT = PropertyX%Evolution%DTInterval
                else !update every time
                    PropertyX%Evolution%NextCompute = Me%ExtVar%Now
                    DT = Me%ExtVar%DT
                endif

                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    !Change property units for Module Sediment Quality
                    !Gases O2, CO2, CH4
                    if ((PropertyX%ID%IDNumber == Oxygen_)        .OR.   &
                        (PropertyX%ID%IDNumber == CarbonDioxide_) .OR.   &
                        (PropertyX%ID%IDNumber == Methane_)            ) then
                        
                        WarningString = "Entering"
                        call ChangePropertyUnits (PropertyX, WarningString)
                    endif

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%Concentration,          &
                                          WaterPoints3D = Me%ExtVar%WaterPoints3D,          &
                                          DTProp        = DT,                               &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'SoilQuality_Processes - ModulePorousMediaProperties - ERR060'
                    
                    !Change property units from Module Sediment Quality
                    !Gases O2, CO2, CH4
                    if ((PropertyX%ID%IDNumber == Oxygen_)        .OR.   &
                        (PropertyX%ID%IDNumber == CarbonDioxide_) .OR.   &
                        (PropertyX%ID%IDNumber == Methane_)            ) then
                        
                        WarningString = "Exiting"
                        call ChangePropertyUnits (PropertyX, WarningString)
                    endif
                    
                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        SedimentRateX => Me%FirstSedimentRate

        do while (associated(SedimentRateX))

            call GetRateFlux(InterfaceID    = Me%ObjInterface,                          &
                             FirstProp      = SedimentRateX%FirstProp%IDNumber,         &
                             SecondProp     = SedimentRateX%SecondProp%IDNumber,        &
                             RateFlux3D     = SedimentRateX%Field,                      &
                             WaterPoints3D  = Me%ExtVar%WaterPoints3D,                  &
                             STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'SedimentQuality_Processes - ModulePorousMediaProperties - ERR03'
            
            !mg/L to g/s - from SedimentQuality
            !rate units are always mg/L for dissolved or particulate ->  * [1g/1000mg] * [1000L/m3] * [m3] / [s]
            where (Me%ExtVar%WaterPoints3D == WaterPoint)                                                &
            SedimentRateX%Field = SedimentRateX%Field * Me%ExtVar%WaterContent * Me%ExtVar%Cellvolume    &
            / Me%Coupled%SoilQuality_DT            
           
            call BoxDif(Me%ObjBoxDif,                                       &
                        SedimentRateX%Field,                                &
                        "soil_"//trim(adjustl(SedimentRateX%ID%Name)),       &
                        Me%ExtVar%WaterPoints3D,                            &
                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                    &
                stop 'SedimentQuality_Processes - ModulePorousMediaProperties - ERR05'

            SedimentRateX => SedimentRateX%Next
        enddo

        if (Me%Coupled%MinConcentration)     call SetLimitsConcentration 
        if (Me%Coupled%WarnOnNegativeValues) call WarnOnNegativeValues   ('After Soil Quality')

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SoilQualityProcesses")

    end subroutine SoilQualityProcesses
    
    !-----------------------------------------------------------------------------    

    subroutine ChangePropertyUnits (Property, Condition)

        !Arguments-------------------------------------------------------------
        type (T_Property),          pointer     :: Property
        character (Len = StringLength)          :: Condition
        !External--------------------------------------------------------------
        integer                                      :: i, j, k, CHUNK
        
        !Begin-----------------------------------------------------------------
        
        !Convert Sediment Quality units for O2
        if (Property%ID%IDNumber == Oxygen_) then
            
            !from mg/L in PMP to mol/L in SedimentQuality
            if (Condition == "Entering") then
                
                CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                !$OMP PARALLEL PRIVATE(I,J,K)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k=Me%WorkSize%KLB, Me%WorkSize%KUB
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        
                        ! [mol/L] = [mg/L] * 1g/1000mg / [g/mol] 
                        !O2 molecular weight is 31.9989 g/mol  
                        Property%Concentration(i,j,k) = Property%Concentration(i,j,k) /(1000. * 31.9989)
                    endif
                
                enddo
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL
                
                
            !from mol/L in SedimentQuality to mg/L in PMP
            elseif (Condition == "Exiting") then

                CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                !$OMP PARALLEL PRIVATE(I,J,K)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k=Me%WorkSize%KLB, Me%WorkSize%KUB
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        
                        ! [mg/L] = [mol/L] * [g/mol] * 1000mg/g
                        !O2 molecular weight is 31.9989 g/mol  
                        Property%Concentration(i,j,k) = Property%Concentration(i,j,k) * 1000. * 31.9989
                    endif
                
                enddo
                enddo
                enddo

                !$OMP END DO
                !$OMP END PARALLEL
            
            endif
        
        elseif ((Property%ID%IDNumber == CarbonDioxide_) .or. (Property%ID%IDNumber == Methane_)) then
            
            !Sediment Quality concentrations are the same as solid phases (mg/kgsoil)
            !assuming conversion to dissolved using dissolved to particulate 

            !from mg/L in PMP to mg/kgsoil in SedimentQuality
            if (Condition == "Entering") then

                CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                !$OMP PARALLEL PRIVATE(I,J,K)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k=Me%WorkSize%KLB, Me%WorkSize%KUB
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        
                        ! [mg/kgsoil] = [mg/L] * [L/kgsoil] 
                        Property%Concentration(i,j,k) = Property%Concentration(i,j,k) * Me%DissolvedToParticulate3D(i,j,k)
                    endif
                
                enddo
                enddo
                enddo

                !$OMP END DO
                !$OMP END PARALLEL
            
            !from mg/kgsoil in SedimentQuality to mg/L in PMP
            elseif (Condition == "Exiting") then

                CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                !$OMP PARALLEL PRIVATE(I,J,K)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do k=Me%WorkSize%KLB, Me%WorkSize%KUB
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        
                        ! [mg/L] = [mg/kgsoil] / [L/kgsoil] 
                        Property%Concentration(i,j,k) = Property%Concentration(i,j,k) / Me%DissolvedToParticulate3D(i,j,k)

                    endif
                
                enddo
                enddo
                enddo

                !$OMP END DO
                !$OMP END PARALLEL
            
            endif
        
        endif
            
    end subroutine ChangePropertyUnits

    !-----------------------------------------------------------------------------  

    subroutine ComputeWindVelocity

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                      :: i, j, k
        
        !Begin-----------------------------------------------------------------

        call SetMatrixValue(Me%ExtVar%WindVelocity3D, Me%Size, 0.0, Me%ExtVar%WaterPoints3D)

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            k = Me%WorkSize%KUB
        
            if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                
                ! km/day                        =  m/s  * 1E-3km/m * 86400s/day 
                Me%ExtVar%WindVelocity3D(i,j,k) = Me%ExtVar%WindVelocity2D(i,j) * 1E-3 * 86400
            endif

        enddo
        enddo

    end subroutine ComputeWindVelocity

    !-----------------------------------------------------------------------------    

    subroutine Partition_Processes

        !Local----------------------------------------------------------------- 
        type (T_Property), pointer              :: PartPropX
        type (T_Property), pointer              :: PropertyX
!        type (T_Property), pointer              :: Salinity
!        type (T_Property), pointer              :: CohesiveSediment
        real                                    :: DT, MassTransfer
        real                                    :: DissolvedFraction
        real                                    :: ParticulateFraction
        real                                    :: TransferRate
!        real                                    :: ReferencePartitionCoef
!        real                                    :: PartitionCoef
!        real                                    :: b, SalinityConcentration
!        real                                    :: RefSedFactor
        integer                                 :: STAT_CALL
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k
        integer                                 :: Couple_ID


        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "Partition_Processes")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KLB = Me%WorkSize%KLB 
        KUB = Me%WorkSize%KUB 

        PropertyX => Me%FirstProperty

do0:    do while(associated(PropertyX))
cd0:        if(.not. PropertyX%Particulate .and. PropertyX%Evolution%Partitioning) then

                !days
                !if DTInterval, only update at given time
                if (PropertyX%Evolution%DTIntervalAssociated) then
                    DT = PropertyX%Evolution%DTInterval / 86400.
                else !update every time
                    PropertyX%Evolution%NextCompute = Me%ExtVar%Now
                    DT = Me%ExtVar%DT /86400.
                endif

cd1:            if(Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    Couple_ID = PropertyX%Evolution%Partition%Couple_ID

                    call Search_Property(PartPropX, PropertyXID = Couple_ID, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Partition_Processes - ModulePorousMediaProperties - ERR02'
                    
!                    if(PropertyX%Evolution%Partition%UseSedimentRefConc)then
!                    
!                        call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) then                                      
!                            stop 'Partition_Processes - ModulePorousMediaProperties - ERR03'
!                        endif
!                    endif
                    
    !                if (PropertyX%Evolution%Partition%SalinityEffect) then
    !
    !                    call Search_Property(Salinity, PropertyXID = Salinity_, STAT = STAT_CALL)
    !                    if (STAT_CALL /= SUCCESS_)                                      &
    !                        stop 'Partition_Processes - ModuleRunoffProperties - ERR03'
    !                endif
                    
                    
    do3:            do k = KLB, KUB
    do1:            do j = JLB, JUB
    do2:            do i = ILB, IUB
        
    cd2:                if ((Me%ExtVar%WaterPoints3D(i,j,k) == WaterPoint)) then

!                            if (PropertyX%Evolution%Partition%SalinityEffect) then
!
!                                ReferencePartitionCoef = PartPropX%Evolution%Partition%Fraction / &
!                                                         PropertyX%Evolution%Partition%Fraction
!
!                                b                      = PropertyX%Evolution%Partition%EmpiricCoef
!
!                                if (Salinity%Concentration(i, j, k) > 36.) then
!
!                                    SalinityConcentration = 0.036
!
!                                else
!
!                                    SalinityConcentration = Salinity%Concentration(i, j, k) / 1000.
!
!                                endif
!
!                                PartitionCoef       = (1. + SalinityConcentration)**b *  &
!                                                       ReferencePartitionCoef
!
!                                DissolvedFraction   = 1. / (1. + PartitionCoef)
!
!                                ParticulateFraction = PartitionCoef / (1. + PartitionCoef)
!
!
!                            else

                                DissolvedFraction   = PropertyX%Evolution%Partition%Fraction

                                ParticulateFraction = PartPropX%Evolution%Partition%Fraction

!                            endif

!                            if(PropertyX%Evolution%Partition%UseSedimentRefConc)then
!
!                                RefSedFactor = CohesiveSediment%Concentration(i,j,k) / &
!                                               PropertyX%Evolution%Partition%SedimentRefConc
!
!                                if(RefSedFactor < 1.)then
!
!                                    TransferRate = PropertyX%Evolution%Partition%Rate * RefSedFactor
!                                
!                                else
!
!                                    TransferRate = PropertyX%Evolution%Partition%Rate
!
!                                end if
!
!                            else

                                TransferRate = PropertyX%Evolution%Partition%Rate

!                            end if

                            ! [g/m3]       =          [s]         * [s^-1]        * [g/m3]
                            MassTransfer    =         DT * TransferRate *          &
                            (DissolvedFraction   * PartPropX%Concentration(i, j,k) -        &                  
                             ParticulateFraction * PropertyX%Concentration(i, j,k))

                            
                            if ((MassTransfer .gt. 0.0) .and. (MassTransfer .gt. PartPropX%Concentration(i, j,k))) then
                                MassTransfer = PartPropX%Concentration(i, j,k)
                            elseif ((MassTransfer .lt. 0.0) .and. (-MassTransfer .gt. PropertyX%Concentration(i, j,k))) then
                                MassTransfer = -PropertyX%Concentration(i, j,k)
                            endif
                            
                            PartPropX%Concentration(i, j,k) =                               &
                                               PartPropX%Concentration(i, j,k) - MassTransfer 

                            PropertyX%Concentration(i, j,k) =                               &
                                               PropertyX%Concentration(i, j,k) + MassTransfer

                        endif cd2
                    enddo do2
                    enddo do1
                    enddo do3
                endif cd1
            endif cd0

            PropertyX=>PropertyX%Next

        enddo do0


        nullify (PropertyX, PartPropX)
!        nullify (Salinity) 

        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "Partition_Processes")


    end subroutine Partition_Processes

    !--------------------------------------------------------------------------    

    subroutine SetLimitsConcentration    
!    subroutine SetLimitsConcentration (Message)

        !Arguments-------------------------------------------------------------
!        character(len=*)                            :: Message
        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k!, CHUNK
!        character(len=5)                            :: char_i, char_j, char_k
!        character(len=20)                           :: char_conc        
!        character (len = StringLength)              :: StrWarning        
        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "SetLimitsConcentration")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KLB = Me%WorkSize%KLB 
        KUB = Me%WorkSize%KUB 


        Property => Me%FirstProperty  

do1 :   do while (associated(Property))
cd1 :       if (Property%Evolution%MinConcentration) then
                
!                CHUNK = ChunkK !CHUNK_K(Me%Size%KLB, Me%Size%KUB)
                
!                !$OMP PARALLEL SHARED(CHUNK, Property) PRIVATE(I,J,K)
!                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)

                do k=Me%WorkSize%KLB, Me%WorkSize%KUB
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
    
                        if (Property%Concentration(i, j, k) < Property%MinValue) then
                            
                            ! mass created
                            Property%Mass_created(i, j, k) = Property%Mass_Created(i, j, k)   +  &
                                                   (Property%MinValue                -  &
                                                    Property%Concentration(i, j, k)) *  (Me%ExtVar%WaterContent(i,j,k) * &
                                                    Me%ExtVar%CellVolume (i, j, k))

                            Property%Concentration(i, j, k) = Property%MinValue

!                            write(char_i, '(i4)')i
!                            write(char_j, '(i4)')j
!                            write(char_k, '(i4)')k
!                            write(char_conc, '(f20.8)') Property%Concentration(i,j,k)
!
!                            StrWarning = trim(Property%ID%Name)//' was modified to its MinValue in cell(i,j,k)'// &
!                                                               char_i//','//char_j//','//char_k//' '//char_conc
!
!                            call SetError(WARNING_, INTERNAL_, StrWarning, OFF)
                            
                        endif

                    endif

                enddo
                enddo
                enddo
                
!                !$OMP END DO NOWAIT
!                !$OMP END PARALLEL
                
            endif cd1
                
        Property => Property%Next
        end do do1

        nullify(Property)

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SetLimitsConcentration")


    end subroutine SetLimitsConcentration

    !--------------------------------------------------------------------------
    
    subroutine WarnOnNegativeValues (Message)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Message
        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k!, CHUNK
        character(len=5)                            :: char_i, char_j, char_k
        character(len=15)                           :: char_conc        
        character (len = StringLength)              :: StrWarning
        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "WarnOnNegativeValues")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KLB = Me%WorkSize%KLB 
        KUB = Me%WorkSize%KUB 


        Property => Me%FirstProperty  

do1 :   do while (associated(Property))

cd1 :       if (Property%Evolution%WarnOnNegativeValues) then
                
!                CHUNK = ChunkK !CHUNK_K(Me%Size%KLB, Me%Size%KUB)
                
!                !$OMP PARALLEL SHARED(CHUNK, Property) PRIVATE(I,J,K,char_i,char_j,char_k,char_conc)
!                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)

                do k=Me%WorkSize%KLB, Me%WorkSize%KUB
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
    
                        if (Property%Concentration(i, j, k) < 0.0) then
                            
                            write(char_i, '(i4)')i
                            write(char_j, '(i4)')j
                            write(char_k, '(i4)')k
                            write(char_conc, '(ES10.3)') Property%Concentration(i,j,k) 

                            StrWarning = trim(Property%ID%Name)//' has a negative concentration in cell(i,j,k)'// &
                                                               char_i//','//char_j//','//char_k//' '//char_conc//' '//Message

                            call SetError(WARNING_, INTERNAL_, StrWarning, OFF)
                            
                        endif

                    endif

                enddo
                enddo
                enddo
                
!                !$OMP END DO NOWAIT
!                !$OMP END PARALLEL
                
            endif cd1
                
        Property => Property%Next
        end do do1

        nullify(Property)

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "WarnOnNegativeValues")


    end subroutine WarnOnNegativeValues

    !--------------------------------------------------------------------------    
    

    subroutine ComputeDissolvedToParticulate3D

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
         
        !Local----------------------------------------------------------------- 
        integer                                 :: i, j, k, CHUNK
        real                                    :: DT, InstantValue, ResidualValue
        type(T_Property), pointer               :: SoilDryDensity!, DrySedimentVolume
        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ComputeDissolvedToParticulate3D")


        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR01'

        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR02'

        !Conversion factor to pass from dissolved concentrations (mg/LH20) to particulate (mg/kgsed) - LH20/kgsed

        CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        !$OMP PARALLEL PRIVATE(I,J,K, InstantValue,ResidualValue)
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%ExtVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then

                ! [Lwater/kgsed] = [m3water/m3cell]*[m3cell]*[1000L/m3] / ([kgsed/m3cell] * [m3cell]) 
                InstantValue                       = Me%ExtVar%WaterContent(i,j,k) *  Me%ExtVar%CellVolume(i,j,k) * 1000. / &
                                                    (SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k))

                ResidualValue                      = Me%DissolvedToParticulate3D(i,j,k)

                Me%DissolvedToParticulate3D(i,j,k) = (ResidualValue * Me%ResidualTime +                                     &
                                                      InstantValue * DT) / (Me%ResidualTime + DT)
                                                       
            end if
        end do
        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL 
        
        Me%ResidualTime = Me%ResidualTime + DT
        
        nullify (SoilDryDensity)
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ComputeDissolvedToParticulate3D")


    end subroutine ComputeDissolvedToParticulate3D
    !--------------------------------------------------------------------------

#ifdef _PHREEQC_
    !--------------------------------------------------------------------------
    subroutine SoilChemistryProcesses

        !External--------------------------------------------------------------
        integer :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property), pointer                   :: PropertyX , pHProp, pEProp, TemperatureProp
        integer                                      :: I, J, K
        real             , pointer, dimension(:,:,:) :: CellTheta
        real             , pointer, dimension(:,:,:) :: CellThetaS
        real(8)          , pointer, dimension(:,:,:) :: CellVolume   
        integer                                      :: index     
        real                                         :: DT
        
        !Begin-----------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "SoilChemistryProcesses")

        CellTheta  => Me%ExtVar%WaterContent
        CellVolume => Me%ExtVar%CellVolume
        CellThetaS => Me%ExtVar%ThetaS
        
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                        
                !            m3                =          %         *        m3
                Me%CellWaterVolume(I, J, K) = CellTheta(I, J, K) * CellVolume(I, J, K)
                !          kg             =         kg/m3         *            m3
                Me%PhreeqC%CellWaterMass(I, J, K) = WaterReferenceDensity * Me%CellWaterVolume(I, J, K)
                
                if (associated(Me%PhreeqC%SoilDryDensity)) then
                    !          kg            =                      kg/m3                       *         m3
                    Me%PhreeqC%CellSoilMass(I, J, K) = Me%PhreeqC%SoilDryDensity%Concentration(I, J, K) * CellVolume(I, J, K)
                else
                    !          kg            = 1000 *         m3
                    Me%PhreeqC%CellSoilMass(I, J, K) = 1000 * CellVolume(I, J, K)
                end if
                        
            end if
                    
        end do
        end do
        end do
        
        if (Me%ExtVar%Now .GE. Me%Coupled%SoilChemistry_NextCompute) then
            
            PropertyX => Me%FirstProperty
            do while(associated(PropertyX))

                if (PropertyX%Evolution%SoilChemistry) then

                    select case (PropertyX%ID%IDNumber)
                    case (Temperature_)
                        TemperatureProp => PropertyX
                    case (pH_)
                        pHProp => PropertyX
                    case (pE_)
                        pEProp => PropertyX
                    case (SoilDryDensity_)
                    end select

                end if

                PropertyX => PropertyX%Next
                
            end do
            
            
            PropertyX => Me%FirstProperty
            do while(associated(PropertyX))

                if (PropertyX%Evolution%SoilChemistry) then

                    select case (PropertyX%ID%IDNumber)
                    case (Temperature_, pH_, pE_, SoilDryDensity_)
                        !Do nothing
                    case default
                    
                        do index = 1, Me%PhreeqC%NumberOfModels
                            call Modify_Interface(InterfaceID        = Me%PhreeqC%Models(index)%ObjInterface, &
                                                  PropertyID         = PropertyX%ID%IDNumber,                 &
                                                  WaterVolume        = Me%CellWaterVolume,                    &
                                                  WaterMass          = Me%PhreeqC%CellWaterMass,              &
                                                  SolidMass          = Me%PhreeqC%CellSoilMass,               &
                                                  Temperature        = TemperatureProp%Concentration,         &
                                                  pH                 = pHProp%Concentration,                  &
                                                  pE                 = pEProp%Concentration,                  &
                                                  Concentration      = PropertyX%Concentration,               &
                                                  WaterPoints3D      = Me%ExtVar%WaterPoints3D,               &
                                                  OpenPoints3D       = Me%ExtVar%OpenPoints3D,                &
                                                  PhreeqCID          = Me%PhreeqC%Filter,                     &
                                                  STAT               = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) &
                                stop 'SoilChemistryProcesses - ModulePorousMediaProperties - ERR010'                                
                        enddo
                        
                    end select

                end if
                                
                PropertyX => PropertyX%Next
                
            end do
            
            Me%Coupled%SoilChemistry_NextCompute = Me%Coupled%SoilChemistry_NextCompute + &
                                                   Me%Coupled%SoilChemistry_DT

        end if


        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if (PropertyX%Evolution%SoilChemistry) then

                !if DTInterval, only update at given time
                if (PropertyX%Evolution%DTIntervalAssociated) then
                    DT = PropertyX%Evolution%DTInterval
                else !update every time
                    PropertyX%Evolution%NextCompute = Me%ExtVar%Now
                    DT = Me%ExtVar%DT
                endif

                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then
                
                    select case (PropertyX%ID%IDNumber)
                    case (Temperature_, pH_, pE_, SoilDryDensity_)
                    case default
                    
                        do index = 1, Me%PhreeqC%NumberOfModels                        
                            call Modify_Interface(InterfaceID   = Me%PhreeqC%Models(index)%ObjInterface, &
                                                  PropertyID    = PropertyX%ID%IDNumber,                 &
                                                  Concentration = PropertyX%Concentration,               &
                                                  WaterPoints3D = Me%ExtVar%WaterPoints3D,               &
                                                  DTProp        = DT,                                    &
                                                  STAT          = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) &
                                stop 'SoilChemistryProcesses - ModulePorousMediaProperties - ERR020'
                        enddo
                        
                    end select

                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        if (Me%Coupled%MinConcentration)     call SetLimitsConcentration !('After Soil Chemistry')
        if (Me%Coupled%WarnOnNegativeValues) call WarnOnNegativeValues   ('After Soil Chemistry')
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SoilChemistryProcesses")
        
        !End-------------------------------------------------------------------
        
    end subroutine SoilChemistryProcesses   
    !-----------------------------------------------------------------------------    
#endif

    !-----------------------------------------------------------------------------    
    subroutine ChainReactionsProcesses
    
        !Local-------------------------------------------------------------------- 
        integer                            :: STAT
        real, dimension(:,:,:), pointer    :: SoilDensity 
        type (T_Property), pointer         :: PropertyX 
        logical                            :: UseSoilMass
        
        !Begin-----------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ChainReactionsProcesses")

        call Search_Property(PropertyX, SoilVolumetricDensity_, STAT) 
        if (STAT .NE. SUCCESS_) then
            UseSoilMass = .false.            
        else
            UseSoilMass = .true.
            SoilDensity => PropertyX%Concentration
        endif
        
        call SetChainReactionsProperties
        
        if (UseSoilMass) then
            call ModifyChainReactions (Me%ObjChainReactions,       &
                                       Me%ExtVar%WaterContent,     &
                                       Me%ExtVar%DT,               &                                   
                                       SoilDensity = SoilDensity,  &
                                       STAT = STAT)                                       
            if (STAT .NE. SUCCESS_) &
                stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR010'
        else
            call ModifyChainReactions (Me%ObjChainReactions,       &
                                       Me%ExtVar%WaterContent,     &                                       
                                       Me%ExtVar%DT,               &                                   
                                       STAT = STAT)                                       
            if (STAT .NE. SUCCESS_) &
                stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR020'
        endif        
                                    
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ChainReactionsProcesses")    
        !-------------------------------------------------------------------------
    
    end subroutine ChainReactionsProcesses
    !-----------------------------------------------------------------------------    
    

    !-----------------------------------------------------------------------------    
    subroutine SetChainReactionsProperties
    
        !Local-------------------------------------------------------------------- 
        type (T_Property), pointer         :: PropertyX
        integer, dimension(:), pointer     :: CRPropertiesList
        integer                            :: PropertiesCount
        integer                            :: STAT
        integer                            :: Index
        
        !Begin-----------------------------------------------------------------
        call GetCRPropertiesList(Me%ObjChainReactions, CRPropertiesList, PropertiesCount, STAT = STAT)
        if (STAT .NE. SUCCESS_) &
            stop 'SetChainReactionsProperties - ModulePororusMediaProperties - ERR010'
    
        do Index = 1, PropertiesCount
    
            call Search_Property(PropertyX, CRPropertiesList(Index), STAT)    
            if (STAT .NE. SUCCESS_) &                            
                stop 'SetChainReactionsProperties - ModulePororusMediaProperties - ERR020'
                    
            call SetCRPropertyConcentration (Me%ObjChainReactions,       &
                                             CRPropertiesList(Index),    &
                                             PropertyX%Concentration,    &
                                             STAT)
            if (STAT .NE. SUCCESS_) &                            
                stop 'SetChainReactionsProperties - ModulePororusMediaProperties - ERR030'
                                    
        end do                
        
        call UnGetChainReactions(Me%ObjChainReactions, CRPropertiesList, STAT)
        if (STAT .NE. SUCCESS_) &
            stop 'SetChainReactionsProperties - ModulePororusMediaProperties - ERR040'        
        !-------------------------------------------------------------------------
    
    end subroutine SetChainReactionsProperties
    !-----------------------------------------------------------------------------    
    
    subroutine DecayProcesses
    
        !Local--------------------------------------------------------------------
        type (T_Property), pointer                         :: Property
        !Begin--------------------------------------------------------------------
        
        
        Property => Me%FirstProperty  

do1 :   do while (associated(Property))    

            if (Property%Evolution%Decay) then
            
                if (Property%Evolution%DecayEquationCoupled) then
                    if (Property%Evolution%DecayEquation == Peyrard_) then
                        !Decay from Tolouse team adapted from Peyrard et al. 2011
                        call DecayPeyrard(Property)
                    !elseif (Property%%Evolution%DecayEquation == SomeEquation_) then
                    !    call DecaySomeEquation
                    endif
                else
                    if (Property%Evolution%DecayMass) then
                        call DecayFirstOrder_Mass(Property)
                    else
                        call DecayFirstOrder_Conc(Property)
                    endif
                endif
                              
            endif     

            Property => Property%Next
        end do do1   

        nullify(Property)
        
            
    end subroutine DecayProcesses

    !-----------------------------------------------------------------------------
    !Decay applied to mass
    subroutine DecayFirstOrder_Mass(Property)
        
        !Arguments----------------------------------------------------------------
        type (T_Property), pointer                         :: Property
        
        !Local--------------------------------------------------------------------
        type (T_Property), pointer                         :: SoilDryDensity
        real                                               :: DT
        real(8)                                            :: WaterVolume
        real(8)                                            :: OldMass
        real(8)                                            :: NewMass
        real(8)                                            :: MassSink, SoilWeight
        integer                                            :: i,j,k, CHUNK, STAT_CALL
        !Begin--------------------------------------------------------------------
        
        
        !days
        !if DTInterval, only update at given time
        if (Property%Evolution%DTIntervalAssociated) then
            DT = Property%Evolution%DTInterval / 86400.
        else !update every time
            Property%Evolution%NextCompute = Me%ExtVar%Now
            DT = Me%ExtVar%DT /86400.
        endif
                  
        if(Me%ExtVar%Now .GE. Property%Evolution%NextCompute) then            
    
            if (Check_Particulate_Property(Property%ID%IDNumber)) then

                call SearchProperty(SoilDryDensity, SoilDryDensity_, .false., STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'DecayFirstOrder - ModulePorousMediaProperties - ERR01'
                
                CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                !$OMP PARALLEL PRIVATE(I,J,K, SoilWeight, OldMass, MassSink, NewMass)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                        
                        !WQ process without volume change. only mass change
                        !kgsoil = kgsoil/m3cell * m3cell
                        SoilWeight = SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k)
                        
                        !g =  mg/kgsoil  * kgsoil * 1E-3 g/mg
                        OldMass = Property%Concentration(i,j,k) * SoilWeight * 1E-3

                        !Particulate mg/kgsoil.day-1 = mg/kgsoil * day-1
                        Property%PropertyDecay(i,j,k) = Property%Concentration(I,J,K) * Property%Evolution%DecayRate
                                                                       
                        !P = P0*exp(-kt)
                        !g
                        MassSink = min (OldMass - OldMass * exp(-Property%Evolution%DecayRate * DT),  OldMass)
                        
                        NewMass = OldMass - MassSink
                        
                        !mg/kgsoil = g * 1E3 mg/g / kgsoil 
                        Property%Concentration(I,J,K) = NewMass * 1E3 / SoilWeight 
                        
                    endif
                
                enddo
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL 
    
            else
   
                CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                !$OMP PARALLEL PRIVATE(I,J,K, WaterVolume, OldMass, MassSink, NewMass)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                        
                        !WQ process without volume change. only mass change
                        !m3 H20
                        WaterVolume = Me%ExtVar%WaterContent(i,j,k)* Me%ExtVar%CellVolume(i,j,k)
                        
                        !g  = g/m3 * m3
                        OldMass = Property%Concentration(I,J,K)  * WaterVolume

                        !Dissolved g/m3.day-1 = g/m3 * day-1
                        Property%PropertyDecay(i,j,k) = Property%Concentration(I,J,K) * Property%Evolution%DecayRate
                        
                        !P = P0*exp(-kt)
                        !g limit the sink to the existing one (old)
                        MassSink = min (OldMass - OldMass * exp(-Property%Evolution%DecayRate * DT),  OldMass)
                        
                        NewMass = OldMass - MassSink
                        
                        !g/m3 = g / m3 
                        Property%Concentration(I,J,K) = NewMass / WaterVolume
                        
                    endif
                
                enddo
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL                    
            
            endif
        
        endif
                      
            
    end subroutine DecayFirstOrder_Mass

    !-----------------------------------------------------------------------------    
    !Decay applied to Concentration. Organisms usually repond to concentration
    subroutine DecayFirstOrder_Conc(Property)
        
        !Arguments----------------------------------------------------------------
        type (T_Property), pointer                         :: Property
        
        !Local--------------------------------------------------------------------
        real                                               :: DT
        integer                                            :: i,j,k, CHUNK
        !Begin--------------------------------------------------------------------
        
        
        !days
        !if DTInterval, only update at given time
        if (Property%Evolution%DTIntervalAssociated) then
            DT = Property%Evolution%DTInterval / 86400.
        else !update every time
            Property%Evolution%NextCompute = Me%ExtVar%Now
            DT = Me%ExtVar%DT /86400.
        endif
                  
        if(Me%ExtVar%Now .GE. Property%Evolution%NextCompute) then            
    
            CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
            !$OMP PARALLEL PRIVATE(I,J,K)
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                    
                    !Particulate mg/kgsoil.day-1 = mg/kgsoil * day-1
                    !Dissolved g/m3.day-1 = g/m3 * day-1
                    Property%PropertyDecay(i,j,k) = Property%Concentration(i,j,k) * Property%Evolution%DecayRate
                    
                    !P = P0*exp(-kt). k is day-1
                    !mg/kgsoil or g/m3
                    Property%Concentration(i,j,k) = max(Property%Concentration(i,j,k)*exp(-Property%Evolution%DecayRate*DT),0.0)
                    
                endif
            
            enddo
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL 
    
        endif
                      
            
    end subroutine DecayFirstOrder_Conc

    !-----------------------------------------------------------------------------    
    !Tolouse team adapted equation from Peyrard et al. 2011 - Longitudinal transformation of nitrogen and carbon 
    !in the hyporheic zone of an N-rich stream: A combined modelling and field study. Physics and Chemistry of 
    !the Earth 36 (2011) 599–611. 
    !It is needed DOC, POC, and Nitrate
    !This property (for now used for NO3) and DOC and POC are used in the equation in uM for DOC and NO3 and mg.g-1 for POC
    !so conversions were needed from mg/L (=g/m3) and mg/kg
    subroutine DecayPeyrard(Property)
        
        !Arguments----------------------------------------------------------------
        type (T_Property), pointer                         :: Property
        
        !Local--------------------------------------------------------------------
        type (T_Property), pointer                         :: SoilDryDensity, DOC, POC, POCr, POCl
        real             , pointer, dimension(:,:,:)       :: Porosity, ThetaF, auxPOC
        real                                               :: DT
        real                                               :: AdjustedDensity, PropFactor
        real                                               :: AnaerobioseFactor, Decay
        real(8)                                            :: NO3Conc, DOCConc, POCConc, NO3HalfSaturation
        integer                                            :: i,j,k, CHUNK, STAT_CALL
        logical                                            :: POCAvailable = .true.
        !Begin--------------------------------------------------------------------
        
        !now is used Peyrard only for nitrate but in the future may be added other properties to transform
        if (Property%ID%IDNumber == Nitrate_) then
        
            !days
            !if DTInterval, only update at given time
            if (Property%Evolution%DTIntervalAssociated) then
                DT = Property%Evolution%DTInterval / 86400.
            else !update every time
                Property%Evolution%NextCompute = Me%ExtVar%Now
                DT = Me%ExtVar%DT /86400.
            endif
                      
            if(Me%ExtVar%Now .GE. Property%Evolution%NextCompute) then  
            
                call SearchProperty(SoilDryDensity, SoilDryDensity_, .false., STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'DecayPeyrard - ModulePorousMediaProperties - ERR01'

                call SearchProperty(DOC, DOC_, .false., STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'DecayPeyrard - ModulePorousMediaProperties - ERR010'
                
                POCAvailable = .true.
                
                call SearchProperty(POC, POC_, .false., STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) then
                                      
                    !try with POCr and POCl (the normal properties used with sediment quality and vegetation)
                    call SearchProperty(POCr, RefreactaryOrganicC_, .false., STAT = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DecayPeyrard - ModulePorousMediaProperties - ERR011'
                    
                    call SearchProperty(POCl, LabileOrganicC_, .false., STAT = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DecayPeyrard - ModulePorousMediaProperties - ERR012'
                    
                    POCAvailable = .false.
                    
                    allocate(auxPOC(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                    auxPOC(:,:,:) = null_real
                    
                    !Compute sum of POC
                    CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                    !$OMP PARALLEL PRIVATE(I,J,K)
                    
                    !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                    do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
                    do I = Me%WorkSize%ILB, Me%WorkSize%IUB                    
                        if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then 
                            auxPOC(i,j,k) = POCr%Concentration(i,j,k) + POCl%Concentration(i,j,k)
                        endif
                    enddo
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL   
                    
                    !point to the computed poc
                    POC%Concentration => auxPOC
                                                         
                endif
                 
                
                
                !shorten variables
                Porosity => Me%ExtVar%ThetaS
                ThetaF   => Me%ExtVar%ThetaF
                
                CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                !$OMP PARALLEL PRIVATE(I,J,K,AdjustedDensity,NO3Conc,NO3HalfSaturation,  &
                !$OMP& PropFactor,AnaerobioseFactor,POCConc,DOCConc,Decay)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                        
                        !density for equation is in kg/dm3 and then is adjusted to porosity
                        !kg/dm3 = kg/m3 * 1-3 m3/dm3 * -
                        AdjustedDensity = SoilDryDensity%Concentration(i,j,k) * 1E-3 * ((1 - Porosity(i,j,k)) / Porosity(i,j,k))
                        
                        !Nitrate concentration for the equation is in uM
                        !uM = 1e6uM/M * mol/L (M) = 1E6 * g/m3 * 1E-3m3/L / 14g/mol(molecular mass of N (NO3-N))
                        NO3Conc = 1E6 * Property%Concentration(i,j,k) * 1E-3 / 14.
                        NO3HalfSaturation = 1E6 * Property%Evolution%DecayHalfSaturation * 1E-3 / 14.
                        
                        ! adimensional
                        PropFactor = NO3Conc / (NO3HalfSaturation + NO3Conc)
                        
                        !From SedimentQuality model but without lower limit
                        !adimensional. factor goes > 1 so it is limited
                        AnaerobioseFactor = min (0.000304*exp(0.0815*ThetaF(i,j,k)*100), 1.0)
                        
                        !POC concentration for the equation is mg/g
                        !mg/g = mg/kg * 1E-3kg/g
                        POCConc = POC%Concentration(i,j,k) * 1E-3
                        
                        !DOC concentration for the equation is uM
                        !uM = 1E6uM/M * mol/L = 1E-6 * g/m3 * 1E-3 m3/L / 12g/mol(molecular mass of Carbon)
                        DOCConc = 1E6 * DOC%Concentration(i,j,k) * 1E-3 / 12.
                        
                        !the  nitrate decay equation adapted from Peyrard et al. 2011 uses concentration in uM for DOC and NO3
                        !and mg.g-1 for POC
                        !uM.d-1 = (kg/dm3 * 1E3g/kg * d-1 * mg/g * 1E-3g/mg * 1E6uM/M / 12 g/mol) + uM.d-1
                        !uM.d-1 =               uM.d-1                                            + uM.d-1
                        Decay = 0.8 * (AdjustedDensity * POC%Evolution%DecayRate * POCConc * 1E6 / 12.  +       &
                                       DOC%Evolution%DecayRate * DOCConc                                 ) *    &
                                       PropFactor * AnaerobioseFactor
                                                
                        !Convert DecayRate in g/m3.d-1 (convert uM to g/m3 in Nitrate)
                        !g/m3.d-1 = 1E-6M/uM * uM.d-1 * 1E3L/m3 * 14g/mol(molecular mass of N (NO3-N))
                        Property%PropertyDecay(i,j,k) = 1E-6 * Decay * 1E3 * 14.
                        
                        !Zero order decay
                        !P = P0 -kt. k is g/m3.day-1                      
                        !g/m3 = g/m3 - g/m3.d-1  * d
                        !decay can not be higher that OldConc
                        Property%Concentration(I,J,K) = max (Property%Concentration(I,J,K)-(Property%PropertyDecay(i,j,k)*DT),0.0)
                                               
                    endif
                
                enddo
                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL                    
                
                if (.not. POCAvailable) then
                    deallocate (auxPOC)
                endif
                
            endif
        
        !elseif (Property%ID%IDNumber == OtherProperty) then
        
        
        endif              
            
    end subroutine DecayPeyrard

    !-----------------------------------------------------------------------------    
    
    ! This subroutine is responsable for defining       
    ! the next time to actualize the value of each      
    ! property                                          
    subroutine Actualize_Time_Evolution

        !Local--------------------------------------------------------------
        type (T_Property), pointer :: Property
        type (T_Time    )          :: Actual

        !----------------------------------------------------------------------

        Property => Me%FirstProperty  

        Actual = Me%ExtVar%Now

do1 :   do while (associated(Property))

cd1:        if (Property%Evolution%Variable .and. Property%Evolution%DTIntervalAssociated) then
cd2 :       if (Actual.GE.Property%Evolution%NextCompute) then
                    Property%Evolution%LastCompute = Property%Evolution%NextCompute
                    Property%Evolution%NextCompute = Property%Evolution%NextCompute &
                                                   + Property%Evolution%DTInterval
            end if cd2
            end if cd1


            Property => Property%Next
        end do do1   

        nullify(Property)


    end subroutine Actualize_Time_Evolution

    
    !--------------------------------------------------------------------------


    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX
!        real, dimension(:,:), pointer           :: SurfaceDiffusivity

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "OutPut_TimeSeries")


        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data3D = PropertyX%Concentration,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR010'

                if (PropertyX%Evolution%AdvectionDiffusion) then

                    call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                        Data2D = PropertyX%ConcentrationOnInfColumn,   &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR020'
                    
!                    call WriteTimeSerie(Me%ObjTimeSerie,                    &
!                                        Data3D = PropertyX%Diffusivity,     &
!                                        STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_)                              &
!                        stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR030'
!                    
!                    allocate(SurfaceDiffusivity(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
!                    
!                    !Only For Debug
!                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!                        SurfaceDiffusivity(i, j) =  PropertyX%Diffusivity(i, j, Me%WorkSize%KUB +1)
!                    enddo            
!                    enddo         
!
!                    call WriteTimeSerie(Me%ObjTimeSerie,                    &
!                                        Data2D = SurfaceDiffusivity,        &
!                                        STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_)                              &
!                        stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR040'
!                        
!                    deallocate(SurfaceDiffusivity)
                endif
                !output decay rate in case of decay (g/m3.day-1 for dissolved and mg/kg.day-1 in particulated)
                if (PropertyX%Evolution%Decay) then
                    call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                        Data3D = PropertyX%PropertyDecay,   &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                              &
                        stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR045' 
                        
                    if (PropertyX%OutputAverageDecay) then
                        call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                            Data2D = PropertyX%AverageAquiferDecay,   &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                    &
                            stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR046'   

                        call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                            Data2D = PropertyX%AverageVadozeDecay,    &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                    &
                            stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR047'                              
                                             
                    endif  
                    
                    if (PropertyX%OutputIntegratedDecay) then
                        call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                            Data2D = PropertyX%IntegratedDecay,       &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                    &
                            stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR047.5'                     
                    endif
                                         
                endif 

                if (PropertyX%OutputAverageConc) then

                    call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                        Data2D = PropertyX%AverageAquiferConc,    &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                    &
                        stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR048'   

                    call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                        Data2D = PropertyX%AverageVadozeConc,     &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                    &
                        stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR049' 
                    
                endif
                                               
            endif
            PropertyX=>PropertyX%Next
        enddo
        
        if (Me%CalculateECw) then
            call WriteTimeSerie(Me%ObjTimeSerie,     &
                                Data3D = Me%ECw,     &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR050'        
        endif

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "OutPut_TimeSeries")

    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------


    subroutine OutPut_HDF

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        type(T_Time)                                :: Actual, LastTime, EndTime
        integer                                     :: OutPutNumber
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
             
        !Begin----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "OutPut_HDF")

        Actual   = Me%ExtVar%Now
        EndTime  = Me%ExtVar%EndTime
        LastTime = Me%LastOutPutHDF5
         
        OutPutNumber = Me%OutPut%NextOutput

TNum:   if (OutPutNumber <= Me%OutPut%Number)            then 
TOut:       if (Actual .GE. Me%OutPut%OutTime(OutPutNumber)) then 
                
                call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))

First:          if (LastTime.LT.Actual) then 
                    
                    !Writes Time
                    TimePtr => AuxTime
                    call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR010'

                    call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                         Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR020'
           
                    Me%LastOutPutHDF5 = Actual
       
                endif First

                !Limits 
                call HDF5SetLimits   (Me%ObjHDF5,        &
                                      Me%WorkSize%ILB,   &
                                      Me%WorkSize%IUB,   &
                                      Me%WorkSize%JLB,   &
                                      Me%WorkSize%JUB,   &
                                      Me%WorkSize%KLB-1, &
                                      Me%WorkSize%KUB,   &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaProperties - ERR030'


                !Vertical 
                call HDF5WriteData  ( Me%ObjHDF5,  "/Grid/VerticalZ", & 
                                     "Vertical",   "m",               & 
                                      Array3D      = Me%ExtVar%SZZ,   &
                                      OutputNumber = OutPutNumber,    &
                                      STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaProperties - ERR040'
                
                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5,                                &
                                      Me%WorkSize%ILB,                           &
                                      Me%WorkSize%IUB,                           &
                                      Me%WorkSize%JLB,                           &
                                      Me%WorkSize%JUB,                           &
                                      Me%WorkSize%KLB,                           &
                                      Me%WorkSize%KUB,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR050'

                !Writes the Open Points
                call HDF5WriteData   (Me%ObjHDF5, "//Grid/OpenPoints",              &
                                      "OpenPoints", "-",                            &
                                      Array3D = Me%ExtVar%OpenPoints3D,             &
                                      OutputNumber = OutPutNumber,                  &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR060'


                PropertyX => Me%FirstProperty
                do while (associated(PropertyX))

                    if (PropertyX%OutputHDF) then
 
                        call HDF5WriteData   (Me%ObjHDF5,                                    &
                                              "/Results/"//trim(PropertyX%ID%Name),          &
                                              trim(PropertyX%ID%Name),                       &
                                              trim(PropertyX%ID%Units),                      &
                                              Array3D = PropertyX%Concentration,             &
                                              OutputNumber = OutPutNumber,                   &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR070'
                        
                        if (PropertyX%Evolution%Decay) then
                            call HDF5WriteData   (Me%ObjHDF5,                                          &
                                                  "/Results/"//trim(PropertyX%ID%Name)//"_DecayRate",  &
                                                  trim(PropertyX%ID%Name)//"_DecayRate",               &
                                                  trim(PropertyX%ID%Units)//".day-1",                  &
                                                  Array3D = PropertyX%PropertyDecay,                   &
                                                  OutputNumber = OutPutNumber,                         &
                                                  STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR075'
                            
                            if (PropertyX%OutputAverageDecay) then
                                
                                !call ComputeAverageDecay (PropertyX)

                                !Sets limits for next write operations
                                call HDF5SetLimits   (Me%ObjHDF5,                                &
                                                      Me%WorkSize%ILB,                           &
                                                      Me%WorkSize%IUB,                           &
                                                      Me%WorkSize%JLB,                           &
                                                      Me%WorkSize%JUB,                           &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR080'
                                
                                call HDF5WriteData   (Me%ObjHDF5,                                                 &
                                                      "/Results/"//trim(PropertyX%ID%Name)//"_AvrgAquiferDecay",  &
                                                      trim(PropertyX%ID%Name)//"_AvrgAquiferDecay",               &
                                                      trim(PropertyX%ID%Units)//".day-1",                         &
                                                      Array2D = PropertyX%AverageAquiferDecay,                    &
                                                      OutputNumber = OutPutNumber,                                &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR090' 

                                call HDF5WriteData   (Me%ObjHDF5,                                                 &
                                                      "/Results/"//trim(PropertyX%ID%Name)//"_AvrgVadozeDecay",   &
                                                      trim(PropertyX%ID%Name)//"_AvrgVadozeDecay",                &
                                                      trim(PropertyX%ID%Units)//".day-1",                         &
                                                      Array2D = PropertyX%AverageVadozeDecay,                     &
                                                      OutputNumber = OutPutNumber,                                &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR0100'  
                            
                            endif    
                            
                            if (PropertyX%OutputIntegratedDecay) then        
                                !Sets limits for next write operations
                                call HDF5SetLimits   (Me%ObjHDF5,                                &
                                                      Me%WorkSize%ILB,                           &
                                                      Me%WorkSize%IUB,                           &
                                                      Me%WorkSize%JLB,                           &
                                                      Me%WorkSize%JUB,                           &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR0101'
                                
                                call HDF5WriteData   (Me%ObjHDF5,                                                 &
                                                      "/Results/"//trim(PropertyX%ID%Name)//"_IntegratedDecay",   &
                                                      trim(PropertyX%ID%Name)//"_IntegratedDecay",                &
                                                      "kg/ha",                                                    &
                                                      Array2D = PropertyX%IntegratedDecay,                        &
                                                      OutputNumber = OutPutNumber,                                &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR0102'
                            endif                                            
                            
                        endif
                    
                        if (PropertyX%OutputAverageConc) then
                            
                            !call ComputeAverageConc (PropertyX)

                            !Sets limits for next write operations
                            call HDF5SetLimits   (Me%ObjHDF5,                                &
                                                  Me%WorkSize%ILB,                           &
                                                  Me%WorkSize%IUB,                           &
                                                  Me%WorkSize%JLB,                           &
                                                  Me%WorkSize%JUB,                           &
                                                  STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR0110'
                            
                            call HDF5WriteData   (Me%ObjHDF5,                                                 &
                                                  "/Results/"//trim(PropertyX%ID%Name)//"_AvrgAquiferConc",   &
                                                  trim(PropertyX%ID%Name)//"_AvrgAquiferConc",                &
                                                  trim(PropertyX%ID%Units),                                   &
                                                  Array2D = PropertyX%AverageAquiferConc,                     &
                                                  OutputNumber = OutPutNumber,                                &
                                                  STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR0120' 

                            call HDF5WriteData   (Me%ObjHDF5,                                                 &
                                                  "/Results/"//trim(PropertyX%ID%Name)//"_AvrgVadozeConc",    &
                                                  trim(PropertyX%ID%Name)//"_AvrgVadozeConc",                 &
                                                  trim(PropertyX%ID%Units),                                   &
                                                  Array2D = PropertyX%AverageVadozeConc,                      &
                                                  OutputNumber = OutPutNumber,                                &
                                                  STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR0130'  
                        
                        endif                    
                    
                    endif                    
                    
                    PropertyX => PropertyX%Next

                enddo

                Me%OutPut%NextOutput = OutPutNumber + 1

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR0140'
            
            endif  TOut
        endif  TNum

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "OutPut_HDF")


    end subroutine OutPut_HDF

    !----------------------------------------------------------------------------
    !the routimne goes vertically (K) from bottom to top summing prop mass and volume to compute
    !average concentration at end of aquifer ot top of soil. if aquifer or vadoze zone do not exist
    !the value is null_real
    subroutine ComputeAverageConc(PropertyX)
    
        !Arguments---------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        !Local-------------------------------------------------------------------
        type (T_Property), pointer                  :: SoilDryDensity       
        integer, dimension(:,:), pointer            :: GWlayer
        integer                                     :: STAT_CALL, CHUNK, i, j, k
        real                                        :: SoilMass, PropMass, SumPropMass, SumSoilMass
        real                                        :: WaterVolume, SumWaterVolume
        !Begin-------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ComputeAverageConc")
        
        !this is the first non saturated layer
        call GetGWLayer   (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeAverageConc - ModulePorousMediaProperties - ERR001'

!        call GetGeometryKFloor(Me%ObjGeometry,                                          &
!                               Z    = Me%ExtVar%KFloor,                                 &
!                               STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                      &
!            call SetError(FATAL_, INTERNAL_, "ComputeAverageConc - ModulePorousMediaProperties. ERR002")
        
        !if not computed will be stupid number (e.g. no vadoze or no aquifer)
        PropertyX%AverageAquiferConc = null_real
        PropertyX%AverageVadozeConc  = null_real
        
        if (Check_Particulate_Property(PropertyX%ID%IDNumber)) then

            call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ComputeAverageConc - ModulePorousMediaProperties - ERR003'

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(I,J,K,SumPropMass,SumSoilMass,PropMass,SoilMass)
                        
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            SumPropMass = 0.
            SumSoilMass = 0.
            
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                
                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then        
                    
                    ! [kgsoil] = [kgsoil/m3cell] * [m3cell]) 
                    SoilMass  = SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
                                                                    
                    ![mg]     =  [mg/kgsoil] * [kgsoil]
                    PropMass  = PropertyX%Concentration(i,j,k) * SoilMass                   
                    
                    !mg
                    SumPropMass   = SumPropMass + PropMass
                    !kgsoil
                    SumSoilMass   = SumSoilMass + SoilMass 
                    
                    !reached aquifer top (can be top of soil if all soil saturated) - compute aquifer conc
                    if (K == GWlayer(i,j)) then
                        
                        ![mg/kgsoil]                  = [mg] / [kgsoil]
                        PropertyX%AverageAquiferConc(i,j) = SumPropMass / SumSoilMass
                        
                        !reset
                        SumPropMass = 0.
                        SumSoilMass = 0.
                    
                    !reached top of soil - compute vadoze conc. In case of all soil saturated condition is trapped
                    !above and no vadoze concentration is computed here                       
                    elseif (K == Me%WorkSize%KUB) then
                        
                        ![mg/kgsoil]                  = [mg] / [kgsoil]
                        PropertyX%AverageVadozeConc(i,j) = SumPropMass / SumSoilMass                 

                    endif
                        
                endif
                
            enddo
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL 
        
        else

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(I,J,K,SumPropMass,SumWaterVolume,PropMass,WaterVolume)
                        
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            SumPropMass    = 0.
            SumWaterVolume = 0.
            
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                
                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then        
                    
                    ! [m3] = [m3H20/m3cell] * [m3cell]) 
                    WaterVolume  = Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k) 
                                                                    
                    ![g]     =  [g/m3] * [m3]
                    PropMass  = PropertyX%Concentration(i,j,k) * WaterVolume                   
                    
                    !g
                    SumPropMass    = SumPropMass + PropMass
                    !m3
                    SumWaterVolume = SumWaterVolume + WaterVolume 
                    
                    !reached aquifer top (can be top of soil if all soil saturated) - compute aquifer conc
                    if (K == GWlayer(i,j)) then
                        
                        ![mg/kgsoil]                  = [mg] / [m3]
                        PropertyX%AverageAquiferConc(i,j) = SumPropMass / SumWaterVolume
                        
                        !reset
                        SumPropMass    = 0.
                        SumWaterVolume = 0.
                    
                    !reached top of soil - compute vadoze conc. In case of all soil saturated condition is trapped
                    !above and no vadoze concentration is computed here
                    elseif (K == Me%WorkSize%KUB) then
                        
                        ![mg/kgsoil]                  = [mg] / [m3]
                        PropertyX%AverageVadozeConc(i,j) = SumPropMass / SumWaterVolume                 

                    endif
                        
                endif
                
            enddo
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL         
        
        endif

!        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%KFloor,       STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                      &
!            call SetError(FATAL_, INTERNAL_, "ComputeAverageConc - ModulePorousMediaProperties. ERR10")

        call UnGetPorousMedia (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeAverageConc - ModulePorousMediaProperties - ERR020'      

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ComputeAverageConc")
    
    end subroutine ComputeAverageConc
    
    !----------------------------------------------------------------------------

    !the routimne goes vertically (K) from bottom to top averaging decay 
    !at end of aquifer ot top of soil. if aquifer or vadoze zone do not exist
    !the value is null_real
    subroutine ComputeAverageDecay(PropertyX)
    
        !Arguments---------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        !Local-------------------------------------------------------------------
        integer, dimension(:,:), pointer            :: GWlayer
        integer                                     :: STAT_CALL, CHUNK, i, j, k
        real                                        :: SumDecay, SumVolume
        !Begin-------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ComputeAverageDecay")
        
        !this is the first non saturated layer
        call GetGWLayer   (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeAverageDecay - ModulePorousMediaProperties - ERR001'

!        call GetGeometryKFloor(Me%ObjGeometry,                                          &
!                               Z    = Me%ExtVar%KFloor,                                 &
!                               STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                      &
!            call SetError(FATAL_, INTERNAL_, "ComputeAverageDecay - ModulePorousMediaProperties. ERR002")
        
        !if not computed will be stupid number (e.g. no vadoze or no aquifer)
        PropertyX%AverageAquiferDecay = null_real
        PropertyX%AverageVadozeDecay  = null_real
        

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(I,J,K,SumDecay,SumVolume)
                    
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
        SumDecay     = 0.   
        SumVolume    = 0.
        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            
            if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then        
                
                ![g/m3.day-1 or mg/kg.day-1].[m3]
                SumDecay = SumDecay + (PropertyX%PropertyDecay(i,j,k) * Me%ExtVar%CellVolume(i,j,k))
                SumVolume = SumVolume + Me%ExtVar%CellVolume(i,j,k)
                               
                !reached aquifer top (can be top of soil if all soil saturated) - compute aquifer conc
                if (K == GWlayer(i,j)) then
                    
                    ![g/m3.day-1 or mg/kg.day-1]               
                    PropertyX%AverageAquiferDecay(i,j) = SumDecay / SumVolume
                    
                    !reset
                    SumDecay     = 0.
                    SumVolume    = 0.
                
                !reached top of soil - compute vadoze conc. In case of all soil saturated condition is trapped
                !above and no vadoze concentration is computed here                      
                elseif (K == Me%WorkSize%KUB) then
                    
                    ![g/m3.day-1 or mg/kg.day-1]                 
                    PropertyX%AverageVadozeDecay(i,j) = SumDecay / SumVolume                 

                endif
                    
            endif
            
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL         
        

!        call UnGetGeometry( Me%ObjGeometry, Me%ExtVar%KFloor,       STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                      &
!            call SetError(FATAL_, INTERNAL_, "ComputeAverageConc - ModulePorousMediaProperties. ERR10")

        call UnGetPorousMedia (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeAverageDecay - ModulePorousMediaProperties - ERR020'      

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ComputeAverageDecay")
    
    end subroutine ComputeAverageDecay
    
    !----------------------------------------------------------------------------
    
    !the routimne goes vertically (K) from bottom to top integrating decay per cell column for all simulation
    subroutine ComputeIntegratedDecay(PropertyX)
    
        !Arguments---------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        !Local-------------------------------------------------------------------
        type (T_Property), pointer                  :: SoilDryDensity
        integer                                     :: STAT_CALL, CHUNK, i, j, k
        real(8)                                     :: SumDecay, SoilMass, WaterVolume, DT
        !Begin-------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ComputeIntegratedDecay")
                
        DT = Me%ExtVar%DT/86400.

        if (Check_Particulate_Property(PropertyX%ID%IDNumber)) then

            call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ComputeIntegratedDecay - ModulePorousMediaProperties - ERR003'

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(I,J,K,SumDecay,SoilMass)
                        
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            SumDecay = 0.
            
            !integrate in vertical
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                
                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then        
                    
                    ! [kgsoil] = [kgsoil/m3cell] * [m3cell]) 
                    SoilMass  = SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
                                                                    
                    ![mg]     =            [mg/kgsoil.day] * [kgsoil] * [day]
                    SumDecay  = SumDecay + PropertyX%PropertyDecay(i,j,k) * SoilMass * DT    
                   
                endif
                
            enddo
            
            !integrate in time
            ![kg/ha]                       =                                    [mg] * 1E-6 [kg/mg] / ( [m2] * 1E-4 [ha/m2] )
            PropertyX%IntegratedDecay(i,j) = PropertyX%IntegratedDecay(i,j) + (SumDecay * 1E-6 / (Me%ExtVar%Area(i,j) * 1E-4))
            
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL                       

        else
        
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(I,J,K,WaterVolume,SumDecay)
                        
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            SumDecay     = 0.   
            
            !integrate in vertical
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                
                if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then        
 
                    ![m3]        = [m3H20/m3cell] * [m3cell]) 
                    WaterVolume  = Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k) 
                   
                    ![mg]    =            [g/m3.day-1].[m3].[day]
                    SumDecay = SumDecay + (PropertyX%PropertyDecay(i,j,k) * WaterVolume * DT)
                       
                endif
                
            enddo
            
            !integrate in time
            ![kg/ha]                       =                                  [mg] * 1E-6 [kg/mg] / ( [m2] * 1E-4 [ha/m2] )
            PropertyX%IntegratedDecay(i,j) = PropertyX%IntegratedDecay(i,j) + (SumDecay * 1E-6 / (Me%ExtVar%Area(i,j) * 1E-4))
            
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL         
            
        endif
    

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ComputeIntegratedDecay")
    
    end subroutine ComputeIntegratedDecay
    
    !----------------------------------------------------------------------------    
    
    
    subroutine Output_Boxes_Mass

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL, CHUNK
        type (T_Property), pointer                  :: CurrProperty, SoilDryDensity
        real                                        :: ConversionFactor
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "Output_Boxes")

       
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Output_Boxes - ModulePorousMediaProperties - ERR10'        
        
        call GetWaterContent    (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Output_Boxes - ModulePorousMediaProperties - ERR020'

        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ    = Me%ExtVar%CellVolume,                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Output_Boxes - ModulePorousMediaProperties - ERR030'
        
        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR02'

        CurrProperty => Me%FirstProperty
        do while (associated(CurrProperty)) 
            
            if (CurrProperty%BoxTimeSerie) then
               
                Me%CellMass(:,:,:) = 0.
                
                if (Check_Particulate_Property(CurrProperty%ID%IDNumber)) then

                    CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                    !$OMP PARALLEL PRIVATE(I,J,K, ConversionFactor)
                    
                    !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                                
                        if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        
                            ! [Lwater/kgsed] = [m3water/m3cell]*[m3cell]*[1000L/m3] / ([kgsed/m3cell] * [m3cell]) 
                            ConversionFactor  = Me%ExtVar%WaterContent(i,j,k) *  Me%ExtVar%CellVolume(i,j,k) * 1000. / &
                                            (SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k))
                                                                            
                            !g =  mg/kgsoil / L/kgsoil * m3H20/m3cell * m3cell * 1E3L/m3 * 1E-3 g/mg
                            Me%CellMass(i,j,k) = CurrProperty%Concentration(i,j,k) / ConversionFactor *  &
                                                  Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)    

                        endif

                    enddo
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL                    
                
                else

                    CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
                    !$OMP PARALLEL PRIVATE(I,J,K, ConversionFactor)
                    
                    !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                                
                        if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        
                            !g =  mg/L  * m3H20/m3cell * m3cell * 1E3L/m3 * 1E-3 mg/g
                            Me%CellMass(i,j,k) = CurrProperty%Concentration(i,j,k)  *  &
                                                  Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)    

                        endif

                    enddo
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL                  
                
                endif
                
                call BoxDif(Me%ObjBoxDif, Me%CellMass,                      &
                            "soil_"//trim(adjustl(CurrProperty%ID%Name)),    &
                            Me%ExtVar%WaterPoints3D,                        &
                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  &
                    stop 'Output_Boxes - ModulePorousMediaProperties - ERR01'        
                
                Me%CellMass(:,:,:) = null_real                                    

            endif
            
            CurrProperty => CurrProperty%Next
        
        end do 
        
        nullify (CurrProperty)
        nullify (SoilDryDensity)
        
        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Output_Boxes - ModulePorousMediaProperties - ERR040'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Output_Boxes - ModulePorousMediaProperties - ERR050'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%CellVolume,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'Output_Boxes - ModulePorousMediaProperties - ERR060'

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "Output_Boxes")    
    
    end subroutine Output_Boxes_Mass

    !----------------------------------------------------------------------------
    
    subroutine Output_Boxes_Fluxes (Property)
    
        !Arguments---------------------------------------------------------------
        type (T_Property)                         :: Property
        !Local-------------------------------------------------------------------
        integer                                   :: CHUNK, STAT_CALL, i, j, k
        !Begin-------------------------------------------------------------------

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                                
        !$OMP PARALLEL PRIVATE(I,J,K)
do2 :   do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do4 :   do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%Fluxes%MassFluxesX (I,J,K) = Me%Fluxes%AdvFluxX(I,J,K) + Me%Fluxes%DifFluxX (I,J,K)
            Me%Fluxes%MassFluxesY (I,J,K) = Me%Fluxes%AdvFluxY(I,J,K) + Me%Fluxes%DifFluxY (I,J,K)
        end do do4
        end do do3
        !$OMP END DO NOWAIT
        end do do2
        !$OMP END PARALLEL

        if (Me%WorkSize%KUB > Me%WorkSize%KLB) then
            !$OMP PARALLEL PRIVATE(I,J,K)
do5 :       do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do6 :       do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do7 :       do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%Fluxes%MassFluxesZ (I,J,K) = Me%Fluxes%AdvFluxZ(I,J,K) + Me%Fluxes%DifFluxZ (I,J,K)
            end do do7
            end do do6
            !$OMP END DO NOWAIT
            end do do5
            !$OMP END PARALLEL
        endif

        !Integration of fluxes
        call BoxDif(Me%ObjBoxDif,                        &
                    Me%Fluxes%MassFluxesX,               &
                    Me%Fluxes%MassFluxesY,               &
                    Me%Fluxes%MassFluxesZ,               &
                    "soil_"//trim(Property%ID%Name),     &
                    Me%ExtVar%WaterPoints3D,             &
                    STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                     &
        stop 'Output_Boxes_Fluxes - ModulePorousMediaProperties - ERR300'
                
    end subroutine Output_Boxes_Fluxes
    
    !--------------------------------------------------------------------------

    subroutine WriteFinalFile

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        integer                                     :: OutPutNumber
        integer                                     :: HDF5_CREATE
        character(LEN = PathLength)                 :: FileName
        integer                                     :: ObjHDF5
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        type (T_Time)                               :: Actual            
        real                                        :: Total_Mass_Created
        character (Len = StringLength)              :: str_mass_created, string_to_be_written          
        !Begin----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "WriteFinalFile")

        !Gets a pointer to Topography
        call GetGridData        (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR00'

        call GetGeometryDistances (Me%ObjGeometry,                                      &
                                  SZZ         = Me%ExtVar%SZZ,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR01'

        !OpenPoints3D
        call GetOpenPoints3D    (Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR02'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Checks if it's at the end of the run 
        !or !if it's supposed to overwrite the final HDF file
        if ((Me%ExtVar%Now == Me%ExtVar%EndTime) .or. Me%Output%RestartOverwrite) then

            filename = trim(Me%Files%FinalFile)

        else

            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%ExtVar%Now))//".fin")

        endif


        ObjHDF5 = 0
        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5,                                                     &
                            trim(filename),                                              &
                            HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalFile - ModulePorousMediaProperties - ERR10'


        Actual   = Me%ExtVar%Now
         
        call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
        !Writes Time
        TimePtr => AuxTime
        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR11'

        call HDF5WriteData  (ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                             Array1D = TimePtr, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR12'
                    
        !Limits 
        call HDF5SetLimits   (ObjHDF5,        &
                              Me%WorkSize%ILB,   &
                              Me%WorkSize%IUB,   &
                              Me%WorkSize%JLB,   &
                              Me%WorkSize%JUB,   &
                              Me%WorkSize%KLB-1, &
                              Me%WorkSize%KUB,   &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR030'

        !Vertical 
        call HDF5WriteData  ( ObjHDF5,  "/Grid/VerticalZ", & 
                             "Vertical",   "m",               & 
                              Array3D      = Me%ExtVar%SZZ,   &
                              OutputNumber = OutPutNumber,    &
                              STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR040'
        
        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5,                                &
                              Me%WorkSize%ILB,                           &
                              Me%WorkSize%IUB,                           &
                              Me%WorkSize%JLB,                           &
                              Me%WorkSize%JUB,                           &
                              Me%WorkSize%KLB,                           &
                              Me%WorkSize%KUB,                           &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaproperties - ERR050'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR25'
        
        !Writes the Grid
        call HDF5WriteData      (ObjHDF5, "/Grid", "Topography", "m",                    &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR035'


        !Writes the Open Points
        call HDF5WriteData   (ObjHDF5, "//Grid/OpenPoints",              &
                              "OpenPoints", "-",                            &
                              Array3D = Me%ExtVar%OpenPoints3D,             &
                              OutputNumber = OutPutNumber,                  &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaproperties - ERR060'


        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))

            call HDF5SetLimits   (ObjHDF5,                                &
                                  Me%WorkSize%ILB,                           &
                                  Me%WorkSize%IUB,                           &
                                  Me%WorkSize%JLB,                           &
                                  Me%WorkSize%JUB,                           &
                                  Me%WorkSize%KLB,                           &
                                  Me%WorkSize%KUB,                           &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaproperties - ERR065'

            call HDF5WriteData   (ObjHDF5,                                    &
                                  "/Results/"//trim(PropertyX%ID%Name),          &
                                  trim(PropertyX%ID%Name),                       &
                                  trim(PropertyX%ID%Units),                      &
                                  Array3D = PropertyX%Concentration,             &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaproperties - ERR070'


            if (PropertyX%Evolution%MinConcentration .and. Me%ExtVar%Now == Me%ExtVar%EndTime) then

                call HDF5WriteData   (ObjHDF5,                                        &
                                      "/Results/"//trim(PropertyX%ID%Name)//" Mass Created",& 
                                      "Property Mass Created",                        &
                                      "g",                                            &
                                      Array3D = PropertyX%Mass_Created,               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR10.5'

                !g/1000 = kg, to avoid big numbers
                Total_Mass_Created = SUM(PropertyX%Mass_Created)/1000

                write(str_mass_created, '(f20.8)') Total_Mass_Created
      
                string_to_be_written = 'Due to MinConcentration PMP Total mass (kg) created on property ' //&
                                        trim(adjustl(adjustr(PropertyX%ID%Name)))//' = ' //&
                                        trim(adjustl(adjustr(str_mass_created))) 
            
                !Writes total mass created to "Error_and_Messages.log" file
                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

            endif

            PropertyX => PropertyX%Next

        enddo

        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaproperties - ERR080'
            

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR085'
        
        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%SZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR100'

        call UnGetGridData              (Me%ObjGridData, Me%ExtVar%Topography,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModulePorousMediaProperties - ERR120'


        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "WriteFinalFile")


    end subroutine WriteFinalFile

    !----------------------------------------------------------------------------

    real function Tortuosity(WC, Porosity)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: WC
        real, intent(IN)                            :: Porosity
        !Begin----------------------------------------------------------------

        !tortuosity = (WC**(10/3))/(porosity)
        Tortuosity = (WC**(7./3.))/(Porosity**2)
        

    end function Tortuosity   


    !----------------------------------------------------------------------------

    subroutine ProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: nProperties
        integer                                             :: n                                    

        nProperties = Me%PropertiesNumber 
    
        n=1
        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))

            call WriteProfile(Me%ObjProfile,                                        &
                              Data3D = PropertyX%Concentration,                              &
                              SZZ    = Me%ExtVar%SZZ,                               &
                              STAT   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR01'

            PropertyX=>PropertyX%Next

        enddo
        
    end subroutine ProfileOutput

    !------------------------------------------------------------------------------ 

    subroutine CalculateTotalStoredMass

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL
        type (T_Property), pointer                  :: CurrProperty, SoilDryDensity
        real                                        :: ConversionFactor
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "CalculateTotalStoredMass")

        
        CurrProperty => Me%FirstProperty
        
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR10'        
        
        call GetWaterContent    (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR020'

        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ    = Me%ExtVar%CellVolume,                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR030'

        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR02'

        
        do while (associated(CurrProperty)) 
            
            if (Check_Particulate_Property(CurrProperty%ID%IDNumber)) then
               
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            
                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                    
                        ! [Lwater/kgsed] = [m3water/m3cell]*[m3cell]*[1000L/m3] / ([kgsed/m3cell] * [m3cell]) 
                        ConversionFactor  = Me%ExtVar%WaterContent(i,j,k) *  Me%ExtVar%CellVolume(i,j,k) * 1000. / &
                                        (SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k))
                                                                
                        !kg = kg + m3H20/m3cell * m3cell * mg/kgsed / L/kgsed * 1E3 L/m3 * 1E-3 g/mg * 1E-3 kg/g
                        CurrProperty%MB%TotalStoredMass = CurrProperty%MB%TotalStoredMass                                  &
                                                          + ( Me%ExtVar%WaterContent(i,j,k)* Me%ExtVar%CellVolume(i,j,k)   &
                                                           * CurrProperty%Concentration(i,j,k) / ConversionFactor * 1E-3    )
                            
                    endif

                enddo
                enddo
                enddo               
               
            else
            
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            
                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                        !kg = kg + (m3H20/m3cell * m3cell * g/m3 * 1E-3 kg/g
                        CurrProperty%MB%TotalStoredMass = CurrProperty%MB%TotalStoredMass                                  &
                                                          + ( Me%ExtVar%WaterContent(i,j,k)* Me%ExtVar%CellVolume(i,j,k)   &
                                                           * CurrProperty%Concentration(i,j,k) * 1E-3    )
                            
                    endif

                enddo
                enddo
                enddo
            
            endif
                        
            CurrProperty => CurrProperty%Next
        end do 
        
        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR040'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR050'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%CellVolume,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR060'

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "CalculateTotalStoredMass")


    end subroutine CalculateTotalStoredMass


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillPorousMediaProperties(ObjPorousMediaPropertiesID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjPorousMediaPropertiesID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers,STAT_CALL  
        type(T_property), pointer           :: PropertyX
        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mPorousMediaProperties_,  Me%InstanceID)
            
            call WriteFinalFile

            if (Me%Output%IntegratedDecay_ON) then
            
                PropertyX => Me%FirstProperty

                do while(associated(PropertyX))
                    

                    if (PropertyX%OutputIntegratedDecay) then
                    
                        call WriteGridData  (PropertyX%IntegratedDecayFile,        &
                             COMENT1          = "IntegratedDecayFile",             &
                             COMENT2          = "IntegratedDecayFile",             &
                             HorizontalGridID = Me%ObjHorizontalGrid,              &
                             FillValue        = -99.0,                             &
                             OverWrite        = .true.,                            &
                             GridData2D_Real  = PropertyX%IntegratedDecay,         &
                             STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR001'
                    
                    endif
                    
                    PropertyX => PropertyX%Next
                    
                end do                       
            

            endif

            PropertyX => Me%FirstProperty
            
            do while (associated(PropertyX)) 
                if(PropertyX%ID%SolutionFromFile)then

                    call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR00'
                end if
                
                PropertyX => PropertyX%Next
            end do 
            
            if (associated(Me%Disper_Longi%Field))then
                deallocate(Me%Disper_Longi%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR01'
                nullify   (Me%Disper_Longi%Field)
            end if

            if (associated(Me%Disper_Trans%Field))then
                deallocate(Me%Disper_Trans%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR02'
                nullify   (Me%Disper_Trans%Field)
            end if
            
            if (nUsers == 0) then

                !Kills the TimeSerie
                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR05'
                endif

               
                if (Me%OutPut%HDF_ON) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillVegetation - ModulePorousMediaProperties  - ERR08'
                endif
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR07'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjGridData)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR07.5'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR08'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR11'
                
                nUsers = DeassociateInstance (mPOROUSMEDIA_,  Me%ObjPorousMedia)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR12'

                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR13'

                nUsers = DeassociateInstance (mMAP_,  Me%ObjMap)
                if (nUsers == 0) stop 'KillPorousMedia - ModulePorousMediaProperties - ERR14'

                if (Me%Coupled%ChainReactions) then
                    call KillChainReactions (Me%ObjChainReactions, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillPorousMedia - ModulePorousMediaProperties - ERR140'
                endif


                if (Me%Output%RateFluxes .or. Me%Output%Boxes_ON) then
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR09'
                endif
                
                call DeallocateVariables

#ifdef _ENABLE_CUDA                
                !Kills ModuleCuda. Do this after DeallocateVariables, since ModuleCuda is needed in DeallocateVariables
                call KillCuda (Me%ObjCuda, STAT = STAT_CALL)
                ! No need to give error yet, Module still has users
#endif _ENABLE_CUDA

                !Deallocates Instance
                call DeallocateInstance ()

                ObjPorousMediaPropertiesID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillPorousMediaProperties
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_PorousMediaProperties), pointer          :: AuxObjPorousMediaProperties
        type (T_PorousMediaProperties), pointer          :: PreviousObjPorousMediaProp

        !Updates pointers
        if (Me%InstanceID == FirstObjPorousMediaProperties%InstanceID) then
            FirstObjPorousMediaProperties => FirstObjPorousMediaProperties%Next
        else
            PreviousObjPorousMediaProp => FirstObjPorousMediaProperties
            AuxObjPorousMediaProperties      => FirstObjPorousMediaProperties%Next
            do while (AuxObjPorousMediaProperties%InstanceID /= Me%InstanceID)
                PreviousObjPorousMediaProp => AuxObjPorousMediaProperties
                AuxObjPorousMediaProperties      => AuxObjPorousMediaProperties%Next
            enddo

            !Now update linked list
            PreviousObjPorousMediaProp%Next => AuxObjPorousMediaProperties%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance
    
    !--------------------------------------------------------------------------

    subroutine DeallocateVariables

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        !griflet
        integer                 :: p
        type(T_VECGW), pointer  :: VECGW

        !Water Content---------------------------------------------------------
        deallocate (Me%ExtVar%WindVelocity3D   )         
        
        if (Me%Coupled%AdvectionDiffusion) then
            deallocate (Me%WaterVolume)
            deallocate (Me%FluxWCorr)
#ifdef _USE_PAGELOCKED
            ! FreePageLocked will also nullify the pointers and arrays
            call FreePageLocked(Me%ObjCuda, Me%TICOEF3Ptr, Me%TICOEF3)
#else
            deallocate(Me%TICOEF3)
#endif _USE_PAGELOCKED            
            
            deallocate (Me%ThetaAtFaces%ThetaW)
            deallocate (Me%ThetaAtFaces%ThetaU)
            deallocate (Me%ThetaAtFaces%ThetaV)

            deallocate(Me%COEF3_VertAdv%C_Flux)
            deallocate(Me%COEF3_VertAdv%D_Flux)
            deallocate(Me%COEF3_VertAdv%E_Flux)
            deallocate(Me%COEF3_VertAdv%F_Flux)

            if (.not. Me%Vertical1D) then
                deallocate(Me%COEF3_HorAdvXX%C_Flux)
                deallocate(Me%COEF3_HorAdvXX%D_Flux)
                deallocate(Me%COEF3_HorAdvXX%E_Flux)
                deallocate(Me%COEF3_HorAdvXX%F_Flux)
                
                deallocate(Me%COEF3_HorAdvYY%C_Flux)
                deallocate(Me%COEF3_HorAdvYY%D_Flux)
                deallocate(Me%COEF3_HorAdvYY%E_Flux)
                deallocate(Me%COEF3_HorAdvYY%F_Flux)
                
            endif
            
            if (.not. Me%AdvDiff_Explicit) then
#ifdef _USE_PAGELOCKED
                ! FreePageLocked will also nullify the pointers and arrays
                call FreePageLocked(Me%ObjCuda, Me%COEF3%DPtr, Me%COEF3%D)
                call FreePageLocked(Me%ObjCuda, Me%COEF3%EPtr, Me%COEF3%E)
                call FreePageLocked(Me%ObjCuda, Me%COEF3%FPtr, Me%COEF3%F)
#else
                deallocate(Me%COEF3%D)
                deallocate(Me%COEF3%E)
                deallocate(Me%COEF3%F)
#endif _USE_PAGELOCKED

!                deallocate(Me%VECG)
!                deallocate(Me%VECW)
              
            endif
            
            !griflet
            do p = 1, Me%MaxThreads                
                VECGW => Me%THOMAS%VEC(p)
                deallocate(VECGW%G)
                deallocate(VECGW%W)
            enddo 
            deallocate(Me%THOMAS%VEC)
            deallocate(Me%THOMAS%COEF3)
            deallocate(Me%THOMAS)

        endif
        
        deallocate (Me%CellWaterVolume)
        
#ifdef _PHREEQC_   
        if (associated(Me%PhreeqC%CellSoilMass)) then   
            deallocate (Me%PhreeqC%CellSoilMass)    
        endif
#endif      
    
        if (associated(Me%PropertiesList)) then
            deallocate (Me%PropertiesList)
        endif
    
    end subroutine DeallocateVariables 

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjPorousMediaProperties_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaProperties_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjPorousMediaProperties_ID > 0) then
            call LocateObjPorousMediaProperties (ObjPorousMediaProperties_ID)
            ready_ = VerifyReadLock (mPorousMediaProperties_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjPorousMediaProperties (ObjPorousMediaPropertiesID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaPropertiesID

        !Local-----------------------------------------------------------------

        Me => FirstObjPorousMediaProperties
        do while (associated (Me))
            if (Me%InstanceID == ObjPorousMediaPropertiesID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModulePorousMediaProperties - LocateObjPorousMediaProperties - ERR01'

    end subroutine LocateObjPorousMediaProperties

    !--------------------------------------------------------------------------

    subroutine SearchProperty(PropertyX, PropertyXIDNumber, PrintWarning, STAT)


        !Arguments-------------------------------------------------------------
        type(T_Property), optional, pointer         :: PropertyX
        integer         , optional, intent (IN)     :: PropertyXIDNumber
        logical,          optional, intent (IN)     :: PrintWarning
        integer         , optional, intent (OUT)    :: STAT

        !Local-----------------------------------------------------------------

        integer                                     :: STAT_ 
        
        !----------------------------------------------------------------------

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
                if (PrintWarning) write (*,*)'Property Not Found in Module PorousMediaProperties ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SearchProperty

    !--------------------------------------------------------------------------


    subroutine ReadLockExternalVar                

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, GWFlowLink
        !Begin-----------------------------------------------------------------

        call GetOldWaterContent (Me%ObjPorousMedia, Me%ExtVar%WaterContentOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR010'

        call GetWaterContent    (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR020'

        call GetFluxU           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR030'

        call GetFluxV           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR040'
        
        call GetUnsatV          (Me%ObjPorousMedia, Me%ExtVar%UnsatV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR061'

        call GetUnsatU          (Me%ObjPorousMedia, Me%ExtVar%UnsatU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR062'   
        
        call GetFluxW           (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR050'

!        call GetUnsatW          (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'
        
        call GetUnsatWFinal     (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'
        
        call GetThetaS          (Me%ObjPorousMedia, Me%ExtVar%ThetaS, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR65'
        
        call GetThetaF          (Me%ObjPorousMedia, Me%ExtVar%ThetaF, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR66'
        
        call GetPotentialInfiltration (Me%ObjPorousMedia, Me%ExtVar%InfiltrationColumn, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleModulePorousMediaProperties. ERR10.'


        call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                 GridCellArea = Me%ExtVar%Area,                            & 
                                 STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR070'

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR080'

        !LandPoints3D
        call GetLandPoints3D    (Me%ObjMap, Me%ExtVar%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSoilProperties - ERR11'

        !OpenPoints3D
        call GetOpenPoints3D    (Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR085'

        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ    = Me%ExtVar%CellVolume,                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR090'
                       
        call GetGeometryDistances (Me%ObjGeometry,                                      &
                                  SZZ         = Me%ExtVar%SZZ,                          &
                                  DWZ         = Me%ExtVar%DWZ,                          &
                                  DZZ         = Me%ExtVar%DZZ,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR100'

        !AreaU, AreaV
        call GetGeometryAreas(Me%ObjGeometry,                     &
                              AreaU = Me%ExtVar%AreaU,            &
                              AreaV = Me%ExtVar%AreaV,            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR105'

       
        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                    &
                                  DUX         = Me%ExtVar%DUX,                          &
                                  DZY         = Me%ExtVar%DZY,                          &
                                  DVY         = Me%ExtVar%DVY,                          &
                                  DZX         = Me%ExtVar%DZX,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR101'

        call GetGridData  (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR110'


        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesW3D = Me%ExtVar%ComputeFacesW3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR120'


        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesU3D = Me%ExtVar%ComputeFacesU3D,             &
                               ComputeFacesV3D = Me%ExtVar%ComputeFacesV3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR130'
        
        if (Me%ExtVar%CoupledDN) then
            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR0140'

            call GetGWFlowOption (Me%ObjPorousMedia, GWFlowLink, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop ' ModifyAdvectionDiffusion - ModulePorousMediaProperties - ERR01'
            
            if (GWFlowLink /= Layer_) then
                call GetGWFlowToChannels   (Me%ObjPorousMedia, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR150'
            else
                call GetGWFlowToChannelsByLayer   (Me%ObjPorousMedia, Me%ExtVar%FlowToChannelsLayer, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR151'
            endif
            
            call GetGWToChannelsLayers   (Me%ObjPorousMedia, Me%ExtVar%GWFlowBottomLayer, Me%ExtVar%GWFlowTopLayer, &
                                          STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR152'

            call GetGWLayerOld   (Me%ObjPorousMedia, Me%ExtVar%GWlayerOld, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR160'
        endif
        
        if (Me%ExtVar%BoundaryImposed) then
            if (Me%ExtVar%BoundaryImposedWalls) then
                call GetBoundaryFluxWalls (Me%ObjPorousMedia, Me%ExtVar%BoundaryFluxWalls, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR170'
            endif
            if (Me%ExtVar%BoundaryImposedBottom) then
                call GetBoundaryFluxBottom (Me%ObjPorousMedia, Me%ExtVar%BoundaryFluxBottom, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR180'
            endif            
        endif
        
    end subroutine ReadLockExternalVar

    !-----------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, GWFlowLink
        !Begin-----------------------------------------------------------------

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContentOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR010'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR020'


        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR030'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR040'
        
        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%UnsatV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR061'
        
        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%UnsatU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR064'
        
        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR050'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR060'

        call UnGetPorousMedia           (Me%ObjPorousMedia,Me%ExtVar%ThetaS, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR061'

        call UnGetPorousMedia           (Me%ObjPorousMedia,Me%ExtVar%ThetaF, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR062'        
       
        call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR070'
        

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%InfiltrationColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR063'


        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR080'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%LandPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR081'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR085'
        
        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%CellVolume,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR090'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%SZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR100'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%DWZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR110'
        
        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%DZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR111'

        !AreaU
        call UnGetGeometry(Me%ObjGeometry,                                              &
                           Me%ExtVar%AreaU, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR111.1'

        !AreaV
        call UnGetGeometry(Me%ObjGeometry,                                              &
                           Me%ExtVar%AreaV, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR111.2'


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR113'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR114'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR115'


        call UnGetGridData              (Me%ObjGridData, Me%ExtVar%Topography,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR120'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesW3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR130'
        
        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesU3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR140'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesV3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR150'

        if (Me%ExtVar%CoupledDN) then
            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0160'               

            call GetGWFlowOption (Me%ObjPorousMedia, GWFlowLink, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop ' ModifyAdvectionDiffusion - ModulePorousMediaProperties - ERR01'
            
            if (GWFlowLink /= Layer_) then
                call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0170'               
            else
                call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%FlowToChannelsLayer, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0171'               
            endif
            
            call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%GWFlowBottomLayer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0172'               

            call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%GWFlowTopLayer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0173'               

            call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%GWlayerOld, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0180'      
        
        endif
        
        if (Me%ExtVar%BoundaryImposed) then
            if (Me%ExtVar%BoundaryImposedWalls) then        
                call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%BoundaryFluxWalls, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0190'      
            endif
            if (Me%ExtVar%BoundaryImposedBottom) then        
                call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%BoundaryFluxBottom, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0200'      
            endif            
        endif
        
    endsubroutine ReadUnlockExternalVar

    !--------------------------------------------------------------------------
    
end module ModulePorousMediaProperties

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 








