!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land
! MODULE        : PorousMediaProperties
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Fev 2010
! REVISION      : 
! DESCRIPTION   : Module to serve as PorousMediaProperties
!
!------------------------------------------------------------------------------
!
!
!Units in porous media properties
!   Transported properties (soluble)  : g/m3 (or mg/l)                     (needs to convert concentrations to SedimentQuality and 
!                                                                           PREEQC at entrance and exit)
!   Adsorbed properties (non soluble) : mg/kgsoil                          (needs to convert concentrations to PREEQC )
!   Solid Phases (non soluble)        : mg/kgsoil                          (used in PhreeqC to make equilibrium with soil solution)
!   Gas Phases (non soluble)          : mol (mass); atm (partial pressure) (used in PhreeqC to make equilibrium with soil solution)
!   Soil Dry Density                  : kg/m3                              
!
! ADVDIFF_EXPLICIT              : 0/1                [1]        !REMARK: Horizontal diffusion is always explicit
!                                                               !(1 - adv and diff are explicit in all directions; 0 - adv and diff 
!                                                               !are implicit in vertical, horizontal adv may be explicit or impl) 
!   ADVDIFF_ADVECTION_H_IMP_EXP : 0/1                [1]        !(read if ADVDIFF_EXPLICIT : 0; 0 - horiz adv implicit; 
!                                                                1 - horiz adv explicit)
!
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
                                         KillGridData
    use ModuleTimeSerie,          only : StartTimeSerie, WriteTimeSerie, KillTimeSerie         
    use ModuleHorizontalGrid,     only : GetHorizontalGrid, GetGridCellArea,               &
                                         WriteHorizontalGrid, UnGetHorizontalGrid
    use ModuleBasinGeometry,      only : GetBasinPoints, GetRiverPoints,  UnGetBasin 
                                       
    use ModuleFillMatrix,         only : ConstructFillMatrix, GetDefaultValue,             &
                                         KillFillMatrix, ModifyFillMatrix
    use ModuleGeometry
    use ModuleMap
    use ModulePorousMedia,        only : GetOldWaterContent, GetWaterContent, GetFluxU,    &
                                         GetFluxV, GetFluxW, GetUnsatW, GetUnsatV,         &
                                         GetUnsatU, UnGetPorousMedia,                      &
                                         GetThetaS, GetGWFlowToChannels,                   &
                                         GetGWLayerOld, GetPotentialInfiltration,          &
                                         GetGWFlowToChannelsByLayer, GetGWToChannelsLayers,  &
                                         GetGWFlowOption, GetTranspiration
                                         
    use ModuleChainReactions,     only : StartChainReactions, SetCRPropertyConcentration,  &
                                         GetCRPropertiesList, UnGetChainReactions,         &
                                         ModifyChainReactions, KillChainReactions
                                         
#ifdef _PHREEQC_                                          
    use ModuleInterface,          only : ConstructInterface, Modify_Interface, GetPhreeqCID 
#else
    use ModuleInterface,          only : ConstructInterface, Modify_Interface     
#endif
                                         
    use ModuleAdvectionDiffusion, only : StartAdvectionDiffusion, AdvectionDiffusion,      &
                                         GetAdvFlux, GetDifFlux, GetBoundaryConditionList, &
                                         UngetAdvectionDiffusion, KillAdvectionDiffusion

#ifdef _PHREEQC_     
    use ModulePhreeqC
#endif

   implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !PorousMediaProperties
    public  :: ConstructPorousMediaProperties
    private ::      AllocateInstance
    private ::      ReadFileNames
    private ::      ReadGlobalOptions
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
    public  :: GetPMPCoupled
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
    private ::      SoilQualityProcesses
#ifdef _PHREEQC_    
    private ::      SoilChemistryProcesses
#endif        
    private ::      ChainReactionsProcesses
    private ::      Partition_Processes
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
 
    !Types---------------------------------------------------------------------
    
    private :: T_PorousMediaProperties

    type T_RelatedID
        integer                       :: IDNumber = -1
        character(LEN = StringLength) :: name     = ''   
    end type T_RelatedID

    type T_ID
        integer                       :: IDNumber
        character(LEN = StringLength) :: name
        character(LEN = StringLength) :: description
        character(LEN = StringLength) :: units
    end type T_ID

    type T_Property_3D
        type(T_PropertyID)               :: ID
        real, pointer, dimension (:,:,:) :: Field
        real                             :: Scalar
    end type T_Property_3D


    type T_ExtVar
        !Map
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D
        integer, pointer, dimension(:,:,:)          :: OpenPoints3D
        integer, pointer, dimension(:,:,:)          :: LandPoints3D
        integer, dimension(:,:), pointer            :: BasinPoints
        integer, dimension(:,:), pointer            :: RiverPoints
        integer, pointer, dimension(:,:,:)          :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:)          :: ComputeFacesV3D
        integer, pointer, dimension(:,:,:)          :: ComputeFacesW3D         !from basin
        real                                        :: PorousMediapropDT
        type(T_Time)                                :: Now
        type(T_Time)                                :: BeginTime
        type(T_Time)                                :: EndTime
   
        ! from porousMedia
        real,    dimension(:,:,:), pointer          :: UnSatW
        real,    dimension(:,:,:), pointer          :: UnSatV
        real,    dimension(:,:,:), pointer          :: UnSatU
        real,    pointer, dimension(:,:,:)          :: WaterContent
        real,    pointer, dimension(:,:,:)          :: WaterContentOld
        real(8), pointer, dimension(:,:)            :: WaterColumn
        real(8), pointer, dimension(:,:)            :: InfiltrationColumn
        real(8), pointer, dimension(:,:,:)          :: CellVolume        
        real(8), dimension(:,:,:), pointer          :: FluxU
        real(8), dimension(:,:,:), pointer          :: FluxV
        real(8), dimension(:,:,:), pointer          :: FluxW
        real,    pointer, dimension(:,:,:)          :: ThetaS        
        real   , pointer, dimension(:,:  )          :: Area
        real   , pointer, dimension(:,:,:)          :: AreaU
        real   , pointer, dimension(:,:,:)          :: AreaV
        real                                        :: DT
        real   , pointer, dimension(:,:,:)          :: DWZ
        real   , pointer, dimension(:,:,:)          :: DZZ
        real   , pointer, dimension(:,:  )          :: DVY
        real   , pointer, dimension(:,:  )          :: DZX
        real   , pointer, dimension(:,:  )          :: DUX
        real   , pointer, dimension(:,:  )          :: DZY
        real   , pointer, dimension(:,:  )          :: Topography  
        real ,   pointer, dimension(:,:,:)          :: SZZ
        real,    pointer, dimension(:,:)            :: FlowToChannels
        real,    pointer, dimension(:,:,:)          :: FlowToChannelsLayer
        integer, pointer, dimension(:,:)            :: GWFlowBottomLayer
        integer, pointer, dimension(:,:)            :: GWFlowTopLayer
        integer, pointer, dimension(:,:)            :: GWLayerOld
        
        !from vegetation
        logical                                     :: ComputeVegInterfaceFluxes 
        logical, dimension(:,:  ), pointer          :: SoilFluxesActive
        real,    dimension(:,:  ), pointer          :: GrazingBiomass
        real,    dimension(:,:  ), pointer          :: GrazingNitrogen
        real,    dimension(:,:  ), pointer          :: GrazingPhosphorus
        real,    dimension(:,:  ), pointer          :: ManagementAerialBiomass
        real,    dimension(:,:  ), pointer          :: ManagementNitrogen
        real,    dimension(:,:  ), pointer          :: ManagementPhosphorus
        real,    dimension(:,:  ), pointer          :: ManagementRootBiomass
        real,    dimension(:,:  ), pointer          :: DormancyBiomass
        real,    dimension(:,:  ), pointer          :: DormancyNitrogen
        real,    dimension(:,:  ), pointer          :: DormancyPhosphorus
        real,    dimension(:,:  ), pointer          :: FertilNitrateSurface
        real,    dimension(:,:  ), pointer          :: FertilNitrateSubSurface
        real,    dimension(:,:  ), pointer          :: FertilAmmoniaSurface
        real,    dimension(:,:  ), pointer          :: FertilAmmoniaSubSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicNSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicNSubSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicPSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicPSubSurface
        real,    dimension(:,:  ), pointer          :: FertilMineralPSurface
        real,    dimension(:,:  ), pointer          :: FertilMineralPSubSurface
        real,    dimension(:,:  ), pointer          :: NitrogenFraction
        real,    dimension(:,:  ), pointer          :: PhosphorusFraction
        real,    dimension(:,:,:), pointer          :: NitrogenUptake
        real,    dimension(:,:,:), pointer          :: PhosphorusUptake
        real,    dimension(:,:  ), pointer          :: RootDepth
        integer, dimension(:,:  ), pointer          :: TranspirationBottomLayer
        logical                                     :: Grazing
        logical                                     :: Management
        logical                                     :: Dormancy
        logical                                     :: Fertilization
        logical                                     :: ModelWater
        logical                                     :: ModelNitrogen
        logical                                     :: ModelPhosphorus
        logical                                     :: GrowthModel
        logical                                     :: CoupledVegetation
        real                                        :: VegetationDT
        
        logical                                     :: CoupledDN  = .false.
        !from basin
        real,    dimension(:,:  ), pointer          :: WindVelocity2D  !m/s
        real,    dimension(:,:,:), pointer          :: WindVelocity3D  !km/day
     end type T_ExtVar

    type T_OutPut
        type (T_Time), pointer, dimension(:)    :: OutTime
        integer                                 :: NextOutPut
        integer                                 :: Number
        logical                                 :: Yes = .false.
        logical                                 :: TimeSerie_ON
        logical                                 :: HDF_ON
        logical                                 :: Profile_ON
    end type T_OutPut

    type T_AdvectionDiffusion   
        !--For AdvectionDiffusion module use
        integer                                :: BoundaryCondition
        real                                   :: SchmidtNumberH
        real                                   :: SchmidtCoefV
        real                                   :: SchmidtBackgroundV
        real                                   :: DiffusionH_imp_exp
        real                                   :: ImplicitH_direction
        logical                                :: Nulldif          = .false.
        logical                                :: NumericStability = .false.
        real                                   :: VolumeRelMax
        integer                                :: AdvMethodH, TVDLimitationH
        integer                                :: AdvMethodV, TVDLimitationV
        logical                                :: Upwind2H, Upwind2V  
        logical                                :: Adv_Dif_Explicit                      
        !--For both models use
        real                                   :: Molecular_Diff_Coef 
                           
    end type T_AdvectionDiffusion

    type       T_Partition
        character (LEN = StringLength)          :: Couple
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

    type T_Evolution
        logical                                 :: Variable = .false.
        real                                    :: DTInterval
        type(T_Time)                            :: LastCompute
        type(T_Time)                            :: NextCompute
        logical                                 :: SoilQuality
        logical                                 :: SoilChemistry
        logical                                 :: Partitioning
        logical                                 :: CationExchangeProcess
        logical                                 :: ChemEquilibriumProcess
        logical                                 :: AdvectionDiffusion
        logical                                 :: SoilWaterFluxes
        logical                                 :: Macropores
        logical                                 :: MinConcentration
        type (T_AdvectionDiffusion)             :: AdvDiff
        type (T_Partition                    )  :: Partition
    end type T_Evolution
    
    type T_ThetaAtFaces
        real, dimension(:,:,:), pointer         :: ThetaW
        real, dimension(:,:,:), pointer         :: ThetaV
        real, dimension(:,:,:), pointer         :: ThetaU    
    end type T_ThetaAtFaces
    
    type T_MassBalance
        real(8)                                 :: TotalStoredMass
        real(8)                                 :: TranspiredMass
        real(8)                                 :: DNExchangeMass
        real(8)                                 :: RPExchangeMass
    end type T_MassBalance
    
    type T_Files
        character(PathLength)                   :: InitialFile
        character(PathLength)                   :: DataFile
        character(PathLength)                   :: FinalFile
        character(PathLength)                   :: TransientHDF
        character(PathLength)                   :: DataSedimentQualityFile
        character(PathLength)                   :: ChainReactionsDataFile
        integer                                 :: AsciiUnit        
    end type T_Files    

    type T_Property
        type (T_PropertyID)                     :: ID
        real, dimension(:,:,:), pointer         :: Concentration            => null()
        real, dimension(:,:,:), pointer         :: ConcentrationOld         => null()
        real, dimension(:,:), pointer           :: ConcentrationOnInfColumn      => null()
        real, dimension(:,:), pointer           :: ConcentrationDN               => null()
        real, dimension(:,:,:), pointer         :: ConcInInterfaceDN             => null()
        real, pointer, dimension(:,:,:)         :: Mass_Created
        real(8)                                 :: TotalStoredMass
        real, pointer, dimension(:,:,:)         :: ViscosityU
        real, pointer, dimension(:,:,:)         :: ViscosityV
        type (T_Property), pointer              :: Next, Prev                     => null()
        logical                                 :: Particulate
        type (T_Evolution)                      :: Evolution
        type (T_MassBalance)                    :: MB
        real, pointer, dimension(:,:,:)         :: Diffusivity
        real, pointer, dimension(:,:,:)         :: Diff_Turbulence_H
        real, pointer, dimension(:,:,:)         :: Diff_Turbulence_V
        real, pointer, dimension(:,:,:)         :: Viscosity
        logical                                 :: Old     = .false.
        real                                    :: MinValue        = FillValueReal
        logical                                 :: TimeSerie        = .false.
        logical                                 :: BoxTimeSerie     = .false.
        logical                                 :: BoxTimeSerie2D   = .false.
        logical                                 :: OutputHDF        = .false.       
    end type T_Property

    type T_Coupled
        logical                                 :: SoilQuality          = .false. !Sediment source/sink model (Sediment Quality)
        real                                    :: SoilQuality_DT
        type (T_Time)                           :: SoilQuality_NextCompute

#ifdef _PHREEQC_        
        logical                                 :: SoilChemistry        = .false.  !Chemical reactions model (PhreeqC)
        real                                    :: SoilChemistry_DT
        type (T_Time)                           :: SoilChemistry_NextCompute
#endif        

        logical                                 :: ChainReactions       = .false. !Simple generic reactions ZERO & FIRST order rate

        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: Partition            = .false.
        logical                                 :: MinConcentration     = .false.
    end type T_Coupled
    
    !Implicit coef for thomas matrix
    type       T_D_E_F
        real   , pointer, dimension(: , : , :)  :: D
        real(8), pointer, dimension(: , : , :)  :: E
        real   , pointer, dimension(: , : , :)  :: F
    end type T_D_E_F
    
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
        real, pointer, dimension(: , : , :)  :: CoefInterfRunoff
        real, pointer, dimension(: , : , :)  :: CoefInterfDN
        real, pointer, dimension(: , : , :)  :: CoefInterfVeg
    end type T_A_B_C_Explicit

    type       T_FluxCoef
        real   , pointer, dimension(: , : , :)  :: C_flux    !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: D_flux    !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: E_flux    !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: F_flux    !Coeficient to calculate AdvFlux and DifFlux
    end type T_FluxCoef

    type T_PorousMediaProperties
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjAdvectionDiffusion = 0
        integer                                     :: ObjBasinGeometry     = 0
        integer                                     :: ObjPorousMedia       = 0
        integer                                     :: ObjGeometry          = 0
        integer                                     :: ObjMap               = 0
        integer                                     :: ObjGridData          = 0
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjtimeSerie         = 0
        integer                                     :: ObjSedimentQuality   = 0
        integer                                     :: ObjHDF5              = 0
        integer                                     :: ObjBottomTopography  = 0
        integer                                     :: ObjProfile           = 0
        integer                                     :: ObjInterface         = 0
        integer                                     :: ObjChainReactions         = 0 !ChainReactionsID
        
#ifdef _PHREEQC_        
        integer                                     :: ObjPhreeqC                = 0
        integer                                     :: ObjInterfaceSoilChemistry = 0 
        real,    pointer, dimension(:,:,:)          :: CellSoilMass
        real,    pointer, dimension(:,:,:)          :: CellWaterMass          
        type (T_Property), pointer                  :: PropertySoilDryDensity
#endif        
        real,    pointer, dimension(:,:,:)          :: CellWaterVolume            !Used by SoilChemistry & ChainReactions modules        
        integer, pointer, dimension(:)              :: PropertiesList             !List with ID of all properties used 
 
        type (T_ExtVar)                             :: ExtVar
        type (T_Files)                              :: Files
        type (T_OutPut)                             :: OutPut
        type (T_Property), pointer                  :: FirstProperty    => null() 
        type (T_Property), pointer                  :: LastProperty        
        type (T_PorousMediaProperties), pointer     :: Next             => null()
        type (T_Coupled)                            :: Coupled
        type (T_Time)                               :: LastOutputHDF5
        type(T_D_E_F)                               :: COEF3
        type(T_A_B_C_Explicit)                      :: COEFExpl 
        type(T_FluxCoef)                            :: COEF3_VertAdv            !Vertical advection coeficients
        type(T_FluxCoef)                            :: COEF3_HorAdvXX           !Horizont advection coeficients
        type(T_FluxCoef)                            :: COEF3_HorAdvYY           !Horizont advection coeficients
        real, pointer, dimension(: , : , :)         :: TICOEF3               
        real(8), pointer, dimension(:)              :: VECG                     !Auxiliar thomas arrays 
        real(8), pointer, dimension(:)              :: VECW                     !Auxiliar thomas arrays     

        logical                                     :: PorousMediaProperties
        real,    pointer, dimension(:,:,:)          :: Volume   
        integer                                     :: PropertiesNumber    = 0
        real   , pointer, dimension(:,:,:)          :: DissolvedToParticulate3D
        real                                        :: ResidualTime
        
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D

        type(T_Property_3D)                         :: Disper_Trans
        type(T_Property_3D)                         :: Disper_Longi
               
        logical                                     :: AdvDiff_Explicit            ! True - All explicit; False - Vertical Implicit
        logical                                     :: AdvDiff_AdvectionH_ImpExp   ! True - Horiz Adv Expl; False - Horiz Adv Impl
        logical                                     :: AdvDiff_CheckCoefs    !
        logical                                     :: Vertical1D
        logical                                     :: XZFlow
        
        type(T_ThetaAtFaces)                        :: ThetaAtFaces
        
        !--For PorousMediaProperties Advection-Diffusion Method
        real,    pointer, dimension(:,:,:)          :: DifusionNumber
        real,    pointer, dimension(:,:,:)          :: ReynoldsMNumber    
        
        logical                                     :: CheckGlobalMass       
        
        real(8), pointer, dimension(:,:,:)          :: WaterVolume
        real(8), pointer, dimension(:,:,:)          :: FluxWCorr             

    end type  T_PorousMediaProperties

    !Global Module Variables
    type (T_PorousMediaProperties), pointer         :: FirstObjPorousMediaProperties
    type (T_PorousMediaProperties), pointer         :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructPorousMediaProperties(ObjPorousMediaPropertiesID,                 &
                                              ComputeTimeID,                              &
                                              HorizontalGridID,                           &
                                              HorizontalMapID,                            &
                                              BasinGeometryID,                            &
                                              PorousMediaID,                              &
                                              GeometryID,                                 &
                                              MapID,                                      &
                                              CoupledDN,                                  &
                                              CheckGlobalMass,                            &
                                              STAT)
     
        !Arguments---------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: BasinGeometryID
        integer                                         :: PorousMediaID
        integer                                         :: GeometryID
        integer                                         :: MapID
        logical, optional                               :: CoupledDN
        logical                                         :: CheckGlobalMass
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
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjPorousMedia    = AssociateInstance (mPOROUSMEDIA_,    PorousMediaID   )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMap_,            MapID           )
        

            if (present(CoupledDN)) then
                Me%ExtVar%CoupledDN = CoupledDN
            endif
            Me%CheckGlobalMass       = CheckGlobalMass
            
            
            call ReadFileNames


            !Constructs the DataFile
            call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR01'                
           
            call ReadGlobalOptions

            call Construct_PropertyList

            if (Me%Coupled%Partition) then
                call ConstructPartition
            end if     
            
            call AllocateVariables
       
            call ConstructHDF    
    
            call ConstructTimeSerie
            

            
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
    
    subroutine CheckFlowDirections

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: ComputeFacesX, ComputeFacesY  
        logical                                     :: YFlow, XFlow, ZFlow  
        !Begin-----------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesW3D = Me%ExtVar%ComputeFacesW3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckDirections - ModulePorousMediaProperties - ERR120'


        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesU3D = Me%ExtVar%ComputeFacesU3D,             &
                               ComputeFacesV3D = Me%ExtVar%ComputeFacesV3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckDirections - ModulePorousMediaProperties - ERR130'

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckDirections - ModulePorousMediaProperties - ERR0140'

        
        !1D flow in Z direction
        Me%Vertical1D = .false.
        !2D flow in XZ direction
        Me%XZFlow     = .false.
        
        
        do I = ILB, IUB
        do J = JLB, JUB
        do K = KLB, KUB
            
            if (Me%ExtVar%WaterPoints3D(i, j, k) == WaterPoint) then     
                
                ComputeFacesX = 0
                ComputeFacesY = 0    
                
                if (Me%ExtVar%ComputeFacesU3D(i,j,k) == 1) then
                    ComputeFacesX = ComputeFacesX + 1
                    XFlow         = .true.
                endif
                if (Me%ExtVar%ComputeFacesV3D(i,j,k) == 1) then
                    ComputeFacesY = ComputeFacesY + 1
                    YFlow         = .true.
                endif
                if (Me%ExtVar%ComputeFacesW3D(i,j,k) == 1) then
                    ZFlow         = .true.
                endif

                !If any cell has flow in two directions than is 3D solution
                if (ComputeFacesX + ComputeFacesY .eq. 2) then
                    Me%Vertical1D = .false.
                    Me%XZFlow     = .false.
                    return
                endif
                
            endif
        
        enddo
        enddo
        enddo   
        
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

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesW3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CheckDirections - ModulePorousMediaProperties - ERR170'
        
        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesU3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CheckDirections - ModulePorousMediaProperties - ERR180'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesV3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CheckDirections - ModulePorousMediaProperties - ERR190'

        
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
 !           write(*, *)"---Discharges       : ", CurrentProperty%Evolution%Discharges
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
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('POROUS_PROP_HDF', Me%Files%TransientHDF, "PorousMedia HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01b'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('POROUS_PROP_FIN', Me%Files%FinalFile, "PorousMedia Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01c'
   
        !Reads the name of the file where to store final data
        call ReadFileName ('POROUS_PROP_INI', Me%Files%InitialFile, "PorousMedia Initial File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01d'

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
            stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR90'
        
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
                stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR95'
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

    subroutine AllocateVariables        
        
        !Local-----------------------------------------------------------------        
        integer                                         :: ILB, IUB, JLB,  JUB 
        integer                                         :: KLB, KUB, IJKLB, IJKUB 

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
        allocate (Me%CellWaterMass (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%CellSoilMass (ILB:IUB,JLB:JUB,KLB:KUB))
        
        Me%CellWaterMass = 0.
        Me%CellSoilMass  = 0.                
#endif

        if (Me%Coupled%AdvectionDiffusion) then
            allocate (Me%WaterVolume          (ILB:IUB,JLB:JUB,KLB:KUB))
            allocate (Me%FluxWCorr                (ILB:IUB,JLB:JUB,KLB:KUB))
            
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
            allocate (Me%COEFExpl%CoefInterfVeg    (ILB:IUB,JLB:JUB,KLB:KUB))

            allocate(Me%TICOEF3                 (ILB:IUB, JLB:JUB, KLB:KUB))

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
                
                allocate(Me%COEF3%D                 (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3%E                 (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(Me%COEF3%F                 (ILB:IUB, JLB:JUB, KLB:KUB))
                
                IJKLB = min (ILB, JLB, KLB)
                IJKUB = max (IUB, JUB, KUB)
                allocate(Me%VECG                    (IJKLB:IJKUB))
                allocate(Me%VECW                    (IJKLB:IJKUB))

                Me%COEF3%D                  = 0.0
                Me%COEF3%E                  = 1.0
                Me%COEF3%F                  = 0.0
                Me%TICOEF3                  = 0.0
                Me%VECG                     = Null_real 
                Me%VECW                     = Null_real 
           
            endif
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
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR01'

        if (NewProperty%Particulate)then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR03'
            end if
        endif
        

    end subroutine Construct_PropertyState

    !--------------------------------------------------------------------------

    subroutine Construct_PropertyEvolution(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        real                                        :: ErrorAux, AuxFactor, DTAux
        real                                        :: ModelDT
        !----------------------------------------------------------------------

        call GetData(NewProperty%Evolution%AdvectionDiffusion,                           &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'ADVECTION_DIFFUSION',                               &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR10'

        if (NewProperty%Evolution%AdvectionDiffusion .and. NewProperty%Particulate) then
            write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is Particulate'
            write(*,*) 'and can not have ADVECTION_DIFFUSION on'
            stop 'Construct_PropertyEvolution - ModulePorousMediaProeprties - ERR15'      
        endif  

        if (NewProperty%Evolution%AdvectionDiffusion) then
            Me%Coupled%AdvectionDiffusion = .true.
            NewProperty%Evolution%Variable = .true.
        endif
        
        if (NewProperty%Evolution%AdvectionDiffusion) then

            call ReadAdvectionDiffusionParam (NewProperty)
        
        end if

        call ConstructPropertyDiffusivity (NewProperty)


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
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR20'

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
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR30'

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
            stop 'Subroutine Construct_PropertyEvolution - ModuleRunoffProperties - ERR50'

        if(NewProperty%Evolution%Partitioning) then
            NewProperty%Evolution%Variable = .true.
            Me%Coupled%Partition           = .true.
        endif
        
        if(NewProperty%Evolution%Partitioning)                                           &
            call Read_Partition_Parameters(NewProperty)

        !Property time step
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR40'

            ModelDT = Me%ExtVar%DT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DTINTERVAL',                                    &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModulePorousMediaProperties',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR050'
                                       
            
            if (NewProperty%Evolution%DTInterval < ModelDT) then
                write(*,*) 
                write(*,*) 'Property time step is smaller then model time step'
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR60'

            elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                !Property time step must be a multiple of the model time step
                auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) 'Property time step must be a multiple of model time step.'
                    write(*,*) 'Please review your input data.'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR70'
                endif

                !Run period in seconds
                DTaux = Me%ExtVar%EndTime - Me%ExtVar%Now

                !The run period   must be a multiple of the Property DT
                auxFactor = DTaux / NewProperty%Evolution%DTInterval

                ErrorAux = auxFactor - int(auxFactor)
                if (ErrorAux /= 0) then

                    write(*,*) 
                    write(*,*) 'Property time step is not a multiple of model time step.'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR80'
                end if
            endif

            NewProperty%Evolution%NextCompute = Me%ExtVar%Now + NewProperty%Evolution%DTInterval

        else

            call null_time(NewProperty%Evolution%NextCompute)

            NewProperty%Evolution%DTInterval = FillValueReal

        endif

    end subroutine Construct_PropertyEvolution     

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
        integer                   :: iflag, BoundaryCondition

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

        call GetData(BoundaryCondition,                            &
                     Me%ObjEnterData,  iflag,                      &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_BOUNDARY_CONDITION',  &
                     Default      = MassConservation,              &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR70'

        ! By default it's imposed a value dependent only from the exterior
        ! value and of the decay time. However this method doesn't conserve mass
        ! when the water fluxes near the frontier are dominant

        if (BoundaryCondition /= MassConservation     .and. &
            BoundaryCondition /= ImposedValue         .and. &
            BoundaryCondition /= NullGradient         .and. &
            BoundaryCondition /= CyclicBoundary       .and. &
            BoundaryCondition /= Orlanski             .and. &
            BoundaryCondition /= MassConservNullGrad) &
            stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR80'

        NewProperty%Evolution%AdvDiff%BoundaryCondition = BoundaryCondition

        !By default the horizontal Diffusion discretization is explicit
        NewProperty%Evolution%AdvDiff%DiffusionH_imp_exp  = ExplicitScheme

        NewProperty%Evolution%AdvDiff%ImplicitH_Direction = DirectionX
        
        !Spatial Method for horizontal advection
        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodH,     &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ADVDIFF_METHOD_H',            &
                     Default      = UpwindOrder1,                  &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR90'

        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationH, &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ADVDIFF_TVD_LIMIT_H',         &
                     Default      = Superbee,                      &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR100'

        !Spatial Method for horizontal advection
        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodV,     &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ADVDIFF_METHOD_V',            &
                     Default      = UpwindOrder1,                  &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR110'

        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationV, &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ADVDIFF_TVD_LIMIT_V',         &
                     Default      = Superbee,                      &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR120'


        call GetData(NewProperty%Evolution%AdvDiff%VolumeRelMax,   &
                     Me%ObjEnterData, iflag,                       &
                     Keyword      = 'ADVDIFF_VOLUME_RELATION_MAX', &
                     Default      = 5.,                            &
                     SearchType   = FromBlock,                      &
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
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR10'
        NewProperty%Diffusivity       = 0. 
        
           
        allocate (NewProperty%Viscosity (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR20'
        NewProperty%Viscosity         = 0.
        
!        allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR30'
!        NewProperty%Diff_Turbulence_H = 0.
!
!        allocate (NewProperty%Diff_Turbulence_V (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR40'
!        NewProperty%Diff_Turbulence_V = 0.                   
        
        allocate (NewProperty%ViscosityU (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR60'
        
        NewProperty%ViscosityU  = 0.0
        
        allocate (NewProperty%ViscosityV (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR70'
    
        NewProperty%ViscosityV  = 0.0    
            

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

        allocate(NewProperty%Concentration(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR10'
        NewProperty%Concentration(:,:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOld(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR20'
        NewProperty%ConcentrationOld(:,:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOnInfColumn(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR40'
        NewProperty%ConcentrationOnInfColumn(:,:) = 0.

        
        if (Me%ExtVar%CoupledDN) then
            allocate(NewProperty%ConcentrationDN(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR50'
            NewProperty%ConcentrationDN(:,:) = FillValueReal
            

        endif
        
        !it has to be always allocated
        allocate(NewProperty%ConcInInterfaceDN(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR50'
        NewProperty%ConcInInterfaceDN(:,:,:) = FillValueReal            
        
        call GetData(NewProperty%MinValue,                                                  &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'MIN_VALUE',                                            &
                     ClientModule = 'ModuleSoilProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR100'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
            Me%Coupled%MinConcentration = .true.
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        if(NewProperty%Evolution%MinConcentration)then
            allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)&
                stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR110'
            NewProperty%Mass_Created(:,:,:) = 0.
        endif

        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleVegetation',                                 &
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
                                       PointsToFill3D       = Me%ExtVar%WaterPoints3D,     &
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
                    write(*,*) 'The property name is not recognised by the model'
                    stop       'ConstructPartition - ModuleRunoffProperties - ERR00'
                else
                    
                    DissolvedProperty%Evolution%Partition%Couple_ID = Couple_ID

                end if
                
                !Search for the particulate
                call Search_Property(ParticulateProperty, PropertyXID = Couple_ID, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ConstructPartition - ModuleWaterProperties - ERR10'
                
                !Check if coupled names of each other match
                if(trim(ParticulateProperty%Evolution%Partition%Couple) .ne.                &
                   trim(DissolvedProperty%ID%Name))                                         &
                   stop 'ConstructPartition - ModuleWaterProperties - ERR20'

                if(ParticulateProperty%Evolution%DTInterval .ne.                            &
                   DissolvedProperty%Evolution%DTInterval)                                  &
                    stop 'ConstructPartition - ModuleWaterProperties - ERR30'

                if(DissolvedProperty%Evolution%Partition%Rate /=                            &
                   ParticulateProperty%Evolution%Partition%Rate   )then
                    write(*,*)'Particulate and dissolved phases must have equal partition rates'
                    stop 'ConstructPartition - ModuleWaterProperties - ERR40'
                end if

                TotalPartition = DissolvedProperty%Evolution%Partition%Fraction  +          &
                                 ParticulateProperty%Evolution%Partition%Fraction

                Error = abs(1. - TotalPartition)
                     
                if(Error > 0.001)then
                    write(*,*)'Particulate and dissolved phases fractions must sum iqual to 1.'
                    stop 'ConstructPartition - ModuleWaterProperties - ERR50'
                end if

                ! .EQV. = Logical equivalence: the expression is true if both A and B 
                !  are true, or both are false.  
!                if (.NOT. (DissolvedProperty%Evolution%Partition%SalinityEffect .EQV.       &
!                           ParticulateProperty%Evolution%Partition%SalinityEffect))         &
!                    stop 'ConstructPartition - ModuleWaterProperties - ERR60'

!                if (DissolvedProperty%Evolution%Partition%EmpiricCoef .ne.                  &
!                    ParticulateProperty%Evolution%Partition%EmpiricCoef )                   &
!                    stop 'ConstructPartition - ModuleWaterProperties - ERR70'

                if (.NOT. (DissolvedProperty%Evolution%Partition%UseSedimentRefConc .EQV.   &
                           ParticulateProperty%Evolution%Partition%UseSedimentRefConc))     &
                    stop 'ConstructPartition - ModuleWaterProperties - ERR80'

                if(DissolvedProperty%Evolution%Partition%SedimentRefConc /=                 &
                   ParticulateProperty%Evolution%Partition%SedimentRefConc   )then
                    write(*,*)'Particulate and dissolved phases must have equal cohesive sediment'
                    write(*,*)'reference concentration'
                    stop 'ConstructPartition - ModuleWaterProperties - ERR90'
                end if


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
                            nProperties     = Me%PropertiesNumber ,                                        &
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
        !Begin------------------------------------------------------------------
        
        !Counts the number of Properties which has timeserie option set to true (2x for inf col concentration)
        PropertyX => Me%FirstProperty
        nProperties = 0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                nProperties = nProperties + 4
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates PropertyList
        allocate(PropertyList(nProperties))
        
        !Property names
        n=1
        PropertyX  => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                PropertyList(n)  = trim(PropertyX%ID%Name) //' [mg/l]'
                n=n+1
                PropertyList(n)  = trim(PropertyX%ID%Name) //'_in_InfilColumn [mg/l]'
                n=n+1
                PropertyList(n)  = trim(PropertyX%ID%Name) //'_Diffusivity [m2/s]'
                n=n+1
                PropertyList(n)  = trim(PropertyX%ID%Name) //'_DiffusivitySurface [m2/s]'
                n=n+1                
            endif
            PropertyX=>PropertyX%Next
        enddo

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMediaProperties',                      &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMediaProperties - ERR01' 

        if (iflag == 1) then
            Me%OutPut%TimeSerie_ON = .true.
        else
            Me%OutPut%TimeSerie_ON = .false.
        endif
        
        !Get water points
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR03.7'


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "spp",                                        &
                            WaterPoints3D = Me%ExtVar%WaterPoints3D,                    &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMediaProperties - ERR02' 

        !Unget
        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR085'

        !Deallocates PropertyList
        deallocate(PropertyList)
       
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

        endif

    end subroutine ConstructHDF


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
        
        call GetGridData (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
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

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
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


        inquire (FILE=trim(Me%Files%InitialFile)//"5", EXIST = Exist)

cd0:    if (Exist) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%InitialFile)//"5",&
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

            call HDF5ReadData   (ObjHDF5, "/Concentration/"//NewProperty%ID%Name,        &
                                 NewProperty%ID%Name,                                    &
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
        type(T_Property), pointer                           :: PropertyX
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
        allocate(Me%DissolvedToParticulate3D(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB, &
                 Me%WorkSize%KLB:Me%WorkSize%KUB))
        Me%DissolvedToParticulate3D(:,:,:) = null_real

        Me%ResidualTime = 0.
    
    end subroutine CoupleSoilQuality

    !--------------------------------------------------------------------------


#ifdef _PHREEQC_
    !--------------------------------------------------------------------------
    subroutine CoupleSoilChemistry        

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        integer, pointer, dimension(:)                      :: SoilChemistryPropertyList
        integer                                             :: STAT_CALL
        real                                                :: SoilChemistryDT
        integer                                             :: nProp = 0 

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
        
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR01'

        !Start Interface
        call ConstructInterface(InterfaceID         = Me%ObjInterfaceSoilChemistry,  &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = PhreeqCModel,                  &
                                DT                  = SoilChemistryDT,               &
                                PropertiesList      = SoilChemistryPropertyList,     &
                                WaterPoints3D       = Me%ExtVar%WaterPoints3D,       &
                                Size3D              = Me%WorkSize,                   &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR02'

        call SetupSoilChemistry  

        call UnGetMap (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR03'

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
                Me%PropertySoilDryDensity => PropertyX  
                
                !Done to ensure that this property will not be sended to PhreeqC
                Me%PropertySoilDryDensity%Evolution%SoilChemistry = .false.
                
                found = .true.
                    
            endif

            PropertyX => PropertyX%Next
        enddo

        if (.not. found) then !stop 'PhreeqC needs property soil dry density - PorousMediaProperties - SetupSoilChemistry - ERR04'
        
            Me%PropertySoilDryDensity => null()            
        
        end if
        
!        end if
    
    end subroutine SetupSoilChemistry
   
    !--------------------------------------------------------------------------
    
#endif

    !--------------------------------------------------------------------------
    subroutine CoupleChainReactions
        
        !Local-----------------------------------------------------------------
        integer :: STAT
        integer :: PropertiesCount

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
            stop 'PorousMediaProperties - CoupleChainReaction - Err010'
            
        if (PropertiesCount .EQ. 0) then
            
            call KillChainReactions(Me%ObjChainReactions, STAT = STAT)
            if (STAT /= SUCCESS_) &
                stop 'PorousMediaProperties - CoupleChainReaction - Err020'
                
            Me%Coupled%ChainReactions = .false.
            
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

    subroutine SetVegetationPMProperties(PorousMediaPropertiesID,                  &
                                         SoilFluxesActive,                         &
                                         GrazingBiomass,                           &
                                         GrazingNitrogen,                          &
                                         GrazingPhosphorus,                        &
                                         ManagementAerialBiomass,                  &
                                         ManagementNitrogen,                       &
                                         ManagementPhosphorus,                     &
                                         ManagementRootBiomass,                    &
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
                                         Management,                               &
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
                                         STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        logical, dimension(:,:), pointer, optional  :: SoilFluxesActive
        real, dimension(:,:), pointer, optional     :: GrazingBiomass
        real, dimension(:,:), pointer, optional     :: GrazingNitrogen
        real, dimension(:,:), pointer, optional     :: GrazingPhosphorus
        real, dimension(:,:), pointer, optional     :: ManagementAerialBiomass
        real, dimension(:,:), pointer, optional     :: ManagementNitrogen
        real, dimension(:,:), pointer, optional     :: ManagementPhosphorus
        real, dimension(:,:), pointer, optional     :: ManagementRootBiomass
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
        logical,  optional                          :: Management
        logical,  optional                          :: Dormancy
        logical,  optional                          :: Fertilization
        logical,  optional                          :: NutrientFluxesWithSoil
        logical,  optional                          :: ModelWater
        logical,  optional                          :: ModelNitrogen
        logical,  optional                          :: ModelPhosphorus
        logical,  optional                          :: GrowthModel
        logical,  optional                          :: CoupledVegetation
        real,     optional                          :: VegetationDT
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
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
            if (present(Management)) then
                Me%ExtVar%Management               = Management
                if (present(ManagementAerialBiomass )) Me%ExtVar%ManagementAerialBiomass  => ManagementAerialBiomass
                if (present(ManagementNitrogen      )) Me%ExtVar%ManagementNitrogen       => ManagementNitrogen
                if (present(ManagementPhosphorus    )) Me%ExtVar%ManagementPhosphorus     => ManagementPhosphorus
                if (present(ManagementRootBiomass   )) Me%ExtVar%ManagementRootBiomass    => ManagementRootBiomass
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
                                           WaterColumn,                         &
                                           STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaPropertiesID
        real(8), dimension(:,:), pointer            :: WaterColumn
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_,STAT_CALL
        type(T_Property), pointer                   :: PropertyX
        !----------------------------------------------------------------------

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

            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))
                
                call SetMatrixValue (PropertyX%ConcentrationOld, Me%Size, PropertyX%Concentration, Me%ExtVar%WaterPoints3D)
                
                if (Me%CheckGlobalMass) then
                    PropertyX%MB%TotalStoredMass     = 0.0
                    PropertyX%MB%TranspiredMass      = 0.0
                    PropertyX%MB%DNExchangeMass      = 0.0
                    PropertyX%MB%RPExchangeMass      = 0.0
                endif
                
                PropertyX => PropertyX%Next
                
            end do                        
            
            !Actualize properties if evolution from file
            call ActualizePropertiesFromFile
            
            !Nutrient sources from vegetation. Sinks are explicit or implict cared in transport
            call InterfaceFluxes

            if (Me%Coupled%AdvectionDiffusion) then
                call AdvectionDiffusionProcesses
            endif

            if (Me%Coupled%SoilQuality) then
                call SoilQualityProcesses
            endif

#ifdef _PHREEQC_
            if (Me%Coupled%SoilChemistry) then
                call SoilChemistryProcesses
            endif
#endif            

            if (Me%Coupled%ChainReactions) then
                call ChainReactionsProcesses
            endif

            if (Me%Coupled%Partition)then
                call Partition_Processes 
            endif           

            if (Me%Coupled%MinConcentration) then
                call SetLimitsConcentration 
            endif

            if (Me%Output%Timeserie_ON) then
                call OutPut_TimeSeries
            endif

            if (Me%Output%HDF_ON) then
                call OutPut_HDF
            endif

    !       call ProfileOutput    em teste no construct e no kill

            call Actualize_Time_Evolution

            if (Me%CheckGlobalMass) then
                call CalculateTotalStoredMass
            endif
        
            call ReadUnlockExternalVar

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyPorousMediaProperties

    !-----------------------------------------------------------------------------

    subroutine ComputeVolumes

        !Local-----------------------------------------------------------------
        integer :: i, j, k, CHUNK        

        !----------------------------------------------------------------------
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        !$OMP PARALLEL PRIVATE(I,J,K)
        
        !Compute volumes and correct top flux taking FluxW(KUB+1) because it would be interpreted by module advection diffusion
        !as an additional water flux with the conc of C(i,j,k)
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,k) == WaterPoint) then             
                 
                Me%WaterVolume(i,j,k)        = Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)

            endif
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL 
               
        !Correct fluxw - take FluxW(KUB+1) because it would be interpreted by module advection diffusion
        !as an additional water flux with the conc of C(i,j,k)
         
         CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(I,J,K)

        k = Me%WorkSize%KUB
        call SetMatrixValue (Me%FluxWCorr, Me%Size, Me%ExtVar%FluxW)
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
    
        !Local-----------------------------------------------------------------
        integer                         :: I,J,K, CHUNK
        real, pointer, dimension(:,:,:) :: ThetaOld, Theta
        !Begin-----------------------------------------------------------------
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K)

        ThetaOld  => Me%ExtVar%WaterContentOld
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
            
    end subroutine ComputeThetaAtFaces

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
            if (Me%ExtVar%ComputeVegInterfaceFluxes) then
                call VegetationInterfaceFluxes
            endif
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
    ! This routine solves mass sources due to vegetation. 
    
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
        real                                        :: ManagementNotCarbon, ManagementNitrogen, ManagementPhosphorus
        real                                        :: ManagementAerialBiomass, ManagementCarbon
        real                                        :: ManagementRootNotCarbon, ManagementRootNitrogen, ManagementRootPhosphorus
        real                                        :: ManagementRootBiomass, ManagementRootCarbon
        real                                        :: NitrogenFraction, PhosphorusFraction, RootDistribution
        real                                        :: FertilizationAmmonia, FertilizationNitrate
        real                                        :: FertilizationOrganicN, FertilizationOrganicP
        real                                        :: FertilizationMineralP
!        real                                        :: NitrogenUptake, PhosphorusUptake
        real                                        :: ModelDT, VegDT, CellWaterVolume, CellSoilMass
        type (T_Property), pointer                  :: Property
        integer                                     :: STAT_CALL

        !Begin--------------------------------------------------------------------

                
        !!CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR01'  


        !!!$OMP PARALLEL PRIVATE(I,J,K)
        !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
        if (Me%ExtVar%BasinPoints(i,j) == BasinPoint .and. Me%ExtVar%SoilFluxesActive(i,j)) then
            
            Area         = Me%ExtVar%Area(i,j)
            RootDepth    = Me%ExtVar%RootDepth(i,j)
            FoundEnd     = .false.
            BottomDepth  = 0.
      
do3:        do K = Me%WorkSize%KUB, Me%WorkSize%KLB, -1                
            
                
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

                GrazingNotCarbon  = 0.0
                GrazingNitrogen   = 0.0
                GrazingPhosphorus = 0.0
                DormancyNotCarbon  = 0.0
                DormancyNitrogen   = 0.0
                DormancyPhosphorus = 0.0
                ManagementNotCarbon  = 0.0
                ManagementNitrogen   = 0.0
                ManagementPhosphorus = 0.0                
                FertilizationAmmonia    = 0.0
                FertilizationNitrate    = 0.0
                FertilizationOrganicN   = 0.0
                FertilizationOrganicP   = 0.0
                FertilizationMineralP   = 0.0
                ManagementRootCarbon     = 0.0
                ManagementRootNotCarbon  = 0.0
                ManagementRootNitrogen   = 0.0
                ManagementRootPhosphorus = 0.0                
                
                if (Me%ExtVar%GrowthModel) then

                    !Fluxes only occuring in surface (aerial biomass residue from grazing, dormancy and management; 
                    !surface fertilization)
                    if (k == Me%WorkSize%KUB) then
                
                        !Grazing
                        GrazingNotCarbon  = 0.0
                        GrazingNitrogen   = 0.0
                        GrazingPhosphorus = 0.0
                        if (Me%ExtVar%Grazing) then
                
                            if (Me%ExtVar%ModelNitrogen) then
                            
                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                GrazingNitrogen  = Me%ExtVar%GrazingNitrogen(i,j) * 1e9 * Area / 10000.
                            
                                GrazingNotCarbon = GrazingNotCarbon + GrazingNitrogen
                
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                GrazingPhosphorus = Me%ExtVar%GrazingPhosphorus(i,j) * 1e9 * Area / 10000.

                                GrazingNotCarbon  = GrazingNotCarbon + GrazingPhosphorus
                
                            endif                          
                
                            !      ug       = Kg/ha * 1E9mg/kg * (m2) * 1ha/10000m2                     
                            GrazingBiomass = Me%ExtVar%GrazingBiomass(i,j) * 1e9 * Area / 10000.
                
                            GrazingCarbon  = GrazingBiomass - GrazingNotCarbon

                        endif

        !                !Dormancy
                        DormancyNotCarbon  = 0.0
                        DormancyNitrogen   = 0.0
                        DormancyPhosphorus = 0.0
                        if (Me%ExtVar%Dormancy) then
                
                            if (Me%ExtVar%ModelNitrogen) then

                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                DormancyNitrogen  = Me%ExtVar%DormancyNitrogen(i,j) * 1e9 * Area / 10000.

                                DormancyNotCarbon = DormancyNotCarbon + DormancyNitrogen
                
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                DormancyPhosphorus = Me%ExtVar%DormancyPhosphorus(i,j) * 1e9 * Area / 10000.

                                DormancyNotCarbon  = DormancyNotCarbon + DormancyPhosphorus
                
                            endif                          
                
                            !      ug       = Kg/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                            DormancyBiomass = Me%ExtVar%DormancyBiomass(i,j) * 1e9 * Area / 10000.
                
                            DormancyCarbon  = DormancyBiomass - DormancyNotCarbon

                        endif


        !                !Management
                        ManagementNotCarbon  = 0.0
                        ManagementNitrogen   = 0.0
                        ManagementPhosphorus = 0.0
                        if (Me%ExtVar%Management) then
                
                            if (Me%ExtVar%ModelNitrogen) then

                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                ManagementNitrogen  = Me%ExtVar%ManagementNitrogen(i,j) * 1e9 * Area / 10000.

                                ManagementNotCarbon = ManagementNotCarbon + ManagementNitrogen
                
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                ManagementPhosphorus = Me%ExtVar%ManagementPhosphorus(i,j) * 1e9 * Area / 10000.

                                ManagementNotCarbon  = ManagementNotCarbon + ManagementPhosphorus
                
                            endif                          
                
                            ManagementAerialBiomass = Me%ExtVar%ManagementAerialBiomass(i,j) * 1e9 * Area / 10000.
                
                            !      ug       = Kg/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                            ManagementCarbon  = ManagementAerialBiomass - ManagementNotCarbon

                        endif

        !                !Fertilization in Surface
                        FertilizationAmmonia    = 0.0
                        FertilizationNitrate    = 0.0
                        FertilizationOrganicN   = 0.0
                        FertilizationOrganicP   = 0.0
                        FertilizationMineralP   = 0.0
                        if (Me%ExtVar%Fertilization) then
                

                            if (Me%ExtVar%ModelNitrogen) then
                    
                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                FertilizationNitrate  = Me%ExtVar%FertilNitrateSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationAmmonia  = Me%ExtVar%FertilAmmoniaSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationOrganicN = Me%ExtVar%FertilOrganicNSurface(i,j) * 1e9 * Area / 10000.
                    
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                FertilizationOrganicP = Me%ExtVar%FertilOrganicPSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationMineralP = Me%ExtVar%FertilMineralPSurface(i,j) * 1e9 * Area / 10000.
                    
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
                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2 
                                FertilizationOrganicN = Me%ExtVar%FertilOrganicNSubSurface(i,j) * 1e9 * Area / 10000.
                    
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                FertilizationOrganicP = Me%ExtVar%FertilOrganicPSubSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationMineralP = Me%ExtVar%FertilMineralPSubSurface(i,j) * 1e9 * Area / 10000.
                    
                            endif                          
                
                        endif

                    endif

                    !Root death to soil (occurrs in all layers until root end )
                    if (RootDepth .gt. 0.0) then
                              !      ug       = Kg/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                        ManagementRootBiomass = Me%ExtVar%ManagementRootBiomass(i,j) * 1e9 * Area / 10000.
                        !Root distribution (based on SWAT formulation)
                        RootDistribution = (1.0 - exp(-10.0 * BottomDepth/RootDepth))/(1.0 - exp(-10.0))                   &
                                           - (1.0 - exp(-10.0 * TopDepth/RootDepth))/(1.0 - exp(-10.0))

                        ManagementRootNotCarbon  = 0.0
                        ManagementRootNitrogen   = 0.0
                        ManagementRootPhosphorus = 0.0
                        if (Me%ExtVar%ModelNitrogen) then
                            NitrogenFraction        = Me%ExtVar%NitrogenFraction (i,j)
                            ManagementRootNitrogen  = ManagementRootBiomass * NitrogenFraction * RootDistribution
                            ManagementRootNotCarbon = ManagementRootNotCarbon + ManagementRootNitrogen
                        endif

                        if (Me%ExtVar%ModelPhosphorus) then
                            PhosphorusFraction        = Me%ExtVar%PhosphorusFraction (i,j)
                            ManagementRootPhosphorus  = ManagementRootBiomass * PhosphorusFraction * RootDistribution
                            ManagementRootNotCarbon   = ManagementRootNotCarbon + ManagementRootPhosphorus
                        endif
                
                        !     ug
                        ManagementRootCarbon = ManagementRootBiomass - ManagementRootNotCarbon

                    endif
                endif


!                !Plant Uptake (occurrs in all layers until root end )
!                if (Me%ExtVar%ModelNitrogen) then
!                    
!                    !      gN       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2                     
!                    NitrogenUptake = Me%ExtVar%NitrogenUptake(i,j,k) * 1e3 * Area / 10000.
!                    
!                endif
!
!                if (Me%ExtVar%ModelPhosphorus) then
!                    
!                    !      gP       = KgP/ha * 1E3g/kg * (m2) * 1ha/10000m2                     
!                    PhosphorusUptake = Me%ExtVar%PhosphorusUptake(i,j,k) * 1e3 * Area / 10000.
!
!                endif

                
                !!NOW UPDATE CONCENTRATIONS WITH THE MASS FLUXES COMPUTED. ConcentrationOld used because done 
                !previous to transport that changes to new concentration
                
                !s
                ModelDT         = Me%ExtVar%DT
                VegDT           = Me%ExtVar%VegetationDT
                !m3             = m3H20/m3cell * m3cell
                CellWaterVolume = Me%ExtVar%WaterContentOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 

                !Soil mass to compute organic and microorganisms pools
                call SearchProperty(Property, SoilDryDensity_        , .false., STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
                !kgsoil         = kg/m3  * m3cell
                CellSoilMass    = Property%ConcentrationOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
                
                if (Me%ExtVar%GrowthModel) then

                    ! Property Calculation
                    !!Carbon
                    call SearchProperty (Property, RefreactaryOrganicC_        , .false., STAT = STAT_CALL)    
                    if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR80'
                
                    !         ug/kgsoil            = ug/kgsoil + ug / kgsoil
                    Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingCarbon + DormancyCarbon    &
                                                        + ManagementCarbon + ManagementRootCarbon) * ModelDT / VegDT)            &
                                                        / CellSoilMass)

    !                call SearchProperty (Property, LabileOrganicC_        , .false., STAT = STAT_CALL)    
    !                if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
    !
    !                Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicC)               &
    !                                                * ModelDT / VegDT) / CellSoilMass)
                endif

                !!Nitrogen
                if (Me%ExtVar%ModelNitrogen) then
                    
                    if (Me%ExtVar%GrowthModel) then
                        call SearchProperty (Property, RefreactaryOrganicN_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR100'
                    
                        !         ug/kgsoil 
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingNitrogen              &
                                                           + DormancyNitrogen + ManagementNitrogen                              &
                                                           + ManagementRootNitrogen) * ModelDT / VegDT)                         &
                                                           / CellSoilMass)


                        call SearchProperty (Property, PON_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR110'
                    
                        !         ug/kgsoil 
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicN)        &
                                                           * ModelDT / VegDT) / CellSoilMass)


                        call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR120'

                        !         g/m3                = g/m3 + g / m3H20
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationNitrate)          &
!                                                            - NitrogenUptake) * ModelDT / VegDT) / CellWaterVolume)
                                                             * ModelDT / VegDT) / CellWaterVolume)


!                        if (Me%CheckGlobalMass) then
!                            !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
!                            Property%MB%TranspiredMass = Property%MB%TranspiredMass + (NitrogenUptake                            &
!                                                                 * 1E-3 * ModelDT / VegDT)
!                        endif

                        call SearchProperty (Property, Ammonia_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR130'

                        !         g/m3                = g/m3 + g / m3H20 
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationAmmonia)         &
                                                           * ModelDT / VegDT) / CellWaterVolume)
                    else

!                        call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
!                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR120'
!
!                        !         g/m3                = g/m3 + g / m3H20 
!                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) - (((NitrogenUptake)               &
!                                                           * ModelDT / VegDT) / CellWaterVolume)

!                        if (Me%CheckGlobalMass) then
!                            !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
!                            Property%MB%TranspiredMass = Property%MB%TranspiredMass + (NitrogenUptake                            &
!                                                                *  1E-3 * ModelDT / VegDT)
!                        endif
                                                           
                    endif

                endif
                
                !Phosphorus
                if (Me%ExtVar%ModelPhosphorus) then

                    if (Me%ExtVar%GrowthModel) then

                        call SearchProperty (Property, RefreactaryOrganicP_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR140'

                        !         ug/kgsoil  
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingPhosphorus             &
                                                           + DormancyPhosphorus + ManagementPhosphorus                           &
                                                           + ManagementRootNitrogen) * ModelDT / VegDT)                          &
                                                           / CellSoilMass)


                        call SearchProperty (Property, POP_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR150'

                        !         ug/kgsoil  
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicP)        &
                                                           * ModelDT / VegDT) / CellSoilMass)


                        call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR160'

                        !         g/m3                    = g/m3 + g / m3H20  
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationMineralP)         &
!                                                            - PhosphorusUptake) * ModelDT / VegDT) / CellWaterVolume)   
                                                             * ModelDT / VegDT) / CellWaterVolume)   

!                        if (Me%CheckGlobalMass) then
!                            !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
!                            Property%MB%TranspiredMass = Property%MB%TranspiredMass + (PhosphorusUptake                          &
!                                                                *  1E-3 * ModelDT / VegDT)
!                        endif                                                            
                                                            
                    else
                        
                        call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR160'

!                        !         g/m3                    = g/m3 + g / m3H20  
!                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) - (((PhosphorusUptake)             &
!                                                           * ModelDT / VegDT) / CellWaterVolume) 
                                                           
                                                           
!                        if (Me%CheckGlobalMass) then
!                            !kg = [kg] + [g] * [1e-3kg/g] * [s model/s veg]
!                            Property%MB%TranspiredMass = Property%MB%TranspiredMass + (PhosphorusUptake                          &
!                                                                * 1E-3 * ModelDT / VegDT)
!                        endif                                                                                      
                     
                    endif             

                endif

            enddo do3

        endif
        
        enddo
        enddo

        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR1070'  



    end subroutine VegetationInterfaceFluxes

    !-----------------------------------------------------------------------------
    
    subroutine AdvectionDiffusionProcesses

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        
        !begin--------------------------------------------------------------------
        

       !Compute water volume and remove infiltration from fluxZ 
       !(it would create error in AdvectionDiffusion routines)
        call ComputeVolumes
        
        call ComputeThetaAtFaces

             
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
        call SetMatrixValue (Me%COEFExpl%CoefInterfVeg, Me%Size, 0.0)
        call SetMatrixValue (CurrProperty%ConcInInterfaceDN, Me%Size, 0.0)
        
        call SetMatrixValue (Me%TICOEF3, Me%Size, 0.0)

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
        
        if (.not. Me%AdvDiff_Explicit) then
            call SetMatrixValue (Me%COEF3%D, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3%E, Me%Size, dble(1.0))
            call SetMatrixValue (Me%COEF3%F, Me%Size, 0.0)

        endif
    
    end subroutine RestartVariables

    !-----------------------------------------------------------------------------

    subroutine ModifyCoefs(PropertyX)

        !Argument-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
       
        !begin--------------------------------------------------------------------
        
        !In this subroutine are computed all the sources/sinks of water and mass
        !that exist in the Porous Media Module (for new theta computation)
               
        !Diffusivity will be used in Infiltration coef and in AdvectionDiffusion
        call ModifyDiffusivity_New(PropertyX)

        !Fluxes in X, Y and Z direction explicit or implicit. 
        call ModifyAdvectionDiffusionCoefs(PropertyX)

        !Evaporation
        !Evaporation fluxes do not take mass and do not need to be accounted
        
        !Transpiration fluxes - in cells along roots - take all dissolved properties
        if (Me%ExtVar%CoupledVegetation .and. Me%ExtVar%ModelWater ) then
            call ModifyVegetationCoefs(PropertyX)
        endif
        
        !Infiltration - in surface cells
        call ModifyInfiltrationCoefs(PropertyX)
        
        !Fluxes with Drainage network - in the cells that link with river
        if (Me%ExtVar%CoupledDN) then
            call ModifyDrainageNetworkCoefs(PropertyX)
        endif

        
    end subroutine ModifyCoefs
    
    !-----------------------------------------------------------------------------
    
    subroutine ModifyAdvectionDiffusionCoefs(PropertyX)

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                       :: PropertyX
        real                                             :: ImpExp_AdvXX, ImpExp_AdvYY           
        real                                             :: AdvectionV_Imp_Exp  
        real                                             :: DiffusionV_Imp_Exp  
        real                                             :: AdvectionH_Imp_Exp                  
        !begin--------------------------------------------------------------------
        
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

        CHUNK = Chunk_K(Me%WorkSize%KLB,Me%WorkSize%KUB)
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

        CHUNK = Chunk_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
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

                    call THOMAS_3D(ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS, di, dj,    &
                         Me%COEF3%D,                                                    &
                         Me%COEF3%E,                                                    &
                         Me%COEF3%F,                                                    &
                         Me%TICOEF3,                                                    &
                         CurrProp%Concentration,                                        &
                         Me%VECG,                                                       &
                         Me%VECW)      
                
                else if (ImpExp_AdvYY == ImplicitScheme) then cd2

                    di = 1
                    dj = 0

                    call THOMAS_3D(JLBWS, JUBWS, ILBWS, IUBWS, KLBWS, KUBWS, di, dj,    &
                         Me%COEF3%D,                                                    &
                         Me%COEF3%E,                                                    &
                         Me%COEF3%F,                                                    &
                         Me%TICOEF3,                                                    &
                         CurrProp%Concentration,                                        &
                         Me%VECG,                                                       &
                         Me%VECW)      

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

        CHUNK = CHUNK_I(ILB, IUB)
        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxX,DT2,DT1)

            
        !CHUNK = CHUNK_K(Me%Size%KLB, Me%Size%KUB)

        !!$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,K)
        !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)

k1:     do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
i1:     do i = ILB, IUB


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
        
        !!$OMP END PARALLEL


        CHUNK = CHUNK_K(KLB, KUB)

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
        
        
        !CHUNK = CHUNK_K(Me%Size%KLB, Me%Size%KUB)
        CHUNK = CHUNK_J(JLB, JUB)
      
        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxY,DT2,DT1)

       
        !!$OMP PARALLEL SHARED(CHUNK) PRIVATE(J,K)
        !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)

k1:     do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
j1:     do j = JLB, JUB


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
        
        !!$OMP END DO NOWAIT
        !!$OMP END PARALLEL


cd6:    if (ImpExp_AdvYY == ExplicitScheme)  then !ExplicitScheme = 0

dok3 :      do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
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
            !$OMP END DO NOWAIT
            end do dok3

        else if (ImpExp_AdvYY == ImplicitScheme) then cd6 !ImplicitScheme = 1

dok4 :      do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
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
            !$OMP END DO NOWAIT
            end do dok4

        else cd6

            stop 'sub. HorizontalAdvectionYY - ModuleAdvectionDiffusion - ERR01'
        
        endif cd6

        !$OMP END PARALLEL

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
        
            CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
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

            CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
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

       
        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxZ,DT1,DT2)

        !CHUNK = CHUNK_J(Me%Size%JLB, Me%Size%JUB)
    
        !!$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)

j1:      do j = Me%WorkSize%JLB, Me%WorkSize%JUB
i1:      do i = Me%WorkSize%ILB, Me%WorkSize%IUB

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
        !!$OMP END PARALLEL

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


        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "VerticalAdvection")


    end subroutine VerticalAdvection

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
!        !!CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
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
!        !!!$OMP PARALLEL PRIVATE(I,J,K)
!        !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
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
!        !!!$OMP END DO
!        !!!$OMP END PARALLEL
!
!
!    end subroutine ModifyExplicitCoefs
    
    !---------------------------------------------------------------------------

    subroutine ModifyDrainageNetworkCoefs (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k, STAT_CALL, CHUNK
        real(8)                                     :: aux 
        integer                                     :: GWFlowLink
        !Begin-----------------------------------------------------------------
   
        call GetGWFlowOption (Me%ObjPorousMedia, GWFlowLink, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop ' ModifyDrainageNetworkCoefs - ModulePorousMediaProperties - ERR01'
        
        CurrProperty => PropertyX
        
       !Flux between river and runoff in layers
       
        if (GWFlowLink == Layer_) then

            CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,aux)
            
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if ((Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) .and. (Me%ExtVar%RiverPoints(I,J) == BasinPoint)) then   
                           
                    if((K .ge. Me%ExtVar%GWFlowBottomLayer(i,j)) .and. (K .le. Me%ExtVar%GWFlowTopLayer(i,j))) then
                        
                        !Auxuliar value for transport - units of flow^-1
                        !s/m3
                        aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                        
                        ! Positive flow -> looses mass
                        Me%COEFExpl%CoefInterfDN(i,j,k) = - aux * Me%ExtVar%FlowToChannelsLayer(i,j,k)

                      
                        ! mass going to channel -> conc from soil
                        if (Me%ExtVar%FlowToChannelsLayer(i,j,k) .gt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) =  CurrProperty%ConcentrationOld(i,j,k)
                            
                       
                        !mass coming from channel -> conc from DN
                        elseif (Me%ExtVar%FlowToChannelsLayer(i,j,k) .lt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) = CurrProperty%ConcentrationDN(i,j)
                            
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
            
            CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,aux)
            
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint .and. Me%ExtVar%RiverPoints(I,J) == BasinPoint) then   

                    if (K == Me%ExtVar%GWlayerOld(i,j)) then

                        !Auxuliar value for transport - units of flow^-1
                        !s/m3
                        aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))

                        ! Positive flow -> looses mass
                        Me%COEFExpl%CoefInterfDN(i,j,k) = - aux * Me%ExtVar%FlowToChannels(i,j)

                        
                        ! mass going to channel -> conc from soil
                        if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) =  CurrProperty%ConcentrationOld(i,j,k)

                        !mass coming from channel -> conc from DN
                        elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
                            
                            CurrProperty%ConcInInterfaceDN(i,j,k) = CurrProperty%ConcentrationDN(i,j)

                        endif
                    endif
                endif
            enddo
            enddo
            !$OMP END DO
            enddo  
            !$OMP END PARALLEL
       
        endif
   
   
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
        
       
        CurrProperty => PropertyX
        
        FluxW => Me%ExtVar%FluxW
       
        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
      
        !$OMP PARALLEL PRIVATE(i,j,k,aux,TopDiffusion)
        
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

   
   end subroutine ModifyInfiltrationCoefs

    !---------------------------------------------------------------------------

    subroutine ModifyVegetationCoefs (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        real, pointer, dimension(:,:,:)             :: TranspirationLayer
        integer                                     :: i, j, k, STAT_CALL, CHUNK
!        real                                        :: NitrogenUptake, PhosphorusUptake
!        real                                        :: CorrectedTranspirationFlux
        real                                        :: aux

        !Begin-----------------------------------------------------------------
        
        !all dissolved properties with advection diffusion will be transported with the transpiration
        !if in vegetation module is used swat uptake, then nitrate and phosphorus uptake has to be changed.
        !see this in the future
        
        CurrProperty => PropertyX
        
        call GetTranspiration(Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop ' ModifyVegetationCoefs - ModulePorousMediaProperties - ERR01'

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
      
        !$OMP PARALLEL PRIVATE(i,j,k,aux)
        
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)        
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then   
                
                if(K .ge. Me%ExtVar%TranspirationBottomLayer(i,j)) then
                    
                    !Auxuliar value for transport - units of flow-1
                    !s/m3
                    aux             = (Me%ExtVar%DT/(Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)))
                    
                    !Positive flow - looses mass
                    Me%COEFExpl%CoefInterfVeg(i,j,k) = - aux * TranspirationLayer(i,j,k)
                    
                    
!                    !Change nitrate and diss phosphorus flux because it may be disconnected from flow (if SWAT model original
!                    ! formulation is used)
!                    !in this case the transpiration flux * conc does not match the mass flux so here is changed
!                    !Concentration old was the concentration used in module vegetation to the uptake so it is maintained
!                    ! for consistency
!                    if (Me%ExtVar%ModelNitrogen .and. CurrProperty%ID%IDNumber == Nitrate_) then
!                        !      gN/s       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2 / s                     
!                        NitrogenUptake = Me%ExtVar%NitrogenUptake(i,j,k) * 1e3 * Me%ExtVar%Area(i,j) * 1e-4      &
!                                         / Me%ExtVar%VegetationDT
!                        
!                        !AssociatedFlux - g/s / g/m3 = g/s * m3/g = m3/s
!                        CorrectedTranspirationFlux        = NitrogenUptake * CurrProperty%ConcentrationOld(i,j,k)
!                        Me%COEFExpl%CoefInterfVeg(i,j,k) = - aux * CorrectedTranspirationFlux
! 
!                    elseif (Me%ExtVar%ModelPhosphorus .and. CurrProperty%ID%IDNumber == Inorganic_Phosphorus_) then
!                        !      gN/s       = KgN/ha * 1E3g/kg * (m2) * 1ha/10000m2 / s                     
!                        PhosphorusUptake = Me%ExtVar%PhosphorusUptake(i,j,k) * 1e3 * Me%ExtVar%Area(i,j) * 1e-4    &
!                                           / Me%ExtVar%VegetationDT
!                        
!                        !AssociatedFlux - g/s / g/m3 = g/s * m3/g = m3/s
!                        CorrectedTranspirationFlux        = PhosphorusUptake * CurrProperty%ConcentrationOld(i,j,k)
!                        Me%COEFExpl%CoefInterfVeg(i,j,k) = - aux * CorrectedTranspirationFlux
!                            
!                    endif

                endif
           endif
        
        enddo
        enddo
        !$OMP END DO
        enddo
        
        !$OMP END PARALLEL        
        
        call UngetPorousMedia (Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyVegetationCoefs - ModulePorousMediaProperties - ERR050'
                           
   
   end subroutine ModifyVegetationCoefs

    !---------------------------------------------------------------------------

    subroutine ModifyPropertyValues(PropertyX)
        
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k, CHUNK
        real(8)                                     :: coefB, CoefInterfDN, CoefInterfVeg
        real(8)                                     :: CoefInterfRunoff
        real(8), pointer, dimension(:,:,:)          :: FluxW
        !Begin-----------------------------------------------------------------    
        
        CurrProperty => PropertyX
        
        if (Me%AdvDiff_Explicit) then

            CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,CoefInterfRunoff,CoefInterfDN,CoefInterfVeg, coefB)
               
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)              
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then   
                    
                    CoefInterfRunoff = Me%COEFExpl%CoefInterfRunoff(i,j,k)
                    CoefInterfDN     = Me%COEFExpl%CoefInterfDN(i,j,k)
                    CoefInterfVeg    = Me%COEFExpl%CoefInterfVeg(i,j,k)
                   
                    CoefB = Me%ExtVar%WaterContentOld(i,j,k)/Me%ExtVar%WaterContent(i,j,k)
                    
                    Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k) + CoefB * CurrProperty%ConcentrationOld(i,j,k  )                  &
                                                      + CoefInterfRunoff * CurrProperty%ConcentrationOnInfColumn(i,j)         &
                                                      + CoefInterfDN     * CurrProperty%ConcInInterfaceDN(i,j,k)              &
                                                      + CoefInterfVeg    * CurrProperty%ConcentrationOld(i,j,k  )     
                                                                     
                    CurrProperty%Concentration(i,j,k) = Me%TICOEF3(i,j,k  )

                endif
           enddo
           enddo
           enddo
           !$OMP END DO
           !$OMP END PARALLEL  
           
        else !implicit
            
            FluxW => Me%ExtVar%FluxW

            CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
          
            !$OMP PARALLEL PRIVATE(i,j,k,CoefInterfRunoff,CoefInterfDN,CoefInterfVeg, coefB)
               
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)                          
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then           
                    
                    CoefInterfRunoff = Me%COEFExpl%CoefInterfRunoff(i,j,k)
                    CoefInterfDN     = Me%COEFExpl%CoefInterfDN(i,j,k)
                    CoefInterfVeg    = Me%COEFExpl%CoefInterfVeg(i,j,k)
                    
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
                    if (CoefInterfVeg .lt. 0.0) then
                        Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k) - CoefInterfVeg
                        CoefInterfVeg       = 0.0
                    endif

                    CoefB = Me%ExtVar%WaterContentOld(i,j,k)/Me%ExtVar%WaterContent(i,j,k)

                    Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k)                                                      &
                                        + coefB * CurrProperty%ConcentrationOld(i,j,k)                           &
                                        + CoefInterfRunoff * CurrProperty%ConcentrationOnInfColumn(i,j)          &
                                        + CoefInterfDN     * CurrProperty%ConcInInterfaceDN(i,j,k)               & 
                                        + CoefInterfVeg    * CurrProperty%ConcentrationOld(i,j,k  )

                endif
            enddo
            enddo
            enddo

           !$OMP END DO
           !$OMP END PARALLEL  
           
           !3D model or 1D vertical
            call THOMASZ(Me%WorkSize%ILB, Me%WorkSize%IUB,                              &
                         Me%WorkSize%JLB, Me%WorkSize%JUB,                              &
                         Me%WorkSize%KLB, Me%WorkSize%KUB,                              &
                         Me%COEF3%D,                                                    &
                         Me%COEF3%E,                                                    &
                         Me%COEF3%F,                                                    &
                         Me%TICOEF3,                                                    &
                         CurrProperty%Concentration,                                    &
                         Me%VECG,                                                       &
                         Me%VECW)      
                 
           
        endif
    
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
        real, pointer, dimension(:,:,:)             :: TranspirationLayer
        !Begin-----------------------------------------------------------------    
        
        CurrProperty => PropertyX
        
        !!Drainage network interface mass balance 
        if (Me%ExtVar%CoupledDN) then
            call GetGWFlowOption (Me%ObjPorousMedia, GWFlowLink, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop ' ModifyAdvectionDiffusion - ModulePorousMediaProperties - ERR01'
            
            CurrProperty => PropertyX
            
           !Flux between river and runoff in layers
           
            if (GWFlowLink == Layer_) then

                !!!$OMP PARALLEL PRIVATE(I,J,K)
                !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
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
                
                !!!$OMP PARALLEL PRIVATE(I,J,K)
                !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                    if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint .and. Me%ExtVar%RiverPoints(I,J) == BasinPoint) then   

                        if (K == Me%ExtVar%GWlayerOld(i,j)) then
                           
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
                enddo  
           
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
        
        !!Vegetation mass balance
        if (Me%ExtVar%CoupledVegetation .and. Me%ExtVar%ModelWater) then
    
            call GetTranspiration(Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop ' ModifyInterfaceMassFluxes - ModulePorousMediaProperties - ERR01'
    
    
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then   
                           
                    if(K .ge. Me%ExtVar%TranspirationBottomLayer(i,j) ) then
                        
                        !kg = [kg] + [m3/s] * [g/m3] * [1e-3kg/g] * [s] 
                        CurrProperty%MB%TranspiredMass = CurrProperty%MB%TranspiredMass +                                   &
                                                          (TranspirationLayer(i,j,k) * CurrProperty%Concentration(i,j,k)    &
                                                            *  1E-3 * Me%ExtVar%DT)
                        
    
                    endif
               endif
            
            enddo
            enddo
            enddo

            call UngetPorousMedia (Me%ObjPorousMedia, TranspirationLayer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyInterfaceMassFluxes - ModulePorousMediaProperties - ERR050'

        
        endif

    
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
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
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
    
    end subroutine ModifyDiffusivity_New

    !---------------------------------------------------------------------------

    subroutine SoilQualityProcesses

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property),          pointer     :: PropertyX
        type (T_Property),          pointer     :: SoilDryDensity, Salinity, pH
        type (T_Property),          pointer     :: IonicStrength, PhosphorusAdsortionIndex

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
            write(*,*) 'property phosphorus sdsortion index not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR40'
        endif

        call ComputeWindVelocity
 
        
        if (Me%ExtVar%Now .GE. Me%Coupled%SoilQuality_NextCompute) then
            
            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))
                

                call Modify_Interface(InterfaceID               = Me%ObjInterface,                         &
                                      PropertyID                = PropertyX%ID%IDNumber,                   &
                                      Concentration             = PropertyX%Concentration,                 &
                                      WaterPoints3D             = Me%ExtVar%WaterPoints3D,                 &
                                      OpenPoints3D              = Me%ExtVar%OpenPoints3D,                  &
                                      WaterPercentage           = Me%ExtVar%WaterContent,                  &
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
                

                PropertyX => PropertyX%Next
                

            end do
            
            Me%Coupled%SoilQuality_NextCompute = Me%Coupled%SoilQuality_NextCompute +       &
                                                     Me%Coupled%SoilQuality_DT

        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))


            if (PropertyX%Evolution%SoilQuality) then

                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%Concentration,          &
                                          WaterPoints3D = Me%ExtVar%WaterPoints3D,          &
                                          DTProp        = PropertyX%Evolution%DTInterval,   &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'SoilQuality_Processes - ModulePorousMediaProperties - ERR060'

                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SoilQualityProcesses")

    end subroutine SoilQualityProcesses
    
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
        type (T_Property), pointer              :: CohesiveSediment
        real                                    :: DT, MassTransfer
        real                                    :: DissolvedFraction
        real                                    :: ParticulateFraction
        real                                    :: TransferRate
!        real                                    :: ReferencePartitionCoef
!        real                                    :: PartitionCoef
!        real                                    :: b, SalinityConcentration
        real                                    :: RefSedFactor
        integer                                 :: STAT_CALL
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k
        integer                                 :: Couple_ID


        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "Partition_Processes")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KLB = Me%WorkSize%KLB 
        KUB = Me%WorkSize%KUB 

        PropertyX => Me%FirstProperty

do0:    do while(associated(PropertyX))
cd0:        if(.not. PropertyX%Particulate .and. PropertyX%Evolution%Partitioning) then

cd1:            if(Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    Couple_ID = PropertyX%Evolution%Partition%Couple_ID

                    call Search_Property(PartPropX, PropertyXID = Couple_ID, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'Partition_Processes - ModuleRunoffProperties - ERR02'

                    call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'Partition_Processes - ModuleRunoffProperties - ERR03'


    !                if (PropertyX%Evolution%Partition%SalinityEffect) then
    !
    !                    call Search_Property(Salinity, PropertyXID = Salinity_, STAT = STAT_CALL)
    !                    if (STAT_CALL /= SUCCESS_)                                      &
    !                        stop 'Partition_Processes - ModuleRunoffProperties - ERR03'
    !                endif

                
                    DT = PropertyX%Evolution%DTInterval
                    
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

                            if(PropertyX%Evolution%Partition%UseSedimentRefConc)then

                                RefSedFactor = CohesiveSediment%Concentration(i,j,k) / &
                                               PropertyX%Evolution%Partition%SedimentRefConc

                                if(RefSedFactor < 1.)then

                                    TransferRate = PropertyX%Evolution%Partition%Rate * RefSedFactor
                                
                                else

                                    TransferRate = PropertyX%Evolution%Partition%Rate

                                end if

                            else

                                TransferRate = PropertyX%Evolution%Partition%Rate

                            end if

                            ! [g/m3]       =          [s]         * [s^-1]        * [g/m3]
                            MassTransfer    =         DT * TransferRate *          &
                            (DissolvedFraction   * PartPropX%Concentration(i, j,k) -        &                  
                             ParticulateFraction * PropertyX%Concentration(i, j,k))

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

        
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "Partition_Processes")


    end subroutine Partition_Processes

    !--------------------------------------------------------------------------    
    
    subroutine SetLimitsConcentration

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k!, CHUNK
        
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
                
!                CHUNK = CHUNK_K(Me%Size%KLB, Me%Size%KUB)
                
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
    

    subroutine ComputeDissolvedToParticulate3D

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
         
        !Local----------------------------------------------------------------- 
        integer                                 :: i, j, k
        real                                    :: DT, InstantValue, ResidualValue
        type(T_Property), pointer               :: SoilDryDensity!, DrySedimentVolume
        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ComputeDissolvedToParticulate3D")


        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR01'

        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR02'


        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%ExtVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then

                ! [m3water/kgsed] = [m3water/m3cell]*[m3cell] / ([kgsed/m3cell] * [m3cell]) 
                InstantValue                       = Me%ExtVar%WaterContent(i,j,k) *  Me%ExtVar%CellVolume(i,j,k) / &
                                                    (SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k))

                ResidualValue                      = Me%DissolvedToParticulate3D(i,j,k)

                Me%DissolvedToParticulate3D(i,j,k) = (ResidualValue * Me%ResidualTime +     &
                                                      InstantValue * DT) / (Me%ResidualTime + DT)
                                                       
            end if
        end do
        end do
        end do

        Me%ResidualTime = Me%ResidualTime + DT
        
        
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
                Me%CellWaterMass(I, J, K) = WaterReferenceDensity * Me%CellWaterVolume(I, J, K)
                
                if (associated(Me%PropertySoilDryDensity)) then
                    !          kg            =                      kg/m3                       *         m3
                    Me%CellSoilMass(I, J, K) = Me%PropertySoilDryDensity%Concentration(I, J, K) * CellVolume(I, J, K)
                else
                    !          kg            = 1000 *         m3
                    Me%CellSoilMass(I, J, K) = 1000 * CellVolume(I, J, K)
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
                        call Modify_Interface(InterfaceID        = Me%ObjInterfaceSoilChemistry,  &
                                              PropertyID         = PropertyX%ID%IDNumber,         &
                                              WaterVolume        = Me%CellWaterVolume,            &
                                              WaterMass          = Me%CellWaterMass,              &
                                              SolidMass          = Me%CellSoilMass,               &
                                              Temperature        = TemperatureProp%Concentration, &
                                              pH                 = pHProp%Concentration,          &
                                              pE                 = pEProp%Concentration,          &
                                              Concentration      = PropertyX%Concentration,       &
                                              WaterPoints3D      = Me%ExtVar%WaterPoints3D,       &
                                              OpenPoints3D       = Me%ExtVar%OpenPoints3D,        &
                                              STAT               = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) &
                            stop 'SoilChemistryProcesses - ModulePorousMediaProperties - ERR001'
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

                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then
                
                    select case (PropertyX%ID%IDNumber)
                    case (Temperature_, pH_, pE_, SoilDryDensity_)
                    case default
                        call Modify_Interface(InterfaceID   = Me%ObjInterfaceSoilChemistry,     &
                                              PropertyID    = PropertyX%ID%IDNumber,            &
                                              Concentration = PropertyX%Concentration,          &
                                              WaterPoints3D = Me%ExtVar%WaterPoints3D,          &
                                              DTProp        = PropertyX%Evolution%DTInterval,   &
                                              STAT          = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                            &
                            stop 'SoilChemistryProcesses - ModulePorousMediaProperties - ERR002'
                    end select

                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SoilChemistryProcesses")
        
        !End-------------------------------------------------------------------
        
    end subroutine SoilChemistryProcesses   
    !-----------------------------------------------------------------------------    
#endif

    !-----------------------------------------------------------------------------    
    subroutine ChainReactionsProcesses
    
        !Local-------------------------------------------------------------------- 
        type (T_Property), pointer         :: PropertyX
        integer, dimension(:), pointer     :: CRPropertiesList
        integer                            :: PropertiesCount
        integer                            :: STAT
        integer                            :: Index
        real, pointer, dimension(:,:,:)    :: WaterVolume !L
        real, pointer, dimension(:,:,:)    :: SoilMass    !kg
        type (T_Property), pointer         :: SoilDensity !kg/m3
        integer                            :: I, J, K
        
        !Begin-----------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ChainReactionsProcesses")

        !Change this so the allocation/deallocation will be done out of this function if ChainReactions is used
        allocate (WaterVolume(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
        allocate (SoilMass(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
        
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                WaterVolume(I, J, K) = Me%ExtVar%WaterContent(I, J, K) * WaterReferenceDensity
            endif
        
        enddo
        enddo
        enddo
        
        !Change this code because in some situations will not exist "soil properties"
        call Search_Property(SoilDensity, SoilVolumetricDensity_, STAT)        
        
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                SoilMass(I,J,K) = Me%ExtVar%CellVolume(I,J,K) * SoilDensity%Concentration(I,J,K)
            endif
        
        enddo
        enddo
        enddo
        
        call GetCRPropertiesList(Me%ObjChainReactions, CRPropertiesList, PropertiesCount, STAT = STAT)
        if (STAT .NE. SUCCESS_) &
            stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR010'
    
        do Index = 1, PropertiesCount
    
            call Search_Property(PropertyX, CRPropertiesList(Index), STAT)    
            if (STAT .NE. SUCCESS_) &                            
                stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR020'
                    
            call SetCRPropertyConcentration (Me%ObjChainReactions,      &
                                             CRPropertiesList(Index),   &
                                             PropertyX%Concentration,   &
                                             STAT)
            if (STAT .NE. SUCCESS_) &                            
                stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR030'
                                    
        end do                
        
        call UnGetChainReactions(Me%ObjChainReactions, CRPropertiesList, STAT)
        if (STAT .NE. SUCCESS_) &
            stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR040'        
               
        if (STAT .NE. SUCCESS_) &
            stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR050'
        
        call ModifyChainReactions (Me%ObjChainReactions, &
                                   WaterVolume,          &
                                   SoilMass,             &
                                   Me%ExtVar%DT,         &                                   
                                   STAT)                                       
        if (STAT .NE. SUCCESS_) &
            stop 'ChainReactionsProcesses - ModulePororusMediaProperties - ERR060'
            
        deallocate (SoilMass)
        deallocate (WaterVolume)
                                    
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ChainReactionsProcesses")    
        !-------------------------------------------------------------------------
    
    end subroutine ChainReactionsProcesses
    !-----------------------------------------------------------------------------    
    

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

cd1:        if (Property%Evolution%Variable) then
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
        integer                                 :: STAT_CALL, i, j

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX
        real, dimension(:,:), pointer           :: SurfaceDiffusivity

        !----------------------------------------------------------------------

        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data3D = PropertyX%Concentration,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR01'

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data2D = PropertyX%ConcentrationOnInfColumn,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR02'
                
                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data3D = PropertyX%Diffusivity,     &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR03'
                
                allocate(SurfaceDiffusivity(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
                
                !Only For Debug
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    SurfaceDiffusivity(i, j) =  PropertyX%Diffusivity(i, j, Me%WorkSize%KUB +1)
                enddo            
                enddo         

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data2D = SurfaceDiffusivity,        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR04'
                    
                deallocate(SurfaceDiffusivity)
            endif
            PropertyX=>PropertyX%Next
        enddo

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

                    endif

                    PropertyX => PropertyX%Next

                enddo

                Me%OutPut%NextOutput = OutPutNumber + 1

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR080'
            
            endif  TOut
        endif  TNum

    end subroutine OutPut_HDF

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
        type (T_Property), pointer                  :: CurrProperty
        !Begin-----------------------------------------------------------------
        
        CurrProperty => Me%FirstProperty
        
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR10'        
        
        call GetWaterContent    (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR020'

        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ    = Me%ExtVar%CellVolume,                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR030'
        
        do while (associated(CurrProperty)) 

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
                        
            CurrProperty => CurrProperty%Next
        end do 
        
        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR040'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR050'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%CellVolume,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModulePorousMediaProperties - ERR060'


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
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - PorousmediaProperties - ERR05'
                endif

               
                if (Me%OutPut%HDF_ON) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillVegetation - PorousmediaProperties  - ERR08'
                endif
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR07'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR08'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR11'
                
                nUsers = DeassociateInstance (mPOROUSMEDIA_,  Me%ObjPorousMedia)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR12'

                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR13'

                nUsers = DeassociateInstance (mMAP_,  Me%ObjMap)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR14'

                if (Me%Coupled%ChainReactions) then
                    call KillChainReactions (Me%ObjChainReactions, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillPorousMedia - PorousmediaProperties - ERR140'
                endif
                
                call DeallocateVariables

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

        !Water Content---------------------------------------------------------
        deallocate (Me%ExtVar%WindVelocity3D   )         
        
        if (Me%Coupled%AdvectionDiffusion) then
            deallocate (Me%WaterVolume)
            deallocate (Me%FluxWCorr)
            deallocate(Me%TICOEF3)
            
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
                deallocate(Me%COEF3%D)
                deallocate(Me%COEF3%E)
                deallocate(Me%COEF3%F)
                deallocate(Me%VECG)
                deallocate(Me%VECW)
              
            endif
        endif
        
        deallocate (Me%CellWaterVolume)
        
#ifdef _PHREEQC_        
        deallocate (Me%CellSoilMass)    
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

        call GetUnsatW          (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'
        
        call GetThetaS          (Me%ObjPorousMedia, Me%ExtVar%ThetaS, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR65'

        call GetPotentialInfiltration (Me%ObjPorousMedia, Me%ExtVar%InfiltrationColumn, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - Module ModulePorousMediaProperties. ERR10.'


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

        call GetGridData  (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
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


        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
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
        
    endsubroutine ReadUnlockExternalVar


end module ModulePorousMediaProperties

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 








