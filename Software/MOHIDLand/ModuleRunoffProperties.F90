!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land
! MODULE        : RunoffProperties
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Fev 2010 David Brito
! REVISION      : 
! DESCRIPTION   : Module to handle Properties in Runoff
!
!------------------------------------------------------------------------------
!
!
!Units in runoff properties
!   Transported properties (soluble)             : g/m3 (or mg/l)  
!   Adsorbed properties (non soluble)            : ug/kgsoil   
!   Bottom layer Properties (bottom transition)  : kg/m2    
!
! ADVDIFF_EXPLICIT              : 0/1               [1]        !REMARK: Horizontal diffusion is always explicit
!                                                               (1 - horiz adv is explicit; 0 - horiz adv is implicit 
! NEW_FORMULATION               : 0/1                [0]        !if 1 then spatial methods will be the same for all properties
!      ADVDIFF_METHOD_H         : integer      [UpwindOrder1]   !Spatial methods for horizontal advection
!                                                               !UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4,
!                                                                CentralDif = 5, LeapFrog = 6    !
! <beginproperty>
!   PARTICULATE                 : 0/1               [0]         !Property physical state: 0 - Dissolved ; 1 - Particulate
!     EROSION                   : 0/1               [0]         !Compute erosion (source/sink term) - only read if PARTICULATE : 1
!     DEPOSITION                : 0/1               [0]         !Compute deposition (source/sink) - only read if PARTICULATE : 1
!       WS_TYPE                 : integer           [1]         !1 -constant;2 -concentration function - only read if DEPOSITION : 1
!         WS_VALUE              : real                          !Fall velocity value - only read if WS_TYPE : 1
!   ADVECTION_DIFFUSION         : 0/1               [0]         !Property advection - diffusion
!       ADVDIFF_METHOD_H        : integer      [UpwindOrder1]   !Spatial methods for horizontal advection
!                                                               !UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4,
!                                                                CentralDif = 5, LeapFrog = 6    
!       ADVDIFF_TVD_LIMIT_H     : integer        [Superbee]     !Horizontal advection non-linear stability conditions
!                                                                MinMod = 1, VanLeer = 2, Muscl = 3, Superbee = 4, PDM = 5
!       ADVDIFF_VOLUME_RELATION_MAX : real          5.          !The relation between adjacent volumes above which 
!                                                               !the advection is upwind
!   SOIL_CHEMISTRY              : 0/1               [0]         !Use PREEQC model to change property (source/sink model)
!   SOIL_QUALITY                : 0/1               [0]         !Use SedimentQuality model to change property (source/sink model)
!   PARTITION                   : 0/1               [0]         !Compute partition between dissolved-particulate phases
!       PARTITION_COUPLE        : char               +          !Name of the property (oposite phase) to compute partition
!       PARTITION_FRACTION      : real               -          !Percentage of mass of a property in a determined phase 
!       PARTITION_RATE          : real            [1 s-1]       !Kinetic rate of partition to reach equilibrium
!       USE_SED_REF_CONC        : 0/1               [0]         !Use cohesive sediment concentration as a reference
!           SED_REF_CONC        : real              [1]         !Reference cohesive sediment concentration to partition
! <endproperty>
!
!------------------------------------------------------------------------------

Module ModuleRunoffProperties

    use ModuleGlobalData
    use ModuleStopWatch
    use ModuleFunctions
    use ModuleTime
    use ModuleHDF5
    use ModuleEnterData
    use ModuleProfile,          only : StartProfile, WriteProfile, KillProfile
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData,    &
                                       KillGridData
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie,                   &
                                       GetNumberOfTimeSeries, GetTimeSerieLocation,      &
                                       GetTimeSerieName, TryIgnoreTimeSerie,             &
                                       CorrectsCellsTimeSerie, KillTimeSerie         
    use ModuleHorizontalGrid,   only : GetHorizontalGridSize, GetHorizontalGrid,         &
                                       GetGridCellArea, GetXYCellZ,                      &
                                       WriteHorizontalGrid, UnGetHorizontalGrid
    use ModuleBasinGeometry,    only : GetBasinPoints, GetRiverPoints,  UnGetBasin 
                                       
    use ModuleFillMatrix,         only : ConstructFillMatrix, GetDefaultValue,             &
                                         KillFillMatrix, ModifyFillMatrix,                 &
                                         GetIfMatrixRemainsConstant
    use ModuleGeometry
    use ModuleHorizontalMap,      only : GetComputeFaces2D, UngetHorizontalMap
    use ModuleRunoff,             only : GetOverLandFlow, UnGetRunoff, GetRunoffWaterColumn,  &
                                         GetFlowToChannels, GetRunoffCenterVelocity,          &
                                         GetRunoffWaterColumnOld, GetRunoffWaterColumn,       &
                                         GetManning, GetManningDelta, GetRunoffCenterVelocity, &
                                         GetRunoffWaterColumnAT, GetBoundaryImposed,          &
                                         GetBoundaryCells, GetBoundaryFlux, GetFlowDischarge, &
                                         GetRouteDFour, GetRouteDFourCells, GetRouteDFourNeighbours, &
                                         GetRouteDFourFlux
                                           
!    use ModuleInterface,          only : ConstructInterface, Modify_Interface
    use ModuleAdvectionDiffusion, only : StartAdvectionDiffusion, AdvectionDiffusion,      &
                                         GetAdvFlux, GetDifFlux, GetBoundaryConditionList, &
                                         UngetAdvectionDiffusion, KillAdvectionDiffusion

    use ModuleDischarges        ,only : Construct_Discharges, GetDischargesNumber,       &
                                        GetDischargesGridLocalization,                   &
                                        GetDischargeWaterFlow, GetDischargesIDName,      &
                                        TryIgnoreDischarge, GetDischargeSpatialEmission, &
                                        CorrectsCellsDischarges, Kill_Discharges,        &
                                        GetDischargeConcentration, GetDischargeFlowDistribuiton, &
                                        UnGetDischarges
    
    use ModuleBoxDif,               only: StartBoxDif, GetBoxes, GetNumberOfBoxes, UngetBoxDif,     &
                                          BoxDif, KillBoxDif        

   implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !RunoffProperties
    public  :: ConstructRunoffProperties
    private ::      AllocateInstance
    private ::      ReadFileNames
    private ::      ReadGlobalOptions
    private ::      Construct_PropertyList
    private ::      ConstructPartition
    private ::      ConstructHDF
    private ::      ConstructTimeSerie
!    private ::      CoupleSoilQuality
!#ifdef _PHREEQC_       
!    private ::      CoupleSoilChemistry
!#endif    

    !Selector
    public  :: GetRPMassBalance
    public  :: GetRPnProperties
    public  :: GetRPPropertiesIDByIdx 
    public  :: GetRPOptions   
    public  :: GetRPConcentration     !New Conc at the end (after link with DN)
    public  :: GetRPConcentrationAT   !Conc After Transtport (needed in ModuleDrainageNetwork)
    public  :: GetRPConcentrationOld  !Old Conc (before transport)
    public  :: GetRPDecayRate
    public  :: CheckRPProperty
    public  :: SetDNConcRP        !RP gets DN conc 
    public  :: SetBasinConcRP     !RP gets Basin WC conc (updated each time a module changes water column)
    public  :: SetBasinToRPSplash !RP gets through fall and canopy height to compute splash erosion
    public  :: SetVegetationRP    !RP gets flux from vegetation to fluffy - fertilization organic particulate material
                                  ! and pesticide to water column if exists

!    public  :: SetWindVelocity
    public  :: UnGetRunoffProperties
    
    !Modifier
    public  :: ModifyRunoffProperties
    private ::      ActualizePropertiesFromFile
    private ::      AdvectionDiffusionProcesses 
!    private ::      SoilQualityProcesses
!#ifdef _PHREEQC_    
!    private ::      SoilChemistryProcesses
!#endif        
    private ::      ModifyBottomFluxes
    private ::      Partition_Processes
    private ::      SetLimitsConcentration
    private ::      OutPut_TimeSeries
    private ::      Output_HDF
    
    !Destructor
    public  :: KillRunoffProperties                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjRunoffProperties 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetRunoffProperties2D_I
    private :: UnGetRunoffProperties2D_R8
    private :: UnGetRunoffProperties2D_R4
    interface  UnGetRunoffProperties
        module procedure UnGetRunoffProperties2D_I
        module procedure UnGetRunoffProperties2D_R8
        module procedure UnGetRunoffProperties2D_R4
    end interface  UnGetRunoffProperties


    !Parameters----------------------------------------------------------------

    real,    parameter                :: WaterReferenceDensity = 1000. ![kg/m3]
    
    integer, parameter                :: DirectionX = 1
    integer, parameter                :: DirectionY = 2

    character(LEN = StringLength), parameter :: prop_block_begin     = '<beginproperty>'
    character(LEN = StringLength), parameter :: prop_block_end       = '<endproperty>'
    integer, parameter                :: AdvDif_ModuleRP_  = 1
    integer, parameter                :: AdvDif_ModuleAD_   = 2
    integer, parameter                :: AdvDif_Upwind_     = 1
    integer, parameter                :: AdvDif_CentralDif_ = 2    
    integer, parameter                :: AdvDif_Diff_Jury_  = 1
    integer, parameter                :: AdvDif_Diff_Old_   = 2        
    
    !Boundary condition
    integer, parameter                      :: ImposedValue_  = 1
    integer, parameter                      :: NullGradient_  = 2

    
    !Types---------------------------------------------------------------------
    
    private :: T_RunoffProperties

    type T_RelatedID
        integer                       :: IDNumber = -1
        character(LEN = StringLength) :: name     = ''   
    end type T_RelatedID

    type T_ID
        integer                       :: IDNumber     = null_int
        character(LEN = StringLength) :: name         = null_str
        character(LEN = StringLength) :: description  = null_str
        character(LEN = StringLength) :: units        = null_str
    end type T_ID

    type T_Property_2D
        type(T_PropertyID)               :: ID
        real, pointer, dimension (:,:)   :: Field     => null()
        real                             :: Scalar    = null_real
        logical                          :: Constant  = .false.
    end type T_Property_2D


    type T_ExtVar
        !Map
        integer, pointer, dimension(:,:,:)      :: LandPoints3D   => null()
        integer, dimension(:,:), pointer        :: BasinPoints    => null()
        integer, dimension(:,:), pointer        :: RiverPoints    => null()
        real                                    :: RunoffpropDT   = null_real
        type(T_Time)                            :: Now
        type(T_Time)                            :: BeginTime
        type(T_Time)                            :: EndTime
   
        ! from Runoff
        real,    dimension(:,:), pointer           :: CenterVelV   => null()
        real,    dimension(:,:), pointer           :: CenterVelU   => null()
        real,    dimension(:,:), pointer           :: FlowToChannels => null()
        real(8), pointer, dimension(:,:)           :: WaterColumn    => null()       !Final WaterColumn
        real(8), pointer, dimension(:,:)           :: WaterColumnOld => null()       !Initial Water Column
        real(8), pointer, dimension(:,:)           :: WaterColumnAT  => null()       !Water Column After Transport
        real(8), pointer, dimension(:,:)           :: CellVolume     => null()
        real(8), pointer, dimension(:,:)           :: CellWaterMass  => null()
        real(8),    dimension(:,:), pointer        :: FluxU          => null()
        real(8),    dimension(:,:), pointer        :: FluxV          => null()
        real,    pointer, dimension(:,:  )         :: Area           => null()
        real                                       :: DT             = null_real
        real,   pointer, dimension(:,:  )          :: DZY            => null()
        real,   pointer, dimension(:,:  )          :: DZX            => null()
        real,   pointer, dimension(:,:  )          :: DXX            => null()
        real,   pointer, dimension(:,:  )          :: DYY            => null()
        real,   pointer, dimension(:,:  )          :: DUX            => null()
        real,   pointer, dimension(:,:  )          :: DVY            => null()
        real   , pointer, dimension(:,:  )          :: Topography    => null()
        
       
        real(8), pointer, dimension(:,:)            :: InfiltrationFlux => null()

        logical                                     :: CoupledDN        = .false.

        !from basin
!        real,    dimension(:,:  ), pointer          :: WindVelocity2D  !m/s
!        real,    dimension(:,:  ), pointer          :: WindVelocity  !km/day
        real(8),    dimension(:,:  ), pointer       :: ThroughFall    => null()
        real,       dimension(:,:  ), pointer       :: CanopyHeight   => null()
        real(8),    dimension(:,:  ), pointer       :: CanopyDrainage => null()
        
        !from vegetation to fluff layer
        !If vegetation organic particulated material in fertilization
        logical                                     :: CoupledVegetation         = .false.
        logical                                     :: VegParticFertilization    = .false.
        logical                                     :: Pesticide                 = .false.
        real,       dimension(:,:  ), pointer       :: FertilOrganicNParticFluff => null()
        real,       dimension(:,:  ), pointer       :: FertilOrganicPParticFluff => null()
        integer                                     :: VegetationDT              = null_int
        
        !From Runoff Boundary
        logical                                     :: BoundaryImposed     = .false.
        integer, pointer, dimension(:,:)            :: BoundaryCells       => null()
        real   , pointer, dimension(:,:)            :: BoundaryFlux        => null()

        !From Runoff DFour
        logical                                     :: RouteDFour          = .false.
        integer, pointer, dimension(:,:)            :: RouteDFourCells     => null()
        integer, pointer, dimension(:,:)            :: RouteDFourLowerI    => null()
        integer, pointer, dimension(:,:)            :: RouteDFourLowerJ    => null()
        real   , pointer, dimension(:,:)            :: RouteDFourFlux      => null()
        
     end type T_ExtVar

    type T_OutPut
        type (T_Time), pointer, dimension(:)    :: OutTime                       => null()
        type (T_Time), dimension(:), pointer    :: RestartOutTime                => null()
        integer                                 :: NextOutPut                    = null_int
        integer                                 :: Number                        = null_int
        logical                                 :: Yes                           = .false.
        logical                                 :: TimeSerie_ON                  = .false.
        logical                                 :: HDF_ON                        = .false.
        logical                                 :: Profile_ON                    = .false.
        logical                                 :: WriteRestartFile              = .false.        
        logical                                 :: RestartOverwrite              = .false.
        logical                                 :: Boxes_ON                      = .false.
        integer                                 :: NextRestartOutput             = 1        
    end type T_OutPut

    type T_AdvectionDiffusion   
        !--For AdvectionDiffusion module use
        integer                                :: BoundaryCondition            = null_int
        real                                   :: SchmidtNumberH               = null_real
        real                                   :: SchmidtCoefV                 = null_real
        real                                   :: SchmidtBackgroundV           = null_real
        real                                   :: DiffusionH_imp_exp           = null_real
        real                                   :: ImplicitH_direction          = null_real
        logical                                :: Nulldif                      = .false.
        logical                                :: NumericStability             = .false.
        real                                   :: VolumeRelMax                 = null_real
        integer                                :: AdvMethodH                   = null_int
        integer                                :: TVDLimitationH               = null_int
        integer                                :: AdvMethodV                   = null_int
        integer                                :: TVDLimitationV               = null_int
        logical                                :: Upwind2H                     = .false.
        logical                                :: Upwind2V                     = .false.
        real                                   :: Molecular_Diff_Coef          = null_real
    end type T_AdvectionDiffusion

    type       T_Partition
        character (LEN = StringLength)          :: Couple               = null_str
        integer                                 :: Couple_ID            = null_int
        real                                    :: Fraction             = null_real
        real                                    :: Rate                 = null_real
        real                                    :: EmpiricCoef          = null_real
        real                                    :: SedimentRefConc      = null_real
        logical                                 :: UseSedimentRefConc   = .false.
        logical                                 :: SalinityEffect       = .false.
        
        !Isothermic Adsortion
        logical                                 :: NonLinear            = .false.
        character(LEN = StringLength)           :: NonLinear_ks_Units   = null_str
        type(T_Property_2D)                     :: Nu            
        type(T_Property_2D)                     :: Be          
        type(T_Property_2D)                     :: ks
        type(T_Property_2D)                     :: PartitionRate
        type(T_Property_2D)                     :: Fraction2D 
        character (LEN = StringLength)          :: Partition_Couple     = null_str   
    end type T_Partition

    type T_Boundary
        integer                                 :: BoundaryCondition   = null_int
        real                                    :: DefaultBoundary     = null_real
    end type T_Boundary

    type T_Evolution
        logical                                 :: Partitioning         = .false.
        logical                                 :: BottomFluxes         = .false.
        logical                                 :: Variable             = .false.
        logical                                 :: Pesticide            = .false.
        real                                    :: DTInterval           = null_real
        logical                                 :: DTIntervalAssociated = .false.
        type(T_Time)                            :: LastCompute
        type(T_Time)                            :: NextCompute
!        logical                                 :: SoilQuality
!        logical                                 :: SoilChemistry
        logical                                 :: CationExchangeProcess  = .false.
        logical                                 :: ChemEquilibriumProcess = .false.
        logical                                 :: AdvectionDiffusion     = .false.
        logical                                 :: Erosion                = .false.
        logical                                 :: Deposition             = .false.
        logical                                 :: SplashErosion          = .false.
        logical                                 :: SoilWaterFluxes        = .false.
        logical                                 :: Macropores             = .false.
        logical                                 :: MinConcentration       = .false.
        logical                                 :: WarnOnNegativeValues   = .false.
        logical                                 :: Decay                  = .false.
        real                                    :: DecayRate              = null_real
        logical                                 :: Discharges             = .false.        
        type (T_AdvectionDiffusion)             :: AdvDiff
        type (T_Partition                    )  :: Partition
        type (T_Boundary)                       :: Boundary                
    end type T_Evolution

    type T_MassBalance
        real(8)                                 :: TotalStoredMass        = null_real
        real(8)                                 :: DNExchangeMass         = null_real
        real(8)                                 :: TotalDischargeMass     = null_real
    end type T_MassBalance
    
    type T_Files
        character(PathLength)                   :: InitialFile            = null_str
        character(PathLength)                   :: DataFile               = null_str
        character(PathLength)                   :: FinalFile              = null_str
        character(PathLength)                   :: TransientHDF           = null_str
        character(PathLength)                   :: DataSedimentQualityFile = null_str
        character(PathLength)                   :: BoxesFile = null_str
        integer                                 :: AsciiUnit              = null_int
    end type T_Files    

    type T_Property
        type (T_PropertyID)                     :: ID
        real, dimension(:,:), pointer           :: Concentration            => null()
        real, dimension(:,:), pointer           :: ConcentrationAT          => null()
!        real, dimension(:,:), pointer           :: Mass                     => null()
        real, dimension(:,:), pointer           :: ConcentrationOld         => null()
        real, dimension(:,:), pointer           :: BottomConcentration      => null()
        real                                    :: BottomMinConc            = null_real
!        real, dimension(:,:), pointer           :: BottomConcentrationOld   => null()
        real, dimension(:,:), pointer           :: ConcentrationDN          => null()
        real, dimension(:,:), pointer           :: TotalConcentration       => null()
        real, dimension(:,:), pointer           :: ConcInInterfaceDN        => null()
        real, dimension(:,:), pointer           :: ConcInBoundary           => null()
        real, dimension(:,:), pointer           :: PesticideFlux            => null()
        logical                                 :: FallVelocity             = .false.
        real, dimension(:,:), pointer           :: ErosionRate              => null()
        real, dimension(:,:), pointer           :: DepositionRate           => null()
        real, dimension(:,:), pointer           :: Ws                       => null()
        real, dimension(:,:), pointer           :: SplashRate               => null() 
        integer                                 :: Ws_Type                  = null_int
        real                                    :: Ws_Value                 = null_real
!        type (T_Property_2D)                    :: CHS
!        type(T_Property_2D)                     :: KL
!        type(T_Property_2D)                     :: KL1
!        type(T_Property_2D)                     :: ML
!        type(T_Property_2D)                     :: M
!        type(T_Property_2D)                     :: IScoefficient       
        real                                    :: CHS                      = null_real
        real                                    :: KL                       = null_real
        real                                    :: KL1                      = null_real
        real                                    :: ML                       = null_real
        real                                    :: M                        = null_real
        real                                    :: IScoefficient            = null_real
        
        real, pointer, dimension(:,:)           :: Mass_Created             => null()
        real,    pointer, dimension(:,:)        :: ViscosityU               => null()
        real,    pointer, dimension(:,:)        :: ViscosityV               => null()
        type (T_Property), pointer              :: Next                     => null()
        type (T_Property), pointer              :: Prev                     => null()
        logical                                 :: Particulate              = .false.
        logical                                 :: Pesticide                = .false.
        type (T_Evolution)                      :: Evolution
        type (T_MassBalance)                    :: MB
        real, pointer, dimension(:,:)           :: Diff_Turbulence_H        => null()
        real, pointer, dimension(:,:)           :: Viscosity                => null()
        real, pointer, dimension(:,:)           :: Diffusivity              => null()

        logical                                 :: Old              = .false.
        real                                    :: MinValue         = null_real
        logical                                 :: WarnOnNegativeValues = .false.
        logical                                 :: TimeSerie        = .false.
        logical                                 :: BoxTimeSerie     = .false.
        logical                                 :: BoxTimeSerie2D   = .false.
        logical                                 :: OutputHDF        = .false.
        
    end type T_Property

    type T_Coupled
!        logical                                 :: SoilQuality          = .false. !Sediment source/sink model (Sediment Quality)
!        real                                    :: SoilQuality_DT
!        type (T_Time)                           :: SoilQuality_NextCompute

!#ifdef _PHREEQC_        
!        logical                                 :: SoilChemistry        = .false.  !Chemical reactions model (PhreeqC)
!        real                                    :: SoilChemistry_DT
!        type (T_Time)                           :: SoilChemistry_NextCompute
!#endif        
        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: MinConcentration     = .false.
        logical                                 :: WarnOnNegativeValues = .false.
        logical                                 :: BottomFluxes         = .false.
        logical                                 :: ErosionFluxes        = .false.
        logical                                 :: DepositionFluxes     = .false.
        logical                                 :: SplashErosionFluxes  = .false.
        logical                                 :: Partition            = .false.       
        logical                                 :: Decay                = .false.
    end type T_Coupled

    type T_ComputeOptions

        logical                                     :: AdvDiff_Explicit         = .false.   ! 0 - Implicit; 1 explicit
        logical                                     :: Splash_CriticalHeight    = .false.
        integer                                     :: Splash_ErosiveRainMethod = null_int   ! 1-constant, 2 - real rain
        logical                                     :: DTIntervalAssociated     = .false.
        logical                                     :: Discharges               = .false.
            
    end type T_ComputeOptions

    !Implicit coef for thomas matrix
    type       T_DEF
        real   , pointer, dimension(: , : )  :: D    => null()
        real(8), pointer, dimension(: , : )  :: E    => null()
        real   , pointer, dimension(: , : )  :: F    => null()
    end type T_DEF
    
    !Explicit coefs
    type       T_A_B_C_Explicit
        real, pointer, dimension(: , : )  :: CoefInterfDN         => null()
        real, pointer, dimension(: , : )  :: CoefInterfBoundary   => null()
    end type T_A_B_C_Explicit

    type       T_FluxCoef
        real   , pointer, dimension(: , : )  :: C_flux  => null()    !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : )  :: D_flux  => null()  !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : )  :: E_flux  => null()  !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : )  :: F_flux  => null()  !Coeficient to calculate AdvFlux and DifFlux
    end type T_FluxCoef

   type  T_Fluxes
        real, pointer, dimension(:,:)         :: AdvFluxX   => null()
        real, pointer, dimension(:,:)         :: AdvFluxY   => null()
        real, pointer, dimension(:,:)         :: AdvFluxZ   => null()
        
        real, pointer, dimension(:,:)         :: DifFluxX   => null()
        real, pointer, dimension(:,:)         :: DifFluxY   => null()
        real, pointer, dimension(:,:)         :: DifFluxZ   => null()
        
        real, pointer, dimension(:,:)         :: MassFluxesX   => null()
        real, pointer, dimension(:,:)         :: MassFluxesY   => null()
        real, pointer, dimension(:,:)         :: MassFluxesZ   => null()
    end type T_Fluxes
    
    type T_RunoffProperties
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjBasinGeometry     = 0
        integer                                     :: ObjRunoff            = 0
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjGridData          = 0
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjtimeSerie         = 0
        integer                                     :: ObjHDF5              = 0
        integer                                     :: ObjProfile           = 0
        integer                                     :: ObjInterface         = 0
        integer                                     :: ObjDischarges        = 0
        integer                                     :: ObjBoxDif            = 0
!#ifdef _PHREEQC_        
!        integer                                     :: ObjPhreeqC                = 0
!        integer                                     :: ObjInterfaceSoilChemistry = 0 
!#endif        
        type (T_ExtVar)                             :: ExtVar
        logical                                     :: CheckGlobalMass   = .false.    
        type (T_Files)                              :: Files
        type (T_OutPut)                             :: OutPut
        type (T_Property), pointer                  :: FirstProperty    => null() 
        type (T_Property), pointer                  :: LastProperty     => null()   
        type (T_RunoffProperties), pointer          :: Next             => null() 
        type (T_Coupled)                            :: Coupled
        type (T_Time)                               :: LastOutputHDF5
        type (T_ComputeOptions)                     :: ComputeOptions

        type(T_DEF)                                 :: COEF3
        type(T_A_B_C_Explicit)                      :: COEFExpl 
        type(T_FluxCoef)                            :: COEF3_HorAdvXX           !Horizont advection coeficients
        type(T_FluxCoef)                            :: COEF3_HorAdvYY           !Horizont advection coeficients
        type(T_Fluxes)                              :: Fluxes
        real, pointer, dimension(: , :    )         :: TICOEF3         => null()      
        real(8), pointer, dimension(:)              :: VECG            => null()         !Auxiliar thomas arrays 
        real(8), pointer, dimension(:)              :: VECW            => null()         !Auxiliar thomas arrays  
        
        integer                                     :: di              = null_int         !auxiliar direction   
        integer                                     :: dj              = null_int

        logical                                     :: RunoffProperties = .false.
        integer                                     :: PropertiesNumber = 0
        integer                                     :: NumberPropForBoxes     = 0
        real   , pointer, dimension(:,:)            :: DissolvedToParticulate2D => null()
        real                                        :: ResidualTime     = null_real
        
        integer                                     :: InstanceID       = null_int
        type (T_Size2D)                             :: Size, WorkSize

        type(T_Property_2D)                         :: Disper_Trans
        type(T_Property_2D)                         :: Disper_Longi
        
        real, dimension(:,:), pointer               :: ShearStress      => null()
        real, dimension(:,:), pointer               :: RainKineticRate  => null()
        type(T_Property_2D)                         :: ErosionCriticalShear
        type(T_Property_2D)                         :: ErosionCoefficient
        type(T_Property_2D)                         :: DepositionCriticalShear 
!        type(T_Property_2D)                         :: StoneFraction
!        type(T_Property_2D)                         :: ClayFraction    
        type(T_Property_2D)                         :: Kclay
!        type(T_Property_2D)                         :: ErosiveRain
         
!        real, dimension(:,:), pointer               :: ShearStressY
               
        real(8), pointer, dimension(:,:)            :: WaterVolume     => null()
        real(8), pointer, dimension(:,:)            :: CellMass        => null()
        integer, pointer, dimension(:,:)            :: DummyOpenPoints => null()

        real(8), pointer, dimension(:,:)            :: WaterColumnBT  => null()       !Water Column Before Transport
        
        integer                                     :: nPropWithDischarge   = 0
       
        logical                                     :: NewFormulation       = .false.   !New formulation for advection 
        integer                                     :: AdvDiff_AdvMethodH   = null_int  !methods are general for all the properties
        
        real                                        :: HminChezy            = null_real    !for shear stress computation
        real                                        :: HcriticSplash        = null_real    !for splash erosion
        real                                        :: Splash_ErosiveRainValue = null_real !for splash erosion

    end type  T_RunoffProperties

    !Global Module Variables
    type (T_RunoffProperties), pointer              :: FirstObjRunoffProperties => null()
    type (T_RunoffProperties), pointer              :: Me                       => null()


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructRunoffProperties(ObjRunoffPropertiesID,                           &
                                              ComputeTimeID,                              &
                                              HorizontalGridID,                           &
                                              HorizontalMapID,                            &
                                              BasinGeometryID,                            &
                                              RunoffID,                                   &
                                              GridDataID,                                 &
                                              DischargesID,                               &
                                              CoupledDN,                                  &
                                              CheckGlobalMass,                            &
                                              CoupledVegetation,                          &
                                              VegParticFertilization,                     &
                                              Pesticide,                                  &
                                              STAT)
     
        !Arguments---------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: BasinGeometryID
        integer                                         :: RunoffID
        integer                                         :: GridDataID
        integer                                         :: DischargesID
        logical, optional                               :: CoupledDN 
        logical                                         :: CheckGlobalMass
        logical                                         :: CoupledVegetation
        logical                                         :: VegParticFertilization
        logical                                         :: Pesticide
        integer, optional, intent(OUT)                  :: STAT 
        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_,STAT_CALL
        !------------------------------------------------------------------------
                                    

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mRunoffProperties_)) then
            nullify (FirstObjRunoffProperties)
            call RegisterModule (mRunoffProperties_) 
        endif

        call Ready(ObjRunoffPropertiesID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            

            call AllocateInstance

            !Associate External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           ComputeTimeID   )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjRunoff         = AssociateInstance (mRUNOFF_,         RunoffID        )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID      )
            
            Me%CheckGlobalMass = CheckGlobalMass

            if (present(CoupledDN)) then
                Me%ExtVar%CoupledDN  = CoupledDN
            endif            
            
            Me%ExtVar%CoupledVegetation      = CoupledVegetation
            Me%ExtVar%VegParticFertilization = VegParticFertilization
            Me%ExtVar%Pesticide              = Pesticide

            call CheckBoundary
            
            call CheckRouteDFour
            
            call ReadFileNames


            !Constructs the DataFile
            call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunoffProperties - ModuleRunoffProperties - ERR01'                
           
            call ReadGlobalOptions

            call Construct_PropertyList

            !call VerifyOptions
                        
            if (Me%Coupled%Partition) then
                call ConstructPartition
            end if            
            
            if (Me%ComputeOptions%Discharges) then
                call ConstructDischarges(DischargesID)
            endif
            
            call ConstructData2D
            
            call AllocateVariables
            
       
            call ConstructHDF    
    
            call ConstructTimeSerie

            if (Me%Output%Boxes_ON) then
                call StartOutputBoxFluxes
            endif            
            
            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunoffProperties - ModuleRunoffProperties - ERR010'



            !Couple nutrient, carbon and oxygen sources and sinks model
!            if (Me%Coupled%SoilQuality) then
!                call CoupleSoilQuality
!            endif
            
!#ifdef _PHREEQC_            
!            !Couple soil chemical model
!            if (Me%Coupled%SoilChemistry) then
!                call CoupleSoilChemistry
!            endif
!#endif
            
            !Message to user
            call ConstructLog

            if (Me%CheckGlobalMass) then
                call CalculateTotalStoredMass
            endif
            
            !First Output
            if (Me%Output%HDF_ON) then
                call OutPut_HDF
            endif


            !Returns ID
            ObjRunoffPropertiesID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructRunoffProperties - ModuleRunoffProperties - ERR040' 

        end if cd0


        if (present(STAT)) STAT = STAT_


    end subroutine ConstructRunoffProperties
 
    !--------------------------------------------------------------------------

    subroutine ConstructDischarges(DischargesID)

        !Arguments--------------------------------------------------------------
        integer                                     :: DischargesID 
        !Local------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: nDischarges, iDis
        character(len=StringLength)                 :: DischargeName
        real                                        :: CoordinateX, CoordinateY
        logical                                     :: CoordinatesON, IgnoreOK
        integer                                     :: Id, Jd
         

        !ObjDischarges comes from ModueRunoff
!        call Construct_Discharges(Me%ObjDischarges,                              &
!                                  Me%ObjTime,                                    &
!                                  STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR01'  

        if (DischargesID == 0)  then                                                
            write(*,*)'You need to define DISCHARGES : 1 in the RUNOFF input file' 
            stop      'ModuleRunoffProperties - ConstructDischarges - ERR01'
        else            
            Me%ObjDischarges     = AssociateInstance (mDISCHARGES_,     DischargesID    )
        endif
        
        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR010'

        do iDis = 1, nDischarges

            call GetDischargesGridLocalization(Me%ObjDischarges,                &
                                               DischargeIDNumber = iDis,        &
                                               CoordinateX   = CoordinateX,     &
                                               CoordinateY   = CoordinateY,     & 
                                               CoordinatesON = CoordinatesON,   &
                                               STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR020' 
                    
            call GetDischargesIDName (Me%ObjDischarges, iDis, DischargeName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR030' 

            if (CoordinatesON) then
                
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordinateX, CoordinateY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR040' 

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreDischarge(Me%ObjDischarges, iDis, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR050' 

                    if (IgnoreOK) then
                        write(*,*) 'Discharge outside the domain - ',trim(DischargeName)
                        cycle
                    else
                        stop 'ModuleRunoffProperties - ConstructDischarges - ERR060' 
                    endif

                endif

                call CorrectsCellsDischarges(Me%ObjDischarges, iDis, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModuleRunoffProperties - ConstructDischarges - ERR070' 
                    
            endif
            
        enddo
    
    
    endsubroutine ConstructDischarges            
    
    !--------------------------------------------------------------------------

    subroutine VerifyOptions
!        !Arguments-------------------------------------------------------------
!                                                    
!        !Local-----------------------------------------------------------------
!        logical                                        :: RouteD4
!        integer                                        :: STAT_CALL
!        !Begin-----------------------------------------------------------------
!        
!        call GetRouteD4 (Me%ObjRunoff, RouteD4, STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleRunoffProperties - ERR01'
!        
!        !RouteD4 is not yet accounted in runoff properties computation 
!        if (RouteD4) then
!            write(*,*)
!            write(*,*) 'Using ROUTE_D4 in ModuleRunoff'
!            write(*,*) 'This process is not yet accounted in RunoffProperties'
!            write(*,*) 'Disconnect ROUTE_D4 from ModuleRunoff'
!            stop 'VerifyOptions - ModuleRunoffProperties - ERR010'
!        endif
            
    end subroutine VerifyOptions
    
    !--------------------------------------------------------------------------
    subroutine CheckBoundary

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                        :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        !check if boundary is imposed
        call GetBoundaryImposed (Me%ObjRunoff, Me%ExtVar%BoundaryImposed, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckBoundary - ModuleRunoffProperties - ERR01'
        
        if (Me%ExtVar%BoundaryImposed) then
            call GetBoundaryCells (Me%ObjRunoff, Me%ExtVar%BoundaryCells, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckBoundary - ModuleRunoffProperties - ERR010'
        endif
        
    end subroutine CheckBoundary

    !--------------------------------------------------------------------------

    subroutine CheckRouteDFour

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                        :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        !check if RouteDFour is used
        call GetRouteDFour (Me%ObjRunoff, Me%ExtVar%RouteDFour, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRouteDFour - ModuleRunoffProperties - ERR01'
        
        if (Me%ExtVar%RouteDFour) then
            call GetRouteDFourCells (Me%ObjRunoff, Me%ExtVar%RouteDFourCells, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckRouteDFour - ModuleRunoffProperties - ERR010'
            
            call GetRouteDFourNeighbours (Me%ObjRunoff, Me%ExtVar%RouteDFourLowerI,  &
                                          Me%ExtVar%RouteDFourLowerJ, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckRouteDFour - ModuleRunoffProperties - ERR020'
        endif
        
    end subroutine CheckRouteDFour

    !--------------------------------------------------------------------------

    subroutine ConstructLog


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        integer                                     :: STAT_CALL


        write(*, *)"-------------------- RUNOFF PROPERTIES -------------------"
        write(*, *)
        write(*, *)"Num of Properties : ", Me%PropertiesNumber
        write(*, *)

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))

            write(*, *)"Property            : ", trim(CurrentProperty%ID%Name)
            write(*, *)"---Adv. Diff.       : ", CurrentProperty%Evolution%AdvectionDiffusion
            write(*, *)"---Particulate      : ", CurrentProperty%Particulate
            write(*, *)"---Bottom Fluxes    : ", CurrentProperty%Evolution%BottomFluxes                 
            write(*, *)"---Hydro. Erosion   : ", CurrentProperty%Evolution%Erosion
            write(*, *)"---Splash Erosion   : ", CurrentProperty%Evolution%SplashErosion
            write(*, *)"---Settl. Deposition: ", CurrentProperty%Evolution%Deposition
            write(*, *)"---Infil. Deposition:  T"                  !by default
            write(*, *)"---Partitioning     : ", CurrentProperty%Evolution%Partitioning
            write(*, *)"---Discharges       : ", CurrentProperty%Evolution%Discharges
            write(*, *)

            CurrentProperty=>CurrentProperty%Next
        enddo

        if(Me%Coupled%BottomFluxes) then
            call SearchProperty(CurrentProperty, PropertyXIDNumber = TSS_, STAT = STAT_CALL)
            !give a warning, if TSS is not chosen as property
            !cohesive sediment concentration must be taken in this case instead!!!
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'Bottom Fluxes are activated, but TSS is not chosen as property'
                write(*,*) 'Cohesive sediment will be taken to calculate erosion rates!'
                write(*,*)
            endif 
        endif            

    end subroutine ConstructLog
    
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_RunoffProperties), pointer                         :: NewObjRunoffProperties
        type (T_RunoffProperties), pointer                         :: PreviousObjRunoffProp


        !Allocates new instance
        allocate (NewObjRunoffProperties)
        nullify  (NewObjRunoffProperties%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjRunoffProperties)) then
            FirstObjRunoffProperties         => NewObjRunoffProperties
            Me                    => NewObjRunoffProperties
        else
            PreviousObjRunoffProp      => FirstObjRunoffProperties
            Me                    => FirstObjRunoffProperties%Next
            do while (associated(Me))
                PreviousObjRunoffProp  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjRunoffProperties
            PreviousObjRunoffProp%Next => NewObjRunoffProperties
        endif

        Me%InstanceID = RegisterNewInstance (mRunoffProperties_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    subroutine ReadFileNames

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
!        integer                                     :: iflag

        !Reads the name of the data file from nomfich
        call ReadFileName ('RUNOFF_PROP_DATA', Me%Files%DataFile, "Runoff Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoffProperties - ERR010'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('RUNOFF_PROP_HDF', Me%Files%TransientHDF, "Runoff HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoffProperties - ERR020'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('RUNOFF_PROP_FIN', Me%Files%FinalFile, "Runoff Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoffProperties - ERR030'


    end subroutine ReadFileNames
    
    !--------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        !Begin-----------------------------------------------------------------

        !Geometry Size
        call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                    Size     = Me%Size,                              &
                                    WorkSize = Me%WorkSize,                          &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR01'


        call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR20'

        call GetComputeTimeLimits(Me%ObjTime,                      &
                                  EndTime   = Me%ExtVar%EndTime,   &
                                  BeginTime = Me%ExtVar%BeginTime, &
                                  STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR30'

        ! Sets the last output equal to zero 
        call SetDate(Me%LastOutPutHDF5, 0, 0, 0, 0, 0, 0)

!        call GetData(Me%ComputeOptions%AdvDiff_DiffMethod,                                               &   
!                     Me%ObjEnterData, iflag,                                              &
!                     SearchType = FromFile,                                               &
!                     keyword    = 'ADVDIFF_DIFF_METHOD',                                  &
!                     Default    = AdvDif_Diff_Jury_,                                      &
!                     ClientModule ='ModuleRunoffProperties',                              &
!                     STAT       = STAT_CALL)            
!        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProeprties - ERR085'            
    

        call GetData(Me%ComputeOptions%AdvDiff_Explicit,           &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'ADVDIFF_EXPLICIT',            &
                     Default      =  .true.,                       & 
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)  &
            stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR90'

        !Want new formulation with advection methods not depending on property
        call GetData(Me%NewFormulation,                            &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'NEW_FORMULATION',             &
                     Default      = .false.,                       &
                     ClientModule = 'ModuleRunoffProperties',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)  &
            stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR100'

        if (Me%NewFormulation) then
            !horizontal advection method that will be general for all props
            call GetData(Me%AdvDiff_AdvMethodH,                        &
                         Me%ObjEnterData, iflag,                       &
                         SearchType   = FromFile,                      &
                         keyword      = 'ADVDIFF_METHOD_H',            &
                         Default      = UpwindOrder1,                  &
                         ClientModule = 'ModuleRunoffProperties', &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)  &
                stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR110'
            
        endif    

        !Min water column for chezy computation - used if erosion active
        call GetData(Me%HminChezy,                                                  &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'HMIN_CHEZY',                                   &
                     Default      =  AlmostZero,                                    & 
                     ClientModule = 'ModuleRunoffProperties',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR120' 
        if (Me%HminChezy .lt. 0.0) then
            write(*,*)'Minimum water column height for chezy computation HMIN_CHEZY can not be negative'
            stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR115'
        endif

!        !Splash rain model used 
!        call GetData(Me%ComputeOptions%Splash_Model,                                &
!                     Me%ObjEnterData, iflag,                                        &
!                     SearchType   = FromFile,                                       &
!                     keyword      = 'SPLASH_MODEL',                                 &
!                     Default      =  Default,                                       & 
!                     ClientModule = 'ModuleRunoffProperties',                       &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR1150' 
        
        !Critical water column for splash - used if splash erosion active
        call GetData(Me%ComputeOptions%Splash_CriticalHeight,                       &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'SPLASH_HCRITIC_COMPUTE',                       &
                     Default      =  .false.,                                       & 
                     ClientModule = 'ModuleRunoffProperties',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR120' 
        
        if (Me%ComputeOptions%Splash_CriticalHeight) then
            !Critical water column for splash - inflexion point of exponential curve - used if splash erosion active
            call GetData(Me%HcriticSplash,                                              &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'SPLASH_HCRITIC',                               &
                         Default      =  0.1,                                           & 
                         ClientModule = 'ModuleRunoffProperties',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR120' 
            if (Me%HcriticSplash .lt. 0.0) then
                write(*,*)'Critical water column height for splash computation SPLASH_HCRITIC can not be negative'
                stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR125'
            endif
        endif
        
        !in splash erosion use erosive rain (case 1 - constant) or real rain (case 2) - used if splash erosion active
        call GetData(Me%ComputeOptions%Splash_ErosiveRainMethod,                    &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'SPLASH_EROSIVERAIN_METHOD',                    &
                     Default      =  2,                                             & 
                     ClientModule = 'ModuleRunoffProperties',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR130' 
        
        !in formulation of splash erosion it is advised 30mm/h as erosiverain for mediterranean region
        if (Me%ComputeOptions%Splash_ErosiveRainMethod == 1) then
            call GetData(Me%Splash_ErosiveRainValue,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'SPLASH_EROSIVERAIN_VALUE',                     &
                         Default      =  30.,                                           & 
                         ClientModule = 'ModuleRunoffProperties',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ReadGlobalOptions - ERR140' 
        endif


    end subroutine ReadGlobalOptions        

    !--------------------------------------------------------------------------

    subroutine AllocateVariables        
        
        !Local-----------------------------------------------------------------        
        integer                                         :: ILB, IUB, JLB,  JUB, IJLB, IJUB 
        integer                                         :: STAT_CALL
        type (T_Property), pointer                      :: NewProperty, ParticulateProperty 
        
        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        
        if (Me%Coupled%AdvectionDiffusion) then
            allocate(Me%WaterVolume(ILB:IUB, JLB:JUB))
            Me%WaterVolume = 0.0
            
            allocate(Me%DummyOpenPoints(ILB:IUB, JLB:JUB))
            Me%DummyOpenPoints = 0
            
            allocate(Me%WaterColumnBT(ILB:IUB, JLB:JUB))
            Me%WaterColumnBT = 0.0
            
            allocate (Me%COEFExpl%CoefInterfDN     (ILB:IUB,JLB:JUB))
            Me%COEFExpl%CoefInterfDN        = 0.0
            allocate (Me%COEFExpl%CoefInterfBoundary (ILB:IUB,JLB:JUB))
            Me%COEFExpl%CoefInterfBoundary  = 0.0

            allocate(Me%TICOEF3                 (ILB:IUB, JLB:JUB))
            Me%TICOEF3                   = 0.0

           
            allocate(Me%COEF3_HorAdvXX%C_Flux    (ILB:IUB, JLB:JUB))
            allocate(Me%COEF3_HorAdvXX%D_Flux    (ILB:IUB, JLB:JUB))
            allocate(Me%COEF3_HorAdvXX%E_Flux    (ILB:IUB, JLB:JUB))
            allocate(Me%COEF3_HorAdvXX%F_Flux    (ILB:IUB, JLB:JUB))
            
            allocate(Me%COEF3_HorAdvYY%C_Flux    (ILB:IUB, JLB:JUB))
            allocate(Me%COEF3_HorAdvYY%D_Flux    (ILB:IUB, JLB:JUB))
            allocate(Me%COEF3_HorAdvYY%E_Flux    (ILB:IUB, JLB:JUB))
            allocate(Me%COEF3_HorAdvYY%F_Flux    (ILB:IUB, JLB:JUB)) 
                
            Me%COEF3_HorAdvXX%C_Flux     = Null_real
            Me%COEF3_HorAdvXX%D_Flux     = Null_real
            Me%COEF3_HorAdvXX%E_Flux     = Null_real
            Me%COEF3_HorAdvXX%F_Flux     = Null_real 

            Me%COEF3_HorAdvYY%C_Flux     = Null_real
            Me%COEF3_HorAdvYY%D_Flux     = Null_real
            Me%COEF3_HorAdvYY%E_Flux     = Null_real
            Me%COEF3_HorAdvYY%F_Flux     = Null_real                                
            
         
            if (.not. Me%ComputeOptions%AdvDiff_Explicit) then
                
                allocate(Me%COEF3%D                 (ILB:IUB, JLB:JUB))
                allocate(Me%COEF3%E                 (ILB:IUB, JLB:JUB))
                allocate(Me%COEF3%F                 (ILB:IUB, JLB:JUB))
                
                IJLB = min (ILB, JLB)
                IJUB = max (IUB, JUB)
                allocate(Me%VECG                    (IJLB:IJUB))
                allocate(Me%VECW                    (IJLB:IJUB))

                Me%COEF3%D                  = 0.0
                Me%COEF3%E                  = 1.0
                Me%COEF3%F                  = 0.0
                Me%VECG                     = Null_real 
                Me%VECW                     = Null_real 
           
            endif

            if (Me%Output%Boxes_ON) then
                allocate(Me%Fluxes%AdvFluxX    (ILB:IUB, JLB:JUB))
                allocate(Me%Fluxes%AdvFluxY    (ILB:IUB, JLB:JUB))
                allocate(Me%Fluxes%AdvFluxZ    (ILB:IUB, JLB:JUB))
                
                allocate(Me%Fluxes%DifFluxX    (ILB:IUB, JLB:JUB))
                allocate(Me%Fluxes%DifFluxY    (ILB:IUB, JLB:JUB))
                allocate(Me%Fluxes%DifFluxZ    (ILB:IUB, JLB:JUB))
                
                allocate(Me%Fluxes%MassFluxesX (ILB:IUB, JLB:JUB))
                allocate(Me%Fluxes%MassFluxesY (ILB:IUB, JLB:JUB))
                allocate(Me%Fluxes%MassFluxesZ (ILB:IUB, JLB:JUB))
            
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
            
        endif
        
        if(Me%Coupled%BottomFluxes) then
            !This test must be here and not inside Construct property because it looks for the desired properties 
!            call SearchProperty(NewProperty, PropertyXIDNumber = TSS_, STAT = STAT_CALL)
!            !give a warning, if TSS is not chosen as property
!            !cohesive sediment concentration must be taken in this case instead!!!
!            if (STAT_CALL /= SUCCESS_) then
!                write(*,*) 'Bottom Fluxes are activated, but TSS is not chosen as property'
!                write(*,*) 'Cohesive sediment will be taken to calculate erosion rates!'
!            endif             
!            if (Me%ComputeOptions%CalcFractionSediment)then
!                call SearchProperty(PropertyX, PropertyXIDNumber = COHSED_FINE_, STAT = STAT_CALL)
!                if (STAT_CALL == SUCCESS_) then
!                    allocate (Me%ShearStress (Me%TotalReaches))            
!                    Me%ShearStress = 0.0
!                else
!                    call SearchProperty(PropertyX, PropertyXIDNumber = COHSED_MEDIUM_, STAT = STAT_CALL)
!                    if (STAT_CALL == SUCCESS_) then
!                        allocate (Me%ShearStress (Me%TotalReaches))            
!                        Me%ShearStress = 0.0
!                    else
!                        call SearchProperty(PropertyX, PropertyXIDNumber = COHSED_COARSE_, STAT = STAT_CALL)
!                        if (STAT_CALL == SUCCESS_) then
!                            allocate (Me%ShearStress (Me%TotalReaches))            
!                            Me%ShearStress = 0.0
!                        else
!                            write (*,*)
!                            write (*,*) 'Bottom Fluxes needs at least one Cohesive Sediment Fraction'
!                            stop        'ConstructSubModules - ModuleDrainageNetwork - ERR02_Wassim'
!                        end if
!                    end if
!                end if
!            else
            call SearchProperty(NewProperty, PropertyXIDNumber = Cohesive_Sediment_, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                allocate (Me%ShearStress(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                Me%ShearStress = 0.0
    !                allocate (Me%ShearStressY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
    !                Me%ShearStressY = 0.0
                
            else
                write (*,*)
                write (*,*) 'Bottom Fluxes needs Cohesive_Sediment_'
                stop        'AllocateVariables - ModuleRunoffProperties - ERR0100'
            end if
    !            end if  

            if (Me%Coupled%SplashErosionFluxes) then
                
                allocate(Me%RainKineticRate(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR130'
                Me%RainKineticRate = 0.0

                allocate(Me%ExtVar%ThroughFall(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR150'
                Me%ExtVar%ThroughFall = FillValueReal

                allocate(Me%ExtVar%CanopyHeight(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR160'
                Me%ExtVar%CanopyHeight = FillValueReal

                allocate(Me%ExtVar%CanopyDrainage(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR170'
                Me%ExtVar%CanopyDrainage = FillValueReal
               
            endif
            
        endif
        
        if (Me%ExtVar%CoupledVegetation .and. Me%ExtVar%VegParticFertilization) then

            !Check if Property PON exists
            call Search_Property(ParticulateProperty, PropertyXID = PON_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'Need property particulate organic nitrogen'
                write(*,*) 'in module runoff properties because vegetation is fertilizing'
                write(*,*) 'particulate organic fertilizer that goes to fluff layer'
                stop 'AllocateVariables - ModuleRunoffProperties - ERR0180'
            endif
            
            !Check if Property POP exists
            call Search_Property(ParticulateProperty, PropertyXID = POP_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) 'Need property particulate organic phosphorus'
                write(*,*) 'in module runoff properties because vegetation is fertilizing'
                write(*,*) 'particulate organic fertilizer that goes to fluff layer'
                stop 'AllocateVariables - ModuleRunoffProperties - ERR0190'
            endif
            
            !flux in kg/ha
            allocate(Me%ExtVar%FertilOrganicNParticFluff(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleRunoffProperties - ERR200'
            Me%ExtVar%FertilOrganicNParticFluff = 0.

            allocate(Me%ExtVar%FertilOrganicPParticFluff(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleRunoffProperties - ERR210'
            Me%ExtVar%FertilOrganicPParticFluff = 0.
            
        endif

    end subroutine AllocateVariables


    !--------------------------------------------------------------------------
    
!    !--------------------------------------------------------------------------
!  
!    subroutine ConstructScalar2D(Scalar2D, ExtractType, ClientNumber, block_begin, block_end)
!
!        !Arguments-------------------------------------------------------------
!        type(T_Property_2D), pointer        :: Scalar2D
!        integer, intent(in)                 :: ExtractType
!        integer, intent(in), optional       :: ClientNumber
!        character(len=*)                    :: block_begin, block_end
!
!        !External--------------------------------------------------------------
!        integer                             :: STAT_CALL
!        logical                             :: BlockFound
!        integer                             :: BlockClientNumber
!
!        !Local-----------------------------------------------------------------
!        integer                             :: ILB, IUB, JLB, JUB
!        
!        !----------------------------------------------------------------------
!        ILB = Me%Size%ILB
!        IUB = Me%Size%IUB
!        JLB = Me%Size%JLB
!        JUB = Me%Size%JUB
!
!        select case(ExtractType)
!
!            case(FromBlock)
!
!                call ExtractBlockFromBuffer(Me%ObjEnterData, BlockClientNumber, block_begin, block_end,  &
!                                            BlockFound, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleRunoffProperties - ERR01'
!
!
!            case(FromBlockInBlock)
!
!                if(.not. present(ClientNumber))then
!                    stop 'ConstructScalar2D - ModuleSoilProperties - ERR02'
!                end if
!                
!                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber, block_begin, block_end,  &
!                                           BlockFound, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR20'
!
!        end select
!
!        if(BlockFound)then
!
!            allocate(Scalar2D%Field(ILB:IUB, JLB:JUB))
!
!            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR030'
!
!            call ConstructFillMatrix  (PropertyID           = Scalar2D%ID,                      &
!                                       EnterDataID          = Me%ObjEnterData,                  &
!                                       TimeID               = Me%ObjTime,                       &
!                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
!                                       ExtractType          = ExtractType,                      &
!                                       PointsToFill2D       = Me%ExtVar%BasinPoints,            &
!                                       Matrix2D             = Scalar2D%Field,                   &
!                                       TypeZUV              = TypeZ_,                           &
!                                       STAT                 = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR040'
!
!
!            call GetDefaultValue(Scalar2D%ID%ObjFillMatrix, Scalar2D%Scalar, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR050'
!
!            call UnGetBasin(Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL) 
!            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR60'
!
!            call KillFillMatrix(Scalar2D%ID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR065'
!
!            if(ExtractType == FromBlockInBlock)then
!                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR070'
!            end if
!
!
!            if(ExtractType == FromBlock)then
!                call Block_Unlock(Me%ObjEnterData, BlockClientNumber, STAT = STAT_CALL) 
!                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR080'
!
!                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR090'
!            end if
!        
!        else
!            write(*,*) 'Block not present:', block_begin, block_end
!            stop 'ConstructScalar2D - ModuleRunoffProperties - ERR100'
!        end if
!
!   
!    end subroutine ConstructScalar2D
!
!    !--------------------------------------------------------------------------
  
    subroutine Construct_PropertyList

        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_Property), pointer          :: NewProperty

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
                        stop 'Construct_PropertyList - ModuleRunoffProperties - ERR01'
                    exit do1    !No more blocks
                
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_PropertyList - ModuleRunoffProperties - ERR02'
            
            else cd1
                
                stop 'Construct_PropertyList - ModuleRunoffProperties - ERR03'
            
            end if cd1
        
        end do do1

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------------    
    
    subroutine Construct_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer           :: NewProperty

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_)stop 'Construct_Property - ModuleRunoffProperties - ERR00'
        
        nullify(NewProperty%Prev, NewProperty%Next)
        nullify(NewProperty%Concentration)
!        nullify(NewProperty%MassOnWaterColumn)

        call ConstructPropertyID            (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call Construct_PropertyState        (NewProperty)
        
        call Construct_PropertyEvolution    (NewProperty)

        call Construct_PropertyValues       (NewProperty)

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
                     ClientModule = 'ModuleRunoffProperties',                       &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModuleRunoffProperties - ERR01'
        
        if (NewProperty%Particulate)then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'Construct_PropertyState - ModuleRunoffProperties - ERR03'
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
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR10'

        if (NewProperty%Evolution%AdvectionDiffusion) then
            Me%Coupled%AdvectionDiffusion = .true.
            NewProperty%Evolution%Variable = .true.
        endif
        
        if (NewProperty%Evolution%AdvectionDiffusion) then

            call ReadAdvectionDiffusionParam (NewProperty)
        
        end if

        call ConstructPropertyDiffusivity (NewProperty)

        !By default the boundary value is imposed - currently water can only exit but
        !useful in the future
        call GetData(NewProperty%Evolution%Boundary%BoundaryCondition,                   &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'BOUNDARY_CONDITION',                                &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = ImposedValue_,                                       &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR020'

        if (NewProperty%Evolution%Boundary%BoundaryCondition /= ImposedValue_   .and.    &
            NewProperty%Evolution%Boundary%BoundaryCondition /= NullGradient_        ) then 
            write(*,*) ' Boundary Condition can only be ImposedValue = 1 or'
            write(*,*) ' NullGradient = 2. Check BOUNDARY_CONDITION keyword'
            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR25'
        endif
        
        !By default if not given the property enters with zero concentration.
        call GetData(NewProperty%Evolution%Boundary%DefaultBoundary,                     &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DEFAULTBOUNDARY',                                   &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = 0.0,                                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR030'
        

        !Checks for Bottom Fluxes
        call GetData(NewProperty%Evolution%BottomFluxes,                            &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'BOTTOM_FLUXES',                                &
                     Default      = .false.,                                        & 
                     ClientModule = 'ModuleRunoffProperties',                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR31' 

        if (NewProperty%Evolution%BottomFluxes) then
        
            Me%Coupled%BottomFluxes = .true.
            
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR40' 
            end if            
        endif

       !in Runoff properties all properties recognized by the model as particulate need to
       !have BOTTOM FLUXES ON because if all water infiltrates the mass needs to go somewhere
       !and so needs the bottom concentration (process controlled in Basin)
        if(Check_Particulate_Property(NewProperty%ID%IDNumber) .and.  &
           .not. NewProperty%Evolution%BottomFluxes) then 
            write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' has not BOTTOM_FLUXES ON'
            write(*,*) 'but is recognised by the model as particulate.'
            write(*,*) 'Particulated recognized properties can accumulate in bottom and'
            write(*,*) 'need BOTTOM_FLUXES to be active for the propery in Runoff Properties'
            stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR50'
        end if 
        
        if (NewProperty%Evolution%BottomFluxes) then
            
            !Compute SplashErosion
            call GetData(NewProperty%Evolution%SplashErosion,                           &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'SPLASH_EROSION',                               &
                         Default      =  .false.,                                       & 
                         ClientModule = 'ModuleRunoffProperties',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR35.6' 
            
            !Compute erosion fluxes
            call GetData(NewProperty%Evolution%Erosion,                                 &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'EROSION',                                      &
                         Default      =  .false.,                                       & 
                         ClientModule = 'ModuleRunoffProperties',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR36' 
        

            !Compute deposition fluxes
            call GetData(NewProperty%Evolution%Deposition,                              &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'DEPOSITION',                                   &
                         Default      =  .false.,                                        & 
                         ClientModule = 'ModuleRunoffProperties',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR39' 
            
        endif


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

!        call GetData(NewProperty%Evolution%SoilQuality,                                  &
!                     Me%ObjEnterData,iflag,                                              &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'SOIL_QUALITY',                                      &
!                     ClientModule = 'ModuleRunoffProperties',                            &
!                     default      = OFF,                                                 &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR20'
!
!        if (NewProperty%Evolution%SoilQuality) then
!            Me%Coupled%SoilQuality     = .true.
!            NewProperty%Evolution%Variable = .true.
!        endif


!#ifdef _PHREEQC_
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

!        call GetData(NewProperty%Evolution%SoilChemistry,                                &
!                     Me%ObjEnterData,iflag,                                              &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'SOIL_CHEMISTRY',                                    &
!                     ClientModule = 'ModuleRunoffProperties',                            &
!                     default      = OFF,                                                 &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR30'
!
!        if (NewProperty%Evolution%SoilChemistry) then
!            Me%Coupled%SoilChemistry       = .true.
!            NewProperty%Evolution%Variable = .true.
!        endif
!#endif

!        call GetData(NewProperty%Evolution%Pesticide,                                   &
!                     Me%ObjEnterData, iflag,                                             &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'PESTICIDE',                                         &
!                     Default      = .false.,                                             &
!                     ClientModule = 'ModuleRunoffProperties',                            &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Subroutine Construct_PropertyEvolution - ModuleRunoffProperties - ERR30'

        call GetData(NewProperty%Evolution%Decay,                                        &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DECAY',                                             &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR60'

        if (NewProperty%Evolution%Decay) then
            Me%Coupled%Decay       = .true.
        endif


        if (NewProperty%Evolution%Decay) then
            
            !Decay rate k (day-1) in P = Po*exp(-kt)
            call GetData(NewProperty%Evolution%DecayRate,                                    &
                         Me%ObjEnterData,iflag,                                              &
                         SearchType   = FromBlock,                                           &
                         keyword      = 'DECAY_RATE',                                        &
                         ClientModule = 'ModuleRunoffProperties',                            &
                         default      = 0.,                                                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                     &
                stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR70'            
            
        endif

        call GetData(NewProperty%Evolution%Discharges,                              &
                     Me%ObjEnterData, iflag,                                        &
                     Keyword        = 'DISCHARGES',                                 &
                     ClientModule   = 'ModuleRunoffProperties',                      &
                     SearchType     = FromBlock,                                    &
                     Default        = OFF,                                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR80' 
        
        if (NewProperty%Evolution%Discharges) then
            Me%ComputeOptions%Discharges = .true.
            Me%nPropWithDischarge = Me%nPropWithDischarge + 1
        endif
        
        !Property time step
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR40'

            ModelDT = Me%ExtVar%DT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DTINTERVAL',                                    &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModuleRunoffProperties',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR050'

            if (iflag == 1) then
                NewProperty%Evolution%DTIntervalAssociated = .true.
                Me%ComputeOptions%DTIntervalAssociated     = .true.                           
                                       
                if (NewProperty%Evolution%DTInterval < ModelDT) then
                    write(*,*) 
                    write(*,*) 'Property time step is smaller then model time step'
                    stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR60'

                elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                    !Property time step must be a multiple of the model time step
                    auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                    Erroraux = auxFactor - int(auxFactor)
                    if (Erroraux /= 0) then
                        write(*,*) 
                        write(*,*) 'Property time step must be a multiple of model time step.'
                        write(*,*) 'Please review your input data.'
                        stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR70'
                    endif

                    !Run period in seconds
                    DTaux = Me%ExtVar%EndTime - Me%ExtVar%Now

                    !The run period   must be a multiple of the Property DT
                    auxFactor = DTaux / NewProperty%Evolution%DTInterval

                    ErrorAux = auxFactor - int(auxFactor)
                    if (ErrorAux /= 0) then

                        write(*,*) 
                        write(*,*) 'Property time step is not a multiple of model time step.'
                        stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR80'
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
                     ClientModule = 'ModuleRunoffProperties',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR10'


        call GetData(NewProperty%Evolution%AdvDiff%NumericStability, &
                     Me%ObjEnterData, iflag,                         &
                     SearchType   = FromBlock,                       &
                     keyword      = 'ADVDIFF_NUM_STABILITY',         &
                     Default      = .FALSE.,                         &
                     ClientModule = 'ModuleRunoffProperties',   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR20'

        call GetData(NewProperty%Evolution%AdvDiff%SchmidtNumberH, &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_SCHMIDT_NUMBER_H',    &
                     Default      = 1.0,                           &
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR30'

        call GetData(NewProperty%Evolution%AdvDiff%SchmidtCoefV,   &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_SCHMIDT_COEF_V',      &
                     Default      = 1.0,                           &
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR40'

        call GetData(NewProperty%Evolution%AdvDiff%SchmidtBackgroundV, &
                     Me%ObjEnterData, iflag,                           &
                     SearchType   = FromBlock,                         &
                     keyword      = 'ADVDIFF_SCHMIDT_BACKGROUND_V',    &
                     Default      = 0.,                                &
                     ClientModule = 'ModuleRunoffProperties',     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR50'

        call GetData(NewProperty%Evolution%AdvDiff%NullDif,        &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                     &
                     keyword      = 'ADVDIFF_NULLDIF',             &
                     Default      = .false.,                       &
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR60'

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
!                     ClientModule = 'ModuleRunoffProperties', &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR70'
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
!            stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR80'
!
!        NewProperty%Evolution%AdvDiff%BoundaryCondition = BoundaryCondition

        !By default the horizontal Diffusion discretization is explicit
        NewProperty%Evolution%AdvDiff%DiffusionH_imp_exp  = ExplicitScheme

        NewProperty%Evolution%AdvDiff%ImplicitH_Direction = DirectionX

        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodH,     &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ADVDIFF_METHOD_H',            &
                     Default      = UpwindOrder1,                  &
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR90'

        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationH, &
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ADVDIFF_TVD_LIMIT_H',         &
                     Default      = Superbee,                      &
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR100'

        call GetData(NewProperty%Evolution%AdvDiff%VolumeRelMax,   &
                     Me%ObjEnterData, iflag,                       &
                     Keyword      = 'ADVDIFF_VOLUME_RELATION_MAX', &
                     Default      = 5.,                            &
                     SearchType   = FromBlock,                      &
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR130'


        if (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodH == P2_TVD) then
            NewProperty%Evolution%AdvDiff%Upwind2H = .true.
        else
            NewProperty%Evolution%AdvDiff%Upwind2H = .false.
        endif


        if (.not. Me%ComputeOptions%AdvDiff_Explicit .and.&
           (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3)) then

            write(*,*) 'If the advection of mass in the horizontal is implicit'
            write(*,*) 'the advection method can not be a second or third order upwind'
            stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR140.'

        endif


    end subroutine ReadAdvectionDiffusionParam

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyDiffusivity (NewProperty)

        !Arguments---------------------------------------------------------
        type(T_Property),    pointer                :: NewProperty

        !Local-------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: STAT_CALL
        
        !Begin-------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        
       
!            allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
!            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR30'
!            NewProperty%Diff_Turbulence_H = 0.

            
        allocate (NewProperty%ViscosityU (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR60'
        
        NewProperty%ViscosityU  = 0.0
        
        allocate (NewProperty%ViscosityV (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR70'
    
        NewProperty%ViscosityV  = 0.0   
        

!        allocate (NewProperty%Viscosity (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR20'
!        NewProperty%Viscosity         = 0.
!        
!        allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR30'
!        NewProperty%Diff_Turbulence_H = 0.

              


    end subroutine ConstructPropertyDiffusivity

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
        
    subroutine Construct_PropertyValues(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),              pointer      :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, i, j

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        integer                                     :: ILB,IUB
        integer                                     :: JLB,JUB
        integer                                     :: WorkSizeILB, WorkSizeIUB
        integer                                     :: WorkSizeJLB, WorkSizeJUB
        real                                        :: BottomInitialConc
        !Begin-----------------------------------------------------------------
        
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        WorkSizeILB = Me%WorkSize%ILB
        WorkSizeIUB = Me%WorkSize%IUB
        WorkSizeJLB = Me%WorkSize%JLB
        WorkSizeJUB = Me%WorkSize%JUB

        allocate(NewProperty%Concentration(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR00'
        NewProperty%Concentration(:,:) = FillValueReal

!        allocate(NewProperty%Mass(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR15'
!        NewProperty%Mass(:,:) = FillValueReal
        
        !Concentration After Transport. is the concentration that DN needs because is the one
        !used to remove flux from runoff
        allocate(NewProperty%ConcentrationAT(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR10'
        NewProperty%Concentration(:,:) = FillValueReal


        allocate(NewProperty%ConcentrationOld(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR20'
        NewProperty%ConcentrationOld(:,:) = FillValueReal

      
        if (Me%ExtVar%CoupledDN) then
            allocate(NewProperty%ConcentrationDN(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR22'
            NewProperty%ConcentrationDN(:,:) = FillValueReal
        endif

        call GetData(NewProperty%Pesticide,                                                 &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'PESTICIDE',                                            &
                     default      = .false.,                                                &
                     ClientModule = 'ModuleRunoffProperties',                               &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR23'

        !by default zero because user may copy property blocks and leave pesticide ON in a property that is not pesticide.
        if (NewProperty%Pesticide) then
            allocate(NewProperty%PesticideFlux(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR24'
            NewProperty%PesticideFlux(:,:) = 0.
        endif


        !it has to be always allocated
        allocate(NewProperty%ConcInInterfaceDN(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR25'
        NewProperty%ConcInInterfaceDN(:,:) = FillValueReal         
        
        !Currently is not needed since water only leaves runoff and concentration can only be from runoff
        !Here for future use
        allocate(NewProperty%ConcInBoundary(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR25'
        NewProperty%ConcInBoundary(:,:) = FillValueReal  
      
        call GetData(NewProperty%MinValue,                                                  &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'MIN_VALUE',                                            &
                     ClientModule = 'ModuleRunoffProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR26'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
            Me%Coupled%MinConcentration = .true.
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        if(NewProperty%Evolution%MinConcentration)then
            allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)&
                stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR27'
            NewProperty%Mass_Created(:,:) = 0.
        endif

        call GetData(NewProperty%WarnOnNegativeValues,                                   &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType     = FromBlock,                                         &
                     keyword        = 'WARN_ON_NEGATIVE_VALUES',                         &
                     Default        = .false.,                                           &                     
                     ClientModule   = 'ModuleRunoffProperties',                          &
                     STAT           = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR28' 
        
        if (NewProperty%WarnOnNegativeValues) Me%Coupled%WarnOnNegativeValues = .true.

        if (NewProperty%Evolution%BottomFluxes) then
            
            !in Runoff Properties if a propertie is particulated needs to have bottom processes
            !because if all the water infiltrates it has to go somewhere and can not enter soil
            !so needs to be deposited
            
            allocate(NewProperty%BottomConcentration(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR29'
            NewProperty%BottomConcentration   (:,:) = FillValueReal 

            allocate(NewProperty%TotalConcentration(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR30'
            NewProperty%TotalConcentration   (:,:) = FillValueReal
            
            !Bottom Initial Concentration
            call GetData(BottomInitialConc,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOTTOM_CONC',                                  &
                         Default      =  0.0,                                           & 
                         ClientModule = 'ModuleRunoffProperties',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR31' 
            
            NewProperty%BottomConcentration   (:,:) = BottomInitialConc 

            !Bottom Initial Concentration
            call GetData(NewProperty%BottomMinConc,                                     &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BOTTOM_MIN_CONC',                              &
                         Default      =  0.0,                                           & 
                         ClientModule = 'ModuleRunoffProperties',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR32' 

            if (NewProperty%Evolution%Erosion) then
                
                Me%Coupled%ErosionFluxes = .true.
                
                allocate (NewProperty%ErosionRate    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR33'
                NewProperty%ErosionRate = 0.0
                
                !Critial Erosion Shear Stress [Pa] - read for a grid in rourine ConstructData2D
!                call GetData(NewProperty%ErosionCriticalShear,                              &
!                             Me%ObjEnterData, iflag,                                        &
!                             SearchType   = FromBlock,                                      &
!                             keyword      = 'CRIT_SS_EROSION',                              &
!                             Default      =  0.2,                                           & 
!                             ClientModule = 'ModuleRunoffProperties',                        &
!                             STAT         = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR37' 

                !Erosion Coefficient [kg m-2 s-1] - in routine ConstrucData2D
!                call GetData(NewProperty%ErosionCoefficient,                                &
!                             Me%ObjEnterData, iflag,                                        &
!                             SearchType   = FromBlock,                                      &
!                             keyword      = 'EROSION_COEF',                                 &
!                             Default      =  5.0E-4,                                        & 
!                             ClientModule = 'ModuleRunoffProperties',                        &
!                             STAT         = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR38' 

            end if

            if (NewProperty%Evolution%Deposition) then
                
                Me%Coupled%DepositionFluxes = .true.

                allocate (NewProperty%DepositionRate (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR34'
                NewProperty%DepositionRate= 0.0
                
                allocate (NewProperty%Ws             (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR35'
                NewProperty%Ws = 0.0            
           
                !Critial Deposition Shear Stress [Pa] - in routine ConstrucData2D
!                call GetData(NewProperty%DepositionCriticalShear,                           &
!                             Me%ObjEnterData, iflag,                                        &
!                             SearchType   = FromBlock,                                      &
!                             keyword      = 'CRIT_SS_DEPOSITION',                           &
!                             Default      =  0.1,                                           & 
!                             ClientModule = 'ModuleRunoffProperties',                       &
!                             STAT         = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - Construct_PropertyEvolution - ERR50' 
        
                !See ModuleFreeVerticalMovement - Hindered settling  - CHS - kg m-3
                call GetData(NewProperty%CHS,                                               &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'CHS',                                          &
                             Default      =  4.0,                                           & 
                             ClientModule = 'ModuleRunoffProperties',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ConstructPropertyValues - ERR30'                 
                
                !Settling type: WSConstant = 1, SPMFunction = 2
                !Compute FallVelocity 
                call GetData(NewProperty%Ws_Type,                                           &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'WS_TYPE',                                      &
                             Default      =  WSConstant,                                    & 
                             ClientModule = 'ModuleRunoffProperties',                        &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ConstructPropertyValues - ERR40' 

                if (NewProperty%Ws_Type.EQ.WSConstant) then
                    !See ModuleFreeVerticalMovement - Constant settling velocity [m s-1]
                    call GetData(NewProperty%Ws_Value,                                          &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromBlock,                                      &
                                 keyword      = 'WS_VALUE',                                     &
                                 Default      =  0.0001,                                        & 
                                 ClientModule = 'ModuleRunoffProperties',                        &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ConstructPropertyValues - ERR50' 

                else
                    !See ModuleFreeVerticalMovement
                    call GetData(NewProperty%KL,                                                &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromBlock,                                      &
                                 keyword      = 'KL',                                           &
                                 Default      =  0.1,                                           & 
                                 ClientModule = 'ModuleRunoffProperties',                        &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ConstructPropertyValues - ERR60' 
                   
                    !See ModuleFreeVerticalMovement
                    call GetData(NewProperty%KL1,                                               &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromBlock,                                      &
                                 keyword      = 'KL1',                                          &
                                 Default      =  0.1,                                           & 
                                 ClientModule = 'ModuleRunoffProperties',                        &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ConstructPropertyValues - ERR70' 

                    !See ModuleFreeVerticalMovement
                    call GetData(NewProperty%ML,                                                &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromBlock,                                      &
                                 keyword      = 'ML',                                           &
                                 Default      =  4.62,                                          & 
                                 ClientModule = 'ModuleRunoffProperties',                        &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ConstructPropertyValues - ERR80' 

                    !See ModuleFreeVerticalMovement
                    call GetData(NewProperty%M,                                                 &
                                 Me%ObjEnterData, iflag,                                        &
                                 SearchType   = FromBlock,                                      &
                                 keyword      = 'M',                                            &
                                 Default      =  1.0,                                           & 
                                 ClientModule = 'ModuleRunoffProperties',                        &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ModuleRunoffProperties - ConstructPropertyValues - ERR90'
                endif
            endif
            
            if (NewProperty%Evolution%SplashErosion) then
                Me%Coupled%SplashErosionFluxes = .true.
                
                allocate(NewProperty%SplashRate(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR110'
                NewProperty%SplashRate = 0.0

               
            endif
            
            
        endif
            

        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleRunoffProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR120'
          
        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not. NewProperty%Old) then

            !Get water points
            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR130'

            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExtVar%BasinPoints,            &
                                       Matrix2D             = NewProperty%Concentration,        &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR140'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR0150'
            end if

            call GetRunoffWaterColumn     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR151'            

            call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                     GridCellArea = Me%ExtVar%Area,                            & 
                                     STAT = STAT_CALL )    
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR0151.5'
            
            !initial concentration based on initial water column
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then 
                    if (Me%ExtVar%WaterColumn(i,j) .gt. 0.0) then
                        NewProperty%ConcentrationOld(i,j) = NewProperty%Concentration(i,j)
                        !g = g/m3 * m * m2
!                        NewProperty%Mass(i,j)             = NewProperty%Concentration(i,j) * Me%ExtVar%InitialWaterColumn &
!                                                            * Me%ExtVar%Area(i,j)
                    else
                        NewProperty%Concentration(i,j)    = 0.0
                        NewProperty%ConcentrationOld(i,j) = NewProperty%Concentration(i,j)
!                        NewProperty%Mass(i,j)             = 0.0 
                    endif
                    if (NewProperty%Evolution%BottomFluxes) then
                        !kg/m2 = (g/m3 * (m * m2) * 1E-3 kg/g) / m2 + kg/m2
                        NewProperty%TotalConcentration(i,j) = NewProperty%Concentration(i,j)          &
                                                              * Me%ExtVar%WaterColumn(i,j) * 1E-3   &
                                                              + NewProperty%BottomConcentration(i,j)
                    endif
                endif
            enddo
            enddo
            
!            call SetMatrixValue(NewProperty%ConcentrationOld, Me%Size, NewProperty%Concentration,Me%ExtVar%BasinPoints)

            call UnGetRunoff     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR152'            

            call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR0153'

            call CheckFieldConsistence (NewProperty)

            call UnGetBasin(Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR160'

        else

            ! If the property is old then the program is going to try to find a property
            ! with the same name in the Water properties initial file written in HDF format  
            call ReadOldConcBoundariesHDF(NewProperty)

        end if   

    end subroutine Construct_PropertyValues

      !--------------------------------------------------------------------------


    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),    pointer        :: NewProperty

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------


        call GetData(NewProperty%TimeSerie,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'TIME_SERIE',                                        &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR01'
        

        call GetData(NewProperty%BoxTimeSerie,                                           &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'BOX_TIME_SERIE',                                    &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR02'

        if (NewProperty%BoxTimeSerie) then
            Me%Output%Boxes_ON = .true.
            Me%NumberPropForBoxes = Me%NumberPropForBoxes + 1
        endif

        call GetData(NewProperty%BoxTimeSerie2D,                                           &
                     Me%ObjEnterData, iflag,                                               &
                     Keyword      = 'BOX_TIME_SERIE2D',                                    &
                     Default      = .false.,                                               &
                     SearchType   = FromBlock,                                             &
                     ClientModule = 'ModuleRunoffProperties',                              &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR03'



        call GetData(NewProperty%OutputHDF,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'OUTPUT_HDF',                                        &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR04'
        
    end subroutine Construct_PropertyOutPut
   
   !---------------------------------------------------------------------------

    subroutine ConstructData2D

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: i,j, STAT_CALL
        !Begin-----------------------------------------------------------------
        

!        call ConstructScalar2D(Me%Disper_Longi, ExtractType = FromBlock,        &
!                              block_begin = '<begin_dispersion_long>',          &
!                              block_end   = '<end_dispersion_long>')
!        
        call Read_Property_2D(Me%Disper_Trans, FromBlock,  '<begin_dispersion_trans>', '<end_dispersion_trans>')            
        
        if (Me%Coupled%BottomFluxes .and. Me%Coupled%ErosionFluxes) then
            
            call Read_Property_2D (Me%ErosionCriticalShear, FromBlock, "<begin_critical_shear_erosion>",     &
                                    "<end_critical_shear_erosion>")

            call Read_Property_2D (Me%ErosionCoefficient, FromBlock, "<begin_erosion_coefficient>", "<end_erosion_coefficient>")
        
        endif


        if (Me%Coupled%BottomFluxes .and. Me%Coupled%DepositionFluxes) then
            
            call Read_Property_2D (Me%DepositionCriticalShear, FromBlock, "<begin_critical_shear_deposition>", &
                                    "<end_critical_shear_deposition>")
       
        endif  
        
        
        if (Me%Coupled%BottomFluxes .and. Me%Coupled%SplashErosionFluxes) then
            
!            call Read_Property_2D (Me%StoneFraction, FromBlock, "<begin_Soil_StoneFraction>", "<end_Soil_StoneFraction>")
!            call Read_Property_2D (Me%ClayFraction, FromBlock, "<begin_Soil_ClayFraction>", "<end_Soil_ClayFraction>")
            call Read_Property_2D (Me%KClay, FromBlock, "<begin_soil_detach>", "<end_soil_detach>") 
                 
        endif 
        
        if (Me%Coupled%BottomFluxes .and. Me%Coupled%ErosionFluxes .and. Me%Coupled%DepositionFluxes) then

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructData2D - ModuleRunoffProperties - ERR020'
        
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
                    if (Me%DepositionCriticalShear%Field(i,j) >= Me%ErosionCriticalShear%Field(i,j)) then
                        write (*,*) 'critical shear for erosion must be higher than critical shear for deposition in cell ', i , j
                        stop 'ModuleRunoffProperties - ConstructData2D - ERR10' 
                    endif
                endif
            enddo
            enddo

            call UnGetBasin                   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructData2D - ModuleRunoffProperties - ERR040'
        
        endif
    
    
    end subroutine ConstructData2D

   !---------------------------------------------------------------------------

    subroutine ConstructPartition

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: Property
        type(T_Property), pointer                   :: DissolvedProperty
        type(T_Property), pointer                   :: CohesiveSediment
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
                if (STAT_CALL /= SUCCESS_) then
                    write(*,*)
                    write(*,*) 'Particulate property was not found'
                    stop 'ConstructPartition - ModuleRunoffProperties - ERR10'
                else
                    if(.not. Check_Particulate_Property (Couple_ID)) then
                        write(*,*)
                        write(*,*) 'Couple property', trim(GetPropertyName(Couple_ID))
                        write(*,*) 'is not recognized by the model as being particulate'
                        stop 'ConstructPartition - ModuleRunoffProperties - ERR15'                        
                    endif
                endif
                
                !Check if coupled names of each other match
                if(trim(ParticulateProperty%Evolution%Partition%Couple) .ne.                &
                   trim(DissolvedProperty%ID%Name)) then
                    write(*,*)
                    write(*,*) 'Particulate property pair and dissolved property names do not match'
                    stop 'ConstructPartition - ModuleRunoffProperties - ERR20'
                endif
                
                if(ParticulateProperty%Evolution%DTInterval .ne.                            &
                   DissolvedProperty%Evolution%DTInterval) then
                    write(*,*)
                    write(*,*) 'Particulate property and dissolved property DTs do not match'
                    stop 'ConstructPartition - ModuleRunoffProperties - ERR30'
                endif
                
                if(DissolvedProperty%Evolution%Partition%Rate /=                            &
                   ParticulateProperty%Evolution%Partition%Rate   )then
                    write(*,*)'Particulate and dissolved phases must have equal partition rates'
                    stop 'ConstructPartition - ModuleRunoffProperties - ERR40'
                end if

                TotalPartition = DissolvedProperty%Evolution%Partition%Fraction  +          &
                                 ParticulateProperty%Evolution%Partition%Fraction

                Error = abs(1. - TotalPartition)
                     
                if(Error > 0.001)then
                    write(*,*)'Particulate and dissolved phases fractions must sum iqual to 1.'
                    stop 'ConstructPartition - ModuleRunoffProperties - ERR50'
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
                           ParticulateProperty%Evolution%Partition%UseSedimentRefConc)) then
                    write(*,*)
                    write(*,*) 'Both particulate property and dissolved property must use or not sediment ref conc'
                    stop 'ConstructPartition - ModuleRunoffProperties - ERR80'
                endif

                if (DissolvedProperty%Evolution%Partition%UseSedimentRefConc) then
                    call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then
                        write(*,*) 'Using sediment ref conc needs cohesive sediment property'
                        stop 'ConstructPartition - ModuleRunoffProperties - ERR03'
                    endif
                endif
                
                if(DissolvedProperty%Evolution%Partition%SedimentRefConc /=                 &
                   ParticulateProperty%Evolution%Partition%SedimentRefConc   )then
                    write(*,*)'Particulate and dissolved phases must have equal cohesive sediment'
                    write(*,*)'reference concentration'
                    stop 'ConstructPartition - ModuleRunoffProperties - ERR90'
                end if


            end if

            nullify(DissolvedProperty, ParticulateProperty)

            Property => Property%Next


        end do do1
          
        nullify(Property)


    end subroutine ConstructPartition

    !--------------------------------------------------------------------------

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
        
        !Counts the number of Properties which has timeserie option set to true
        PropertyX => Me%FirstProperty
        nProperties = 0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                nProperties = nProperties + 1
                if(PropertyX%Evolution%BottomFluxes) then
                    nProperties = nProperties + 2
                    if(PropertyX%Evolution%Erosion) then
                        nProperties = nProperties + 1
                    endif
                    if(PropertyX%Evolution%Deposition) then
                        nProperties = nProperties + 2
                    endif 
                    if(PropertyX%Evolution%SplashErosion) then
                        nProperties = nProperties + 1
                    endif                          
                endif
            endif
            PropertyX => PropertyX%Next
        enddo
        if (Me%Coupled%BottomFluxes .and. (Me%Coupled%ErosionFluxes .or. Me%Coupled%DepositionFluxes)) then
            nProperties = nProperties +1
        endif
        if(Me%Coupled%SplashErosionFluxes) then
            nProperties = nProperties +1
        endif
        
        !Allocates PropertyList
        allocate(PropertyList(nProperties))
        
        !Property names
        n=1
        PropertyX  => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                PropertyList(n)  = trim(PropertyX%ID%Name)//" mg/l"
                n=n+1
                if(PropertyX%Evolution%BottomFluxes) then
                    PropertyList(n)  = trim(PropertyX%ID%Name)//" Bottom kg/m2"
                    n=n+1
                    PropertyList(n)  = trim(PropertyX%ID%Name)//" Total kg/m2"
                    n=n+1
                    if(PropertyX%Evolution%Erosion) then
                        PropertyList(n)  = trim(PropertyX%ID%Name)//" ErosionRate kg.m-2.s-1"
                        n=n+1
                    endif
                    if(PropertyX%Evolution%Deposition) then
                        PropertyList(n)  = trim(PropertyX%ID%Name)//" DepositionRate kg.m-2.s-1"
                        n=n+1
                        PropertyList(n)  = trim(PropertyX%ID%Name)//" FallVelocity m/s"
                        n=n+1
                    endif  
                    if(PropertyX%Evolution%SplashErosion) then
                        PropertyList(n)  = trim(PropertyX%ID%Name)//" SplashRate kg.m-2.s-1"
                        n=n+1
                    endif                                             
                endif
            endif
            PropertyX=>PropertyX%Next
        enddo
        if (Me%Coupled%BottomFluxes .and. (Me%Coupled%ErosionFluxes .or. Me%Coupled%DepositionFluxes)) then
            PropertyList(n)  = "Shear Stress N/m2"
            n=n+1
        endif
        if(Me%Coupled%SplashErosionFluxes) then
            PropertyList(n)  = "Rain Kinetic Work W/m2"
        endif
        
        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleRunoffProperties',                           &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR01' 

        if (iflag == 1) then
            Me%OutPut%TimeSerie_ON = .true.
        else
            Me%OutPut%TimeSerie_ON = .false.
        endif
        
        !Get water points
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR020'


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srrp",                                       &
                            WaterPoints2D = Me%ExtVar%BasinPoints,                      &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR030' 

        !Unget
        call UnGetBasin                   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR040'

        !Deallocates PropertyList
        deallocate(PropertyList)

        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR03'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR04'
            
            call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR04'
            
i1:         if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR05'

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR06'

                    if (IgnoreOK) then
                        write(*,*) 'Time Serie outside the domain - ',trim(TimeSerieName)
                        cycle
                    else
                        stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR07'
                    endif

                endif


                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR08'

            endif i1

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      LocalizationI   = Id,                             &
                                      LocalizationJ   = Jd,                             & 
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR09'

            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR090'

            if (Me%ExtVar%BasinPoints(Id, Jd) /= WaterPoint) then
                 write(*,*) 'Time Serie in a cell outside basin - ',trim(TimeSerieName)
            endif
            
            call UnGetBasin                (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR0140'


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
                stop 'ConstructHDF - ModuleRunoffProperties - ERR01' 

            if (Me%OutPut%HDF_ON) then

                Me%OutPut%NextOutPut = 1

                call Open_HDF5_OutPut_File

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructHDF - ModuleRunoffProperties - ERR02'
            endif 

            !Output for restart
            call GetOutPutTime(Me%ObjEnterData,                                             &
                               CurrentTime  = Me%ExtVar%Now,                                &
                               EndTime      = Me%ExtVar%EndTime,                                   &
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
        integer                                             :: ILB, IUB, JLB, JUB
        logical                                             :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        type(T_Property    ),                       pointer :: PropertyX
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList
        character(len=StringLength), dimension(:),  pointer :: FluxesOutputList
        integer                                             :: nScalars, n, nFluxes

        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        ! This keyword have two functions if exist fluxes between boxes are compute 
        ! and the value read is the name file where the boxes are defined
        call GetData(Me%Files%BoxesFile,                                            &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'BOXFLUXES',                                    &
                     ClientModule = 'ModuleRunoffProperties',                       &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR01'
        if (iflag .EQ. 0)                                                           &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR02'    
        
        inquire(File = Me%Files%BoxesFile, Exist = exist)
        if (exist) then
            inquire(File = Me%Files%BoxesFile, Opened  = Opened)
            if (opened) then
                write(*,*    ) 
                write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Files%BoxesFile))
                write(*,*    ) 'Already opened.'
                stop           'StartOutputBoxFluxes - ModuleRunoffProperties - ERR03'    
            end if
        else
            write(*,*) 
            write(*,*)     'Could not find the boxes file.'
            write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
            stop           'StartOutputBoxFluxes - ModuleRunoffProperties - ERR04'    
        end if
        
        !Output Rates and Properties inside box (.bxm)
        nScalars = Me%NumberPropForBoxes
        !Output Properties fluxes between boxes
        nFluxes  = Me%NumberPropForBoxes
            
        allocate(ScalarOutputList(nScalars), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR05'

        allocate(FluxesOutputList(nFluxes), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR06'

        n = 0
        PropertyX => Me%FirstProperty
        do while(associated(PropertyX))
            if (PropertyX%BoxTimeSerie) then
                n = n + 1
                ScalarOutputList(n) = "runoff_"//trim(PropertyX%ID%Name)
                FluxesOutputList(n) = "runoff_"//trim(PropertyX%ID%Name)
            endif
            PropertyX => PropertyX%Next
        end do        
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR01'

        call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                 &
                         TimeID             = Me%ObjTime,                   &
                         HorizontalGridID   = Me%ObjHorizontalGrid,         &
                         BoxesFilePath      = Me%Files%BoxesFile,           &
                         FluxesOutputList   = FluxesOutputList,             &
                         ScalarOutputList   = ScalarOutputList,             &
                         WaterPoints2D      = Me%ExtVar%BasinPoints,        &
                         STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR07'

        !Unget
        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR08'  


        deallocate(ScalarOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR09'

        deallocate(FluxesOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR10'


        allocate(Me%CellMass(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleRunoffProperties - ERR170'
        Me%CellMass(:,:) = 0.            
        
        
        
    end subroutine StartOutputBoxFluxes

    !--------------------------------------------------------------------------

     subroutine Open_HDF5_OutPut_File        

        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        !Begin-----------------------------------------------------------------

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Gets a pointer to Topography
        call GetGridData        (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR00'

        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR01'



        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR02'

      
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR02'


        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR05'
        
        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR05'

        !WriteBasinPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR07'

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR08'  


        !Unget
        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR90'  

        !UnGets Topography
        call UnGetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR100'


    end subroutine Open_HDF5_OutPut_File   
   
    !--------------------------------------------------------------------------

    subroutine ConstructAsciiFile            

        !Local-----------------------------------------------------------------
        integer               :: status
        integer               :: STAT_CALL                
        integer               :: Counter
        character(LEN=4)      :: Number

        call UnitsManager(Me%Files%AsciiUnit, OPEN_FILE, STAT = status) 
        if (status /= SUCCESS_) stop "ConstructAsciiOutPut - ModulePorousMediaProperties - ERR01"

        Counter  = 1
do1:     do
            Number = '    '
            write(Number, fmt='(i4)')Counter
            open(UNIT   = Me%Files%AsciiUnit,                                      &
                 FILE   = '..\res\RP_ADCoefs_'//trim(adjustl(Number))//'.log', &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do1
            else
                Counter = Counter + 1
            end if
        enddo do1
        
        write (Me%Files%AsciiUnit, FMT=*) 'YY     MM   DD   HH   MM     SS     i   j    cofA_U cofB cofC_U cofA_V cofC_V'
    
    end subroutine ConstructAsciiFile
    
    !--------------------------------------------------------------------------
   
    subroutine ReadOldConcBoundariesHDF(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        character (Len=StringLength)                :: PropertyName
        logical                                     :: EXIST
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ
        logical, save                               :: FileNameRead = .false.                                     

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        !----------------------------------------------------------------------

        !Reads the name of the file where to read restart options
        if (.not. FileNameRead) then

            call ReadFileName ('RUNOFF_PROP_INI', Me%Files%InitialFile, "Runoff Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR00'

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
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR01'


            PropertyName = trim(adjustl(NewProperty%ID%name))

            NewProperty%Concentration(:,:) = FillValueReal
!            NewProperty%Mass(:,:) = FillValueReal

            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR02'

            call HDF5ReadData   (ObjHDF5, "/Results/"//trim(adjustl(NewProperty%ID%Name)), &
                                 trim(adjustl(NewProperty%ID%Name)),                       &
                                 Array2D = NewProperty%Concentration,                    &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR03'

            
            if (NewProperty%Evolution%BottomFluxes) then
                call HDF5ReadData   (ObjHDF5, "/Results/"//trim(adjustl(NewProperty%ID%Name))//" Bottom",   &
                                     trim(adjustl(NewProperty%ID%Name))//" Bottom",                         &
                                     Array2D = NewProperty%Concentration,                    &
                                     STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                   &
                    stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR03'
            endif


!            call HDF5ReadData   (ObjHDF5, "/Concentration/"//NewProperty%ID%Name,        &
!                                 NewProperty%ID%Name,                                    &
!                                 Array2D = NewProperty%Mass,                             &
!                                 STAT    = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_)                                                   &
!                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR03'

            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR06'

        else
            
            write(*,*)
            stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR07'

        end if cd0

    end subroutine ReadOldConcBoundariesHDF


    !--------------------------------------------------------------------------


    subroutine CheckFieldConsistence(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer               :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: i,j
        integer                                 :: ILB, IUB, JLB, JUB
        logical                                 :: StopSubroutine = .false.
        integer                                 :: UnitAux
        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB


            
        !Verification if the values read are lower than zero in water points
        do I = ILB, IUB
        do J = JLB, JUB
            
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                               
                if (NewProperty%Concentration(i, j) < 0.) then
                    
                    StopSubroutine = .true.
                    NewProperty%Concentration(i, j) = 0.0

                endif

            else

                NewProperty%Concentration(i, j) = FillValueReal

            endif

        enddo
        enddo

        if (StopSubroutine) then                                                   
            
            call UnitsManager(UnitAux, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR02' 

            open(UnitAux, FILE = trim(NewProperty%ID%name)//'.new',                 &
                 FORM = 'FORMATTED', STATUS = 'UNKNOWN', IOSTAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR03' 

            write(UnitAux,*) '<ConcentrationBegin>'
           
            do I = ILB, IUB
            do J = JLB, JUB
            
                write(UnitAux,*) NewProperty%Concentration(i, j)

            enddo
            enddo

            write(UnitAux,*) '<ConcentrationEnd>'

            call UnitsManager(UnitAux, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR04' 

            write(*,*) 'A new concentration file was created for property: ', trim(NewProperty%ID%Name)
            write(*,*) 'Run again with this new file ', trim(NewProperty%ID%name)//'.new'
            stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR05'  

        endif

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------


!    subroutine CoupleSoilQuality        
!
!        !Local-----------------------------------------------------------------
!        type(T_Property), pointer                           :: PropertyX
!        integer, pointer, dimension(:)                      :: SoilQualityPropertyList
!        integer                                             :: STAT_CALL
!        real                                                :: SoilQualityDT
!        integer                                             :: nProp = 0 
!
!        !Begin------------------------------------------------------------------
!
!        !Counts the number of Properties which has WaterQuality option set to true
!        PropertyX => Me%FirstProperty
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilQuality) then
!                nProp = nProp + 1
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        !Allocates Array to hold IDs
!        allocate (SoilQualityPropertyList(1:nProp))
!
!        !Fills Array
!        PropertyX => Me%FirstProperty
!        nProp = 0
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilQuality) then
!                nProp = nProp + 1
!                SoilQualityPropertyList(nProp) = PropertyX%ID%IDNumber
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModuleRunoffProperties - ERR01'
!
!        !Start Interface
!        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
!                                TimeID              = Me%ObjTime,                    &
!                                SinksSourcesModel   = SedimentQualityModel,          &
!                                DT                  = SoilQualityDT,                 &
!                                PropertiesList      = SoilQualityPropertyList,       &
!                                WaterPoints2D       = Me%ExtVar%BasinPoints,         &
!                                Size2D              = Me%WorkSize,                   &
!                                STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                   &
!            stop 'CoupleSoilQuality - ModuleRunoffProperties - ERR02'
!
!
!        call UnGetBasin                   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModuleRunoffProperties - ERR03'
!
!
!        deallocate (SoilQualityPropertyList)
!
!        Me%Coupled%SoilQuality_DT          = SoilQualityDT 
!        Me%Coupled%SoilQuality_NextCompute = Me%ExtVar%Now    
!
!        nullify (Me%DissolvedToParticulate2D)
!        allocate(Me%DissolvedToParticulate2D(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
!        Me%DissolvedToParticulate2D(:,:) = null_real
!
!        Me%ResidualTime = 0.
!    
!    end subroutine CoupleSoilQuality
!
!    !--------------------------------------------------------------------------


!#ifdef _PHREEQC_
!    !--------------------------------------------------------------------------
!    subroutine CoupleSoilChemistry        
!
!        !Local-----------------------------------------------------------------
!        type(T_Property), pointer                           :: PropertyX
!        integer, pointer, dimension(:)                      :: SoilChemistryPropertyList
!        integer                                             :: STAT_CALL
!        real                                                :: SoilChemistryDT
!        integer                                             :: nProp = 0 
!
!        !Begin------------------------------------------------------------------
!
!        !Counts the number of Properties which has SoilChemistry option set to true
!        PropertyX => Me%FirstProperty
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilChemistry) then
!                nProp = nProp + 1
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        !Allocates Array to hold IDs
!        allocate (SoilChemistryPropertyList(1:nProp))
!
!        !Fills Array
!        PropertyX => Me%FirstProperty
!        nProp = 0
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilChemistry) then
!                nProp = nProp + 1
!                SoilChemistryPropertyList(nProp) = PropertyX%ID%IDNumber
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        !Question: What does this function? Is it necessary to SoilChemistry process or Interface?
!        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModuleRunoffProperties - ERR01'
!
!        !Start Interface
!        call ConstructInterface(InterfaceID         = Me%ObjInterfaceSoilChemistry,  &
!                                TimeID              = Me%ObjTime,                    &
!                                SinksSourcesModel   = PhreeqCModel,                  &
!                                DT                  = SoilChemistryDT,               &
!                                PropertiesList      = SoilChemistryPropertyList,     &
!                                WaterPoints2D       = Me%ExtVar%BasinPoints,         &
!                                Size3D              = Me%WorkSize,                   &
!                                STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                   &
!            stop 'CoupleSoilChemistry - ModuleRunoffProperties - ERR02'
!
!        !Question: What does this function? 
!        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModuleRunoffProperties - ERR03'
!
!        deallocate (SoilChemistryPropertyList)
!
!        Me%Coupled%SoilChemistry_DT          = SoilChemistryDT 
!        Me%Coupled%SoilChemistry_NextCompute = Me%ExtVar%Now    
!
!            
!    end subroutine CoupleSoilChemistry
!    !--------------------------------------------------------------------------
!#endif


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
            STAT_  = NOT_FOUND_ERR_  
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine Search_Property

    !--------------------------------------------------------------------------

    subroutine Read_Property_2D(Property, ExtractType, BeginBlock, EndBlock, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Property_2D)                :: Property
        character(len=*)                    :: BeginBlock, EndBlock
        integer                             :: ExtractType
        integer, optional                   :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        logical                             :: BlockFound
        integer                             :: LocalClientNumber

        !Begin-----------------------------------------------------------------

        select case(ExtractType)

            case(FromBlock)
                
                call ExtractBlockFromBuffer (Me%ObjEnterData, LocalClientNumber,    &
                                             BeginBlock, EndBlock,                  &
                                             BlockFound, STAT = STAT_CALL)
            case(FromBlockInBlock)
                
                call ExtractBlockFromBlock (Me%ObjEnterData, ClientNumber,          &
                                            BeginBlock, EndBlock,                   &
                                            BlockFound, STAT = STAT_CALL)
            case default

        end select


        if(present(ClientNumber)) LocalClientNumber = ClientNumber
        
        if (BlockFound) then

            !Allocates Variables
            allocate (Property%Field (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            Property%Field (:,:) = null_real

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR01'        


            call ConstructFillMatrix  (PropertyID           = Property%ID,                      &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = ExtractType,                      &
                                       PointsToFill2D       = Me%ExtVar%BasinPoints,            &
                                       Matrix2D             = Property%Field,                   &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR05'

            call UngetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR06'  

            call GetDefaultValue(Property%ID%ObjFillMatrix, Property%Scalar, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR10'

            call GetIfMatrixRemainsConstant(FillMatrixID    = Property%ID%ObjFillMatrix,        &
                                            RemainsConstant = Property%Constant,                &
                                            STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR20'


            if(.not. Property%ID%SolutionFromFile)then
                call KillFillMatrix(Property%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR30'
            end if

            if(ExtractType == FromBlock)then
                call Block_Unlock(Me%ObjEnterData, LocalClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR40'
            end if

        else

            write (*,*)'Block ',trim(BeginBlock),' ',trim(EndBlock),' not found'
            stop 'Read_Property_2D - ModuleRunoffProperties - ERR60'

        endif
        
        
        select case(ExtractType)

            case(FromBlock)

                call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR70'

            case(FromBlockInBlock)
                
                call RewindBlock  (Me%ObjEnterData, LocalClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleRunoffProperties - ERR80'

            case default

        end select


    end subroutine Read_Property_2D
    
    !---------------------------------------------------------------------------
    


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
!
    !---------------------------------------------------------------------------

!    subroutine GetRPCoupled(RunoffPropertiesID, &
!                             SoilQuality,             &
!#ifdef _PHREEQC_
!                             SoilChemistry,           &
!#endif                             
!                             STAT) 
!
!        !Arguments-------------------------------------------------------------
!        integer                        :: RunoffPropertiesID
!        integer, optional, intent(OUT) :: STAT
!        logical, optional, intent(OUT) :: SoilQuality        
!#ifdef _PHREEQC_
!        logical, optional, intent(OUT) :: SoilChemistry
!#endif
!
!        !External--------------------------------------------------------------
!
!        integer :: ready_              
!
!        !Local-----------------------------------------------------------------
!
!        integer :: STAT_              !Auxiliar local variable
!        !----------------------------------------------------------------------
!
!        STAT_ = UNKNOWN_
!
!        call Ready(RunoffPropertiesID, ready_)
!        
!cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
!            (ready_ .EQ. READ_LOCK_ERR_)) then
!
!            if (present(SoilQuality   )) SoilQuality    = Me%Coupled%SoilQuality
!            
!#ifdef _PHREEQC_            
!            if (present(SoilChemistry )) SoilChemistry  = Me%Coupled%SoilChemistry
!#endif
!
!            STAT_ = SUCCESS_
!        else 
!            STAT_ = ready_
!        end if cd1
!
!        if (present(STAT)) STAT = STAT_
!
!        !----------------------------------------------------------------------
!
!    end subroutine GetRPCoupled
!    !--------------------------------------------------------------------------
    
    subroutine GetRPConcentration(RunoffPropertiesID, ConcentrationX, &  !MassX, 
                                PropertyXIDNumber,                    &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: RunoffPropertiesID
        real, pointer, dimension(:,:)               :: ConcentrationX
!        real, pointer, dimension(:,:), optional     :: MassX
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

        call Ready(RunoffPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mRunoffPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                
                ConcentrationX => PropertyX%Concentration
!                if (present(MassX)) MassX => PropertyX%Mass

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
            
    end subroutine GetRPConcentration

    !--------------------------------------------------------------------------------

    subroutine GetRPConcentrationAT(RunoffPropertiesID, ConcentrationXAT, &  !MassX, 
                                PropertyXIDNumber,                    &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: RunoffPropertiesID
        real, pointer, dimension(:,:)               :: ConcentrationXAT
!        real, pointer, dimension(:,:), optional     :: MassX
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

        call Ready(RunoffPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mRunoffPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                
                ConcentrationXAT => PropertyX%ConcentrationAT
!                if (present(MassX)) MassX => PropertyX%Mass

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
            
    end subroutine GetRPConcentrationAT

    !--------------------------------------------------------------------------------

    subroutine GetRPConcentrationOld(RunoffPropertiesID, ConcentrationXOld, &  !MassX, 
                                PropertyXIDNumber,                    &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: RunoffPropertiesID
        real, pointer, dimension(:,:)               :: ConcentrationXOld
!        real, pointer, dimension(:,:), optional     :: MassX
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

        call Ready(RunoffPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mRunoffPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                
                ConcentrationXOld => PropertyX%ConcentrationOld
!                if (present(MassX)) MassX => PropertyX%Mass

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
            
    end subroutine GetRPConcentrationOld

    !--------------------------------------------------------------------------------

    subroutine GetRPDecayRate(RunoffPropertiesID, DecayRate,              &  
                                PropertyXIDNumber,                        &
                                STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: RunoffPropertiesID
        real                                        :: DecayRate          !day-1
        integer,                      intent(IN )   :: PropertyXIDNumber
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_CALL              
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: STAT_    

        !------------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                
                DecayRate = PropertyX%Evolution%DecayRate

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetRPDecayRate

    !--------------------------------------------------------------------------------

!    subroutine SetWindVelocity (RunoffPropertiesID, WindModulus, STAT)
!                                  
!        !Arguments--------------------------------------------------------------
!        integer                                     :: RunoffPropertiesID
!        real, dimension(:,:), pointer               :: WindModulus
!
!        integer, optional, intent(OUT)              :: STAT
!
!        !Local-----------------------------------------------------------------
!        integer                                     :: ready_        
!        integer                                     :: STAT_
!        
!        !----------------------------------------------------------------------
!
!        STAT_ = UNKNOWN_
!
!        call Ready(RunoffPropertiesID, ready_)
!
!        if (ready_ .EQ. IDLE_ERR_)then
!            
!
!            Me%ExtVar%WindVelocity2D   => WindModulus
!
!
!            STAT_ = SUCCESS_
!        else
!            STAT_ = ready_
!        end if
!
!        if (present(STAT))STAT = STAT_
!
!    end subroutine SetWindVelocity 

    !---------------------------------------------------------------------------

   
    subroutine SetDNConcRP (RunoffPropertiesID, PropertyID, DNConcentration, ChannelsID, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer                                         :: PropertyID   
        real, dimension (:), pointer                    :: DNConcentration
        integer, dimension(:, :), pointer               :: ChannelsID
        integer                                         :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
        
            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetDNConcRP - ModuleRunoffProperties - ERR01'
            
           
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
                write(*,*) 'Looking for Drainage Network Property in Runoff Properties', GetPropertyName(PropertyID)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetDNConcRP - ModuleRunoffProperties - ERR010'
            end if

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetDNConcRP - ModuleRunoffProperties - ERR040'               


        else
            STAT_ = ready_
        end if

        STAT = STAT_        
                     
    end subroutine SetDNConcRP

    !---------------------------------------------------------------------------

    subroutine SetBasinConcRP   (RunoffPropertiesID, BasinConcentration, & !BasinMass, &
                                                PropertyXIDNumber, MassToBottom,   & 
                                                STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        real,    dimension(:, :), pointer               :: BasinConcentration
!        real,    dimension(:, :), pointer, optional     :: BasinMass
        integer                                         :: PropertyXIDNumber
        integer, intent(OUT), optional                  :: STAT
        real(8), dimension(:, :), pointer               :: MassToBottom

        !Local------------------------------------------------------------------
        integer                                         :: j ,i
        integer                                         :: STAT_, STAT_CALL, ready_
        type(T_Property), pointer                       :: PropertyX
        real(8)                                         :: WaterVolume

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR01'            
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_CALL)

            if (STAT_CALL == SUCCESS_) then

                call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                         GridCellArea = Me%ExtVar%Area,                            & 
                                         STAT = STAT_CALL )    
                if (STAT_CALL /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR010'

                call GetRunoffWaterColumn     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR20'
            
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
                        
                        PropertyX%Concentration (i,j) = BasinConcentration (i,j)
!                        PropertyX%Mass (i,j)          = BasinMass (i,j)
                        
                        if (PropertyX%Evolution%BottomFluxes) then
                            !Kg/m2 = ((Kg/m2 * m2) + (g * 1E-3kg/g)) / m2 
                            PropertyX%BottomConcentration(i,j) = (PropertyX%BottomConcentration(i,j) * Me%ExtVar%Area(i, j)   &
                                                                  + MassToBottom(i,j) * 1E-3) / Me%ExtVar%Area(i, j)
                            WaterVolume = Me%ExtVar%WaterColumn(i,j) * Me%ExtVar%Area(i, j)
                            ![kg/m2] = [g/m3]* [m3] * [1E-3kg/g] /[m2] + [kg/m2]
                            PropertyX%TotalConcentration (i,j) = PropertyX%Concentration (i,j) * 1E-3 * WaterVolume            &
                                                                 / Me%ExtVar%Area(i, j) &
                                                                 + PropertyX%BottomConcentration (i,j) 
                        endif
                    endif
                enddo
                enddo            

                call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR040'
        
                call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_)   
                if (STAT_ /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR020'
                    
            else
                write(*,*) 'Looking for Runoff Property in Runoff Property ???', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetBasinConcRP - ModuleDrainageNetwork - ERR010'
            end if

            call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR030'
            
            STAT_ = SUCCESS_
            
        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetBasinConcRP 

    !---------------------------------------------------------------------------


    subroutine SetBasinToRPSplash   (RunoffPropertiesID, ThroughFall,         &
                                      CanopyDrainage, CanopyHeight,  STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        real(8),    dimension(:, :), pointer            :: CanopyDrainage, ThroughFall
        real,       dimension(:, :), pointer            :: CanopyHeight
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetBasinToRPSplash - ModuleRunoffProperties - ERR01'            

            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then

                    Me%ExtVar%ThroughFall(i,j)    = ThroughFall(i,j)
                    Me%ExtVar%CanopyHeight(i,j)   = CanopyHeight(i,j)
                    Me%ExtVar%CanopyDrainage(i,j) = CanopyDrainage(i,j)
                endif
            enddo
            enddo

            call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetBasinToRPSplash - ModuleRunoffProperties - ERR010'
                
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetBasinToRPSplash 

    !---------------------------------------------------------------------------

    subroutine SetVegetationRP   (RunoffPropertiesID, FertilOrganicNParticFluff, &
                                      FertilOrganicPParticFluff, PesticideIDNumber, &
                                      PesticideSoilFlux, VegetationDT, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        real,    dimension(:, :), pointer, optional     :: FertilOrganicNParticFluff
        real,    dimension(:, :), pointer, optional     :: FertilOrganicPParticFluff
        integer, optional                               :: PesticideIDNumber
        real,    dimension(:, :), pointer, optional     :: PesticideSoilFlux
        real,    optional                               :: VegetationDT
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j
        type (T_Property), pointer                      :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
            if (present(FertilOrganicNParticFluff)) then
                call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'SetVegetationRP - ModuleRunoffProperties - ERR01'            

                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then

                        Me%ExtVar%FertilOrganicNParticFluff(i,j)   = FertilOrganicNParticFluff(i,j)
                    endif
                enddo
                enddo

                call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'SetVegetationRP - ModuleRunoffProperties - ERR010'
            endif
            
            if (present(FertilOrganicPParticFluff)) then
                call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'SetVegetationRP - ModuleRunoffProperties - ERR020'            

                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then

                        Me%ExtVar%FertilOrganicPParticFluff(i,j)   = FertilOrganicPParticFluff(i,j)
                    endif
                enddo
                enddo

                call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
                if (STAT_ /= SUCCESS_) stop 'SetVegetationRP - ModuleRunoffProperties - ERR030'            
            endif
            
            if (present(VegetationDT)) Me%ExtVar%VegetationDT = VegetationDT

            if (present(PesticideIDNumber       )) then
                call SearchProperty(PropertyX, PropertyXIDNumber = PesticideIDNumber, STAT = STAT_)
                if (STAT_ == SUCCESS_) then
                    if (.not. present(PesticideSoilFlux)) stop 'SetVegetationRP - ModuleRunoffProperties - ERR40'
                    if (.not. associated(PropertyX%PesticideFlux)) then
                        write(*,*) 'Property', GetPropertyName(PesticideIDNumber) 
                        write(*,*) 'is a pesticide and needs PESTICIDE keyword in Runoff Properties'
                        stop 'SetVegetationRP - ModuleRunoffProperties - ERR50'
                    else    
                        PropertyX%PesticideFlux => PesticideSoilFlux
                    endif
                else
                    write(*,*)
                    write(*,*) 'Not found pesticide', GetPropertyName(PesticideIDNumber) 
                    write(*,*) 'in Runoff Properties property list'
                    stop ' SetVegetationRP - ModuleRunoffProperties - ERR60'
                endif    
            endif
                
            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetVegetationRP 

    !---------------------------------------------------------------------------

    subroutine GetRPMassBalance (RunoffPropertiesID,                        &
                                  PropertyID,                               &
                                  TotalStoredMass,                          &
                                  DNExchangeMass,                           &
                                  STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer                                         :: PropertyID
        real(8), optional                               :: TotalStoredMass
        real(8), optional                               :: DNExchangeMass
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type (T_Property), pointer                      :: PropertyX
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then 
                if (present (TotalStoredMass)) TotalStoredMass = PropertyX%MB%TotalStoredMass
                if (present (DNExchangeMass))  DNExchangeMass  = PropertyX%MB%DNExchangeMass
                
                STAT_CALL = SUCCESS_
            else
                stop 'GetRPMassBalance - ModuleRunoffProperties - ERR01'
            endif
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetRPMassBalance

    !---------------------------------------------------------------------------

    subroutine GetRPnProperties (RunoffPropertiesID, nProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer                                         :: nProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            nProperties       = Me%PropertiesNumber
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetRPnProperties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetRPPropertiesIDByIdx (RunoffPropertiesID, Idx, ID,PropAdvDiff, Particulate,Decay, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer, intent(IN)                             :: Idx
        integer, intent(OUT)                            :: ID
        logical, intent(OUT), optional                  :: PropAdvDiff
        logical, intent(OUT), optional                  :: Particulate
        logical, intent(OUT), optional                  :: Decay
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, i
        type (T_Property), pointer                      :: CurrProp

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrProp => Me%FirstProperty
            do i = 1, idx - 1
                CurrProp => CurrProp%Next
            enddo

            ID        = CurrProp%ID%IDNumber
            if (present(PropAdvDiff)) PropAdvDiff = CurrProp%Evolution%AdvectionDiffusion
            if (present(Particulate)) Particulate = CurrProp%Particulate
            if (present(Decay      )) Decay       = CurrProp%Evolution%Decay

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetRPPropertiesIDByIdx

    !---------------------------------------------------------------------------

    subroutine CheckRPProperty (RunoffPropertiesID,                         &
                                  PropertyID,                               &
                                  STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer                                         :: PropertyID
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, STAT_
        type (T_Property), pointer                      :: PropertyX
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then 
                if (.not. PropertyX%Evolution%AdvectionDiffusion) then
                    write(*,*)
                    write(*,*)'Property', GetPropertyName(PropertyID)
                    write(*,*)'has advection diffusion inactive in RunoffProperties Module'
                    write(*,*)'and it is unconsistent with activation in other Modules'
                    stop 'CheckRPProperty - ModuleRunoffProperties - ERR01' 
                else               
                    STAT_CALL = SUCCESS_
                endif
            else
                write(*,*)
                write(*,*)'Could not find property', GetPropertyName(PropertyID)
                write(*,*)'in RunoffProperties Module'
                stop 'CheckRPProperty - ModuleRunoffProperties - ERR010'
            endif
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine CheckRPProperty

    !---------------------------------------------------------------------------

    subroutine GetRPOptions (RunoffPropertiesID, SplashErosion, Erosion, Deposition, Decay, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        logical, optional                               :: SplashErosion, Erosion, Deposition
        logical, optional                               :: Decay
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            if (present(SplashErosion)) then
                SplashErosion       = Me%Coupled%SplashErosionFluxes
            endif
            if (present(Erosion)) then
                Erosion             = Me%Coupled%ErosionFluxes
            endif
            if (present(Deposition)) then
                Deposition          = Me%Coupled%DepositionFluxes
            endif

            if (present(Decay)) then
                Decay          = Me%Coupled%Decay
            endif
            
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetRPOptions

    !---------------------------------------------------------------------------

    subroutine UnGetRunoffProperties2D_I(ObjRunoffPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID
        integer, dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRunoffProperties_, Me%InstanceID, "UnGetRunoffProperties2D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunoffProperties2D_I

    !--------------------------------------------------------------------------

    subroutine UnGetRunoffProperties2D_R4(ObjRunoffPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRunoffProperties_, Me%InstanceID,  "UnGetRunoffProperties2D_R4")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunoffProperties2D_R4


    !--------------------------------------------------------------------------

    subroutine UnGetRunoffProperties2D_R8(ObjRunoffPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID
        real(8), dimension(:, :), pointer                :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRunoffProperties_, Me%InstanceID,  "UnGetRunoffProperties2D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunoffProperties2D_R8


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyRunoffProperties(ObjRunoffPropertiesID,          &
                                           STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunoffPropertiesID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_,STAT_CALL
        type (T_Property), pointer                  :: PropertyX

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunoffProperties - ModuleRunoffProperties - ERR02'
            
            !Actualize the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunoffProperties - ModuleRunoffProperties - ERR03'
                      
            call ReadLockExternalVar

            !Actualize property matrix
            PropertyX => Me%FirstProperty

            do while (associated(PropertyX))

                call SetMatrixValue (PropertyX%ConcentrationOld, Me%Size, PropertyX%Concentration, Me%ExtVar%BasinPoints)

                if (Me%CheckGlobalMass) then
                    PropertyX%MB%TotalStoredMass     = 0.0
                    PropertyX%MB%DNExchangeMass      = 0.0
                    PropertyX%MB%TotalDischargeMass  = 0.0
                endif
                
                PropertyX => PropertyX%Next
                
            enddo

            
            !Eduardo Jauch
            !Actualize properties if evolution from file
            call ActualizePropertiesFromFile

            !Nutrient sources from vegetation - fertilization particulate to fluff layer 
            call InterfaceFluxes
            
            !Water discharges fluxes - explicit because may exist more than one discharge in the same cell
            if (Me%ComputeOptions%Discharges) call ModifyDischarges
            
            !Now adevection diffusion has only the transport between cells after discharge took place
            !Discharge flow is accounted for "Old" volume and advection/diffusion new volume is the one after transport
            if (Me%Coupled%AdvectionDiffusion) then
                call AdvectionDiffusionProcesses
            endif

            if (Me%Coupled%MinConcentration)     call SetLimitsConcentration  
            if (Me%Coupled%WarnOnNegativeValues) call WarnOnNegativeValues    ('After AdectionDiffusion')
            
            !Actualize property matrix after transport
            PropertyX => Me%FirstProperty

            do while (associated(PropertyX))
                !Concentration after transport that is the one used for DN flux, and needed in ModuleDrainageNetwork
                !Becuase DN flux could remove all water and concentration would be zero
                !but DN flux occurs with concentration after transport and not with conc after transport and DN flux
                call SetMatrixValue (PropertyX%ConcentrationAT, Me%Size, PropertyX%Concentration, Me%ExtVar%BasinPoints)
                PropertyX => PropertyX%Next
            enddo
            
            !!!WARNING!!!
            !The processes associated to Storm Water Model and Storm Water Drainage are not accounted in properties
            !In the future need to add them
            !!!WARNING!!!
            
            !DN fluxes had to be separated from Advection Diffusion (explained above)
            !routine to compute new explicit concentration
            if (Me%ExtVar%CoupledDN) then
                call ModifyDrainageNetworkInterface
            endif
            
            !Route DFourPoints fluxes need to be updated
            if (Me%ExtVar%RouteDFour) then
                call ModifyRouteDFourInterface
            endif
            
            !Boundary fluxes need to be updated (e.g. for mass balance in the future)
            if (Me%ExtVar%BoundaryImposed) then
                call ModifyBoundaryInterface
            endif
            
!            if (Me%Coupled%SoilQuality) then
!                call SoilQualityProcesses
!            endif
!
!#ifdef _PHREEQC_
!            if (Me%Coupled%SoilChemistry) then
!                call SoilChemistryProcesses
!            endif
!#endif            
            if (Me%Coupled%BottomFluxes) then
                call ModifyBottomFluxes
            endif

            if (Me%Coupled%Decay) then
                call DecayProcesses
            endif
            
            if (Me%Coupled%Partition)then
                call Partition_Processes 
            endif           
            
            if (Me%Coupled%MinConcentration)     call SetLimitsConcentration 
            if (Me%Coupled%WarnOnNegativeValues) call WarnOnNegativeValues    ('After Quality Processes')
            
            if (Me%Output%Timeserie_ON) then
                call OutPut_TimeSeries
            endif

            if (Me%Output%HDF_ON) then
                call OutPut_HDF
            endif

            if (Me%Output%Boxes_ON) then
                call Output_Boxes_Mass
            endif

            !Restart Output
            if (Me%Output%WriteRestartFile .and. .not. (Me%ExtVar%Now == Me%ExtVar%EndTime)) then
                if(Me%ExtVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    call WriteFinalFile
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                endif
            endif
            
            if (Me%ComputeOptions%DTIntervalAssociated) then
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

    end subroutine ModifyRunoffProperties

    !-----------------------------------------------------------------------------

    subroutine ModifyDischarges ()

        !Arguments--------------------------------------------------------------
        !Local------------------------------------------------------------------
        integer                                 :: iDis, nDischarges
        integer                                 :: i, j, k, n
        integer                                 :: STAT_CALL
        real, dimension(:, :), pointer          :: FlowDischarge
        real,    dimension(:,:), pointer        :: DischargesConc
        type (T_Property), pointer              :: Property
        integer                                 :: iProp
        real(8)                                 :: VolumeOld
        integer, dimension(:),   pointer        :: VectorI, VectorJ
        integer                                 :: nCells
        integer                                 :: FlowDistribution

        !Get integrated flow from runoff to be sure using same values
        !ATTENTION - DIFFERENT DISCHARGES CAN NOT OCCUR AT THE SAME CELL
        !Test this here in RP because is only problem if properties are used
        !or create duplications of integrated discharge (per discharge) in runoff        
        call GetFlowDischarge (Me%ObjRunoff, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOffProperties - ModifyDischarges - ERR01'
        
        !Get concentration
        !Gets the number of discharges
        call GetDischargesNumber(Me%ObjDischarges, nDischarges, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOffProperties - ModifyDischarges - ERR010'
        
        allocate(DischargesConc(nDischarges, Me%nPropWithDischarge))
        DischargesConc = null_real
        
        do iDis = 1, nDischarges

            call GetDischargesGridLocalization(Me%ObjDischarges,                        &
                                               DischargeIDNumber = iDis,                &
                                               Igrid = i,                               &
                                               JGrid = j,                               &
                                               KGrid = k,                               &
                                               STAT = STAT_CALL)
            if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOffProperties - ModifyDischarges - ERR020'
            
            if (k == 0) then

                call GetDischargeFlowDistribuiton(Me%ObjDischarges, iDis, nCells, FlowDistribution, &
                                                  VectorI, VectorJ, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_) stop 'ModulePorousMediaProperties - ModifyDischarges - ERR040'             
               
 dn:            do n=1, nCells
 
                    if (nCells > 1) then
                        i         = VectorI(n)
                        j         = VectorJ(n)                      
                    endif

                    if (Me%ExtVar%BasinPoints(i, j) /= WaterPoint)  Cycle                
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
                                    stop 'ModuleRunOffProperties - ModifyDischarges - ERR030'
                                endif
                            endif
                        
                            !In case of negative discharge flux for mass balance is done using old concentration
                            !and before concentration is updated in routine DischargeProperty
                            !Do not move this computation to after DischargeProperty
                            !In case of positive use dicharge concentration
                            if (Me%CheckGlobalMass) then
                                if (FlowDischarge(i,j) .lt. 0.0) then                        
                                    !kg = kg + m3/s * s * g/m3 * 1e-3kg/g
                                    Property%MB%TotalDischargeMass = Property%MB%TotalDischargeMass + (FlowDischarge(i,j)   &
                                                                     * Me%ExtVar%DT * Property%ConcentrationOld(i,j))
                                else
                                    !kg = kg + m3/s * s * g/m3 * 1e-3kg/g
                                    Property%MB%TotalDischargeMass = Property%MB%TotalDischargeMass + (FlowDischarge(i,j)  &
                                                                     * Me%ExtVar%DT * DischargesConc(iDis, iProp))
                            
                                endif
                            endif
                        
                            !initial volume in runoff - no discarge
                            VolumeOld = Me%ExtVar%WaterColumnOld(i,j) * Me%ExtVar%Area(i,j)
                        
                            !Update old concentration (same as doing new concentrarion and then old = new before other processes)
                            call DischargeProperty (FlowDischarge(i,j), DischargesConc(iDis, iProp),        &
                                                    i, j, VolumeOld,   Property, Me%ExtVar%DT, .false.)
                        
                        end if
                                    
                        Property => Property%Next

                    enddo
                    
                enddo dn
 
                call UnGetDischarges(Me%ObjDischarges, VectorI, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModuleRunOffProperties - ModifyDischarges - ERR070'

                call UnGetDischarges(Me%ObjDischarges, VectorJ, STAT = STAT_CALL)             
                if (STAT_CALL/=SUCCESS_)                                                    &
                    stop 'ModuleRunOffProperties - ModifyDischarges - ERR080'
           
           endif
           
        enddo
        
        deallocate (DischargesConc)
        
        call UngetRunoff (Me%ObjRunoff, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOffProperties - ModifyDischarges - ERR040'

    end subroutine ModifyDischarges  
    
    !--------------------------------------------------------------------------

    subroutine DischargeProperty (DischargeFlow, DischargeConc, Igrid, Jgrid, VolumeOld,         &
                                   Property, LocalDT, Accumulate)
        !Arguments--------------------------------------------------------------
        real                                        :: DischargeFlow, DischargeConc
        real(8)                                     :: VolumeOld
        integer                                     :: Igrid, Jgrid
        type (T_Property), pointer                  :: Property
        real                                        :: LocalDT
        logical                                     :: Accumulate

        !Local------------------------------------------------------------------
        real(8)                                     :: DischargeVolume
        real(8)                                     :: OldMass, NewMass
        real                                        :: Concentration
        
        Concentration = Property%ConcentrationOld(Igrid,Jgrid) 

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

                !If the property acculumlates in the water column 
                !(e.g particulate properties during infiltration) then the concentration will increase

                if (Accumulate) then
                    NewMass          = OldMass
                else
                    NewMass          = OldMass * (1.0 + DischargeVolume / VolumeOld)
                endif
                
                !if water remains
                if (abs(DischargeVolume) < VolumeOld) then
                   
                   ![g/m3] = [g] / [m3]
                    Concentration    = NewMass / (VolumeOld + DischargeVolume)
                
                else   !if all water exits node than accumulated mass needs to be accounted in bottom!
                    
                    Concentration  = 0.0
                    
                    if (Accumulate) then
                        ![kg/m2] = [kg/m2] + [g] * 1e-3 [kg/g] / m2
                        Property%BottomConcentration(Igrid,Jgrid) = Property%BottomConcentration(Igrid,Jgrid) +    &
                                                                    (NewMass * 1e-3 / Me%ExtVar%Area(Igrid,Jgrid))
                    endif
                endif
            endif

        else
            
            !Do Nothing            

        endif
                    
        !Update concOld with discharge
        Property%ConcentrationOld(Igrid,Jgrid) = Concentration


    end subroutine DischargeProperty

    !---------------------------------------------------------------------------

    subroutine InterfaceFluxes !routine for mass sources - from vegetation (OM fertilization to fluff layer)
                               !and pesticide if water column exists
        !Local--------------------------------------------------------------------
        !Begin--------------------------------------------------------------------

        if (Me%ExtVar%CoupledVegetation) then
            if(Me%ExtVar%VegParticFertilization .or. Me%ExtVar%Pesticide) then
                call VegetationInterfaceFluxes
            endif
        endif


    end subroutine InterfaceFluxes
 
    !-----------------------------------------------------------------------------

    subroutine VegetationInterfaceFluxes !routine for mass sources - from vegetation (OM fertilization to fluff layer)
        !Local--------------------------------------------------------------------
        type (T_Property), pointer                :: Property
        integer                                   :: STAT_CALL, i,j
        real                                      :: ModelDT, VegDT
        real                                      :: CellWaterVolume
        !Begin--------------------------------------------------------------------

        
        call GetComputeTimeStep (Me%ObjTime, ModelDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VegetationInterfaceFluxes - ModuleRunoffProperties - ERR010'
        
        VegDT = Me%ExtVar%VegetationDT
        
        if(Me%ExtVar%VegParticFertilization) then
        
            call Search_Property(Property, PON_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VegetationInterfaceFluxes - ModuleRunoffProperties - ERR01'
                
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
                
                    !m3             = m  * m2
                    !CellWaterVolume = Me%ExtVar%WaterColumnOld(i,j) * Me%ExtVar%Area(i,j) 
                    !kgN/m2                          = kg/m2  + (kg/ha * dt/vegdt * 1ha/10000m2)  
                    Property%BottomConcentration(i,j) = Property%BottomConcentration(i,j)                          &
                                                       + (Me%ExtVar%FertilOrganicNParticFluff(i,j)                     &
                                                        * (ModelDT / VegDT)  / 10000.)   
                                                    
                endif
            enddo
            enddo
     
            call Search_Property(Property, POP_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VegetationInterfaceFluxes - ModuleRunoffProperties - ERR01'
                
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
                
                    !m3             = m  * m2
                    !CellWaterVolume = Me%ExtVar%WaterColumnOld(i,j) * Me%ExtVar%Area(i,j) 
                    !kgN/m2                          = kg/m2  + (kg/ha * dt/vegdt * 1ha/10000m2)  
                    Property%BottomConcentration(i,j) = Property%BottomConcentration(i,j)                          &
                                                       + (Me%ExtVar%FertilOrganicPParticFluff(i,j)                     &
                                                        * (ModelDT / VegDT)  / 10000.)   
                                                    
                endif
            enddo
            enddo
        endif
        
        if (Me%ExtVar%Pesticide) then
            !Pesticide flux for runoff wate if water column exists
            Property => Me%FirstProperty

            do while (associated(Property)) 
                
                if (Property%Pesticide) then
                
                    !if (Property%Evolution%AdvectionDiffusion) then

                        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                            
                            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint .and. Me%ExtVar%WaterColumnOld(i,j) .gt. AlmostZero) then
                            
                                !m3             = m  * m2
                                CellWaterVolume = Me%ExtVar%WaterColumnOld(i,j) * Me%ExtVar%Area(i,j) 
                                !g/m3                          = g/m3  + (kg/ha * dt/vegdt * 1E3g/kg * (m2) * 1ha/10000m2) / m3 
                                Property%ConcentrationOld(i,j) = Property%ConcentrationOld(i,j) + (Property%PesticideFlux(i,j) &
                                                                    * ModelDT / VegDT * 1e3 * Me%ExtVar%Area(i,j) / 10000.)        &
                                                                     / CellWaterVolume
                                                                
                            endif
                        enddo
                        enddo
    !                    else
                      
    !                    endif         
                endif
            
                Property => Property%Next
                
            enddo  
        endif

    end subroutine VegetationInterfaceFluxes
 
    !-----------------------------------------------------------------------------
    
    subroutine ModifyDiffusivity_Old(PropertyX)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer           :: PropertyX

        !External--------------------------------------------------------------
        integer                              :: i, j,  CHUNK

        !Begin----------------------------------------------------------------------

        call ModifyTurbulence(PropertyX)

        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J)

               
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                !Vertical Diffusivity
                !m2/s = m2/s * [-]  + m2/s 
                PropertyX%Diffusivity(i,j) = 0.0
                
                !Horizontal Diffusivity is called viscosity to maintain the format of when is called in Module Advection Diffusion
                !m2/s = m2/s * [-]  + m2/s
                PropertyX%Viscosity(i,j)   = PropertyX%Evolution%AdvDiff%Molecular_Diff_Coef &
                                             + PropertyX%Diff_Turbulence_H (i,j)
                
                    PropertyX%ViscosityU(i,j) = PropertyX%Viscosity(i,j) 
                    PropertyX%ViscosityV(i,j) = PropertyX%Viscosity(i,j)
                    
                                                                       

            endif

        enddo
        enddo

        !$OMP END DO
        !$OMP END PARALLEL


    end subroutine ModifyDiffusivity_Old

    !--------------------------------------------------------------------------

    subroutine ModifyTurbulence (PropertyX)

        !Arguments------------------------------------------------------------------
        type (T_Property), pointer :: PropertyX

        !Local----------------------------------------------------------------------
        integer                    :: i, j, CHUNK
!        real                       :: VelMedU, VelMedV, VelMedW, VelU, VelV
        real                       :: VelU, VelV

        !Begin----------------------------------------------------------------------

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J, VelU, VelV)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
!                VelMedU = 0.5 * (Me%ExtVar%UnsatU(i,j,k) * Me%ExtVar%ComputeFacesU3D(I,J,K) +  Me%ExtVar%UnsatU(i,j+1,k)  &
!                           * Me%ExtVar%ComputeFacesU3D(I,J+1,K))
!                VelMedV = 0.5 * (Me%ExtVar%UnsatV(i,j,k) * Me%ExtVar%ComputeFacesV3D(I,J,K) +  Me%ExtVar%UnsatV(i+1,j,k)   &
!                            * Me%ExtVar%ComputeFacesV3D(I+1,J+,K))
! 

!                !m2/s = m * m/s
!                PropertyX%Diff_Turbulence_H(i,j,k) = Me%Disper_Longi%Field(i,j,k) *  abs(0.5 * (VelMedU + VelMedV)                &
!                                                    + Me%Disper_Trans%Field(i,j,k) * abs(VelMedW)
                VelU = 0.5 * (Me%ExtVar%CenterVelU(i,j)+Me%ExtVar%CenterVelU(i,j-1))
                VelV = 0.5 * (Me%ExtVar%CenterVelV(i,j)+Me%ExtVar%CenterVelV(i-1,j))
                !m2/s = m * m/s
                PropertyX%Diff_Turbulence_H(i,j) = Me%Disper_Longi%Field(i,j) *                              &
                                                         abs(VelV + VelU  / 2.)
               
            endif

        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        

    end subroutine ModifyTurbulence

    !---------------------------------------------------------------------------

    subroutine ActualizePropertiesFromFile
    
        !Local--------------------------------------------------------------------        
        type (T_Property), pointer :: PropertyX 
        integer                    :: STAT_CALL   
        
        !-------------------------------------------------------------------------    

        !Properties can not be defined from file - RunoffProperties is not a forcing module
        !The concentrations have to be computed or it will not be consistent with transport
!        PropertyX => Me%FirstProperty
!
!        do while (associated(PropertyX))
!
!            if (PropertyX%ID%SolutionFromFile) then
!            
!                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
!                                       Matrix2D       = PropertyX%Concentration,    &
!                                       PointsToFill2D = Me%ExtVar%BasinPoints,    &
!                                       STAT           = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR01'
!            
!            endif
!            
!            PropertyX => PropertyX%Next
!            
!        enddo
        
        
        if (Me%Coupled%BottomFluxes .and. Me%Coupled%ErosionFluxes) then
             
             if (Me%ErosionCriticalShear%ID%SolutionFromFile) then

                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,       &
                                       Matrix2D       = Me%ErosionCriticalShear%Field,    &
                                       PointsToFill2D = Me%ExtVar%BasinPoints,            &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR10'
                
             endif


             if (Me%ErosionCoefficient%ID%SolutionFromFile) then

                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,       &
                                       Matrix2D       = Me%ErosionCoefficient%Field,    &
                                       PointsToFill2D = Me%ExtVar%BasinPoints,            &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR20'
                
             endif
         
         endif
        
        if (Me%Coupled%BottomFluxes .and. Me%Coupled%DepositionFluxes) then
             
             if (Me%DepositionCriticalShear%ID%SolutionFromFile) then

                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,       &
                                       Matrix2D       = Me%DepositionCriticalShear%Field,    &
                                       PointsToFill2D = Me%ExtVar%BasinPoints,            &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR30'
                
             endif
         
         endif    
         
        if (Me%Coupled%BottomFluxes .and. Me%Coupled%SplashErosionFluxes) then
             
!             if (Me%StoneFraction%ID%SolutionFromFile) then
!
!                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,       &
!                                       Matrix2D       = Me%StoneFraction%Field,           &
!                                       PointsToFill2D = Me%ExtVar%BasinPoints,            &
!                                       STAT           = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR40'
!                
!             endif
             
!             if (Me%ClayFraction%ID%SolutionFromFile) then
!
!                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,       &
!                                       Matrix2D       = Me%ClayFraction%Field,    &
!                                       PointsToFill2D = Me%ExtVar%BasinPoints,            &
!                                       STAT           = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR50'
!                
!             endif
             
             if (Me%KClay%ID%SolutionFromFile) then

                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,       &
                                       Matrix2D       = Me%KClay%Field,    &
                                       PointsToFill2D = Me%ExtVar%BasinPoints,            &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR50'
                
             endif            
         
         endif      
        !-------------------------------------------------------------------------    
    
    end subroutine ActualizePropertiesFromFile
    
    !-------------------------------------------------------------------------    

    subroutine AdvectionDiffusionProcesses

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        
        !begin--------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "AdvectionDiffusionProcesses")
        

       !Compute water volume. Remove discharges because they were explicity computed
        call ComputeVolumes
        
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
               
!                !Update property mass fluxes between modules
!                if (Me%CheckGlobalMass) then
!                    call ModifyInterfaceMassFluxes (PropertyX)
!                endif
            
            endif
                        
            PropertyX => PropertyX%Next

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "AdvectionDiffusionProcesses")


    end subroutine AdvectionDiffusionProcesses

    !-----------------------------------------------------------------------------

    subroutine ComputeVolumes

        !Local-----------------------------------------------------------------
        integer                                   :: i, j, CHUNK, STAT_CALL        
        real, dimension(:,:), pointer             :: FlowDischarge

        !----------------------------------------------------------------------
        !The discharge flows have to be added because they are accounted separately 
        !see ModifyDischarges. Advection/Diffusion will onlye act between WaterColumnBT 
        !(Before Transport) and WaterColumnAT (After Transport) volumes

        !Get integrated flow from runoff to be sure using same values
        call GetFlowDischarge (Me%ObjRunoff, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOffProperties - ComputeVolumes - ERR01'
        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J)     
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             
                 
                !Me%WaterVolume(i,j)        = Me%ExtVar%WaterColumn(i,j) * Me%ExtVar%Area(i,j)
                !m3/s = m * m2 + (m3/s * s)
                Me%WaterVolume(i,j)        = Me%ExtVar%WaterColumnAT(i,j) * Me%ExtVar%Area(i,j)
                
                !Adding discharge flow (because is separated) for computing "Old" for tansport
                !m = m + (m3/s * s / m2)
                Me%WaterColumnBT(i,j)      = Me%ExtVar%WaterColumnOld(i,j) +   &
                                            (FlowDischarge(i,j) * Me%ExtVar%DT / Me%ExtVar%Area(i,j))

                
                if (Me%WaterVolume(i, j) > AlmostZero) then
                    Me%DummyOpenPoints(i,j) = 1
                else
                    Me%DummyOpenPoints(i,j) = 0
                endif

            endif
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL 

        call UngetRunoff (Me%ObjRunoff, FlowDischarge, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ModuleRunOffProperties - ComputeVolumes - ERR040'

        
    end subroutine ComputeVolumes

    !-----------------------------------------------------------------------------

    subroutine RestartVariables(PropertyX)
        
        !Argument-----------------------------------------------------------------
         type (T_Property), pointer                  :: PropertyX

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        !begin--------------------------------------------------------------------
        
        CurrProperty => PropertyX
        

        call SetMatrixValue (Me%COEFExpl%CoefInterfDN, Me%Size, 0.0)
        call SetMatrixValue (CurrProperty%ConcInInterfaceDN, Me%Size, 0.0)
        call SetMatrixValue (Me%COEFExpl%CoefInterfBoundary, Me%Size, 0.0)
        call SetMatrixValue (CurrProperty%ConcInBoundary, Me%Size, 0.0)
       
        call SetMatrixValue (Me%TICOEF3, Me%Size, 0.0)

        if (.not. Me%NewFormulation) then
       
            call SetMatrixValue (Me%COEF3_HorAdvXX%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%F_flux, Me%Size, 0.0)
            
            call SetMatrixValue (Me%COEF3_HorAdvYY%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%F_flux, Me%Size, 0.0)                
        
        endif    
        
        if (.not. Me%ComputeOptions%AdvDiff_Explicit) then
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

           

        Direction = 'Horizontal'
        
        select case (Me%AdvDiff_AdvMethodH)
        
            case (UpwindOrder1)
                
                call AdvectionUpwindOrder1 (Direction)
            
            case (CentralDif)
            
                call AdvectionCentralDifferences(Direction)
                
            case default
                
                write(*,*) 'Undefined method for horizontal advection check ADVDIFF_METHOD_H'
                stop 'AdvectionDiffusionProcesses - Runoff Properties - ERR01'
            
         end select
             
    
    end subroutine ComputeAdvectionTerms

    !-----------------------------------------------------------------------------

    subroutine AdvectionUpwindOrder1(Direction)
        
        !Argument-----------------------------------------------------------------
        character(len=StringLength)         :: Direction
        !Local--------------------------------------------------------------------
        integer                             :: i,j, CHUNK                             
       
        !begin-------------------------------------------------------------------- 

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "AdvectionUpwindOrder1")

        
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

            
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(i,j) 
     
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%BasinPoints(i, j) == 1) then
                    
                    !x direction
!                    if (Me%ExtVar%ComputeFaces2D(i,j) == 1) then
                        if (Me%ExtVar%FluxU(i,j) .gt. 0.0) then
                            Me%COEF3_HorAdvXX%D_flux(i,j) = Me%ExtVar%FluxU(i,j)
                        else
                            Me%COEF3_HorAdvXX%E_flux(i,j) = Me%ExtVar%FluxU(i,j)
                        endif    
!                    endif
                    
                    !y direction
!                    if (Me%ExtVar%ComputeFaces2D(i,j) == 1) then                
                        if (Me%ExtVar%FluxV(i,j) .gt. 0.0) then
                            Me%COEF3_HorAdvYY%D_flux(i,j) = Me%ExtVar%FluxV(i,j)
                        else
                            Me%COEF3_HorAdvYY%E_flux(i,j) = Me%ExtVar%FluxV(i,j)
                        endif   
 !                   endif                 
                    
                endif

            end do do1
            end do do2
            !$OMP END DO
                
            !$OMP END PARALLEL
            
        
        else
            stop 'AdvectionUpwindOrder1 - Runoff Properties - ERR010'
        endif

       
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "AdvectionUpwindOrder1")
           
           
    end subroutine AdvectionUpwindOrder1
    
    !-----------------------------------------------------------------------------

    subroutine AdvectionCentralDifferences(Direction)

        !Argument-----------------------------------------------------------------
        character(len=StringLength)         :: Direction
        
        !Local--------------------------------------------------------------------
        integer                             :: i,j, CHUNK   
        real                                :: Aux1, Aux2
        real, dimension(:,:  ), pointer     :: DUX, DVY
        !begin-------------------------------------------------------------------- 

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "AdvectionCentralDifferences")

        
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
            
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
            !$OMP PARALLEL PRIVATE(i,j, Aux1, Aux2) 
     
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :       do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%BasinPoints(i, j) == 1) then
                    
                    !x direction
 !                   if (Me%ExtVar%ComputeFaces2D(i,j) == 1) then
                        Aux1 = DUX(i,j  )/(DUX(i,j) + DUX(i,j-1))
                        Aux2 = DUX(i,j-1)/(DUX(i,j) + DUX(i,j-1))
                        Me%COEF3_HorAdvXX%D_flux(i,j) = Me%ExtVar%FluxU(i,j) * Aux1
                        Me%COEF3_HorAdvXX%E_flux(i,j) = Me%ExtVar%FluxU(i,j) * Aux2
 !                   endif
                    
                    !y direction
 !                   if (Me%ExtVar%ComputeFaces2D(i,j) == 1) then                
                        Aux1 = DVY(i  ,j)/(DVY(i,j) + DVY(i-1,j))
                        Aux2 = DVY(i-1,j)/(DVY(i,j) + DVY(i-1,j))
                        Me%COEF3_HorAdvYY%D_flux(i,j) = Me%ExtVar%FluxV(i,j) * Aux1
                        Me%COEF3_HorAdvYY%E_flux(i,j) = Me%ExtVar%FluxV(i,j) * Aux2
 !                   endif                 
                    
                    
                endif

            end do do1
            end do do2
            !$OMP END DO
                
            !$OMP END PARALLEL
        
        else
            stop 'AdvectionCentralDifferences - Runoff Properties - ERR010'
        
        endif
        
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "AdvectionCentralDifferences")
           
           
    end subroutine AdvectionCentralDifferences
    
    !-----------------------------------------------------------------------------

    subroutine ModifyCoefs(PropertyX)

        !Argument-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
       
        !begin--------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyCoefs")

        
        !In this subroutine are computed all the sources/sinks of water and mass
        !that exist in the Runoff Module (for new Water Column computation)
        
        !Diffusivity will be used in AdvectionDiffusion
        call ModifyDiffusivity_New(PropertyX)

        !Fluxes in X, Y direction explicit or implicit. 
        call ModifyAdvectionDiffusionCoefs(PropertyX)

        !Discharges
        !Discharges not yet accounted
        
        !Fluxes with Drainage network - in the cells that link with river
        !This link was disconnected since in one time step runoff may run out of water
        !to DN, and if it happens in several time steps than flow is different from zero
        !and conc is zero that is wrong. So transport will be evaluated for water columns
        !before and after transport and link to DN is made after (see ModifyDrainageNetworkInterface)
        !if (Me%ExtVar%CoupledDN) then
            !call ModifyDrainageNetworkCoefs(PropertyX)
        !endif

        !BoundaryFluxes - moved. same reason as for DN fluxes
!        if (Me%ExtVar%BoundaryImposed) then
!            call ModifyBoundaryCoefs(PropertyX)
!        endif

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyCoefs")
        
        
    end subroutine ModifyCoefs
    
    !-----------------------------------------------------------------------------
    
    subroutine ModifyAdvectionDiffusionCoefs(PropertyX)

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                       :: PropertyX
        real                                             :: ImpExp_AdvXX, ImpExp_AdvYY           
        real                                             :: AdvectionH_Imp_Exp                  
        !begin--------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyAdvectionDiffusionCoefs")
       
       !check spatial and temporal options
        call CheckTransportOptions(PropertyX, AdvectionH_Imp_Exp, ImpExp_AdvXX, ImpExp_AdvYY)
        
        !Horizontal fluxes 
            
        !Routines from ModuleAdvectionDiffusion
        call HorizontalDiffusion(PropertyX)                              !always explicit
        call HorizontalAdvection(PropertyX, ImpExp_AdvXX, ImpExp_AdvYY)  !explicit or implicit


cd5 :   if (Me%Output%Boxes_ON) then

            if (ImpExp_AdvXX == ImplicitScheme)    &
                call CalcHorizontalAdvFluxXX(PropertyX, ImpExp_AdvXX)

            if (ImpExp_AdvYY == ImplicitScheme)    &
                call CalcHorizontalAdvFluxYY(PropertyX, ImpExp_AdvYY)

            call Output_Boxes_Fluxes (PropertyX)
            
        end if cd5

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyAdvectionDiffusionCoefs")

        
    end subroutine ModifyAdvectionDiffusionCoefs
    
    !-----------------------------------------------------------------------------
    
    subroutine CheckTransportOptions(PropertyX, AdvectionH_Imp_Exp, ImpExp_AdvXX, ImpExp_AdvYY)
        !Local--------------------------------------------------------------------
        type (T_Property), pointer                       :: PropertyX
        real, intent(OUT)                                :: ImpExp_AdvXX, ImpExp_AdvYY           
        real, intent(OUT)                                :: AdvectionH_Imp_Exp          
        !begin--------------------------------------------------------------------

        if (Me%ComputeOptions%AdvDiff_Explicit) then
            
            !if Explicit - all explicit 
            AdvectionH_Imp_Exp = ExplicitScheme
        
        else
            
            AdvectionH_Imp_Exp = ImplicitScheme
            
            !Horizontal diffusion is always explicit
            
        endif
            
        if(AdvectionH_Imp_Exp == ImplicitScheme) then

            if(PropertyX%Evolution%AdvDiff%ImplicitH_Direction == DirectionX)then
                                           
                !Direction X implicit
                ImpExp_AdvXX = ImplicitScheme
                ImpExp_AdvYY = ExplicitScheme
                
                Me%dj = 1
                Me%di = 0 

                PropertyX%Evolution%AdvDiff%ImplicitH_Direction = DirectionY

            else 
            
                !Direction Y implicit
                ImpExp_AdvXX = ExplicitScheme 
                ImpExp_AdvYY = ImplicitScheme 

                Me%dj = 0
                Me%di = 1 
                
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
        
        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "HorizontalDiffusion")

        
        call HorizontalDiffusionXX(CurrProp)
        
        call HorizontalDiffusionYY(CurrProp)
                
        
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "HorizontalDiffusion")


    end subroutine HorizontalDiffusion

    !--------------------------------------------------------------------------

    subroutine HorizontalDiffusionXX(CurrProp) !Routine from ModuleAdvectionDiffusion

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp
        !Local-----------------------------------------------------------------
        real(8)                                     :: DTPropDouble 
        real(8)                                     :: AuxJ, AreaU
        integer                                     :: i, j
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "HorizontalDiffusionXX")


        DTPropDouble = dble(Me%ExtVar%DT) 

        CHUNK = ChunkI !CHUNK_I(Me%WorkSize%ILB,Me%WorkSize%IUB)
        !$OMP PARALLEL PRIVATE(i,j,AuxJ,AreaU) 
 
do2 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            !Only if there will be water at the end of the time step it can have diffusion
            if ((Me%WaterVolume(i, j) .gt. AlmostZero) .and. (Me%WaterVolume(i, j-1) .gt. AlmostZero)) then
                
                AreaU = (0.5 * (Me%ExtVar%WaterColumnOld(i,j) + Me%ExtVar%WaterColumnOld(i,j-1))) * Me%ExtVar%DYY(i,j  )
                
                AuxJ = CurrProp%ViscosityU      (i,j)                                  &
                       * AreaU                                                         &
                       / Me%ExtVar%DZX          (i,j-1  )                    

                Me%TICOEF3(i,j-1) = Me%TICOEF3(i,j-1) + AuxJ * DTPropDouble /          &
                                      Me%WaterVolume(i, j-1) *                         &
                                     (CurrProp%Concentration(i,j) - CurrProp%Concentration(i,j-1))


                Me%TICOEF3(i,j  ) = Me%TICOEF3(i,j  ) - AuxJ * DTPropDouble /       &
                                      Me%WaterVolume(i, j  ) *                       &
                                     (CurrProp%Concentration(i,j) - CurrProp%Concentration(i,j-1))


            endif

        end do do1
        !$OMP END DO
        end do do2
            
        !$OMP END PARALLEL

        if (Me%Output%Boxes_ON) call CalcHorizontalDifFluxXX(CurrProp)

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "HorizontalDiffusionXX")
    
    end subroutine HorizontalDiffusionXX

    !--------------------------------------------------------------------------

    subroutine HorizontalDiffusionYY(CurrProp)  !Routine from ModuleAdvectionDiffusion

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp
        !Local-----------------------------------------------------------------
        real(8)                                     :: DTPropDouble 
        real(8)                                     :: AuxI, AreaV
        integer                                     :: i, j
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "HorizontalDiffusionYY")


        DTPropDouble = dble(Me%ExtVar%DT) 

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(i,j,AuxI,AreaV) 

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            !Only if there will be water at the end of the time step it can have diffusion
            if ((Me%WaterVolume(i, j) .gt. AlmostZero) .and. (Me%WaterVolume(i-1, j) .gt. AlmostZero)) then

                AreaV = (0.5 * (Me%ExtVar%WaterColumnOld(i,j) + Me%ExtVar%WaterColumnOld(i-1,j))) * Me%ExtVar%DXX(i,j  )
                       
                AuxI = CurrProp%ViscosityV      (i  ,j)                               &
                       * AreaV                                                        &
                       / Me%ExtVar%DZY          (i-1,j  )                    

                Me%TICOEF3(i-1,j) = Me%TICOEF3(i-1,j) + AuxI * DTPropDouble /       &
                                      Me%WaterVolume(i-1, j) *                      &
                                     (CurrProp%Concentration(i,j) - CurrProp%Concentration(i-1,j))


                Me%TICOEF3(i,j  ) = Me%TICOEF3(i,j  ) - AuxI * DTPropDouble /       &
                                      Me%WaterVolume(i  , j) *                      &
                                     (CurrProp%Concentration(i,j) - CurrProp%Concentration(i-1,j))
            endif
        end do do1
        end do do2
        !$OMP END DO

        !$OMP END PARALLEL

        if (Me%Output%Boxes_ON) call CalcHorizontalDifFluxYY(CurrProp)

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "HorizontalDiffusionYY")


    end subroutine HorizontalDiffusionYY

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvection(CurrProp, ImpExp_AdvXX, ImpExp_AdvYY)  !Routine adapted from ModuleAdvectionDiffusion

        !Arguments-------------------------------------------------------------
        real                                :: ImpExp_AdvXX, ImpExp_AdvYY
        type (T_Property), pointer          :: CurrProp

        !Local-----------------------------------------------------------------
!        integer                             :: di,    dj    
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: ILBWS, IUBWS, JLBWS, JUBWS

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "HorizontalAdvection")

        ILBWS = Me%WorkSize%ILB
        IUBWS = Me%WorkSize%IUB

        JLBWS = Me%WorkSize%JLB
        JUBWS = Me%WorkSize%JUB

        ILB   = Me%Size%ILB
        IUB   = Me%Size%IUB

        JLB   = Me%Size%JLB 
        JUB   = Me%Size%JUB

        
        call HorizontalAdvectionXX(CurrProp, ImpExp_AdvXX)

        call HorizontalAdvectionYY(CurrProp, ImpExp_AdvYY)

        
        !The implicit method only needs a thomas array at this point in 3D cases
        !in the routine to update properties the thomas will be called if implicit and will go in implicit direction assigned
        
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "HorizontalAdvection")


    end subroutine HorizontalAdvection

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvectionXX(CurrProp, ImpExp_AdvXX)    !Routine adapted from ModuleAdvectionDiffusion

        !Arguments--------------------------------------------------------------
        type (T_Property), pointer          :: CurrProp
        real                                :: ImpExp_AdvXX

        !Local-----------------------------------------------------------------               

        real(8) :: AdvFluxX, DT1, DT2

        integer :: i,     j                             
        integer :: ILB, IUB, JLB, JUB
        integer :: CHUNK
        !----------------------------------------------------------------------


        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "HorizontalAdvectionXX")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !CHUNK = ChunkI !CHUNK_I(ILB, IUB)
        !! $OMP PARALLEL PRIVATE(i,j,AdvFluxX,DT2,DT1)

        if (.not. Me%NewFormulation) then

            !! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
i1:         do i = ILB, IUB


                call ComputeAdvection1D_V2(JLB+1, JUB+1, Me%ExtVar%DT,                  &
                                        Me%ExtVar%DUX               (i,:),              &
                                        CurrProp%Concentration      (i,:),              &
                                        Me%ExtVar%FluxU             (i,:),              &
                                        Me%WaterVolume              (i,:),              & 
                                        Me%DummyOpenPoints          (i,:),              &
                                        Me%COEF3_HorAdvXX%C_flux    (i,:),              &
                                        Me%COEF3_HorAdvXX%D_flux    (i,:),              &
                                        Me%COEF3_HorAdvXX%E_flux    (i,:),              &
                                        Me%COEF3_HorAdvXX%F_flux    (i,:),              &
                                        CurrProp%Evolution%AdvDiff%AdvMethodH,          &
                                        CurrProp%Evolution%AdvDiff%TVDLimitationH,      &
                                        CurrProp%Evolution%AdvDiff%VolumeRelMax,        &
                                        CurrProp%Evolution%AdvDiff%Upwind2H)

            end do i1
            !! $OMP END DO
            
            !! $OMP END PARALLEL

        endif

cd6:    if (ImpExp_AdvXX == ExplicitScheme)  then !ExplicitScheme = 0

            CHUNK = ChunkI !CHUNK_I(ILB, IUB)
            !$OMP PARALLEL PRIVATE(i,j,AdvFluxX)

doj3 :      do j = JLB, JUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doi3 :      do i = ILB, IUB
            
            !computation needs volumes
            if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i, j-1) .gt. AlmostZero) then

                AdvFluxX =    (Me%COEF3_HorAdvXX%C_flux(i,   j)                          &
                            *  CurrProp%Concentration  (i, j-2)                          &
                            +  Me%COEF3_HorAdvXX%D_flux(i,   j)                          &
                            *  CurrProp%Concentration  (i, j-1)                          &
                            +  Me%COEF3_HorAdvXX%E_flux(i,   j)                          &
                            *  CurrProp%Concentration  (i,   j)                          &
                            +  Me%COEF3_HorAdvXX%F_flux(i,   j)                          &
                            *  CurrProp%Concentration  (i, j+1))

                Me%TICOEF3(i,j  ) = Me%TICOEF3(i,j  ) + AdvFluxX * Me%ExtVar%DT / Me%WaterVolume(i,j  )
                Me%TICOEF3(i,j-1) = Me%TICOEF3(i,j-1) - AdvFluxX * Me%ExtVar%DT / Me%WaterVolume(i,j-1)

            endif

            end do doi3
            !$OMP END DO
            end do doj3
            !$OMP END PARALLEL

        else if (ImpExp_AdvXX == ImplicitScheme) then cd6 !ImplicitScheme = 1

            CHUNK = ChunkI !CHUNK_I(ILB, IUB)
            !$OMP PARALLEL PRIVATE(i,j, DT2, DT1)

doj4 :      do j = JLB, JUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doi4 :      do i = ILB, IUB

            !computation needs volumes
            if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i, j-1) .gt. AlmostZero) then

                DT2 = Me%ExtVar%DT / Me%WaterVolume(i,j  )
                DT1 = Me%ExtVar%DT / Me%WaterVolume(i,j-1)

                Me%COEF3%D(i,j  ) = Me%COEF3%D(i,j  ) - Me%COEF3_HorAdvXX%D_flux(i,   j) * DT2
                Me%COEF3%E(i,j  ) = Me%COEF3%E(i,j  ) - Me%COEF3_HorAdvXX%E_flux(i,   j) * DT2

                Me%COEF3%E(i,j-1) = Me%COEF3%E(i,j-1) + Me%COEF3_HorAdvXX%D_flux(i,   j) * DT1
                Me%COEF3%F(i,j-1) = Me%COEF3%F(i,j-1) + Me%COEF3_HorAdvXX%E_flux(i,   j) * DT1


            endif


            end do doi4
            !$OMP END DO
            end do doj4
            !$OMP END PARALLEL

        else cd6

            stop 'sub. ModulePorousMediaProperties - HorizontalAdvectionXX - ERR01'
        
        endif cd6
        
        
        if (Me%Output%Boxes_ON .and. ImpExp_AdvXX == ExplicitScheme) call CalcHorizontalAdvFluxXX(CurrProp, ImpExp_AdvXX)       


        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "HorizontalAdvectionXX")


    end subroutine HorizontalAdvectionXX

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvectionYY(CurrProp, ImpExp_AdvYY)    !Routine adapted from ModuleAdvectionDiffusion

        !Arguments--------------------------------------------------------------
        real                                :: ImpExp_AdvYY
        type (T_Property), pointer          :: CurrProp

        !Local-----------------------------------------------------------------               
        real(8)                             :: AdvFluxY, DT1, DT2
        integer                             :: i,     j                       
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "HorizontalAdvectionYY")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        
        if (.not. Me%NewFormulation) then

            !CHUNK = ChunkJ !CHUNK_J(JLB, JUB)
          
            !! $OMP PARALLEL PRIVATE(i,j,AdvFluxY,DT2,DT1)       
             !! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
j1:         do j = JLB, JUB


                call ComputeAdvection1D_V2(ILB+1, IUB+1, Me%ExtVar%DT,                              &
                                        Me%ExtVar%DVY               (:,j),                          &
                                        CurrProp%Concentration      (:,j),                          &
                                        Me%ExtVar%FluxV             (:,j),                          &
                                        Me%WaterVolume              (:,j),                          & 
                                        Me%DummyOpenPoints          (:,j),                          &
                                        Me%COEF3_HorAdvYY%C_flux    (:,j),                          &
                                        Me%COEF3_HorAdvYY%D_flux    (:,j),                          &
                                        Me%COEF3_HorAdvYY%E_flux    (:,j),                          &
                                        Me%COEF3_HorAdvYY%F_flux    (:,j),                          &
                                        CurrProp%Evolution%AdvDiff%AdvMethodH,                      &
                                        CurrProp%Evolution%AdvDiff%TVDLimitationH,                  &
                                        CurrProp%Evolution%AdvDiff%VolumeRelMax,                    &
                                        CurrProp%Evolution%AdvDiff%Upwind2H)

            end do j1
            !! $OMP END DO
            
            !! $OMP END DO NOWAIT
            !! $OMP END PARALLEL
        endif

cd6:    if (ImpExp_AdvYY == ExplicitScheme)  then !ExplicitScheme = 0
            
            CHUNK = ChunkJ !CHUNK_J(JLB, JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,AdvFluxY)

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj3 :      do j = JLB, JUB
doi3 :      do i = ILB, IUB

            !computation needs volumes
            if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i-1, j) .gt. AlmostZero) then

                AdvFluxY =    (Me%COEF3_HorAdvYY%C_flux(  i, j)                          &
                            *  CurrProp%Concentration  (i-2, j)                          &
                            +  Me%COEF3_HorAdvYY%D_flux(  i, j)                          &
                            *  CurrProp%Concentration  (i-1, j)                          &
                            +  Me%COEF3_HorAdvYY%E_flux(  i, j)                          &
                            *  CurrProp%Concentration  (  i, j)                          &
                            +  Me%COEF3_HorAdvYY%F_flux(  i, j)                          &
                            *  CurrProp%Concentration  (i+1, j))

                Me%TICOEF3(i  ,j) = Me%TICOEF3(i  ,j) + AdvFluxY * Me%ExtVar%DT / Me%WaterVolume(i  ,j)
                Me%TICOEF3(i-1,j) = Me%TICOEF3(i-1,j) - AdvFluxY * Me%ExtVar%DT / Me%WaterVolume(i-1,j)

            endif

            end do doi3
            end do doj3
            !$OMP END DO
            !$OMP END PARALLEL

        else if (ImpExp_AdvYY == ImplicitScheme) then cd6 !ImplicitScheme = 1

            CHUNK = ChunkJ !CHUNK_J(JLB, JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,DT2,DT1)

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj4 :      do j = JLB, JUB
doi4 :      do i = ILB, IUB


            !computation needs volumes
            if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i-1, j) .gt. AlmostZero) then

                DT2 = Me%ExtVar%DT / Me%WaterVolume(i  ,j  )
                DT1 = Me%ExtVar%DT / Me%WaterVolume(i-1,j  )

                Me%COEF3%D(i,j  ) = Me%COEF3%D(i,j  ) - Me%COEF3_HorAdvYY%D_flux(i,   j) * DT2
                Me%COEF3%E(i,j  ) = Me%COEF3%E(i,j  ) - Me%COEF3_HorAdvYY%E_flux(i,   j) * DT2

                Me%COEF3%E(i-1,j) = Me%COEF3%E(i-1,j) + Me%COEF3_HorAdvYY%D_flux(i,   j) * DT1
                Me%COEF3%F(i-1,j) = Me%COEF3%F(i-1,j) + Me%COEF3_HorAdvYY%E_flux(i,   j) * DT1

            endif


            end do doi4
            end do doj4
            !$OMP END DO
            !$OMP END PARALLEL
            
        else cd6

            stop 'sub. HorizontalAdvectionYY - ModuleAdvectionDiffusion - ERR01'
        
        endif cd6

        if (Me%Output%Boxes_ON .and. ImpExp_AdvYY == ExplicitScheme) call CalcHorizontalAdvFluxXX(CurrProp, ImpExp_AdvYY)

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "HorizontalAdvectionYY")

    end subroutine HorizontalAdvectionYY
    
    !--------------------------------------------------------------------------

    subroutine CalcHorizontalAdvFluxXX(CurrProp, Weigth)

        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp    
        real, intent(IN) :: Weigth !Refers to the wigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------

        integer :: i,     j
        integer :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "CalcHorizontalAdvFluxXX")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(i,j)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i, j - 1) .gt. AlmostZero) then

            Me%Fluxes%AdvFluxX(i, j) =                                               &
                          Me%Fluxes%AdvFluxX(i, j)                                   &
                        + Weigth                                                     &
                        * (Me%COEF3_HorAdvXX%C_flux(i,   j)                          &
                        *  CurrProp%Concentration  (i, j-2)                          &
                        +  Me%COEF3_HorAdvXX%D_flux(i,   j)                          &
                        *  CurrProp%Concentration  (i, j-1)                          &
                        +  Me%COEF3_HorAdvXX%E_flux(i,   j)                          &
                        *  CurrProp%Concentration  (i,   j)                          &
                        +  Me%COEF3_HorAdvXX%F_flux(i,   j)                          &
                        *  CurrProp%Concentration  (i, j+1))

        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "CalcHorizontalAdvFluxXX")

        !----------------------------------------------------------------------

    end subroutine CalcHorizontalAdvFluxXX

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalAdvFluxYY(CurrProp, Weigth)

        !External--------------------------------------------------------------
    
        real, intent(IN) :: Weigth !Refers to the wigth of Implicit-Explicit calculations
        type (T_Property), pointer                  :: CurrProp    

        !Local-----------------------------------------------------------------

        integer :: i,     j
        integer :: CHUNK
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "CalcHorizontalAdvFluxYY")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j)
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i-1, j) .gt. AlmostZero) then
            Me%Fluxes%AdvFluxY(i, j) =                                                  &
                          Me%Fluxes%AdvFluxY(i, j)                                      &
                        + Weigth                                                        &
                        * (Me%COEF3_HorAdvYY%C_flux(i  , j)                          &
                        *  CurrProp%Concentration  (i-2, j)                          &
                        +  Me%COEF3_HorAdvYY%D_flux(i,   j)                          &
                        *  CurrProp%Concentration  (i-1, j)                          &
                        +  Me%COEF3_HorAdvYY%E_flux(i,   j)                          &
                        *  CurrProp%Concentration  (i,   j)                          &
                        +  Me%COEF3_HorAdvYY%F_flux(i,   j)                          &
                        *  CurrProp%Concentration  (i+1, j))
        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "CalcHorizontalAdvFluxYY")

        !----------------------------------------------------------------------

    end subroutine CalcHorizontalAdvFluxYY

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalDifFluxXX(CurrProp)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp    

        !Local-----------------------------------------------------------------
        integer                                     :: i, j  
        integer                                     :: CHUNK
        real                                        :: AreaU
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "CalcHorizontalDifFluxXX")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(i,j,AreaU)
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i, j-1) .gt. AlmostZero) then                           
                
                AreaU = (0.5 * (Me%ExtVar%WaterColumnOld(i,j) + Me%ExtVar%WaterColumnOld(i,j-1))) * Me%ExtVar%DYY(i,j)
                
                Me%Fluxes%DifFluxX(i, j) =                                           &
                              Me%Fluxes%DifFluxX      (i,  j)                        &
                            - CurrProp%ViscosityU     (i,  j)                        &
                            * AreaU                                                  &
                            / Me%ExtVar%DZX           (i,j-1)                        &
                            *(CurrProp%Concentration  (i,  j)                        &
                            - CurrProp%Concentration  (i,j-1))

            endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "CalcHorizontalDifFluxXX")


    end subroutine CalcHorizontalDifFluxXX

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalDifFluxYY(CurrProp)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProp    

        !Local-----------------------------------------------------------------
        integer                                     :: i, j 
        integer                                     :: CHUNK
        real                                        :: AreaV

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "CalcHorizontalDifFluxYY")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,AreaV)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%WaterVolume(i, j) .gt. AlmostZero .and. Me%WaterVolume(i-1, j) .gt. AlmostZero) then
            
            AreaV = (0.5 * (Me%ExtVar%WaterColumnOld(i,j) + Me%ExtVar%WaterColumnOld(i-1,j))) * Me%ExtVar%DXX(i,j  )
        
            Me%Fluxes%DifFluxY(i, j) =                                               &
                          Me%Fluxes%DifFluxY      (i  , j)                           &
                        - CurrProp%ViscosityV     (i  , j)                           &
                        * AreaV                                                      &
                        / Me%ExtVar%DZY           (i-1, j)                           &
                        *(CurrProp%Concentration  (i  , j)                           &
                        - CurrProp%Concentration  (i-1, j))
        endif
        end do doi1
        end do doj1
        !$OMP END DO
        
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "CalcHorizontalDifFluxYY")
        

    end subroutine CalcHorizontalDifFluxYY

    !--------------------------------------------------------------------------

!    subroutine ModifyDrainageNetworkCoefs (PropertyX)
!    
!        !Arguments-------------------------------------------------------------
!        type (T_Property), pointer                  :: PropertyX
!
!
!        !Local-----------------------------------------------------------------
!        type (T_Property), pointer                  :: CurrProperty
!        integer                                     :: i, j, CHUNK
!        real(8)                                     :: aux 
!        !Begin-----------------------------------------------------------------
!
!        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyDrainageNetworkCoefs")
!   
!        
!        CurrProperty => PropertyX
!        
!       !Flux between river and runoff in layers
!       
!        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
!      
!        !$OMP PARALLEL PRIVATE(i,j,aux)
!        
!        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!            
!            if ((Me%ExtVar%BasinPoints(I,J) == BasinPoint) .and.  &
!                (Me%ExtVar%RiverPoints(I,J) == BasinPoint) .and.  &
!                (Me%DummyOpenPoints(I,J)    == BasinPoint)) then   
!                       
!                !Auxuliar value for transport - units of flow^-1
!                !s/m3
!                aux             = (Me%ExtVar%DT/Me%WaterVolume(i,j) )
!                
!                ! Positive flow -> looses mass
!                Me%COEFExpl%CoefInterfDN(i,j) = - aux * Me%ExtVar%FlowToChannels(i,j)
!
!              
!                ! mass going to channel -> conc from runoff
!                if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
!                    
!                    CurrProperty%ConcInInterfaceDN(i,j) =  CurrProperty%ConcentrationOld(i,j)
!                    
!               
!                !mass coming from channel -> conc from DN
!                elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
!                    
!                    CurrProperty%ConcInInterfaceDN(i,j) = CurrProperty%ConcentrationDN(i,j)
!                    
!                endif
!                
!           endif
!        
!        enddo
!        enddo
!        !$OMP END DO
!        
!        !$OMP END PARALLEL
!                           
!        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyDrainageNetworkCoefs")
!
!   
!   end subroutine ModifyDrainageNetworkCoefs

    !---------------------------------------------------------------------------

    !DN fluxes had to be separated from Advection diffusion since in one time step
    !runoff could run out of water and the mixing between runoff cels and link to DN had to be
    !separated.
    !Concentration only changes if flow from river .or. (flow to river .and. no water column) but
    !was left the two cases (from and to river) because of mass balance accounting (FlowMass)
    subroutine ModifyDrainageNetworkInterface
    
        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, CHUNK !, STAT_CALL
        real(8)                                     :: Prop, WaterVolumeOld, WaterVolumeNew
        real(8)                                     :: FlowMass !, WaterColumnNew
        !real(8), dimension(:,:), pointer            :: WaterColumnFinal
        character (Len = 5)                         :: str_i, str_j
        character (Len = 15)                        :: str_conc
        character (len = StringLength)              :: string_to_be_written        
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyDrainageNetworkInterface")
   
!        call GetRunoffWaterColumn (Me%ObjRunoff, WaterColumnFinal, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ModifyDrainageNetworkInterface - ModuleRunoffProperties - ERR10'

        CurrProperty => Me%FirstProperty
        
        do while (associated(CurrProperty))

            if (CurrProperty%Evolution%AdvectionDiffusion) then        
               !Flux between river and runoff 
               
                CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
              
                !!$OMP PARALLEL PRIVATE(i,j,Prop, WaterVolumeOld, WaterVolumeNew, FlowMass)
                
                !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                    if ((Me%ExtVar%BasinPoints(I,J) == BasinPoint) .and.  &
                        (Me%ExtVar%RiverPoints(I,J) == BasinPoint)) then   
                        
                        ! mass going to channel -> conc from runoff - concentration does not change
                        if (Me%ExtVar%FlowToChannels(i,j) .ge. 0.0) then
                            
                            !this is conc after transport, the same as used in Drainage network
                            Prop =  CurrProperty%Concentration(i,j)

                            if (CurrProperty%Concentration(i,j) < 0.0) then
                                write(str_i, '(i3)') i 
                                write(str_j, '(i3)') j
                                write(str_conc, '(ES10.3)') CurrProperty%Concentration(i,j)                               
                                string_to_be_written = trim(CurrProperty%ID%Name) //' in RP is < 0 going from runoff to river' &
                                                       //' in cell (i,j) with conc: '   &
                                                      //str_i//','//str_j//' '//str_conc
                                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)
                            endif     
                       
                        !mass coming from channel -> conc from DN
                        elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
                            
                            Prop = CurrProperty%ConcentrationDN(i,j)

                            if (CurrProperty%ConcentrationDN(i,j) < 0.0) then
                                write(str_i, '(i3)') i 
                                write(str_j, '(i3)') j
                                write(str_conc, '(ES10.3)') CurrProperty%ConcentrationDN(i,j)                               
                                string_to_be_written = trim(CurrProperty%ID%Name) //' in RP is < 0 going from river to runoff' &
                                                       //' in cell (i,j) with conc: '   &
                                                      //str_i//','//str_j//' '//str_conc
                                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

                            endif
                            
                        endif

                        !Prop and WaterVolumeOld are after transport and before this update
                        WaterVolumeOld = Me%ExtVar%WaterColumnAT(i,j) * Me%ExtVar%Area(i,j)
                        !New watervolume can not be the final because it may exist boundary flux (See BoundaryInterface).
                        !So new is old (after transport) less flux. Positive flux removes WC
                        !WaterVolumeNew = WaterColumnFinal(i,j) * Me%ExtVar%Area(i,j)
                        WaterVolumeNew = WaterVolumeOld - Me%ExtVar%FlowToChannels(i,j) * Me%ExtVar%DT
                        !WaterColumnNew = WaterVolumeNew / Me%ExtVar%Area(i,j)
                        
                        !g = m3/s * s * g/m3
                        FlowMass       = Me%ExtVar%FlowToChannels(i,j) * Me%ExtVar%DT * Prop                    
                        
                        !if (WaterColumnFinal(i,j) .gt. AlmostZero) then
                        if (WaterVolumeNew .gt. AlmostZero) then
                            !Update New Concentration
                            !g/m3 = ((g/m3 * m3) + g)/ m3
                            CurrProperty%Concentration(i,j) = ((CurrProperty%Concentration(i,j) * WaterVolumeOld)  &
                                                                - FlowMass) / WaterVolumeNew
                        else
                            CurrProperty%Concentration(i,j) = 0.0                        
                        endif
                        
                        if (Me%CheckGlobalMass) then
                            !kg = g * 1e-3 kg/g
                            CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass + (FlowMass * 1e-3)
                        endif
!                        write(*,*) CurrProperty%ID%Name, i,j, Me%ExtVar%FlowToChannels(i,j), &
!                                    Me%ExtVar%WaterColumnOld(i,j),Me%ExtVar%WaterColumnAT(i,j),Me%ExtVar%WaterColumn(i,j) &
!                                  , CurrProperty%ConcentrationOld(i,j),Prop,CurrProperty%Concentration(i,j)
                   endif
                
                enddo
                enddo
                !!$OMP END DO
                
                !!$OMP END PARALLEL

            endif
            
            CurrProperty => CurrProperty%Next
        
        enddo
        
!        call UngetRunoff (Me%ObjRunoff, WaterColumnFinal, STAT_CALL)  
!        if (STAT_CALL /= SUCCESS_)  stop 'ModifyDrainageNetworkInterface - ModuleRunoffProperties - ERR020'     
                           
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyDrainageNetworkInterface")

   
   end subroutine ModifyDrainageNetworkInterface

    !---------------------------------------------------------------------------

!    subroutine ModifyBoundaryCoefs(PropertyX)
!
!        !Arguments----------------------------------------------------------------
!        type (T_Property), pointer                  :: PropertyX
!        
!        !Local--------------------------------------------------------------------
!        integer                                     :: i,j, CHUNK !, STAT_CALL
!        real(8)                                     :: aux
!        
!        !Begin--------------------------------------------------------------------
!        
!        
!        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB) 
!      
!        !$OMP PARALLEL PRIVATE(i,j,aux)
!        
!        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!            
!            if ((Me%ExtVar%BoundaryCells(I,J) == WaterPoint) .and. (Me%WaterVolume(i,j) .gt. AlmostZero)) then
!                
!                !Auxuliar value for transport - units of flow^-1
!                !s/m3
!                aux             = (Me%ExtVar%DT/Me%WaterVolume(i,j) )
!                
!                ! Positive flow -> gains mass
!                !Currently BoundaryFlux is always negative (water can only exit) so also the coef
!                Me%COEFExpl%CoefInterfBoundary(i,j) = aux * Me%ExtVar%BoundaryFlux(i,j)
!                
!                !negative flow - exiting - Prop is from runoff if explicit
!                if (Me%ExtVar%BoundaryFlux(i,j) .le. 0.0) then
!                    
!                    PropertyX%ConcInBoundary(i,j) = PropertyX%ConcentrationOld(i,j)
!                
!                else !entering - Prop is imposed
!                    
!                    if (PropertyX%Evolution%Boundary%BoundaryCondition == ImposedValue_) then
!                        
!                        PropertyX%ConcInBoundary(i,j) = PropertyX%Evolution%Boundary%DefaultBoundary
!                    
!                    elseif (PropertyX%Evolution%Boundary%BoundaryCondition == NullGradient_) then
!                        
!                        PropertyX%ConcInBoundary(i,j) = PropertyX%ConcentrationOld(i,j)
!                    
!                    endif
!                endif
!                
!            endif
!            
!        enddo
!        enddo
!        !$OMP END DO
!        !$OMP END PARALLEL
!
!    
!    end subroutine ModifyBoundaryCoefs
!
!    !-----------------------------------------------------------------------------
    !Concentration changes in recieving cell .or. (origin cell .and. it gets no water at the end)
    subroutine ModifyRouteDFourInterface
    
        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, it, jt, CHUNK !, STAT_CALL
        real(8)                                     :: Prop, WaterVolumeOldDFour, WaterVolumeNewDFour
        real(8)                                     :: WaterVolumeOldLowNeighbour, WaterVolumeNewLowNeighbour
        real(8)                                     :: FlowMass !, WaterColumnNew
        !real(8), dimension(:,:), pointer            :: WaterColumnFinal      
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyRouteDFourInterface")
        
        CurrProperty => Me%FirstProperty
        
        do while (associated(CurrProperty))

            if (CurrProperty%Evolution%AdvectionDiffusion) then        
               !Flux in routeD4 points 
               
                CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
              
                !!$OMP PARALLEL PRIVATE(i,j,Prop, WaterVolumeOld, WaterVolumeNew, FlowMass)
                
                !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                    if (Me%ExtVar%RouteDFourCells(I,J) == BasinPoint .and. Me%ExtVar%RouteDFourLowerI(i,j) /= null_int) then  

                        it = Me%ExtVar%RouteDFourLowerI(i,j)
                        jt = Me%ExtVar%RouteDFourLowerJ(i,j)
                        
                        ! mass from cell i,j to diagonal it, jt
                        if (Me%ExtVar%RouteDFourFlux(i,j) .ge. 0.0) then
                            
                            Prop =  CurrProperty%Concentration(i,j)                            
                       
                        !mass coming from diagonal
                        elseif (Me%ExtVar%RouteDFourFlux(i,j) .lt. 0.0) then

                            Prop = CurrProperty%Concentration(it,jt)                            

                        endif

                        !Change conc in recieving cell. origin cell does not change conc.
                        !WaterVolumeOld is after transport and after DN update in DrainageNetworkInterface
                        WaterVolumeOldDFour        = Me%ExtVar%WaterColumnAT(i,j) * Me%ExtVar%Area(i,j)
                        WaterVolumeOldLowNeighbour = Me%ExtVar%WaterColumnAT(it,jt) * Me%ExtVar%Area(it,jt)
                        if (Me%ExtVar%CoupledDN) then
                            WaterVolumeOldDFour        = WaterVolumeOldDFour -                             &
                                                         Me%ExtVar%FlowToChannels(i,j) * Me%ExtVar%DT
                            WaterVolumeOldLowNeighbour = WaterVolumeOldLowNeighbour -                      &
                                                          Me%ExtVar%FlowToChannels(it,jt) * Me%ExtVar%DT
                        endif
                        
                        !New watervolume after routedfour flux
                        !So new is old less flux. Positive flux adds WC
                        WaterVolumeNewDFour        = WaterVolumeOldDFour                                         &
                                                     - Me%ExtVar%RouteDFourFlux(i,j) * Me%ExtVar%DT
                        WaterVolumeNewLowNeighbour = WaterVolumeOldLowNeighbour                                  &
                                                     + Me%ExtVar%RouteDFourFlux(i,j) * Me%ExtVar%DT
                        
                        !g = m3/s * s * g/m3
                        FlowMass       = Me%ExtVar%RouteDFourFlux(i,j) * Me%ExtVar%DT * Prop                    

                        if (WaterVolumeNewDFour .gt. AlmostZero) then
                            !Update New Concentration
                            !g/m3 = ((g/m3 * m3) + g)/ m3
                            CurrProperty%Concentration(i,j) = ((CurrProperty%Concentration(i,j)    &
                                                                 * WaterVolumeOldDFour)            &
                                                                   - FlowMass) / WaterVolumeNewDFour
                        else
                            CurrProperty%Concentration(i,j) = 0.0                        
                        endif                        
                        if (WaterVolumeNewLowNeighbour .gt. AlmostZero) then
                            !Update New Concentration
                            !g/m3 = ((g/m3 * m3) + g)/ m3
                            CurrProperty%Concentration(it,jt) = ((CurrProperty%Concentration(it,jt) &
                                                                   * WaterVolumeOldLowNeighbour)    &
                                                                   + FlowMass) / WaterVolumeNewLowNeighbour
                        else
                            CurrProperty%Concentration(it,jt) = 0.0                        
                        endif
                        
                   endif
                
                enddo
                enddo
                !!$OMP END DO
                
                !!$OMP END PARALLEL

            endif
            
            CurrProperty => CurrProperty%Next
        
        enddo
        
                           
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyRouteDFourInterface")

   
   end subroutine ModifyRouteDFourInterface

    !---------------------------------------------------------------------------

    !If river cells exist that are boundary than river flow has to be added/removed from
    !volume after transport to get water volume old
    subroutine ModifyBoundaryInterface
    
        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, CHUNK, STAT_CALL
        real(8)                                     :: Prop, WaterVolumeOld, WaterVolumeNew
        real(8)                                     :: FlowMass
        real(8), dimension(:,:), pointer            :: WaterColumnFinal
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyBoundaryInterface")
   
        call GetRunoffWaterColumn (Me%ObjRunoff, WaterColumnFinal, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyBoundaryInterface - ModuleRunoffProperties - ERR10'

        CurrProperty => Me%FirstProperty
        
        do while (associated(CurrProperty))

            if (CurrProperty%Evolution%AdvectionDiffusion) then        
               !Flux between runoff and outside 
               
                CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
              
                !!$OMP PARALLEL PRIVATE(i,j,Prop, WaterVolumeOld, WaterVolumeNew, FlowMass)
                
                !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                    if ((Me%ExtVar%BoundaryCells(I,J) == BasinPoint)) then   
                        
                        ! mass going to outside -> conc from runoff - concentration does not change
                        if (Me%ExtVar%BoundaryFlux(i,j) .le. 0.0) then
                            
                            !this is conc after transport, just before going to outside
                            Prop =  CurrProperty%Concentration(i,j)
                       
                        !mass coming from outsie 
                        elseif (Me%ExtVar%BoundaryFlux(i,j) .gt. 0.0) then
                            
                            if (CurrProperty%Evolution%Boundary%BoundaryCondition == ImposedValue_) then
                                
                                Prop = CurrProperty%Evolution%Boundary%DefaultBoundary
                            
                            elseif (CurrProperty%Evolution%Boundary%BoundaryCondition == NullGradient_) then
                                
                                Prop = CurrProperty%Concentration(i,j)
                            
                            endif                            
                        endif

                        !Prop and WaterVolumeOld are after transport and before this update
                        !Old watervolume can not be the after transport because it may exist river flux (See DN Interface).
                        !So old is new (final) less flux. Negative flux removes WC                         
                        !WaterVolumeOld = Me%ExtVar%WaterColumnAT(i,j) * Me%ExtVar%Area(i,j)
                        WaterVolumeNew = WaterColumnFinal(i,j) * Me%ExtVar%Area(i,j)
                        WaterVolumeOld = WaterVolumeNew - Me%ExtVar%BoundaryFlux(i,j) * Me%ExtVar%DT
                        
                        !g = m3/s * s * g/m3
                        FlowMass       = Me%ExtVar%BoundaryFlux(i,j) * Me%ExtVar%DT * Prop                    
                        
                        if (WaterVolumeNew .gt. AlmostZero) then
                            !Update New Concentration
                            !g/m3 = ((g/m3 * m3) + g)/ m3
                            CurrProperty%Concentration(i,j) = ((CurrProperty%Concentration(i,j) * WaterVolumeOld)  &
                                                                 + FlowMass) / WaterVolumeNew
                        else
                            CurrProperty%Concentration(i,j) = 0.0
                        endif

                   endif
                
                enddo
                enddo
                !!$OMP END DO
                
                !!$OMP END PARALLEL

            endif
            
            CurrProperty => CurrProperty%Next
        
        enddo
        
        call UngetRunoff (Me%ObjRunoff, WaterColumnFinal, STAT_CALL)  
        if (STAT_CALL /= SUCCESS_)  stop 'ModifyBoundaryInterface - ModuleRunoffProperties - ERR020'     
                           
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyBoundaryInterface")

   
   end subroutine ModifyBoundaryInterface

    !---------------------------------------------------------------------------

    subroutine ModifyPropertyValues(PropertyX)
        
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, CHUNK, di, dj
        integer                                     :: IJmin, IJmax, JImin, JImax
        real(8)                                     :: coefB, CoefInterfDN
        real(8)                                     :: CoefInterfBoundary        
        !Begin-----------------------------------------------------------------    

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyPropertyValues")

        
        CurrProperty => PropertyX
        
        if (Me%ComputeOptions%AdvDiff_Explicit) then

            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,CoefInterfDN,CoefInterfBoundary, coefB)
               
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)              
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then   
                    
                    !evaluate if there is water
                    if (Me%WaterVolume(i,j) .gt. AlmostZero) then
                        CoefInterfDN       = Me%COEFExpl%CoefInterfDN(i,j)
                        CoefInterfBoundary = Me%COEFExpl%CoefInterfBoundary(i,j)
                       
                        !CoefB = Me%ExtVar%WaterColumnOld(i,j)/Me%ExtVar%WaterColumn(i,j)
                        !CoefB = Me%ExtVar%WaterColumnOld(i,j)/Me%ExtVar%WaterColumnAT(i,j)
                        CoefB = Me%WaterColumnBT(i,j)/Me%ExtVar%WaterColumnAT(i,j)
                        
                        Me%TICOEF3(i,j  ) = Me%TICOEF3(i,j) + CoefB * CurrProperty%ConcentrationOld(i,j  )                  &
                                                            + CoefInterfDN     * CurrProperty%ConcInInterfaceDN(i,j)        &
                                                            + CoefInterfBoundary * CurrProperty%ConcInBoundary(i,j)
                    else
                        Me%TICOEF3(i,j  ) = 0.0
                    endif                                                 
                    
                    CurrProperty%Concentration(i,j) = Me%TICOEF3(i,j  )

                    if (CurrProperty%Evolution%BottomFluxes) then
                        ![kg/m2] = [g/m3]* [m * m2] * [1E-3kg/g] /[m2] + [kg/m2]
                        CurrProperty%TotalConcentration (i,j) = ((CurrProperty%Concentration (i,j) * 1E-3                  &
!                                                                  * Me%ExtVar%WaterColumn(i,j) * Me%ExtVar%Area(i, j))    &
                                                                  * Me%ExtVar%WaterColumnAT(i,j) * Me%ExtVar%Area(i, j))   &
                                                                 / Me%ExtVar%Area(i, j))                                   &
                                                                 + CurrProperty%BottomConcentration (i,j)   
                    endif      

                endif
           enddo
           enddo
           !$OMP END DO
           !$OMP END PARALLEL  
           
        else !implicit
            
            CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
          
            !$OMP PARALLEL PRIVATE(i,j,CoefInterfDN,CoefInterfBoundary, coefB)
               
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)                          
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then           
                    
                    if (Me%WaterVolume(i,j) .gt. AlmostZero) then
                    
                        CoefInterfDN       = Me%COEFExpl%CoefInterfDN(i,j)
                        CoefInterfBoundary = Me%COEFExpl%CoefInterfBoundary(i,j)
                        
                        !Water exiting runoff to DN - Comput Conc Implicitly
                        if (CoefInterfDN .lt. 0.0) then
                            !Add coef to implicit part
                            Me%COEF3%E(i,j) = Me%COEF3%E(i,j) - CoefInterfDN
                            !Nullify explicit coef that enters in TICoef
                            CoefInterfDN        = 0.0
                        endif
                        
                        !Currently is alwasy negative - water exiting runoff
                        if (CoefInterfBoundary .lt. 0.0) then
                            !Add coef to implicit part
                            Me%COEF3%E(i,j) = Me%COEF3%E(i,j) - CoefInterfBoundary
                            !Nullify explicit coef that enters in TICoef
                            CoefInterfBoundary  = 0.0
                        endif

                        !CoefB = Me%ExtVar%WaterColumnOld(i,j)/Me%ExtVar%WaterColumn(i,j)
                        !CoefB = Me%ExtVar%WaterColumnOld(i,j)/Me%ExtVar%WaterColumnAT(i,j)
                        CoefB = Me%WaterColumnBT(i,j)/Me%ExtVar%WaterColumnAT(i,j)

                        Me%TICOEF3(i,j) = Me%TICOEF3(i,j)                                                      &
                                            + coefB * CurrProperty%ConcentrationOld(i,j)                       &
                                            + CoefInterfDN     * CurrProperty%ConcInInterfaceDN(i,j)           &
                                            + CoefInterfBoundary * CurrProperty%ConcInBoundary(i,j)
                                                      
                    else
                        Me%TICOEF3(i,j) = 0.0
                    endif
                endif
            enddo
            enddo

           !$OMP END DO
           !$OMP END PARALLEL  
            
            dj = Me%dj
            di = Me%di
                        
            IJmin = Me%WorkSize%ILB * dj + Me%WorkSize%JLB * di
            IJmax = Me%WorkSize%IUB * dj + Me%WorkSize%JUB * di

            JImin = Me%WorkSize%ILB * di + Me%WorkSize%JLB * dj
            JImax = Me%WorkSize%IUB * di + Me%WorkSize%JUB * dj  
                     
            call THOMAS_2D(IJmin, IJmax,                                                &
                         JImin, JImax,                                                  &
                         di, dj,                                                        &
                         Me%COEF3%D,                                                    &
                         Me%COEF3%E,                                                    &
                         Me%COEF3%F,                                                    &
                         Me%TICOEF3,                                                    &
                         CurrProperty%Concentration,                                    &
                         Me%VECG,                                                       &
                         Me%VECW)      
            
            !Update total conc
            if (CurrProperty%Evolution%BottomFluxes) then
            
                CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
              
                !$OMP PARALLEL PRIVATE(i,j)
                   
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)                          
                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then           
                        ![kg/m2] = [g/m3]* [m * m2] * [1E-3kg/g] /[m2] + [kg/m2]
                        CurrProperty%TotalConcentration (i,j) = ((CurrProperty%Concentration (i,j) * 1E-3                 &
!                                                                  * Me%ExtVar%WaterColumn(i,j) * Me%ExtVar%Area(i, j))    &
                                                                  * Me%ExtVar%WaterColumnAT(i,j) * Me%ExtVar%Area(i, j))    &
                                                                 / Me%ExtVar%Area(i, j))                                  &
                                                                 + CurrProperty%BottomConcentration (i,j) 
                    endif 
                enddo
                enddo 
               !$OMP END DO
               !$OMP END PARALLEL  
                
            endif               
        endif

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyPropertyValues")

    
    end subroutine ModifyPropertyValues
    
    !---------------------------------------------------------------------------

!    subroutine ModifyInterfaceMassFluxes(PropertyX) 
!
!        !Arguments-------------------------------------------------------------
!        type (T_Property), pointer                  :: PropertyX
!
!        !Local-----------------------------------------------------------------
!        type (T_Property), pointer                  :: CurrProperty
!        integer                                     :: i, j !, STAT_CALL !, CHUNK
!        !Begin-----------------------------------------------------------------    
!        
!        CurrProperty => PropertyX
!        
!        !!Drainage network interface mass balance 
!        if (Me%ExtVar%CoupledDN) then
!           
!            !!! $OMP PARALLEL PRIVATE(I,J,K)
!            !!! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
!            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!                
!                if (Me%ExtVar%BasinPoints(I,J) == BasinPoint .and. Me%ExtVar%RiverPoints(I,J) == BasinPoint) then   
!
!                    ! mass going to channel -> conc from runoff
!                    if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
!                        
!                        !Global Mass Exchange
!                        ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
!                        CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
!                                  + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%ConcentrationOld(i,j)          &
!!                          + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%Concentration(i,j)             &
!                             * 1e-3 * Me%ExtVar%DT)
!                    
!                    !mass coming from channel -> conc from DN
!                    elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
!                        
!                        !Global Mass Exchange
!                        ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
!                        CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
!                          + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%ConcentrationDN(i,j)             &
!                             * 1e-3 * Me%ExtVar%DT)  
!                        
!                    endif
!                endif
!            enddo
!            enddo
!           
!        endif
!        
!    
!    end subroutine ModifyInterfaceMassFluxes
    
    !---------------------------------------------------------------------------

!    subroutine ModifyAdvectionDiffusion_Explicit (PropertyX)
!    
!        !Arguments-------------------------------------------------------------
!        type (T_Property), pointer                  :: PropertyX
!
!
!        !Local-----------------------------------------------------------------
!        type (T_Property), pointer                  :: CurrProperty
!        integer                                     :: i, j!, CHUNK
!        real(8)                                     :: Area_Vertical, Area_Top_U, Area_Top_V
!        real(8)                                     :: Area_Bottom_U, Area_Bottom_V
!!        real(8)                                     :: AdvTermB_Top, AdvTermC
!        real(8)                                     :: AdvTermA_U, AdvTermA_V
!        real(8)                                     :: AdvTermB_Top_U, AdvTermB_Top_V
!        real(8)                                     :: AdvTermB_Bottom_U, AdvTermB_Bottom_V
!        real(8)                                     :: AdvTermC_U, AdvTermC_V
!        real(8)                                     :: DifTerm_Top_U, DifTerm_Top_V
!        real(8)                                     :: DifTerm_Bottom_U, DifTerm_Bottom_V         
!        real(8)                                     :: aux 
!        real(8)                                     :: cofA_U,cofB_U,cofC_U
!        real(8)                                     :: cofA_V,cofB_V,cofC_V, CofB
!        real(8)                                     :: cofInterfaceDN, ConcInInterfaceDN
!        real, pointer, dimension(:,:)               :: FluxU, FluxV
!        real, pointer, dimension(:,:  )             :: DZX, DZY, DXX, DYY
!        real(8), pointer, dimension(:,:)            :: WaterColumn, WaterColumnOld
!        !Begin-----------------------------------------------------------------
!   
!        !!CHUNK = ChunkK !CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
!
!        CurrProperty => PropertyX
!
!        if (Me%ComputeOptions%AdvDiff_DiffMethod == AdvDif_Diff_Jury_) then ! new formulation based on Jury
!            call ModifyDiffusivity_New(CurrProperty)
!        elseif (Me%ComputeOptions%AdvDiff_DiffMethod == AdvDif_Diff_Old_) then !old formulation advection diffusion
!            call ModifyDiffusivity_Old(CurrProperty)
!        else
!            write(*,*)'Diffusion method to be used unrecognized,'
!            write(*,*)'Please check ADVDIFF_DIFF_METHOD keyword'
!            stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR10'                            
!        endif
!        
!        
!        FluxU          => Me%ExtVar%FluxU
!        FluxV          => Me%ExtVar%FluxV        
!        DZX            => Me%ExtVar%DZX
!        DZY            => Me%ExtVar%DZY
!        DXX            => Me%ExtVar%DXX
!        DYY            => Me%ExtVar%DYY
!        WaterColumn    => Me%ExtVar%WaterColumn
!        WaterColumnOld => Me%ExtVar%WaterColumnOld
!
!        !!! $OMP PARALLEL PRIVATE(I,J,K)
!        !!! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then
!                
!                if (WaterColumn(i,j) .gt. AlmostZero) then
!                    Area_Vertical      = Me%ExtVar%Area(i, j)
!                    !m2 = WC m * face m
!                    Area_Top_U         = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i,j+1))) * DYY(i,j+1)
!                    Area_Top_V         = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i+1,j))) * DXX(i+1,j)
!                    Area_Bottom_U      = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i,j-1))) * DYY(i,j  )
!                    Area_Bottom_V      = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i-1,j))) * DXX(i,j  )
!                    
!                    !s/m3 = s / (m * m2)
!                    aux      = (Me%ExtVar%DT/(WaterColumn(i,j)* Area_Vertical))
!                   
!                   !!!FLUXES WITH Porous Media (infiltration/exfiltration) are not treated here (but in basin)
!                   !!!because Module Basin handles the WC height and concentrations changes and ModuleRunoff
!                   !!!and Module Runoff Properties only transport water horizontally (the remaining water)
!                   
!                   !!!FLUXES WITH Drainage Network 
!                   CofInterfaceDN    = 0.0
!                   ConcInInterfaceDN = 0.0
!                    if (Me%ExtVar%CoupledDN) then
!                       !Flux between river and runoff
!                        if (Me%ExtVar%RiverPoints(I,J) == BasinPoint) then                        
!                            
!                            ! Positive flow -> looses mass
!                            cofInterfaceDN = - aux * Me%ExtVar%FlowToChannels(i,j)
!                            
!                            ! mass going to channel -> conc from runoff
!                            if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
!                                
!                                ConcInInterfaceDN =  CurrProperty%ConcentrationOld(i,j)
!                            
!                                if (Me%CheckGlobalMass) then
!                                    !Global Mass Exchange
!                                    ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
!                                    CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
!                                      + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%ConcentrationOld(i,j)            &
!                                         * 1e-3 * Me%ExtVar%DT)
!                                endif                            
!                            
!                            !mass coming from channel -> conc from DN
!                            elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
!                                
!                                ConcInInterfaceDN = CurrProperty%ConcentrationDN(i,j)
!
!                                if (Me%CheckGlobalMass) then
!                                    !Global Mass Exchange
!                                    ![kg] = [kg] + [m3/s] * [g/m3] * [1e-3kg/g]* [s] 
!                                    CurrProperty%MB%DNExchangeMass =  CurrProperty%MB%DNExchangeMass                   &
!                                      + (Me%ExtVar%FlowToChannels(i,j) * CurrProperty%ConcentrationDN(i,j)            &
!                                         * 1e-3 * Me%ExtVar%DT)
!                                endif                                  
!                                
!                            endif
!                        
!                        endif
!                    endif
!                    
!                    !!DISCHARGE FLUXES IN RUNOFF (NOT YET DONE)
!                    
!                    !!BOUNDARY FLUXES IN RUNOFF (NOT YET DONE)
!                    
!                    !!FLUXES IN X AND Y DIRECTION                        
!                    if (Me%ComputeOptions%AdvDiff_SpatialMethod==AdvDif_CentralDif_) then ! diferenças centrais
!
!                        
!                        AdvTermA_U        = (aux * FluxU(i,j  ) / 2.) 
!                        AdvTermA_V        = (aux * FluxV(i,j  ) / 2.)
!                        AdvTermB_Top_U    = (aux * FluxU(i,j+1) / 2.) 
!                        AdvTermB_Bottom_U = (aux * FluxU(i,j  ) / 2.)
!                        AdvTermB_Top_V    = (aux * FluxV(i+1,j) / 2.) 
!                        AdvTermB_Bottom_V = (aux * FluxV(i  ,j) / 2.)
!                        AdvTermC_U        = (aux * FluxU(i,j+1) / 2.) 
!                        AdvTermC_V        = (aux * FluxV(i+1,j) / 2.)
!                        
!
!                    elseif (Me%ComputeOptions%AdvDiff_SpatialMethod==AdvDif_Upwind_) then ! upwind
!
!                        !DirecU face j
!                        if (FluxU(i,j) .lt. 0.0) then !Left face, Negative - exiting
!                            AdvTermA_U        = 0.0
!                            AdvTermB_Bottom_U = aux * FluxU(i,j)
!                        else !Positive - entering or zero.
!                            AdvTermA_U        = aux * FluxU(i,j)
!                            AdvTermB_Bottom_U = 0.0
!                        endif
!
!                        !DirecV face i
!                        if (FluxV(i,j) .lt. 0.0) then !Left face, Negative - exiting
!                            AdvTermA_V        = 0.0
!                            AdvTermB_Bottom_V = aux * FluxV(i,j)
!                        else !Positive - entering or zero.
!                            AdvTermA_V        = aux * FluxV(i,j)
!                            AdvTermB_Bottom_V = 0.0
!                        endif
!                        
!                        !DirecU face j+1
!                        if (FluxU(i,j+1) .lt. 0.0) then !Right face, Negative - entering
!                            AdvTermC_U        = aux * FluxU(i,j+1)
!                            AdvTermB_Top_U    = 0.0
!                        else !Positive - exiting or zero.
!                            AdvTermC_U        = 0.0
!                            AdvTermB_Top_U    = aux * FluxU(i,j+1)
!                        endif                    
!
!                        !DirecV face i+1
!                        if (FluxV(i+1,j) .lt. 0.0) then !Right face, Negative - entering
!                            AdvTermC_V        = aux * FluxV(i+1,j)
!                            AdvTermB_Top_V    = 0.0
!                        else !Positive - exiting or zero.
!                            AdvTermC_V        = 0.0
!                            AdvTermB_Top_V    = aux * FluxV(i+1,j)
!                        endif                    
!                           
!
!                    endif
!                    
!                    DifTerm_Top_U    = CurrProperty%ViscosityU(i,j+1) * Area_Top_U    * aux / DZX(i,j  )
!                    DifTerm_Bottom_U = CurrProperty%ViscosityU(i,j  ) * Area_Bottom_U * aux / DZX(i,j-1)
!                    DifTerm_Top_V    = CurrProperty%ViscosityV(i+1,j) * Area_Top_V    * aux / DZY(i  ,j)
!                    DifTerm_Bottom_V = CurrProperty%ViscosityV(i  ,j) * Area_Bottom_V * aux / DZY(i-1,j)
!                    
!                    cofA_U = AdvTermA_U                                                          &
!                              + DifTerm_Bottom_U 
!                    
!                    cofA_V = AdvTermA_V                                                          &
!                              + DifTerm_Bottom_V 
!                
!                    cofB_U = - AdvTermB_Top_U                                                    &
!                             + AdvTermB_Bottom_U                                                 &
!                             - DifTerm_Bottom_U                                                  &
!                             - DifTerm_Top_U
!                    
!                    cofB_V = - AdvTermB_Top_V                                                    &
!                             + AdvTermB_Bottom_V                                                 &
!                             - DifTerm_Bottom_V                                                  &
!                             - DifTerm_Top_V        
!                     
!                    cofC_U = - AdvTermC_U                                                        &
!                             + DifTerm_Top_U          
!
!                    cofC_V = - AdvTermC_V                                                        &
!                             + DifTerm_Top_V    
!                    
!                    CofB = ((WaterColumnOld(i,j)*Area_Vertical) / (WaterColumn(i,j)*Area_Vertical)) + cofB_U + cofB_V
!                    
!                    CurrProperty%Concentration(i,j)=  cofA_U * CurrProperty%ConcentrationOld(i,j-1)  &
!                                                     + cofC_U * CurrProperty%ConcentrationOld(i,j+1) &
!                                                     + cofA_V * CurrProperty%ConcentrationOld(i-1,j) &
!                                                     + cofC_V * CurrProperty%ConcentrationOld(i+1,j) &
!                                                     + cofB   * CurrProperty%ConcentrationOld(i,j)   &
!                                                     + CofInterfaceDN * ConcInInterfaceDN  
!!                    if (CurrProperty%Particulate) then
!!                        CurrProperty%BottomConcentration(i,j)=  cofA_U * CurrProperty%BottomConcentrationOld(i,j-1)  &
!!                                                         + cofC_U * CurrProperty%BottomConcentrationOld(i,j+1) &
!!                                                         + cofA_V * CurrProperty%BottomConcentrationOld(i-1,j) &
!!                                                         + cofC_V * CurrProperty%BottomConcentrationOld(i+1,j) &
!!                                                         + cofB   * CurrProperty%BottomConcentrationOld(i,j)   
!!                    endif
!                                                     
!                else !No Volume
!                    CurrProperty%Concentration(i,j) = 0.0                                                         
!                endif
!
!                if (CurrProperty%Particulate) then
!                    ![kg/m2] = [g/m3]* [m * m2] * [1E-3kg/g] /[m2] + [kg/m2]
!                    CurrProperty%TotalConcentration (i,j) = ((CurrProperty%Concentration (i,j) * 1E-3 * WaterColumn(i,j)       &
!                                                             * Me%ExtVar%Area(i, j))  &
!                                                             / Me%ExtVar%Area(i, j)) + CurrProperty%BottomConcentration (i,j)   
!                endif      
!                
!                !Check if any of the coeffs get a negative value. If true, stop program
!                if ((Me%ComputeOptions%AdvDiff_CheckCoefs) .AND. ((cofA_U < 0.0) .OR. (cofB < 0.0) .OR. (cofC_U < 0.0) .OR.    &
!                                                   (cofA_V < 0.0) .OR. (cofC_V < 0.0) )) then
!                                                   
!                    call LogCoefs(i,j,cofA_U,cofB,cofC_U,cofA_V,cofC_V)
!                
!                endif
!              
!            endif
!        enddo
!        enddo
!        !!! $OMP END DO
!        !!! $OMP END PARALLEL
!
!
!    end subroutine ModifyAdvectionDiffusion_Explicit
!    
!    !---------------------------------------------------------------------------
    
    subroutine ModifyDiffusivity_New(CurrProperty)
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        !Local-----------------------------------------------------------------
        integer                                     :: I,J, CHUNK
        real                                        :: DiffCoef, velU, VelV        
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyDiffusivity_New")

        
        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J,DiffCoef,VelU,VelV)

        DiffCoef  = CurrProperty%Evolution%AdvDiff%Molecular_Diff_Coef
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then
                
                if ((Me%ExtVar%WaterColumnOld(i,j) .gt. AlmostZero) .and. (Me%ExtVar%WaterColumnOld(i,j-1) .gt. AlmostZero)) then
                    VelU = 0.5 * (Me%ExtVar%CenterVelU(i,j)+Me%ExtVar%CenterVelU(i,j-1))
                
                    CurrProperty%ViscosityU(i,j)  = (DiffCoef +(abs(VelU) * Me%Disper_Trans%Field(i,j)))
                else
                    CurrProperty%ViscosityU(i,j)  = 0.0
                endif                                                   
                
                if ((Me%ExtVar%WaterColumnOld(i,j) .gt. AlmostZero) .and. (Me%ExtVar%WaterColumnOld(i-1,j) .gt. AlmostZero)) then
                    VelV = 0.5 * (Me%ExtVar%CenterVelV(i,j)+Me%ExtVar%CenterVelV(i-1,j))
                    
                    CurrProperty%ViscosityV(i,j)  = (DiffCoef +(abs(VelV) * Me%Disper_Trans%Field(i,j)))
                else
                    CurrProperty%ViscosityV(i,j)  = 0.0
                endif                                                    

            endif                                                                        
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyDiffusivity_New")

    
    end subroutine ModifyDiffusivity_New

    !---------------------------------------------------------------------------
    
!    subroutine SoilQualityProcesses
!
!        !External--------------------------------------------------------------
!        integer                                 :: STAT_CALL        
!        
!        !Local----------------------------------------------------------------- 
!        type (T_Property),          pointer     :: PropertyX
!        type (T_Property),          pointer     :: SoilDryDensity, Salinity, pH
!        type (T_Property),          pointer     :: IonicStrength, PhosphorusAdsortionIndex
!
!!        type (T_SoilRate),      pointer     :: SoilRateX
!!        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB 
!!        integer                                 :: i, j, k
!        
!        !Begin-----------------------------------------------------------------
!        
!!        WIUB = Me%WorkSize%IUB
!!        WJUB = Me%WorkSize%JUB
!!        WILB = Me%WorkSize%ILB
!!        WJLB = Me%WorkSize%JLB
!!        WKUB = Me%WorkSize%KUB
!!        WKLB = Me%WorkSize%KLB
!        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "SoilQualityProcesses")
!
!        call ComputeDissolvedToParticulate2D
!        
!        !Properties not modified by sediment quality (not state variables) but needed in argument
!        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property soil dry density not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR00'
!        endif
!
!        call SearchProperty(Salinity, Salinity_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property salinity not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR10'
!        endif
!
!        call SearchProperty(pH, pH_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property pH not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR20'
!        endif
!
!        call SearchProperty(IonicStrength, IonicStrength_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property ionic strength not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR30'
!
!        endif
!
!        call SearchProperty(PhosphorusAdsortionIndex, PhosphorusAdsortionIndex_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property phosphorus sdsortion index not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR40'
!        endif
!
!        call ComputeWindVelocity
! 
!        
!        if (Me%ExtVar%Now .GE. Me%Coupled%SoilQuality_NextCompute) then
!            
!            PropertyX => Me%FirstProperty
!
!            do while(associated(PropertyX))
!                
!
!!                call Modify_Interface(InterfaceID               = Me%ObjInterface,                         &
!!                                      PropertyID                = PropertyX%ID%IDNumber,                   &
!!                                      Concentration             = PropertyX%Concentration,                 &
!!                                      WaterPoints2D             = Me%ExtVar%BasinPoints,                   &
!!                                      OpenPoints2D              = Me%ExtVar%BasinPoints,                   &
!!                                      WaterPercentage           = Me%ExtVar%WaterContent,                  &
!!                                      DissolvedToParticulate3D  = Me%DissolvedToParticulate2D,             &
!!                                      SoilDryDensity            = SoilDryDensity%Concentration,            &
!!                                      Salinity                  = Salinity%Concentration,                  &
!!                                      pH                        = pH%Concentration,                        &
!!                                      IonicStrength             = IonicStrength%Concentration,             &
!!                                      PhosphorusAdsortionIndex  = PhosphorusAdsortionIndex%Concentration,  &
!!                                      WindVelocity              = Me%ExtVar%WindVelocity,                  &
!!                                      STAT                      = STAT_CALL)
!!                if (STAT_CALL .NE. SUCCESS_)                                                               &
!!                    stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR01'
!                
!
!                PropertyX => PropertyX%Next
!                
!
!            end do
!            
!            Me%Coupled%SoilQuality_NextCompute = Me%Coupled%SoilQuality_NextCompute +       &
!                                                     Me%Coupled%SoilQuality_DT
!
!        end if
!
!        PropertyX => Me%FirstProperty
!
!        do while(associated(PropertyX))
!
!
!            if (PropertyX%Evolution%SoilQuality) then
!
!                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then
!
!                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
!                                          PropertyID    = PropertyX%ID%IDNumber,            &
!                                          Concentration = PropertyX%Concentration,          &
!                                          WaterPoints2D = Me%ExtVar%BasinPoints,          &
!                                          DTProp        = PropertyX%Evolution%DTInterval,   &
!                                          STAT          = STAT_CALL)
!                    if (STAT_CALL .NE. SUCCESS_)                                            &
!                        stop 'SoilQuality_Processes - ModuleRunoffProperties - ERR02'
!
!                end if
!
!            end if
!
!            PropertyX => PropertyX%Next
!            
!        end do
!
!        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "SoilQualityProcesses")
!
!    end subroutine SoilQualityProcesses
!    
!    !-----------------------------------------------------------------------------    
!
!    subroutine ComputeWindVelocity
!
!        !Arguments-------------------------------------------------------------
!
!        !External--------------------------------------------------------------
!        integer                                      :: i, j, k
!        
!        !Begin-----------------------------------------------------------------
!
!        call SetMatrixValue(Me%ExtVar%WindVelocity, Me%Size, 0.0, Me%ExtVar%BasinPoints)
!
!        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
!
!            if (Me%ExtVar%BasinPoints(i, j) == 1) then
!                
!                ! km/day                        =  m/s  * 1E-3km/m * 86400s/day 
!                Me%ExtVar%WindVelocity(i,j) = Me%ExtVar%WindVelocity2D(i,j) * 1E-3 * 86400
!            endif
!
!        enddo
!        enddo
!
!    end subroutine ComputeWindVelocity
!
!    !----------------------------------------------------------------------------- 

    subroutine ModifyBottomFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ModifyBottomFluxes")

        if (Me%Coupled%SplashErosionFluxes) then
        
            call ModifySplashErosionFluxes

        endif

            
        call ModifyShearStress

!        if (Me%ComputeOptions%ErosionDepositionModel == Krone_) then
            if (Me%Coupled%ErosionFluxes) then
            
                call ComputeErosionFluxes
            
            endif
            
            if (Me%Coupled%DepositionFluxes) then
                
                call ComputeDepositionFluxes
            
            endif

!        elseif (Me%ComputeOptions%ErosionDepositionModel == SomethingElse_) then
!        
!        endif

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ModifyBottomFluxes")

    
    end subroutine ModifyBottomFluxes
    
    !---------------------------------------------------------------------------

    subroutine ModifySplashErosionFluxes
       
       !Local-----------------------------------------------------------------
        type (T_Property) , pointer                 :: Property
        real(8), dimension(:,:), pointer            :: ThroughFall, WaterColumn
        real, dimension(:,:), pointer               :: CanopyHeight, BottomSedimentConc
        real(8), dimension(:,:), pointer            :: CanopyDrainage
        real                                        :: KE_Leaf_Drainage, KE_ThroughFall
        real                                        :: SplashRate, CanopyDrain, DirectRain, DirectRainRate        
        real                                        :: SplashErodedMass, RainKineticEnergy, SplashConc
        real                                        :: BottomArea,WaterVolume
        integer                                     :: STAT_CALL, i ,j 
        !---------------------------------------------------------------------
    
        call SearchProperty(Property, PropertyXIDNumber = TSS_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call SearchProperty(Property, PropertyXIDNumber = Cohesive_Sediment_)
        endif
        
        ThroughFall        => Me%ExtVar%ThroughFall
        CanopyHeight       => Me%ExtVar%CanopyHeight
        CanopyDrainage     => Me%ExtVar%CanopyDrainage
        WaterColumn        => Me%ExtVar%WaterColumn
        BottomSedimentConc => Property%BottomConcentration        
        
        nullify (Property)
        Property => Me%FirstProperty                                                    
        do while (associated (Property)) 
 
if1:        if (Property%Evolution%BottomFluxes   .and. Property%Evolution%SplashErosion  &
                .AND. Property%ID%IDNumber /= VSS_ .AND. Property%ID%IDNumber /= TSS_) then
         
                KE_Leaf_Drainage = 0.0
                KE_ThroughFall   = 0.0
                SplashErodedMass = 0.0
                SplashRate       = 0.0
                SplashConc       = 0.0
                
                Property%SplashRate       = 0.0
                Me%RainKineticRate        = 0.0
                
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
if2:                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then   
                        
                        !Only compute if there is water to recieve flux or if there is sediment to take
                        if (WaterColumn(i,j) > AllmostZero .and. BottomSedimentConc (i,j) > AllmostZero) then
                        
                            if (CanopyHeight(i,j) .LT. 0.15) then
                           
                                KE_Leaf_Drainage = 0.0
                               
                            else
                                ![mm] = [m] * 1E3 mm/m
                                CanopyDrain = CanopyDrainage(i,j) * 1E3
                               ! [J.m-2] = [mm] . [J.m-2.mm-1] ! units not clearly defined using canopy height in [m]) 
                               !! see RPC Morgan 2007 in Earth Surface Processes and Landforms and RPC Morgan 1998 in Earth 
                               !Surface Processes and Landforms
                                KE_Leaf_Drainage = CanopyDrain * (( 15.8 * CanopyHeight(i,j)**0.5 ) - 5.87)
                               
                            endif 
                            
                            !Direct Throughfall is Troughfall - Canopy Drainage (see basin formulation)
                            ! [mm] =  ( [m] - [m]) * 1E+3 mm/m 
                            DirectRain       = (ThroughFall(i,j) - CanopyDrainage(i,j)) * 1E+3
                            
                            if (Me%ComputeOptions%Splash_ErosiveRainMethod == 1) then
                                ! [J.m-2]  =  [mm] . [J.m-2.mm-1] ! units not clearly defined using erosive rain in [mm.h-1]) 
                                !! see RPC Morgan 2007 in Earth Surface Processes and Landforms and RPC Morgan 1998 in Earth 
                                !Surface Processes and Landforms                            
                                KE_ThroughFall = DirectRain * ( 8.95 + 8.44 * log10(Me%Splash_ErosiveRainValue))
                            else
                                ![mm.h-1] = [mm] * [s] * 3600 s/h
                                DirectRainRate   = DirectRain * Me%ExtVar%DT * 3600.
                               
                                ! [J.m-2]  =  [mm] . [J.m-2.mm-1] ! units not clearly defined using rain in [mm.h-1]) 
                                !! see RPC Morgan 2007 in Earth Surface Processes and Landforms and RPC Morgan 1998 in Earth 
                                !Surface Processes and Landforms
                                !KE_ThroughFall = ThroughFall * ( 8.95 + 8.44 * log10(Me%ErosiveRain%Field(i,j)) )
                                if (DirectRainRate .gt. 0.0) then
                                    KE_ThroughFall = DirectRain * ( 8.95 + 8.44 * log10(DirectRainRate))
                                    !For small rains expression gets negative
                                    if (KE_ThroughFall .lt. 0.0) then
                                        KE_ThroughFall = 0.0
                                    endif
                                else
                                    KE_ThroughFall = 0.0
                                endif
                            endif
                           
                            ! [J.m-2] or [kg.s-2]
                            RainKineticEnergy = KE_Leaf_Drainage + KE_Throughfall
                            
                            ! [J.s-1.m-2] or [W/m2] or [kg.s-3] =  [J.m-2] / [s]
                            Me%RainKineticRate  (i,j) =  RainKineticEnergy / Me%ExtVar%DT

                            ! [kg.m-2] =  [g.J-1] . [] . [1E-3kg/g] . [] * [J.m-2] 
    !                        SplashFlux = (Me%Kclay%Field(i,j) * Me%ClayFraction%Field(i,j) * 1E-3 *   &
!                            (1 - Me%StoneFraction%Field(i,j)) * Me%RainKineticEnergy(i,j))
                            
                            ! [kg.m-2.s-1] =  [g.J-1] * [1E-3kg/g] * [J.s-1.m-2] * [kg/m2 Prop]/[kg/m2 Sed] 
                            SplashRate = (Me%Kclay%Field(i,j) * 1E-3 * Me%RainKineticRate(i,j))   &
                                          * Property%BottomConcentration (i,j)                    &
                                          / BottomSedimentConc (i,j)                              
                            
                            if (Me%ComputeOptions%Splash_CriticalHeight) then
                                ! Correction of the splash flux by the water colomn level. Exponential decay gets inflexion point
                                ! when water column is height HcriticSplash
                                ! [kg.m-2.s-1] =  [kg.m-2.s-1] * [] 
                                SplashRate = SplashRate * exp(- WaterColumn(i,j) / Me%HcriticSplash)
                            endif
                            
                            ![kg.m-2] = [kg.m-2.s-1] * [s]
                            SplashConc = SplashRate * Me%ExtVar%DT
                           
                            if (SplashConc < Property%BottomConcentration (i,j)) then

                                Property%BottomConcentration (i,j) = Property%BottomConcentration (i,j)       &
                                                                     - SplashConc

                            else

                                SplashConc = Property%BottomConcentration (i,j) - Property%BottomMinConc 
                                Property%BottomConcentration (i,j) = Property%BottomMinConc
                                ![kg m-2 s-1] = [kg m-2] / [s]                      
                                SplashRate = SplashConc / Me%ExtVar%DT
                            
                            end if                            

                            BottomArea = Me%ExtVar%Area(i,j)
                                
                            ![g]  = [kg m-2] * [m2] * [1000g/kg]
                            SplashErodedMass  = SplashConc * BottomArea * 1E3  
                                                               
                           ![m3] = [m] * [m2]
                            WaterVolume = WaterColumn(i,j) * BottomArea 
                            
                            ![g/m3] = (([g/m3] * [m3]) + [g]) / [m3]
                            Property%Concentration (i,j) = ((Property%Concentration (i,j) * WaterVolume)  &
                                                                + SplashErodedMass) / WaterVolume 

                            Property%SplashRate(i,j) = SplashRate

                            ![kg/m2] = [g/m3]* [m3] * [1E-3kg/g] /[m2] + [kg/m2]
                            Property%TotalConcentration (i,j) = Property%Concentration (i,j) * 1E-3 * WaterVolume / BottomArea &
                                                                + Property%BottomConcentration (i,j)            
                                                                
                        endif
 
                    end if if2
 
                enddo
                enddo
            
            endif if1
 
            Property => Property%Next
        enddo        
                
    
    end subroutine ModifySplashErosionFluxes   
    
    !---------------------------------------------------------------------------

    subroutine ModifyShearStress

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i, j
        real                                        :: Chezy, CenterVelocity, FinalManning
!        real                                        :: Area_U, Area_V
!        real                                        :: WetPerimeter_U, WetPerimeter_V
!        real                                        :: HydraulicRadius
        real, dimension(:,:), pointer               :: CenterVelocityX, CenterVelocityY
        real, dimension(:,:), pointer               :: DUX, DVY
        real(8), dimension(:,:), pointer            :: WaterColumn
        real, dimension(:,:), pointer               :: Manning
        real, dimension(:,:), pointer               :: ManningDelta
        !----------------------------------------------------------------------
        
!        if (Me%ComputeOptions%ShearModel == Chezy_) then
            DUX            => Me%ExtVar%DUX
            DVY            => Me%ExtVar%DVY
            !The new water column computed in module runoff will be used - water column where center velocity field was computed
            !In cells that loose all water in the time step shear will be zero (and center velocity from runoff) 
            !but in fact they had flux field
            !and a shear 
            WaterColumn => Me%ExtVar%WaterColumn
            
            call GetManning (ObjRunOffID = Me%ObjRunoff, Manning = Manning, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_)  stop 'ModifyShearStress - ModuleRunoffProperties - ERR01'        
  
  
            call GetManningDelta (ObjRunOffID = Me%ObjRunoff, ManningDelta = ManningDelta, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_)  stop 'ModifyShearStress - ModuleRunoffProperties - ERR05' 
            
            call GetRunoffCenterVelocity (Me%ObjRunoff, CenterVelocityX, CenterVelocityY, STAT_CALL)  
            if (STAT_CALL /= SUCCESS_)  stop 'ModifyShearStress - ModuleRunoffProperties - ERR010'   
            
            Me%ShearStress(:,:) = 0.0
                      
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB

                if ((Me%ExtVar%BasinPoints(I,J) == BasinPoint) .and. (WaterColumn(i,j) .gt. Me%HminChezy)) then  
                    
                   
                    !To compute overland resistance in bottom for shear computation (erosion/deposition).
                    !This process was created to remove from manning the resistance given by 
                    !aerial vegetation parts that affect flow but do not affect bottom shear. Without that,
                    !a manning increase (e.g. forestation) in one cell increases water depth (and reduces velocity)
                    !but may increase shear stress (because water height increase is transformed in bottom resistance 
                    !using manning - chezy)
                    !Manning change if user loaded values in runoff else ManningDelta will be zero and no change                    
                    FinalManning = Manning(i,j) - ManningDelta(i,j)
                    
                    !Hydraulic radius is water column in runoff
                    if (WaterColumn(i,j) > Me%HminChezy) then
                        ![-]  = [m.s-2] * [(s.m(-1/3))^2 == s2.m(-2/3)] / [m(-1/3)] 
                        Chezy = Gravity * FinalManning**2.0 / (WaterColumn(i,j)** (1./3.))
                    else 
                        Chezy = 0.0
                    end if
                    
                    !m/s - Attention CenterVelocities in Runoff are computed with fluxes and water column new
                    CenterVelocity       = sqrt(CenterVelocityX(i,j)**2.0 + CenterVelocityY(i,j)**2.0)
                    ![Pa] or [kg/m.s2]   = [kg/m3] * [-] * [m2/s2]
                    Me%ShearStress (i,j) = SigmaDensityReference * Chezy * CenterVelocity**2.0
                
                else
                
                    Me%ShearStress (i,j) = 0.0
                                
                end if

            enddo
            enddo


            call UnGetRunoff (Me%ObjRunoff, Manning, STAT_CALL)  
            if (STAT_CALL /= SUCCESS_)  stop 'ModifyShearStress - ModuleRunoffProperties - ERR020'        
 
            call UnGetRunoff (Me%ObjRunoff, ManningDelta, STAT_CALL)  
            if (STAT_CALL /= SUCCESS_)  stop 'ModifyShearStress - ModuleRunoffProperties - ERR025'        
            
            call UnGetRunoff (Me%ObjRunoff, CenterVelocityX, STAT_CALL)  
            if (STAT_CALL /= SUCCESS_)  stop 'ModifyShearStress - ModuleRunoffProperties - ERR030'        

            call UnGetRunoff (Me%ObjRunoff, CenterVelocityY, STAT_CALL)  
            if (STAT_CALL /= SUCCESS_)  stop 'ModifyShearStress - ModuleRunoffProperties - ERR040'        

       
!        elseif (Me%ComputeOptions%ShearModel == SomethingElse_) then
!
!        endif
        
    end subroutine ModifyShearStress

    !---------------------------------------------------------------------------

    subroutine ComputeErosionFluxes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property) , pointer                 :: Property
        real, dimension(:,:), pointer               :: BottomSedimentConc
        real(8), dimension(:,:), pointer            :: WaterColumn
        real                                        :: ErosionRate, WaterVolume
        real                                        :: ErodedConc, ErodedMass        
        real                                        :: BottomArea
        real                                        :: aux
        integer                                     :: STAT_CALL, i ,j 
        !---------------------------------------------------------------------

        call SearchProperty(Property, PropertyXIDNumber = TSS_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call SearchProperty(Property, PropertyXIDNumber = Cohesive_Sediment_)
        endif
        
        BottomSedimentConc => Property%BottomConcentration
        WaterColumn => Me%ExtVar%WaterColumn

        nullify (Property)
        Property => Me%FirstProperty                                                    
        do while (associated (Property))

if1:        if (Property%Evolution%BottomFluxes  .and. Property%Evolution%Erosion  &
                .AND. Property%ID%IDNumber /= VSS_ .AND. Property%ID%IDNumber /= TSS_) then

                Property%ErosionRate = 0.0

                do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                do I = Me%WorkSize%ILB, Me%WorkSize%IUB

                    ErodedConc  = 0.0        
                    ErosionRate = 0.0

if2:                if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then
                   
if3:                    if (Me%ShearStress (i,j) > Me%ErosionCriticalShear%Field(i,j) .and.  &
                            WaterColumn(i,j) > Me%HminChezy .and.                            &
                            BottomSedimentConc (i,j) > AllmostZero) then

                            aux = Me%ShearStress (i,j) / Me%ErosionCriticalShear%Field(i,j) - 1.0

                            ![kg m-2 s-1] = [kg m-2 s-1] * [(Kg m-2 Prop / Kg m-2 Cohesive)] * [-]
                            ErosionRate = Me%ErosionCoefficient%Field(i,j)                  &
                                        * Property%BottomConcentration (i,j)                &
                                        / BottomSedimentConc (i,j)                          &
                                        * aux 
                            
                            ![kg m-2] = [kg m-2 s-1] * [s]
                            ErodedConc = ErosionRate * Me%ExtVar%DT
                                                  
                            if (ErodedConc < Property%BottomConcentration (i,j) - Property%BottomMinConc) then

                                Property%BottomConcentration (i,j) = Property%BottomConcentration (i,j)       &
                                                                     - ErodedConc

                            else

                                ErodedConc = Property%BottomConcentration (i,j) - Property%BottomMinConc 
                                Property%BottomConcentration (i,j) = Property%BottomMinConc
                                ![kg m-2 s-1] = [kg m-2] / [s]                      
                                ErosionRate = ErodedConc / Me%ExtVar%DT
                            
                            end if

                            BottomArea = Me%ExtVar%Area(i,j)
                            
                            ![g]  = [kg m-2] * [m2] * [1000g/kg]
                            ErodedMass  = ErodedConc * BottomArea * 1E3
                            
                            !m3 = m * m2
                            WaterVolume = WaterColumn(i,j) * BottomArea
                            
                            ![g/m3] = (([g/m3] * [m3]) + [g]) / [m3]
                            Property%Concentration (i,j) = ((Property%Concentration (i,j) * WaterVolume) &
                                                            + ErodedMass) / WaterVolume 
                            ![kg m-2 s-1]
                            Property%ErosionRate (i,j) = ErosionRate
                            
                            ![kg/m2] = [g/m3]* [m3] * [1E-3kg/g] /[m2] + [kg/m2]
                            Property%TotalConcentration (i,j) = Property%Concentration (i,j) * 1E-3 * WaterVolume / BottomArea &
                                                                + Property%BottomConcentration (i,j) 
                                                              
                        end if if3
                    end if if2

                enddo
                enddo
            
            endif if1

            Property => Property%Next
        enddo

    end subroutine ComputeErosionFluxes
    !---------------------------------------------------------------------------
    
    subroutine ComputeDepositionFluxes
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Property) , pointer                 :: Property
        real(8), dimension(:,:), pointer            :: WaterColumn
        real, dimension(:,:), pointer               :: WaterSedimentConc
        real                                        :: WaterVolume, BottomArea
        real                                        :: DepositionRate, DepositedMass        
        real                                        :: Mass, SPMConc
        real                                        :: aux, MinimumMass
        integer                                     :: STAT_CALL, i ,j
        !Begin------------------------------------------------------------------

        call SearchProperty(Property, PropertyXIDNumber = TSS_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call SearchProperty(Property, PropertyXIDNumber = Cohesive_Sediment_)
        endif
        WaterSedimentConc  => Property%Concentration
        
        WaterColumn        => Me%ExtVar%WaterColumn
 
        nullify (Property)
        Property => Me%FirstProperty                                                    
        
        do while (associated (Property))
if1:        if (Property%Evolution%BottomFluxes   .and.   Property%Evolution%Deposition    &
                .AND. Property%ID%IDNumber /= VSS_ .AND. Property%ID%IDNumber /= TSS_) then
                
                Property%DepositionRate = 0.0
                
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    DepositedMass  = 0.0
                    DepositionRate = 0.0

if2:                if ((Me%ExtVar%BasinPoints(i,j) == BasinPoint)) then

if3:                    if (Me%ShearStress (i,j) < Me%DepositionCriticalShear%Field(i,j) .and.                          &
                            WaterColumn(i,j) > AllmostZero) then

                            !ModifySettlingVelocity - [m.s-1]
                            if (Property%Ws_Type == SPMFunction) then
                                
                                SPMConc = Property%Concentration (i,j)
                                
                                Property%Ws (i,j) = SettlingVelocity_RP(SPMConc,                            &
                                                                         Property%CHS, Property%KL,         &
                                                                         Property%KL1, Property%M,          &
                                                                         Property%ML, i, j)
                            else                       
                                Property%Ws (i,j) = Property%Ws_Value
                            end if
                            
                            ![-] = [Pa].[Pa]
                            aux = 1.0 - Me%ShearStress (i,j) / Me%DepositionCriticalShear%Field(i,j) 
                            
                            ![kg.m-2.s-1] = [g.m-3] * [m.s-1] * [1E-3kg/g] * [-] * [(mg.l-1 Prop / mg.l-1 Cohesive)]
                            DepositionRate = Property%Concentration (i,j) * 1E-3 * Property%Ws (i,j) * aux   ! &
                                             !* Property%Concentration (i,j)                            &
                                             !/ WaterSedimentConc (i,j)                                 &
                            
                            BottomArea = Me%ExtVar%Area(i,j)
                            
                            ![g] = [kg.m-2.s-1] * 1E3[g/kg] * [m2] * [s]
                            DepositedMass = DepositionRate * 1E3 * BottomArea *  Me%ExtVar%DT
                            
                            ![m3] = [m] * [m2]
                            WaterVolume = WaterColumn(i,j) * BottomArea                    
                            
                            ![g] = [g.m-3] * [m3]
                            Mass = Property%Concentration (i,j) * WaterVolume
                            
                            if (DepositedMass < Mass ) then
                                ![g.m-3] = [g.m-3] - [g] / [m3]
                                Property%Concentration (i,j) = Property%Concentration (i,j)                 &
                                                                - DepositedMass / WaterVolume               
                            else
                                if(Property%Evolution%MinConcentration) then
                                    ! [g] = [g/m3] * [m3]
                                    MinimumMass = Property%MinValue * WaterVolume
                                    DepositedMass = Mass - MinimumMass
                                    Property%Concentration (i,j) = Property%MinValue
                                else
                                    DepositedMass = Mass
                                    Property%Concentration (i,j) = 0.0
                                endif
                                
                                ![kg/m2.s] = [g] * 1E-3[kg/g]  / ([m2]* [s])
                                DepositionRate = DepositedMass * 1E-3 / (BottomArea * Me%ExtVar%DT)
                                
                            endif
                            
                            ![kg/m2] = [kg/m2] + ([g] * 1E-3 kg/g) / [m2]
                            Property%BottomConcentration (i,j) = Property%BottomConcentration (i,j)       &
                                                                 + (DepositedMass * 1E-3) / BottomArea  
                            
                            Property%DepositionRate (i,j) = DepositionRate
                            
                            ![kg/m2] = [g/m3]* [m3] * [1E-3kg/g] /[m2] + [kg/m2]
                            Property%TotalConcentration (i,j) = Property%Concentration (i,j) * 1E-3 * WaterVolume   &
                                                                / Me%ExtVar%Area(i,j)                               &
                                                                + Property%BottomConcentration (i,j)                          
                        
                        end if if3
                    end if if2                    
                
                end do    
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
    real function SettlingVelocity_RP (ConcentrationIN, CHS, KL, KL1, M, ML, i, j)      
        !Arguments-------------------------------------------------
        real,       intent(IN) :: ConcentrationIN    !g/m3
        real,       intent(IN) :: CHS              !Hindered settling concentration
        real,       intent(IN) :: KL, KL1, M, ML
        integer,    intent(IN) :: i,j 
        real                   :: Concentration       

        !Local-----------------------------------------------------
        real                    :: Aux
        
        !Begin-----------------------------------------------------

        !kg/m3 = g/m3 * 1E-3 kg/g
        Concentration = ConcentrationIN * 1E-3
        
        Aux = KL1 * (Concentration - CHS)

        if (Concentration < CHS .and. Concentration >= 0.) then 

            SettlingVelocity_RP = KL*(Concentration)**M

        elseif(Aux < 1. .and. Aux >= 0.) then

            SettlingVelocity_RP = KL*(CHS)**M*(1.0-Aux)**ML 

        elseif(Aux > 1. .and. Concentration < 100000.) then

            SettlingVelocity_RP = 0. !if concentration is to high settling velocity is null

        else

            write(*,*)'Concentration (g/l)          = ', Concentration
            write(*,*)'KL1 * (Concentration - CHS)  = ', Aux
            write(*,*)'Cell                         = ', i,j
            stop 'Error computing the settling velocity - SettlingVelocity - ModuleRunoffProperties'

       endif

    end function SettlingVelocity_RP

    !---------------------------------------------------------------------------

    subroutine DecayProcesses
    
        !Local--------------------------------------------------------------------
        type (T_Property), pointer                         :: Property
        real                                               :: DT
        real(8)                                            :: WaterVolume
        real(8)                                            :: OldMass
        real(8)                                            :: NewMass
        real(8)                                            :: MassSink
        integer                                            :: i,j
        !Begin--------------------------------------------------------------------
        
        
        Property => Me%FirstProperty  

do1 :   do while (associated(Property))    

            if (Property%Evolution%Decay) then
                
                !days
                !if DTInterval, only update at given time
                if (Property%Evolution%DTIntervalAssociated) then
                    DT = Property%Evolution%DTInterval / 86400.
                else !update every time
                    Property%Evolution%NextCompute = Me%ExtVar%Now
                    DT = Me%ExtVar%DT /86400.
                endif
                    
                if(Me%ExtVar%Now .GE. Property%Evolution%NextCompute) then            
            
                    do J = Me%WorkSize%JLB, Me%WorkSize%JUB       
                    do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                        if ((Me%ExtVar%BasinPoints(I,J) == WaterPoint) .and. (Me%ExtVar%WaterColumn(i,j) .gt. 0.)) then
                            
                            !WQ process without volume change. only mass change
                            !m3 H20
                            WaterVolume = Me%ExtVar%WaterColumn(i,j) * Me%ExtVar%Area(i, j)
                            
                            !g  = g/m3 * m3
                            OldMass = Property%Concentration(I,J) * WaterVolume
                            
                            !P = P0*exp(-kt)
                            !g
                            MassSink = min (OldMass - OldMass * exp(-Property%Evolution%DecayRate * DT),  OldMass)
                            
                            NewMass = OldMass - MassSink
                            
                            Property%Concentration(I,J) = NewMass / WaterVolume
                            
                        endif
                    
                    enddo
                    enddo
                
                endif
                
            endif     

            Property => Property%Next
        end do do1   

        nullify(Property)
        
            
    end subroutine DecayProcesses

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
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: Couple_ID


        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "Partition_Processes")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 

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
                        stop 'Partition_Processes - ModuleRunoffProperties - ERR02'
                    
                    if(PropertyX%Evolution%Partition%UseSedimentRefConc)then
                    
                        call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'Partition_Processes - ModuleRunoffProperties - ERR03'
                    
                    endif

    !                if (PropertyX%Evolution%Partition%SalinityEffect) then
    !
    !                    call Search_Property(Salinity, PropertyXID = Salinity_, STAT = STAT_CALL)
    !                    if (STAT_CALL /= SUCCESS_)                                      &
    !                        stop 'Partition_Processes - ModuleRunoffProperties - ERR03'
    !                endif
                    

    do1:            do j = JLB, JUB
    do2:            do i = ILB, IUB
        
    cd2:                if ((Me%ExtVar%BasinPoints(i,j) == BasinPoint)) then

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

                                RefSedFactor = CohesiveSediment%Concentration(i,j) / &
                                               PropertyX%Evolution%Partition%SedimentRefConc

                                if(RefSedFactor < 1.)then

                                    TransferRate = PropertyX%Evolution%Partition%Rate * RefSedFactor
                                
                                else

                                    TransferRate = PropertyX%Evolution%Partition%Rate

                                end if

                            else

                                TransferRate = PropertyX%Evolution%Partition%Rate

                            end if

                            ! [g/m3]       =          [d]         * [d^-1]        * [g/m3]
                            MassTransfer    =         DT * TransferRate *          &
                            (DissolvedFraction   * PartPropX%Concentration(i, j) -        &                  
                             ParticulateFraction * PropertyX%Concentration(i, j))

                            if ((MassTransfer .gt. 0.0) .and. (MassTransfer .gt. PartPropX%Concentration(i, j))) then
                                MassTransfer = PartPropX%Concentration(i, j)
                            elseif ((MassTransfer .lt. 0.0) .and. (-MassTransfer .gt. PropertyX%Concentration(i, j))) then
                                MassTransfer = -PropertyX%Concentration(i, j)
                            endif

                            PartPropX%Concentration(i, j) =                               &
                                               PartPropX%Concentration(i, j) - MassTransfer 

                            PropertyX%Concentration(i, j) =                               &
                                               PropertyX%Concentration(i, j) + MassTransfer

                        endif cd2
                    enddo do2
                    enddo do1
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
!    subroutine SetLimitsConcentration (Message)

        !Arguments-------------------------------------------------------------
!        character(len=*)                            :: Message
        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: i, j !, CHUNK
!        character(len=5)                            :: char_i, char_j
!        character(len=20)                           :: char_conc
!        character (len = StringLength)              :: StrWarning  
                
        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "SetLimitsConcentration")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 

        Property => Me%FirstProperty  

do1 :   do while (associated(Property))
cd1 :       if (Property%Evolution%MinConcentration) then
                

                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%BasinPoints(i, j) == 1) then
    
                        if (Property%Concentration(i, j) < Property%MinValue) then
                            
                            ! mass created (g) = g + (g/m3)* (m * m2)
                            Property%Mass_created(i, j) = Property%Mass_Created(i, j)   +  &
                                                   (Property%MinValue                -  &
                                                    Property%Concentration(i, j)) *  (Me%ExtVar%WaterColumn(i,j) * &
                                                    Me%ExtVar%Area (i, j))

                            Property%Concentration(i, j) = Property%MinValue
!
!                            write(char_i, '(i4)')i
!                            write(char_j, '(i4)')j
!                            write(char_conc, '(f20.8)') Property%Concentration(i,j) 
!                            StrWarning = trim(Property%ID%Name)//' was modified to its MinValue in cell(i,j)'// &
!                                                               char_i//','//char_j//' '//char_conc//' '//Message
!
!                            call SetError(WARNING_, INTERNAL_, StrWarning, OFF)
                            
                        endif

                    endif

                enddo
                enddo
                
                
            endif cd1
                
        Property => Property%Next
        end do do1

        nullify(Property)

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "SetLimitsConcentration")


    end subroutine SetLimitsConcentration

    !--------------------------------------------------------------------------

    subroutine WarnOnNegativeValues (Message)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Message
        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: i, j !, CHUNK
        character(len=5)                            :: char_i, char_j
        character(len=15)                           :: char_conc
        character (len = StringLength)              :: StrWarning        
        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "WarnOnNegativeValues")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 

        Property => Me%FirstProperty  

do1 :   do while (associated(Property))

cd1 :       if (Property%Evolution%WarnOnNegativeValues) then
                
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%BasinPoints(i, j) == 1) then
    
                        if (Property%Concentration(i, j) < 0.0) then
                            
                            write(char_i, '(i4)')i
                            write(char_j, '(i4)')j
                            write(char_conc, '(ES10.3)') Property%Concentration(i,j) 
                            StrWarning = trim(Property%ID%Name)//' has a negative concentration in cell(i,j)'// &
                                                               char_i//','//char_j//' '//char_conc//' '//Message

                            call SetError(WARNING_, INTERNAL_, StrWarning, OFF)
                            
                        endif

                    endif

                enddo
                enddo
                
                
            endif cd1
                
        Property => Property%Next
        end do do1

        nullify(Property)

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "WarnOnNegativeValues")


    end subroutine WarnOnNegativeValues

    !--------------------------------------------------------------------------
    
!    subroutine ComputeDissolvedToParticulate2D
!
!        !External--------------------------------------------------------------
!        integer                                 :: STAT_CALL        
!         
!        !Local----------------------------------------------------------------- 
!        integer                                 :: i, j, k
!        real                                    :: DT, InstantValue, ResidualValue
!        type(T_Property), pointer               :: SoilDryDensity!, DrySedimentVolume
!        !Begin-----------------------------------------------------------------
!        
!        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ComputeDissolvedToParticulate3D")
!
!
!        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)  &
!            stop 'ComputeDissolvedToParticulate3D - ModuleRunoffProperties - ERR01'
!
!        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) stop 'ComputeDissolvedToParticulate3D - ModuleRunoffProperties - ERR02'
!
!
!        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if(Me%ExtVar%BasinPoints(i,j) .eq. BasinPoint)then
!
!                ! [m3water/kgsed] = [m3water/m3cell]*[m3cell] / ([kgsed/m3cell] * [m3cell]) 
!                InstantValue                       = Me%ExtVar%WaterContent(i,j) *  Me%ExtVar%CellVolume(i,j) / &
!                                                    (SoilDryDensity%Concentration(i,j) * Me%ExtVar%CellVolume(i,j))
!
!                ResidualValue                      = Me%DissolvedToParticulate2D(i,j)
!
!                Me%DissolvedToParticulate2D(i,j) = (ResidualValue * Me%ResidualTime +     &
!                                                      InstantValue * DT) / (Me%ResidualTime + DT)
!                                                       
!            end if
!        end do
!        end do
!
!        Me%ResidualTime = Me%ResidualTime + DT
!        
!        
!        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ComputeDissolvedToParticulate3D")
!
!
!    end subroutine ComputeDissolvedToParticulate2D
!    !--------------------------------------------------------------------------

!#ifdef _PHREEQC_
!    !--------------------------------------------------------------------------
!    subroutine SoilChemistryProcesses
!
!        !External--------------------------------------------------------------
!        integer :: STAT_CALL        
!        
!        !Local----------------------------------------------------------------- 
!        type (T_Property), pointer                   :: PropertyX, pH
!        integer                                      :: I, J, K
!        real             , pointer, dimension(:,:,:) :: Theta
!        real(8)          , pointer, dimension(:,:,:) :: Volume
!        
!        !Begin-----------------------------------------------------------------
!        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "SoilChemistryProcesses")
!
!        Theta  => Me%ExtVar%WaterContent
!        Volume => Me%ExtVar%Cellvolume
!        
!        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!        
!            if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then
!                        
!                Me%ExtVar%CellWaterMass(I, J, K) = Theta(I, J, K) * Volume(I, J, K) * WaterReferenceDensity    
!                        
!            end if
!                    
!        end do
!        end do
!        end do
!        
!
!        call SearchProperty(pH, pH_ , .false., STAT = STAT_CALL)    
!            
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property pH not found in porous media properties'
!            stop 'SoilChemistryProcesses - ModuleRunoffProperties - ERR01'
!        endif
!
!        !Question: Is this necessary to any other process than the SoilQualityProcesses?
!        !call ComputeDissolvedToParticulate3D
!        
!        !Question: Is this necessary to any other process than the SoilQualityProcess?
!        !call ComputeWindVelocity
! 
! 
!        !Question: Why here is used SoilChemistry_NextCompute and after is used NextCompute?
!        !          This is related to the possibility of different DT's beetwen 0D model and MOHID general DT?
!        !          
!        if (Me%ExtVar%Now .GE. Me%Coupled%SoilChemistry_NextCompute) then
!            
!            PropertyX => Me%FirstProperty
!
!            do while(associated(PropertyX))
!                
!                call Modify_Interface(InterfaceID   = Me%ObjInterfaceSoilChemistry , &
!                                      PropertyID    = PropertyX%ID%IDNumber        , &
!                                      WaterMass     = Me%ExtVar%CellWaterMass      , &
!                                      Concentration = PropertyX%Concentration      , &
!                                      WaterPoints2D = Me%ExtVar%BasinPoints      , &
!                                      pH            = pH%Concentration             , &
!                                      OpenPoints2D  = Me%ExtVar%BasinPoints       , &
!                                      STAT          = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) &
!                    stop 'SoilChemistryProcesses - ModuleRunoffProperties - ERR02'
!                
!                PropertyX => PropertyX%Next
!                
!            end do
!            
!            Me%Coupled%SoilChemistry_NextCompute = Me%Coupled%SoilChemistry_NextCompute + &
!                                                   Me%Coupled%SoilChemistry_DT
!
!        end if
!
!        PropertyX => Me%FirstProperty
!
!        do while(associated(PropertyX))
!
!            if (PropertyX%Evolution%SoilChemistry) then
!
!                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then
!
!                    call Modify_Interface(InterfaceID   = Me%ObjInterfaceSoilChemistry,                  &
!                                          PropertyID    = PropertyX%ID%IDNumber,            &
!                                          Concentration = PropertyX%Concentration,          &
!                                          WaterPoints2D = Me%ExtVar%BasinPoints,            &
!                                          DTProp        = PropertyX%Evolution%DTInterval,   &
!                                          STAT          = STAT_CALL)
!                    if (STAT_CALL .NE. SUCCESS_)                                            &
!                        stop 'SoilChemistryProcesses - ModuleRunoffProperties - ERR03'
!
!                end if
!
!            end if
!
!            PropertyX => PropertyX%Next
!            
!        end do
!
!        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "SoilChemistryProcesses")
!        
!        !End-------------------------------------------------------------------
!        
!    end subroutine SoilChemistryProcesses   
!    !-----------------------------------------------------------------------------    
!#endif

    !--------------------------------------------------------------------------
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

        !----------------------------------------------------------------------

        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data2D = PropertyX%Concentration,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR01'
                    
                if(PropertyX%Evolution%BottomFluxes) then
                    call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                        Data2D = PropertyX%BottomConcentration,   &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                    &
                        stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR02'                

                    call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                        Data2D = PropertyX%TotalConcentration,    &
                                        STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                    &
                        stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR03'    
                    
                    if(PropertyX%Evolution%Erosion) then

                        call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                            Data2D = PropertyX%ErosionRate,           &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                    &
                            stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR04'    
                    endif
                    
                    if(PropertyX%Evolution%Deposition) then
                        
                        call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                            Data2D = PropertyX%DepositionRate,        &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                    &
                            stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR05'    
                            
                        call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                            Data2D = PropertyX%Ws,                    &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                    &
                            stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR06'    
                            
                    endif
                    
                     if(PropertyX%Evolution%SplashErosion) then
                        
                        call WriteTimeSerie(Me%ObjTimeSerie,                          & 
                                            Data2D = PropertyX%SplashRate,            &
                                            STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                    &
                            stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR07'    
                            
                    endif               
                endif
                    

            endif
            PropertyX=>PropertyX%Next
        enddo

        if (Me%Coupled%BottomFluxes .and. (Me%Coupled%ErosionFluxes .or. Me%Coupled%DepositionFluxes)) then
            call WriteTimeSerie(Me%ObjTimeSerie,                          &
                                Data2D = Me%ShearStress,                  &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                    &
                stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR08'   
        endif
        if (Me%Coupled%BottomFluxes .and. Me%Coupled%SplashErosionFluxes) then                
             call WriteTimeSerie(Me%ObjTimeSerie,                         &
                                Data2D = Me%RainKineticRate,              &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                    &
                stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR09'                   
                
        endif

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
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR00'

                    call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                         Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR01'
           
                    Me%LastOutPutHDF5 = Actual
       
                endif First

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5,                                &
                                      Me%WorkSize%ILB,                           &
                                      Me%WorkSize%IUB,                           &
                                      Me%WorkSize%JLB,                           &
                                      Me%WorkSize%JUB,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR02'


                PropertyX => Me%FirstProperty
                do while (associated(PropertyX))

                    if (PropertyX%OutputHDF) then
 
                        call HDF5WriteData   (Me%ObjHDF5,                                    &
                                              "/Results/"//trim(PropertyX%ID%Name),          &
                                              trim(PropertyX%ID%Name),                       &
                                              trim(PropertyX%ID%Units),                      &
                                              Array2D = PropertyX%Concentration,             &
                                              OutputNumber = OutPutNumber,                   &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR04'
                        
                        if(PropertyX%Evolution%BottomFluxes) then
                    
                            call HDF5WriteData   (Me%ObjHDF5,                                      &
                                                  "/Results/"//trim(PropertyX%ID%Name)//" Bottom", &
                                                  trim(PropertyX%ID%Name)//" Bottom",              &
                                                  "kg/m2",                                         &
                                                  Array2D = PropertyX%BottomConcentration,         &
                                                  OutputNumber = OutPutNumber,                     &
                                                  STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR05'

                            call HDF5WriteData   (Me%ObjHDF5,                                      &
                                                  "/Results/"//trim(PropertyX%ID%Name)//" Total",  &
                                                  trim(PropertyX%ID%Name)//" Total",               &
                                                  "kg/m2",                                         &
                                                  Array2D = PropertyX%TotalConcentration,          &
                                                  OutputNumber = OutPutNumber,                     &
                                                  STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR06'

                            if(PropertyX%Evolution%Erosion) then

                                call HDF5WriteData   (Me%ObjHDF5,                                      &
                                                      "/Results/"//trim(PropertyX%ID%Name)//" Erosion Rate", &
                                                      trim(PropertyX%ID%Name)//" Erosion Rate",       &
                                                      "kg/m2.s",                                       &
                                                      Array2D = PropertyX%ErosionRate,                 &
                                                      OutputNumber = OutPutNumber,                     &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR07'
                            endif

                            if(PropertyX%Evolution%Deposition) then

                                call HDF5WriteData   (Me%ObjHDF5,                                      &
                                                      "/Results/"//trim(PropertyX%ID%Name)//" Deposition Rate", &
                                                      trim(PropertyX%ID%Name)//" Deposition Rate",     &
                                                      "kg/m2.s",                                       &
                                                      Array2D = PropertyX%DepositionRate,              &
                                                      OutputNumber = OutPutNumber,                     &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR07.5'

                                call HDF5WriteData   (Me%ObjHDF5,                                     &
                                                      "/Results/"//trim(PropertyX%ID%Name)//" Fall Velocity",& 
                                                      trim(PropertyX%ID%Name)//" Fall Velocity",      &
                                                      "m/s",                                          &
                                                      Array2D = PropertyX%Ws,                         &
                                                      OutputNumber = OutPutNumber,                    &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR7.6'
                                
                            endif
                    
                           if(PropertyX%Evolution%SplashErosion) then
                           
                                call HDF5WriteData   (Me%ObjHDF5,                                      &
                                                      "/Results/"//trim(PropertyX%ID%Name)//" Splash Rate", &
                                                      trim(PropertyX%ID%Name)//" Splash Rate",         &
                                                      "kg/m2.s",                                       &
                                                      Array2D = PropertyX%SplashRate,                  &
                                                      OutputNumber = OutPutNumber,                     &
                                                      STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR07.5'   

                                                 
                           endif
                        endif     
                    endif
                    
                    PropertyX => PropertyX%Next

                enddo
                if (Me%Coupled%BottomFluxes .and. (Me%Coupled%ErosionFluxes .or. Me%Coupled%DepositionFluxes)) then
                    call HDF5WriteData   (Me%ObjHDF5,                                     &
                                          "/Results/"//"Shear Stress",                    & 
                                          "Shear Stress",                                 &
                                          "N/m2",                                         &
                                          Array2D = Me%ShearStress,                       &
                                          OutputNumber = OutPutNumber,                    &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR08'
                endif

                if (Me%Coupled%BottomFluxes .and. Me%Coupled%SplashErosionFluxes) then
                    call HDF5WriteData   (Me%ObjHDF5,                                     &
                                          "/Results/"//"Total Rain Kinetic Work",         & 
                                          "Total Rain Kinetic Work",                      &
                                          "W/m2",                                         &
                                          Array2D = Me%RainKineticRate,                   &
                                          OutputNumber = OutPutNumber,                    &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR08'
                endif


                if (Me%Coupled%BottomFluxes .and. Me%Coupled%ErosionFluxes) then

                    call HDF5WriteData   (Me%ObjHDF5,                                     &
                                          "/Results/"//"Critical Shear Stress Erosion",   & 
                                          "Critical Shear Stress Erosion",                &
                                          "N/m2",                                         &
                                          Array2D = Me%ErosionCriticalShear%Field,        &
                                          OutputNumber = OutPutNumber,                    &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR09'

                    call HDF5WriteData   (Me%ObjHDF5,                                     &
                                          "/Results/"//"Erosion Coefficient",             & 
                                          "Erosion Coefficient",                          &
                                          "kg.m-2.s-1",                                   &
                                          Array2D = Me%ErosionCoefficient%Field,          &
                                          OutputNumber = OutPutNumber,                    &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR010'
                    
                endif

                if (Me%Coupled%BottomFluxes .and. Me%Coupled%DepositionFluxes) then

                    call HDF5WriteData   (Me%ObjHDF5,                                     &
                                          "/Results/"//"Critical Shear Stress Deposition",& 
                                          "Critical Shear Stress Deposition",             &
                                          "N/m2",                                         &
                                          Array2D = Me%DepositionCriticalShear%Field,     &
                                          OutputNumber = OutPutNumber,                    &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR10.5'

                endif
                
                Me%OutPut%NextOutput = OutPutNumber + 1

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR011'
            
            endif  TOut
        endif  TNum

    end subroutine OutPut_HDF

    !----------------------------------------------------------------------------

    subroutine Output_Boxes_Mass

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL, CHUNK
        type (T_Property), pointer                  :: CurrProperty
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "Output_Boxes")

        CurrProperty => Me%FirstProperty
        do while (associated(CurrProperty)) 
            
            if (CurrProperty%BoxTimeSerie) then
               
                Me%CellMass(:,:) = 0.
                
                CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                !$OMP PARALLEL PRIVATE(I,J)
                
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            
                    if (Me%ExtVar%BasinPoints(i, j) == 1) then
                    
                        !g =  g/m3 * m * m2
                        Me%CellMass(i,j) = CurrProperty%Concentration(i,j)  *  &
                                              Me%ExtVar%WaterColumn(i,j) * Me%ExtVar%Area(i,j)    

                    endif

                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL                  
                
               
                call BoxDif(Me%ObjBoxDif, Me%CellMass,                         &
                            "runoff_"//trim(adjustl(CurrProperty%ID%Name)),    &
                            Me%ExtVar%BasinPoints,                             &
                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  &
                    stop 'Output_Boxes - ModuleRunoffProperties - ERR01'        
                
                Me%CellMass(:,:) = null_real                                    

            endif
            
            CurrProperty => CurrProperty%Next
        
        end do 
        
        nullify (CurrProperty)
       
        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "Output_Boxes")    
    
    end subroutine Output_Boxes_Mass

    !----------------------------------------------------------------------------
    
    subroutine Output_Boxes_Fluxes (Property)
    
        !Arguments---------------------------------------------------------------
        type (T_Property)                         :: Property
        !Local-------------------------------------------------------------------
        integer                                   :: CHUNK, STAT_CALL, i, j
        !Begin-------------------------------------------------------------------

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
                                
        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do J = Me%WorkSize%JLB, Me%WorkSize%JUB
do4 :   do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            Me%Fluxes%MassFluxesX (I,J) = Me%Fluxes%AdvFluxX(I,J) + Me%Fluxes%DifFluxX (I,J)
            Me%Fluxes%MassFluxesY (I,J) = Me%Fluxes%AdvFluxY(I,J) + Me%Fluxes%DifFluxY (I,J)
        end do do4
        end do do3
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        !Integration of fluxes
        call BoxDif(Me%ObjBoxDif,                        &
                    Me%Fluxes%MassFluxesX,               &
                    Me%Fluxes%MassFluxesY,               &
                    "runoff_"//trim(Property%ID%Name),   &
                    Me%ExtVar%BasinPoints,               &
                    STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                     &
        stop 'Output_Boxes_Fluxes - ModuleRunoffProperties - ERR300'
                
    end subroutine Output_Boxes_Fluxes
    
    !--------------------------------------------------------------------------

    subroutine WriteFinalFile

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        !integer                                     :: OutPutNumber
        integer                                     :: HDF5_CREATE
        character(LEN = PathLength)                 :: FileName
        integer                                     :: ObjHDF5
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
        type (T_Time)                               :: Actual           
        real                                        :: Total_Mass_Created
        character (Len = StringLength)              :: str_mass_created, string_to_be_written         
        !Begin----------------------------------------------------------------

        !Gets a pointer to Topography
        call GetGridData        (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR00'

        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR01'

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
            stop 'WriteFinalFile - ModuleRunoffProperties - ERR10'

        Actual   = Me%ExtVar%Now
         
        call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
        !Writes Time
        TimePtr => AuxTime
        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR11'

        call HDF5WriteData  (ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                             Array1D = TimePtr, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR12'


        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5,                                &
                              Me%WorkSize%ILB,                           &
                              Me%WorkSize%IUB,                           &
                              Me%WorkSize%JLB,                           &
                              Me%WorkSize%JUB,                           &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffproperties - ERR02'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR25'

        !Writes the Grid
        call HDF5WriteData   (ObjHDF5, "//Grid/Topography", "Topography", "m",           &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR05'

        !WriteBasinPoints
        call HDF5WriteData   (ObjHDF5, "//Grid/BasinPoints", "BasinPoints", "-",          &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR07'


        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))

            call HDF5SetLimits   (ObjHDF5,                                &
                                  Me%WorkSize%ILB,                           &
                                  Me%WorkSize%IUB,                           &
                                  Me%WorkSize%JLB,                           &
                                  Me%WorkSize%JUB,                           &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffproperties - ERR10'

            call HDF5WriteData   (ObjHDF5,                                    &
                                  "/Results/"//trim(PropertyX%ID%Name),          &
                                  trim(PropertyX%ID%Name),                       &
                                  trim(PropertyX%ID%Units),                      &
                                  Array2D = PropertyX%Concentration,             &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffproperties - ERR14'
                
            if(PropertyX%Evolution%BottomFluxes) then
        
                call HDF5WriteData   (ObjHDF5,                                      &
                                      "/Results/"//trim(PropertyX%ID%Name)//" Bottom", &
                                      trim(PropertyX%ID%Name)//" Bottom",              &
                                      "kg/m2",                                         &
                                      Array2D = PropertyX%BottomConcentration,         &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffproperties - ERR15'

                   
            endif

            if (PropertyX%Evolution%MinConcentration .and. Me%ExtVar%Now == Me%ExtVar%EndTime) then
  
                call HDF5WriteData   (ObjHDF5,                                        &
                                      "/Results/"//trim(PropertyX%ID%Name)//" Mass Created",& 
                                      "Property Mass Created",                        &
                                      "g",                                            &
                                      Array2D = PropertyX%Mass_Created,               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR10.5'

                !g/1000 = kg, to avoid big numbers
                Total_Mass_Created = SUM(PropertyX%Mass_Created)/1000

                write(str_mass_created, '(f20.8)') Total_Mass_Created
      
                string_to_be_written = 'Due to MinConcentration RP Total mass (kg) created on property ' //&
                                        trim(adjustl(adjustr(PropertyX%ID%Name)))//' = ' //&
                                        trim(adjustl(adjustr(str_mass_created))) 
            
                !Writes total mass created to "Error_and_Messages.log" file
                call SetError(WARNING_, INTERNAL_, string_to_be_written, OFF)

            endif
            
            PropertyX => PropertyX%Next

        enddo

        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffproperties - ERR030'

        !Unget
        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR90'  

        !UnGets Topography
        call UnGetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleRunoffProperties - ERR100'
            

    end subroutine WriteFinalFile

    !----------------------------------------------------------------------------
    
    subroutine CalculateTotalStoredMass

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        type (T_Property), pointer                  :: CurrProperty
        real                                        :: BottomMass
        !Begin-----------------------------------------------------------------
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModuleRunoffProperties - ERR01'        

        call GetRunoffWaterColumn     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModuleRunoffProperties - ERR10'

        call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                 GridCellArea = Me%ExtVar%Area,                            & 
                                 STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModuleRunoffProperties - ERR020'
        
        CurrProperty => Me%FirstProperty
        
        do while (associated(CurrProperty)) 
            
            !This cycle is not paralelizable because it changes the same address in memory
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        
                if (Me%ExtVar%BasinPoints(i, j) == 1) then
                    BottomMass = 0.0
                    if (CurrProperty%Evolution%BottomFluxes) then
                        !kg = kg/m2 * m2
                        BottomMass = CurrProperty%BottomConcentration(i,j) * Me%ExtVar%Area(i,j)
                    endif
                    !kg = kg + kg + [m * m2] * [g/m3] * [1E-3 kg/g])
                    CurrProperty%MB%TotalStoredMass = CurrProperty%MB%TotalStoredMass  + Bottommass           &
                                                      + (Me%ExtVar%WaterColumn(i,j)* Me%ExtVar%Area(i,j)      &
                                                       * CurrProperty%Concentration(i,j) * 1E-3    )
                        
                endif

            enddo
            enddo
            
            CurrProperty => CurrProperty%Next
        end do 

        call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModuleRunoffProperties - ERR030'

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModuleRunoffProperties - ERR040'

        call UngetBasin           (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalculateTotalStoredMass - ModuleRunoffProperties - ERR050'        
        

    end subroutine CalculateTotalStoredMass

    !----------------------------------------------------------------------------

    subroutine LogCoefs(i,j,Coef1, Coef2, Coef3, Coef4, Coef5)

        !Arguments-------------------------------------------------------------
        integer                                     :: i,j
        real(8)                                     :: Coef1, Coef2, Coef3, Coef4, Coef5

        !Local-----------------------------------------------------------------
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second

        call ExtractDate(Me%ExtVar%Now, Year, Month, Day, Hour, Minute, Second)
        
        write (Me%Files%AsciiUnit, fmt=1000) Year, Month, Day, Hour, Minute, Second,    &
                                             i,j, Coef1, Coef2, Coef3, Coef4, Coef5

        1000 format(f5.0, f5.0, f5.0, f5.0, f5.0, f12.5, i3, i3, 5f13.8)

    end subroutine LogCoefs

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillRunoffProperties(ObjRunoffPropertiesID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjRunoffPropertiesID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers,STAT_CALL  
        type(T_property), pointer           :: PropertyX
        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mRunoffProperties_,  Me%InstanceID)
            
            call WriteFinalFile

            PropertyX => Me%FirstProperty
            
            do while (associated(PropertyX)) 
                if(PropertyX%ID%SolutionFromFile)then

                    call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillRunoffProperties - ModuleRunoffProperties - ERR00'
                end if
                
                PropertyX => PropertyX%Next
            end do 


            if (Me%Coupled%BottomFluxes .and. Me%Coupled%ErosionFluxes) then
             
                if (Me%ErosionCriticalShear%ID%SolutionFromFile) then            

                    call KillFillMatrix(Me%ErosionCriticalShear%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillRunoffProperties - ModuleRunoffProperties - ERR01'
                endif


                if (Me%ErosionCriticalShear%ID%SolutionFromFile) then            

                    call KillFillMatrix(Me%ErosionCoefficient%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillRunoffProperties - ModuleRunoffProperties - ERR02'
                endif 
                
            endif           

            if (Me%Output%Boxes_ON) then
                call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillRunOff - RunOff - ERR04'
            endif            
            
            if (associated(Me%Disper_Longi%Field))then
                deallocate(Me%Disper_Longi%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillRunoffProperties - ModuleRunoffProperties - ERR03'
                nullify   (Me%Disper_Longi%Field)
            end if

            if (associated(Me%Disper_Trans%Field))then
                deallocate(Me%Disper_Trans%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillRunoffProperties - ModuleRunoffProperties - ERR04'
                nullify   (Me%Disper_Trans%Field)
            end if
            
            if (nUsers == 0) then

                !Kills the TimeSerie
                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillRunoff - ModuleRunoffProperties - ERR05'
                endif

                
                if (Me%OutPut%HDF_ON) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillVegetation - ModuleRunoffProperties  - ERR08'
                endif
                
                !Runoff constructed, not RunoffProperties
!                if (Me%ComputeOptions%Discharges) then
!                    call Kill_Discharges(Me%ObjDischarges, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'KillRunOff - ModuleRunoffProperties - ERR05a'
!                endif  
                if (Me%ComputeOptions%Discharges) then
                    nUsers = DeassociateInstance(mDISCHARGES_,    Me%ObjDischarges)
                    if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR06'
                endif
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR07'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR08'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR11'
                
                nUsers = DeassociateInstance (mRunoff_,  Me%ObjRunoff)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR12'

!                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
!                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR13'

                nUsers = DeassociateInstance (mGRIDDATA_,  Me%ObjGridData)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR13'

                
                call DeallocateVariables

                !Deallocates Instance
                call DeallocateInstance ()

                ObjRunoffPropertiesID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillRunoffProperties
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_RunoffProperties), pointer          :: AuxObjRunoffProperties
        type (T_RunoffProperties), pointer          :: PreviousObjRunoffProp

        !Updates pointers
        if (Me%InstanceID == FirstObjRunoffProperties%InstanceID) then
            FirstObjRunoffProperties => FirstObjRunoffProperties%Next
        else
            PreviousObjRunoffProp => FirstObjRunoffProperties
            AuxObjRunoffProperties      => FirstObjRunoffProperties%Next
            do while (AuxObjRunoffProperties%InstanceID /= Me%InstanceID)
                PreviousObjRunoffProp => AuxObjRunoffProperties
                AuxObjRunoffProperties      => AuxObjRunoffProperties%Next
            enddo

            !Now update linked list
            PreviousObjRunoffProp%Next => AuxObjRunoffProperties%Next

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
        
!#ifdef _PHREEQC_
!        deallocate (Me%ExtVar%CellWaterMass)
!#endif        

        if (Me%Coupled%AdvectionDiffusion) then
            deallocate (Me%WaterVolume)
            deallocate(Me%TICOEF3)

            deallocate(Me%COEF3_HorAdvXX%C_Flux)
            deallocate(Me%COEF3_HorAdvXX%D_Flux)
            deallocate(Me%COEF3_HorAdvXX%E_Flux)
            deallocate(Me%COEF3_HorAdvXX%F_Flux)
            
            deallocate(Me%COEF3_HorAdvYY%C_Flux)
            deallocate(Me%COEF3_HorAdvYY%D_Flux)
            deallocate(Me%COEF3_HorAdvYY%E_Flux)
            deallocate(Me%COEF3_HorAdvYY%F_Flux)
                
            if (.not. Me%ComputeOptions%AdvDiff_Explicit) then
                deallocate(Me%COEF3%D)
                deallocate(Me%COEF3%E)
                deallocate(Me%COEF3%F)
                deallocate(Me%VECG)
                deallocate(Me%VECW)
              
            endif
        endif
      
   
    end subroutine DeallocateVariables 

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjRunoffProperties_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunoffProperties_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjRunoffProperties_ID > 0) then
            call LocateObjRunoffProperties (ObjRunoffProperties_ID)
            ready_ = VerifyReadLock (mRunoffProperties_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjRunoffProperties (ObjRunoffPropertiesID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunoffPropertiesID

        !Local-----------------------------------------------------------------

        Me => FirstObjRunoffProperties
        do while (associated (Me))
            if (Me%InstanceID == ObjRunoffPropertiesID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleRunoffProperties - LocateObjRunoffProperties - ERR01'

    end subroutine LocateObjRunoffProperties

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
                if (PrintWarning) write (*,*)'Property Not Found in Module RunoffProperties ', &
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
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR01'        
        
        call GetOverLandFlow      (Me%ObjRunoff, Me%ExtVar%FluxU, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR030'

        call GetRunoffCenterVelocity (Me%ObjRunoff, Me%ExtVar%CenterVelU, Me%ExtVar%CenterVelV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR040'

        call GetRunoffWaterColumn     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR64'

        !Water Column After Transport (Advection Diffusion only includes transport)
        call GetRunoffWaterColumnAT     (Me%ObjRunoff, Me%ExtVar%WaterColumnAT, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR64'

        call GetRunoffWaterColumnOld   (Me%ObjRunoff, Me%ExtVar%WaterColumnOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR65'
        

        call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                 GridCellArea = Me%ExtVar%Area,                            & 
                                 STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR070'


        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                    &
                                  DVY         = Me%ExtVar%DVY,                          &
                                  DUY         = Me%ExtVar%DUX,                          &
                                  DXX         = Me%ExtVar%DXX,                          &
                                  DYY         = Me%ExtVar%DYY,                          &
                                  DZY         = Me%ExtVar%DZY,                          &
                                  DZX         = Me%ExtVar%DZX,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR101'

        if (Me%ExtVar%CoupledDN) then

            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR01'

            call GetFlowToChannels   (Me%ObjRunoff, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR10'            
            
        endif    

        if (Me%ExtVar%BoundaryImposed) then
            call GetBoundaryFlux (Me%ObjRunoff, Me%ExtVar%BoundaryFlux, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR170'
        endif

        if (Me%ExtVar%RouteDFour) then
            call GetRouteDFourFlux (Me%ObjRunoff, Me%ExtVar%RouteDFourFlux, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR180'
        endif

    end subroutine ReadLockExternalVar

    !-----------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        call UngetBasin           (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR01'        

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%FluxU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR030'

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR040'
        
        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%CenterVelU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR050'    
        
        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%CenterVelV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR060'             

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR061'
        
        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%WaterColumnAT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR063'

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%WaterColumnOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR064'

        call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleRunoffProperties - ERR070'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR080'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR090'        


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR112'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR113'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DXX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR114'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR115'

        if (Me%ExtVar%CoupledDN) then

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR140'               
            
            call UnGetRunoff (Me%ObjRunoff, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR150'     
  
            
        endif
        
        if (Me%ExtVar%BoundaryImposed) then
            call UnGetRunoff (Me%ObjRunoff, Me%ExtVar%BoundaryFlux, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0190'      
        endif

        if (Me%ExtVar%RouteDFour) then
            call UnGetRunoff (Me%ObjRunoff, Me%ExtVar%RouteDFourFlux, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0200'      
        endif

    endsubroutine ReadUnlockExternalVar
    
    !--------------------------------------------------------------------------

end module ModuleRunoffProperties

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 
