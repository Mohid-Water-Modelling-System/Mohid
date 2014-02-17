    !-----------------------------------------------------------------------------------
    !        IST/MARETEC, Water Modelling Group, Mohid modelling system 
    !-----------------------------------------------------------------------------------
    !
    ! TITLE         : Mohid Model
    ! PROJECT       : Mohid Base 1
    ! MODULE        : Bivalve
    ! URL           : http://www.mohid.com
    ! AFFILIATION   : IST/MARETEC, Marine Modelling Group
    ! DATE          : 30 Set 2013
    ! REVISION      : Sofia Saraiva
    ! DESCRIPTION   : Individual Based Population (or individual) model for one/several bivalve
    !                 species following DED theory 
    !
    !-----------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    !
    !DataFile example
    !
    !DT                   : 1800.
    !DENSITY_UNITS        : 0 !0:m2, 1:m3
    !BIVALVE_OUTPUT_TIME  : 0 1800.
    !PELAGIC_MODEL        : WaterQuality
    !NITROGEN             : 1
    !PHOSPHOR             : 1
    !SIMPLE_FILTRATION    : 0
    !CORRECT_FILTRATION   : 1
    !INDEX_OUTPUT         : 411 ! cell for time serie results
    !MASS_BALANCE         : 1 ! only works with no predation + the same NC ratio food/bivalve (or with life) + complex filtrat
    !MIN_NUMBER           : 1 ! minimum number of organism in a cohort
    !TESTING_PARAMETERS   : 0 ! 0/1, the name of the output file will incude the parameters
    !
    !
    !<begin_species>
    !NAME                 : bivalve1
    !DESCRIPTION          : Mytilus edulis
    !TESTING_FILENAME     : All_ParameterTest.dat
    !POPULATION           : 1
    !FEED_ON_LARVAE       : 0
    !LARVAE_MAXSIZE       : 0.026 !cm
    !NUMBER_OF_COHORTS    : 1
    !MIN_OBS_LENGTH       : 0.1
    !COHORT_OUTPUT        : 1
    !BYSIZE_OUTPUT        : 1
    !!SIZE_STEP            : 1        !1cm classes from 0 to the maximum size of the species
    !!MAX_SIZECLASS        : 10       !maximum size for size distribution representation
    !<<begin_size_classes>>
    !0.0
    !0.026
    !0.03
    !0.05
    !0.1
    !0.15
    !0.2
    !0.5
    !1.2
    !2.6
    !3.
    !4.
    !6.
    !8.
    !10.
    !12.
    !14.
    !16.
    !<<end_size_classes>>
    !
    !RESERVES_nH          : 1.8      !molH/molC, chemical index of hydrogen in bivalve reserves (Kooijman 2010)
    !RESERVES_nO          : 0.53     !molO/molC, chemical index of oxygen in bivalve reserves (Kooijman 2010)
    !RESERVES_nN          : 0.18     !molN/molC, chemical index of nitrogen in bivalve reserves
    !RESERVES_nP          : 0.006    !molP/molC, chemical index of phosphorus in bivalve reserves
    !
    !STRUCTURE_nH         : 1.8      !molH/molC, chemical index of hydrogen in bivalve structure
    !STRUCTURE_nO         : 0.53     !molO/molC, chemical index of oxygen in bivalve structure
    !STRUCTURE_nN         : 0.18     !molN/molC, chemical index of nitrogen in bivalve structure
    !STRUCTURE_nP         : 0.006    !molP/molC, chemical index of phosphorus in bivalve structure
    !
    !Tref                 : 293      !293, K, Rate Temperature reference (van der Veer etal., 2006)
    !TA                   : 7022     !7022, K, Arrhenius temperature (van der Veer etal., 2006)
    !TL                   : 275      !273, K, Lower Boundary tolerance range (van der Veer etal., 2006)
    !TH                   : 296      !290, K, Upper Boundary tolerance range (van der Veer etal., 2006)
    !TAL                  : 45430    !45430, K, Arrhenius temperature for lower boundary (van der Veer etal., 2006)
    !TAH                  : 31376    !31376, K, Arrhenius temperature for upper boundary (van der Veer etal., 2006)
    !
    !F_FIX                : 1        !1, adim, constant food density parameter (only if simple filtration)
    !PAM_FIX              : 94.79    !80.5, Jd-1cm-2, bivalve surface-specific assimilation rate if fix (Saraiva etal., inpress)
    !DELTA_M              : 0.297    !0.297, cm(volumetric)/cm(real), shape coefficient  (Saraiva etal., inpress)
    !LIFE_SPAN            : 24       !24, years, max life span for a mussel under natural conditions (Sukhotin et al. (2007))
    !M_VELOCITY           : 0.0      !/d, fraction of individuals that die due to high velocity
    !MAX_VELOCITY         : 0.5      !m/s, maximum  water velocity tolerable for this species
    !M_NATURAL            : 0.50000E-02
    !M_SPAT               : 0.9
    !DENSITYLIMIT         : 1        ! density limitation?
    !DENSITY_MAXVALUE     : 3000     !3000 #/m2, maxium density found in field observations       
    !V_COND               : 0.056    !0.056, cm/d, energy conductance (Saraiva etal., in press)
    !KAPPA                : 0.67     !0.67, adim, allocation fraction to growth/somatic maintenace (Saraiva etal., in press)
    !KAP_R                : 0.95     !0.95, adim, fraction of flux allocated to reproduction (Kooijman, 2010)
    !pM                   : 11.6     !11.6, J/(d.cm3), volume specific somatic maintenace energy flux (Saraiva etal., inpress)
    !EG                   : 5993     !5993, J/cm3(volumetric), energy costs for structural volume growth (Saraiva etal., inpress)
    !EH_B                 : 2.95e-5  !2.99e-5, J, Maturity threshold for birth (Saraiva etal., inpress)
    !EH_P                 : 1.58e2   !1.58e2, J, Maturity threshold for puberty (Saraiva etal., inpress)
    !CRM                  : 0.096    !0.096, m3/d.cm2, maximum clearance rate (Saraiva etal., 2011)
    !JX1FM                : 4.8e-4   !4.8e-4, molC/(d.cm2), algae maximum surface area-specific filtration rate (Thomas etal., 2011)
    !JX0FM                : 3.5      !3.5, g/(d.cm2), inorganic material maximum surface area-specific filtration rate (Saraiva etal.
    !RO_X1                : 0.4      !0.4, adim, algae binding probability (Saraiva etal., inpress)
    !RO_X0                : 0.4      !0.4, adim, inorganic material binding probability (Saraiva etal., inpress)
    !JX1IM                : 1.3e4    !1.3e-4, molC/(d.cm2), algae maximum surface area-specific ingestion rate (Saraiva etal., 2011)
    !JX0IM                : 0.11     !0.11, g/(d.cm2), inorganic material maximum surface area-specific ingestion rate (Saraiva etal.
    !YEX                  : 0.75     !0.65, molCE/molCV, yield coeficienct of reserves in algae structure
    !GSR_MIN              : 0.1      !0.1, molC(gam)/molC(struc), minimum gonado-somatic ratio in the organism (Cardoso et al., 2007)
    !GSR_SPAWN            : 0.2      !0.2, molC(gam)/molC(struc), gonado-somatic ratio to spawn (Saraiva etal., submited)
    !T_SPAWN              : 9.6      !9.6, C, minimum temperature for spawning (Hummel etal., 1989)
    !MIN_SPAWN_TIME       : 0        !15, d, minimum time between spawning events
    !ME_0                 : 1.49e-10 !1.48e-10, molC(reser), reserves in an embryo at optimal food conditions (Saraiva etal., submite
    !ME_B                 : 6.0e-11  !1.0e-7, molC, reserves in a new born individual at optimal food conditions (Saraiva etal., subm
    !MV_B                 : 7.92e-11 !7.52e-11, molC, structure in a new born individual at optimal food conditions (Saraiva etal., s
    !MH_B                 : 4.24e-11 !4.24e-11, molC, maturity in a new born individual at optimal food conditions (Saraiva etal., su
    !L_B                  : 7.3e-3   !7.3e-3, molC, length in a new born individual at optimal food conditions (Saraiva etal., submit
    !DV                   : 0.2      !0.2, g(dw)/cm3, bivalve structure and reserves specific density (Rosland etal., 2009 and Brey,
    !MU_E                 : 6.97e5   !6.97e5, J/molC(reser), chemical potential of reserves (van der Veer etal., 2006)
    !SIMPLE_ASSI          : 0        !0, 0/1 option to compute simple assimilation
    !SIMPLE_TEMP          : 1        !0, 0/1 option to compute simple temperature
    !
    !<<begin_particle>>
    !NAME                 : phytoplankton
    !DESCRIPTION          : description
    !ORGANIC              : 1       !1/0, is this an organic particle?
    !SILICA_USE           : 0       !1/0, does it have silica?, just important in case of life model (diatoms)
    !RATIO_VARIABLE       : 0       !1/0, the N/C and P/C ratios of this property are variable?
    !RATIOHC              : 0.15    !0.15, mgH/mgC in the food, same as bivalve
    !RATIOOC              : 0.71    !0.71, mgO/mgC in the food, same as bivalve
    !RATIONC              : 0.3396     !0.3, 0.18, Redfield, to define if RATIO_VARIABLE:0 and FOOD:1 and ORGANIC:1, must be equal to wq
    !RATIOPC              : 0.07    !0.07, 0.024, Redfield, to define if RATIO_VARIABLE:0 and FOOD:1 and ORGANIC:1, must be equal to wq
    !RATIOSiC             : 0.89    !0.89, to define if RATIO_VARIABLE:0 and FOOD:1 and ORGANIC:1, must be equal to wq
    !RATIOCHLC            : 0.017   !0.017, to define if RATIO_VARIABLE:0 and FOOD:1 and ORGANIC:1, must be equal to wq
    !SIZE                 : 0.2     !0.2, to define if FOOD:1, cm, mean size of property
    !F_E                  : 0.8     !0.5, to define if FOOD:1, molCReserves/molCTotalFood, fraction of reserves in the food
    !<<end_particle>>
    !
    !<<begin_predator>>
    !NAME                 : shrimp
    !DESCRIPTION          : shrimp
    !SIZE                 : 1.2     !0.0, cm, shrimp average size (Beukema 1992)
    !MINPREYSIZE          : 0.00   !0.0, cm, minimum size of the prey
    !MAXPREYSIZE          : 0.2    !0.2, cm, maximum size of the prey, 10%
    !FEEDING_RATE         : 48.35   !48.35; %J/cm2.d shrimp average ingestion rate (Campos etal., 2009)
    !FEEDING_UNITS        : 3       !1-#/d.ind; 2-AFDW/d.ind, 3-J/cm2.d
    !FEEDING_TIME         : 1       !1-Always; 2-LowTide, 3-HighTide
    !DIET_FRACTION        : 0.00000E+00
    !AFDW_DW              : 0.85    !-, conversion of afdw in dw (this study)
    !DW_C                 : 2.5     !-, conversion of dw in gCarbon (Solobodkin and Richman 1961)
    !
    !CORRECT_TEMP         : 1        !0, 0/1 option to compute simple temperature for predator rate
    !SIMPLE_TEMP          : 1        !0, 0/1 option to compute simple temperature for predator rate
    !P_Tref               : 283.     !293, K, Rate Temperature reference for predator (same as the prey...for now)
    !P_TA                 : 9000     !7022, K, Arrhenius temperature for predator
    !P_TL                 : 273      !273, K, Lower Boundary tolerance range for predator
    !P_TH                 : 303      !290, K, Upper Boundary tolerance range for predator
    !P_TAL                : 6700000  !45430, K, Arrhenius temperature for lower boundary for predator
    !P_TAH                : 49368    !31376, K, Arrhenius temperature for upper boundary for predator
    !<<end_predator>>
    !<end_species>

    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------

    Module ModuleBivalve

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleFunctions, only: Chunk_I
    use ModuleTime
    !use ifport

    implicit none

    private 

    !Subroutines---------------------------------------------------------------------

    !Constructor
    public  ::      ConstructBivalve
    private ::          AllocateInstance
    private ::          ReadDataBivalve
    private ::              ConstructGlobalVariables
    private ::              ConstructSpecies
    private ::                  AddSpecies
    private ::                  ConstructSpeciesParameters
    private ::                      ConstructSpeciesComposition
    private ::                      ConstructIndividualParameters
    private ::                      ConstructParticles
    private ::                          AddParticles
    private ::                          ConstructParticlesParameters
    private ::                      ConstructPredator
    private ::                          AddPredator
    private ::                          ConstructPredatorParameters
    private ::                  ConstructCohort
    private ::                      AddCohort
    private ::                      ConstructCohortOutput
    private ::              ConstructOutputs
    private ::          PropertyIndexNumber
    private ::          ConstructPropertyList
    private ::          KillEnterData    

    !Selector
    public  ::      GetBivalvePropertyList
    public  ::      GetDTBivalve
    public  ::      GetBivalveSize
    public  ::      GetBivalvePropIndex
    public  ::      GetBivalveListDeadIDS
    public  ::      GetBivalveNewborns
    public  ::      GetBivalveNewBornParameters
    public  ::      GetBivalveOtherParameters
    public  ::      UnGetBivalve
    public  ::      SearchPropIndex
    public  ::      SetSettlementOnBivalve

    !Modifier
    public  ::      ModifyBivalve
    public  ::          AllocateAndInitializeByTimeStep
    private ::          ComputeBivalve
    private ::              PrepareRunByIndex
    private ::                  ConvertUnits
    private ::                  RestoreParticlesList
    private ::                  ConstructParticlesFromLarvae
    private ::              ComputeMortalityByVelocity
    private ::              ComputeMortalityByWrongSettlement
    private ::              ComputeIndividualProcesses
    private ::                  ComputeChemicalIndices
    private ::                  ComputeAuxiliarParameters
    private ::                  ComputeBivalveCondition
    private ::                  ComputeFeedingProcesses
    private ::                      ComputeSimpleFiltration
    private ::                      ComputeComplexFiltration
    private ::                          ComputeClearanceRate
    private ::                          ComputeFiltrationRate
    private ::                          ComputeIngestionAssimilationRate
    private ::                          UpdateLarvaeOnMatrixMass
    private ::                          UpdateMatrixMass
    private ::                  ComputeSomaticMaintenance
    private ::                  ComputeMobilization
    private ::                  ComputeReservesDynamics
    private ::                  ComputeStructureDynamics
    private ::                  ComputeMaturity
    private ::                  ComputeLengthAgeDynamics
    private ::                  ComputeSpawning
    private ::                  ComputeReproductionDynamics
    private ::                  ComputeInorganicFluxes
    private ::                  ImposeCohortDeath
    private ::              ComputeExtraStarvationMortality
    private ::              ComputeNaturalMortality
    private ::              ComputePredation
    private ::                  ComputePredationByPredator
    private ::              ComputePopulationVariables
    private ::              WriteOutput
    private ::                  WriteSizeDistribution
    private ::                  WriteMassBalance
    private ::              UpdateCohortState
    private ::              RestoreUnits
    private ::          UpdateListDeadAndNewBornIDs
    public  ::          UpdateBivalvePropertyList
    !Destructor
    public  ::      KillBivalve
    private ::          WriteTestingFile           
    private ::      DeAllocateInstance
    private ::      CloseFiles
    !Management
    private ::      Ready
    private ::      LocateObjBivalve 

    !Interfaces---------------------------------------------------------------------
    private ::       UnGetBivalve1D_I
    private ::       UnGetBivalve2D_R4
    interface        UnGetBivalve
    module procedure UnGetBivalve1D_I
    module procedure UnGetBivalve2D_R4
    end interface    UnGetBivalve

    !Types--------------------------------------------------------------------------
    type     T_ExternalVar
        type(T_Time)                     :: CurrentTime
        real,    pointer, dimension(:  ) :: Salinity       => null()
        real,    pointer, dimension(:  ) :: Temperature    => null()
        real,    pointer, dimension(:,:) :: Mass           => null()    !g/m3
        integer, pointer, dimension(:  ) :: OpenPoints     => null()
        real(8), pointer, dimension(:  ) :: CellVolume     => null()    !m3
        real,    pointer, dimension(:  ) :: CellArea       => null()    !m2
        real,    pointer, dimension(:  ) :: Velocity       => null()
        real,    pointer, dimension(:  ) :: InitialPhyto   => null()    
        real,    pointer, dimension(:  ) :: InitialShrimp  => null()     
    end type T_ExternalVar

    type     T_ComputeOptions
        logical                          :: Nitrogen           = .false.
        logical                          :: Phosphorus         = .false.
        logical                          :: SimpleFiltration   = .false.
        logical                          :: CorrectFiltration  = .true.
        logical                          :: MassBalance        = .false.
        character(len=StringLength)      :: PelagicModel       = null_str
    end type T_ComputeOptions 


    type      T_PropIndex     
        !Sediments             
        integer                          :: sediments             = null_int   
        !Nitrogen                            
        integer                          :: AM                    = null_int
        integer                          :: PON                   = null_int   
        !Phosphorus                          
        integer                          :: IP                    = null_int   
        integer                          :: POP                   = null_int   
        !Carbon 
        integer                          :: POC                   = null_int   
        integer                          :: CarbonDioxide         = null_int   
        !Oxygen 
        integer                          :: Oxygen                = null_int
        !Silica                            
        integer                          :: BioSilica             = null_int
        integer                          :: DissSilica            = null_int
        !Particles  
        integer                          :: phyto                 = null_int
        integer                          :: diatoms               = null_int
        integer                          :: zoo                   = null_int
        integer                          :: ciliate               = null_int
        integer                          :: bacteria              = null_int
        integer                          :: silica                = null_int
        integer                          :: DiatomsC              = null_int
        integer                          :: DiatomsN              = null_int
        integer                          :: DiatomsP              = null_int
        integer                          :: DiatomsChl            = null_int
        integer                          :: DiatomsSi             = null_int
        integer                          :: Mix_FlagellateC       = null_int
        integer                          :: Mix_FlagellateN       = null_int
        integer                          :: Mix_FlagellateP       = null_int
        integer                          :: Mix_FlagellateChl     = null_int
        integer                          :: PicoalgaeC            = null_int
        integer                          :: PicoalgaeN            = null_int
        integer                          :: PicoalgaeP            = null_int
        integer                          :: PicoalgaeChl          = null_int
        integer                          :: FlagellateC           = null_int
        integer                          :: FlagellateN           = null_int
        integer                          :: FlagellateP           = null_int
        integer                          :: FlagellateChl         = null_int
        integer                          :: MicrozooplanktonC     = null_int
        integer                          :: MicrozooplanktonN     = null_int
        integer                          :: MicrozooplanktonP     = null_int
        integer                          :: Het_NanoflagellateC   = null_int
        integer                          :: Het_NanoflagellateN   = null_int
        integer                          :: Het_NanoflagellateP   = null_int
        integer                          :: MesozooplanktonC      = null_int
        integer                          :: MesozooplanktonN      = null_int
        integer                          :: MesozooplanktonP      = null_int
        integer                          :: Het_BacteriaC         = null_int
        integer                          :: Het_BacteriaN         = null_int
        integer                          :: Het_BacteriaP         = null_int
        !Predators  
        integer                          :: Shrimp                = null_int
        integer                          :: Crab                  = null_int
        integer                          :: OysterCatcher         = null_int
        integer                          :: EiderDuck             = null_int
        integer                          :: HerringGull           = null_int
    end type T_PropIndex  

    type     T_ID         
        integer                          :: ID, IDNumber    = null_int
        character(len=StringLength)      :: Name            = null_str
        character(len=StringLength)      :: Description     = null_str
    end type T_ID         

    type     T_Composition
        real                             :: nC   = null_real   !molC/molC, chemical index of carbon = 1
        real                             :: nH   = null_real   !molH/molC, chemical index of hydrogen
        real                             :: nO   = null_real   !molO/molC, chemical index of oxygen
        real                             :: nN   = null_real   !molN/molC, chemical index of nitrogen
        real                             :: nP   = null_real   !molP/molC, chemical index of phosphorus
    end type T_Composition          

    type     T_Ratios  
        real                             :: HC_Ratio   = null_real !Actual Ratio H/C [mgH/mgC]
        real                             :: OC_Ratio   = null_real !Actual Ratio O/C [mgO/mgC]
        real                             :: NC_Ratio   = null_real !Actual Ratio N/C [mgN/mgC]
        real                             :: PC_Ratio   = null_real !Actual Ratio P/C [mgP/mgC]
        real                             :: ChlC_Ratio = null_real !Actual Ratio P/C [mgP/mgC]
        real                             :: SiC_Ratio  = null_real !Actual Ratio P/C [mgP/mgC]
    end type T_Ratios 

    type     T_SpeciesComposition       
        type (T_Composition  )           :: ReservesComposition  
        type (T_Composition  )           :: StructureComposition
    end type T_SpeciesComposition

    type     T_StateIndex      
        integer                          :: M_V      = null_int !molC (struc), structure biomass
        integer                          :: M_E      = null_int !molC (reser), reserves biomass
        integer                          :: M_H      = null_int !molC (reser), maturity
        integer                          :: M_R      = null_int !molC (reser), reproduction buffer
        integer                          :: L        = null_int !cm, bivalve real length
        integer                          :: Age      = null_int !days, bivalve age
        integer                          :: Number   = null_int !#, number of organims in the cohort
    end type T_StateIndex 

    type    T_IndividualParameters
        real                             :: Tref               = null_real !K, Rate Temperature reference
        real                             :: TA                 = null_real !K, Arrhenius Temp
        real                             :: TL                 = null_real !K, Lower Boundary tolerance range
        real                             :: TH                 = null_real !K, Upper Boundary tolerance range
        real                             :: TAL                = null_real !K, Arrhenius lower boundary
        real                             :: TAH                = null_real !K, Arrhenius upper boundary
        real                             :: F_FIX              = null_real !adim, constant food density parameter
        real                             :: PAM_FIX            = null_real !Jd-1cm-2, bivalve sur-spec assimilation rate
        real                             :: delta_M            = null_real !cm(volumetric)/cm(real), shape coefficient
        real                             :: Em                 = null_real !cm(volumetric)/cm(real), shape coefficient        
        real                             :: LifeSpan           = null_real !years, max life span
        real                             :: m_wrongSettle      = null_real !/d, fraction of the ones that settle in the wrong place
        real                             :: m_velocity         = null_real !/d, fraction of individuals that die from high velocity
        real                             :: MAX_velocity       = null_real !m/s, maximum  water velocity tolerable for this species
        real                             :: m_natural          = null_real !/d, constant natural mortality rate
        real                             :: m_spat             = null_real !/d, constant natural mortality rate
        logical                          :: DensityLimOption   = .false.   !density limitation?
        real                             :: MAX_density        = null_real !3000 #/m2, maxium density found in field observations       
        real                             :: v_cond             = null_real !cm/d, energy conductance
        real                             :: kappa              = null_real !adim, allocation fraction to growth/SomMaintenace
        real                             :: kap_R              = null_real !adim, reproduction efficiency
        real                             :: pM                 = null_real !J/(d.cm3), volume specific somatic maintenace
        real                             :: EG                 = null_real !J/cm3(volumetric), energy costs for structure
        real                             :: EHb                = null_real !J, maturity threshold for birth
        real                             :: EHp                = null_real !J, maturity threshold for puberty
        real                             :: Crm                = null_real !m3/d.cm2, maximum clearance rate
        real                             :: JX1Fm              = null_real !molC/(d.cm2),max surf area-specific filt algae
        real                             :: JX0Fm              = null_real !g/(d.cm2),max surf area-specific filt for inorgmat
        real                             :: ro_X1              = null_real !adim, binding  probability for algae
        real                             :: ro_X0              = null_real !adim, binding  probability for inorganic material        
        real                             :: JX1Im              = null_real !molC/d, max surf area-spec ing rate for algae
        real                             :: JX0Im              = null_real !molC/d, max surf area-spec ing rate for inorg mat
        real                             :: YEX                = null_real !molCE/molCX, yield coef of reser in algae struc       
        real                             :: GSR_MIN            = null_real !molCE/molCV, min GSR in the organism
        real                             :: GSR_SPAWN          = null_real !molCE/molCV, gonado-somatic ratio to spawn
        real                             :: T_SPAWN            = null_real !oC, minimum temperature for spawning
        real                             :: MIN_SPAWN_TIME     = null_real !d, minimum time interval between spawning events
        real                             :: ME_0               = null_real !molC, reserves in an embryo at optimal food
        real                             :: MEb                = null_real !molC, reserves in a new born at optimal food
        real                             :: MVb                = null_real !molC, structure in a new born at optimal food
        real                             :: MHb                = null_real !molC, maturity in a new born at optimal food
        real                             :: Lb                 = null_real !molC, length in a new born at optimal food
        real                             :: d_V                = null_real !g(dw)/cm3, density of  structure
        real                             :: mu_E               = null_real !J/molC(reser), chemical potential of reserves          
        integer                          :: SIMPLE_ASSI        = null_int  !1/0,Compute simple assimilation? 
        integer                          :: SIMPLE_TEMP        = null_int  !1/0,Compute simple temperature correction factor?
    end type T_IndividualParameters          

    type     T_AuxiliarParameters       
        real                             :: C_AtomicMass      = 12        !mgC/molC
        real                             :: H_AtomicMass      = 1         !mgH/molH
        real                             :: O_AtomicMass      = 16        !mgO/molO
        real                             :: N_AtomicMass      = 14        !mgN/molN
        real                             :: P_AtomicMass      = 31        !mgP/molP
        real                             :: TempCorrection    = null_real !adim, temperature correction factor
        real                             :: WE                = null_real !gDW/molC, AFDW to carbonV convertion
        real                             :: WV                = null_real !gDW/molC, AFDW to carbonE convertion
        real                             :: Mv                = null_real !molC(struc)/cm3, volume specific struc mass
        real                             :: MHb               = null_real !molC, maturity threshold for birth
        real                             :: MHp               = null_real !molC, maturity threshold for puberty
        real                             :: y_VE              = null_real !molCV/molCE,yield coefficient of struc on reser
        real                             :: kM                = null_real !/d, somatic maintenace rate coefficient
        real                             :: kJ                = null_real !/d, maturity maintenace rate coefficient
        real                             :: Lm                = null_real !cm, maximum length of the species
    end type T_AuxiliarParameters 

    type     T_BivalveCondition       
        real                             :: Vol               = null_real !cm3, volumetric length
        real                             :: E                 = null_real !molC(reser)/cm3, reserve density
        real                             :: ScaledE           = null_real !adim, scaled reserves density
        real                             :: GSR               = null_real !molC(reser)/molC(total), gonado-somatic ratio
        real                             :: DW                = null_real !g(dw), organism total dry weight
        real                             :: TotalmolC         = null_real !molC, organism total molC
        real                             :: TotalmolN         = null_real !molN, organism total molN
        real                             :: TotalmolP         = null_real !molP, organism total molP
    end type T_BivalveCondition 

    type     T_ByElement          
        real                             :: C    = null_real !molC/d
        real                             :: H    = null_real !molH/d
        real                             :: O    = null_real !molO/d
        real                             :: N    = null_real !molN/d
        real                             :: P    = null_real !molP/d
    end type T_ByElement
    
    type   T_Particles      
        type (T_ID           )           :: ID
        type (T_Composition  )           :: Composition
        type (T_Ratios       )           :: Ratios                             
        integer                          :: ParticleIndex          = null_int
        integer                          :: RatioVariable          = null_int   !1/0 ratios NC and PC
        integer                          :: Organic                = null_int   !1/0, is this an organic particle?
        integer                          :: Silica                 = null_int   !1/0, this particle has silica?
        real                             :: Size                   = null_real  !cm, mean length of the food
        real                             :: F_E                    = null_real  !molCReserves/molCTotalFood, fraction of res
        logical                          :: Larvae                 = .false.
        integer                          :: LarvaeSpeciesID        = null_int   !species ID
        real                             :: LarvaeDensity          = null_real  !#/m3
        real                             :: LarvaeBiomass          = null_real  !molC/#
        real                             :: LarvaeSelfPredated     = null_real  !#/m3
        real                             :: LarvaePredatedByOthers = null_real  !#/m3
        real                             :: Total_CR               = null_real  !from all cohorts that use this particle
        type (T_Particles),pointer       :: Next   
    end type T_Particles  

    type   T_Predator      
        type (T_ID           )           :: ID             
        real                             :: PredatorSize       = null_real  !cm, size of the predator
        real                             :: MinPreySize        = null_real  !cm, minimum size of the prey
        real                             :: MaxPreySize        = null_real  !cm, maximum size of the prey
        real                             :: Feeding_Rate       = null_real  !Feeding rate per individual
        integer                          :: Feeding_Units      = null_int   !1-#/d.ind; 2-AFDW/d.ind, 3-J/cm2.d
        integer                          :: Feeding_Time       = null_int   !1-Always; 2-LowTide, 3-HighTide
        logical                          :: Feeding            = OFF        !is the predator feeding?
        real                             :: Diet               = null_real  !fraction of mussels in the predator food
        real                             :: AfdwToC            = null_real  !conversion of afdw to Dw
        real                             :: DwToC              = null_real  !conversion of Dw to Carbon
        integer                          :: SIMPLE_TEMP        = null_int   !1/0,Compute simple temperature correction factor?
        integer                          :: CORRECT_TEMP       = null_int   !1/0,Compute simple temperature correction factor?
        real                             :: P_Tref             = null_real  !K, Rate Temperature reference, for predators
        real                             :: P_TA               = null_real  !K, Arrhenius temperature, for predators
        real                             :: P_TL               = null_real  !K, Lower Boundary tolerance range, for predators
        real                             :: P_TH               = null_real  !K, Upper Boundary tolerance range, for predators
        real                             :: P_TAL              = null_real  !K, Arrhenius lower boundary, for predators
        real                             :: P_TAH              = null_real  !K, Arrhenius upper boundary, for predators
        real                             :: TempCorrection     = null_real  !adim, temperature correction factor
        real                             :: TotalFeeding_Rate  = null_real  !Computed based on the predator abundance
        type (T_Predator),pointer        :: Next   
    end type T_Predator  

    type     T_InorganicFluxes       
        real                             :: CO2  = null_real  !molCO2/d, amount of CO2 produced
        real                             :: H2O  = null_real  !molH2O/d, amount of H2O produced
        real                             :: O2   = null_real  !molO2/d, amount of O2 consumed
        real                             :: NH3  = null_real  !molNH3/d, amount of ammonia produced
        real                             :: PO4  = null_real  !molPO4/d, amount of phosphate produced
    end type T_InorganicFluxes       

    type   T_Processes 
        real                             :: ClearanceRate                  = null_real  !m3/d, ClearanceRate
        real                             :: FilteredInorganic              = null_real  !g/d, filtered inorganic material 
        type(T_ByElement )               :: FilteredFood                      
        real                             :: IngestionInorganic             = null_real  !mg/d, ingested inorganic material 
        type(T_ByElement )               :: IngestionFood                                
        real                             :: PFContributionInorganic        = null_real  !mg/d, pseudofaeces from inorg material
        type(T_ByElement )               :: PFContributionFood                        
        type(T_ByElement )               :: Assimilation                                
        real                             :: FaecesContributionInorganic    = null_real  !mg/d, faeces from inorganic material
        type(T_ByElement )               :: FaecesContributionFood                       
        real                             :: SomaticMaintenance             = null_real  !molC(res)/d, JEM, somatic maintenance
        real                             :: Mobilization                   = null_real  !molC(res)/d, JEC, mobilization flux
        real                             :: ReservesDynamics               = null_real  !molC(res)/d, JE, reserves dynamics
        real                             :: ToGrowthAndSomatic             = null_real  !molC(res)/d, k fraction
        real                             :: ToGrowth                       = null_real  !molC(res)/d,  flux to growth
        real                             :: GametesLoss                    = null_real  !molC(res)/d, loss for maintenance
        real                             :: StructureLoss                  = null_real  !molC(res)/d, to shrink
        real                             :: SomaticMaintenanceNeeds        = null_real  !molC(res)/d, somatic maintenance needs
        real                             :: StructureDynamics              = null_real  !molC(str)/d, JV, flux to growth
        real                             :: ToMaturityAndReproduction      = null_real  !molC(res)/d, flux to mat and reprod     
        real                             :: MaturityMaintenance            = null_real  !molC(res)/d, JEJ, maturity maintenance
        real                             :: FluxToMatORRepr                = null_real  !molC(res)/d, JER, flux to mat or reprod     
        real                             :: MaturityLoss                   = null_real  !molC(res)/d, lose of maturity
        real                             :: FluxToGametes                  = null_real  !molC(res)/d, JER_R, flux to reproduction
        real                             :: FluxToMaturity                 = null_real  !molC(res)/d, JER_M, flux to maturity
        real                             :: MaturityDynamics               = null_real  !molC(str)/d, maturity dynamics
        real                             :: RemainMRReproduction           = null_real  !molC(str)/d, remains in reprod buffer
        real                             :: Spawning                       = null_real  !molC(res)/d, JESpaw, gamets spawned    
        real                             :: SpawningOverhead               = null_real  !molC(res)/d, JESpawOV, before spawning    
        real                             :: GametesToRelease               = null_real  !#/d, JESpaw, number of gametes released    
        real                             :: NewbornsThisCohort             = null_real  !#/d, number of new born individuals  
        real                             :: NONewbornsThisCohort           = null_real  !#/d, number of new born individuals  
        real                             :: ReproductionDynamics           = null_real  !molC(res)/d, JR, flux to reprod
        type(T_InorganicFluxes)          :: InorganicFluxes             
        real                             :: DeathByAge                     = 0.0  !#/d, dead by age             
        real                             :: DeathByOxygen                  = 0.0  !#/d, dead by lack of oxygen             
        real                             :: DeathByStarvation              = 0.0  !#/d, dead by starvation            
        real                             :: DeathByExtraStarvation         = 0.0  !#/d, extra mortality by starvation            
        real                             :: DeathByNatural                 = 0.0  !#/d, dead by natural            
        real                             :: DeathByVelocity                = 0.0  !#/d, dead by velocity            
        real                             :: DeathByWrongSettlement         = 0.0  !#/d, dead by settle in the wrong place  
        real                             :: DeathByDensityLimit            = 0.0  !#/d, dead by settle in the wrong place  
        real                             :: PredationByShrimps             = 0.0  !#/d, dead by shrimp             
        real                             :: PredationByCrabs               = 0.0  !#/d, dead by crab             
        real                             :: PredationByOysterCatchers      = 0.0  !#/d, dead by oystercatcher             
        real                             :: PredationByEiderDucks          = 0.0  !#/d, dead by eider ducks             
        real                             :: PredationByHerringGull         = 0.0  !#/d, dead by HerringGulls             
        real                             :: DeathByLowNumbers              = 0.0  !#/d, dead by low numbers in the cohort             
        real                             :: DeathBySelfPredation           = 0.0  !#/d, dead by self predation           
        real                             :: DeathByLarvaePredationByOthers = 0.0  !#/d, dead by self predation 
    end type T_Processes        


    type   T_PopulationProcesses 
        real                             :: TNStartTimeStep              = 0.0  
        real                             :: TNNonLarvaeStartTimeStep     = 0.0  
        real                             :: TN                           = 0.0  
        real                             :: NCoh                         = 0.0  
        real                             :: TBio                         = 0.0                 
        real                             :: Cr                           = 0.0  
        real                             :: Fil                          = 0.0                 
        real                             :: Ing                          = 0.0  
        real                             :: Ass                          = 0.0                 
        real                             :: CO2                          = 0.0                 
        real                             :: H2O                          = 0.0  
        real                             :: O                            = 0.0                 
        real                             :: NH3                          = 0.0  
        real                             :: PO4                          = 0.0  
        real                             :: m_A                          = 0.0  
        real                             :: m_O                          = 0.0  
        real                             :: m_F                          = 0.0  
        real                             :: m_nat                        = 0.0  
        real                             :: m_shr                        = 0.0  
        real                             :: m_cra                        = 0.0  
        real                             :: m_oys                        = 0.0  
        real                             :: m_duck                       = 0.0  
        real                             :: m_gull                       = 0.0  
        real                             :: m_low                        = 0.0
        real                             :: m_self                       = 0.0  
        real                             :: m_others                     = 0.0  
        real                             :: m_vel                        = 0.0  
        real                             :: m_settle                     = 0.0  
        real                             :: Massm_A                      = 0.0  
        real                             :: Massm_O                      = 0.0  
        real                             :: Massm_F                      = 0.0  
        real                             :: Massm_nat                    = 0.0  
        real                             :: Massm_shr                    = 0.0  
        real                             :: Massm_cra                    = 0.0  
        real                             :: Massm_oys                    = 0.0  
        real                             :: Massm_duck                   = 0.0  
        real                             :: Massm_gull                   = 0.0  
        real                             :: Massm_low                    = 0.0  
        real                             :: Massm_self                   = 0.0  
        real                             :: Massm_others                 = 0.0  
        real                             :: Massm_vel                    = 0.0  
        real                             :: Massm_settle                 = 0.0  
        real                             :: TNField                      = 0.0             
        real                             :: MaxLength                    = 0.0 !#cm,             
        real                             :: LastLength                   = null_real !length of the last cohort before death
        real, pointer, dimension(:)      :: SumLogAllMortalityInNumbers  => null() !Product of mortalities in numbers
        real, pointer, dimension(:)      :: SumAllMortalityInMass        => null() !sum of all death rates of all instants
        real, pointer, dimension(:)      :: AverageMortalityInNumbers    => null() !geometricAverage of all death rates 
        integer                          :: nInstantsForAverage          = null_int
        integer                          :: nSpawning                    = null_int !number of spawning events
        real                             :: nNewborns                    = 0
    end type T_PopulationProcesses        
     
    type     T_Output
        integer, dimension(:), pointer   :: Unit                         => null()
        integer                          :: nParticles
        character(len=StringLength)      :: FileName    = '   '    
        real, dimension(:), pointer      :: Aux
    end type T_Output

    type     T_Cohort
        type(T_ID                    )   :: ID
        type(T_StateIndex            )   :: StateIndex
        type(T_BivalveCondition      )   :: BivalveCondition
        type(T_Processes             )   :: Processes
        real,  pointer, dimension(:,:)   :: FeedingOn   => null()  !to store, Col = Filtered ingested assimilated (molC/g.d.ind)
        type(T_Output                )   :: CohortOutput
        logical                          :: Larvae      = OFF
        integer                          :: Dead        = 0
        integer                          :: GlobalDeath = 1
        !integer,  pointer, dimension(:)  :: LarvaeState => null()  
        logical                          :: AtLeastOneLarvae         = .false.
        logical                          :: PreviousAtLeastOneLarvae = .false.
        type(T_Cohort ), pointer         :: Next
    end type T_Cohort

    type     T_Species
        type(T_ID                   )    :: ID
        type(T_Cohort)    , pointer      :: FirstCohort
        type(T_Particles) , pointer      :: FirstParticles
        type(T_Predator)  , pointer      :: FirstPredator
        type(T_SpeciesComposition   )    :: SpeciesComposition
        type(T_IndividualParameters )    :: IndividualParameters
        type(T_AuxiliarParameters   )    :: AuxiliarParameters
        type(T_Output               )    :: PopulationOutput
        type(T_Output               )    :: SizeDistributionOutput
        type(T_Output               )    :: TestingParametersOutput
        type(T_PopulationProcesses  )    :: PopulationProcesses
        logical                          :: CohortOutput          = OFF 
        logical                          :: Population            = OFF
        logical                          :: FeedOnLarvae          = OFF
        logical                          :: ExtraStarvation       = OFF
        real                             :: LarvaeMaxSize         = null_real
        logical                          :: BySizeOutput          = OFF 
        character(len=StringLength)      :: Testing_File          = null_str
        real                             :: SizeStep              = null_real
        real                             :: Max_SizeClass         = null_real
        real, pointer, dimension(:)      :: SizeClasses           => null()
        real, pointer, dimension(:)      :: SizeFrequency         => null()
        integer                          :: nSizeClasses          = null_int
        character(len=4)                 :: SizeClassesCharFormat = '   '
        integer                          :: Initial_nCohorts      = 0
        real                             :: MinObsLength          = 0
        integer                          :: LastCohortID          = 0
        integer                          :: nCohorts              = 0
        logical                          :: NewbornCohort         = .false.
        integer                          :: nParticles            = 0
        integer                          :: nPredator             = 0
        real                             :: Total_CR              = null_real  !total from all cohorts that filter the particles
        real                             :: Total_CR_Larvae       = null_real  !total including the larvae
        real, pointer, dimension(:)      :: SettlementProbability => null()
        type(T_Species),    pointer      :: Next
    end type T_Species
    
    type     T_RestartSpecies
        character(len = StringLength)    :: Name = " "
        integer                          :: nCohorts
        integer, dimension(:), pointer   :: CohortIDs     => null()
    end type T_RestartSpecies

    type     T_Bivalve
        integer                              :: InstanceID
        integer                              :: ObjTime                  = 0
        integer                              :: ObjEnterData             = 0
        integer                              :: DensityUnits             = 0 ! 0: m2, 1:m3
        type(T_Time)                         :: InitialDate, FinalDate
        logical                              :: Old                      = .false.
        real                                 :: DT                       = null_real
        real                                 :: DTDay                    = null_real
        type (T_Size1D)                      :: Array                    
        type (T_Size1D)                      :: Prop                   
        integer                              :: nSpecies                 = 0
        real                                 :: LackOfFood               = 0
        integer, dimension(:), pointer       :: PropertyList             => null()
        integer                              :: nPropertiesFromBivalve   = 0
        integer                              :: nCohortProperties        = 7          !Each cohort has 7 associated properties
        real                                 :: MinNumber                = null_real
        integer, dimension(:), pointer       :: ListDeadIDs              => null()
        integer                              :: nLastDeadID              = 0
        integer, dimension(:), pointer       :: ListNewbornsIDs          => null()    !List of SpeciesID with newborns  
        integer                              :: nLastNewbornsID          = 0
        real, dimension(:,:)   , pointer     :: MatrixNewborns           => null()    !col = SpeciesID | Index +1  (nNewborns)
        real                                 :: DT_OutputTime            = null_real
        logical                              :: Testing_Parameters       = OFF
        logical                              :: OutputON                 = OFF
        integer                              :: TotalOutputs, NextOutPut = null_int
        type (T_Time), pointer, dimension(:) :: BivalveOutputTimes
        integer, dimension(1:30)             :: IndexOutputs             = null_int 
        integer                              :: nIndexOutputs            = null_int              
        real                                 :: ConvertionFactor         = 1.0 !convertion from m2 to m3
        type (T_Output        )              :: MassOutput
        type (T_PropIndex     )              :: PropIndex
        type (T_ComputeOptions)              :: ComputeOptions
        type (T_Species       ), pointer     :: FirstSpecies
        type (T_ExternalVar   )              :: ExternalVar
        type (T_Bivalve       ), pointer     :: Next
        real                                 :: MassLoss                 = 0.0
        real                                 :: MaxTNField               = 0.0  ! #/m2
        !character(len = PathLength)         :: PathFileName = '/home/saraiva/00_Projects/Parametric/Running/' !biocluster
        character(len = PathLength)          :: PathFileName = ''  
        character(len = PathLength)          :: InitialFileName = ''  
        character(len = PathLength)          :: FinalFileName = '' 
        
        type (T_RestartSpecies), dimension(:), pointer  :: RestartSpecies
        integer                              :: nRestartSpecies
        logical                              :: ComputeThisIndex        = .true.
        logical                              :: OutputThisIndex         = .false.
    end type T_Bivalve


    !Global Module Variables
    type (T_Bivalve       ), pointer         :: FirstObjBivalve
    type (T_Bivalve       ), pointer         :: Me

    !integer            :: mBivalve_ 

    !-------------------------------------------------------------------------------

    contains

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRU

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructBivalve(ObjBivalveID, FileName, BeginTime, EndTime, ArraySize, STAT)

        !Arguments------------------------------------------------------------------
        integer                           :: ObjBivalveID
        character(len=*)                  :: FileName
        type(T_Time)                      :: BeginTime, EndTime
        type(T_Size1D)                    :: ArraySize
        integer, optional, intent(OUT)    :: STAT     

        !External-------------------------------------------------------------------
        integer                           :: ready_, STAT_CALL         

        !Local----------------------------------------------------------------------
        integer                           :: STAT_

        !---------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mBivalve_)) then
            nullify (FirstObjBivalve)
            call RegisterModule (mBivalve_) 
        endif

        call Ready(ObjBivalveID, ready_)    

cd1 :   if (ready_ .EQ. OFF_ERR_) then


            call AllocateInstance

            call ConstructEnterData(Me%ObjEnterData, trim(FileName), STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_)  &         
                stop 'Subroutine ConstructBivalve - ModuleBivalve - ERR01'

            Me%InitialDate  = BeginTime
            Me%FinalDate    = EndTime
            
            Me%Array%ILB = ArraySize%ILB
            Me%Array%IUB = ArraySize%IUB

            call ReadDataBivalve
            
            call PropertyIndexNumber

            call ConstructPropertyList

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 

            if (STAT_CALL .NE. SUCCESS_)  &         
                stop 'Subroutine ConstructBivalve - ModuleBivalve - ERR02'

            !Returns ID
            ObjBivalveID          = Me%InstanceID

            STAT_ = SUCCESS_

        else 

            stop 'Subroutine ConstructBivalve - ModuleBivalve - ERR03' 

        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine ConstructBivalve

    !-------------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments------------------------------------------------------------------

        !Local----------------------------------------------------------------------
        type (T_Bivalve), pointer           :: NewObjBivalve
        type (T_Bivalve), pointer           :: PreviousObjBivalve


        !Allocates new instance
        allocate (NewObjBivalve     )
        nullify  (NewObjBivalve%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjBivalve)) then
            FirstObjBivalve          => NewObjBivalve
            Me         => NewObjBivalve
        else
            PreviousObjBivalve       => FirstObjBivalve
            Me         => FirstObjBivalve%Next
            
            do while (associated(Me))
                PreviousObjBivalve   => Me
                Me     => Me%Next
            enddo
            
            Me         => NewObjBivalve
            PreviousObjBivalve%Next  => NewObjBivalve
        endif

        Me%InstanceID = RegisterNewInstance (mBivalve_)

    end subroutine AllocateInstance

    !-------------------------------------------------------------------------------

    subroutine ReadDataBivalve
    
        !Local----------------------------------------------------------------------
        integer                         :: i
        type(T_Species), pointer        :: Species

        !Begin----------------------------------------------------------------------

        call ConstructGlobalVariables
        
        call ConstructSpecies

        if (Me%OutputON) then
        
            Species => Me%FirstSpecies
            do while(associated(Species))
            
                if(Species%Population)then
                    allocate(Species%PopulationOutput%Unit(1:Me%nIndexOutputs))
                endif
                
                if (Species%BySizeOutput) then
                    allocate(Species%SizeDistributionOutput%Unit(1:Me%nIndexOutputs))
                endif
                
                Species => Species%Next
            enddo
            
            if(Me%ComputeOptions%MassBalance) then
                allocate(Me%MassOutput%Unit(1:Me%nIndexOutputs))
            endif
        
            do i=1, Me%nIndexOutputs
                call ConstructOutputs(i)
            enddo
        end if

    end subroutine ReadDataBivalve

    !-------------------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !External-------------------------------------------------------------------
        integer             :: flag, STAT_CALL

        !Begin----------------------------------------------------------------------
        
        call ReadFileName("ROOT_SRT", Me%PathFileName, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR00'
            
        call GetOutPutTime(Me%ObjEnterData,                                         &
                           CurrentTime      = Me%InitialDate ,                      &
                           EndTime          = Me%FinalDate,                         &
                           keyword          = 'BIVALVE_OUTPUT_TIME',                &
                           SearchType       = FromFile,                             &
                           OutPutsTime      = Me%BivalveOutputTimes,                &
                           OutPutsOn        = Me%OutputON,                          &
                           OutPutsNumber    = Me%TotalOutputs,                      &
                           STAT             = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                          &
            stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR01'

        if (Me%OutputON) then

            Me%NextOutPut = 1

        endif 

        call GetData(Me%DT                                                  , &
                    Me%ObjEnterData, flag                                   , &
                    SearchType   = FromFile                                 , &
                    keyword      ='DT'                                      , &       
                    default      = 3600.0                                   , & 
                    ClientModule = 'ModuleBivalve'                          , &
                    STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                          &
            stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR02'

cd1:    if (flag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Keyword DT not found in ModuleBivalve data file.'
            write(*,*) 'Subroutine ConstructGlobalVariables-ModuleBivalve-WRN01'
            write(*,*) 'Assumed ', Me%DT , &
                       'seconds (',  Me%DT / 3600.0, 'hour).'
            write(*,*) 
        end if cd1

        !Convert DT [seconds] in DT [day]
        Me%DTDay = Me%DT / 24.0 / 60.0 / 60.0

        call GetData(Me%ComputeOptions%PelagicModel                         , &
                    Me%ObjEnterData,  flag                                  , &
                    SearchType   = FromFile                                 , &
                    keyword      = 'PELAGIC_MODEL'                          , &
                    ClientModule = 'ModuleBivalve'                          , &
                    STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                          &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR10'

cd2:    if(flag==0)then
            write(*,*)'Please define the pelagic model to couple with ModuleBivalve'
            stop 'ConstructGlobalVariables - ModuleBivalve - ERR03'
        end if cd2

cd3:    if((Me%ComputeOptions%PelagicModel .ne. WaterQualityModel .and. Me%ComputeOptions%PelagicModel .ne. LifeModel))then
            write(*,*)'Pelagic model to couple with ModuleBivalve must be one of the following:'
            write(*,*)trim(WaterQualityModel)
            write(*,*)trim(LifeModel)
            stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR20'
        end if cd3

        call GetData(Me%ComputeOptions%Nitrogen                             , &
                    Me%ObjEnterData, flag                                   , &
                    SearchType = FromFile                                   , &
                    keyword='NITROGEN'                                      , &
                    ClientModule = 'ModuleBivalve'                          , &
                    STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                          &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR30'


        call GetData(Me%ComputeOptions%Phosphorus                           , &
                    Me%ObjEnterData, flag                                   , &
                    SearchType   = FromFile                                 , &
                    keyword      = 'PHOSPHOR'                               , &
                    Default      = .false.                                  , &
                    ClientModule = 'ModuleBivalve'                          , &
                    STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                           &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR40'


        call GetData(Me%ComputeOptions%SimpleFiltration                    , &
                    Me%ObjEnterData, flag                                  , &
                    SearchType   = FromFile                                , &
                    keyword      = 'SIMPLE_FILTRATION'                     , &
                    Default      = .false.                                 , &
                    ClientModule = 'ModuleBivalve'                         , &
                    STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                          &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR50'

        call GetData(Me%ComputeOptions%CorrectFiltration                   , &
                    Me%ObjEnterData, flag                                  , &
                    SearchType   = FromFile                                , &
                    keyword      = 'CORRECT_FILTRATION'                    , &
                    Default      = .false.                                 , &
                    ClientModule = 'ModuleBivalve'                         , &
                    STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                          &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR60'

        call GetData(Me%ComputeOptions%MassBalance                         , &
                    Me%ObjEnterData, flag                                  , &
                    SearchType   = FromFile                                , &
                    keyword      = 'MASS_BALANCE'                          , &
                    Default      = .false.                                 , &
                    ClientModule = 'ModuleBivalve'                         , &
                    STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                          &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR70'

        call GetData(Me%MinNumber                                         , &
                    Me%ObjEnterData, flag                                 , &
                    SearchType   = FromFile                               , &
                    keyword      = 'MIN_NUMBER'                           , &
                    Default      = 1.0e-1                                 , &
                    ClientModule = 'ModuleBivalve'                        , &
                    STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                         &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR80'

        call GetData(Me%DT_OutputTime                                    , &
                    Me%ObjEnterData, flag                                , &
                    SearchType   = FromFile                              , &
                    keyword      = 'DT_OUTPUT_TIME'                      , &
                    default      = Me%DT                                 , &
                    ClientModule = 'ModuleBivalve'                       , &
                    STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                       &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR90'
        
        call GetData(Me%Testing_Parameters                               , &
                    Me%ObjEnterData, flag                                , &
                    SearchType   = FromFile                              , &
                    keyword      = 'TESTING_PARAMETERS'                  , &
                    default      = .false.                               , &
                    ClientModule = 'ModuleBivalve'                       , &
                    STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                       &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR100'
        
        call GetData(Me%DensityUnits                                     , &
                    Me%ObjEnterData, flag                                , &
                    SearchType   = FromFile                              , &
                    keyword      = 'DENSITY_UNITS'                       , &
                    default      = 0                                     , &
                    ClientModule = 'ModuleBivalve'                       , &
                    STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                       &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR110'
        
        call GetData(Me%IndexOutputs                                     , &
                    Me%ObjEnterData, flag                                , &
                    SearchType   = FromFile                              , &
                    keyword      = 'INDEX_OUTPUTS'                       , &
                    ClientModule = 'ModuleBivalve'                       , &
                    STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)then
            !By default 30 values are read so this error always exist.
            if (STAT_CALL /= SIZE_ERR_) then 
                stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR120'
            else
                Me%nIndexOutputs = flag
            endif
        else
            Me%nIndexOutputs = 0
        endif
        
        call GetData(Me%Old                                              , &
                     Me%ObjEnterData, flag                               , &
                     SearchType   = FromFile                             , &
                     keyword      = 'OLD'                                , &       
                     default      = .false.                              , & 
                     ClientModule = 'ModuleBivalve'                      , &
                     STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                       &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR130'
        
        if(Me%Old)then
            
            call ReadFileName("BIV_INI", Me%InitialFileName, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                       &
            stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR140'
            
            call ReadInitialBivalveFile

        endif

        call ReadFileName("BIV_FIN", Me%FinalFileName, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                       &
        stop 'Subroutine ConstructGlobalVariables - ModuleBivalve - ERR150'
        
        allocate(Me%ExternalVar%InitialPhyto (Me%Array%ILB:Me%Array%IUB))
        allocate(Me%ExternalVar%InitialShrimp(Me%Array%ILB:Me%Array%IUB))

        
    end subroutine ConstructGlobalVariables

    !-------------------------------------------------------------------------------
    
    subroutine ReadInitialBivalveFile
    
        !Local----------------------------------------------------------------------
        integer                                     :: ClientNumber, STAT_CALL, flag
        logical                                     :: BlockFound, BlockCohortsFound
        integer                                     :: iCohort, iSpecies
        integer                                     :: iLine, FirstLine, LastLine, ID
        integer                                     :: ObjEnterData = 0
        
        !Begin----------------------------------------------------------------------

        call ConstructEnterData(ObjEnterData, trim(Me%InitialFileName), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadInitialBivalveFile - ModuleBivalve - ERR01'        
                
        call GetData(Me%nRestartSpecies                         , &
                     ObjEnterData, flag                         , &
                     SearchType   = FromFile                    , &
                     keyword      = 'NUMBER_OF_SPECIES'         , &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadInitialBivalveFile - ModuleBivalve - ERR10'
        
        allocate(Me%RestartSpecies(1:Me%nRestartSpecies))
        
        iSpecies = 0   

do1 :   do
            call ExtractBlockFromBuffer(ObjEnterData                       , &
                                        ClientNumber    = ClientNumber     , &
                                        block_begin     = '<begin_species>', &
                                        block_end       = '<end_species>'  , &
                                        BlockFound      = BlockFound       , &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then

                    iSpecies = iSpecies + 1

                    call GetData(Me%RestartSpecies(iSpecies)%Name           , &
                                 ObjEnterData, flag                         , &
                                 SearchType   = FromBlock                   , &
                                 keyword      = 'NAME'                      , &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInitialBivalveFile - ModuleBivalve - ERR20'

                    call GetData(Me%RestartSpecies(iSpecies)%nCohorts       , &
                                 ObjEnterData, flag                         , &
                                 SearchType   = FromBlock                   , &
                                 keyword      = 'NUMBER_OF_COHORTS'         , &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInitialBivalveFile - ModuleBivalve - ERR30'
                    
                    allocate(Me%RestartSpecies(iSpecies)%CohortIDs(1:Me%RestartSpecies(iSpecies)%nCohorts))
                    
                    do iCohort = 1, Me%RestartSpecies(iSpecies)%nCohorts
                    
                        BlockCohortsFound = .false.
                    
                        call ExtractBlockFromBlock(ObjEnterData, ClientNumber,      &
                                                   '<begin_cohort>', '<end_cohort>',&
                                                   BlockCohortsFound,               &
                                                   FirstLine = FirstLine,           &
                                                   LastLine  = LastLine,            &
                                                   STAT      = STAT_CALL)
                        if (STAT_CALL .EQ. SUCCESS_) then    
                            if (BlockCohortsFound) then                 
                            
                                do iLine = FirstLine+1, LastLine-1
                                
                                    call GetData(ID, ObjEnterData, flag, Buffer_Line  = iLine, STAT = STAT_CALL)
                                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadInitialBivalveFile - ModuleBivalve - ERR40'  
                                    
                                    Me%RestartSpecies(iSpecies)%CohortIDs(iCohort) = ID
                                            
                                enddo

                            else
                                stop 'Subroutine ReadInitialBivalveFile - ModuleBivalve - ERR60'
                            endif      
                        else
                            stop 'Subroutine ReadInitialBivalveFile - ModuleBivalve - ERR70'
                        endif
                    enddo

                else cd2
                
                    call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadInitialBivalveFile - ModuleBivalve - ERR80'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Subroutine ReadInitialBivalveFile - ModuleBivalve - ERR90'
            else cd1
                stop 'Subroutine ReadInitialBivalveFile - ModuleBivalve - ERR100'
            end if cd1
        end do do1


        call KillEnterData(ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadInitialBivalveFile - ModuleBivalve - ERR110'        
    
    
    
    end subroutine ReadInitialBivalveFile

    !-------------------------------------------------------------------------------

    subroutine ConstructSpecies

        !Arguments------------------------------------------------------------------

        !Local----------------------------------------------------------------------
        type (T_Species)           , pointer        :: NewSpecies
        integer                                     :: ClientNumber, STAT_CALL
        logical                                     :: BlockFound
        integer                                     :: iCohort, iSpecies, RestartCohortID

        !Begin----------------------------------------------------------------------

        iSpecies = 0
        
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData                    , &
                                        ClientNumber    = ClientNumber     , &
                                        block_begin     = '<begin_species>', &
                                        block_end       = '<end_species>'  , &
                                        BlockFound      = BlockFound       , &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then        

                    call AddSpecies (NewSpecies)

                    call ConstructSpeciesParameters (NewSpecies, ClientNumber)
                    
                    if(.not. Me%Old)then
                    
                        do iCohort = 1, NewSpecies%Initial_nCohorts

                            call ConstructCohort (NewSpecies)

                        end do
                    
                    else
                        
                        !continuacao de calculo
                        iSpecies = iSpecies + 1
                        
                        if(trim(Me%RestartSpecies(iSpecies)%Name) .NE. trim(NewSpecies%ID%Name))then
                            write(*,*)"Inconsistent species names between Bivalve restart file and input file"
                            stop 'Subroutine ConstructSpecies - ModuleBivalve - ERRO1'
                        end if
                        
                        do iCohort = 1, Me%RestartSpecies(iSpecies)%nCohorts
                        
                            RestartCohortID = Me%RestartSpecies(iSpecies)%CohortIDs(iCohort)
                        
                            call ConstructCohort (NewSpecies, RestartCohortID)
                            
                        enddo
                        
                    endif

                        

                    nullify(NewSpecies)

                else cd2
                
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)        &
                        stop 'Subroutine ConstructSpecies - ModuleBivalve - ERRO1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Subroutine ConstructSpecies - ModuleBivalve - ERR20'
            else cd1
                stop 'Subroutine ConstructSpecies - ModuleBivalve - ERR30'
            end if cd1
        end do do1

    end subroutine ConstructSpecies

    !-------------------------------------------------------------------------------

    subroutine AddSpecies (ObjSpecies)

        !Arguments------------------------------------------------------------------
        type (T_Species),      pointer           :: ObjSpecies
        !Local----------------------------------------------------------------------
        type (T_Species),      pointer           :: PreviousSpecies
        type (T_Species),      pointer           :: NewSpecies
        integer, save                            :: NextSpeciesID = 1

        !Allocates new Species
        allocate (NewSpecies     )
        nullify  (NewSpecies%Next)

        !Insert new Species into list 
cd1:    if (.not. associated(Me%FirstSpecies)) then
            Me%FirstSpecies             => NewSpecies
            ObjSpecies    => NewSpecies
        else
            PreviousSpecies             => Me%FirstSpecies
            ObjSpecies    => Me%FirstSpecies%Next

do1:        do while (associated(ObjSpecies))
                PreviousSpecies         => ObjSpecies
                ObjSpecies=> ObjSpecies%Next
            enddo do1
            ObjSpecies    => NewSpecies
            PreviousSpecies%Next        => NewSpecies
        endif cd1

        !Attributes ID
        ObjSpecies%ID%ID  = NextSpeciesID
        NextSpeciesID     = NextSpeciesID + 1

        Me%nSpecies       = Me%nSpecies + 1

    end subroutine AddSpecies

    !-------------------------------------------------------------------------------

    subroutine ConstructCohort(NewSpecies, CohortID)

        !Arguments------------------------------------------------------------------
        type (T_Species)        , pointer             :: NewSpecies
        integer, optional                             :: CohortID

        !Local----------------------------------------------------------------------
        type(T_Cohort)          , pointer             :: NewCohort
        character(len=5)                              :: CohortIDStr
        integer                                       :: i !, Index

        !Begin----------------------------------------------------------------------

        allocate(NewCohort)
        
        if(present(CohortID))then
            call AddCohort(NewSpecies, NewCohort, CohortID)
        else
            call AddCohort(NewSpecies, NewCohort)
        endif


        write(CohortIDStr, ('(i5)'))NewCohort%ID%ID

        NewCohort%ID%Name = trim(adjustl(NewSpecies%ID%Name))//" cohort "//trim(adjustl(CohortIDStr))

        if (NewSpecies%CohortOutput) then
        
            allocate(NewCohort%CohortOutput%Unit(1:Me%nIndexOutputs))
        
            do i=1, Me%nIndexOutputs
                call ConstructCohortOutput (NewCohort, i)
            enddo
        end if
        
!        allocate(NewCohort%LarvaeState(Me%Array%ILB:Me%Array%IUB))
!        
!        do Index = Me%Array%ILB, Me%Array%IUB 
!            NewCohort%LarvaeState(Index) = -99
!        enddo
        
        nullify(NewCohort)

    end subroutine ConstructCohort

    !-------------------------------------------------------------------------------

    subroutine AddCohort (Species, NewCohort, NewCohortID)

        !Arguments-------------------------------------------------------------
        type (T_Species),            pointer            :: Species
        type (T_Cohort),            pointer             :: NewCohort
        integer, optional                               :: NewCohortID

        !Local-----------------------------------------------------------------
        type (T_Cohort),            pointer             :: ObjCohort

        nullify  (NewCohort%Next)

        !Insert new cohort into list
cd1:    if (.not. associated(Species%FirstCohort)) then
            Species%FirstCohort   => NewCohort
        else

            ObjCohort => Species%FirstCohort

do1:        do while (associated(ObjCohort%Next))

                ObjCohort => ObjCohort%Next
            enddo do1

            ObjCohort%Next => NewCohort
        endif cd1

        !Attributes ID
        if(present(NewCohortID))then
            NewCohort%ID%ID        = NewCohortID
        else
            NewCohort%ID%ID        = Species%LastCohortID + 1
        endif

        Species%LastCohortID       = NewCohort%ID%ID

        Species%nCohorts           = Species%nCohorts + 1

    end subroutine AddCohort

    !-------------------------------------------------------------------------------
    
    subroutine ComputeCohortLarvaeState(Cohort, Index, LarvaeMaxSize)
    
        !Arguments-------------------------------------------------------------
        type (T_Cohort),            pointer             :: Cohort
        integer                                         :: Index
        real                                            :: LarvaeMaxSize
        
        !Local-----------------------------------------------------------------

        if ((Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) .gt. 0.0) .and. &
            (Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      .gt. 0.0) .and. &
            (Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      .le. LarvaeMaxSize)) then
                        
            Cohort%Larvae = .true.
                    
        else
                    
            Cohort%Larvae = .false.
                        
        end if
    
    end subroutine ComputeCohortLarvaeState
    
    !-------------------------------------------------------------------------------

    subroutine ConstructOutputs(iIndexOutput)
        
        !Arguments-------------------------------------------------------------
        integer                                         :: iIndexOutput

        !Local----------------------------------------------------------------------
            integer                                     :: STAT_CALL,i
            type (T_Species)        , pointer           :: Species
            type (T_Predator)       , pointer           :: Predator
            character(len=900)                          :: SizeDistributionHeader
            character(len=900)                          :: OuputHeader
            character(len=500)                          :: OuputFileName
            character(len=16)                           :: IndexOutputStr
            character(len=16)                           :: SizeClassNameStr
            character(len=16)                           :: ParameterValueStr
            !character(len=16)                           :: ArgumentInComand
        !Begin----------------------------------------------------------------------
        
        !call getarg(1,ArgumentInComand) 
        write(IndexOutputStr, ('(I5)')) Me%IndexOutputs(iIndexOutput)

        Species => Me%FirstSpecies
        do while(associated(Species))
        
            if (Species%Population) then

                call UnitsManager(Species%PopulationOutput%Unit(iIndexOutput), OPEN_FILE, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ConstructOutputs - ModuleBivalve - ERR20'

                !OuputFileName = "Output/"//trim(ArgumentInComand)  !biocluster

                OuputFileName = trim(IndexOutputStr)//"_"//trim(Species%ID%Name)

                !Species population output
                if (Me%Testing_Parameters) then
                
                    Predator => Species%FirstPredator
                    do while(associated(Predator))
                    
                        write(ParameterValueStr, ('(E12.5)')) Predator%Diet
                        
                        OuputFileName = trim(OuputFileName)//'_'//trim(ParameterValueStr)
                    
                        Predator => Predator%Next
                    end do
                
                    write(ParameterValueStr, ('(E12.5)')) Me%DT

                    OuputFileName = trim(OuputFileName)//'_'//trim(ParameterValueStr)

                    write(ParameterValueStr, ('(E12.5)')) Species%IndividualParameters%m_spat

                    OuputFileName = trim(OuputFileName)//'_'//trim(ParameterValueStr)

                    write(ParameterValueStr, ('(E12.5)')) Species%IndividualParameters%m_natural

                    OuputFileName = trim(OuputFileName)//'_'//trim(ParameterValueStr)

                end if

                Species%PopulationOutput%FileName = trim(OuputFileName)
                
                OuputFileName = trim(Me%PathFileName)// &
                                trim(Species%PopulationOutput%FileName)//'_'//'population.srw'

                open(Unit = Species%PopulationOutput%Unit(iIndexOutput), File = trim(OuputFileName), Status = 'REPLACE')


                !time serie format
                call WriteDataLine(Species%PopulationOutput%Unit(iIndexOutput), "Time Serie Results File")

                call WriteDataLine(Species%PopulationOutput%Unit(iIndexOutput), "NAME", 'Population File')

                call WriteDataLine(Species%PopulationOutput%Unit(iIndexOutput), 'SERIE_INITIAL_DATA', Me%InitialDate)

                call WriteDataLine(Species%PopulationOutput%Unit(iIndexOutput), 'TIME_UNITS', 'SECONDS')


                102 format(A900)

                OuputHeader =   "!Seconds_1 YY_2 MM_3 DD_4 hh_5 mm_6 ss_7 "                                     // &
                                "#/m2_8 #_9 molC/m2_10 m3/d.m2_11 m3/d.m2_12 m3/d.m2_13 m3/d.m2_14 molC/m3_15 " // &
                                "m3/d.m2_16 m3/d.m2_17 m3/d.m2_18 m3/d.m2_19 -_20 "                             // &
                                "#/d.m2_21 #/d.m2_22 #/d.m2_23 #/d.m2_24 #/d.m2_25 "                            // &
                                "#/d.m2_26 #/d.m2_27 #/d.m2_28 #/d.m2_29 #/d.m2_30 "                            // &
                                "#/d.m2_31 #/d.m2_32 #/d.m2_33 #/d.m2_34 mol_35 "                               // &
                                "mol_36 mol_37 mol_38 mol_39 mol_40 "                                           // &
                                "mol_41 mol_42 mol_43 mol_44 mol_45 "                                           // &
                                "mol_46 mol_47 mol_48 /d_49 /d_50 "                                             // &
                                "/d_51 /d_52 /d_53 /d_54 /d_55 "                                                // &
                                "/d_56 /d_57 /d_58 /d_59 /d_60 "                                                // &
                                "/d_61 /d_62 #/m2_63 mg/l_64 #/m2_65 "                                          // &
                                "cm_66 #_67 #_68"

                write(Species%PopulationOutput%Unit(iIndexOutput), 102) OuputHeader

                OuputHeader = "Seconds_1 YY_2 MM_3 DD_4 hh_5 mm_6 ss_7 "                                        // &
                              "TN_8 NCoh_9 TBio_10 Cr_11 Fil_12 Ing_13 Ass_14 CO_15 "                           // &
                              "H2O_16 O_17 NH3_18 PO4_19 LackOfFood_20 "                                        // &
                              "m_A_21 m_O_22 m_F_23 m_nat_24 m_shr_25 "                                         // &
                              "m_cra_26 m_oys_27 m_duck_28 m_gull_29 m_low_30 "                                 // &
                              "m_self_31 m_others_32 m_vel_33 m_settle_34 TMASSm_A_35 "                         // &
                              "TMASSm_O_36 TMASSm_F_37 TMASSm_nat_38 TMASSm_shr_39 TMASSm_cra_40 "              // &
                              "TMASSm_oys_41 TMASSm_duck_42 TMASSm_gull_43 TMASSm_low_44 TMASSm_self_45 "       // &
                              "TMASSm_others_46 TMASSm_vel_47 TMASSm_settle_48 GEOm_A_49 GEOm_O_50 "            // &
                              "GEOm_F_51 GEOm_nat_52 GEOm_shr_53 GEOm_cra_54 GEOm_oys_55 "                      // &
                              "GEOm_duck_56 GEOm_gull_57 GEOm_low_58 GEOm_self_59 GEOm_others_60 "              // &
                              "GEOm_low_61 GEOm_settle_62 TNField_63 InitialPhyto_64 InitialShrimp_65 "         // &
                              "MaxLength_66 MaxTNField_67 SpawningEvents_68"
   
                write(Species%PopulationOutput%Unit(iIndexOutput), 102) OuputHeader
                
                call WriteDataLine(Species%PopulationOutput%Unit(iIndexOutput), '<BeginTimeSerie>')
                
                if (Species%BySizeOutput) then
                
                    call UnitsManager(Species%SizeDistributionOutput%Unit(iIndexOutput), OPEN_FILE, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ConstructOutputs - ModuleBivalve - ERR21'

                    Species%SizeDistributionOutput%FileName = Species%PopulationOutput%FileName
                
                    OuputFileName = trim(Me%PathFileName)//              &
                                    trim(Species%SizeDistributionOutput%FileName)//  &
                                    '_'//'SizeDistribution.srw'

                    !Size Distribution output
                    open(Unit = Species%SizeDistributionOutput%Unit(iIndexOutput), File = trim(OuputFileName), Status = 'REPLACE')


                    !time serie format
                    call WriteDataLine(Species%SizeDistributionOutput%Unit(iIndexOutput), "Time Serie Results File")

                    call WriteDataLine(Species%SizeDistributionOutput%Unit(iIndexOutput), "NAME", 'Population File')

                    call WriteDataLine(Species%SizeDistributionOutput%Unit(iIndexOutput), 'SERIE_INITIAL_DATA', Me%InitialDate)

                    call WriteDataLine(Species%SizeDistributionOutput%Unit(iIndexOutput), 'TIME_UNITS', 'SECONDS')

                    SizeDistributionHeader = 'Seconds YY MM DD hh mm ss'

                    do i = 1, Species%nSizeClasses

                        write(SizeClassNameStr, ('(F16.3)')) Species%SizeClasses(i)

                        SizeDistributionHeader = trim(SizeDistributionHeader)//' '//trim(SizeClassNameStr)

                    end do

                    write(Species%SizeDistributionOutput%Unit(iIndexOutput),102) SizeDistributionHeader

                    !write a dynamic format
                    if    (Species%nSizeClasses .lt. 10                                        )then
                        write(Species%SizeClassesCharFormat, '(i1)')Species%nSizeClasses
                    elseif(Species%nSizeClasses .ge. 10   .and. Species%nSizeClasses .lt. 100  )then
                        write(Species%SizeClassesCharFormat, '(i2)')Species%nSizeClasses
                    elseif(Species%nSizeClasses .ge. 100  .and. Species%nSizeClasses .lt. 1000 )then
                        write(Species%SizeClassesCharFormat, '(i3)')Species%nSizeClasses
                    elseif(Species%nSizeClasses .ge. 1000 .and. Species%nSizeClasses .lt. 10000)then
                        write(Species%SizeClassesCharFormat, '(i4)')Species%nSizeClasses
                    else
                        stop 'Number of volumes limited to 9999.'
                    endif

                call WriteDataLine(Species%SizeDistributionOutput%Unit(iIndexOutput), '<BeginTimeSerie>')

                end if
            
            end if !population

            Species => Species%Next
        enddo

        103 format(A50)

        if (Me%ComputeOptions%MassBalance) then

            !mass balance results
            call UnitsManager(Me%MassOutput%Unit(iIndexOutput), OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ConstructOutputs - ModuleBivalve - ERR30'

            open(Unit = Me%MassOutput%Unit(iIndexOutput), File = trim(Me%PathFileName)//trim(IndexOutputStr)// &
                            '_MassBalance.dat',&

         Status = 'REPLACE')

            if (Me%ComputeOptions%PelagicModel .eq. LifeModel) then

                write(Me%MassOutput%Unit(iIndexOutput), 103) "YY MM DD hh mm ss SumC_g SumN_g SumP_g"         

            else

                write(Me%MassOutput%Unit(iIndexOutput), 103) "Mass Balance only possible to estimate for N and P"         
                write(Me%MassOutput%Unit(iIndexOutput), 103) "Carbon is not followed when using WaterQuality model"         
                write(Me%MassOutput%Unit(iIndexOutput), 103) "YY MM DD hh mm ss SumN_g SumP_g"         

            end if
        end if 

    end subroutine ConstructOutputs

    !-------------------------------------------------------------------------------

    subroutine ConstructCohortOutput (Cohort, iIndexOutput)
        
        !Arguments------------------------------------------------------------------
        type (T_Cohort),        pointer   :: Cohort
        integer                           :: iIndexOutput


        !Local----------------------------------------------------------------------
        integer                           :: STAT_CALL

        character(len=500)                :: CohortFileName
        character(len=900)                :: OuputHeader
        character(len=16)                 :: IndexOutputStr
        
        !Begin----------------------------------------------------------------------

        call UnitsManager(Cohort%CohortOutput%Unit(iIndexOutput), OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ConstructCohortOutput - ModuleBivalve - ERR01'

        !Bivalve processes, bivalve1.dat
        
        write(IndexOutputStr, ('(I5)')) Me%IndexOutputs(iIndexOutput)
        
        CohortFileName = trim(Me%PathFileName)//trim(IndexOutputStr)//'_'//trim(Cohort%ID%Name)//'.srw'
        
        open(Unit = Cohort%CohortOutput%Unit(iIndexOutput), File = trim(CohortFileName), Status = 'REPLACE')
                                              
        !time serie format
        call WriteDataLine(Cohort%CohortOutput%Unit(iIndexOutput), "Time Serie Results File")

        call WriteDataLine(Cohort%CohortOutput%Unit(iIndexOutput), "NAME", trim(IndexOutputStr)//'_'//trim(Cohort%ID%Name))

        call WriteDataLine(Cohort%CohortOutput%Unit(iIndexOutput), 'SERIE_INITIAL_DATA', Me%InitialDate)
        
        call WriteDataLine(Cohort%CohortOutput%Unit(iIndexOutput), 'TIME_UNITS', 'SECONDS')

        101 format(A800)

        OuputHeader =  "!Seconds_1 YY_2	MM_3 DD_4 hh_5 mm_6	ss_7 "                    // &
                       "#/m2_8 mol_9 mol_10 mol_11 mol_12 cm_13y_14 m3/d.ind_15 "     // &
                       "g/d.ind_16 mol/d.ind_17 g/d.ind_18 mol/d.ind_19 g/d.ind_20 "  // &
                       "mol/d.ind_21 mol/d.ind_22 g/d.ind_23 mol/d_24 mol/d_25 "      // &
                       "mol/dm3_26 mol/dm3_27 mol/d_28 mol/d_29 mol/d_30 "            // &
                       "mol/d_31 mol/d_32 mol/d_33 #_34 mol/d_35 "                    // &
                       "mol/d_36 mol/d_37 mol/d_38 mol/d_39 mol/d_40 "                // &
                       "#/d.m2_41 #/d.m2_42 #/d.m2_43 #/d.m2_44 #/d.m2_45 "           // &
                       "#/d.m2_46 #/d.m2_47 #/d.m2_48 #/d.m2_49 #/d.m2_50 "           // &
                       "#/d.m2_51 #/d.m2_52 #/d.m2_53 #/d.m2_54 adim_55"

        write(Cohort%CohortOutput%Unit(iIndexOutput), 101) OuputHeader

        OuputHeader = "Seconds_1 YY_2 MM_3 DD_4 hh_5 mm_6 ss_7 "            // &
                      "Number_8 ME_9 MV_10 MH_11 MR_12 L_13 A_14 Cr_15 "    // &
                      "FInorg_16 F_17 IInorg_18 I_19 PFInorg_20 "           // &
                      "PF_21 Ass_22 FAEIng_23 FAE_24 JEM_25 "               // &
                      "JE_26 dE_27 GamLoss_28 StruLoss_29 JV_30 "           // &
                      "JH_31 MatLoss_32 JS_33 Gam_34 JR_35 "                // &
                      "CO2_36 H2O_37 O2_38 NH3_39 PO4_40 "                  // &
                      "m_A_41 m_O_42 m_F_43 m_nat_44 m_shr_45 "             // &
                      "m_cra_46 m_oys_47 m_duck_48 m_gull_49 m_low_50 "     // &
                      "m_self_51 m_others_52 m_vel_53 m_settle_54 ScaledE"
        
        write(Cohort%CohortOutput%Unit(iIndexOutput), 101) OuputHeader
        
        call WriteDataLine(Cohort%CohortOutput%Unit(iIndexOutput), '<BeginTimeSerie>')

    end subroutine ConstructCohortOutput

    !-------------------------------------------------------------------------------

    subroutine ConstructSpeciesParameters (NewSpecies, ClientNumber)

        !Arguments------------------------------------------------------------------
        type (T_Species),      pointer              :: NewSpecies
        integer                                     :: ClientNumber

        !External-------------------------------------------------------------------
        integer                                     :: flag, STAT_CALL, i
        integer                                     :: FirstLine, Line, LastLine
        logical                                     :: BlockLayersFound
        character(len = PathLength)                 :: ArgumentInComand 
        !Begin-----------------------------------------------------------------
        
        !Name of Species
        call GetData(NewSpecies%ID%Name             , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'NAME'           , &
                    !default      = ' '             , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)
                    
        if (STAT_CALL .NE. SUCCESS_)                  &
        stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR01'

        if(.not. Checkpropertyname(trim(NewSpecies%ID%Name), NewSpecies%ID%IDNumber))then
            write(*,*)trim(NewSpecies%ID%Name)
            stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR10'
        end if 

        !Description of Species
        call GetData(NewSpecies%ID%Description      , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'DESCRIPTION'    , &
                    !default      = ' '             , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                  &
        stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR20'

        !Description of Species
        call GetData(NewSpecies%Population          , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'POPULATION'     , &
                    default      = .false.          , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)

        !Description of Species
        call GetData(NewSpecies%FeedOnLarvae        , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'FEED_ON_LARVAE' , &
                    default      = .false.          , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                  &
        stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR30'

        !Description of Species
        call GetData(NewSpecies%LarvaeMaxSize       , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'LARVAE_MAXSIZE' , &
                    default      = 0.026            , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)      &
        stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR40'

        !Use an extra starvation mortality?
        call GetData(NewSpecies%ExtraStarvation     , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'M_STARVATION'   , &
                    default      = .false.          , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)      &
        stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR41'
        

        !Number of Cohorts from this species
        call GetData(NewSpecies%Initial_nCohorts       , &
                    Me%ObjEnterData, flag              , &
                    SearchType   = FromBlock           , &
                    keyword      = 'NUMBER_OF_COHORTS' , &
                    !default     = ' '                  , &
                    ClientModule = 'ModuleBivalve'     , &
                    STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)         &
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR50'

        if(flag==0)then
        write(*,*)"NUMBER_OF_COHORTS must be defined for bivalve species : ", trim(adjustl(NewSpecies%ID%Name))
        stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR051'
        endif
        
        !Minimum observed length
        call GetData(NewSpecies%MinObsLength        , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'MIN_OBS_LENGTH' , &
                    default      = 0.1              , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)      &
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR060'

        !Option to write Cohort output file
        call GetData(NewSpecies%CohortOutput        , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'COHORT_OUTPUT'  , &
                    default      = .true.           , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)      &
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR70'

        !Option to write Size Distribution output file
        call GetData(NewSpecies%BySizeOutput        , &
                    Me%ObjEnterData, flag           , &
                    SearchType   = FromBlock        , &
                    keyword      = 'BYSIZE_OUTPUT'  , &
                    default      = .true.           , &
                    ClientModule = 'ModuleBivalve'  , &
                    STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)      &
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR80'

        call GetData(NewSpecies%Testing_File          , &
                    Me%ObjEnterData, flag             , &
                    SearchType   = FromFile           , &
                    keyword      = 'TESTING_FILENAME' , &
                    !default      = .false.            , &
                    ClientModule = 'ModuleBivalve'    , &
                    STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)        &        
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR90'  
                    
        if ((flag .eq. 0) .and. (Me%Testing_Parameters)) then 
        
            !call getarg(1,ArgumentInComand)
            
            NewSpecies%Testing_File = trim(Me%PathFileName)// &
                            trim(ArgumentInComand)//"_ParameterTest.dat"
            
        end if
        
        !Step to define size classes in the Size Distribution output file
        call GetData(NewSpecies%SizeStep              , &
                    Me%ObjEnterData, flag             , &
                    SearchType   = FromBlock          , &
                    keyword      = 'SIZE_STEP'        , &
                    default      = 2.0                , &
                    ClientModule = 'ModuleBivalve'    , &
                    STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)        &
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR100'

        !Step to define size classes in the Size Distribution output file
        call GetData(NewSpecies%Max_SizeClass         , &
                    Me%ObjEnterData, flag             , &
                    SearchType   = FromBlock          , &
                    keyword      = 'MAX_SIZECLASS'    , &
                    default      = 18.0               , &
                    ClientModule = 'ModuleBivalve'    , &
                    STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)        &
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR110'

        if(flag == 0)then

            call ExtractBlockFromBlock( Me%ObjEnterData,            & 
                    ClientNumber,               &
                    '<<begin_size_classes>>',   &
                    '<<end_size_classes>>',     &
                    BlockLayersFound,           &
                    FirstLine = FirstLine,      &
                    LastLine  = LastLine,       &
                    STAT      = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_)then 
               
                if (BlockLayersFound) then

                    NewSpecies%nSizeClasses = (LastLine - FirstLine - 1)
                    allocate(NewSpecies%SizeClasses(1:NewSpecies%nSizeClasses))  !To sort individuals by size class
                    NewSpecies%SizeClasses = 0.0

                    allocate(NewSpecies%SizeFrequency(1:NewSpecies%nSizeClasses))  !To sort individuals by size class
                    NewSpecies%SizeFrequency = 0.0

                    Line = FirstLine + 1

                    do  i = 1, NewSpecies%nSizeClasses

                        call GetData(NewSpecies%SizeClasses(i), Me%ObjEnterData,    &
                        flag, Buffer_Line  = Line, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR120'

                        Line = Line + 1

                    enddo

                else
                    write(*,*)'Block <<begin_size_classes>> <<end_size_classes>> not found'
                    stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR130'
                endif

        else

            stop 'Subroutine ConstructSpeciesParameters - ModuleBivalve - ERR140'

        endif

    else

        NewSpecies%nSizeClasses = int(NewSpecies%Max_SizeClass/NewSpecies%SizeStep) + 1

        allocate(NewSpecies%SizeClasses(1:NewSpecies%nSizeClasses))  !To sort individuals by size class
        NewSpecies%SizeClasses = 0.0

        allocate(NewSpecies%SizeFrequency(1:NewSpecies%nSizeClasses))  !To sort individuals by size class
        NewSpecies%SizeFrequency = 0.0

        do i = 2, NewSpecies%nSizeClasses

            NewSpecies%SizeClasses(i) = NewSpecies%SizeClasses(i-1) + NewSpecies%SizeStep

        end do

    endif
        
        allocate(NewSpecies%PopulationProcesses%SumLogAllMortalityInNumbers(14)) 
        allocate(NewSpecies%PopulationProcesses%SumAllMortalityInMass(14))      
        allocate(NewSpecies%PopulationProcesses%AverageMortalityInNumbers(14))   
        
        NewSpecies%PopulationProcesses%SumLogAllMortalityInNumbers = 0.0
        NewSpecies%PopulationProcesses%SumAllMortalityInMass       = 0.0
        NewSpecies%PopulationProcesses%AverageMortalityInNumbers   = 0.0
        NewSpecies%PopulationProcesses%nInstantsForAverage         = 0.0       
        NewSpecies%PopulationProcesses%nSpawning                   = 0
        
        call ConstructSpeciesComposition    (NewSpecies%SpeciesComposition   )
        call ConstructIndividualParameters  (NewSpecies%IndividualParameters )
        call ConstructParticles             (NewSpecies, ClientNumber        )
        call ConstructPredator              (NewSpecies, ClientNumber        )
        
        allocate(NewSpecies%SettlementProbability(Me%Array%ILB:Me%Array%IUB))

    end subroutine ConstructSpeciesParameters

    !-------------------------------------------------------------------------------

    subroutine ConstructSpeciesComposition (SpeciesComposition)

        !Arguments------------------------------------------------------------------
        type (T_SpeciesComposition)         :: SpeciesComposition

        !External-------------------------------------------------------------------
        integer :: flag, STAT_CALL

        !Local----------------------------------------------------------------------
        integer :: FromBlock 

        !Begin----------------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        !molH/molC, chemical index of hydrogen in bivalve reserves
        call GetData(SpeciesComposition%ReservesComposition%nC       , &
                    Me%ObjEnterData, flag                            , &
                    SearchType   = FromBlock                         , &
                    keyword      = 'RESERVES_nC'                     , &
                    default      = 1.                                , &
                    ClientModule = 'ModuleBivalve'                   , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                       &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR00'


        !molH/molC, chemical index of hydrogen in bivalve reserves, (Kooijman, 2010)
        call GetData(SpeciesComposition%ReservesComposition%nH       , &
                    Me%ObjEnterData, flag                            , &
                    SearchType   = FromBlock                         , &
                    keyword      = 'RESERVES_nH'                     , &
                    default      = 1.8                               , &
                    ClientModule = 'ModuleBivalve'                   , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                       &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR10'

        !molO/molC, chemical index of oxygen in bivalve reserves, (Kooijman, 2010)
        call GetData(SpeciesComposition%ReservesComposition%nO       , &
                    Me%ObjEnterData, flag                            , &
                    SearchType   = FromBlock                         , &
                    keyword      = 'RESERVES_nO'                     , &
                    default      = 0.53                              , &
                    ClientModule = 'ModuleBivalve'                   , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                       &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR20'

        !molN/molC, chemical index of nitrogen in bivalve reserves, (Smaal and Vonck, 1997)
        call GetData(SpeciesComposition%ReservesComposition%nN       , &
                    Me%ObjEnterData, flag                            , &
                    SearchType   = FromBlock                         , &
                    keyword      = 'RESERVES_nN'                     , &
                    default      = 0.15                              , &
                    ClientModule = 'ModuleBivalve'                   , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                       &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR30'

        !molN/molC, chemical index of phosphorus in bivalve reserves, (Smaal and Vonck, 1997)
        call GetData(SpeciesComposition%ReservesComposition%nP       , &
                    Me%ObjEnterData, flag                            , &
                    SearchType   = FromBlock                         , &
                    keyword      = 'RESERVES_nP'                     , &
                    default      = 0.006                             , &
                    ClientModule = 'ModuleBivalve'                   , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                             &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR40'


        !molC/molC, chemical index of hydrogen in bivalve structure
        call GetData(SpeciesComposition%StructureComposition%nC      , &
                    Me%ObjEnterData, flag                            , &
                    SearchType   = FromBlock                         , &
                    keyword      = 'STRUCTURE_nC'                    , &
                    default      =  1.                               , &
                    ClientModule = 'ModuleBivalve'                   , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                       &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR50'


        !molH/molC, chemical index of hydrogen in bivalve structure, (Kooijman, 2010)
        call GetData(SpeciesComposition%StructureComposition%nH     , &
                    Me%ObjEnterData, flag                           , &
                    SearchType   = FromBlock                        , &
                    keyword      = 'STRUCTURE_nH'                   , &
                    default      =  1.8                             , &
                    ClientModule = 'ModuleBivalve'                  , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                      &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR60'

        !molO/molC, chemical index of oxygen in bivalve structure, (Kooijman, 2010)
        call GetData(SpeciesComposition%StructureComposition%nO     , &
                    Me%ObjEnterData, flag                           , &
                    SearchType   = FromBlock                        , &
                    keyword      = 'STRUCTURE_nO'                   , &
                    default      = 0.53                             , &
                    ClientModule = 'ModuleBivalve'                  , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                      &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR70'

        !molN/molC, chemical index of nitrogen in bivalve structure, (Smaal and Vonck, 1997)
        call GetData(SpeciesComposition%StructureComposition%nN     , &
                    Me%ObjEnterData, flag                           , &
                    SearchType   = FromBlock                        , &
                    keyword      = 'STRUCTURE_nN'                   , &
                    default      = 0.15                             , &
                    ClientModule = 'ModuleBivalve'                  , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                      &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR80'

        !molN/molC, chemical index of phosphorus in bivalve structure, (Smaal and Vonck, 1997)
        call GetData(SpeciesComposition%StructureComposition%nP     , &
                    Me%ObjEnterData, flag                           , &
                    SearchType   = FromBlock                        , &
                    keyword      = 'STRUCTURE_nP'                   , &
                    default      = 0.006                            , &
                    ClientModule = 'ModuleBivalve'                  , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)                      &
                    stop 'Subroutine ConstructSpeciesComposition - ModuleBivalve - ERR90'

    end subroutine ConstructSpeciesComposition

    !-------------------------------------------------------------------------------

    subroutine ConstructIndividualParameters (IndividualParameters)

        !Arguments------------------------------------------------------------------
        type (T_IndividualParameters)       :: IndividualParameters

        !External-------------------------------------------------------------------
        integer :: flag, STAT_CALL

        !Local----------------------------------------------------------------------
        integer :: FromBlock 

        !Begin----------------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        !K, Rate Temperature reference
        call GetData(IndividualParameters%Tref              , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'Tref'                   , &
                    default      = 293.0                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR01'


        !K, Arrhenius temperature (van der Veer etal., 2006)
        call GetData(IndividualParameters%TA                , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'TA'                     , &
                    default      = 7022.0                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR02'

        !K, Lower Boundary tolerance range  (van der Veer etal., 2006)
        call GetData(IndividualParameters%TL                , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'TL'                     , &
                    default      = 273.0                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR03'

        !K, Upper Boundary tolerance rang  (van der Veer etal., 2006)
        call GetData(IndividualParameters%TH                , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'TH'                     , &
                    default      = 290.0                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR04'

        !K, Arrhenius temperature for lower boundary  (van der Veer etal., 2006)
        call GetData(IndividualParameters%TAL               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'TAL'                    , &
                    default      = 45430.0                  , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR05'

        !K, Arrhenius temperature for upper boundary  (van der Veer etal., 2006)
        call GetData(IndividualParameters%TAH               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'TAH'                    , &
                    default      = 31376.0                  , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR06'

        !adim, constant food density parameter (only if simple filtration)
        call GetData(IndividualParameters%F_FIX             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'F_FIX'                  , &
                    default      = 1.                       , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR07'

        !Jd-1cm-2, bivalve surface-specific assimilation rate if fix (Saraiva etal., inpress)   
        call GetData(IndividualParameters%PAM_FIX           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'PAM_FIX'                , &
                    default      = 80.5                     , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                if (STAT_CALL .NE. SUCCESS_)                  &
                stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR08'

        !cm(volumetric)/cm(real), shape coefficient  (Saraiva etal., inpress)
        call GetData(IndividualParameters%delta_M           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'DELTA_M'                , &
                    default      = 0.297                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR09'

        !years, life span of a mussel, Sukhotin et al. (2007)
        call GetData(IndividualParameters%LifeSpan          , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'LIFE_SPAN'              , &
                    default      = 24.0                     , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR10'

        !fraction, natural mortality rate
        call GetData(IndividualParameters%m_natural         , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'M_NATURAL'              , &
                    default      = 0.0                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR11'

        !%, spat mortality rate
        call GetData(IndividualParameters%m_spat            , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'M_SPAT'                 , &
                    default      = 0.0                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR12'

        !cm/d, energy conductance (Saraiva etal., in press)    
        call GetData(IndividualParameters%v_cond            , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'V_COND'                 , &
                    default      = 0.056                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR013'

        !adim, allocation fraction to growth/somatic maintenace (Saraiva etal., inpress) 
        call GetData(IndividualParameters%kappa             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'KAPPA'                  , &
                    default      = 0.67                     , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR14'

        !adim, fraction of flux allocated to reproduction (Kooijman, 2010)
        call GetData(IndividualParameters%kap_R             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'KAP_R'                  , &
                    default      = 0.95                     , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR15'

        !J/(d.cm3), volume specific somatic maintenace energy flux (Saraiva etal., inpress) 
        call GetData(IndividualParameters%pM                , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'pM'                     , &
                    default      = 11.6                     , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR16'

        !J/cm3(volumetric), energy costs for structural volume growth (Saraiva etal., inpress)
        call GetData(IndividualParameters%EG                , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'EG'                     , &
                    default      = 5993.                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR17'

        !J, Maturity threshold for birth (Saraiva etal., inpress)
        call GetData(IndividualParameters%EHb               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'EH_B'                   , &
                    default      = 2.99e-5                  , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR18'

        !J, Maturity threshold for puberty (Saraiva etal., inpress)
        call GetData(IndividualParameters%EHp               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'EH_P'                   , &
                    default      = 1.58e2                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR19'

        !m3/d.cm2, maximum clearance rate (Saraiva etal., 2011) 
        call GetData(IndividualParameters%Crm               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'CRM'                    , &
                    default      = 0.096                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR20'

        !molC/(d.cm2), algae maximum surface area-specific filtration rate (Thomas etal., 2011) 
        call GetData(IndividualParameters%JX1Fm             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'JX1FM'                  , &
                    default      = 4.8e-4                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR21'

        !g/(d.cm2), inorganic material maximum surface area-specific filtration rate (Saraiva etal., 2011) 
        call GetData(IndividualParameters%JX0Fm             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'JX0FM'                  , &
                    default      = 3.5                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR22'

        !adim, algae binding probability (Saraiva etal., inpress)
        call GetData(IndividualParameters%ro_X1             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'RO_XI'                  , &
                    default      = 0.4                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR23'

        !adim, inorganic material binding probability (Saraiva etal., inpress)
        call GetData(IndividualParameters%ro_X0             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'RO_X0'                  , &
                    default      = 0.4                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR24'

        !molC/(d.cm2), algae maximum surface area-specific ingestion rate (Saraiva etal., 2011) 
        call GetData(IndividualParameters%JX1Im             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'JX1IM'                  , &
                    default      = 1.3e4                    , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR25'

        !g/(d.cm2), inorganic material maximum surface area-specific ingestion rate (Saraiva etal., 2011) 
        call GetData(IndividualParameters%JX0Im             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'JX0IM'                  , &
                    default      = 0.11                     , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR26'

        !molC(res)/molC(food), yield coeficienct of reserves in algae structure (Kooijman, 2010)        
        call GetData(IndividualParameters%YEX               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'YEX'                    , &
                    default      = 0.75                     , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR27'

        !molC(gam)/molC(struc), minimum gonado-somatic ratio in the organism (Cardoso et al., 2007)
        call GetData(IndividualParameters%GSR_MIN           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'GSR_MIN'                , &
                    default      = 0.1                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR28'

        !molC(gam)/molC(struc), gonado-somatic ratio to spawn (Saraiva etal., submited)  
        call GetData(IndividualParameters%GSR_SPAWN         , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'GSR_SPAWN'              , &
                    default      = 0.2                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR29'

        !oC, minimum temperature for spawning (Hummel etal., 1989)
        call GetData(IndividualParameters%T_SPAWN           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'T_SPAWN'                , &
                    default      = 9.6                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR30'

        !d, minimum time interval between spawning events    
        call GetData(IndividualParameters%MIN_SPAWN_TIME    , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'MIN_SPAWN_TIME'         , &
                    default      = 0.                       , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR31'



        !molC(reser), reserves in an embryo at optimal food conditions (Saraiva etal., submited)    
        call GetData(IndividualParameters%ME_0              , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'ME_0'                   , &
                    default      = 1.48e-10                 , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)               &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR32'



        !molC, reserves in a new born individual at optimal food conditions (Saraiva etal., submited) 
        call GetData(IndividualParameters%MEb               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'ME_B'                   , &
                    default      = 1.0e-7                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR33'


        !molC, structure in a new born individual at optimal food conditions (Saraiva etal., submited) 
        call GetData(IndividualParameters%MVb               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'MV_B'                   , &
                    default      = 7.52e-11                 , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR34'

        !molC, maturity in a new born individual at optimal food conditions (Saraiva etal., submited) 
        call GetData(IndividualParameters%MHb               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'MH_B'                   , &
                    default      = 4.24e-11                 , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR35'

        !molC, length in a new born individual at optimal food conditions (Saraiva etal., submited) 
        call GetData(IndividualParameters%Lb                , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'L_B'                    , &
                    default      = 7.3e-3                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR36'

        !g(dw)/cm3, bivalve structure and reserves specific density (Rosland etal., 2009 and Brey, 2001)
        call GetData(IndividualParameters%d_V               , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'DV'                     , &
                    default      = 0.2                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR37'

        !J/molC(reser), chemical potential of reserves (van der Veer etal., 2006)
        call GetData(IndividualParameters%mu_E              , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'MU_E'                   , &
                    default      = 6.97e5                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR38'

        !option to compute simple assimilation 
        call GetData(IndividualParameters%SIMPLE_ASSI       , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'SIMPLE_ASSI'            , &
                    default      = 0                        , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR39'

        !option to compute simple temperature correction factor 
        call GetData(IndividualParameters%SIMPLE_TEMP       , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'SIMPLE_TEMP'            , &
                    default      = 1                        , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR40'

        !/d, fraction of individuals that die due to high velocity 
        call GetData(IndividualParameters%m_velocity        , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'M_VELOCITY'             , &
                    default      = 0.0                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR41'

        !m/s, maximum  water velocity tolerable for this species
        call GetData(IndividualParameters%MAX_velocity      , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'MAX_VELOCITY'           , &
                    default      = 0.5                      , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR42'
                    
        !Em, maximum  reserves capacity, J/cm2
        call GetData(IndividualParameters%Em                , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'E_M'                    , &
                    default      = 1438.0                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR43'
                    
        !logical, assume dnesity limitation? 
        call GetData(IndividualParameters%DensityLimOption  , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'DENSITYLIMIT'           , &
                    default      = .false.                  , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR44'

        !#/m2, maximum density found in the observations 
        call GetData(IndividualParameters%MAX_density       , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlock                , &
                    keyword      = 'DENSITY_MAXVALUE'       , &
                    default      = 3000.0                   , &
                    ClientModule = 'ModuleBivalve'          , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructIndividualParameters - ModuleBivalve - ERR45'

    end subroutine ConstructIndividualParameters

    !-------------------------------------------------------------------------------

    subroutine ConstructParticles (NewSpecies, ClientNumber)

        !Arguments------------------------------------------------------------------
        type (T_Species)  , pointer                 :: NewSpecies
        integer                                     :: ClientNumber

        !External-------------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local----------------------------------------------------------------------
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockInBlockFound
        type (T_Particles) , pointer                :: NewParticles

        !Begin----------------------------------------------------------------------

do1 :   do
        call ExtractBlockFromBlock(Me%ObjEnterData                       , &
                                ClientNumber      = ClientNumber         , &
                                block_begin       = '<<begin_particle>>' , &
                                block_end         = '<<end_particle>>'   , &
                                BlockInBlockFound = BlockInBlockFound    , &
                                FirstLine         = FirstLine            , &
                                LastLine          = LastLine             , &
                                STAT= STAT_CALL)

cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockInBlockFound) then        

                    allocate(NewParticles)

                    call AddParticles (NewSpecies, NewParticles)

                    call ConstructParticlesParameters (NewParticles)

                    nullify(NewParticles)

                else cd2
                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Subroutine ConstructParticles - ModuleBivalve - ERR02'
            else cd1
                stop 'Subroutine ConstructParticles - ModuleBivalve - ERR03'
            end if cd1
        end do do1

    end subroutine ConstructParticles

    !-------------------------------------------------------------------------------

    subroutine AddParticles (ObjSpecies, NewParticles)

        !Arguments-------------------------------------------------------------
        type (T_Species),            pointer  :: ObjSpecies
        type (T_Particles),          pointer  :: NewParticles

        !Local-----------------------------------------------------------------
        type (T_Particles),          pointer  :: ObjParticles
        integer, save                         :: NextParticlesID = 1

        nullify  (NewParticles%Next)

        !Insert new Food into list
cd1:    if (.not. associated(ObjSpecies%FirstParticles)) then
            NextParticlesID= 1
            ObjSpecies%FirstParticles    => NewParticles
        else

            ObjParticles => ObjSpecies%FirstParticles

do1:        do while (associated(ObjParticles%Next))
                ObjParticles => ObjParticles%Next
            enddo do1

            ObjParticles%Next => NewParticles
        endif cd1

        !Attributes ID
        NewParticles%ID%ID            = NextParticlesID

        NextParticlesID = NextParticlesID + 1

        ObjSpecies%nParticles         = ObjSpecies%nParticles + 1

    end subroutine AddParticles

    !-------------------------------------------------------------------------------

    subroutine ConstructParticlesParameters (NewParticles)

        !Arguments------------------------------------------------------------------
        type (T_Particles),      pointer            :: NewParticles

        !External-------------------------------------------------------------------
        integer                                     :: flag, STAT_CALL
        integer                                     :: FromBlockInBlock
        !Begin----------------------------------------------------------------------


        call GetExtractType    (FromBlockInBlock = FromBlockInBlock)


        call GetData(NewParticles%ID%Name                   , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'NAME'                   , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR01'


        call GetData(NewParticles%ID%Description            , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'DESCRIPTION'            , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR02'

        !Organic?, Is this property organic?  
        call GetData(NewParticles%Organic                   , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'ORGANIC'                , &
                    default      = 0                        , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR03'

        !Silica?, Just if the pelagic model is life 
        call GetData(NewParticles%Silica                    , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'SILICA_USE'             , &
                    default      = 0                        , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR04'


        !The NC and PC ratios of this property are variable?
        call GetData(NewParticles%RatioVariable             , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'RATIO_VARIABLE'         , &
                    default      =  0                       , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR05'

        !Ratio H/C [mgH/mgC]
        call GetData(NewParticles%Ratios%HC_Ratio           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'RATIOHC'                , &
                    default      = 0.15                     , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR06'      

        !Ratio O/C [mgO/mgC]
        call GetData(NewParticles%Ratios%OC_Ratio           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'RATIOOC'                , &
                    default      = 0.67                     , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              & 
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR07'      

        !Ratio N/C [mgN/mgC]
        call GetData(NewParticles%Ratios%NC_Ratio           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'RATIONC'                , &
                    default      = 0.16                     , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR08'      

        !Ratio P/C [mgP/mgC]
        call GetData(NewParticles%Ratios%PC_Ratio           , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'RATIOPC'                , &
                    default      = 0.024                    , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR09'      

        !Ratio Si/C [mgSi/mgC]
        call GetData(NewParticles%Ratios%SiC_Ratio          , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'RATIOSiC'               , &
                    default      = 0.89                     , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR10'      

        !Ratio Chla/C [mgChla/mgC]
        call GetData(NewParticles%Ratios%ChlC_Ratio         , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'RATIOCHLC'              , &
                    default      = 0.017                    , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR11'      

        !cm, mean size of property
        call GetData(NewParticles%Size                      , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'SIZE'                   , &
                    default      = 0.2                      , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR12'      

        !molCReserves/molCTotalFood, fraction of reserves in the food    
        call GetData(NewParticles%f_E                       , &
                    Me%ObjEnterData, flag                   , &
                    SearchType   = FromBlockInBlock         , &
                    keyword      = 'F_E'                    , &
                    default      = 0.5                      , &
                    ClientModule = 'ModuleBivalve'          , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)              &
                    stop 'Subroutine ConstructParticlesParameters - ModuleBivalve - ERR13'


    end subroutine ConstructParticlesParameters

    !-------------------------------------------------------------------------------

    subroutine ConstructPredator (NewSpecies, ClientNumber)

        !Arguments------------------------------------------------------------------
        type (T_Species)             , pointer            :: NewSpecies
        integer                                           :: ClientNumber

        !External-------------------------------------------------------------------
        integer                                           :: STAT_CALL

        !Local----------------------------------------------------------------------
        integer                                           :: FirstLine, LastLine
        logical                                           :: BlockInBlockFound
        type (T_Predator)            , pointer            :: NewPredator

        !Begin----------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBlock(Me%ObjEnterData                       , &
                                    ClientNumber      = ClientNumber         , &
                                    block_begin       = '<<begin_predator>>' , &
                                    block_end         = '<<end_predator>>'   , &
                                    BlockInBlockFound = BlockInBlockFound    , &
                                    FirstLine         = FirstLine            , &
                                    LastLine          = LastLine             , &
                                    STAT= STAT_CALL)

cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockInBlockFound) then        

                    allocate(NewPredator)

                    call AddPredator (NewSpecies, NewPredator)

                    call ConstructPredatorParameters (NewPredator)

                    nullify(NewPredator)

                else cd2
                !
                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Subroutine ConstructPredator - ModuleBivalve - ERR02'
            else cd1
                stop 'Subroutine ConstructPredator - ModuleBivalve - ERR03'
            end if cd1
        end do do1

    end subroutine ConstructPredator

    !-------------------------------------------------------------------------------

    subroutine AddPredator (ObjSpecies, NewPredator)

        !Arguments-------------------------------------------------------------
        type (T_Species),            pointer:: ObjSpecies
        type (T_Predator),          pointer :: NewPredator

        !Local-----------------------------------------------------------------
        type (T_Predator),          pointer :: ObjPredator
        integer, save                       :: NextPredatorID = 1

        nullify  (NewPredator%Next)

        !Insert new Food into list
cd1:    if (.not. associated(ObjSpecies%FirstPredator)) then
            NextPredatorID= 1
            ObjSpecies%FirstPredator    => NewPredator
        else

            ObjPredator => ObjSpecies%FirstPredator

do1:        do while (associated(ObjPredator%Next))
                ObjPredator => ObjPredator%Next
            enddo do1

            ObjPredator%Next => NewPredator
        endif cd1

        !Attributes ID
        NewPredator%ID%ID     = NextPredatorID

        NextPredatorID        = NextPredatorID + 1

        ObjSpecies%nPredator  = ObjSpecies%nPredator + 1

    end subroutine AddPredator

    !-------------------------------------------------------------------------------

    subroutine ConstructPredatorParameters (NewPredator)

        !Arguments------------------------------------------------------------------
        type (T_Predator),      pointer           :: NewPredator

        !External-------------------------------------------------------------------
        integer                                   :: flag, STAT_CALL
        integer                                   :: FromBlockInBlock
        !Begin----------------------------------------------------------------------


        call GetExtractType    (FromBlockInBlock = FromBlockInBlock)

        call GetData(NewPredator%ID%Name                   , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'NAME'                  , &
                    !default      = ' '                    , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR01'


        call GetData(NewPredator%ID%Description            , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'DESCRIPTION'           , &
                    !default      = ' '                    , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR02'

        !cm, predator size
        call GetData(NewPredator%PredatorSize              , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'SIZE'                  , &
                    default      = 0.0                     , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR03'


        !cm, minimum size of the prey
        call GetData(NewPredator%MinPreySize               , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'MINPREYSIZE'           , &
                    default      = 0.                      , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR04'

        !cm, maximum size of the prey
        call GetData(NewPredator%MaxPreySize               , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'MAXPREYSIZE'           , &
                    default      = 0.                      , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR05'

        !the units depend on the predator, Feeding_Units
        call GetData(NewPredator%Feeding_Rate              , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'FEEDING_RATE'          , &
                    default      =  0.                     , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR06'      

        !1-#/d.ind; 2-AFDW/d.ind, 3-J/cm2.d
        call GetData(NewPredator%Feeding_Units             , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'FEEDING_UNITS'         , &
                    default      =  1                      , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR07'      

        !1-Always; 2-LowTide, 3-HighTide
        call GetData(NewPredator%Feeding_Time              , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'FEEDING_TIME'          , &
                    default      =  1                      , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR08'      

        !-, fraction of mussels in the shrimp food (this study)
        call GetData(NewPredator%Diet                      , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'DIET_FRACTION'         , &
                    default      = 0.15                    , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR09'      

        !-, conversion of afdw in dw
        call GetData(NewPredator%AfdwToC                   , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'AFDW_DW'               , &
                    default      = 0.85                    , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR10'      

        !-, conversion of afdw in dw
        call GetData(NewPredator%DwToC                     , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlockInBlock        , &
                    keyword      = 'DW_C'                  , &
                    default      = 2.5                     , &
                    ClientModule = 'ModuleBivalve'         , &
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR11'      

        !option to compute simple temperature correction factor 
        call GetData(NewPredator%SIMPLE_TEMP               , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'SIMPLE_TEMP'           , &
                    default      = 1                       , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR12'

        !option to compute simple temperature correction factor 
        call GetData(NewPredator%CORRECT_TEMP              , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'CORRECT_TEMP'          , &
                    default      = 1                       , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR12'

        !K, Rate Temperature reference for predator
        call GetData(NewPredator%P_Tref                    , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'P_Tref'                , &
                    default      = 293.                    , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR13'

        !K, Arrhenius temperature for predator 
        call GetData(NewPredator%P_TA                      , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'P_TA'                  , &
                    default      = 0.                      , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR14'

        ! K, Lower Boundary tolerance range for predator 
        call GetData(NewPredator%P_TL                      , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'P_TL'                  , &
                    default      = 0.                      , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR15'


        !K, Upper Boundary tolerance range for predator 
        call GetData(NewPredator%P_TH                      , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'P_TH'                  , &
                    default      = 0.                      , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR16'

        !K, Arrhenius temperature for lower boundary for predator
        call GetData(NewPredator%P_TAL                     , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'P_TAL'                 , &
                    default      = 0.                      , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR17'

        !K, Arrhenius temperature for upper boundary for predator
        call GetData(NewPredator%P_TAH                     , &
                    Me%ObjEnterData, flag                  , &
                    SearchType   = FromBlock               , &
                    keyword      = 'P_TAH'                 , &
                    default      = 0.                      , &
                    ClientModule = 'ModuleBivalve'         , & 
                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_)             &
                    stop 'Subroutine ConstructPredatorParameters - ModuleBivalve - ERR18'

    end subroutine ConstructPredatorParameters

    !-------------------------------------------------------------------------------

    subroutine PropertyIndexNumber

        !Local-----------------------------------------------------------------
        type(T_Species)   , pointer                :: Species
        type(T_Particles) , pointer                :: Particles
        type(T_Predator)  , pointer                :: Predator
        type(T_Cohort)    , pointer                :: Cohort

        !Begin-----------------------------------------------------------------

        Me%Prop%ILB = 1
        Me%Prop%IUB = 0

        !Sediments             
        Me%PropIndex%sediments              = null_int  
        !Nitrogen             
        Me%PropIndex%AM                     = null_int
        Me%PropIndex%PON                    = null_int   
        !Phosphorus           
        Me%PropIndex%IP                     = null_int    
        Me%PropIndex%POP                    = null_int   
        !Carbon 
        Me%PropIndex%POC                    = null_int   
        Me%PropIndex%CarbonDioxide          = null_int   
        !Oxygen 
        Me%PropIndex%Oxygen                 = null_int
        !Silica             
        Me%PropIndex%BioSilica              = null_int
        Me%PropIndex%DissSilica             = null_int
        !Food  
        Me%PropIndex%phyto                  = null_int
        Me%PropIndex%diatoms                = null_int
        Me%PropIndex%zoo                    = null_int
        Me%PropIndex%ciliate                = null_int
        Me%PropIndex%bacteria               = null_int
        Me%PropIndex%silica                 = null_int

        Me%PropIndex%DiatomsC               = null_int
        Me%PropIndex%DiatomsN               = null_int
        Me%PropIndex%DiatomsP               = null_int
        Me%PropIndex%DiatomsChl             = null_int
        Me%PropIndex%DiatomsSi              = null_int
        Me%PropIndex%Mix_FlagellateC        = null_int
        Me%PropIndex%Mix_FlagellateN        = null_int
        Me%PropIndex%Mix_FlagellateP        = null_int
        Me%PropIndex%Mix_FlagellateChl      = null_int
        Me%PropIndex%PicoalgaeC             = null_int
        Me%PropIndex%PicoalgaeN             = null_int
        Me%PropIndex%PicoalgaeP             = null_int
        Me%PropIndex%PicoalgaeChl           = null_int
        Me%PropIndex%FlagellateC            = null_int
        Me%PropIndex%FlagellateN            = null_int
        Me%PropIndex%FlagellateP            = null_int
        Me%PropIndex%FlagellateChl          = null_int
        Me%PropIndex%MicrozooplanktonC      = null_int
        Me%PropIndex%MicrozooplanktonN      = null_int
        Me%PropIndex%MicrozooplanktonP      = null_int
        Me%PropIndex%Het_NanoflagellateC    = null_int
        Me%PropIndex%Het_NanoflagellateN    = null_int
        Me%PropIndex%Het_NanoflagellateP    = null_int
        Me%PropIndex%MesozooplanktonC       = null_int
        Me%PropIndex%MesozooplanktonN       = null_int
        Me%PropIndex%MesozooplanktonP       = null_int
        Me%PropIndex%Het_BacteriaC          = null_int
        Me%PropIndex%Het_BacteriaN          = null_int
        Me%PropIndex%Het_BacteriaP          = null_int
        Me%PropIndex%Shrimp                 = null_int
        Me%PropIndex%Crab                   = null_int
        Me%PropIndex%OysterCatcher          = null_int
        Me%PropIndex%EiderDuck              = null_int
        Me%PropIndex%HerringGull            = null_int

        if (Me%ComputeOptions%PelagicModel .eq. LifeModel) then

            if(Me%PropIndex%POC .eq. null_int)then

                Me%Prop%IUB        = Me%Prop%IUB + 1
                Me%PropIndex%POC   = Me%Prop%IUB

                Me%Prop%IUB        = Me%Prop%IUB + 1
                Me%PropIndex%AM    = Me%Prop%IUB

                Me%Prop%IUB        = Me%Prop%IUB + 1
                Me%PropIndex%PON   = Me%Prop%IUB

                Me%Prop%IUB        = Me%Prop%IUB + 1
                Me%PropIndex%IP    = Me%Prop%IUB

                Me%Prop%IUB        = Me%Prop%IUB + 1
                Me%PropIndex%POP   = Me%Prop%IUB

            end if
        end if

        if(Me%ComputeOptions%Nitrogen)then

            if(Me%PropIndex%AM .eq. null_int)then

                Me%Prop%IUB         = Me%Prop%IUB + 1
                Me%PropIndex%AM     = Me%Prop%IUB
            end if

            if(Me%PropIndex%PON .eq. null_int)then

                Me%Prop%IUB         = Me%Prop%IUB + 1
                Me%PropIndex%PON    = Me%Prop%IUB
            end if

        end if

        if(Me%ComputeOptions%Phosphorus)then

            if(Me%PropIndex%IP .eq. null_int)then

                Me%Prop%IUB         = Me%Prop%IUB + 1
                Me%PropIndex%IP     = Me%Prop%IUB
            end if

            if(Me%PropIndex%POP .eq. null_int)then

                Me%Prop%IUB         = Me%Prop%IUB + 1
                Me%PropIndex%POP    = Me%Prop%IUB
            end if

        end if

        Me%Prop%IUB       = Me%Prop%IUB + 1
        Me%PropIndex%CarbonDioxide      = Me%Prop%IUB

        Me%Prop%IUB       = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen             = Me%Prop%IUB

        Species => Me%FirstSpecies
        do while(associated(Species))

            !Particles (food) index number      
            Particles => Species%FirstParticles
            do while(associated(Particles))

                !water quality properties - organisms
                if (Particles%ID%Name .eq. 'phytoplankton') then
                    if(Me%PropIndex%phyto .eq. null_int)then
                        Me%Prop%IUB                 = Me%Prop%IUB + 1
                        Me%PropIndex%phyto          = Me%Prop%IUB
                    end if
                end if

                if (Particles%ID%Name .eq. 'zooplankton') then
                    if(Me%PropIndex%zoo .eq. null_int)then
                        Me%Prop%IUB                 = Me%Prop%IUB + 1
                        Me%PropIndex%zoo            = Me%Prop%IUB
                    end if
                end if

                if (Particles%ID%Name .eq. 'ciliate') then
                    if(Me%PropIndex%ciliate .eq. null_int)then
                        Me%Prop%IUB                 = Me%Prop%IUB + 1
                        Me%PropIndex%ciliate        = Me%Prop%IUB
                    end if
                end if

                if (Particles%ID%name .eq. 'bacteria') then
                    if(Me%PropIndex%bacteria .eq. null_int)then
                        Me%Prop%IUB                 = Me%Prop%IUB + 1
                        Me%PropIndex%bacteria       = Me%Prop%IUB
                    end if
                end if

                if (Particles%ID%name .eq. 'cohesive sediment') then
                    if(Me%PropIndex%sediments .eq. null_int)then
                        Me%Prop%IUB                 = Me%Prop%IUB + 1
                        Me%PropIndex%sediments      = Me%Prop%IUB
                    end if
                end if  

                !life and wq - diatoms
                if (Particles%ID%name .eq. 'diatoms') then

                    if (Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then

                        if(Me%PropIndex%diatoms .eq. null_int)then
                            Me%Prop%IUB             = Me%Prop%IUB + 1
                            Me%PropIndex%diatoms    = Me%Prop%IUB
                        end if

                        if(Particles%Silica .eq. 1)then

                            if(Me%PropIndex%BioSilica .eq. null_int)then
                                Me%Prop%IUB                 = Me%Prop%IUB + 1
                                Me%PropIndex%BioSilica      = Me%Prop%IUB
                            end if 

                            if(Me%PropIndex%DissSilica .eq. null_int)then
                                Me%Prop%IUB                 = Me%Prop%IUB + 1
                                Me%PropIndex%DissSilica     = Me%Prop%IUB
                            end if

                        end if

                    else !if life

                        if(Me%PropIndex%DiatomsC .eq. null_int)then

                            Me%Prop%IUB                     = Me%Prop%IUB + 1
                            Me%PropIndex%DiatomsC           = Me%Prop%IUB

                            Me%Prop%IUB                     = Me%Prop%IUB + 1
                            Me%PropIndex%DiatomsN           = Me%Prop%IUB

                            Me%Prop%IUB                     = Me%Prop%IUB + 1
                            Me%PropIndex%DiatomsP           = Me%Prop%IUB

                            Me%Prop%IUB                     = Me%Prop%IUB + 1
                            Me%PropIndex%DiatomsChl         = Me%Prop%IUB

                            Me%Prop%IUB                     = Me%Prop%IUB + 1
                            Me%PropIndex%DiatomsSi          = Me%Prop%IUB

                        end if

                    end if !model 

                end if !if is diatoms  

                !life properties - organisms
                if (Particles%ID%name .eq. 'autotrophic flagellates') then

                    if(Me%PropIndex%Mix_FlagellateC .eq. null_int)then
                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Mix_FlagellateC    = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Mix_FlagellateN    = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Mix_FlagellateP    = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Mix_FlagellateChl  = Me%Prop%IUB
                    end if

                end if

                if (Particles%ID%name .eq. 'picoalgae') then

                    if(Me%PropIndex%PicoalgaeC .eq. null_int)then
                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%PicoalgaeC         = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%PicoalgaeN         = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%PicoalgaeP         = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%PicoalgaeChl       = Me%Prop%IUB
                    end if

                end if

                if (Particles%ID%name .eq. 'Flagellate') then

                    if(Me%PropIndex%FlagellateC .eq. null_int)then
                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%FlagellateC        = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%FlagellateN        = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%FlagellateP        = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%FlagellateChl      = Me%Prop%IUB
                    end if

                end if

                if (Particles%ID%name .eq. 'microzooplankton') then

                    if(Me%PropIndex%MicrozooplanktonC .eq. null_int)then
                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%MicrozooplanktonC  = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%MicrozooplanktonN  = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%MicrozooplanktonP  = Me%Prop%IUB
                    end if

                    end if

                if (Particles%ID%name .eq. 'heterotrophic nanoflagellate') then

                    if(Me%PropIndex%Het_NanoflagellateC .eq. null_int)then
                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Het_NanoflagellateC= Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Het_NanoflagellateN= Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Het_NanoflagellateP= Me%Prop%IUB
                    end if

                end if

                if (Particles%ID%name .eq. 'mesozooplankton') then

                    if(Me%PropIndex%MesozooplanktonC .eq. null_int)then
                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%MesozooplanktonC   = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%MesozooplanktonN   = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%MesozooplanktonP   = Me%Prop%IUB
                    end if

                end if

                if (Particles%ID%name .eq. 'heterotrophic bacteria') then

                    if(Me%PropIndex%Het_BacteriaC .eq. null_int)then
                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Het_BacteriaC      = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Het_BacteriaN      = Me%Prop%IUB

                        Me%Prop%IUB                     = Me%Prop%IUB + 1
                        Me%PropIndex%Het_BacteriaP      = Me%Prop%IUB
                    end if

                end if

                Particles => Particles%Next
            end do

            !Predator index number      
            Predator => Species%FirstPredator
            do while(associated(Predator))

                if (Predator%ID%name .eq. 'shrimp') then
                    Me%Prop%IUB          = Me%Prop%IUB + 1
                    Me%PropIndex%Shrimp  = Me%Prop%IUB
                end if

                if (Predator%ID%name .eq. 'crab') then
                    Me%Prop%IUB          = Me%Prop%IUB + 1
                    Me%PropIndex%Crab    = Me%Prop%IUB
                end if

                if (Predator%ID%name .eq. 'oystercatcher') then
                    Me%Prop%IUB          = Me%Prop%IUB + 1
                    Me%PropIndex%OysterCatcher         = Me%Prop%IUB
                end if

                if (Predator%ID%name .eq. 'eider duck') then
                    Me%Prop%IUB          = Me%Prop%IUB + 1
                    Me%PropIndex%EiderDuck             = Me%Prop%IUB
                end if

                if (Predator%ID%name .eq. 'herring gull') then
                    Me%Prop%IUB          = Me%Prop%IUB + 1
                    Me%PropIndex%HerringGull           = Me%Prop%IUB
                end if

                Predator => Predator%Next
            end do

            !Cohorts index number      
            Cohort => Species%FirstCohort
            do while(associated(Cohort))

                Me%Prop%IUB             = Me%Prop%IUB + 1
                Cohort%StateIndex%M_V   = Me%Prop%IUB
                
                Me%Prop%IUB             = Me%Prop%IUB + 1
                Cohort%StateIndex%M_E   = Me%Prop%IUB

                Me%Prop%IUB             = Me%Prop%IUB + 1
                Cohort%StateIndex%M_H   = Me%Prop%IUB

                Me%Prop%IUB             = Me%Prop%IUB + 1
                Cohort%StateIndex%M_R   = Me%Prop%IUB

                Me%Prop%IUB             = Me%Prop%IUB + 1
                Cohort%StateIndex%L     = Me%Prop%IUB

                Me%Prop%IUB             = Me%Prop%IUB + 1
                Cohort%StateIndex%Age   = Me%Prop%IUB

                Me%Prop%IUB             = Me%Prop%IUB + 1
                Cohort%StateIndex%Number= Me%Prop%IUB

                Cohort => Cohort%Next
            end do

        Species => Species%Next
        end do

    end subroutine PropertyIndexNumber

    !-------------------------------------------------------------------------------

    subroutine ConstructPropertyList

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Species)   , pointer       :: Species
        type(T_Cohort)    , pointer       :: Cohort
        type(T_Particles) , pointer       :: Particles
        type(T_Predator)  , pointer       :: Predator

        !Begin-----------------------------------------------------------------

        !Allocate new propertylist
        if (associated(Me%PropertyList)) deallocate(Me%PropertyList)

        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))

        !life  properties
        if (Me%ComputeOptions%PelagicModel .eq. LifeModel) then
            Me%PropertyList(Me%PropIndex%POC)                   = POC_
        end if

        if(Me%ComputeOptions%Nitrogen)then
            Me%PropertyList(Me%PropIndex%AM)                    = Ammonia_
            Me%PropertyList(Me%PropIndex%PON)                   = PON_
        end if

        if(Me%ComputeOptions%Phosphorus)then
            Me%PropertyList(Me%PropIndex%IP )                   = Inorganic_Phosphorus_
            Me%PropertyList(Me%PropIndex%POP)                   = POP_
        end if


        Me%PropertyList(Me%PropIndex%CarbonDioxide)             = CarbonDioxide_

        Me%PropertyList(Me%PropIndex%Oxygen)                    = Oxygen_

        !Include particles inside the Species in the property list 
        Me%nPropertiesFromBivalve = 0      
        
        Species => Me%FirstSpecies
        do while(associated(Species))

            Particles => Species%FirstParticles
            do while(associated(Particles))

                !wq properties
                !if (GetPropertyIDNumber(trim(Particles%ID%Name)) .eq. GetPropertyIDNumber('phytoplankton')) then
                if (Particles%ID%name .eq.'phytoplankton') then
                    Me%PropertyList(Me%PropIndex%phyto)     = Phytoplankton_
                end if

                if (Particles%ID%name .eq.'zooplankton') then
                    Me%PropertyList(Me%PropIndex%zoo)       = Zooplankton_
                end if

                if (Particles%ID%name .eq.'ciliate') then
                    Me%PropertyList(Me%PropIndex%ciliate)   = Ciliate_
                end if

                if (Particles%ID%name .eq.'bacteria') then
                    Me%PropertyList(Me%PropIndex%bacteria)  = Bacteria_
                end if

                if (Particles%ID%name .eq.'cohesive sediment') then
                    Me%PropertyList(Me%PropIndex%sediments) = Cohesive_Sediment_
                end if

                !life and wq properties
                if (Particles%ID%Name .eq.'diatoms') then

                    if (Me%ComputeOptions%PelagicModel .eq. LifeModel) then

                        Me%PropertyList(Me%PropIndex%DiatomsC)      = Diatom_C_
                        Me%PropertyList(Me%PropIndex%DiatomsN)      = Diatom_N_
                        Me%PropertyList(Me%PropIndex%DiatomsP)      = Diatom_P_
                        Me%PropertyList(Me%PropIndex%DiatomsChl)    = Diatom_Chl_
                        Me%PropertyList(Me%PropIndex%DiatomsSi)     = Diatom_Si_
                        Me%PropertyList(Me%PropIndex%BioSilica)     = BioSilica_
                        Me%PropertyList(Me%PropIndex%DissSilica)    = Silicate_

                    else

                        Me%PropertyList(Me%PropIndex%diatoms)       = Diatoms_

                        if(Particles%Silica .eq. 1)then
                            Me%PropertyList(Me%PropIndex%BioSilica) = BioSilica_
                            Me%PropertyList(Me%PropIndex%DissSilica)= Silicate_
                        end if

                    end if
                
                end if

                if (Particles%ID%Name .eq.'autotrophic flagellates') then

                    Me%PropertyList(Me%PropIndex%Mix_FlagellateC)   = Mix_Flagellate_C_
                    Me%PropertyList(Me%PropIndex%Mix_FlagellateN)   = Mix_Flagellate_N_
                    Me%PropertyList(Me%PropIndex%Mix_FlagellateP)   = Mix_Flagellate_P_
                    Me%PropertyList(Me%PropIndex%Mix_FlagellateChl) = Mix_Flagellate_Chl_

                end if

                if (Particles%ID%Name .eq.'picoalgae') then

                    Me%PropertyList(Me%PropIndex%PicoalgaeC)        = Picoalgae_C_
                    Me%PropertyList(Me%PropIndex%PicoalgaeN)        = Picoalgae_N_
                    Me%PropertyList(Me%PropIndex%PicoalgaeP)        = Picoalgae_P_
                    Me%PropertyList(Me%PropIndex%PicoalgaeChl)      = Picoalgae_Chl_

                end if

                if (Particles%ID%Name .eq.'flagellates') then

                    Me%PropertyList(Me%PropIndex%FlagellateC)       = Flagellate_C_
                    Me%PropertyList(Me%PropIndex%FlagellateN)       = Flagellate_N_
                    Me%PropertyList(Me%PropIndex%FlagellateP)       = Flagellate_P_
                    Me%PropertyList(Me%PropIndex%FlagellateChl)     = Flagellate_Chl_

                end if


                if (Particles%ID%Name .eq.'microzooplankton') then

                    Me%PropertyList(Me%PropIndex%MicrozooplanktonC) = Microzooplankton_C_
                    Me%PropertyList(Me%PropIndex%MicrozooplanktonN) = Microzooplankton_N_
                    Me%PropertyList(Me%PropIndex%MicrozooplanktonP) = Microzooplankton_P_

                end if

                if (Particles%ID%Name .eq.'heterotrophic nanoflagellate') then

                    Me%PropertyList(Me%PropIndex%Het_NanoflagellateC) = Het_Nanoflagellate_C_
                    Me%PropertyList(Me%PropIndex%Het_NanoflagellateN) = Het_Nanoflagellate_N_
                    Me%PropertyList(Me%PropIndex%Het_NanoflagellateP) = Het_Nanoflagellate_P_

                end if

                if (Particles%ID%Name .eq.'mesozooplankton') then

                    Me%PropertyList(Me%PropIndex%MesozooplanktonC)    = Mesozooplankton_C_
                    Me%PropertyList(Me%PropIndex%MesozooplanktonN)    = Mesozooplankton_N_
                    Me%PropertyList(Me%PropIndex%MesozooplanktonP)    = Mesozooplankton_P_

                end if


                if (Particles%ID%Name .eq.'heterotrophic bacteria') then

                    Me%PropertyList(Me%PropIndex%Het_BacteriaC)       = Het_Bacteria_C_
                    Me%PropertyList(Me%PropIndex%Het_BacteriaN)       = Het_Bacteria_N_
                    Me%PropertyList(Me%PropIndex%Het_BacteriaP)       = Het_Bacteria_P_

                end if

                Particles => Particles%Next
            end do

            !Include Predator in the property list        
            Predator => Species%FirstPredator
            do while(associated(Predator))

                if (Predator%ID%name .eq.'shrimp') then
                    Me%PropertyList(Me%PropIndex%Shrimp)         = Shrimp_
                end if

                if (Predator%ID%name .eq.'crab') then
                    Me%PropertyList(Me%PropIndex%Crab)           = Crab_
                end if

                if (Predator%ID%name .eq.'oystercatcher') then
                    Me%PropertyList(Me%PropIndex%OysterCatcher)  = OysterCatcher_
                end if

                if (Predator%ID%name .eq.'eider duck') then
                    Me%PropertyList(Me%PropIndex%EiderDuck)      = EiderDuck_
                end if

                if (Predator%ID%name .eq.'herring gull') then
                    Me%PropertyList(Me%PropIndex%HerringGull)    = HerringGull_
                end if

                Predator => Predator%Next
            end do

            Cohort => Species%FirstCohort
            do while(associated(Cohort))

                Me%PropertyList(Cohort%StateIndex%M_V)     = &
                GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" structure")  
                 
                Me%PropertyList(Cohort%StateIndex%M_E)     = &
                GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" reserves")
                
                Me%PropertyList(Cohort%StateIndex%M_H)     = &
                GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" maturity")
                
                Me%PropertyList(Cohort%StateIndex%M_R)     = &
                GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" reproduction")
                
                Me%PropertyList(Cohort%StateIndex%L)       = &
                GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" length")
                
                Me%PropertyList(Cohort%StateIndex%Age)     = &
                GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" age")
                
                Me%PropertyList(Cohort%StateIndex%Number)  = &
                GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" number")
                
                Me%nPropertiesFromBivalve = Me%nPropertiesFromBivalve + 7

                Cohort => Cohort%Next
            end do

        Species => Species%Next
        end do

    end subroutine ConstructPropertyList

    !-------------------------------------------------------------------------------    

    subroutine PrepareRunByIndex (Index)

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Species)         , pointer           :: Species
        type(T_Species)         , pointer           :: SpeciesAgain
        type(T_Cohort)          , pointer           :: Cohort
        type(T_Particles)       , pointer           :: NewParticles
        integer                                     :: Index
        integer                                     :: Number
        integer                                     :: Shrimp, Crab, OysterCatcher
        integer                                     :: EiderDuck, HerringGull
        real                                        :: TotalNumberOfIndividuals
        integer                                     :: iIndexOutput
        
        !Begin-----------------------------------------------------------------
        
        Shrimp        = Me%PropIndex%Shrimp
        Crab          = Me%PropIndex%Crab
        OysterCatcher = Me%PropIndex%OysterCatcher
        EiderDuck     = Me%PropIndex%EiderDuck
        HerringGull   = Me%PropIndex%HerringGull
        
        if (Me%PropIndex%Shrimp .eq. null_int) then 
            Me%ExternalVar%InitialShrimp (Index) = 0.0
        else
            Me%ExternalVar%InitialShrimp (Index) = Me%ExternalVar%Mass(Me%PropIndex%Shrimp,Index)
        end if


        if (Me%PropIndex%phyto .eq. null_int) then 
            Me%ExternalVar%InitialPhyto (Index) = 0.0
        else
            Me%ExternalVar%InitialPhyto (Index) = Me%ExternalVar%Mass(Me%PropIndex%phyto,Index)
        end if
        
        call ConvertUnits (Index)
        
        Species => Me%FirstSpecies
        do while(associated(Species))
        
            Cohort => Species%FirstCohort
            do while(associated(Cohort))
                
                call ComputeCohortLarvaeState(Cohort, Index, Species%LarvaeMaxSize)
                    
                Cohort => Cohort%Next
            end do
            Species => Species%Next
        end do

                
        Species => Me%FirstSpecies
        do while(associated(Species))
        
            Species%Total_CR        = 0.0
            Species%Total_CR_Larvae = 0.0
        
            call RestoreParticlesList (Species)
        
            SpeciesAgain => Me%FirstSpecies
            do while(associated(SpeciesAgain))

                Cohort => SpeciesAgain%FirstCohort
                do while(associated(Cohort))
                
                    if (Species%FeedOnLarvae .and. Cohort%Larvae) then 
                        
                        !Add this cohort to the list of properties
                        allocate(NewParticles)

                        call AddParticles (Species, NewParticles)

                        call ConstructParticlesFromLarvae (SpeciesAgain, Cohort, NewParticles, Index)

                        nullify(NewParticles)
                        
                    end if
                    
                    Cohort => Cohort%Next

                end do
                
                SpeciesAgain => SpeciesAgain%Next

            end do
           
            Species => Species%Next

        end do
        
        TotalNumberOfIndividuals = 0.0

        Species => Me%FirstSpecies
        do while(associated(Species))

            !/m3
            Species%PopulationProcesses%TNStartTimeStep          = 0.0
            Species%PopulationProcesses%TNNonLarvaeStartTimeStep = 0.0
            
            Species%PopulationProcesses%TN          = 0.0
            Species%PopulationProcesses%NCoh        = 0.0
            Species%PopulationProcesses%TBio        = 0.0
            Species%PopulationProcesses%Cr          = 0.0
            Species%PopulationProcesses%Fil         = 0.0
            Species%PopulationProcesses%Ing         = 0.0
            Species%PopulationProcesses%Ass         = 0.0
            Species%PopulationProcesses%CO2         = 0.0
            Species%PopulationProcesses%H2O         = 0.0
            Species%PopulationProcesses%O           = 0.0
            Species%PopulationProcesses%NH3         = 0.0
            Species%PopulationProcesses%PO4         = 0.0
            Species%PopulationProcesses%m_A         = 0.0
            Species%PopulationProcesses%m_O         = 0.0
            Species%PopulationProcesses%m_F         = 0.0
            Species%PopulationProcesses%m_nat       = 0.0
            Species%PopulationProcesses%m_shr       = 0.0
            Species%PopulationProcesses%m_cra       = 0.0
            Species%PopulationProcesses%m_oys       = 0.0
            Species%PopulationProcesses%m_duck      = 0.0
            Species%PopulationProcesses%m_gull      = 0.0
            Species%PopulationProcesses%m_low       = 0.0
            Species%PopulationProcesses%m_self      = 0.0
            Species%PopulationProcesses%m_others    = 0.0
            Species%PopulationProcesses%m_vel       = 0.0
            Species%PopulationProcesses%m_settle    = 0.0
            Species%PopulationProcesses%Massm_A     = 0.0
            Species%PopulationProcesses%Massm_O     = 0.0
            Species%PopulationProcesses%Massm_F     = 0.0
            Species%PopulationProcesses%Massm_nat   = 0.0
            Species%PopulationProcesses%Massm_shr   = 0.0
            Species%PopulationProcesses%Massm_cra   = 0.0
            Species%PopulationProcesses%Massm_oys   = 0.0
            Species%PopulationProcesses%Massm_duck  = 0.0
            Species%PopulationProcesses%Massm_gull  = 0.0
            Species%PopulationProcesses%Massm_low   = 0.0
            Species%PopulationProcesses%Massm_self  = 0.0
            Species%PopulationProcesses%Massm_others= 0.0
            Species%PopulationProcesses%Massm_vel   = 0.0
            Species%PopulationProcesses%Massm_settle= 0.0
            Species%PopulationProcesses%TNField     = 0.0
            Species%PopulationProcesses%MaxLength   = 0.0
            
            Cohort => Species%FirstCohort
            do while(associated(Cohort))
            
                Number = Cohort%StateIndex%Number 
            
                !in the begining of the simulation, the cohort is always alive
                if(Me%ExternalVar%Mass(Number, Index) == 0.0)then
                    Cohort%Dead = 1
                else
                    Cohort%Dead = 0
                endif

                Cohort%Processes%DeathByAge                   = 0.0
                Cohort%Processes%DeathByOxygen                = 0.0
                Cohort%Processes%DeathByStarvation            = 0.0
                Cohort%Processes%DeathByNatural               = 0.0
                Cohort%Processes%PredationByShrimps           = 0.0
                Cohort%Processes%PredationByCrabs             = 0.0
                Cohort%Processes%PredationByOysterCatchers    = 0.0
                Cohort%Processes%PredationByEiderDucks        = 0.0
                Cohort%Processes%PredationByHerringGull       = 0.0
                Cohort%Processes%DeathByLowNumbers            = 0.0
                Cohort%Processes%DeathBySelfPredation         = 0.0
                Cohort%Processes%DeathByLarvaePredationByOthers = 0.0
                Cohort%Processes%DeathByVelocity              = 0.0
                Cohort%Processes%DeathByWrongSettlement       = 0.0
                Cohort%Processes%DeathByDensityLimit          = 0.0
                
                Species%PopulationProcesses%TNStartTimeStep =  Species%PopulationProcesses%TNStartTimeStep        + &
                                                               Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index)
                if(.not. Cohort%Larvae)then                           
                    Species%PopulationProcesses%TNNonLarvaeStartTimeStep =  &
                                                            Species%PopulationProcesses%TNNonLarvaeStartTimeStep + &
                                                            Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index)
                endif
                                                              
                if(associated(Cohort%FeedingOn)) deallocate (Cohort%FeedingOn)
                
                allocate(Cohort%FeedingOn(Species%nParticles, 3)) !store values, Col = Filtered ingested assimilated 
                Cohort%FeedingOn = 0.0  
                      
                Cohort => Cohort%Next

            end do
            
            TotalNumberOfIndividuals = TotalNumberOfIndividuals + Species%PopulationProcesses%TNStartTimeStep
                                    
            Species => Species%Next

        end do
        
        if(TotalNumberOfIndividuals .gt. 0.0)then
            Me%ComputeThisIndex = .true.
        else
            Me%ComputeThisIndex = .false.
        endif
        
        if(Me%OutputON)then
            Me%OutputThisIndex = .false.
            do iIndexOutput = 1, Me%nIndexOutputs
                if(Index == Me%IndexOutputs(iIndexOutput))then
                    Me%OutputThisIndex = .true.
                endif
            enddo
        endif
        
    end subroutine PrepareRunByIndex
    
    !-------------------------------------------------------------------------------    

    subroutine ConvertUnits (Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        type(T_Species)         , pointer           :: Species
        type(T_Cohort)          , pointer           :: Cohort
        integer                                     :: Number
        integer                                     :: Shrimp, Crab, OysterCatcher
        integer                                     :: EiderDuck, HerringGull

        !Begin-----------------------------------------------------------------
        
        Shrimp        = Me%PropIndex%Shrimp
        Crab          = Me%PropIndex%Crab
        OysterCatcher = Me%PropIndex%OysterCatcher
        EiderDuck     = Me%PropIndex%EiderDuck
        HerringGull   = Me%PropIndex%HerringGull
        
        if (Me%DensityUnits .eq. 0) then

            !From m2  to m3
            Me%ConvertionFactor = Me%ExternalVar%CellArea(Index)/Me%ExternalVar%CellVolume(Index)
                    
        else
        
            !the bivalve units are in m3
            Me%ConvertionFactor = 1.0
        
        end if
                      
        Species => Me%FirstSpecies
        do while(associated(Species))

            Cohort => Species%FirstCohort
            do while(associated(Cohort))
            
                Number = Cohort%StateIndex%Number 

                Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index) * Me%ConvertionFactor
                
                              
                Cohort => Cohort%Next

            end do
                                    
            Species => Species%Next

        end do
        
        !Convert also the predators numbers of organisms density to the same units 
        if(Shrimp .ne. null_int)then

            Me%ExternalVar%Mass(Shrimp,Index) = Me%ExternalVar%Mass(Shrimp,Index) * Me%ConvertionFactor
            
        end if

        if(Crab .ne. null_int)then

            Me%ExternalVar%Mass(Crab,Index) = Me%ExternalVar%Mass(Crab,Index) * Me%ConvertionFactor
            
        end if
        
        if(Me%PropIndex%OysterCatcher .ne. null_int)then

            Me%ExternalVar%Mass(OysterCatcher,Index) = Me%ExternalVar%Mass(OysterCatcher,Index) * Me%ConvertionFactor
            
        end if

        if(Me%PropIndex%EiderDuck .ne. null_int)then

            Me%ExternalVar%Mass(EiderDuck,Index) = Me%ExternalVar%Mass(EiderDuck,Index) * Me%ConvertionFactor
            
        end if
        
        if(Me%PropIndex%HerringGull .ne. null_int)then

            Me%ExternalVar%Mass(HerringGull,Index) = Me%ExternalVar%Mass(HerringGull,Index) * Me%ConvertionFactor
            
        end if

    end subroutine ConvertUnits

    !-------------------------------------------------------------------------------    

    subroutine RestoreParticlesList (ObjSpecies)

        !Arguments-----------------------------------------------------------------
        type (T_Species), pointer        :: ObjSpecies

        !Local-----------------------------------------------------------------
        type (T_Particles), pointer      :: Particles           => null()
        type (T_Particles), pointer      :: PreviousParticles   => null()
        !----------------------------------------------------------------------

        nullify(Particles, PreviousParticles)

        Particles  => ObjSpecies%FirstParticles
        do while (associated(Particles))
        
            Particles%Total_CR = 0.0

            if (Particles%Larvae) then  !this property will be removed             

                if(associated(PreviousParticles))then
                    PreviousParticles%Next      => Particles%Next
                else
                    ObjSpecies%FirstParticles   => Particles%Next
                end if

                ObjSpecies%nParticles  = ObjSpecies%nParticles - 1
                
                Particles  => Particles%Next


            else

                PreviousParticles => Particles
                Particles  => Particles%Next
            
            endif
        enddo 

    end subroutine RestoreParticlesList

    !-------------------------------------------------------------------------------

    subroutine ConstructParticlesFromLarvae (Species, Cohort, NewParticles, Index)

        !Arguments------------------------------------------------------------------
        type (T_Species),      pointer          :: Species
        type (T_Cohort),       pointer          :: Cohort
        type (T_Particles),    pointer          :: NewParticles
        integer                                 :: Index

        !External-------------------------------------------------------------------
        
        !Begin----------------------------------------------------------------------

        NewParticles%ID%IDNumber      = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" number")

        NewParticles%ID%Name          = Cohort%ID%Name
        
        NewParticles%ID%Description   = Cohort%ID%Description
        
        NewParticles%Composition%nC   = Species%SpeciesComposition%ReservesComposition%nC

        NewParticles%Composition%nH   = Species%SpeciesComposition%ReservesComposition%nH

        NewParticles%Composition%nO   = Species%SpeciesComposition%ReservesComposition%nO

        NewParticles%Composition%nN   = Species%SpeciesComposition%ReservesComposition%nN

        NewParticles%Composition%nP   = Species%SpeciesComposition%ReservesComposition%nP
        
        NewParticles%Organic          = 1

        NewParticles%Size             = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
        
        NewParticles%Larvae           = .true.
        
        NewParticles%LarvaeSpeciesID  = Species%ID%IDNumber

        NewParticles%LarvaeDensity    = Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index)

        NewParticles%LarvaeBiomass    = Me%ExternalVar%Mass(Cohort%StateIndex%M_E,Index) +  &
                                        Me%ExternalVar%Mass(Cohort%StateIndex%M_V,Index) +  &
                                        Me%ExternalVar%Mass(Cohort%StateIndex%M_R,Index)
                                        
        NewParticles%LarvaeSelfPredated     = 0.0
        
        NewParticles%LarvaePredatedByOthers = 0.0
        
        if(NewParticles%LarvaeBiomass .gt. 0.0)then
            NewParticles%F_E          = Me%ExternalVar%Mass(Cohort%StateIndex%M_E,Index) / NewParticles%LarvaeBiomass
        else
            NewParticles%F_E          = 0.0
        endif
        
        NewParticles%Total_CR         = 0.0
                  
    end subroutine ConstructParticlesFromLarvae    

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTO

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine GetBivalvePropertyList(Bivalve_ID, PropertyList, STAT)

        !Arguments-------------------------------------------------------------
        integer                                   :: Bivalve_ID
        integer, dimension(:), pointer            :: PropertyList
        integer, optional, intent(OUT)            :: STAT

        !External--------------------------------------------------------------
        integer                                   :: ready_

        !Local-----------------------------------------------------------------
        integer                                   :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Bivalve_ID, ready_)    

if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBivalve_, Me%InstanceID)

            PropertyList => Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetBivalvePropertyList

    !--------------------------------------------------------------------------

    subroutine GetDTBivalve(Bivalve_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Bivalve_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Bivalve_ID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetDTBivalve

    !--------------------------------------------------------------------------

    subroutine GetBivalveSize(Bivalve_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer :: Bivalve_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer :: ready_

        !Local-----------------------------------------------------------------
        integer :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Bivalve_ID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBivalveSize

    !--------------------------------------------------------------------------

    subroutine GetBivalvePropIndex (Bivalve_ID, PropertyIDNumber, PropertyIndex, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Bivalve_ID
        integer,           intent(IN )      :: PropertyIDNumber
        integer,           intent(OUT)      :: PropertyIndex
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, CurrentIndex
        logical                             :: found

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Bivalve_ID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            found = .false.
            do CurrentIndex = Me%Prop%ILB,Me%Prop%IUB

                if (PropertyIDNumber.eq. Me%PropertyList(CurrentIndex))then
                    PropertyIndex = CurrentIndex
                    found = .true.
                    exit
                end if

            end do

            if(.not. found)then
                STAT_ = NOT_FOUND_ERR_
            else
                STAT_ = SUCCESS_
            endif

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBivalvePropIndex

    !--------------------------------------------------------------------------

    subroutine GetBivalveListDeadIDS(BivalveID, ListDeadIDS, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: BivalveID
        integer, dimension(:), pointer      :: ListDeadIDS
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer :: STAT_
        integer :: ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BivalveID, ready_)    

if1 :       if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

                call Read_Lock(mBivalve_, Me%InstanceID)

                ListDeadIDS => Me%ListDeadIDs(:)

                STAT_ = SUCCESS_
            else 
                STAT_ = ready_
            end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBivalveListDeadIDS

    !--------------------------------------------------------------------------

    subroutine GetBivalveNewborns(BivalveID, ListNewbornsIDs, MatrixNewborns, STAT)

        !Arguments-------------------------------------------------------------
        integer                               :: BivalveID
        integer, optional, dimension(:), pointer        :: ListNewbornsIDs
        real,    dimension(:,:), pointer      :: MatrixNewborns
        integer, optional, intent(OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                               :: STAT_
        integer                               :: ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BivalveID, ready_)    

if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.  (ready_ .EQ. READ_LOCK_ERR_)) then

            if(present(ListNewbornsIDs)) then
                call Read_Lock(mBivalve_, Me%InstanceID)
                ListNewbornsIDs => Me%ListNewbornsIDs(:)
            endif
            
            call Read_Lock(mBivalve_, Me%InstanceID)
            MatrixNewborns => Me%MatrixNewborns(:,:)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBivalveNewborns
    
    !--------------------------------------------------------------------------

    subroutine GetBivalveNewBornParameters(Bivalve_ID, SpeciesIDNumber, M_V0, M_E0, M_H0, L_0, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Bivalve_ID
        integer,           intent(IN)       :: SpeciesIDNumber
        real,    optional, intent(OUT)      :: M_V0, M_E0, M_H0 , L_0
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer :: STAT_
        integer :: ready_      
        type(T_Species),  pointer           :: Species

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Bivalve_ID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            Species => Me%FirstSpecies
            do while(associated(Species))

                if(Species%ID%IDNumber == SpeciesIDNumber)then

                    if(present(M_V0)) M_V0 = Species%IndividualParameters%MVb
                    if(present(M_E0)) M_E0 = Species%IndividualParameters%MEb
                    if(present(M_H0)) M_H0 = Species%IndividualParameters%MHb
                    if(present(L_0 )) L_0  = Species%IndividualParameters%Lb

                endif

                Species => Species%Next
            end do 

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBivalveNewBornParameters

     !--------------------------------------------------------------------------

    subroutine GetBivalveOtherParameters(Bivalve_ID, SpeciesIDNumber, MinObsLength, LarvaeMaxSize, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Bivalve_ID
        integer,           intent(IN)       :: SpeciesIDNumber
        real,    optional, intent(OUT)      :: MinObsLength, LarvaeMaxSize
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_      
        type(T_Species),  pointer           :: Species

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Bivalve_ID, ready_)    

cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
        
            Species => Me%FirstSpecies
            do while(associated(Species))

                if(Species%ID%IDNumber == SpeciesIDNumber)then

                    if(present(MinObsLength)) MinObsLength   = Species%MinObsLength
                    if(present(LarvaeMaxSize)) LarvaeMaxSize = Species%LarvaeMaxSize

                endif

                Species => Species%Next
            end do 

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBivalveOtherParameters

    !--------------------------------------------------------------------------
   
    subroutine SetSettlementOnBivalve(Bivalve_ID, SpeciesIDNumber, SettlementProbability, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Bivalve_ID
        integer,           intent(IN)       :: SpeciesIDNumber
        real, dimension(:), pointer         :: SettlementProbability
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_      
        type(T_Species),  pointer           :: Species

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Bivalve_ID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_)then
        
            Species => Me%FirstSpecies
            do while(associated(Species))

                if(Species%ID%IDNumber == SpeciesIDNumber)then

                    Species%SettlementProbability = SettlementProbability

                endif

                Species => Species%Next
            end do 

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine SetSettlementOnBivalve

    !--------------------------------------------------------------------------
   
    subroutine UnGetBivalve1D_I(BivalveID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: BivalveID
        integer, dimension(:), pointer      :: Array
        integer, intent(OUT), optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BivalveID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBivalve_, Me%InstanceID, "UnGetBivalve1D_I")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBivalve1D_I

    !--------------------------------------------------------------------------

    subroutine UnGetBivalve2D_R4(BivalveID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                           :: BivalveID
        real, dimension(:,:), pointer     :: Array
        integer, intent(OUT), optional    :: STAT

        !Local-----------------------------------------------------------------
        integer             :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BivalveID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBivalve_, Me%InstanceID, "UnGetBivalve2D_R4")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBivalve2D_R4

    !--------------------------------------------------------------------------

    integer function SearchPropIndex (PropIDNumber)

        !Arguments-------------------------------------------------------------
        integer,  intent(IN )           :: PropIDNumber
        !Local-----------------------------------------------------------------
        integer           :: CurrentIndex

        !----------------------------------------------------------------------

        SearchPropIndex = UNKNOWN_

        do CurrentIndex = Me%Prop%ILB, Me%Prop%IUB

            if (PropIDNumber == Me%PropertyList(CurrentIndex))then
                SearchPropIndex = CurrentIndex
                exit
            end if

        end do

    end function SearchPropIndex


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIE

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyBivalve(ObjBivalveID         , &
                             Temperature          , &
                             Salinity             , & 
                             Mass                 , &
                             OpenPoints           , &
                             WaterVolumeIN        , &
                             CellAreaIN           , &
                             VelocityIN           , &
                             TimeIN               , &
                             STAT)
        !Arguments------------------------------------------------------------------
        integer                                         :: ObjBivalveID
        real,    pointer, dimension(:  )                :: Temperature
        real,    pointer, dimension(:  )                :: Salinity
        real,    pointer, dimension(:,:)                :: Mass
        integer, pointer, dimension(:  )                :: OpenPoints
        real(8), pointer, dimension(:  )                :: WaterVolumeIN
        real,    pointer, dimension(:  )                :: CellAreaIN
        real,    pointer, dimension(:  )                :: VelocityIN
        type(T_Time)                                    :: TimeIN
        integer, intent(OUT), optional                  :: STAT
       
        !Local----------------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: Index
        !integer                                         :: STAT_CALL
        integer                                     :: Chunk
        !---------------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        CHUNK = CHUNK_I(Me%Array%ILB, Me%Array%IUB)

        
        call Ready(ObjBivalveID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%Temperature  => Temperature
            if (.not. associated(Me%ExternalVar%Temperature))       &
            stop 'ModifyBivalve - ModuleBivalve - ERR01'

            Me%ExternalVar%Salinity     => Salinity
            if (.not. associated(Me%ExternalVar%Salinity))          &   
            stop 'ModifyBivalve - ModuleBivalve - ERR02'

            Me%ExternalVar%Mass         => Mass   !mg/l (g/m3)
            if (.not. associated(Me%ExternalVar%Mass))              &
            stop 'ModifyBivalve - ModuleBivalve - ERR03'
                        
            Me%ExternalVar%CellVolume  => WaterVolumeIN
            if (.not. associated(Me%ExternalVar%CellVolume))        &
            stop 'ModifyBivalve - ModuleBivalve - ERR04'

            Me%ExternalVar%CellArea  => CellAreaIN
            if (.not. associated(Me%ExternalVar%CellArea))          &
            stop 'ModifyBivalve - ModuleBivalve - ERR05'
            
            Me%ExternalVar%Velocity  => VelocityIN
            if (.not. associated(Me%ExternalVar%Velocity))          &
            stop 'ModifyBivalve - ModuleBivalve - ERR06'
            
            Me%ExternalVar%CurrentTime  = TimeIN

            call AllocateAndInitializeByTimeStep
            

            do Index = Me%Array%ILB, Me%Array%IUB

                call ComputeBivalve (Index, OpenPoints(Index))

            enddo
            
            if(Me%OutputON) then
                if (Me%ExternalVar%CurrentTime .ge. Me%BivalveOutputTimes(Me%NextOutPut)) then
                    Me%NextOutPut = Me%NextOutPut + 1
                endif
            end if
            
            call UpdateListDeadAndNewBornIDs
                        
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyBivalve

    !-------------------------------------------------------------------------------    

    subroutine AllocateAndInitializeByTimeStep

        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Species)         , pointer           :: Species
        type(T_Cohort)          , pointer           :: Cohort
        integer                                     :: TotalNumberIndex, Index

        !Begin-----------------------------------------------------------------

        Species  => Me%FirstSpecies
d1:     do while (associated(Species))

            Species%NewbornCohort = .false. 
            
            Cohort  => Species%FirstCohort
d2:         do while (associated(Cohort))

                Cohort%GlobalDeath = 1
                
                Cohort%PreviousAtLeastOneLarvae = Cohort%AtLeastOneLarvae
                
                Cohort%AtLeastOneLarvae = .false.
                
!                if ((Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) .le. Me%ConvertionFactor * 1e-8))then
!                    
!                    Me%ExternalVar%Mass(Cohort%StateIndex%L,    Index) = 0.0
!                    Me%ExternalVar%Mass(Cohort%StateIndex%M_E,  Index) = 0.0
!                    Me%ExternalVar%Mass(Cohort%StateIndex%M_V,  Index) = 0.0
!                    Me%ExternalVar%Mass(Cohort%StateIndex%M_H,  Index) = 0.0
!                    Me%ExternalVar%Mass(Cohort%StateIndex%M_R,  Index) = 0.0
!                    
!                endif
                
            
                do Index = Me%Array%ILB, Me%Array%IUB

                    if ((Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) .gt. 0.0) .and. &
                        (Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      .gt. 0.0) .and. &
                        (Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      .le. Species%LarvaeMaxSize)) then
                                    
                        Cohort%AtLeastOneLarvae = .true.
                                
                    end if

                enddo
                
                if(Cohort%PreviousAtLeastOneLarvae .and. .not. Cohort%AtLeastOneLarvae)then 
                    write(*,*)"BIVALVE : ", trim(Cohort%ID%Name), " is no longer larvae"
                endif

                Cohort  => Cohort%Next  
            end do  d2        

            Species  => Species%Next  
        enddo  d1

        TotalNumberIndex = Me%Array%IUB
                
        !List of Cohorts/Properties IDs from cohorts that died in this time step (died in all cells)
        if(associated(Me%ListDeadIDs)) deallocate (Me%ListDeadIDs)
        allocate(Me%ListDeadIDs(Me%nPropertiesFromBivalve))
        Me%ListDeadIDs = 0
        
        Me%nLastDeadID = 0
        
        !List of Species IDs that will have Newborns in the next time step
        if(associated(Me%ListNewbornsIDs)) deallocate (Me%ListNewbornsIDs)
        allocate(Me%ListNewbornsIDs(Me%nSpecies))
        Me%ListNewbornsIDs = 0
        
        Me%nLastNewbornsID = 0

        !newborns for each species in each index
        if(associated(Me%MatrixNewborns )) deallocate (Me%MatrixNewborns)
        allocate(Me%MatrixNewborns(Me%nSpecies, TotalNumberIndex+1))
        
        Me%MatrixNewborns = 0
               
    end subroutine AllocateAndInitializeByTimeStep

    !--------------------------------------------------------------------------

    subroutine ComputeBivalve(Index, CheckIfOpenPoint)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)           :: Index
        integer, intent(IN)           :: CheckIfOpenPoint
        integer                       :: i, iIndexOutput

        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------

        call PrepareRunByIndex (Index)
        
        if(Me%ComputeThisIndex)then
        
            call ComputeMortalityByVelocity (Index)
            
            call ComputeMortalityByWrongSettlement (Index)

            call ComputeIndividualProcesses (Index, CheckIfOpenPoint)

            call ComputeExtraStarvationMortality (Index)
                    
            call ComputeNaturalMortality (Index)

            call ComputePredation (Index, CheckIfOpenPoint)
            
            call ComputePopulationVariables (Index)
            
            call UpdateCohortState

        end if
          
        
        if(Me%OutputThisIndex)then

            if (Me%ExternalVar%CurrentTime >= Me%BivalveOutputTimes(Me%NextOutPut)) then
            
                do i = 1, Me%nIndexOutputs
                    if(Index == Me%IndexOutputs(i))then
                        iIndexOutput = i
                    endif
                enddo
            
        
                call WriteOutput (Index, iIndexOutput)
            
            endif
            
        endif
        
        call RestoreUnits (Index)
        
   end subroutine ComputeBivalve
   
   !--------------------------------------------------------------------------

   subroutine ComputeMortalityByVelocity (Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)           :: Index

        !Local-----------------------------------------------------------------
        type(T_Species)    , pointer :: Species
        type(T_Cohort)     , pointer :: Cohort
        real                         :: Velocity
        integer                      :: Number, M_V, M_E, M_R
        integer                      :: POC, PON, POP

        !Begin-----------------------------------------------------------------

        Velocity = Me%ExternalVar%Velocity(index)
        
        Species => Me%FirstSpecies
        do while(associated(Species))

        if (Velocity .ge. Species%IndividualParameters%MAX_velocity) then 
            
            Cohort => Species%FirstCohort
            do while(associated(Cohort))
            
                if(Cohort%Dead .eq. 0)then 
            
                    if (.not. Cohort%Larvae) then

                        M_V     = Cohort%StateIndex%M_V
                        M_E     = Cohort%StateIndex%M_E
                        M_R     = Cohort%StateIndex%M_R
                        Number  = Cohort%StateIndex%Number
                        
                        !Mortality by velocity, #/d.m3
                        Cohort%Processes%DeathByVelocity = (Me%ExternalVar%Mass(Number,Index)             * &
                                                          Species%IndividualParameters%m_velocity) / Me%DTDay
                                                          

                        !update the number of organisms in mass matrix
                        Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index)           - &
                                                            (Cohort%Processes%DeathByVelocity * Me%DTDay)
                                                            
                        if (Me%ExternalVar%Mass(Number,Index) .eq. 0.0) then 
                        
                            call ImposeCohortDeath (Index, Species, Cohort) !sets all proc to zero, convert mass to OM, Deadlist
                        
                        else
                        
                            POC     = Me%PropIndex%POC
                            PON     = Me%PropIndex%PON
                            POP     = Me%PropIndex%POP

                            !update POM again
                            Me%ExternalVar%Mass(PON,Index) = Me%ExternalVar%Mass(PON,Index)                                   + &
                                                            ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                            Species%SpeciesComposition%StructureComposition%nN                + &
                                                            (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                            Species%SpeciesComposition%ReservesComposition%nN )               * &
                                                            Species%AuxiliarParameters%N_AtomicMass                           * &
                                                            Cohort%Processes%DeathByVelocity * Me%DTDay           

                            if(Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then


                                if (Me%ComputeOptions%Phosphorus) then

                                    Me%ExternalVar%Mass(POP,Index) = Me%ExternalVar%Mass(POP,Index)                           + &
                                                            ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                            Species%SpeciesComposition%StructureComposition%nP                + &
                                                            (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                            Species%SpeciesComposition%ReservesComposition%nP )               * &
                                                            Species%AuxiliarParameters%P_AtomicMass                           * &
                                                            Cohort%Processes%DeathByVelocity * Me%DTDay           

                                end if

                            else !(if life)


                                Me%ExternalVar%Mass(POC,Index) = Me%ExternalVar%Mass(POC,Index)              + &
                                                                ( Me%ExternalVar%Mass(M_V,Index)             + &
                                                                Me%ExternalVar%Mass(M_E,Index)               + &
                                                                Me%ExternalVar%Mass(M_R,Index) )             * &
                                                                Cohort%Processes%DeathByVelocity * Me%DTDay

                            end if !pelagic model
                            
                        end if ! not (Me%ExternalVar%Mass(Number,Index) .eq. 0.0 
                    
                    else !if larvae 
                        !Mortality by velocity, #/d.m3
                        Cohort%Processes%DeathByVelocity = 0.0
                    end if             
                endif !Cohort%Dead = 0                      

            Cohort => Cohort%Next
            end do
            
        end if 

        Species => Species%Next
        end do

    end subroutine ComputeMortalityByVelocity
   !--------------------------------------------------------------------------

   subroutine ComputeMortalityByWrongSettlement (Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)           :: Index

        !Local-----------------------------------------------------------------
        type(T_Species)    , pointer :: Species
        type(T_Cohort)     , pointer :: Cohort
        integer                      :: Number, M_V, M_E, M_R
        integer                      :: POC, PON, POP
        real                         :: MaxNumberToSettle, NumberToDieFromDensityLimit
        

        !Begin-----------------------------------------------------------------
        
        Species => Me%FirstSpecies
        do while(associated(Species))
            
            Cohort => Species%FirstCohort
            do while(associated(Cohort))
            
                if(Cohort%Dead .eq. 0)then 
                
                    !if the whole has just stopped being a larvae
                    if(Cohort%PreviousAtLeastOneLarvae .and. .not. Cohort%AtLeastOneLarvae)then 
                    
                        M_V     = Cohort%StateIndex%M_V
                        M_E     = Cohort%StateIndex%M_E
                        M_R     = Cohort%StateIndex%M_R
                        Number  = Cohort%StateIndex%Number
                        
                        !exclude current cohort number from the total number non larvae of individuals
                        Species%PopulationProcesses%TNNonLarvaeStartTimeStep = &
                                            Species%PopulationProcesses%TNNonLarvaeStartTimeStep - &
                                            Me%ExternalVar%Mass(Number,Index)

                        !Mortality by settlement, #/d.m3
                        Cohort%Processes%DeathByWrongSettlement = Me%ExternalVar%Mass(Number,Index)          * &
                                                                  (1 - Species%SettlementProbability(Index)) / &
                                                                  Me%DTDay
                                                                  
                        !update the number of organisms in mass matrix
                        Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index)           - &
                                                            (Cohort%Processes%DeathByWrongSettlement * Me%DTDay)
                                   

                        !Mortality by DensityLimitations, #/d.m3
                        if(Species%IndividualParameters%DensityLimOption) then
                            
                            !#/m3
                            MaxNumberToSettle = Species%IndividualParameters%MAX_density * Me%ConvertionFactor - &
                                                Species%PopulationProcesses%TNNonLarvaeStartTimeStep
                        
                            !#/m3
                            NumberToDieFromDensityLimit = Me%ExternalVar%Mass(Number,Index) - MaxNumberToSettle

                            if(NumberToDieFromDensityLimit .lt. 0)then
                                NumberToDieFromDensityLimit = 0.0
                            endif
                        
                            Cohort%Processes%DeathByDensityLimit = NumberToDieFromDensityLimit / Me%DTDay
                                                                      
                            !update the number of organisms in mass matrix
                            Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index)           - &
                                                                (Cohort%Processes%DeathByDensityLimit * Me%DTDay)
                                                                
                            !sum wrongsettlement with densitylimit just to simplify                                   
                            Cohort%Processes%DeathByWrongSettlement = Cohort%Processes%DeathByWrongSettlement + &
                                                                      Cohort%Processes%DeathByDensityLimit            

                        end if
                                                          

                                             
                        if (Me%ExternalVar%Mass(Number,Index) .eq. 0.0) then 
                        !they all died
                        
                            call ImposeCohortDeath (Index, Species, Cohort) !sets all proc to zero, convert mass to OM, Deadlist
                        
                        else
                        
                            POC     = Me%PropIndex%POC
                            PON     = Me%PropIndex%PON
                            POP     = Me%PropIndex%POP

                            !update POM again
                            Me%ExternalVar%Mass(PON,Index) = Me%ExternalVar%Mass(PON,Index)                                   + &
                                                            ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                            Species%SpeciesComposition%StructureComposition%nN                + &
                                                            (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                            Species%SpeciesComposition%ReservesComposition%nN )               * &
                                                            Species%AuxiliarParameters%N_AtomicMass                           * &
                                                            Cohort%Processes%DeathByWrongSettlement * Me%DTDay           

                            if(Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then


                                if (Me%ComputeOptions%Phosphorus) then

                                    Me%ExternalVar%Mass(POP,Index) = Me%ExternalVar%Mass(POP,Index)                           + &
                                                            ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                            Species%SpeciesComposition%StructureComposition%nP                + &
                                                            (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                            Species%SpeciesComposition%ReservesComposition%nP )               * &
                                                            Species%AuxiliarParameters%P_AtomicMass                           * &
                                                            Cohort%Processes%DeathByWrongSettlement * Me%DTDay           

                                end if

                            else !(if life)


                                Me%ExternalVar%Mass(POC,Index) = Me%ExternalVar%Mass(POC,Index)              + &
                                                                ( Me%ExternalVar%Mass(M_V,Index)             + &
                                                                Me%ExternalVar%Mass(M_E,Index)               + &
                                                                Me%ExternalVar%Mass(M_R,Index) )             * &
                                                                Cohort%Processes%DeathByWrongSettlement * Me%DTDay

                            end if !pelagic model
                            
                        end if ! not (Me%ExternalVar%Mass(Number,Index) .eq. 0.0
                    
                    else ! if larvae 
                    
                     Cohort%Processes%DeathByWrongSettlement  = 0.0
                    
                    end if  
                                    
                endif !Cohort%Dead .eq. 0

                Cohort => Cohort%Next
            end do
            
            Species => Species%Next
        end do

    end subroutine ComputeMortalityByWrongSettlement

    !--------------------------------------------------------------------------

    subroutine ComputeIndividualProcesses (Index, CheckIfOpenPoint)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)           :: Index
        integer, intent(IN)           :: CheckIfOpenPoint

        !Local-----------------------------------------------------------------
        type(T_Species)    , pointer :: Species
        type(T_Cohort)     , pointer :: Cohort
        integer                      :: O

        !Begin-----------------------------------------------------------------

        O  = Me%PropIndex%Oxygen

        !if (Me%ExternalVar%Mass(O, Index) .gt. 0.0) then !Bivalve can survive         

            Species => Me%FirstSpecies
d1:         do while(associated(Species))


                call ComputeChemicalIndices    (Index, Species)

                call ComputeAuxiliarParameters (Index, Species)

                call ComputeBivalveCondition   (Index, Species)

                Species => Species%Next

            end do d1

            call ComputeFeedingProcesses (Index, CheckIfOpenPoint) !for all species

            Species => Me%FirstSpecies
    d2:     do while(associated(Species))

                Species%PopulationProcesses%nNewborns = 0

                Cohort => Species%FirstCohort
    d3:         do while(associated(Cohort))
    
                    if (Cohort%Dead .eq. 0) then

                        call ComputeSomaticMaintenance (Species, Cohort)

                        call ComputeMobilization       (Species, Cohort)

                        call ComputeReservesDynamics   (Index, Cohort)

                        call ComputeStructureDynamics  (Index, Species, Cohort)

                        call ComputeMaturity           (Index, Species, Cohort)

                        call ComputeLengthAgeDynamics  (Index, Species, Cohort)

                        call ComputeSpawning (Index, Species, Cohort)

                        call ComputeReproductionDynamics (Index, Cohort)

                        call ComputeInorganicFluxes (Index, Species, Cohort)

                    end  if

                    Cohort => Cohort%Next
                end do d3

                if ((Species%Population) .and. (Species%PopulationProcesses%nNewborns .gt. 0.0 )) then 
                !if it is not population then new borns will be ignored
                
                    Species%NewbornCohort = .true.
                    
                    !Add the SpeciesID to the list of Species ID with new borns in the next time step/m2   
                    Me%MatrixNewborns(Species%ID%ID,1) = Species%ID%IDNumber
                    
                    Me%MatrixNewborns(Species%ID%ID,Index+1) = Species%PopulationProcesses%nNewborns
                    
                    Me%MatrixNewborns(Species%ID%ID,Index+1) = Me%MatrixNewborns(Species%ID%ID,Index+1) * &
                                                                   1./Me%ConvertionFactor                                  
                end if 

            Species => Species%Next
            end do d2

!        else !no oxygen in the system
!
!            Species => Me%FirstSpecies
!d4:         do while(associated(Species))
!
!                Cohort => Species%FirstCohort
!d5:             do while(associated(Cohort))
!
!                    if (Cohort%Dead .eq. 0 ) then
!
!                        Cohort%Processes%DeathByOxygen = Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) / &
!                                                         Me%DTDay ! They will all die
!                                                
!                        !update the number of organisms in mass matrix
!                        Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index)           - &
!                                                            (Cohort%Processes%DeathByOxygen * Me%DTDay)
!
!
!                        call ImposeCohortDeath (Index, Species, Cohort) 
!
!!                        if (Species%nCohorts .eq. Me%nListDeadIDs) then
!                        if ((Species%nCohorts .eq. 1) .and. (Cohort%Dead .eq. 1 )) then                        
!                        !if the last cohort in the population                          
!
!                            Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
!                            
!                        end if
!                        
!                    end if
!
!                    Cohort => Cohort%Next
!                end do d5
!
!            Species => Species%Next
!            end do d4
!
!        end if

    end subroutine ComputeIndividualProcesses

    !--------------------------------------------------------------------------

    subroutine ComputeChemicalIndices (Index, Species)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: Index
        type(T_Species)      , pointer      :: Species

        !Local-----------------------------------------------------------------
        type(T_Particles)    , pointer      :: Particles
        integer                             :: PropertyIndexC, PropertyIndexN
        integer                             :: PropertyIndexP, PropertyIndexChl
        integer                             :: PropertyIndexSi
        real                                :: C_AtomicMass, H_AtomicMass,O_AtomicMass 
        real                                :: P_AtomicMass, N_AtomicMass

        !Begin-----------------------------------------------------------------


        C_AtomicMass  = Species%AuxiliarParameters%C_AtomicMass        
        H_AtomicMass  = Species%AuxiliarParameters%H_AtomicMass        
        O_AtomicMass  = Species%AuxiliarParameters%O_AtomicMass        
        P_AtomicMass  = Species%AuxiliarParameters%P_AtomicMass        
        N_AtomicMass  = Species%AuxiliarParameters%N_AtomicMass        

        Particles => Species%FirstParticles
        d1:     do while(associated(Particles))
        
            if (.not. Particles%Larvae) then 

                if (Particles%Organic .eq. 1) then

                    if (Particles%RatioVariable .eq. 1) then
                    !compute actual ratios based on mass values

                            if (Particles%ID%Name .eq. 'particulate organic matter') then
                                PropertyIndexC   = SearchPropIndex(GetPropertyIDNumber("particulate organic carbon"))
                                PropertyIndexN   = SearchPropIndex(GetPropertyIDNumber("particulate organic nitrogen"))
                                PropertyIndexP   = SearchPropIndex(GetPropertyIDNumber("particulate organic phosphorus"))
                            else
                                PropertyIndexC   = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" carbon")))
                                PropertyIndexN   = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" nitrogen")))
                                PropertyIndexP   = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" phosphorus")))
                            end if

                        !g/g
                        Particles%Ratios%NC_Ratio   = Me%ExternalVar%Mass(PropertyIndexN,Index) /             &
                        Me%ExternalVar%Mass(PropertyIndexC,Index)

                        Particles%Ratios%PC_Ratio   = Me%ExternalVar%Mass(PropertyIndexP,Index) /             &
                        Me%ExternalVar%Mass(PropertyIndexC,Index)

                        if(Particles%Silica .eq. 1)then          

                            PropertyIndexSi  = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" silica")))
                            Particles%Ratios%SiC_Ratio   = Me%ExternalVar%Mass(PropertyIndexSi,Index) /       &
                            Me%ExternalVar%Mass(PropertyIndexC,Index)

                        end if

                        if ((Particles%ID%Name .eq. 'diatoms') .or. (Particles%ID%Name .eq. 'autotrophic flagellates') .or. &
                        (Particles%ID%Name .eq. 'picoalgae') .or. (Particles%ID%Name .eq. 'flagellates')) then 

                            PropertyIndexChl = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" chlorophyll")))

                            Particles%Ratios%ChlC_Ratio = Me%ExternalVar%Mass(PropertyIndexChl,Index) /       &
                            Me%ExternalVar%Mass(PropertyIndexC,Index)

                        end if

                    end if 

                    Particles%Composition%nC = 1

                    Particles%Composition%nH = (Particles%Ratios%HC_Ratio / H_AtomicMass) * C_AtomicMass

                    Particles%Composition%nO = (Particles%Ratios%OC_Ratio / O_AtomicMass) * C_AtomicMass

                    Particles%Composition%nN = (Particles%Ratios%NC_Ratio / N_AtomicMass) * C_AtomicMass

                    Particles%Composition%nP = (Particles%Ratios%PC_Ratio / P_AtomicMass) * C_AtomicMass

                else

                    Particles%Composition%nC = 0.0

                    Particles%Composition%nH = 0.0

                    Particles%Composition%nO = 0.0

                    Particles%Composition%nN = 0.0

                    Particles%Composition%nP = 0.0

                end if
            
            end if !(if not larvae)

            Particles => Particles%Next

        end do d1

    end subroutine ComputeChemicalIndices

    !--------------------------------------------------------------------------

    subroutine ComputeAuxiliarParameters (Index, Species)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: Index
        type(T_Species),    pointer     :: Species

        !Local-----------------------------------------------------------------
        real                            :: Tref, TA, TL, TH, TAL, TAH, T
        real                            :: pM, EG, mu_E, d_V
        real                            :: PAM_FIX, delta_M, kappa

        !Begin-----------------------------------------------------------------

        Tref        = Species%IndividualParameters%Tref
        TA          = Species%IndividualParameters%TA
        TL          = Species%IndividualParameters%TL
        TH          = Species%IndividualParameters%TH
        TAL         = Species%IndividualParameters%TAL
        TAH         = Species%IndividualParameters%TAH
        pM          = Species%IndividualParameters%pM   
        EG          = Species%IndividualParameters%EG     
        mu_E        = Species%IndividualParameters%mu_E
        d_V         = Species%IndividualParameters%d_V
        PAM_FIX     = Species%IndividualParameters%PAM_FIX
        delta_M     = Species%IndividualParameters%delta_M
        kappa       = Species%IndividualParameters%kappa

        !Temperature Correction factor

        !K, Actual Temperature, oC to K
        T = Me%ExternalVar%Temperature(index) + 273.


        if (Species%IndividualParameters%SIMPLE_TEMP .eq. 0.0) then
            Species%AuxiliarParameters%TempCorrection  = exp(TA/Tref-TA/T) *     &
            (1.0 + exp(TAL/Tref - TAL/TL) + exp(TAH/TH-TAH/Tref))/   &
            (1.0 + exp(TAL/T - TAL/TL) + exp(TAH/TH-TAH/T))
        else            
            Species%AuxiliarParameters%TempCorrection  = exp(TA/Tref-TA/T) 
        end if


        !WE, gDW/molC , ash free dry weight to carbon convertion factor for bivalve reserve  
        Species%AuxiliarParameters%WE     = Species%SpeciesComposition%ReservesComposition%nC *  &
                                            Species%AuxiliarParameters%C_AtomicMass            + &
                                            Species%SpeciesComposition%ReservesComposition%nH *  &
                                            Species%AuxiliarParameters%H_AtomicMass            + &
                                            Species%SpeciesComposition%ReservesComposition%nO *  &
                                            Species%AuxiliarParameters%O_AtomicMass            + &
                                            Species%SpeciesComposition%ReservesComposition%nN *  &
                                            Species%AuxiliarParameters%N_AtomicMass            + &
                                            Species%SpeciesComposition%ReservesComposition%nP *  &
                                            Species%AuxiliarParameters%P_AtomicMass   

        !WV, gDW/molC , ash free dry weight to carbon convertion factor for bivalve structure  
        Species%AuxiliarParameters%WV   = Species%AuxiliarParameters%WE

        !Mv, molC(struc)/cm3, volume specific structural mass
        Species%AuxiliarParameters%Mv   = d_V / Species%AuxiliarParameters%WV

        !MHb, molC, Maturity threshold for birth
        Species%AuxiliarParameters%MHb  = Species%IndividualParameters%EHb / mu_E 

        !MHp, molC, Maturity threshold for puberty
        Species%AuxiliarParameters%MHp  = Species%IndividualParameters%EHp / mu_E 

        !y_VE, molC(struc)/molC(reser), yield of structure on reserves
        Species%AuxiliarParameters%y_VE = Species%AuxiliarParameters%Mv * mu_E / EG 

        !kM, d-1, somatic maintenance rate coefficient
        Species%AuxiliarParameters%kM   = pM / EG         

        !kJ, d-1, maturity maintenance rate coefficient
        Species%AuxiliarParameters%kJ   = Species%AuxiliarParameters%kM

        !Lm, cm, maximum length of the species
        if (pM .gt. 0.0) then
            Species%AuxiliarParameters%Lm   = ((kappa * PAM_FIX) / pM) / delta_M 
        else
            Species%AuxiliarParameters%Lm   = 0.0
        end if
        
    end subroutine ComputeAuxiliarParameters

    !--------------------------------------------------------------------------

    subroutine ComputeBivalveCondition (Index, Species)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)           :: Index
        type(T_Species),    pointer   :: Species

        !Local-----------------------------------------------------------------
        type(T_Cohort),     pointer   :: Cohort
        integer         :: L, M_V, M_E, M_R

        !Begin-----------------------------------------------------------------

        Cohort => Species%FirstCohort
        d1:     do while(associated(Cohort))

            L   = Cohort%StateIndex%L
            M_V = Cohort%StateIndex%M_V
            M_E = Cohort%StateIndex%M_E
            M_R = Cohort%StateIndex%M_R

            !Vol, cm3, Structural volume of the organism
            Cohort%BivalveCondition%Vol  = Me%ExternalVar%Mass(M_V, Index) / Species%AuxiliarParameters%Mv

            !E, molC(reserves)/cm3(deb), reserve density
            if (Cohort%BivalveCondition%Vol .eq. 0.0 ) then
            
                Cohort%BivalveCondition%E       = 0.0
                
                Cohort%BivalveCondition%ScaledE = 0.0
                
                Cohort%BivalveCondition%GSR     = 0.0
                
                
            else
                Cohort%BivalveCondition%E = Me%ExternalVar%Mass(M_E, Index) / Cohort%BivalveCondition%Vol
                
                Cohort%BivalveCondition%ScaledE = Cohort%BivalveCondition%E          / &
                                                  (Species%IndividualParameters%Em   / &
                                                   Species%IndividualParameters%mu_E) 

                !GSR, molC(gam)/molC(total), fraction of gametes in the organism
                Cohort%BivalveCondition%GSR = Me%ExternalVar%Mass(M_R, Index)                / &
                                            ( Me%ExternalVar%Mass(M_V, Index)                + &
                                              Me%ExternalVar%Mass(M_E, Index)                + &
                                              Me%ExternalVar%Mass(M_R, Index) )
            end if


            !DW, mgDW, organism total dry weight, 1e3 * (molC * gdw/molC)
            Cohort%BivalveCondition%DW = 1e3                                                                    * &
                                        ( Me%ExternalVar%Mass(M_V, Index) * Species%AuxiliarParameters%WV       + &
                                        (Me%ExternalVar%Mass(M_E, Index) + Me%ExternalVar%Mass(M_R, Index))     * &
                                        Species%AuxiliarParameters%WE )    

            !TotalmolC, molC, organism total molC
            Cohort%BivalveCondition%TotalmolC = Me%ExternalVar%Mass(M_V, Index)                                 + &
                                            Me%ExternalVar%Mass(M_E, Index)                                     + &
                                            Me%ExternalVar%Mass(M_R, Index)

            !TotalmolN, molN, organism total molN
            Cohort%BivalveCondition%TotalmolN = Me%ExternalVar%Mass(M_V, Index)                                 * &
                                            Species%SpeciesComposition%StructureComposition%nN                  + &
                                            (Me%ExternalVar%Mass(M_E, Index) + Me%ExternalVar%Mass(M_R, Index)) * &
                                            Species%SpeciesComposition%ReservesComposition%nN     


            !TotalmolP, molP, organism total molP
            Cohort%BivalveCondition%TotalmolP = Me%ExternalVar%Mass(M_V, Index)                                 * &
                                            Species%SpeciesComposition%StructureComposition%nP                  + &
                                            (Me%ExternalVar%Mass(M_E, Index) + Me%ExternalVar%Mass(M_R, Index)) * &
                                            Species%SpeciesComposition%ReservesComposition%nP     


            Cohort => Cohort%Next
        end do d1

    end subroutine ComputeBivalveCondition

    !--------------------------------------------------------------------------

    subroutine ComputeFeedingProcesses (Index, CheckIfOpenPoint)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)           :: Index
        integer, intent(IN)           :: CheckIfOpenPoint

        !Local-----------------------------------------------------------------
        type(T_Species),    pointer   :: Species
        type(T_Cohort),     pointer   :: Cohort

        !Begin-----------------------------------------------------------------

        !Choose feeding processes model
        if (CheckIfOpenPoint == OpenPoint) then

            if (Me%ComputeOptions%SimpleFiltration) then ! all are simple filtration model

                call ComputeSimpleFiltration (Index)

            else ! complex filtration model

                call ComputeComplexFiltration (Index) ! all will be computed as complex

            end if

        else

            Species => Me%Firstspecies
            d1:     do while(associated(Species))

                Cohort => Species%FirstCohort
                d2:     do while(associated(Cohort))

                    call ImposeNoFiltrationProcess (Cohort)

                    Cohort => Cohort%Next
                end do d2

                Species => Species%Next
            end do d1

        end if

    end subroutine ComputeFeedingProcesses

    !--------------------------------------------------------------------------

    subroutine ComputeSimpleFiltration (Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: Index

        !Local-----------------------------------------------------------------
        type(T_Species),    pointer         :: Species
        type(T_Cohort),     pointer         :: Cohort
        integer                             :: M_H 
        real                                :: F_FIX, PAM_FIX, mu_E, YEX
        real                                :: C_AtomicMass, H_AtomicMass, O_AtomicMass
        real                                :: P_AtomicMass, N_AtomicMass
        real                                :: RATIOHC, RATIOOC, RATIONC, RATIOPC 
        real                                :: Vol, TempCorrection,MHb

        !Begin-----------------------------------------------------------------

        Species   => Me%FirstSpecies
        
d1:     do while(associated(Species))

            F_FIX           = Species%IndividualParameters%F_FIX 
            PAM_FIX         = Species%IndividualParameters%PAM_FIX 
            mu_E            = Species%IndividualParameters%mu_E
            YEX             = Species%IndividualParameters%YEX

            C_AtomicMass    = Species%AuxiliarParameters%C_AtomicMass        
            H_AtomicMass    = Species%AuxiliarParameters%H_AtomicMass        
            O_AtomicMass    = Species%AuxiliarParameters%O_AtomicMass        
            P_AtomicMass    = Species%AuxiliarParameters%P_AtomicMass        
            N_AtomicMass    = Species%AuxiliarParameters%N_AtomicMass        
            TempCorrection  = Species%AuxiliarParameters%TempCorrection         
            MHb             = Species%AuxiliarParameters%MHb  

            !Assumed the samecomposition as the bivalve 
            RATIOHC         = 0.15         
            RATIOOC         = 0.71         
            RATIONC         = 0.3 !0.3395653308501444         
            RATIOPC         = 0.07  

            Cohort => Species%FirstCohort
d2:         do while(associated(Cohort))

                M_H         = Cohort%StateIndex%M_H    
                Vol         = Cohort%BivalveCondition%Vol

                if (Me%ExternalVar%Mass(M_H,Index) .gt. MHb) then !feeding            

                    !Filtration, molC/d
                    Cohort%Processes%FilteredInorganic = 0.0
                    Cohort%Processes%FilteredFood%C    = F_FIX * PAM_FIX / mu_E * TempCorrection * Vol**(2.0/3.0) / YEX
                    Cohort%Processes%FilteredFood%H    = Cohort%Processes%FilteredFood%C * RATIOHC / H_AtomicMass * C_AtomicMass
                    Cohort%Processes%FilteredFood%O    = Cohort%Processes%FilteredFood%C * RATIOOC / O_AtomicMass * C_AtomicMass
                    Cohort%Processes%FilteredFood%N    = Cohort%Processes%FilteredFood%C * RATIONC / N_AtomicMass * C_AtomicMass
                    Cohort%Processes%FilteredFood%P    = Cohort%Processes%FilteredFood%C * RATIOPC / P_AtomicMass * C_AtomicMass

                    !Clearance rate, m3/d        
                    Cohort%Processes%ClearanceRate      = 0.0  

                    !Ingestion rate, molC/d        
                    Cohort%Processes%IngestionInorganic = 0.0  
                    Cohort%Processes%IngestionFood%C    = Cohort%Processes%FilteredFood%C
                    Cohort%Processes%IngestionFood%H    = Cohort%Processes%FilteredFood%H
                    Cohort%Processes%IngestionFood%O    = Cohort%Processes%FilteredFood%O
                    Cohort%Processes%IngestionFood%N    = Cohort%Processes%FilteredFood%N
                    Cohort%Processes%IngestionFood%P    = Cohort%Processes%FilteredFood%P

                    !Pseudo-faeces contribution rate, molC/d        
                    Cohort%Processes%PFContributionInorganic  = 0.0
                    Cohort%Processes%PFContributionFood%C     = 0.0
                    Cohort%Processes%PFContributionFood%H     = 0.0
                    Cohort%Processes%PFContributionFood%O     = 0.0
                    Cohort%Processes%PFContributionFood%N     = 0.0
                    Cohort%Processes%PFContributionFood%P     = 0.0

                    !Compute Assimilation, molC/d
                    Cohort%Processes%Assimilation%C    = F_FIX * PAM_FIX / mu_E * TempCorrection * Vol**(2.0/3.0)
                    Cohort%Processes%Assimilation%H    = Cohort%Processes%Assimilation%C * RATIOHC / H_AtomicMass * C_AtomicMass
                    Cohort%Processes%Assimilation%O    = Cohort%Processes%Assimilation%C * RATIOOC / O_AtomicMass * C_AtomicMass
                    Cohort%Processes%Assimilation%N    = Cohort%Processes%Assimilation%C * RATIONC / N_AtomicMass * C_AtomicMass
                    Cohort%Processes%Assimilation%P    = Cohort%Processes%Assimilation%C * RATIOPC / P_AtomicMass * C_AtomicMass

                    !Compute Faecesproduction, molC/d
                    Cohort%Processes%FaecesContributionInorganic = 0.0                                
                    Cohort%Processes%FaecesContributionFood%C    = Cohort%Processes%IngestionFood%C - &
                                                                   Cohort%Processes%Assimilation%C
                    Cohort%Processes%FaecesContributionFood%H    = Cohort%Processes%IngestionFood%H - &
                                                                   Cohort%Processes%Assimilation%H
                    Cohort%Processes%FaecesContributionFood%O    = Cohort%Processes%IngestionFood%O - &
                                                                   Cohort%Processes%Assimilation%O
                    Cohort%Processes%FaecesContributionFood%N    = Cohort%Processes%IngestionFood%N - &
                                                                   Cohort%Processes%Assimilation%N
                    Cohort%Processes%FaecesContributionFood%P    = Cohort%Processes%IngestionFood%P - &
                                                                   Cohort%Processes%Assimilation%P

                else ! if not (M_H .gt. MHb), dont feed

                    call ImposeNoFiltrationProcess (Cohort)

                end if !(M_H .gt. MHb) feeding

                Cohort%FeedingOn = 0.0

                Cohort => Cohort%Next

            end do d2

            Species => Species%Next

        end do d1

    end subroutine ComputeSimpleFiltration

    !--------------------------------------------------------------------------

    subroutine ImposeNoFiltrationProcess (Cohort)

        !Arguments-------------------------------------------------------------
        type(T_Cohort),       pointer :: Cohort

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        Cohort%Processes%ClearanceRate                 = 0.0
        Cohort%Processes%FilteredInorganic             = 0.0
        Cohort%Processes%FilteredFood%C                = 0.0
        Cohort%Processes%FilteredFood%H                = 0.0
        Cohort%Processes%FilteredFood%O                = 0.0
        Cohort%Processes%FilteredFood%N                = 0.0
        Cohort%Processes%FilteredFood%P                = 0.0

        Cohort%Processes%IngestionInorganic            = 0.0
        Cohort%Processes%IngestionFood%C               = 0.0
        Cohort%Processes%IngestionFood%H               = 0.0
        Cohort%Processes%IngestionFood%O               = 0.0
        Cohort%Processes%IngestionFood%N               = 0.0
        Cohort%Processes%IngestionFood%P               = 0.0

        Cohort%Processes%PFContributionInorganic       = 0.0
        Cohort%Processes%PFContributionFood%C          = 0.0
        Cohort%Processes%PFContributionFood%H          = 0.0
        Cohort%Processes%PFContributionFood%O          = 0.0
        Cohort%Processes%PFContributionFood%N          = 0.0
        Cohort%Processes%PFContributionFood%P          = 0.0
        Cohort%Processes%Assimilation%C                = 0.0
        Cohort%Processes%Assimilation%H                = 0.0
        Cohort%Processes%Assimilation%O                = 0.0
        Cohort%Processes%Assimilation%N                = 0.0
        Cohort%Processes%Assimilation%P                = 0.0

        Cohort%Processes%FaecesContributionInorganic   = 0.0
        Cohort%Processes%FaecesContributionFood%C      = 0.0
        Cohort%Processes%FaecesContributionFood%H      = 0.0
        Cohort%Processes%FaecesContributionFood%O      = 0.0
        Cohort%Processes%FaecesContributionFood%N      = 0.0
        Cohort%Processes%FaecesContributionFood%P      = 0.0

    end subroutine ImposeNoFiltrationProcess

    !--------------------------------------------------------------------------

    subroutine ComputeComplexFiltration(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: Index

        !Local-----------------------------------------------------------------
        type(T_Species),    pointer     :: Species

        !Begin-----------------------------------------------------------------


        Me%LackOfFood = 0.0      
        
        call ComputeClearanceRate (Index)  
        
        Species => Me%FirstSpecies
        do while(associated(Species))
        
            call ComputeFiltrationRate (Species, Index)

            call ComputeIngestionAssimilationRate (Species)

            Species => Species%Next
        end do 
        
        call UpdateLarvaeOnMatrixMass (Index)
                
        Species => Me%FirstSpecies
        do while(associated(Species))
        
            call UpdateMatrixMass (Species, Index)

            Species => Species%Next
        end do 

    end subroutine ComputeComplexFiltration

    !--------------------------------------------------------------------------

    subroutine ComputeClearanceRate(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: Index

        !Local-----------------------------------------------------------------
        type(T_Species)  ,    pointer   :: Species
        type(T_Species)  ,    pointer   :: SpeciesAgain
        type(T_Particles),    pointer   :: Particles
        type(T_Particles),    pointer   :: ParticlesAgain
        type(T_Cohort)   ,    pointer   :: Cohort
        integer                         :: M_H, L, Number, ParticlesIndex
        integer                         :: PON, POP, POC, BioSilica, diatoms, DiatomsSi 
        real                            :: Crm, JX1Fm, JX0Fm, ro_X1, ro_X0, JX1Im, JX0Im
        real                            :: TempCorrection, MHb, Vol
        real                            :: C_AtomicMass, H_AtomicMass, O_AtomicMass, N_AtomicMass, P_AtomicMass
        real                            :: CrDenominator!, ComputedYEX
        real                            :: CrDenominator_Withlarvae
        !Begin-----------------------------------------------------------------

        !Cycle to compute CRdenominator = 1 + sum(Xi*CRM/JXim), Xi in molC/m3
        !base on all the properties that are food for each cohort in each species
        
        POC           = Me%PropIndex%POC
        PON           = Me%PropIndex%PON
        POP           = Me%PropIndex%POP
        BioSilica     = Me%PropIndex%BioSilica
        diatoms       = Me%PropIndex%diatoms
        DiatomsSi     = Me%PropIndex%DiatomsSi


        Species => Me%FirstSpecies
d1:     do while(associated(Species))

            Crm            = Species%IndividualParameters%Crm 
            JX1Fm          = Species%IndividualParameters%JX1Fm      
            JX0Fm          = Species%IndividualParameters%JX0Fm      
            ro_X1          = Species%IndividualParameters%ro_X1      
            ro_X0          = Species%IndividualParameters%ro_X0      
            JX1Im          = Species%IndividualParameters%JX1Im      
            JX0Im          = Species%IndividualParameters%JX0Im

            TempCorrection = Species%AuxiliarParameters%TempCorrection         
            C_AtomicMass   = Species%AuxiliarParameters%C_AtomicMass        
            H_AtomicMass   = Species%AuxiliarParameters%H_AtomicMass        
            O_AtomicMass   = Species%AuxiliarParameters%O_AtomicMass        
            P_AtomicMass   = Species%AuxiliarParameters%P_AtomicMass        
            N_AtomicMass   = Species%AuxiliarParameters%N_AtomicMass        
            MHb            = Species%AuxiliarParameters%MHb  

            !At this point we dont know if there are cohorts filtering larvae, so we compute the 
            !denominator with and without larvae
            CrDenominator            = 1.0
            CrDenominator_Withlarvae = 1.0

            Particles => Species%FirstParticles
d3:         do while(associated(Particles))

                if (Particles%Larvae) then
                
                    CrDenominator_Withlarvae = CrDenominator_Withlarvae                           + &
                                               Particles%LarvaeBiomass * Particles%LarvaeDensity  * &
                                               Crm / JX1Fm
                                               
                else
                                
                    if (Particles%Organic .eq. 1) then

                        if (Particles%ID%Name .eq.'particulate organic matter') then

                            if (Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then

                                if  (Me%ComputeOptions%Nitrogen) then

                                    !Amount of POC based on the PON value, using the NCratio and converted to molC
                                    CrDenominator = CrDenominator                         + &
                                                    (Me%ExternalVar%Mass(PON,Index)       / &
                                                    Particles%Ratios%NC_Ratio)            / &
                                                    C_AtomicMass * Crm / JX1Fm

                                    CrDenominator_Withlarvae = CrDenominator_Withlarvae   + &
                                                    (Me%ExternalVar%Mass(PON,Index)       / &
                                                    Particles%Ratios%NC_Ratio)            / &
                                                    C_AtomicMass * Crm / JX1Fm

                                else

                                    !Only if the model is running just with P
                                    !Amount of POC based on the POP value, using the PCratio and converted to molC
                                    CrDenominator = CrDenominator                         + &
                                                    (Me%ExternalVar%Mass(POP,Index)       / &
                                                    Particles%Ratios%PC_Ratio)            / &
                                                    C_AtomicMass * Crm / JX1Fm

                                    CrDenominator_Withlarvae = CrDenominator_Withlarvae   + &
                                                     (Me%ExternalVar%Mass(POP,Index)      / &
                                                     Particles%Ratios%PC_Ratio) / &
                                                     C_AtomicMass * Crm / JX1Fm

                                end if

                            else
                            ! if Life

                                !For filtration
                                !Amount of POC 
                                CrDenominator = CrDenominator                               + & 
                                                Me%ExternalVar%Mass(POC,Index)/C_AtomicMass * &
                                                Crm / JX1Fm  
                                
                                CrDenominator_Withlarvae = CrDenominator_Withlarvae                    + &
                                                           Me%ExternalVar%Mass(POC,Index)/C_AtomicMass * &
                                                           Crm / JX1Fm  

                            end if ! Me%ComputeOptions%PelagicModel

                        else
                        !if not 'particulate organic matter' 

                            if (Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then

                                ParticlesIndex = SearchPropIndex(GetPropertyIDNumber(Particles%ID%Name))

                            else
                            !if Life

                                ParticlesIndex = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" carbon")))

                            end if

                            CrDenominator = CrDenominator                                          + &
                                            Me%ExternalVar%Mass(ParticlesIndex,Index)/C_AtomicMass * &
                                            Crm / JX1Fm

                            CrDenominator_Withlarvae = CrDenominator_Withlarvae                               + &
                                                       Me%ExternalVar%Mass(ParticlesIndex,Index)/C_AtomicMass * &
                                                       Crm / JX1Fm

                        end if !'particulate organic matter'   

                    else !if not organic

                        ParticlesIndex = SearchPropIndex(GetPropertyIDNumber(Particles%ID%Name))

                        CrDenominator = CrDenominator                             + &
                                        Me%ExternalVar%Mass(ParticlesIndex,Index) * &
                                        Crm / JX0Fm

                        CrDenominator_Withlarvae = CrDenominator_Withlarvae                 + &
                                                   Me%ExternalVar%Mass(ParticlesIndex,Index)* &
                                                   Crm / JX0Fm

                    end if !(Particles%Organic .eq. 1)
                
                end if !(Particles%Larvae) 

                Particles => Particles%Next

            end do d3

            !Cycle to compute the Clearance rate, m3/d/individual for each cohort in each species
            !And the total CR foe each species, with and without larvae
            Cohort => Species%FirstCohort
d2:         do while(associated(Cohort))

                M_H    = Cohort%StateIndex%M_H    
                L      = Cohort%StateIndex%L    
                Number = Cohort%StateIndex%Number    
                Vol    = Cohort%BivalveCondition%Vol        

                if (Me%ExternalVar%Mass(M_H,Index) .gt. MHb) then !feeding
                
                    if ((Species%FeedOnLarvae) .and. (.not.(Cohort%Larvae))) then
                    !this cohort will feed on larvae and clearance rate computed using the denominator with larvae
                        
                        Cohort%Processes%ClearanceRate = Crm * Tempcorrection * Vol**(2.0/3.0) / CrDenominator_Withlarvae 

                        !Total Clearance rate, m3/d (.m3), sum all individuals in the system that are able to filter larvae
                        !m3/d (.m3) = m3/d.ind * ind/m3                           
                        Species%Total_CR_Larvae = Species%Total_CR_Larvae             + &
                                                  Cohort%Processes%ClearanceRate      * &
                                                  Me%ExternalVar%Mass(Number,Index)

                    else 
                    
                        Cohort%Processes%ClearanceRate = Crm * Tempcorrection * Vol**(2.0/3.0) / CrDenominator 
                    
                    end if !(Species%FeedOnLarvae)            

                    !Total Clearance rate, m3/d (.m3), sum all individuals in the system, for the other particles
                    !m3/d (.m3) = m3/d.ind * ind/m3
                    Species%Total_CR = Species%Total_CR               + &
                                       Cohort%Processes%ClearanceRate * &
                                       Me%ExternalVar%Mass(Number,Index)

                else ! if not (M_H .gt. MHb)

                    Cohort%Processes%ClearanceRate = 0.0

                end if !(M_H .gt. MHb) feeding

                !Initialize to proceed with computation       
                Cohort%Processes%FilteredInorganic           = 0.0
                Cohort%Processes%FilteredFood%C              = 0.0
                Cohort%Processes%FilteredFood%H              = 0.0
                Cohort%Processes%FilteredFood%O              = 0.0
                Cohort%Processes%FilteredFood%N              = 0.0
                Cohort%Processes%FilteredFood%P              = 0.0

                Cohort%Processes%IngestionInorganic          = 0.0
                Cohort%Processes%IngestionFood%C             = 0.0
                Cohort%Processes%IngestionFood%H             = 0.0
                Cohort%Processes%IngestionFood%O             = 0.0
                Cohort%Processes%IngestionFood%N             = 0.0
                Cohort%Processes%IngestionFood%P             = 0.0

                Cohort%Processes%PFContributionInorganic     = 0.0
                Cohort%Processes%PFContributionFood%C        = 0.0
                Cohort%Processes%PFContributionFood%H        = 0.0
                Cohort%Processes%PFContributionFood%O        = 0.0
                Cohort%Processes%PFContributionFood%N        = 0.0
                Cohort%Processes%PFContributionFood%P        = 0.0

                Cohort%Processes%Assimilation%C              = 0.0
                Cohort%Processes%Assimilation%H              = 0.0
                Cohort%Processes%Assimilation%O              = 0.0
                Cohort%Processes%Assimilation%N              = 0.0
                Cohort%Processes%Assimilation%P              = 0.0

                Cohort%Processes%FaecesContributionInorganic = 0.0
                Cohort%Processes%FaecesContributionFood%C    = 0.0
                Cohort%Processes%FaecesContributionFood%H    = 0.0
                Cohort%Processes%FaecesContributionFood%O    = 0.0
                Cohort%Processes%FaecesContributionFood%N    = 0.0
                Cohort%Processes%FaecesContributionFood%P    = 0.0

                Cohort => Cohort%Next

            end do d2
            
            Species => Species%Next

        end do d1
        
        !Cycle to compute the total CR by particle (all species) for future compute of filtration
        Species => Me%FirstSpecies
        do while(associated(Species))
        
            Particles => Species%FirstParticles
            do while(associated(Particles))
            
                SpeciesAgain => Me%FirstSpecies
                do while(associated(SpeciesAgain))
                
                    ParticlesAgain => SpeciesAgain%FirstParticles
                    do while(associated(ParticlesAgain))
                    
                        if (ParticlesAgain%ID%Name .eq. Particles%ID%Name) then
                        
                            if (ParticlesAgain%Larvae) then
                                Particles%Total_CR = Particles%Total_CR + SpeciesAgain%Total_CR_Larvae
                            else
                                Particles%Total_CR = Particles%Total_CR + SpeciesAgain%Total_CR
                            end if
                        
                        end if 
                        
                        ParticlesAgain => ParticlesAgain%Next

                    end do
                
                    SpeciesAgain => SpeciesAgain%Next

                end do
            
                Particles => Particles%Next

            end do        
        
            Species => Species%Next

        end do 
        
        
    end subroutine ComputeClearanceRate

    !--------------------------------------------------------------------------

    subroutine ComputeFiltrationRate(Species, Index)

        !Arguments-------------------------------------------------------------
        type(T_Species)    , pointer    :: Species
        integer, intent(IN)             :: Index

        !Local-----------------------------------------------------------------
        type(T_Particles)  , pointer    :: Particles
        type(T_Cohort)     , pointer    :: Cohort
        integer                         :: L, Number, ParticlesIndex
        integer                         :: par
        integer                         :: PON, POP, POC, BioSilica, diatoms, DiatomsSi 
        real                            :: ParticleConcentration, Total_ParticlePotFil
        real                            :: ParticleTempMass, Total_PossibleParticleFil
        real                            :: FilteredByCohort
        real                            :: C_AtomicMass

        !Begin-----------------------------------------------------------------

        POC           = Me%PropIndex%POC
        PON           = Me%PropIndex%PON
        POP           = Me%PropIndex%POP
        BioSilica     = Me%PropIndex%BioSilica
        diatoms       = Me%PropIndex%diatoms
        DiatomsSi     = Me%PropIndex%DiatomsSi
        
        C_AtomicMass  = Species%AuxiliarParameters%C_AtomicMass        

        par = 0

        !search for particle concentration, in molC/d if organic, g/d if inorganic
        Particles => Species%FirstParticles
        do while(associated(Particles))

            par = par + 1 !count which property

            ParticleConcentration = 0.0
            
            if (Particles%Larvae) then
            
                !total molC/m3 = molC/# * #/m3 
                ParticleConcentration = Particles%LarvaeBiomass * Particles%LarvaeDensity
                
            else

                if (Particles%Organic .eq. 1) then

                    if (Particles%ID%Name .eq.'particulate organic matter') then

                        !POMcheck = 1

                        if (Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then

                            if (Me%ComputeOptions%Nitrogen) then

                                !amount of carbon associated with PON property, em molC/m3
                                !mgN/l / mgN/mgC / g/mol = mol/m3 . 
                                ParticleConcentration = (Me%ExternalVar%Mass(PON,Index)             / &
                                                         Particles%Ratios%NC_Ratio)                 / &
                                                         C_AtomicMass
                            else 
                                !Only if the model is running just with P
                                !carbon associated with POP property, mol/d            
                                ParticleConcentration = (Me%ExternalVar%Mass(POP,Index)             / &
                                                         Particles%Ratios%PC_Ratio)                 / &
                                                         C_AtomicMass
                            end if

                        else
                        !Life

                            !amount of carbon associated with PON property, em mol/d
                            ParticleConcentration = Me%ExternalVar%Mass(POC,Index) / C_AtomicMass             
                        
                        end if 

                    else         
                    !if not 'particulate organic matter'            

                        if (Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then
                            ParticlesIndex = SearchPropIndex(GetPropertyIDNumber(Particles%ID%Name))
                        else
                            !if Life
                            ParticlesIndex = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" carbon")))
                        end if

                        !molC/m3 = g/m3 * molC/g 
                        ParticleConcentration = Me%ExternalVar%Mass(ParticlesIndex,Index) / C_AtomicMass 
                        
                    end if  !if 'particulate organic matter'     

                else !not organic

                    ParticlesIndex = SearchPropIndex(GetPropertyIDNumber(Particles%ID%Name))

                    !g/m3
                    ParticleConcentration = Me%ExternalVar%Mass(ParticlesIndex,Index)

                end if !organic
                
                !Potential Total filtration of this particle molC(g)/d(.m3) = m3/d(.m3) * molC/m3
                !Total_CR -  soma de todos os cohorts que filtram cada propriedade
                
            end if !larvae
            
            Total_ParticlePotFil = Particles%Total_CR * ParticleConcentration

            !check if it is possible, in mass
            !molC.m3 = molC/m3 - molC/d(.m3) * d
            ParticleTempMass = ParticleConcentration - Total_ParticlePotFil * Me%DTDay 

            if ((ParticleTempMass .lt. 0.0) .or. (ParticleConcentration .eq. 0.0))then

                Me%LackOfFood = 1.0

                !molC/d (.m3) = molC/m3 / d
                Total_PossibleParticleFil = ParticleConcentration / Me%DTDay

            else

                !molC/d (.m3)
                Total_PossibleParticleFil = Total_ParticlePotFil

            endif

            !Compute filtration based on the concentration and Total_PossibleParticleFil found before
            Cohort => Species%FirstCohort
            do while(associated(Cohort))
            
                if (Cohort%Dead .eq. 0) then


                    Number = Cohort%StateIndex%Number    
                    L      = Cohort%StateIndex%L 
                                        
                    if (Particles%Larvae .and. Cohort%Larvae) then 
                    !this cohort does not feed on larvae if it is a larvae
                        
                        FilteredByCohort = 0.0
                        
                    else
                    
                        if (Particles%Total_CR .ne. 0.0) then

                            if (Me%ComputeOptions%CorrectFiltration) then
                            
                                if (Me%ExternalVar%Mass(Number,Index) .ne. 0.0) then
                                
                                    !molC/d.ind
                                    FilteredByCohort = Cohort%Processes%ClearanceRate                         * &
                                                       Me%ExternalVar%Mass(Number,Index) / Particles%Total_CR * &
                                                       Total_PossibleParticleFil                              / &
                                                       Me%ExternalVar%Mass(Number,Index)
                                
                                else
                                
                                    FilteredByCohort = 0.0
                                    
                                end if
                                
                            else 
                            !molC/d.ind
                                FilteredByCohort = Cohort%Processes%ClearanceRate * ParticleConcentration 

                            end if ! if (Me%ComputeOptions%CorrectFiltration)
                            
                        else !(Total_CR .ne. 0.0) 
                        
                            FilteredByCohort = 0.0
                        
                        end if !(Total_CR_Forlarvae .ne. 0.0)

                    end if

                    !store value in the Filtered column, in molC or molC/d.ind
                    Cohort%FeedingOn(par,1) = FilteredByCohort

                    !Compute the sum of properties for each cohort
                    if (Particles%Organic .eq. 1) then

                        !Total filtration from organic material, molC/d
                        Cohort%Processes%FilteredFood%C = Cohort%Processes%FilteredFood%C + FilteredByCohort

                        Cohort%Processes%FilteredFood%H = Cohort%Processes%FilteredFood%H +                  &
                                                          FilteredByCohort * Particles%Composition%nH
                                                          
                        Cohort%Processes%FilteredFood%O = Cohort%Processes%FilteredFood%O +                  &
                                                          FilteredByCohort * Particles%Composition%nO
                                                          
                        Cohort%Processes%FilteredFood%N = Cohort%Processes%FilteredFood%N +                  &
                                                          FilteredByCohort * Particles%Composition%nN
                                                          
                        Cohort%Processes%FilteredFood%P = Cohort%Processes%FilteredFood%P +                  &
                                                          FilteredByCohort * Particles%Composition%nP

                    else !inorganic

                        !Total filtration from inorganic material, g/d
                        Cohort%Processes%FilteredInorganic = Cohort%Processes%FilteredInorganic + FilteredByCohort

                    end if
                    
                endif !Cohort%Dead

                Cohort => Cohort%Next
            end do 
            Particles => Particles%Next
        end do  !to have the total carbon Filtered by each cohort from all properties
        
    end subroutine ComputeFiltrationRate
    
    !--------------------------------------------------------------------------

    subroutine ComputeIngestionAssimilationRate(Species)

        !Arguments-------------------------------------------------------------
        type(T_Species),      pointer        :: Species

        !Local-----------------------------------------------------------------
        type(T_Cohort)   ,       pointer     :: Cohort
        type(T_Particles),       pointer     :: Particles
        real                                 :: ro_X1, ro_X0, JX1Im, JX0Im, YEX
        integer                              :: par
        real                                 :: IngDenominator !, ComputedYEX
        real                                 :: FilteredByCohort,IngestedByCohort, AssimilatedByCohort
        real                                 :: r_C, r_N, r_P
        real                                 :: AssimilatedStructure, AssimilatedReserves

        !Begin-----------------------------------------------------------------

        ro_X1          = Species%IndividualParameters%ro_X1      
        ro_X0          = Species%IndividualParameters%ro_X0      
        JX1Im          = Species%IndividualParameters%JX1Im      
        JX0Im          = Species%IndividualParameters%JX0Im
        YEX            = Species%IndividualParameters%YEX

        !Compute ingestion and assimilation, molC/d/individual
        Cohort => Species%FirstCohort
        do while(associated(Cohort))
        
            if (Cohort%Dead .eq. 0) then

                !Because ro_Xi and JXiIm are the same for all organic particles
                IngDenominator = 1 + (ro_X1 * Cohort%Processes%FilteredFood%C)    / JX1Im    + &
                                     (ro_X0 * Cohort%Processes%FilteredInorganic) / JX0Im

                par = 0

                Particles => Species%FirstParticles
                do while(associated(Particles)) !round to compute IngestedByCohort (second column of FeedingOn)

                    par = par + 1

                    FilteredByCohort = Cohort%FeedingOn(par,1)

                    if (Particles%Organic .eq. 1) then

                        !Ingestion from each food item, molC/d(.m3) 
                        IngestedByCohort = ro_X1 * FilteredByCohort / IngDenominator 

                        !Ingestion: sum of all ingested carbon form the different food itens, mol/d
                        Cohort%Processes%IngestionFood%C = Cohort%Processes%IngestionFood%C + IngestedByCohort

                        Cohort%Processes%IngestionFood%H = Cohort%Processes%IngestionFood%H                 + &
                                                            IngestedByCohort * Particles%Composition%nH

                        Cohort%Processes%IngestionFood%O = Cohort%Processes%IngestionFood%O                 + &
                                                            IngestedByCohort * Particles%Composition%nO

                        Cohort%Processes%IngestionFood%N = Cohort%Processes%IngestionFood%N                 + &
                                                            IngestedByCohort * Particles%Composition%nN

                        Cohort%Processes%IngestionFood%P = Cohort%Processes%IngestionFood%P                 + &
                                                            IngestedByCohort * Particles%Composition%nP

                        !Pseudofaeces: sum of all contributions form all the different food itens, mol/d(.ind)
                        Cohort%Processes%PFContributionFood%C = Cohort%Processes%PFContributionFood%C       + &
                                                            (FilteredByCohort - IngestedByCohort)

                        Cohort%Processes%PFContributionFood%H = Cohort%Processes%PFContributionFood%H       + &
                                                            (FilteredByCohort - IngestedByCohort) * Particles%Composition%nH

                        Cohort%Processes%PFContributionFood%O = Cohort%Processes%PFContributionFood%O       + &
                                                            (FilteredByCohort - IngestedByCohort) * Particles%Composition%nO

                        Cohort%Processes%PFContributionFood%N = Cohort%Processes%PFContributionFood%N       + &
                                                            (FilteredByCohort - IngestedByCohort) * Particles%Composition%nN 

                        Cohort%Processes%PFContributionFood%P = Cohort%Processes%PFContributionFood%P       + &
                                                            (FilteredByCohort - IngestedByCohort) * Particles%Composition%nP

                        if (Species%IndividualParameters%SIMPLE_ASSI .eq. 1) then !simple assimilation model

                            !Assimilation, molC/d
                            AssimilatedByCohort = IngestedByCohort * YEX

                        else !not SIMPLE_ASSI

                            if (IngestedByCohort .gt. 1e-20) then

                                if ((Particles%ID%Name .eq. 'particulate organic matter')) then 

                                    Particles%F_E = 1

                                end if 

                                !Assimilation: Structure, simple relation using the yield coefficient
                                AssimilatedStructure = IngestedByCohort * (1-Particles%F_E) * YEX

                                !Assimilation: Reserves, Parallel and complementary SU's for three elements, molC/d
                                r_C = IngestedByCohort * Particles%F_E

                                r_N = IngestedByCohort * Particles%F_E                    * &
                                        Particles%Composition%nN                          / &
                                        Species%SpeciesComposition%ReservesComposition%nN

                                r_P = IngestedByCohort * Particles%F_E                     * &
                                        Particles%Composition%nP                           / &
                                        Species%SpeciesComposition%ReservesComposition%nP

                                AssimilatedReserves = ( 1/r_C                       + &
                                                        1/r_N                       + &
                                                        1/r_P                       - &
                                                        1/(r_C + r_N)               - &
                                                        1/(r_C + r_P)               - &
                                                        1/(r_N + r_P)               + &
                                                        1/(r_C + r_N + r_P))          &
                                                        **(-1)

                                !Assimilation: Reserves + Structure, molC/d       
                                AssimilatedByCohort = AssimilatedStructure + AssimilatedReserves       

                            else
                            !Assimilation: Reserves + Structure, molC/d(.m3)       
                            AssimilatedByCohort = 0.0       

                            end if !if ingestion !=0.0

                        end if !if SIMPLE_ASSI

                        !ComputedYEX = Cohort%Processes%Assimilation%C / Cohort%Processes%IngestionFood%C 
                        
                        !Total Assimilation in C, mol/d(.m3)
                        Cohort%Processes%Assimilation%C = Cohort%Processes%Assimilation%C + AssimilatedByCohort

                        Cohort%Processes%Assimilation%H = Cohort%Processes%Assimilation%C                  * &
                                                          Species%SpeciesComposition%ReservesComposition%nH

                        Cohort%Processes%Assimilation%O = Cohort%Processes%Assimilation%C                  * &
                                                          Species%SpeciesComposition%ReservesComposition%nO

                        Cohort%Processes%Assimilation%N = Cohort%Processes%Assimilation%C                  * &
                                                          Species%SpeciesComposition%ReservesComposition%nN

                        Cohort%Processes%Assimilation%P = Cohort%Processes%Assimilation%C                  * &
                                                          Species%SpeciesComposition%ReservesComposition%nP

                        !Total Faeces Contribution in C, mol/d
                        Cohort%Processes%FaecesContributionFood%C = Cohort%Processes%FaecesContributionFood%C + &
                                                                    (IngestedByCohort - AssimilatedByCohort) 
                        
                        !Total Faeces production for the other elements, mol/d
                        Cohort%Processes%FaecesContributionFood%H = Cohort%Processes%FaecesContributionFood%H  + &
                                                               (IngestedByCohort * Particles%Composition%nN    - & 
                                                                AssimilatedByCohort                            * &
                                                                Species%SpeciesComposition%ReservesComposition%nN)
                                                                    
                        Cohort%Processes%FaecesContributionFood%O = Cohort%Processes%FaecesContributionFood%O  + &
                                                               (IngestedByCohort * Particles%Composition%nO    - & 
                                                                AssimilatedByCohort                            * &
                                                                Species%SpeciesComposition%ReservesComposition%nO)

                        Cohort%Processes%FaecesContributionFood%N = Cohort%Processes%FaecesContributionFood%N  + &
                                                               (IngestedByCohort * Particles%Composition%nN    - & 
                                                                AssimilatedByCohort                            * &
                                                                Species%SpeciesComposition%ReservesComposition%nN)

                        Cohort%Processes%FaecesContributionFood%P = Cohort%Processes%FaecesContributionFood%P  + &
                                                               (IngestedByCohort * Particles%Composition%nP    - & 
                                                                AssimilatedByCohort                            * &
                                                                Species%SpeciesComposition%ReservesComposition%nP)

                    else !not organic

                        !Ingestion, g/d
                        IngestedByCohort = ro_X0 * FilteredByCohort / IngDenominator  

                        !Ingestion: sum of all ingested mass of inorganic material that the mussel is able to ingest, g/d
                        Cohort%Processes%IngestionInorganic = Cohort%Processes%IngestionInorganic + IngestedByCohort

                        !Pseudofaeces: sum of all pseudofaeces from inorganic material, g/d
                        Cohort%Processes%PFContributionInorganic = Cohort%Processes%PFContributionInorganic  + &
                                                                   (FilteredByCohort - IngestedByCohort)

                        !Assimilation: inorganic material can not be assimilated
                        AssimilatedByCohort = 0.0

                        !Faeces: sum of all faeces from inorganic material, g/d
                        Cohort%Processes%FaecesContributionInorganic = Cohort%Processes%FaecesContributionInorganic +   &
                                                                       (IngestedByCohort - AssimilatedByCohort)

                    end if !organic

                    !store value in the ingested column
                    Cohort%FeedingOn(par,2) = IngestedByCohort

                    !store value in the assimilated column
                    Cohort%FeedingOn(par,3) = AssimilatedByCohort

                    Particles => Particles%Next

                end do     
            endif !Cohort%Dead
            Cohort => Cohort%Next
        end do 
        
    end subroutine ComputeIngestionAssimilationRate
    
    !--------------------------------------------------------------------------

    subroutine UpdateLarvaeOnMatrixMass(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: Index

        !Local-----------------------------------------------------------------
        type(T_Species)      , pointer      :: Species
        type(T_Species)      , pointer      :: SpeciesAgain
        type(T_Particles)   , pointer       :: Particles
        type(T_Cohort)      , pointer       :: Cohort
        type(T_Cohort)      , pointer       :: CohortAgain
        integer                             :: PON, POP, POC 
        integer                             :: Number, par,DynamicParticlesIndex, NumberAgain
        real                                :: FilteredByCohort,IngestedByCohort, AssimilatedByCohort
        real                                :: PseudoFaecesByCohort,FaecesByCohort,FaecesByCohortN, FaecesByCohortP
        real                                :: C_AtomicMass, H_AtomicMass, O_AtomicMass, N_AtomicMass, P_AtomicMass

        
        !Begin-----------------------------------------------------------------

        POC           = Me%PropIndex%POC
        PON           = Me%PropIndex%PON
        POP           = Me%PropIndex%POP
               
        Species => Me%FirstSpecies
        do while(associated(Species))
        
            C_AtomicMass   = Species%AuxiliarParameters%C_AtomicMass        
            H_AtomicMass   = Species%AuxiliarParameters%H_AtomicMass        
            O_AtomicMass   = Species%AuxiliarParameters%O_AtomicMass        
            P_AtomicMass   = Species%AuxiliarParameters%P_AtomicMass        
            N_AtomicMass   = Species%AuxiliarParameters%N_AtomicMass        

            if (Species%FeedOnLarvae) then
            
                par = 0
                !Update the number of larvae first
                Particles => Species%FirstParticles
                do while(associated(Particles))

                    par = par + 1

                    FilteredByCohort    = 0.0
                    IngestedByCohort    = 0.0
                    AssimilatedByCohort = 0.0
                    
                    if (Particles%Larvae) then

                        Cohort => Species%FirstCohort
                        do while(associated(Cohort))

                            Number = Cohort%StateIndex%Number

                            !molC/d
                            FilteredByCohort     = Cohort%FeedingOn(par,1)

                            if (FilteredByCohort .gt. 0.0) then

                                !molC/d
                                IngestedByCohort     = Cohort%FeedingOn(par,2)

                                !molC/d
                                AssimilatedByCohort  = Cohort%FeedingOn(par,3)

                                !Pseudofaeces and Faeces form this species on this particles
                                PseudoFaecesByCohort = FilteredByCohort - IngestedByCohort

                                FaecesByCohort = IngestedByCohort - AssimilatedByCohort

                                FaecesByCohortN = IngestedByCohort * Particles%Composition%nN                            - & 
                                                  AssimilatedByCohort * Species%SpeciesComposition%ReservesComposition%nN

                                FaecesByCohortP = IngestedByCohort * Particles%Composition%nP                            - & 
                                                  AssimilatedByCohort * Species%SpeciesComposition%ReservesComposition%nP

                                !find the particle cohort to update the number of larvae
                                SpeciesAgain => Me%FirstSpecies
                                do while(associated(SpeciesAgain))
                                
                                    if (Particles%LarvaeSpeciesID .eq. SpeciesAgain%ID%IDNumber) then 
                                    
                                        DynamicParticlesIndex = SearchPropIndex(Particles%ID%IDNumber)
                               
                                        CohortAgain => SpeciesAgain%FirstCohort
                                        do while(associated(CohortAgain))
                                        
                                            NumberAgain = CohortAgain%StateIndex%Number

                                            if(DynamicParticlesIndex == CohortAgain%StateIndex%Number) then
                                            !cohort found
                                            
                                            
                                                !Update number in the matrix mass
                                                Me%ExternalVar%Mass(NumberAgain, Index) = Me%ExternalVar%Mass(NumberAgain, Index)-&
                                                                                    (FilteredByCohort                           * &
                                                                                     Me%ExternalVar%Mass(Number, Index)         / &
                                                                                     Particles%LarvaeBiomass)                   * &
                                                                                    Me%DTDay
                                                                                                     
                                                if (Particles%LarvaeSpeciesID .eq. Species%ID%IDNumber) then 
                                                
                                                    ! #/m3.d, store the total from all the cohorts
                                                    CohortAgain%Processes%DeathBySelfPredation =                             &
                                                                                CohortAgain%Processes%DeathBySelfPredation + &
                                                                                (FilteredByCohort                          * &
                                                                                Me%ExternalVar%Mass(Number, Index)         / &
                                                                                Particles%LarvaeBiomass) 

                                                else    
                                                    
                                                    CohortAgain%Processes%DeathByLarvaePredationByOthers =                     &
                                                                        CohortAgain%Processes%DeathByLarvaePredationByOthers + &
                                                                                (FilteredByCohort                            * &
                                                                                Me%ExternalVar%Mass(Number, Index)           / &
                                                                                Particles%LarvaeBiomass) 
                                                    
                                                end if
                                                
                                             end if
                                                  
                                             CohortAgain => CohortAgain%Next
                                        end do
                                        
                                    end if 
                                   
                                    SpeciesAgain => SpeciesAgain%Next
                                end do 
                                !update PON, POC and POC 

                                if(Me%ComputeOptions%Nitrogen)then
                                !g

                                    Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)                   + &
                                                                      ( PseudoFaecesByCohort * Particles%Composition%nN   + &
                                                                        FaecesByCohortN )                                 * &
                                                                       N_AtomicMass                                       * &
                                                                       Me%ExternalVar%Mass(Number, Index) * Me%DTDay
                               
                                end if   

                                if(Me%ComputeOptions%Phosphorus) then       

                                    Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)                      + &
                                                                      ( PseudoFaecesByCohort * Particles%Composition%nP      + &
                                                                        FaecesByCohortP )                                    * &
                                                                       P_AtomicMass                                          * &
                                                                       Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                                end if   

                                if(Me%ComputeOptions%PelagicModel .eq. LifeModel) then    

                                    Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)              + &
                                                                       ( PseudoFaecesByCohort                        + &
                                                                         FaecesByCohort      )                       * &
                                                                       C_AtomicMass                                  * &
                                                                       Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                                end if
                                
                            end if ! if (FilteredByCohort .gt. 0.0)

                            Cohort => Cohort%Next
                        end do
                        
                    end if !(if larvae)
                    Particles => Particles%Next
                end do             
            end if 

        Species => Species%Next
        end do 

    end subroutine UpdateLarvaeOnMatrixMass
    
    !--------------------------------------------------------------------------

    subroutine UpdateMatrixMass(Species, Index)

        !Arguments-------------------------------------------------------------
        type(T_Species)      , pointer      :: Species
        integer, intent(IN)                 :: Index

        !Local-----------------------------------------------------------------
        type(T_Particles)   , pointer       :: Particles
        type(T_Cohort)      , pointer       :: Cohort
        integer                             :: PON, POP, POC, BioSilica, diatoms, DiatomsSi 
        integer                             :: Number, par, POMcheck
        integer                             :: ParticlesIndex
        integer                             :: PropertyIndexC,PropertyIndexN, PropertyIndexP, PropertyIndexChl  
        real                                :: FilteredByCohort,IngestedByCohort, AssimilatedByCohort
        real                                :: PseudoFaecesByCohort,FaecesByCohort,FaecesByCohortN, FaecesByCohortP
        real                                :: C_AtomicMass, H_AtomicMass, O_AtomicMass, N_AtomicMass, P_AtomicMass

        
        !Begin-----------------------------------------------------------------

        POC           = Me%PropIndex%POC
        PON           = Me%PropIndex%PON
        POP           = Me%PropIndex%POP
        BioSilica     = Me%PropIndex%BioSilica
        diatoms       = Me%PropIndex%diatoms
        DiatomsSi     = Me%PropIndex%DiatomsSi
        
        C_AtomicMass   = Species%AuxiliarParameters%C_AtomicMass        
        H_AtomicMass   = Species%AuxiliarParameters%H_AtomicMass        
        O_AtomicMass   = Species%AuxiliarParameters%O_AtomicMass        
        P_AtomicMass   = Species%AuxiliarParameters%P_AtomicMass        
        N_AtomicMass   = Species%AuxiliarParameters%N_AtomicMass        

        par = 0
        
        POMcheck = 0
        
        Particles => Species%FirstParticles
        do while(associated(Particles))
        
            par = par + 1
            
            if (.not.(Particles%Larvae)) then
            !Larvae have been uptdated before

                FilteredByCohort    = 0.0
                IngestedByCohort    = 0.0
                AssimilatedByCohort = 0.0

                Cohort => Species%FirstCohort
                do while(associated(Cohort))

                    Number = Cohort%StateIndex%Number

                    !molC/d
                    FilteredByCohort     = Cohort%FeedingOn(par,1)

                    !molC/d
                    IngestedByCohort     = Cohort%FeedingOn(par,2)

                    !molC/d
                    AssimilatedByCohort  = Cohort%FeedingOn(par,3)

                    !Pseudofaeces and Faeces form this species on this particles
                    PseudoFaecesByCohort = FilteredByCohort - IngestedByCohort

                    FaecesByCohort = IngestedByCohort - AssimilatedByCohort

                    FaecesByCohortN = IngestedByCohort * Particles%Composition%nN                            - & 
                                      AssimilatedByCohort * Species%SpeciesComposition%ReservesComposition%nN

                    FaecesByCohortP = IngestedByCohort * Particles%Composition%nP                            - & 
                                      AssimilatedByCohort * Species%SpeciesComposition%ReservesComposition%nP

                    if (Particles%Organic .eq. 0) then !inorganic particle

                        ParticlesIndex = SearchPropIndex(GetPropertyIDNumber(Particles%ID%Name))

                        Me%ExternalVar%Mass (ParticlesIndex, Index) = Me%ExternalVar%Mass (ParticlesIndex, Index)   + &
                                                                     ( PseudoFaecesByCohort                         + &
                                                                       FaecesByCohort                               - &
                                                                       FilteredByCohort)                            * &
                                                                       Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                    else !(Particles%Organic .eq. 0)
                    
                        if (Particles%ID%Name .eq.'particulate organic matter') then
                        
                            POMcheck = 1

                            if(Me%ComputeOptions%Nitrogen)then

                                Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)                   + &
                                                                  ( PseudoFaecesByCohort * Particles%Composition%nN   + &
                                                                    FaecesByCohortN                                   - &
                                                                    FilteredByCohort     * Particles%Composition%nN ) * &
                                                                   N_AtomicMass                                       * &
                                                                   Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                            end if   

                            if(Me%ComputeOptions%Phosphorus) then       

                                Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)                   + &
                                                                  ( PseudoFaecesByCohort * Particles%Composition%nP   + &
                                                                    FaecesByCohortP                                   - &
                                                                    FilteredByCohort     * Particles%Composition%nP ) * &
                                                                   P_AtomicMass                                       * &
                                                                   Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                            end if   

                            if(Me%ComputeOptions%PelagicModel .eq. LifeModel) then    

                                Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)                + &
                                                            ( PseudoFaecesByCohort                                 + &
                                                              FaecesByCohort                                       - &
                                                              FilteredByCohort      )                              * &
                                                              C_AtomicMass                                         * &
                                                              Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                            end if

                        else
                        !if not 'particulate organic matter'

                            if(Me%ComputeOptions%PelagicModel .eq. WaterQualityModel ) then

                                if (Particles%ID%Name .eq.'diatoms') then

                                    if(Particles%Silica .eq. 1)then
                                    
                                        !All the silica ingested is returned to the water in the form of biosilica
                                        Me%ExternalVar%Mass (BioSilica, Index) = Me%ExternalVar%Mass(BioSilica, Index)+ &
                                                                        FilteredByCohort                              * &
                                                                        C_AtomicMass * Particles%Ratios%SiC_Ratio     * &
                                                                        Me%ExternalVar%Mass(Number, Index) * Me%DTDay
                                        
                                    end if

                                end if                          
                            
                                ParticlesIndex = SearchPropIndex(GetPropertyIDNumber(Particles%ID%Name))
                                
                                Me%ExternalVar%Mass (ParticlesIndex, Index) = Me%ExternalVar%Mass (ParticlesIndex, Index) - &
                                                                            FilteredByCohort                              * &
                                                                            Me%ExternalVar%Mass(Number, Index)            * &
                                                                            C_AtomicMass * Me%DTDay

                                !check if there is mass loss
                                if (Me%ExternalVar%Mass (ParticlesIndex, Index) .lt. 0.0) then
                                
                                    Me%MassLoss = Me%MassLoss + Me%ExternalVar%Mass (ParticlesIndex, Index) * (-1.)
                                    
                                    Me%ExternalVar%Mass (ParticlesIndex, Index) = 0.0
                                end if

                            else     
                            !if life 

                                if (Particles%ID%Name .eq.'diatoms') then             

                                    Me%ExternalVar%Mass (DiatomsSi, Index) = Me%ExternalVar%Mass (DiatomsSi, Index)     - &
                                                                            FilteredByCohort * C_AtomicMass             * &
                                                                            Particles%Ratios%SiC_Ratio                  * &
                                                                            Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                                    Me%ExternalVar%Mass (BioSilica, Index) = Me%ExternalVar%Mass (BioSilica, Index)     + &
                                                                            FilteredByCohort * C_AtomicMass             * &
                                                                            Particles%Ratios%SiC_Ratio                  * &
                                                                            Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                                end if

                                PropertyIndexC = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" carbon")))
                                PropertyIndexN = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" nitrogen")))
                                PropertyIndexP = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//" phosphorus")))

                                Me%ExternalVar%Mass (PropertyIndexC, Index) = Me%ExternalVar%Mass (PropertyIndexC, Index) - &
                                                                            FilteredByCohort * C_AtomicMass               * &
                                                                            Me%ExternalVar%Mass(Number, Index) * Me%DTDay 

                                Me%ExternalVar%Mass (PropertyIndexN, Index) = Me%ExternalVar%Mass (PropertyIndexN, Index) - &
                                                                            FilteredByCohort * Particles%Composition%nN   * &
                                                                            N_AtomicMass                                  * &
                                                                            Me%ExternalVar%Mass(Number, Index) * Me%DTDay 

                                Me%ExternalVar%Mass (PropertyIndexP, Index) = Me%ExternalVar%Mass (PropertyIndexP, Index) - &
                                                                            FilteredByCohort * Particles%Composition%nP   * & 
                                                                            P_AtomicMass                                  * &
                                                                            Me%ExternalVar%Mass(Number, Index) * Me%DTDay 


                                if ((Particles%ID%Name .eq. 'diatoms')                 .or. &
                                    (Particles%ID%Name .eq. 'autotrophic flagellates') .or. &
                                    (Particles%ID%Name .eq. 'picoalgae')               .or. &
                                    (Particles%ID%Name .eq. 'flagellates')) then 

                                    PropertyIndexChl = SearchPropIndex(GetPropertyIDNumber((trim(Particles%ID%Name)//   &
                                                        " chlorophyll")))

                                    Me%ExternalVar%Mass(PropertyIndexChl,Index)=Me%ExternalVar%Mass(PropertyIndexChl,Index) - &
                                                                                FilteredByCohort *  C_AtomicMass            * &
                                                                                Particles%Ratios%ChlC_Ratio                 * &
                                                                                Me%ExternalVar%Mass(Number, Index) * Me%DTDay  

                                end if

                            !if life
                            end if

                        end if !if Organic matter
                        
                    end if !(Particles%Organic .eq. 0)

                    !update PON, POC and POC when it is not a property from the properties list (does not have filtration)
                    if (POMcheck .eq. 0) then

                        if(Me%ComputeOptions%Nitrogen)then
                        !g

                            Me%ExternalVar%Mass (PON, Index) = Me%ExternalVar%Mass (PON, Index)                   + &
                                                              ( PseudoFaecesByCohort * Particles%Composition%nN   + &
                                                                FaecesByCohortN )                                 * &
                                                               N_AtomicMass                                       * &
                                                               Me%ExternalVar%Mass(Number, Index) * Me%DTDay
                       
                        end if   

                        if(Me%ComputeOptions%Phosphorus) then       

                            Me%ExternalVar%Mass (POP, Index) = Me%ExternalVar%Mass (POP, Index)                      + &
                                                              ( PseudoFaecesByCohort * Particles%Composition%nP      + &
                                                                FaecesByCohortP )                                    * &
                                                               P_AtomicMass                                          * &
                                                               Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                        end if   

                        if(Me%ComputeOptions%PelagicModel .eq. LifeModel) then    

                            Me%ExternalVar%Mass (POC, Index) = Me%ExternalVar%Mass (POC, Index)             + &
                                                              ( PseudoFaecesByCohort       + &
                                                                FaecesByCohort      )    * &
                                                             C_AtomicMass                                   * &
                                                             Me%ExternalVar%Mass(Number, Index) * Me%DTDay

                        end if

                    end if !PON, POP and POC are not in the list of properties        

                    Cohort => Cohort%Next
                end do
                
            end if !.not.(Particles%Larvae)

            Particles => Particles%Next
        end do 
        
    
    end subroutine UpdateMatrixMass
    
    !--------------------------------------------------------------------------

    subroutine ComputeSomaticMaintenance(Species, Cohort)

        !Arguments-------------------------------------------------------------
        type(T_Species),      pointer           :: Species
        type(T_Cohort),       pointer           :: Cohort

        !Local-----------------------------------------------------------------
        real                                    :: pM, mu_E
        real                                    :: TempCorrection, Vol  
        !Begin-----------------------------------------------------------------


        pM               = Species%IndividualParameters%pM   
        mu_E             = Species%IndividualParameters%mu_E

        TempCorrection   = Species%AuxiliarParameters%TempCorrection         
        Vol              = Cohort%BivalveCondition%Vol

        !J_ES, molC (reserves)/d, Somatic maintenance, Fraccion proporcional to the body volume     
        Cohort%Processes%SomaticMaintenance = pM / mu_E * Tempcorrection * Vol

    end subroutine ComputeSomaticMaintenance

    !--------------------------------------------------------------------------

    subroutine ComputeMobilization(Species, Cohort)

        !Arguments-------------------------------------------------------------
        type(T_Species),      pointer       :: Species
        type(T_Cohort),       pointer       :: Cohort


        !Local-----------------------------------------------------------------
        real                                :: Vol,TempCorrection,E
        real                                :: pM, mu_E, EG,v_cond,kappa  
        !Begin-----------------------------------------------------------------

        Vol            = Cohort%BivalveCondition%Vol
        E              = Cohort%BivalveCondition%E 

        TempCorrection = Species%AuxiliarParameters%TempCorrection         

        pM             = Species%IndividualParameters%pM   
        mu_E           = Species%IndividualParameters%mu_E
        EG             = Species%IndividualParameters%EG     
        v_cond         = Species%IndividualParameters%v_cond 
        kappa          = Species%IndividualParameters%kappa

        if (Vol .eq. 0) then

            Cohort%Processes%Mobilization = 0.0

        else

            !molC (res)/d, Mobilization
            Cohort%Processes%Mobilization = E / (EG/mu_E + kappa * E)                                *  &
                                            (EG/mu_E * v_cond * TempCorrection * Vol**(2./3.) +         &
                                            Cohort%Processes%SomaticMaintenance)

        end if

    end subroutine ComputeMobilization

    !--------------------------------------------------------------------------

    subroutine ComputeReservesDynamics(Index, Cohort)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)           :: Index
        type(T_Cohort),       pointer :: Cohort


        !Local-----------------------------------------------------------------
        integer         :: M_E

        !Begin-----------------------------------------------------------------

        M_E = Cohort%StateIndex%M_E

        !molC (res)/d, Reserves dynamics, J_E
        Cohort%Processes%ReservesDynamics = (Cohort%Processes%Assimilation%C -         &  
                                             Cohort%Processes%Mobilization)

        !Matrix Mass update
        Me%ExternalVar%Mass(M_E, Index) = Me%ExternalVar%Mass(M_E, Index) +            &
                                          Cohort%Processes%ReservesDynamics * Me%DTDay

    end subroutine ComputeReservesDynamics

    !--------------------------------------------------------------------------

    subroutine ComputeStructureDynamics(Index, Species, Cohort)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                 :: Index
        type(T_Species)            , pointer                :: Species
        type(T_Cohort)             , pointer                :: Cohort


        !Local-----------------------------------------------------------------
        integer                                             :: M_V, M_R
        real                                                :: kappa, y_VE

        !Begin-----------------------------------------------------------------

        M_V     = Cohort%StateIndex%M_V
        M_R     = Cohort%StateIndex%M_R

        kappa   = Species%IndividualParameters%kappa

        y_VE    = Species%AuxiliarParameters%y_VE 

        !Compute allocation fraction to growth and somatic maintenance(kJC), molC(reserves)/d
        Cohort%Processes%ToGrowthAndSomatic = kappa * Cohort%Processes%Mobilization

        !Compute flux to growth (JVG), molC(structure)/d
        Cohort%Processes%ToGrowth = (Cohort%Processes%ToGrowthAndSomatic - Cohort%Processes%SomaticMaintenance) * y_VE

        if (Cohort%Processes%ToGrowth .gt. 0) then

            Cohort%Processes%GametesLoss       = 0.0
            Cohort%Processes%StructureLoss     = 0.0

        else
        !(ToGrowth < 0)

            !molC(res)/d
            Cohort%Processes%SomaticMaintenanceNeeds =  - Cohort%Processes%ToGrowth

            Cohort%Processes%ToGrowth  = 0.0

            if (Me%ExternalVar%Mass(M_R, Index) .gt. 0) then

                if (Me%ExternalVar%Mass(M_R, Index) .gt. (Cohort%Processes%SomaticMaintenanceNeeds * Me%DTDay)) then

                    !from gamets
                    Cohort%Processes%GametesLoss   = Cohort%Processes%SomaticMaintenanceNeeds
                    Cohort%Processes%StructureLoss = 0.0

                else
                    !from gamets, molC(res)/d
                    Cohort%Processes%GametesLoss   = Me%ExternalVar%Mass(M_R, Index) / Me%DTDay

                    !and from structure, molC(res)/d 
                    Cohort%Processes%StructureLoss = (Cohort%Processes%SomaticMaintenanceNeeds - Cohort%Processes%GametesLoss) &
                    / y_VE

                end if

            else
            !(MRReproduction = 0)

                Cohort%Processes%GametesLoss   = 0.0

                !and from structure
                Cohort%Processes%StructureLoss = Cohort%Processes%SomaticMaintenanceNeeds / y_VE

            end if

        end if

        !Structure Dynamics, molC(struc)/d
        Cohort%Processes%StructureDynamics = Cohort%Processes%ToGrowth - Cohort%Processes%StructureLoss   

        !Matrix Mass update
        Me%ExternalVar%Mass(M_V,Index) = Me%ExternalVar%Mass(M_V,Index) +             &
        Cohort%Processes%StructureDynamics * Me%DTDay

        if (Me%ExternalVar%Mass(M_V,Index) .lt. Species%IndividualParameters%MVb) then

            Cohort%Processes%DeathByStarvation = Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) /  &
            Me%DTDay !All die from starvation

            if ((Cohort%Dead .eq. 0 )) then 
            
                call ImposeCohortDeath (Index, Species, Cohort) !sets all proc to zero, convert mass to OM, Deadlist

                if ((Species%nCohorts .eq. 1) .and. ((Cohort%Dead .eq. 1 ))) then !this was the last cohort of the population...                        

                    Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                    
                end if
                
                !set the state of the dead cohort to zero
                !Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      = 0.0
                !Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) = 0.0

            end if

        end if    

    end subroutine ComputeStructureDynamics

    !--------------------------------------------------------------------------

    subroutine ComputeMaturity(Index, Species, Cohort)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        type(T_Species)         , pointer           :: Species
        type(T_Cohort)          , pointer           :: Cohort


        !Local-----------------------------------------------------------------
        integer                                     :: L, M_V, M_H, M_R
        real                                        :: kappa 
        real                                        :: kJ, Vol, MHp, MHb, Lb, TempCorrection

        !Begin-----------------------------------------------------------------

        L       = Cohort%StateIndex%L
        M_V     = Cohort%StateIndex%M_V
        M_H     = Cohort%StateIndex%M_H
        M_R     = Cohort%StateIndex%M_R

        kappa   = Species%IndividualParameters%kappa  
        Lb      = Species%IndividualParameters%Lb  

        MHp     = Species%AuxiliarParameters%MHp  
        MHb     = Species%AuxiliarParameters%MHb  
        kJ      = Species%AuxiliarParameters%kJ
        TempCorrection = Species%AuxiliarParameters%TempCorrection         


        Vol     = Cohort%BivalveCondition%Vol


        !Compute allocation fraction to maturity and reproduction (1-k), molC(reserves)/d
        Cohort%Processes%ToMaturityAndReproduction = (1 - kappa) * Cohort%Processes%Mobilization 

        !Compute maturity maintenance, molC(reserves)/d, first estimation
        Cohort%Processes%MaturityMaintenance = kJ * Tempcorrection * Me%ExternalVar%Mass(M_H,Index)

        !Compute flux to reproductionOrmaturity, JER, molC(reserves)/d
        Cohort%Processes%FluxToMatORRepr = Cohort%Processes%ToMaturityAndReproduction - Cohort%Processes%MaturityMaintenance

        if (Cohort%Processes%FluxToMatORRepr .le. 0.0) then

            Cohort%Processes%MaturityMaintenance = Cohort%Processes%ToMaturityAndReproduction

            !The organism losses maturity, molCreserves/d
            Cohort%Processes%MaturityLoss = (Me%ExternalVar%Mass(M_H,Index) - Cohort%Processes%MaturityMaintenance/    &
            (kJ* Tempcorrection)) / Me%DTDay

            !New JER
            Cohort%Processes%FluxToMatORRepr = 0.0

            Cohort%Processes%FluxToGametes   = 0.0

            Cohort%Processes%FluxToMaturity  = 0.0

        else

            Cohort%Processes%MaturityLoss = 0.0

            if (Me%ExternalVar%Mass(M_H,Index) .ge. MHp) then

                Cohort%Processes%FluxToMaturity = 0.0

                Cohort%Processes%FluxToGametes  = Cohort%Processes%FluxToMatORRepr

            else

                Cohort%Processes%FluxToMaturity = Cohort%Processes%FluxToMatORRepr

                Cohort%Processes%FluxToGametes  = 0.0

            end if

        end if 

        !Maturity Dynamics, molC(structure)/d
        Cohort%Processes%MaturityDynamics = Cohort%Processes%FluxToMaturity - Cohort%Processes%MaturityLoss

        !Matrix Mass update
        Me%ExternalVar%Mass(M_H,Index) = Me%ExternalVar%Mass(M_H,Index) +             &
        Cohort%Processes%MaturityDynamics * Me%DTDay

        if ((Me%ExternalVar%Mass(M_H,Index) .lt. MHb) .and. (Me%ExternalVar%Mass(L,Index) .gt. Lb)) then

            Cohort%Processes%DeathByStarvation = Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) /   &
            Me%DTDay !All will die from starvation

            if (Cohort%Dead .eq. 0) then 
                
                call ImposeCohortDeath (Index, Species, Cohort) !sets all proc to zero, convert mass to OM, Deadlist

                if ((Species%nCohorts .eq. 1) .and. (Cohort%Dead .eq. 1 )) then !this was the last cohort of the population...                        

                    Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                    
                end if
                
                !set the state of the dead cohort to zero
                !Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      = 0.0
                !Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) = 0.0

            end if

        end if

    end subroutine ComputeMaturity

    !--------------------------------------------------------------------------

    subroutine ComputeLengthAgeDynamics(Index, Species, Cohort)

        !Arguments-------------------------------------------------------------
        type(T_Species)         , pointer   :: Species
        type(T_Cohort)          , pointer   :: Cohort
        integer, intent(IN)                 :: Index

        !Local-----------------------------------------------------------------
        integer                             :: L, Age, M_V
        real                                :: L_aux, delta_M, LifeSpan, Mv

        !Begin-----------------------------------------------------------------

        L         = Cohort%StateIndex%L
        Age       = Cohort%StateIndex%Age
        M_V       = Cohort%StateIndex%M_V

        delta_M   = Species%IndividualParameters%delta_M
        LifeSpan  = Species%IndividualParameters%LifeSpan

        Mv        = Species%AuxiliarParameters%Mv

        !compute the potential new length value       
        L_aux = (Me%ExternalVar%Mass(M_V,Index)/Mv)**(1.0/3.0) / delta_M 

        Me%ExternalVar%Mass(L,Index) = max(Me%ExternalVar%Mass(L,Index),L_aux)
        Me%ExternalVar%Mass(Age,Index) = Me%ExternalVar%Mass(Age,Index) + Me%DTDay

        if (Me%ExternalVar%Mass(Age,Index) .gt. (LifeSpan * 365)) then

            Cohort%Processes%DeathByAge = Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) /   &
            Me%DTDay !All the bivalves in this cohort died from age

            if (Cohort%Dead .eq. 0 ) then 
                
                call ImposeCohortDeath (Index, Species, Cohort) !sets all proc to zero, convert mass to OM, Deadlist

                if ((Species%nCohorts .eq. 1) .and. (Cohort%Dead .eq. 1 )) then !this was the last cohort of the population...                        

                    Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                    
                end if
                
                !set the state of the dead cohort to zero
                !Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      = 0.0
                !Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) = 0.0

            end if

        end if

    end subroutine ComputeLengthAgeDynamics

    !--------------------------------------------------------------------------

    subroutine ComputeSpawning (Index, Species, Cohort)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: Index
        type(T_Species),      pointer   :: Species
        type(T_Cohort),       pointer   :: Cohort

        !Local-----------------------------------------------------------------
        integer                         :: M_R, M_V, M_E,Number 
        integer                         :: POC, PON, POP 
        real                            :: T 
        real                            :: kap_R, GSR_MIN, GSR_SPAWN
        real                            :: T_SPAWN,MIN_SPAWN_TIME, MEb, MVb
        real                            :: kJ, GSR 
        !Begin-----------------------------------------------------------------

        POC            = Me%PropIndex%POC
        PON            = Me%PropIndex%PON
        POP            = Me%PropIndex%POP


        M_R            = Cohort%StateIndex%M_R
        M_V            = Cohort%StateIndex%M_V
        M_E            = Cohort%StateIndex%M_E
        Number         = Cohort%StateIndex%Number

        T              = Me%ExternalVar%Temperature(index)

        kap_R          = Species%IndividualParameters%kap_R
        GSR_MIN        = Species%IndividualParameters%GSR_MIN
        GSR_SPAWN      = Species%IndividualParameters%GSR_SPAWN 
        T_SPAWN        = Species%IndividualParameters%T_SPAWN 
        MIN_SPAWN_TIME = Species%IndividualParameters%MIN_SPAWN_TIME 
        
        MEb            = Species%IndividualParameters%MEb
        MVb            = Species%IndividualParameters%MVb

        kJ             = Species%AuxiliarParameters%kJ 
        GSR            = Cohort%BivalveCondition%GSR 

        !Flux of gametes in the spawning event, molCreserves/d
        if ((GSR .gt. GSR_SPAWN) .and. (T .ge. T_SPAWN)) then 

            !amount of gametes that should remain in the reproduction buffer, molCreserves
            Cohort%Processes%RemainMRReproduction = GSR_MIN/(1-GSR_MIN) *         &
                                        (Me%ExternalVar%Mass(M_V,Index) + Me%ExternalVar%Mass(M_E,Index))

            !amount of gametes released in the spawning event, molCreserves/d
            Cohort%Processes%Spawning = kap_R *      &
                                        (Me%ExternalVar%Mass(M_R,Index) - Cohort%Processes%RemainMRReproduction) &
                                        / Me%DTDay

            !losses in gametes production, just before a spawning event, molCreserves/d
            Cohort%Processes%SpawningOverhead = (1-kap_R) *        &
                                        (Me%ExternalVar%Mass(M_R,Index) - Cohort%Processes%RemainMRReproduction)  & 
                                        / Me%DTDay

            !number of gametes to be released/d.ind
            Cohort%Processes%GametesToRelease = Cohort%Processes%Spawning/(MEb + MVb)
            
            if (Me%OutputThisIndex) then
                       
                Species%PopulationProcesses%nSpawning = Species%PopulationProcesses%nSpawning + 1

            end if
            
            !number of new born/d.ind
            Cohort%Processes%NewbornsThisCohort = Cohort%Processes%GametesToRelease * & 
                                                  (1 - Species%IndividualParameters%m_spat)  

            !number of new borns that will not survive/d.m3
            Cohort%Processes%NONewbornsThisCohort = Cohort%Processes%GametesToRelease * & 
                                                    Species%IndividualParameters%m_spat 

            !ind/m3
            Species%PopulationProcesses%nNewborns = Species%PopulationProcesses%nNewborns + &
                                                    Cohort%Processes%NewbornsThisCohort * Me%DTDay  * &
                                                    Me%ExternalVar%Mass(Number, Index)
            

            !update mass, gametes that dont survive are converted into POM g/m3
            Me%ExternalVar%Mass(PON,Index) = Me%ExternalVar%Mass(PON,Index)                              + &
                                        Cohort%Processes%NONewbornsThisCohort                            * &
                                        (MEb * Species%SpeciesComposition%ReservesComposition%nN         + &
                                         MVb * Species%SpeciesComposition%StructureComposition%nN)       * &
                                        Me%ExternalVar%Mass(Number, Index) * Me%DTDay                    * &
                                        Species%AuxiliarParameters%N_AtomicMass

            if(Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then

                if (Me%ComputeOptions%Phosphorus) then

                Me%ExternalVar%Mass(POP,Index) = Me%ExternalVar%Mass(POP,Index)                          + &
                                        Cohort%Processes%NONewbornsThisCohort                            * &
                                        (MEb * Species%SpeciesComposition%ReservesComposition%nP         + &
                                         MVb * Species%SpeciesComposition%StructureComposition%nP)       * &
                                        Me%ExternalVar%Mass(Number, Index) * Me%DTDay                    * &
                                        Species%AuxiliarParameters%P_AtomicMass

                end if

            else !(if life)

                Me%ExternalVar%Mass(POC,Index) = Me%ExternalVar%Mass(POC,Index)                          + &
                                        Cohort%Processes%NONewbornsThisCohort * (MEb + MVb)              * &
                                        Me%ExternalVar%Mass(Number, Index) * Me%DTDay                    * &
                                        Species%AuxiliarParameters%C_AtomicMass
            end if !pelagic model


        else

            Cohort%Processes%Spawning             = 0.0
            Cohort%Processes%SpawningOverhead     = 0.0
            Cohort%Processes%GametesToRelease     = 0.0     
            Cohort%Processes%NewbornsThisCohort   = 0.0     
            Cohort%Processes%NONewbornsThisCohort = 0.0     

        end if
        
        

    end subroutine ComputeSpawning

    !--------------------------------------------------------------------------

    subroutine ComputeReproductionDynamics (Index, Cohort)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: Index
        type(T_Cohort)             , pointer            :: Cohort

        !Local-----------------------------------------------------------------
        integer                                         :: M_R

        !Begin-----------------------------------------------------------------

        M_R      = Cohort%StateIndex%M_R

        !Reproduction Buffer Dynamics, molC(structure)/d
        Cohort%Processes%ReproductionDynamics = Cohort%Processes%FluxToGametes              - &
                                                Cohort%Processes%GametesLoss                - &
                                                Cohort%Processes%Spawning                   - &
                                                Cohort%Processes%SpawningOverhead

        !Matrix Mass update
        Me%ExternalVar%Mass(M_R,Index) = Me%ExternalVar%Mass(M_R,Index)                     + &
                                         Cohort%Processes%ReproductionDynamics * Me%DTDay

    end subroutine ComputeReproductionDynamics

    !--------------------------------------------------------------------------

    subroutine ImposeCohortDeath(Index, Species, Cohort)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: Index
        type(T_Species)            , pointer            :: Species
        type(T_Cohort   )          , pointer            :: Cohort

        !Local-----------------------------------------------------------------
        integer                                         :: L, M_V, M_E, M_H, M_R    
        integer                                         :: PON, POP, POC    
        integer                                         :: Age, Number
        !Begin-----------------------------------------------------------------

        POC     = Me%PropIndex%POC
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP

        L       = Cohort%StateIndex%L
        M_V     = Cohort%StateIndex%M_V
        M_E     = Cohort%StateIndex%M_E
        M_H     = Cohort%StateIndex%M_H
        M_R     = Cohort%StateIndex%M_R
        Age     = Cohort%StateIndex%Age
        Number  = Cohort%StateIndex%Number

        !the cohort is dead
        Cohort%Dead = 1
        
        !bivalve biomass and what it had assimilated is converted into POM
        Me%ExternalVar%Mass(PON,Index) = Me%ExternalVar%Mass(PON,Index)                                            + &
                                        ( Me%ExternalVar%Mass(M_V,Index)                                           * &
                                        Species%SpeciesComposition%StructureComposition%nN                         + &
                                        (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index))          * &
                                        Species%SpeciesComposition%ReservesComposition%nN                          + &
                                        Cohort%Processes%Assimilation%N * Me%DTDay )                               * &
                                        Species%AuxiliarParameters%N_AtomicMass                                    * &
                                        Me%ExternalVar%Mass(Number,Index)          

        if(Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then

            if (Me%ComputeOptions%Phosphorus) then

                !bivalve biomass and what it had assimilated is converted into POM
                Me%ExternalVar%Mass(POP,Index) = Me%ExternalVar%Mass(POP,Index)                                     + &
                                                ( Me%ExternalVar%Mass(M_V,Index)                                    * &
                                                Species%SpeciesComposition%StructureComposition%nP                  + &
                                                (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index))   * &
                                                Species%SpeciesComposition%ReservesComposition%nP                   + &
                                                Cohort%Processes%Assimilation%P * Me%DTDay )                        * &
                                                Species%AuxiliarParameters%P_AtomicMass                             * &
                                                Me%ExternalVar%Mass(Number,Index)          

            end if

        else !(if life)

            !bivalve biomass and what it had assimilated is converted into POM
            Me%ExternalVar%Mass(POC,Index) = Me%ExternalVar%Mass(POC,Index)                     + &
                                            ( Me%ExternalVar%Mass(M_V,Index)                    + &
                                            Me%ExternalVar%Mass(M_E,Index)                      + &
                                            Me%ExternalVar%Mass(M_R,Index)                      + &
                                            Cohort%Processes%Assimilation%C  * Me%DTDay )       * &
                                            Species%AuxiliarParameters%C_AtomicMass             * &
                                            Me%ExternalVar%Mass(Number,Index)          

        end if !pelagic model


        
        !Me%ExternalVar%Mass(L,Index)                   = 0.0
        !Me%ExternalVar%Mass(M_V,Index)                 = 0.0
        !Me%ExternalVar%Mass(M_E,Index)                 = 0.0
        !Me%ExternalVar%Mass(M_H,Index)                 = 0.0
        !Me%ExternalVar%Mass(M_R,Index)                 = 0.0
        !Me%ExternalVar%Mass(Age,Index)                 = 0.0
        Me%ExternalVar%Mass(Number,Index)              = 0.0        

        !Make all the processes zero       
        Cohort%Processes%ClearanceRate                 = 0.0
        Cohort%Processes%FilteredInorganic             = 0.0
        Cohort%Processes%FilteredFood%C                = 0.0
        Cohort%Processes%FilteredFood%H                = 0.0
        Cohort%Processes%FilteredFood%O                = 0.0
        Cohort%Processes%FilteredFood%N                = 0.0
        Cohort%Processes%FilteredFood%P                = 0.0

        Cohort%Processes%IngestionInorganic            = 0.0
        Cohort%Processes%IngestionFood%C               = 0.0
        Cohort%Processes%IngestionFood%H               = 0.0
        Cohort%Processes%IngestionFood%O               = 0.0
        Cohort%Processes%IngestionFood%N               = 0.0
        Cohort%Processes%IngestionFood%P               = 0.0

        Cohort%Processes%PFContributionInorganic       = 0.0
        Cohort%Processes%PFContributionFood%C          = 0.0
        Cohort%Processes%PFContributionFood%H          = 0.0
        Cohort%Processes%PFContributionFood%O          = 0.0
        Cohort%Processes%PFContributionFood%N          = 0.0
        Cohort%Processes%PFContributionFood%P          = 0.0

        Cohort%Processes%Assimilation%C                = 0.0
        Cohort%Processes%Assimilation%H                = 0.0
        Cohort%Processes%Assimilation%O                = 0.0
        Cohort%Processes%Assimilation%N                = 0.0
        Cohort%Processes%Assimilation%P                = 0.0

        Cohort%Processes%FaecesContributionInorganic   = 0.0
        Cohort%Processes%FaecesContributionFood%C      = 0.0
        Cohort%Processes%FaecesContributionFood%H      = 0.0
        Cohort%Processes%FaecesContributionFood%O      = 0.0
        Cohort%Processes%FaecesContributionFood%N      = 0.0
        Cohort%Processes%FaecesContributionFood%P      = 0.0

        Cohort%Processes%SomaticMaintenance            = 0.0
        Cohort%Processes%Mobilization                  = 0.0
        Cohort%Processes%ReservesDynamics              = 0.0
        Cohort%Processes%ToGrowthAndSomatic            = 0.0
        Cohort%Processes%ToGrowth                      = 0.0
        Cohort%Processes%GametesLoss                   = 0.0
        Cohort%Processes%StructureLoss                 = 0.0
        Cohort%Processes%SomaticMaintenanceNeeds       = 0.0
        Cohort%Processes%StructureDynamics             = 0.0
        Cohort%Processes%ToMaturityAndReproduction     = 0.0
        Cohort%Processes%MaturityMaintenance           = 0.0
        Cohort%Processes%FluxToMatORRepr               = 0.0

        Cohort%Processes%MaturityLoss                  = 0.0
        Cohort%Processes%FluxToGametes                 = 0.0
        Cohort%Processes%FluxToMaturity                = 0.0
        Cohort%Processes%MaturityDynamics              = 0.0
        Cohort%Processes%RemainMRReproduction          = 0.0
        Cohort%Processes%Spawning                      = 0.0
        Cohort%Processes%SpawningOverhead              = 0.0
        Cohort%Processes%GametesToRelease              = 0.0
        Cohort%Processes%ReproductionDynamics          = 0.0

    end subroutine ImposeCohortDeath

    !--------------------------------------------------------------------------

    subroutine ComputeInorganicFluxes (Index, Species, Cohort)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: Index
        type(T_Species),      pointer   :: Species
        type(T_Cohort),       pointer   :: Cohort


        !Local-----------------------------------------------------------------
        integer                         :: AM, IP, CarbonDioxide, Oxygen 
        integer                         :: Number 
        real                            :: nC_CO2, nH_CO2, nO_CO2,nN_CO2, nP_CO2 
        real                            :: nC_H2O,nH_H2O,nO_H2O, nN_H2O,nP_H2O
        real                            :: nC_O2,nH_O2,nO_O2,nN_O2,nP_O2
        real                            :: nC_NH3,nH_NH3,nO_NH3,nN_NH3,nP_NH3
        real                            :: nC_PO4,nH_PO4,nO_PO4,nN_PO4,nP_PO4
        real                            :: nC_Stru,nH_Stru,nO_Stru,nN_Stru,nP_Stru
        real                            :: nC_Rese, nH_Rese, nO_Rese, nN_Rese,nP_Rese
        !Begin-----------------------------------------------------------------


        AM             = Me%PropIndex%AM
        IP             = Me%PropIndex%IP
        CarbonDioxide  = Me%PropIndex%CarbonDioxide
        Oxygen         = Me%PropIndex%Oxygen
        Number         = Cohort%StateIndex%Number

        nC_CO2         = 1 
        nC_H2O         = 0
        nC_O2          = 0
        nC_NH3         = 0
        nC_PO4         = 0
        nC_Stru        = 1
        nC_Rese        = 1

        nH_CO2         = 0 
        nH_H2O         = 2
        nH_O2          = 0
        nH_NH3         = 3
        nH_PO4         = 0
        nH_Stru        = Species%SpeciesComposition%StructureComposition%nH
        nH_Rese        = Species%SpeciesComposition%ReservesComposition%nH

        nO_CO2         = 2 
        nO_H2O         = 1
        nO_O2          = 2
        nO_NH3         = 0
        nO_PO4         = 4
        nO_Stru        = Species%SpeciesComposition%StructureComposition%nO
        nO_Rese        = Species%SpeciesComposition%ReservesComposition%nO

        nN_CO2         = 0 
        nN_H2O         = 0
        nN_O2          = 0
        nN_NH3         = 1
        nN_PO4         = 0
        nN_Stru        = Species%SpeciesComposition%StructureComposition%nN
        nN_Rese        = Species%SpeciesComposition%ReservesComposition%nN

        nP_CO2         = 0 
        nP_H2O         = 0
        nP_O2          = 0
        nP_NH3         = 0
        nP_PO4         = 1
        nP_Stru        = Species%SpeciesComposition%StructureComposition%nP
        nP_Rese        = Species%SpeciesComposition%ReservesComposition%nP

        !Compute Ammonia Fluxes, molN/d
        Cohort%Processes%InorganicFluxes%NH3 = -1/nN_NH3 * ( -Cohort%Processes%FilteredFood%N                 + &
                                                            Cohort%Processes%PFContributionFood%N             + &
                                                            Cohort%Processes%FaecesContributionFood%N         + &
                                                            Cohort%Processes%StructureDynamics     * nN_Stru  + &
                                                            Cohort%Processes%ReservesDynamics      * nN_Rese  + &
                                                            Cohort%Processes%ReproductionDynamics  * nN_Rese  + &
                                                            Cohort%Processes%GametesToRelease                 * &
                                                            (Species%IndividualParameters%MEb      * nN_Stru  + &
                                                             Species%IndividualParameters%MVb      * nN_Rese))

        !Compute Water Fluxes, molH2O/d
        Cohort%Processes%InorganicFluxes%H2O = -1/nH_H2O * ( -Cohort%Processes%FilteredFood%H                 + &
                                                            Cohort%Processes%InorganicFluxes%NH3   * nH_NH3   + &
                                                            Cohort%Processes%PFContributionFood%H             + &
                                                            Cohort%Processes%FaecesContributionFood%H         + &
                                                            Cohort%Processes%StructureDynamics     * nH_Stru  + &
                                                            Cohort%Processes%ReservesDynamics      * nH_Rese  + &
                                                            Cohort%Processes%ReproductionDynamics  * nH_Rese  + &
                                                            Cohort%Processes%GametesToRelease                 * &
                                                            (Species%IndividualParameters%MEb      * nH_Stru  + &
                                                             Species%IndividualParameters%MVb      * nH_Rese))

        !Compute CO2, molCO2/d
        Cohort%Processes%InorganicFluxes%CO2 = -1/nC_CO2 * ( -Cohort%Processes%FilteredFood%C                 + &
                                                            Cohort%Processes%PFContributionFood%C             + &
                                                            Cohort%Processes%FaecesContributionFood%C         + &
                                                            Cohort%Processes%StructureDynamics     * nC_Stru  + &
                                                            Cohort%Processes%ReservesDynamics      * nC_Rese  + &
                                                            Cohort%Processes%ReproductionDynamics  * nC_Rese  + &
                                                            Cohort%Processes%GametesToRelease                 * &
                                                            (Species%IndividualParameters%MEb      * nC_Stru  + &
                                                             Species%IndividualParameters%MVb      * nC_Rese))

        if (Me%ComputeOptions%Phosphorus) then

            !Compute Phosphate Fluxes, molP/d
            Cohort%Processes%InorganicFluxes%PO4 = -1/nP_PO4 * ( -Cohort%Processes%FilteredFood%P                 + &
                                                                Cohort%Processes%PFContributionFood%P             + &
                                                                Cohort%Processes%FaecesContributionFood%P         + &
                                                                Cohort%Processes%StructureDynamics     * nP_Stru  + &
                                                                Cohort%Processes%ReservesDynamics      * nP_Rese  + &
                                                                Cohort%Processes%ReproductionDynamics  * nP_Rese  + &
                                                                Cohort%Processes%GametesToRelease                 * &
                                                                (Species%IndividualParameters%MEb      * nP_Stru  + &
                                                                Species%IndividualParameters%MVb       * nP_Rese))

            !Compute Oxygen Fluxes, molO2/d
            Cohort%Processes%InorganicFluxes%O2  = -1/nO_O2 * ( -Cohort%Processes%FilteredFood%O                     + &
                                                                Cohort%Processes%InorganicFluxes%PO4      * nO_PO4   + &
                                                                Cohort%Processes%InorganicFluxes%CO2      * nO_H2O   + &
                                                                Cohort%Processes%InorganicFluxes%H2O      * nO_CO2   + &
                                                                Cohort%Processes%PFContributionFood%O                + &
                                                                Cohort%Processes%FaecesContributionFood%O            + &
                                                                Cohort%Processes%StructureDynamics        * nO_Stru  + &
                                                                Cohort%Processes%ReservesDynamics         * nO_Rese  + &
                                                                Cohort%Processes%ReproductionDynamics     * nO_Rese  + &
                                                                Cohort%Processes%GametesToRelease                    * &
                                                                (Species%IndividualParameters%MEb         * nO_Stru  + &
                                                                Species%IndividualParameters%MVb          * nO_Rese))

        else !no P

            !Compute Oxygen Fluxes, molO2/d
            Cohort%Processes%InorganicFluxes%O2 = -1/nO_O2 * ( -Cohort%Processes%FilteredFood%O                      + &
                                                                Cohort%Processes%InorganicFluxes%CO2      * nO_H2O   + &
                                                                Cohort%Processes%InorganicFluxes%H2O      * nO_CO2   + &
                                                                Cohort%Processes%PFContributionFood%O                + &
                                                                Cohort%Processes%FaecesContributionFood%O            + &
                                                                Cohort%Processes%StructureDynamics         * nO_Stru + &
                                                                Cohort%Processes%ReservesDynamics          * nO_Rese + &
                                                                Cohort%Processes%ReproductionDynamics      * nO_Rese + &
                                                                Cohort%Processes%GametesToRelease                    * &
                                                                (Species%IndividualParameters%MEb          * nO_Stru + &
                                                                Species%IndividualParameters%MVb           * nO_Rese))
            
        end if


        !update mass, g/m3    
        Me%ExternalVar%Mass(AM,Index) = Me%ExternalVar%Mass(AM,Index)                      + &
                                        (Cohort%Processes%InorganicFluxes%NH3              * &
                                        Species%AuxiliarParameters%N_AtomicMass)           * &
                                        Me%ExternalVar%Mass(Number, Index) * Me%DTDay  


        Me%ExternalVar%Mass(CarbonDioxide,Index) = Me%ExternalVar%Mass(CarbonDioxide,Index)   + &
                                                (Cohort%Processes%InorganicFluxes%CO2         * &
                                                Species%AuxiliarParameters%C_AtomicMass)      * &
                                                Me%ExternalVar%Mass(Number, Index) * Me%DTDay  

        Me%ExternalVar%Mass(Oxygen,Index) = Me%ExternalVar%Mass(Oxygen,Index)                 + &
                                            (Cohort%Processes%InorganicFluxes%O2              * &
                                            Species%AuxiliarParameters%O_AtomicMass)          * &
                                            Me%ExternalVar%Mass(Number, Index) * Me%DTDay  

        if (Me%ComputeOptions%Phosphorus) then

            !Mass Uptade
            Me%ExternalVar%Mass(IP,Index) = Me%ExternalVar%Mass(IP,Index)                     + &
                                            (Cohort%Processes%InorganicFluxes%PO4             * &
                                            Species%AuxiliarParameters%P_AtomicMass)          * &
                                            Me%ExternalVar%Mass(Number, Index) * Me%DTDay  
        end if

    end subroutine ComputeInorganicFluxes

    !--------------------------------------------------------------------------

    subroutine ComputeExtraStarvationMortality (Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        type(T_Species)             , pointer       :: Species
        type(T_Cohort)              , pointer       :: Cohort
        integer                                     :: Number, M_V, M_E, M_R   
        integer                                     :: PON, POC, POP
        real                                        :: StarvationMortality   

        !Begin-----------------------------------------------------------------

        POC     = Me%PropIndex%POC
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP

        Species => Me%FirstSpecies
d1:     do while(associated(Species))

            if (Species%ExtraStarvation) then
            
                 Cohort => Species%FirstCohort
                 do while(associated(Cohort))

                        if ((Cohort%Dead .eq. 0 ) .and. (Cohort%BivalveCondition%ScaledE .le. 0.1))then

                            M_V     = Cohort%StateIndex%M_V
                            M_E     = Cohort%StateIndex%M_E
                            M_R     = Cohort%StateIndex%M_R
                            Number  = Cohort%StateIndex%Number
                        
                            StarvationMortality = 1.-(1.+100.*exp(-70.*Cohort%BivalveCondition%ScaledE))**(-1.)

                            !Extra mortality depending on the bivalve condition, #/d.m3
                            Cohort%Processes%DeathByExtraStarvation = Me%ExternalVar%Mass(Number,Index)     * &
                                                                      StarvationMortality
                                                              

                            !update the number of organisms in mass matrix
                            Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index)           - &
                                                                (Cohort%Processes%DeathByExtraStarvation * Me%DTDay) 


                            !update POM again
                            Me%ExternalVar%Mass(PON,Index) = Me%ExternalVar%Mass(PON,Index)                                   + &
                                                            ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                            Species%SpeciesComposition%StructureComposition%nN                + &
                                                            (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                            Species%SpeciesComposition%ReservesComposition%nN )               * &
                                                            Species%AuxiliarParameters%N_AtomicMass                           * &
                                                            Cohort%Processes%DeathByExtraStarvation * Me%DTDay           

                            if(Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then


                                if (Me%ComputeOptions%Phosphorus) then

                                    Me%ExternalVar%Mass(POP,Index) = Me%ExternalVar%Mass(POP,Index)                           + &
                                                            ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                            Species%SpeciesComposition%StructureComposition%nP                + &
                                                            (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                            Species%SpeciesComposition%ReservesComposition%nP )               * &
                                                            Species%AuxiliarParameters%P_AtomicMass                           * &
                                                            Cohort%Processes%DeathByExtraStarvation * Me%DTDay           

                                end if

                            else !(if life)


                                Me%ExternalVar%Mass(POC,Index) = Me%ExternalVar%Mass(POC,Index)              + &
                                                                ( Me%ExternalVar%Mass(M_V,Index)             + &
                                                                Me%ExternalVar%Mass(M_E,Index)               + &
                                                                Me%ExternalVar%Mass(M_R,Index) )             * &
                                                                Cohort%Processes%DeathByExtraStarvation * Me%DTDay

                            end if !pelagic model

                        end if !its alive             

                        Cohort => Cohort%Next
                    end do
            end if 

            Species => Species%Next
        end do d1

    end subroutine ComputeExtraStarvationMortality

    !--------------------------------------------------------------------------

    subroutine ComputeNaturalMortality (Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        type(T_Species)             , pointer       :: Species
        type(T_Cohort)              , pointer       :: Cohort
        integer                                     :: Number, M_V, M_E, M_R   
        integer                                     :: PON, POC, POP   

        !Begin-----------------------------------------------------------------


        Species => Me%FirstSpecies
d1:     do while(associated(Species))

            Cohort => Species%FirstCohort
d2:         do while(associated(Cohort))

                if (Cohort%Dead .eq. 0 ) then
                
                    POC     = Me%PropIndex%POC
                    PON     = Me%PropIndex%PON
                    POP     = Me%PropIndex%POP
                
                    M_V     = Cohort%StateIndex%M_V
                    M_E     = Cohort%StateIndex%M_E
                    M_R     = Cohort%StateIndex%M_R
                    Number  = Cohort%StateIndex%Number

                    !Natural mortality, #/d.m3
                    Cohort%Processes%DeathByNatural = Me%ExternalVar%Mass(Number,Index)             * &
                                                      Species%IndividualParameters%m_natural
                                                      

                    !update the number of organisms in mass matrix
                    Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index)           - &
                                                        (Cohort%Processes%DeathByNatural * Me%DTDay) 


                    !update POM again
                    Me%ExternalVar%Mass(PON,Index) = Me%ExternalVar%Mass(PON,Index)                                   + &
                                                    ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                    Species%SpeciesComposition%StructureComposition%nN                + &
                                                    (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                    Species%SpeciesComposition%ReservesComposition%nN )               * &
                                                    Species%AuxiliarParameters%N_AtomicMass                           * &
                                                    Cohort%Processes%DeathByNatural * Me%DTDay           

                    if(Me%ComputeOptions%PelagicModel .eq. WaterQualityModel) then


                        if (Me%ComputeOptions%Phosphorus) then

                            Me%ExternalVar%Mass(POP,Index) = Me%ExternalVar%Mass(POP,Index)                           + &
                                                    ( Me%ExternalVar%Mass(M_V,Index)                                  * &
                                                    Species%SpeciesComposition%StructureComposition%nP                + &
                                                    (Me%ExternalVar%Mass(M_E,Index) + Me%ExternalVar%Mass(M_R,Index)) * &
                                                    Species%SpeciesComposition%ReservesComposition%nP )               * &
                                                    Species%AuxiliarParameters%P_AtomicMass                           * &
                                                    Cohort%Processes%DeathByNatural * Me%DTDay           

                        end if

                    else !(if life)


                        Me%ExternalVar%Mass(POC,Index) = Me%ExternalVar%Mass(POC,Index)              + &
                                                        ( Me%ExternalVar%Mass(M_V,Index)             + &
                                                        Me%ExternalVar%Mass(M_E,Index)               + &
                                                        Me%ExternalVar%Mass(M_R,Index) )             * &
                                                        Cohort%Processes%DeathByNatural * Me%DTDay

                    end if !pelagic model

                end if !its alive             

                Cohort => Cohort%Next
            end do d2

            Species => Species%Next
        end do d1

    end subroutine ComputeNaturalMortality

    !--------------------------------------------------------------------------

    subroutine ComputePredation (Index, CheckIfOpenPoint)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: Index
        integer, intent(IN)                 :: CheckIfOpenPoint


        !Local-----------------------------------------------------------------
        type(T_Species)          , pointer  :: Species
        type(T_Cohort)           , pointer  :: Cohort
        character(StringLength)             :: PredatorName


        !Begin-----------------------------------------------------------------

        call CheckFeedingByPredator (CheckIfOpenPoint)

        PredatorName = 'shrimp'
        call ComputePredationByPredator (Index, PredatorName, Me%PropIndex%Shrimp)

        PredatorName = 'crab'
        call ComputePredationByPredator (Index, PredatorName, Me%PropIndex%Crab)

        PredatorName = 'oystercatcher'
        call ComputePredationByPredator (Index, PredatorName, Me%PropIndex%OysterCatcher)

        PredatorName = 'eider duck'
        call ComputePredationByPredator (Index, PredatorName, Me%PropIndex%EiderDuck)

        PredatorName = 'herring gull'
        call ComputePredationByPredator (Index, PredatorName, Me%PropIndex%HerringGull)

        Species => Me%FirstSpecies
d1:     do while(associated(Species))

            Cohort => Species%FirstCohort
d2:         do while(associated(Cohort))

                    if ((.not. Cohort%Larvae ) .and. & 
                       (Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) .lt. Me%MinNumber)) then

                    if (Cohort%Dead .eq. 0 ) then 
                    
                        call ImposeCohortDeath (Index, Species, Cohort) !sets all proc to zero, convert mass to OM, Deadlist

                        if ((Species%nCohorts .eq. 1) .and. (Cohort%Dead .eq. 1 )) then 
                        !if the last cohort in the population                          

                            Cohort%Processes%DeathByLowNumbers = Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) /   &
                                                   Me%DTDay !All die from 'low numbers'
                            
                            Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                            

                        end if
                        
                        !set the state of the dead cohort to zero
                        !Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)      = 0.0
                        !Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index) = 0.0
                        
                    end if


                
                end if

                Cohort => Cohort%Next
            end do d2

        Species => Species%Next
        end do d1

    end subroutine ComputePredation

    !--------------------------------------------------------------------------

    subroutine CheckFeedingByPredator (CheckIfOpenPoint)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: CheckIfOpenPoint

        !Local-----------------------------------------------------------------
        type(T_Species)  ,      pointer :: Species
        type(T_Predator) ,      pointer :: Predator

        !Begin-----------------------------------------------------------------

        Species => Me%FirstSpecies
d1:     do while(associated(Species))

            Predator => Species%FirstPredator
d2:         do while(associated(Predator))

                Predator%Feeding = .false.

                !always
                if (Predator%Feeding_Time .eq. 1) Predator%Feeding = .true. 

                !low tide
                if ((Predator%Feeding_Time .eq. 2) .and. (CheckIfOpenPoint .ne. OpenPoint)) Predator%Feeding = .true.

                !high tide
                if ((Predator%Feeding_Time .eq. 3) .and. (CheckIfOpenPoint .eq. OpenPoint)) Predator%Feeding = .true.

                Predator => Predator%Next
            end do d2

        Species => Species%Next
        end do d1

    end subroutine CheckFeedingByPredator

    !--------------------------------------------------------------------------

    subroutine ComputePredationByPredator (Index, PredatorName, PredatorIndex)

        !Arguments-------------------------------------------------------------
        integer                , intent(IN)   :: Index
        character(StringLength), intent(IN)   :: PredatorName
        integer                , intent(IN)   :: PredatorIndex

        !Local-----------------------------------------------------------------
        type(T_Species),      pointer         :: Species
        type(T_Cohort),       pointer         :: Cohort
        type(T_Predator),     pointer         :: Predator
        real                                  :: Tref, TA, TL, TH, TAL, TAH, T
        real                                  :: TotalPreyAvailable
        real                                  :: PotentialCohortPredation, ThisCohortPredation
        integer                               :: Number, L

        !Begin-----------------------------------------------------------------

        !Compute the total number of organisms available for this predator        
        TotalPreyAvailable = 0

        Species => Me%FirstSpecies
d1:     do while(associated(Species))

            Predator => Species%FirstPredator
d2:         do while(associated(Predator))

                if ((Predator%ID%Name .eq. PredatorName) .and. Predator%Feeding) then

                    !Compute the number of prey available
                    Cohort => Species%FirstCohort
d3:                 do while(associated(Cohort))

                        if (Cohort%Dead .eq. 0 ) then

                            if ((Me%ExternalVar%Mass(Cohort%StateIndex%L,Index) .ge. Predator%MinPreySize) .and.  &
                                (Me%ExternalVar%Mass(Cohort%StateIndex%L,Index) .le. Predator%MaxPreySize)  ) then

                                !#/m3
                                TotalPreyAvailable = TotalPreyAvailable + Me%ExternalVar%Mass(Cohort%StateIndex%Number,Index)

                            end if
                        
                        end if

                        Cohort => Cohort%Next
                    end do d3

                end if

                Predator => Predator%Next
            end do d2

            Species => Species%Next
        end do d1

        !Compute how much will be predated by each cohort, depending on its relative abundance, in numbers
        Species => Me%FirstSpecies
d4:     do while(associated(Species))

            Predator => Species%FirstPredator
d5:         do while(associated(Predator))

                if ((Predator%ID%Name .eq. PredatorName) .and. Predator%Feeding) then

                    !Temperature Correction factor for the predator
                    Tref  = Predator%P_Tref
                    TA    = Predator%P_TA
                    TL    = Predator%P_TL
                    TH    = Predator%P_TH
                    TAL   = Predator%P_TAL
                    TAH   = Predator%P_TAH

                    !K, Actual Temperature, oC to K
                    T = Me%ExternalVar%Temperature(Index) + 273.

                    if (Predator%CORRECT_TEMP .eq. 0.0) then

                        Predator%TempCorrection  = 1.0

                    else
                        if (Predator%SIMPLE_TEMP .eq. 0.0) then
                        
                            Predator%TempCorrection  = exp(TA/Tref-TA/T)                                   * &
                                                    (1.0 + exp(TAL/Tref - TAL/TL) + exp(TAH/TH-TAH/Tref))  / &
                                                    (1.0 + exp(TAL/T - TAL/TL) + exp(TAH/TH-TAH/T))
                                                
                        else 
                                    
                            Predator%TempCorrection  = exp(TA/Tref-TA/T)  
                        
                        end if
                        
                    end if

                    !Compute the predation on this cohort
                    Cohort => Species%FirstCohort
d6:                 do while(associated(Cohort))

                        if (Cohort%Dead .eq. 0) then

                            Number = Cohort%StateIndex%Number
                            L      = Cohort%StateIndex%L

                            if ((Me%ExternalVar%Mass(L,Index) .ge. Predator%MinPreySize) .and.  &
                                (Me%ExternalVar%Mass(L,Index) .le. Predator%MaxPreySize) .and.  &
                                 Me%ExternalVar%Mass(Number,Index) .gt. 0.0 ) then

                                if (TotalPreyAvailable .eq. 0.0) then
                                    Predator%TotalFeeding_Rate = 0.0
                                    PotentialCohortPredation  = 0.0
                                else

                                    !Compute how much will be predated from this cohort
                                    select case (Predator%Feeding_Units)

                                    case (1) !#/d.ind, crab

                                        !#prey/d.m3 = #prey/d.ind * ind/m3
                                        Predator%TotalFeeding_Rate = Predator%Feeding_Rate                                     * &
                                                                     Me%ExternalVar%Mass(PredatorIndex,Index)

                                        !#prey/d.m3
                                        PotentialCohortPredation = Me%ExternalVar%Mass(Number,Index) / TotalPreyAvailable      * &
                                                                   Predator%TotalFeeding_Rate * Predator%TempCorrection        * &
                                                                   Predator%Diet

                                    case (2) !AFDW/d.ind, birds

                                        !molC/d.m3 = gAFDW/d.ind * gdw/gafdw * gC/gdw * molC/gC * ind/m3  
                                        Predator%TotalFeeding_Rate = Predator%Feeding_Rate/Predator%AfdwToC/Predator%AfdwToC/12 * &
                                                                     Me%ExternalVar%Mass(PredatorIndex,Index)

                                        !molC/d.m3
                                        PotentialCohortPredation = Me%ExternalVar%Mass(Number,Index) / TotalPreyAvailable       * &
                                                                   Predator%TotalFeeding_Rate * Predator%TempCorrection         * &
                                                                   Predator%Diet

                                        !#prey/d.m3
                                        PotentialCohortPredation = PotentialCohortPredation / Cohort%BivalveCondition%TotalmolC 

                                    case (3) !J/cm2.d.ind, shrimps

                                        !J/d.m3 = J/cm2.d.ind * cm2 * ind/m3 
                                        Predator%TotalFeeding_Rate = Predator%Feeding_Rate * Predator%PredatorSize ** 2         * &
                                                                     Me%ExternalVar%Mass(PredatorIndex,Index)

                                        !J/d.m3
                                        PotentialCohortPredation = Me%ExternalVar%Mass(Number,Index) / TotalPreyAvailable       * &
                                                                   Predator%TotalFeeding_Rate * Predator%TempCorrection         * &
                                                                   Predator%Diet

                                        !#prey/d.m3, mu_E = mu_V 
                                        PotentialCohortPredation = PotentialCohortPredation / Species%IndividualParameters%mu_E / &
                                                                   Cohort%BivalveCondition%TotalmolC
                                    
                                    end select
                                end if !TotalPreyAvailable==0

                                !Restrict ThisCohortPredation based on the number of individuals in this cohort, #/m3
                                ThisCohortPredation = min(Me%ExternalVar%Mass(Number,Index),PotentialCohortPredation * Me%DTDay) 

                                !Store how much is predated by each predator
                                if (Predator%ID%Name .eq. 'shrimp') then  
                                    
                                    Cohort%Processes%PredationByShrimps =  ThisCohortPredation/Me%DTday
                                    
                                 Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                                    
                                   
                                end if

                                if (Predator%ID%Name .eq. 'crab') then  
                                    
                                    Cohort%Processes%PredationByCrabs =  ThisCohortPredation/Me%DTday
                                    
                                    Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)

                                end if

                                if (Predator%ID%Name .eq. 'oystercatcher') then  
                                    
                                    Cohort%Processes%PredationByOysterCatchers =  ThisCohortPredation/Me%DTday
                                    
                                    Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                                
                                end if

                                if (Predator%ID%Name .eq. 'eider duck') then  
                                    
                                    Cohort%Processes%PredationByEiderDucks =  ThisCohortPredation/Me%DTday
                                    
                                    Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                                    
                                end if

                                if (Predator%ID%Name .eq. 'herring gull') then 
                                 
                                    Cohort%Processes%PredationByHerringGull =  ThisCohortPredation/Me%DTday
                                    
                                    Species%PopulationProcesses%LastLength          = Me%ExternalVar%Mass(Cohort%StateIndex%L,Index)
                                    
                                end if

                                !Update Mass, #/m3
                                Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index) - ThisCohortPredation    

                            end if !if the cohort have the right size for this predator

                        end if
                        Cohort => Cohort%Next
                    end do d6
                end if !it is this predator

                Predator => Predator%Next
            end do d5

            Species => Species%Next
        end do d4

    end subroutine ComputePredationByPredator


   !--------------------------------------------------------------------------

    subroutine ComputePopulationVariables (Index)

        !Arguments-------------------------------------------------------------
        type(T_Species),               pointer  :: Species
        type(T_Cohort),                pointer  :: Cohort
        integer, intent(IN)                     :: Index

        !Local-----------------------------------------------------------------
        type(T_PopulationProcesses),   pointer  :: PopulationProcesses
        integer                                 :: L, Number
        !Begin-----------------------------------------------------------------

        Species => Me%FirstSpecies
d1:     do while(associated(Species))
            
            PopulationProcesses => Species%PopulationProcesses

            Cohort => Species%FirstCohort
d2:         do while(associated(Cohort))

                    Number = Cohort%StateIndex%Number    
                    L      = Cohort%StateIndex%L
                            
                    !#/m3
                    PopulationProcesses%TN        = PopulationProcesses%TN             + &
                                                  Me%ExternalVar%Mass(Number,Index)
                                                    
                    !#/m3
                    PopulationProcesses%NCoh      = PopulationProcesses%NCoh           + 1
                          
                    !molC/m3
                    PopulationProcesses%TBio      = PopulationProcesses%TBio           + &
                                                  Cohort%BivalveCondition%TotalmolC    * &
                                                  Me%ExternalVar%Mass(Number,Index)
                                                    
                    !m3/d.m3
                    PopulationProcesses%Cr        = PopulationProcesses%Cr             + &
                                                  Cohort%Processes%ClearanceRate       * &
                                                  Me%ExternalVar%Mass(Number,Index)
                                                    
                    !molC/m3
                    PopulationProcesses%Fil       = PopulationProcesses%Fil            + &
                                                  Cohort%Processes%FilteredFood%C      * &
                                                  Me%ExternalVar%Mass(Number,Index) 
                                                       
                    !molC//m3
                    PopulationProcesses%Ing       = PopulationProcesses%Ing            + &
                                                  Cohort%Processes%IngestionFood%C     * &
                                                  Me%ExternalVar%Mass(Number,Index)
                                                        
                    !molC/m3
                    PopulationProcesses%Ass       = PopulationProcesses%Ass            + &
                                                  Cohort%Processes%Assimilation%C      * &
                                                  Me%ExternalVar%Mass(Number,Index) 
                                                       
                    !molCO2/m3
                    PopulationProcesses%CO2       = PopulationProcesses%CO2            + &
                                                  Cohort%Processes%InorganicFluxes%CO2 * &
                                                  Me%ExternalVar%Mass(Number,Index) 
                                                       
                    !molH2O/m3
                    PopulationProcesses%H2O       = PopulationProcesses%H2O            + &
                                                  Cohort%Processes%InorganicFluxes%H2O * &
                                                  Me%ExternalVar%Mass(Number,Index) 
                                                         
                    !O2/m3
                    PopulationProcesses%O         = PopulationProcesses%O              + &
                                                  Cohort%Processes%InorganicFluxes%O2  * &
                                                  Me%ExternalVar%Mass(Number,Index)  
                                                        
                    !NH3/m3
                    PopulationProcesses%NH3       = PopulationProcesses%NH3            + &
                                                  Cohort%Processes%InorganicFluxes%NH3 * &
                                                  Me%ExternalVar%Mass(Number,Index)
                                                          
                    !PO4/m3
                    PopulationProcesses%PO4       = PopulationProcesses%PO4            + &
                                                  Cohort%Processes%InorganicFluxes%PO4 * &
                                                  Me%ExternalVar%Mass(Number,Index) 
                                                             
                    !#/m3
                    PopulationProcesses%m_A       = PopulationProcesses%m_A                            + &
                                                  Cohort%Processes%DeathByAge * Me%DTday
                                                             
                    !#/m3
                    PopulationProcesses%m_O       = PopulationProcesses%m_O                            + &
                                                  Cohort%Processes%DeathByOxygen * Me%DTday
                                                      
                    PopulationProcesses%m_F       = PopulationProcesses%m_F                             + &
                                                  Cohort%Processes%DeathByStarvation * Me%DTday         + &
                                                  Cohort%Processes%DeathByExtraStarvation * Me%DTday 
                                                      
                    PopulationProcesses%m_nat     = PopulationProcesses%m_nat                           + &
                                                  Cohort%Processes%DeathByNatural * Me%DTday
                                                      
                    PopulationProcesses%m_shr     = PopulationProcesses%m_shr                           + &
                                                  Cohort%Processes%PredationByShrimps * Me%DTday
                                                      
                    PopulationProcesses%m_cra     = PopulationProcesses%m_cra                           + &
                                                  Cohort%Processes%PredationByCrabs * Me%DTday
                                                      
                    PopulationProcesses%m_oys     = PopulationProcesses%m_oys                           + &
                                                  Cohort%Processes%PredationByOysterCatchers * Me%DTday
                                                      
                    PopulationProcesses%m_duck    = PopulationProcesses%m_duck                          + &
                                                  Cohort%Processes%PredationByEiderDucks  * Me%DTday 
                                                      
                    PopulationProcesses%m_gull    = PopulationProcesses%m_gull                          + &
                                                  Cohort%Processes%PredationByHerringGull * Me%DTday
                                                      
                    PopulationProcesses%m_low     = PopulationProcesses%m_low                           + &
                                                  Cohort%Processes%DeathByLowNumbers * Me%DTday
                                                                                                   
                    PopulationProcesses%m_self    = PopulationProcesses%m_self                           + &
                                                  Cohort%Processes%DeathBySelfPredation * Me%DTday

                    PopulationProcesses%m_others  = PopulationProcesses%m_others                         + &
                                                  Cohort%Processes%DeathByLarvaePredationByOthers * Me%DTday

                    PopulationProcesses%m_vel    = PopulationProcesses%m_vel                             + &
                                                  Cohort%Processes%DeathByVelocity * Me%DTday

                    PopulationProcesses%m_settle    = PopulationProcesses%m_settle                       + &
                                                  Cohort%Processes%DeathByWrongSettlement * Me%DTday

                    !molC/m3
                    PopulationProcesses%Massm_A   = PopulationProcesses%Massm_A                         + &
                                                  Cohort%Processes%DeathByAge * Me%DTday                * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                              
                    PopulationProcesses%Massm_O   = PopulationProcesses%Massm_O                         + &
                                                  Cohort%Processes%DeathByOxygen * Me%DTday             * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                              
                    PopulationProcesses%Massm_F   = PopulationProcesses%Massm_F                         + &
                                                  (Cohort%Processes%DeathByStarvation                   + &
                                                   Cohort%Processes%DeathByExtraStarvation)* Me%DTday   * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                              
                    PopulationProcesses%Massm_nat = PopulationProcesses%Massm_nat                       + &
                                                  Cohort%Processes%DeathByNatural * Me%DTday            * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                              
                    PopulationProcesses%Massm_shr = PopulationProcesses%Massm_shr                       + &
                                                  Cohort%Processes%PredationByShrimps * Me%DTday        * &
                                                  Cohort%BivalveCondition%TotalmolC 
                                                             
                    PopulationProcesses%Massm_cra = PopulationProcesses%Massm_cra                       + &
                                                  Cohort%Processes%PredationByCrabs * Me%DTday          * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                              
                    PopulationProcesses%Massm_oys = PopulationProcesses%Massm_oys                       + &
                                                  Cohort%Processes%PredationByOysterCatchers * Me%DTday * &
                                                  Cohort%BivalveCondition%TotalmolC 
                                                             
                    PopulationProcesses%Massm_duck= PopulationProcesses%Massm_duck                      + &
                                                  Cohort%Processes%PredationByEiderDucks * Me%DTday     * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                              
                    PopulationProcesses%Massm_gull= PopulationProcesses%Massm_gull                      + &
                                                  Cohort%Processes%PredationByHerringGull * Me%DTday    * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                              
                    PopulationProcesses%Massm_low = PopulationProcesses%Massm_low                       + &
                                                  Cohort%Processes%DeathByLowNumbers * Me%DTday         * &
                                                  Cohort%BivalveCondition%TotalmolC

                    PopulationProcesses%Massm_self = PopulationProcesses%Massm_self                     + &
                                                  Cohort%Processes%DeathBySelfPredation * Me%DTday      * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                  
                    PopulationProcesses%Massm_vel = PopulationProcesses%Massm_vel                       + &
                                                  Cohort%Processes%DeathByVelocity * Me%DTday           * &
                                                  Cohort%BivalveCondition%TotalmolC
                                                  
                    PopulationProcesses%Massm_others = PopulationProcesses%Massm_others                 + &
                                                  Cohort%Processes%DeathByLarvaePredationByOthers * Me%DTday   * &
                                                  Cohort%BivalveCondition%TotalmolC

                    PopulationProcesses%Massm_settle = PopulationProcesses%Massm_settle                 + &
                                                  Cohort%Processes%DeathByWrongSettlement * Me%DTday    * &
                                                  Cohort%BivalveCondition%TotalmolC
                    
                    if (Me%ExternalVar%Mass(L,Index) .ge. 0.1) then
                        !#/m3               
                        PopulationProcesses%TNField = PopulationProcesses%TNField + Me%ExternalVar%Mass(Number,Index)
                    end if
                    
                    PopulationProcesses%MaxLength = max(PopulationProcesses%MaxLength,Me%ExternalVar%Mass(L,Index))

                    if (Me%OutputThisIndex) then
                               
                        Me%MaxTNField = max(Me%MaxTNField,PopulationProcesses%TNField)

                    end if
                
                Cohort => Cohort%Next
            end do d2            

            if (PopulationProcesses%TN .ne. 0.0) then
            
                PopulationProcesses%nInstantsForAverage = PopulationProcesses%nInstantsForAverage + 1
                
                !Product of mortalities in numbers
                if (1.-PopulationProcesses%m_A /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(1)  = PopulationProcesses%SumLogAllMortalityInNumbers(1)   + &
                                                            log10(1.-PopulationProcesses%m_A /PopulationProcesses%TNStartTimeStep) 
                end if
                
                if (1.-PopulationProcesses%m_O /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(2)  = PopulationProcesses%SumLogAllMortalityInNumbers(2)   + &
                                                            log10(1.-PopulationProcesses%m_O/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_F /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(3)  = PopulationProcesses%SumLogAllMortalityInNumbers(3)   + &
                                                            log10(1.-PopulationProcesses%m_F/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_nat /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(4)  = PopulationProcesses%SumLogAllMortalityInNumbers(4)   + &
                                                            log10(1.-PopulationProcesses%m_nat/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_shr /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(5)  = PopulationProcesses%SumLogAllMortalityInNumbers(5)   + &
                                                            log10(1.-PopulationProcesses%m_shr/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_cra /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(6)  = PopulationProcesses%SumLogAllMortalityInNumbers(6)   + &
                                                            log10(1.-PopulationProcesses%m_cra/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_oys /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(7)  = PopulationProcesses%SumLogAllMortalityInNumbers(7)   + &
                                                            log10(1.-PopulationProcesses%m_oys/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_duck /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(8)  = PopulationProcesses%SumLogAllMortalityInNumbers(8)   + &
                                                            log10(1.-PopulationProcesses%m_duck/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_gull /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(9)  = PopulationProcesses%SumLogAllMortalityInNumbers(9)   + &
                                                            log10(1.-PopulationProcesses%m_gull/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_low /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(10) = PopulationProcesses%SumLogAllMortalityInNumbers(10)  + &
                                                            log10(1.-PopulationProcesses%m_low/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_self /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(11) = PopulationProcesses%SumLogAllMortalityInNumbers(11)  + &
                                                            log10(1.-PopulationProcesses%m_self/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_others /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(12) = PopulationProcesses%SumLogAllMortalityInNumbers(12)  + &
                                                        log10(1.-PopulationProcesses%m_others/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_vel /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(13) = PopulationProcesses%SumLogAllMortalityInNumbers(13)  + &
                                                            log10(1.-PopulationProcesses%m_vel/PopulationProcesses%TNStartTimeStep)
                end if

                if (1.-PopulationProcesses%m_settle /PopulationProcesses%TNStartTimeStep .gt. 0.0) then 
                    PopulationProcesses%SumLogAllMortalityInNumbers(14) = PopulationProcesses%SumLogAllMortalityInNumbers(14)  + &
                                                         log10(1.-PopulationProcesses%m_settle/PopulationProcesses%TNStartTimeStep)
                end if

                                                        
                !average mortality in numbers /d.m3
                PopulationProcesses%AverageMortalityInNumbers(1 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(1)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(2 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(2)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(3 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(3)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(4 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(4)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(5 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(5)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(6 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(6)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(7 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(7)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(8 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(8)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(9 )  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(9)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday 
                PopulationProcesses%AverageMortalityInNumbers(10)  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(10)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday
                PopulationProcesses%AverageMortalityInNumbers(11)  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(11)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday
            
                PopulationProcesses%AverageMortalityInNumbers(12)  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(12)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday

                PopulationProcesses%AverageMortalityInNumbers(13)  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(13)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday

                PopulationProcesses%AverageMortalityInNumbers(14)  = (1.-10**(PopulationProcesses%SumLogAllMortalityInNumbers(14)*&
                                                                    (1./PopulationProcesses%nInstantsForAverage))) / Me%DTday
                !Sum of mortalities in mass
                PopulationProcesses%SumAllMortalityInMass(1)  = PopulationProcesses%SumAllMortalityInMass(1) + &
                                                            PopulationProcesses%Massm_A    
                PopulationProcesses%SumAllMortalityInMass(2)  = PopulationProcesses%SumAllMortalityInMass(2) + &
                                                            PopulationProcesses%Massm_O    
                PopulationProcesses%SumAllMortalityInMass(3)  = PopulationProcesses%SumAllMortalityInMass(3) + &
                                                            PopulationProcesses%Massm_F    
                PopulationProcesses%SumAllMortalityInMass(4)  = PopulationProcesses%SumAllMortalityInMass(4) + &
                                                            PopulationProcesses%Massm_nat  
                PopulationProcesses%SumAllMortalityInMass(5)  = PopulationProcesses%SumAllMortalityInMass(5) + &
                                                            PopulationProcesses%Massm_shr  
                PopulationProcesses%SumAllMortalityInMass(6)  = PopulationProcesses%SumAllMortalityInMass(6) + &
                                                            PopulationProcesses%Massm_cra  
                PopulationProcesses%SumAllMortalityInMass(7)  = PopulationProcesses%SumAllMortalityInMass(7) + &
                                                            PopulationProcesses%Massm_oys  
                PopulationProcesses%SumAllMortalityInMass(8)  = PopulationProcesses%SumAllMortalityInMass(8) + &
                                                            PopulationProcesses%Massm_duck 
                PopulationProcesses%SumAllMortalityInMass(9)  = PopulationProcesses%SumAllMortalityInMass(9) + &
                                                            PopulationProcesses%Massm_gull 
                PopulationProcesses%SumAllMortalityInMass(10) = PopulationProcesses%SumAllMortalityInMass(10)+ &
                                                            PopulationProcesses%Massm_low  
                PopulationProcesses%SumAllMortalityInMass(11) = PopulationProcesses%SumAllMortalityInMass(11)+ &
                                                            PopulationProcesses%Massm_self  
                PopulationProcesses%SumAllMortalityInMass(12) = PopulationProcesses%SumAllMortalityInMass(12)+ &
                                                            PopulationProcesses%Massm_others  
                PopulationProcesses%SumAllMortalityInMass(13) = PopulationProcesses%SumAllMortalityInMass(13)+ &
                                                            PopulationProcesses%Massm_vel  
                PopulationProcesses%SumAllMortalityInMass(14) = PopulationProcesses%SumAllMortalityInMass(14)+ &
                                                            PopulationProcesses%Massm_settle  

            end if !tn>0
            
            Species => Species%Next
        end do d1

    end subroutine ComputePopulationVariables

    !--------------------------------------------------------------------------

    subroutine WriteOutput (Index, iIndexOutput)

        !Arguments-------------------------------------------------------------
        type(T_Species),      pointer   :: Species
        type(T_Cohort),       pointer   :: Cohort
        integer, intent(IN)             :: Index, iIndexOutput

        !Local-----------------------------------------------------------------
        integer                         :: M_V, M_E, M_H, M_R
        integer                         :: L, Age, Number 
        real                            :: Year, Month, Day, hour, minute, second
        real                            :: TotalSeconds

        !Begin-----------------------------------------------------------------

        call ExtractDate(Me%ExternalVar%CurrentTime, Year, Month, Day, hour, minute, second)
        
        TotalSeconds = Me%ExternalVar%CurrentTime - Me%InitialDate
                                                                          
        Species => Me%FirstSpecies
d1:     do while(associated(Species))

            if (Species%CohortOutput) then

                Cohort => Species%FirstCohort
d2:             do while(associated(Cohort))

                    M_V     = Cohort%StateIndex%M_V
                    M_E     = Cohort%StateIndex%M_E
                    M_H     = Cohort%StateIndex%M_H
                    M_R     = Cohort%StateIndex%M_R
                    L       = Cohort%StateIndex%L
                    Age     = Cohort%StateIndex%Age
                    Number  = Cohort%StateIndex%Number

                    !  header
                    !        write(Cohort%CohortOutput%Unit, 101)trim(adjustl(" Seconds YY1 MM2 DD3 hh4 mm5 ss6 "   // &
                    !           " Number7 ME8 MV9 MH10 MR11 L12 A13 Cr14 FInorg15 F16"        // &
                    !           " IInorg17 I18 PFInorg19 PF20 Ass21 FAEIng22 FAE23 JEM24 JE25 dE26"         // &
                    !           " GamLoss27 StruLoss28 JV29 JH30 MatLoss31 JS32 Gam33 JR34"   // &
                    !           " CO235 H2O36 O237 NH338 PO439"   // &
                    !           " m_A40 m_O41 m_F42 m_nat43 m_shr44 m_cra45 m_oys46 m_duck47 m_gull48  m_low49 m_self50
                    !           " m_ others m_vel49 m_settle50"))

                    102 format( F16.2  , 1x, I4   , 1x, I2,    1x, I2,    1x, I2,    1x, I2,    1x, I2, 1x , &  !7
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                                 , &  !10
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                      , &  !15
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                      , &  !20
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                      , &  !25
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                      , &  !30
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                      , &  !35
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                      , &  !40
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                      , &  !45
                                E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, &
                                E20.13, 1x, E20.13, 1x)!50

                    write(Cohort%CohortOutput%Unit(iIndexOutput), 102) TotalSeconds, int(Year), int(Month), int(Day),& !1,2,3
                                                         int(hour), int(minute), int(second)                        ,& !4,5,6
                                    Me%ExternalVar%Mass(Number,Index)    /Me%ConvertionFactor                       ,& !7, #/m2
                                    Me%ExternalVar%Mass(M_E,Index)                                                  ,& !8
                                    Me%ExternalVar%Mass(M_V,Index)       , Me%ExternalVar%Mass(M_H,Index)           ,& !9,10    
                                    Me%ExternalVar%Mass(M_R,Index)       , Me%ExternalVar%Mass(L,Index)             ,& !11,12    
                                    Me%ExternalVar%Mass(Age,Index)                                                  ,& !13
                                    Cohort%Processes%ClearanceRate       , Cohort%Processes%FilteredInorganic       ,& !14,15    
                                    Cohort%Processes%FilteredFood%C      , Cohort%Processes%IngestionInorganic      ,& !16,17
                                    Cohort%Processes%IngestionFood%C     , Cohort%Processes%PFContributionInorganic ,& !18,19
                                    Cohort%Processes%PFContributionFood%C, Cohort%Processes%Assimilation%C          ,& !20,21
                                    Cohort%Processes%FaecesContributionInorganic                                    ,& !22
                                    Cohort%Processes%FaecesContributionFood%C                                       ,& !23
                                    Cohort%Processes%SomaticMaintenance  , Cohort%Processes%Mobilization            ,& !24,25    
                                    Cohort%Processes%ReservesDynamics    , Cohort%Processes%GametesLoss             ,& !26,27    
                                    Cohort%Processes%StructureLoss       , Cohort%Processes%StructureDynamics       ,& !28,29    
                                    Cohort%Processes%MaturityDynamics    , Cohort%Processes%MaturityLoss            ,& !30,31    
                                    Cohort%Processes%Spawning            , Cohort%Processes%GametesToRelease        ,& !32,33    
                                    Cohort%Processes%ReproductionDynamics, Cohort%Processes%InorganicFluxes%CO2     ,& !34,35
                                    Cohort%Processes%InorganicFluxes%H2O , Cohort%Processes%InorganicFluxes%O2      ,& !36,37
                                    Cohort%Processes%InorganicFluxes%NH3 , Cohort%Processes%InorganicFluxes%PO4     ,& !38,39
                                    Cohort%Processes%DeathByAge                /Me%ConvertionFactor            ,& !40, #.d/m2
                                    Cohort%Processes%DeathByOxygen             /Me%ConvertionFactor            ,& !41, #.d/m2   
                                    Cohort%Processes%DeathByStarvation         /Me%ConvertionFactor            ,& !42, #.d/m2
                                    Cohort%Processes%DeathByNatural            /Me%ConvertionFactor            ,& !43, #.d/m2
                                    Cohort%Processes%PredationByShrimps        /Me%ConvertionFactor            ,& !44, #.d/m2
                                    Cohort%Processes%PredationByCrabs          /Me%ConvertionFactor            ,& !45, #.d/m2
                                    Cohort%Processes%PredationByOysterCatchers /Me%ConvertionFactor            ,& !46, #.d/m2
                                    Cohort%Processes%PredationByEiderDucks     /Me%ConvertionFactor            ,& !47, #.d/m2
                                    Cohort%Processes%PredationByHerringGull    /Me%ConvertionFactor            ,& !48, #.d/m2
                                    Cohort%Processes%DeathByLowNumbers         /Me%ConvertionFactor            ,& !49, #.d/m2
                                    Cohort%Processes%DeathBySelfPredation      /Me%ConvertionFactor            ,& !50, #.d/m2
                                    Cohort%Processes%DeathByLarvaePredationByOthers      /Me%ConvertionFactor   ,& !50, #.d/m2
                                    Cohort%Processes%DeathByVelocity           /Me%ConvertionFactor            ,& !49, #.d/m2
                                    Cohort%Processes%DeathByWrongSettlement    /Me%ConvertionFactor            ,& !50, #.d/m2
                                    Cohort%BivalveCondition%ScaledE                 !50, #.d/m2

                    Cohort => Cohort%Next
                end do d2
            end if
            
            if (Species%Population) then
 
!                OuputHeader =   "Seconds YY MM DD hh mm ss "                                        // &
!                                " TN NCoh TBio Cr Fil Ing Ass"                                      // &
!                                " CO H2O O NH3 PO4 LackOfFood"                                      // &
!                                " m_A m_O m_F m_nat m_shr m_cra m_oys m_duck m_gull  m_low m_self"  // &
!                                " m_vel m_settle"                                                   // &
!                                " TMASSm_A TMASSm_O TMASSm_F TMASSm_nat TMASSm_shr TMASSm_cra"      // &
!                                " TMASSm_oys TMASSm_duck TMASSm_gull TMASSm_low TMASSm_self"        // &
!                                "  TMASSm_low TMASSm_self
!                                " GEOm_A GEOm_O GEOm_F GEOm_nat GEOm_shr GEOm_cra GEOm_oys"         // &
!                                " GEOm_duck GEOm_gull GEOm_low GEOm_self"                           // &
!                                " TNField"                                                          // &
!                                " InitialPhyto InitialShrimp"                                       // &
!                                " MaxLength MaxTNFiels SpawningEvents"

                103 format( F16.2,    1x, I4,    1x, I2,    1x, I2,    1x, I2,    1x, I2,    1x, I2, 1x  , &  !6
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                     , &  !10
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x          , &  !15
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x          , &  !20
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x          , &  !25
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x          , &  !30
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x   , &  !35
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x          , &  !40
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x          , &  !45
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x          , &  !50
                            E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x                     , &  !53
                            E20.13, 1x, E20.13, 1x                                           , &  !Phyto E Shrimp
                            E20.13, 1x, E20.13, 1x, I4, 1x)                                               !54,55
                                            
                    write(Species%PopulationOutput%Unit(iIndexOutput), 103) TotalSeconds                 ,& !1
                        int(Year), int(Month), int(Day)                                                  ,& !2,3,4
                        int(hour), int(minute), int(second)                                              ,& !5,6,7
                        Species%PopulationProcesses%TN                            /Me%ConvertionFactor ,& !7, #/m2
                        Species%PopulationProcesses%NCoh                                                 ,& !8
                        Species%PopulationProcesses%TBio                          /Me%ConvertionFactor ,& !9, molC/m2
                        Species%PopulationProcesses%Cr                            /Me%ConvertionFactor ,& !10, molC/m2
                        Species%PopulationProcesses%Fil,   Species%PopulationProcesses%Ing               ,& !11,12, molC/m3
                        Species%PopulationProcesses%Ass                                                  ,& !13, molC/m3
                        Species%PopulationProcesses%CO2,   Species%PopulationProcesses%H2O               ,& !14,15, mol/m3
                        Species%PopulationProcesses%O  ,   Species%PopulationProcesses%NH3               ,& !16,17, mol/m3
                        Species%PopulationProcesses%PO4                                                  ,& !18, mol/m3
                        Me%LackOfFood                                                                    ,& !19
                        Species%PopulationProcesses%m_A                           /Me%ConvertionFactor ,& !20,#/d.m2
                        Species%PopulationProcesses%m_O                           /Me%ConvertionFactor ,& !21,#/d.m2
                        Species%PopulationProcesses%m_F                           /Me%ConvertionFactor ,& !22,#/d.m2
                        Species%PopulationProcesses%m_nat                         /Me%ConvertionFactor ,& !23,#/d.m2
                        Species%PopulationProcesses%m_shr                         /Me%ConvertionFactor ,& !24,#/d.m2
                        Species%PopulationProcesses%m_cra                         /Me%ConvertionFactor ,& !25,#/d.m2
                        Species%PopulationProcesses%m_oys                         /Me%ConvertionFactor ,& !26,#/d.m2
                        Species%PopulationProcesses%m_duck                        /Me%ConvertionFactor ,& !27,#/d.m2
                        Species%PopulationProcesses%m_gull                        /Me%ConvertionFactor ,& !28,#/d.m2
                        Species%PopulationProcesses%m_low                         /Me%ConvertionFactor ,& !29,#/d.m2    
                        Species%PopulationProcesses%m_self                        /Me%ConvertionFactor ,& !29,#/d.m2    
                        Species%PopulationProcesses%m_vel                         /Me%ConvertionFactor ,& !29,#/d.m2    
                        Species%PopulationProcesses%m_settle                      /Me%ConvertionFactor ,& !29,#/d.m2    
                        Species%PopulationProcesses%m_others                      /Me%ConvertionFactor ,& !29,#/d.m2    
                        Species%PopulationProcesses%SumAllMortalityInMass(1)      /Me%ConvertionFactor ,& !30,#/d.m2    
                        Species%PopulationProcesses%SumAllMortalityInMass(2)      /Me%ConvertionFactor ,& !31,#/d.m2    
                        Species%PopulationProcesses%SumAllMortalityInMass(3)      /Me%ConvertionFactor ,& !32,#/d.m2    
                        Species%PopulationProcesses%SumAllMortalityInMass(4)      /Me%ConvertionFactor ,& !33
                        Species%PopulationProcesses%SumAllMortalityInMass(5)      /Me%ConvertionFactor ,& !34    
                        Species%PopulationProcesses%SumAllMortalityInMass(6)      /Me%ConvertionFactor ,& !35    
                        Species%PopulationProcesses%SumAllMortalityInMass(7)      /Me%ConvertionFactor ,& !36
                        Species%PopulationProcesses%SumAllMortalityInMass(8)      /Me%ConvertionFactor ,& !37
                        Species%PopulationProcesses%SumAllMortalityInMass(9)      /Me%ConvertionFactor ,& !38
                        Species%PopulationProcesses%SumAllMortalityInMass(10)     /Me%ConvertionFactor ,& !39
                        Species%PopulationProcesses%SumAllMortalityInMass(11)     /Me%ConvertionFactor ,& !39
                        Species%PopulationProcesses%SumAllMortalityInMass(12)     /Me%ConvertionFactor ,& !39
                        Species%PopulationProcesses%SumAllMortalityInMass(13)     /Me%ConvertionFactor ,& !39
                        Species%PopulationProcesses%SumAllMortalityInMass(14)     /Me%ConvertionFactor ,& !39
                        Species%PopulationProcesses%AverageMortalityInNumbers(1)  /Me%ConvertionFactor ,& !40    
                        Species%PopulationProcesses%AverageMortalityInNumbers(2)  /Me%ConvertionFactor ,& !41    
                        Species%PopulationProcesses%AverageMortalityInNumbers(3)  /Me%ConvertionFactor ,& !42    
                        Species%PopulationProcesses%AverageMortalityInNumbers(4)  /Me%ConvertionFactor ,& !43
                        Species%PopulationProcesses%AverageMortalityInNumbers(5)  /Me%ConvertionFactor ,& !44    
                        Species%PopulationProcesses%AverageMortalityInNumbers(6)  /Me%ConvertionFactor ,& !45    
                        Species%PopulationProcesses%AverageMortalityInNumbers(7)  /Me%ConvertionFactor ,& !46
                        Species%PopulationProcesses%AverageMortalityInNumbers(8)  /Me%ConvertionFactor ,& !47
                        Species%PopulationProcesses%AverageMortalityInNumbers(9)  /Me%ConvertionFactor ,& !48
                        Species%PopulationProcesses%AverageMortalityInNumbers(10) /Me%ConvertionFactor ,& !49
                        Species%PopulationProcesses%AverageMortalityInNumbers(11) /Me%ConvertionFactor ,& !49
                        Species%PopulationProcesses%AverageMortalityInNumbers(12) /Me%ConvertionFactor ,& !49
                        Species%PopulationProcesses%AverageMortalityInNumbers(13) /Me%ConvertionFactor ,& !49
                        Species%PopulationProcesses%AverageMortalityInNumbers(14) /Me%ConvertionFactor ,& !49
                        Species%PopulationProcesses%TNField                       /Me%ConvertionFactor ,& !50,#/d.m2
                        Me%ExternalVar%InitialPhyto(Index), Me%ExternalVar%InitialShrimp(Index)             ,& !Phyto e Shrimp
                        Species%PopulationProcesses%MaxLength                                               ,& !51
                        Me%MaxTNField                                             /Me%ConvertionFactor ,& !51
                        Species%PopulationProcesses%nSpawning                                                !52
                                                                    
                if (Species%BySizeOutput) then

                    call WriteSizeDistribution (Index, Species, iIndexOutput)

                end if
            
            end if !population
            
            Species => Species%Next
        end do d1

        if (Me%ComputeOptions%MassBalance) call WriteMassBalance (Index, iIndexOutput)

    end subroutine WriteOutput

    !--------------------------------------------------------------------------

    subroutine WriteMassBalance (Index, iIndexOutput)

        !Notes-----------------------------------------------------------------

        !Arguments-------------------------------------------------------------
        integer, intent(IN)             :: Index, iIndexOutput

        !Local-----------------------------------------------------------------
        type(T_Species),      pointer   :: Species
        type(T_Cohort),       pointer   :: Cohort
        real                            :: RATIOHC, RATIOOC, RATIONC, RATIOPC 
        real                            :: SumC, SumN, SumP  
        real                            :: C_AtomicMass, H_AtomicMass, O_AtomicMass, N_AtomicMass, P_AtomicMass
        integer                         :: M_V, M_E, M_R, Number    
        real                            :: Year, Month, Day, hour, minute, second

        !Begin-----------------------------------------------------------------

        call ExtractDate(Me%ExternalVar%CurrentTime, Year, Month, Day, hour, minute, second)

        !Assumed the samecomposition as the default values for bivalves if wq
        !The verification of the mass is only possible if the ratios are the same for bivalve and phytoplankton
        ! or running with life because wq does not follow POC and the ratio PON/POC is not known
        ! 0.3395653308501444 is an approximate value... 
        !Should be tested with life
        RATIOHC  = 0.15         
        RATIOOC  = 0.71         
        RATIONC  = 0.3        
        RATIOPC  = 0.07  

        !in g
        SumC = 0.0
        SumN = 0.0
        SumP = 0.0
        
        !folowing the constructed property list
        if (Me%ComputeOptions%PelagicModel .eq. LifeModel) then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%POC,Index) 

        else 

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%PON,Index) / RATIONC  

        end if

        SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%CarbonDioxide,Index)

        if(Me%ComputeOptions%Nitrogen)then   
            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%AM,Index)
            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%PON,Index)  
        end if

        if(Me%ComputeOptions%Phosphorus)then
            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%IP,Index)
            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%POP,Index)
        end if

        !Include All the cohorts from all the species      
        Species => Me%FirstSpecies
        do while(associated(Species))

            C_AtomicMass  = Species%AuxiliarParameters%C_AtomicMass        
            H_AtomicMass  = Species%AuxiliarParameters%H_AtomicMass        
            O_AtomicMass  = Species%AuxiliarParameters%O_AtomicMass        
            P_AtomicMass  = Species%AuxiliarParameters%P_AtomicMass        
            N_AtomicMass  = Species%AuxiliarParameters%N_AtomicMass        

            Cohort => Species%FirstCohort
            do while(associated(Cohort))

                M_V     = Cohort%StateIndex%M_V
                M_E     = Cohort%StateIndex%M_E
                M_R     = Cohort%StateIndex%M_R
                Number  = Cohort%StateIndex%Number

                SumC     = SumC + Me%ExternalVar%Mass(M_V,Index)                    * &
                                  Me%ExternalVar%Mass(Number,Index) * C_AtomicMass
                                  
                SumC     = SumC + Me%ExternalVar%Mass(M_E,Index)                    * &
                                  Me%ExternalVar%Mass(Number,Index) * C_AtomicMass
                                  
                SumC     = SumC + Me%ExternalVar%Mass(M_R,Index)                    * &
                                  Me%ExternalVar%Mass(Number,Index) * C_AtomicMass

                if(Me%ComputeOptions%Nitrogen)then   

                    SumN = SumN + Me%ExternalVar%Mass(M_V,Index)                     * &
                                  Species%SpeciesComposition%StructureComposition%nN * &
                                  Me%ExternalVar%Mass(Number,Index) * N_AtomicMass

                    SumN = SumN + Me%ExternalVar%Mass(M_E,Index)                     * &
                                  Species%SpeciesComposition%ReservesComposition%nN  * &
                                  Me%ExternalVar%Mass(Number,Index) * N_AtomicMass

                    SumN = SumN + Me%ExternalVar%Mass(M_R,Index)                     * &
                                  Species%SpeciesComposition%ReservesComposition%nN  * &
                                  Me%ExternalVar%Mass(Number,Index) * N_AtomicMass

                end if

                if(Me%ComputeOptions%Phosphorus)then

                    SumP = SumP + Me%ExternalVar%Mass(M_V,Index)                     * &
                                  Species%SpeciesComposition%StructureComposition%nP * &
                                  Me%ExternalVar%Mass(Number,Index) * P_AtomicMass

                    SumP = SumP + Me%ExternalVar%Mass(M_E,Index)                     * &
                                  Species%SpeciesComposition%ReservesComposition%nP  * &
                                  Me%ExternalVar%Mass(Number,Index) * P_AtomicMass

                    SumP = SumP + Me%ExternalVar%Mass(M_R,Index)                     * &
                                  Species%SpeciesComposition%ReservesComposition%nP  * &
                                  Me%ExternalVar%Mass(Number,Index) * P_AtomicMass

                end if

                Cohort => Cohort%Next
            end do

            SumC = SumC + Species%PopulationProcesses%nNewborns                     * &
                            (Species%IndividualParameters%MEb                       * &
                            Species%SpeciesComposition%ReservesComposition%nC       + &
                            Species%IndividualParameters%MVb                        * &
                            Species%SpeciesComposition%StructureComposition%nC) * C_AtomicMass

            SumN = SumN + Species%PopulationProcesses%nNewborns                     * &
                            (Species%IndividualParameters%MEb                       * &
                            Species%SpeciesComposition%ReservesComposition%nN       + &
                            Species%IndividualParameters%MVb                        * &
                            Species%SpeciesComposition%StructureComposition%nN) * N_AtomicMass

            SumP = SumP + Species%PopulationProcesses%nNewborns                     * &
                            (Species%IndividualParameters%MEb                       * &
                             Species%SpeciesComposition%ReservesComposition%nP      + &
                             Species%IndividualParameters%MVb                       * &
                             Species%SpeciesComposition%StructureComposition%nP) * P_AtomicMass


            Species => Species%Next
        end do

        if(Me%PropIndex%phyto .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%phyto,Index)

            if(Me%ComputeOptions%Nitrogen)then   
                SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%phyto,Index) * RATIONC  
            end if

            if(Me%ComputeOptions%Phosphorus)then
                SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%phyto,Index) * RATIOPC
            end if

        end if

        if(Me%PropIndex%zoo .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%zoo,Index)

            if(Me%ComputeOptions%Nitrogen)then   
                SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%zoo,Index) * RATIONC 
            end if

            if(Me%ComputeOptions%Phosphorus)then
                SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%zoo,Index) * RATIOPC 
            end if

        end if   

        if(Me%PropIndex%ciliate .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%ciliate,Index)

            if(Me%ComputeOptions%Nitrogen)then   
                SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%ciliate,Index) * RATIONC 
            end if

            if(Me%ComputeOptions%Phosphorus)then
                SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%ciliate,Index) * RATIOPC 
            end if

        end if      

        if(Me%PropIndex%bacteria .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%bacteria,Index)

            if(Me%ComputeOptions%Nitrogen)then   
                SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%bacteria,Index) * RATIONC 
            end if

            if(Me%ComputeOptions%Phosphorus)then
                SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%bacteria,Index) * RATIOPC 
            end if

        end if      

        if(Me%PropIndex%DiatomsC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%DiatomsC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%DiatomsN,Index) 

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%DiatomsP,Index) 

        end if      

        if(Me%PropIndex%diatoms .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%diatoms,Index)

            if(Me%ComputeOptions%Nitrogen)then   
                SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%diatoms,Index) * RATIONC 
            end if

            if(Me%ComputeOptions%Phosphorus)then
                SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%diatoms,Index) * RATIOPC 
            end if

        end if      

        if(Me%PropIndex%Mix_FlagellateC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%Mix_FlagellateC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%Mix_FlagellateN,Index)  

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%Mix_FlagellateP,Index) 

        end if      


        if(Me%PropIndex%PicoalgaeC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%PicoalgaeC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%PicoalgaeN,Index)  

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%PicoalgaeP,Index) 

        end if      

        if(Me%PropIndex%FlagellateC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%FlagellateC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%FlagellateN,Index) 

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%FlagellateP,Index)

        end if      

        if(Me%PropIndex%MicrozooplanktonC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%MicrozooplanktonC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%MicrozooplanktonN,Index) 

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%MicrozooplanktonP,Index) 

        end if      

        if(Me%PropIndex%Het_NanoflagellateC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%Het_NanoflagellateC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%Het_NanoflagellateN,Index)  

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%Het_NanoflagellateP,Index)  

        end if      

        if(Me%PropIndex%MesozooplanktonC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%MesozooplanktonC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%MesozooplanktonN,Index) 

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%MesozooplanktonP,Index) 

        end if      


        if(Me%PropIndex%Het_BacteriaC .ne. null_int)then

            SumC = SumC + Me%ExternalVar%Mass(Me%PropIndex%Het_BacteriaC,Index)

            SumN = SumN + Me%ExternalVar%Mass(Me%PropIndex%Het_BacteriaN,Index)  

            SumP = SumP + Me%ExternalVar%Mass(Me%PropIndex%Het_BacteriaP,Index)  

        end if      


        101 format(I4,    1x, I2,    1x, I2,    1x, I2,    1x, I2, 1x, I2, 1x, &  !6 
                   E16.8, 1x, E16.8, 1x, E16.8, 1x, E16.8, 1x)                    !10

        102 format(I4,    1x, I2,    1x, I2,    1x, I2,    1x, I2, 1x, I2, 1x  ,& !6 
                   E16.8, 1x, E16.8, 1x, E16.8, 1x)                               !9

        if (Me%ComputeOptions%PelagicModel .eq. LifeModel) then

            write(Me%MassOutput%Unit(iIndexOutput), 101)  int(Year), int(Month), int(Day)   ,& !1,2,3
                                            int(hour), int(minute), int(second)             ,& !4,5,6
                                            SumC, SumN, SumP, Me%MassLoss
        else

            write(Me%MassOutput%Unit(iIndexOutput), 102)  int(Year), int(Month), int(Day)   ,& !1,2,3
                                            int(hour), int(minute), int(second)             ,& !4,5,6
                                            SumN, SumP, Me%MassLoss
        end if

    end subroutine WriteMassBalance

    !--------------------------------------------------------------------------

    subroutine WriteSizeDistribution(Index, Species, iIndexOutput)

        !Notes-----------------------------------------------------------------

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: Index, iIndexOutput
        type(T_Species)          , pointer  :: Species

        !Local-----------------------------------------------------------------            
        type(T_Cohort)           , pointer  :: Cohort
        integer                             :: L, Number, iSizeClass    
        real                                :: lowsize, highsize 
        real                                :: Year, Month, Day, hour, minute, second
        real                                :: TotalSeconds


        !Begin-----------------------------------------------------------------

        call ExtractDate(Me%ExternalVar%CurrentTime, Year, Month, Day, hour, minute, second)

        TotalSeconds = Me%ExternalVar%CurrentTime - Me%InitialDate

        !Initialize
        Species%SizeFrequency = 0.0

        do iSizeClass = 1, Species%nSizeClasses

            lowsize  = Species%SizeClasses(iSizeClass)

            if (iSizeClass .lt. Species%nSizeClasses) then
                highsize = Species%SizeClasses(iSizeClass+1)
            else
                highsize = 1000.           
            end if

            Cohort => Species%FirstCohort
            do while(associated(Cohort))

                L       = Cohort%StateIndex%L
                Number  = Cohort%StateIndex%Number

                if ((Me%ExternalVar%Mass(L,Index) .ge. lowsize) .and. (Me%ExternalVar%Mass(L,Index) .lt. highsize)) then

                    !m2
                    Species%SizeFrequency(iSizeClass) = Species%SizeFrequency(iSizeClass) + Me%ExternalVar%Mass(Number,Index) &
                                                                                            /Me%ConvertionFactor

                end if

                Cohort => Cohort%Next
            end do

        end do

        write(Species%SizeDistributionOutput%Unit(iIndexOutput), '(F16.2, 1x, I4,1x,I2,1x, I2,1x, I2,1x, I2,1x, I2, 1x ,'//&
                                                                Species%SizeClassesCharFormat//'E16.8, 1x)')  &
                                                                TotalSeconds, int(Year), int(Month), int(Day) , & !1,2,3
                                                                int(hour), int(minute), int(second), & !4,5,6
                                                                Species%SizeFrequency

    end subroutine WriteSizeDistribution

    !--------------------------------------------------------------------------

    subroutine WriteTestingFile(Index, iIndexOutput)

        !Notes-----------------------------------------------------------------

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                 :: Index, iIndexOutput

        !Local-----------------------------------------------------------------            
        integer                             :: STAT_CALL
        type(T_Cohort)           , pointer  :: Cohort
        type(T_Species)          , pointer  :: Species
        type(T_Predator)         , pointer  :: Predator
        integer                             :: L, Number    
        real                                :: TN,  Maxlength
        real                                :: m_A, m_O, m_F, m_nat
        real                                :: m_shr, m_cra, m_oys, m_duck, m_gull, m_low, m_self
        real                                :: m_others, m_vel, m_settle        
        integer                             :: TC,Evaluation 
        !real                                :: Year, Month, Day, hour, minute, second
        character(len=1000)                  :: ParameterValueStr
        character(len=1000)                  :: CompleteLineStr

        
        !Begin-----------------------------------------------------------------
        
        Species => Me%FirstSpecies
        do while(associated(Species))

            call UnitsManager(Species%TestingParametersOutput%Unit(iIndexOutput), OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine WriteTestingFile - ModuleBivalve - ERR01'
            
            
            open(Unit = Species%TestingParametersOutput%Unit(iIndexOutput), File = trim(Species%Testing_File)   , &
                                                              Status = 'UNKNOWN', Access = 'APPEND' )
            
            TN = 0.0
            TC = 0.0
            
            if (Species%nCohorts .eq. 0) then
            

                TN = 0.0
                TC = 0.0
                
                Maxlength = Species%PopulationProcesses%LastLength
                
            else
            
                Cohort => Species%FirstCohort
                do while(associated(Cohort))

                    L       = Cohort%StateIndex%L
                    Number  = Cohort%StateIndex%Number

                    TN = TN + Me%ExternalVar%Mass(Number,Index)
                    TC = TC + 1
                    
                    Maxlength = max(Maxlength,Me%ExternalVar%Mass(L,Index))
                     
                    Cohort => Cohort%Next
                end do
                
                if (Maxlength .lt. 1.0e-15) then           
                    Maxlength = 0.0
                end if  
                
                m_A    = 0.0
                m_O    = 0.0
                m_F    = 0.0
                m_nat  = 0.0
                m_shr  = 0.0
                m_cra  = 0.0
                m_oys  = 0.0
                m_duck = 0.0
                m_gull = 0.0
                m_low  = 0.0
                m_self  = 0.0
                m_others  = 0.0
                m_vel  = 0.0
                m_settle  = 0.0
                
            end if
            
            CompleteLineStr = trim(Species%ID%Name)//' '
            
            Predator => Species%FirstPredator
            do while(associated(Predator))
            
                write(ParameterValueStr, ('(F7.4)')) Predator%Diet
                
                CompleteLineStr = trim(CompleteLineStr)//' '//trim(ParameterValueStr)
            
                Predator => Predator%Next
            end do
                
            150 format(E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, I5, 1x     , &
                       E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x             , & !5
                       E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x,   &
                       E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x,E20.13, 1x,   &
                       E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x, E20.13, 1x,E20.13, 1x,   &
                       E20.13, 1x, E20.13, 1x,E20.13, 1x,I5, 1x,E20.13, 1x)        !6
            
            write(ParameterValueStr, 150) Me%DT                                                                          , &   
                                          Species%IndividualParameters%m_spat                      /Me%ConvertionFactor  , &
                                          Species%IndividualParameters%m_natural                   /Me%ConvertionFactor  , &
                                          TN                                                       /Me%ConvertionFactor  , &
                                          Maxlength                                                                      , &
                                          TC                                                                             , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(1)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(2)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(3)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(4)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(5)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(6)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(7)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(8)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(9)      /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(10)     /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(11)     /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(12)     /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(13)     /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%SumAllMortalityInMass(14)     /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(1)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(2)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(3)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(4)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(5)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(6)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(7)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(8)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(9)  /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(10) /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(11) /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(12) /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(13) /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%AverageMortalityInNumbers(14) /Me%ConvertionFactor , &
                                          Species%PopulationProcesses%nSpawning                                          , &
                                          Me%MaxTNField                                             /Me%ConvertionFactor
                                         
            CompleteLineStr = trim(CompleteLineStr)//' '//trim(ParameterValueStr)

            !Compute Evaluation
            if (TN .gt. 0.0) then
                Evaluation = 1 !Alive
            else
                if ( (m_shr+m_cra+m_oys+m_duck+m_gull+m_self+m_others) .gt. 0.0 ) then
                    Evaluation = -1 !Died from some predation
                else
                    if (m_F .gt. 0.0 ) then
                        Evaluation = 0  !Died from starvation
                    else
                        if (m_nat .gt. 0.0 ) then
                            Evaluation = -2  !Died from natural causes
                        else
                            if (m_A .gt. 0.0 ) then
                                Evaluation = -3  !Died from age
                            else
                                if (m_low .gt. 0.0 ) then
                                    Evaluation = -4  !Died from low numbers
                                else
                                    Evaluation = -5  !Died from oxygen depletion
                                end if
                            end if
                        end if
                    end if
                end if
            end if
    
            write(ParameterValueStr, ('(I5)')) Evaluation

            CompleteLineStr = trim(CompleteLineStr)//' '//trim(ParameterValueStr)

            101 format(A1000)             
            
            write(Species%TestingParametersOutput%Unit(iIndexOutput), 101) CompleteLineStr
                 
            Species => Species%Next
        end do

    end subroutine WriteTestingFile

    !--------------------------------------------------------------------------

    subroutine UpdateCohortState
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Cohort)           , pointer          :: Cohort
        type(T_Species)          , pointer          :: Species
        !----------------------------------------------------------------------
  
        Species  => Me%FirstSpecies
d1:     do while (associated(Species))

            Cohort  => Species%FirstCohort
d2:         do while (associated(Cohort))

                Cohort%GlobalDeath = Cohort%GlobalDeath * Cohort%Dead

                Cohort  => Cohort%Next  
            end do  d2        

            Species  => Species%Next  
        enddo  d1

    end subroutine UpdateCohortState

    !--------------------------------------------------------------------------

    subroutine RestoreUnits (Index)
        !Arguments-------------------------------------------------------------
        integer                                     :: Index
        
        !Local-----------------------------------------------------------------
        type(T_Species)         , pointer           :: Species
        type(T_Cohort)          , pointer           :: Cohort
        integer                                     :: Number
        integer                                     :: Shrimp, Crab, OysterCatcher
        integer                                     :: EiderDuck, HerringGull

        !Begin-----------------------------------------------------------------
        
        Shrimp        = Me%PropIndex%Shrimp
        Crab          = Me%PropIndex%Crab
        OysterCatcher = Me%PropIndex%OysterCatcher
        EiderDuck     = Me%PropIndex%EiderDuck
        HerringGull   = Me%PropIndex%HerringGull

        Species => Me%FirstSpecies
        do while(associated(Species))
            
            Cohort => Species%FirstCohort
            do while(associated(Cohort))
            
                Number = Cohort%StateIndex%Number 
                         
                !Reconvert the mussels organisms density by applying 1.0/ConvertionFactor
                
                    Me%ExternalVar%Mass(Number,Index) = Me%ExternalVar%Mass(Number,Index)/Me%ConvertionFactor
                
                Cohort => Cohort%Next

            end do
                                    
            Species => Species%Next

        end do
        
        !Convert the predators organisms density to the units in the begining
        if(Shrimp .ne. null_int)then

            Me%ExternalVar%Mass(Shrimp,Index) = Me%ExternalVar%Mass(Shrimp,Index)               * &
                                                1./Me%ConvertionFactor
        end if

        if(Crab .ne. null_int)then

            Me%ExternalVar%Mass(Crab,Index) = Me%ExternalVar%Mass(Crab,Index)                   * &
                                              1./Me%ConvertionFactor
        end if
        
        if(Me%PropIndex%OysterCatcher .ne. null_int)then

            Me%ExternalVar%Mass(OysterCatcher,Index) = Me%ExternalVar%Mass(OysterCatcher,Index) * &
                                                       1./Me%ConvertionFactor
        end if

        if(Me%PropIndex%EiderDuck .ne. null_int)then

            Me%ExternalVar%Mass(EiderDuck,Index) = Me%ExternalVar%Mass(OysterCatcher,Index)     * &
                                                   1./Me%ConvertionFactor
        end if
        
        
        if(Me%PropIndex%HerringGull .ne. null_int)then

            Me%ExternalVar%Mass(HerringGull,Index) = Me%ExternalVar%Mass(HerringGull,Index)     * &
                                                     1./Me%ConvertionFactor
                                                        
        end if

    end subroutine RestoreUnits

    !--------------------------------------------------------------------------

    subroutine UpdateListDeadAndNewBornIDs
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: DynamicCohortPropertyID
        type(T_Cohort)           , pointer          :: Cohort
        type(T_Species)          , pointer          :: Species
        !----------------------------------------------------------------------
  
        Species  => Me%FirstSpecies
d1:     do while (associated(Species))

            if (Species%NewbornCohort) then
            
                !Add the PopulationProcesses%nNewborns of this species to the ListNewbornsIDs   
                Me%nLastNewbornsID = Me%nLastNewbornsID + 1
                Me%ListNewbornsIDs(Me%nLastNewbornsID) = Species%ID%IDNumber
                
            end if
                
            Cohort  => Species%FirstCohort
d2:         do while (associated(Cohort))

                if (Cohort%GlobalDeath .eq. 1) then
                 
                    Me%nLastDeadID = Me%nLastDeadID + 1
                    DynamicCohortPropertyID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" structure")
                    Me%ListDeadIDs(Me%nLastDeadID) = DynamicCohortPropertyID

                    Me%nLastDeadID = Me%nLastDeadID + 1
                    DynamicCohortPropertyID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" reserves")
                    Me%ListDeadIDs(Me%nLastDeadID) = DynamicCohortPropertyID

                    Me%nLastDeadID = Me%nLastDeadID + 1
                    DynamicCohortPropertyID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" maturity")
                    Me%ListDeadIDs(Me%nLastDeadID) = DynamicCohortPropertyID

                    Me%nLastDeadID = Me%nLastDeadID + 1
                    DynamicCohortPropertyID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" reproduction")
                    Me%ListDeadIDs(Me%nLastDeadID) = DynamicCohortPropertyID

                    Me%nLastDeadID = Me%nLastDeadID + 1
                    DynamicCohortPropertyID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" length")
                    Me%ListDeadIDs(Me%nLastDeadID) = DynamicCohortPropertyID

                    Me%nLastDeadID = Me%nLastDeadID + 1
                    DynamicCohortPropertyID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" age")
                    Me%ListDeadIDs(Me%nLastDeadID) = DynamicCohortPropertyID

                    Me%nLastDeadID = Me%nLastDeadID + 1
                    DynamicCohortPropertyID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" number")
                    Me%ListDeadIDs(Me%nLastDeadID) = DynamicCohortPropertyID
                                
                end if 

                Cohort  => Cohort%Next  

            end do  d2        

            Species  => Species%Next  
        enddo  d1
        
    end subroutine UpdateListDeadAndNewBornIDs

    !--------------------------------------------------------------------------

    subroutine UpdateBivalvePropertyList
        
        !Arguments------------------------------------------------------------------

        !Local----------------------------------------------------------------------

        !---------------------------------------------------------------------------

        call RemoveDeadIDsFromSpecies

        call AddNewbornsToSpecies

        call PropertyIndexNumber

        call ConstructPropertyList

    end subroutine UpdateBivalvePropertyList

    !--------------------------------------------------------------------------

    subroutine RemoveDeadIDsFromSpecies 

        !Local-----------------------------------------------------------------
        type (T_Species)           , pointer      :: Species
        type (T_Cohort)            , pointer      :: Cohort
        type (T_Cohort)            , pointer      :: CohortToContinue
        integer                                   :: iDeadIDs = 0,DynamicCohortID 
        !----------------------------------------------------------------------
        
        nullify(CohortToContinue)

        Species  => Me%FirstSpecies
d1:     do while (associated(Species))

            Cohort  => Species%FirstCohort
d2:         do while (associated(Cohort))

                iDeadIDs = 1
                DynamicCohortID = GetDynamicPropertyIDNumber(trim(adjustl(Cohort%ID%Name))//" structure")

                do iDeadIDs = 1, Me%nLastDeadID

                    if (Me%ListDeadIDs(iDeadIDs) .eq. DynamicCohortID) then
                    
                        nullify(CohortToContinue)

                        CohortToContinue => Cohort%Next

                        call RemoveCohortFromList (Species, Cohort)

                        exit
                                            
                    endif

                end do

                if(associated(CohortToContinue)) then 
                    Cohort => CohortToContinue
                    nullify(CohortToContinue)
                else
                    if (associated(Cohort)) Cohort  => Cohort%Next
                end if

            end do  d2        

            Species  => Species%Next  
        enddo  d1

    end subroutine RemoveDeadIDsFromSpecies

    !--------------------------------------------------------------------------       

    subroutine RemoveCohortFromList (ObjSpecies, ObjCohort)

        !Arguments-----------------------------------------------------------------
        type (T_Species), pointer     :: ObjSpecies
        type (T_Cohort), pointer      :: ObjCohort         

        !Local-----------------------------------------------------------------
        type (T_Cohort), pointer      :: Cohort           => null()
        type (T_Cohort), pointer      :: PreviousCohort   => null()
        integer                       :: STAT_CALL, iIndexOutput
        !----------------------------------------------------------------------

        nullify(Cohort, PreviousCohort)

        Cohort  => ObjSpecies%FirstCohort
        do while (associated(Cohort))

            if (Cohort%ID%ID .eq. ObjCohort%ID%ID) then  !this is the cohort to be removed            

                if(associated(PreviousCohort))then
                    PreviousCohort%Next      => Cohort%Next
                else
                    ObjSpecies%FirstCohort   => Cohort%Next
                end if

                ObjSpecies%nCohorts  = ObjSpecies%nCohorts - 1
                
                if (ObjSpecies%CohortOutput) then
                
                    do iIndexOutput = 1, Me%nIndexOutputs
                
                        call WriteDataLine(Cohort%CohortOutput%Unit(iIndexOutput), '<EndTimeSerie>')

                        call UnitsManager(Cohort%CohortOutput%Unit(iIndexOutput), CLOSE_FILE, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine RemoveCohortFromList - ModuleBivalve - ERR01'
                    
                    enddo
                    
                    deallocate(Cohort%CohortOutput%Unit)

                end if

                if(associated(Cohort%FeedingOn)) deallocate (Cohort%FeedingOn)
                
                deallocate    (Cohort)
                nullify       (Cohort)

                cycle

            else

                PreviousCohort => Cohort
                Cohort  => Cohort%Next

            endif

        enddo 

        nullify(ObjCohort)

    end subroutine RemoveCohortFromList

    !--------------------------------------------------------------------------

    subroutine AddNewbornsToSpecies 

        !Local-----------------------------------------------------------------
        type (T_Species), pointer                   :: Species
        
        Me%nLastNewbornsID = 0

        Species  => Me%FirstSpecies
d1:     do while (associated(Species))

            if (Species%NewbornCohort) then
            
                !Add the PopulationProcesses%nNewborns of this species to the ListNewbornsIDs   
                Me%nLastNewbornsID = Me%nLastNewbornsID + 1
                Me%ListNewbornsIDs(Me%nLastNewbornsID) = Species%ID%IDNumber

                call ConstructCohort (Species)
            
            end if
            
            Species  => Species%Next  
        enddo  d1

    end subroutine AddNewbornsToSpecies  

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DEST

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillBivalve(ObjBivalveID, STAT)

        !Arguments--------------------------------------------------------------------
        integer :: ObjBivalveID
        integer, optional, intent(OUT)      :: STAT

        !External---------------------------------------------------------------------
        integer :: ready_

        !Local------------------------------------------------------------------------
        integer :: STAT_, nUsers           
        type (T_Species)  , pointer         :: Species
        type (T_Cohort)   , pointer         :: Cohort
        integer                             :: STAT_CALL, iIndexOutput
        !-----------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBivalveID, ready_)
        
        if (Me%Testing_Parameters) then
        
            Species => Me%FirstSpecies
            do while(associated(Species))
        
                allocate(Species%TestingParametersOutput%Unit(1:Me%nIndexOutputs))
                
                Species => Species%Next
                
            enddo
            
            do iIndexOutput = 1, Me%nIndexOutputs
        
                call WriteTestingFile (Me%IndexOutputs(iIndexOutput), iIndexOutput)

            enddo
            
        end if

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mBivalve_,  Me%InstanceID)

            if (nUsers == 0) then
            
                call WriteFinalBivalveFile
            
                !deallocate everything inside the cohorts and species
                Species => Me%FirstSpecies
                do while(associated(Species))

                    Cohort => Species%FirstCohort
                    do while(associated(Cohort))

                        !soffs
                        if(associated(Cohort%FeedingOn)) deallocate (Cohort%FeedingOn)
                        
                        !deallocate(Cohort%LarvaeState)

                        !deallocate(Cohort%FeedingOn)

                    Cohort  => Cohort%Next  
                    enddo
                    
                    if (Me%OutputON) call CloseFiles (Species)

                    deallocate (Species%FirstCohort)
                    nullify    (Species%FirstCohort)
                    
                    deallocate (Species%FirstParticles)
                    nullify    (Species%FirstParticles)

                    deallocate (Species%FirstPredator)
                    nullify    (Species%FirstPredator)

                    deallocate (Species%SizeClasses)
                    deallocate (Species%SizeFrequency)

                    deallocate (Species%PopulationProcesses%SumLogAllMortalityInNumbers)
                    deallocate (Species%PopulationProcesses%SumAllMortalityInMass)
                    deallocate (Species%PopulationProcesses%AverageMortalityInNumbers)
                    
                    deallocate (Species%SettlementProbability)
                    nullify    (Species%SettlementProbability)

                   
                Species  => Species%Next  
                enddo
                
                if (Me%ComputeOptions%MassBalance) then
                
                    do iIndexOutput = 1, Me%nIndexOutputs

                        call UnitsManager(Me%MassOutput%Unit(iIndexOutput), CLOSE_FILE, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillBivalve - ModuleBivalve - ERR00'
                    
                    enddo
                    
                end if
                
                deallocate (Me%FirstSpecies)

                if(associated(Me%ListDeadIDs )) deallocate (Me%ListDeadIDs)
                if(associated(Me%MatrixNewborns)) deallocate (Me%MatrixNewborns)
                
                deallocate(Me%ExternalVar%InitialPhyto )
                deallocate(Me%ExternalVar%InitialShrimp)

                if (Me%ObjTime /= 0) then
                    nUsers = DeassociateInstance(mTIME_, Me%ObjTime)
                    if (nUsers == 0) stop 'KillBivalve - ModuleBivalve - ERR01'
                endif    

                !Deallocates Instance
                call DeallocateInstance ()

                ObjBivalveID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    !---------------------------------------------------------------------------

    end subroutine KillBivalve

    !-------------------------------------------------------------------------------
    
    subroutine WriteFinalBivalveFile 
       
        !Local----------------------------------------------------------------------
        integer                                   :: STAT_CALL, FinalUnit
        type (T_Species)         , pointer        :: Species
        type (T_Cohort )         , pointer        :: Cohort

        !Begin----------------------------------------------------------------------

        call UnitsManager(FinalUnit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteFinalBivalveFile - ModuleBivalve - ERR01'

        open(Unit = FinalUnit, File = trim(Me%FinalFileName), status = 'UNKNOWN', IOSTAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteFinalBivalveFile - ModuleBivalve - ERR10'

        call WriteDataLine(FinalUnit, 'NUMBER_OF_SPECIES', Me%nSpecies)

        Species => Me%FirstSpecies
        do while(associated(Species))
        
            call WriteDataLine(FinalUnit, '<begin_species>') 
            call WriteDataLine(FinalUnit, 'NAME', trim(Species%ID%Name))
            call WriteDataLine(FinalUnit, 'NUMBER_OF_COHORTS', Species%nCohorts)

            Cohort => Species%FirstCohort
            do while(associated(Cohort))
            
                call WriteDataLine(FinalUnit, '<begin_cohort>') 
                write(FinalUnit, *) Cohort%ID%ID
                call WriteDataLine(FinalUnit, '<end_cohort>') 

                Cohort  => Cohort%Next  
            enddo

            call WriteDataLine(FinalUnit, '<end_species>') 

            Species  => Species%Next  
        enddo

        call UnitsManager(FinalUnit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteFinalBivalveFile - ModuleBivalve - ERR20'
    
    end subroutine WriteFinalBivalveFile
    
    !-------------------------------------------------------------------------------

    subroutine CloseFiles (Species)

        !Local----------------------------------------------------------------------
        integer                                   :: STAT_CALL
        type (T_Species)         , pointer        :: Species
        type (T_Cohort )         , pointer        :: Cohort
        integer                                   :: iIndexOutput

        !Begin----------------------------------------------------------------------

                
        Cohort  => Species%FirstCohort
        do while (associated(Cohort))

            if (Species%CohortOutput) then
            
                do iIndexOutput = 1, Me%nIndexOutputs
            
                    call WriteDataLine(Cohort%CohortOutput%Unit(iIndexOutput), '<EndTimeSerie>')

                    call UnitsManager(Cohort%CohortOutput%Unit(iIndexOutput), CLOSE_FILE, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine CloseFiles - ModuleBivalve - ERR01'
                    
                enddo
                
                deallocate(Cohort%CohortOutput%Unit)
                
            end if
            Cohort  => Cohort%Next

        enddo 
                        
        if (Species%Population) then
             
            do iIndexOutput = 1, Me%nIndexOutputs

                call WriteDataLine(Species%PopulationOutput%Unit(iIndexOutput), '<EndTimeSerie>')
                
                call UnitsManager(Species%PopulationOutput%Unit(iIndexOutput), CLOSE_FILE, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine CloseFiles - ModuleBivalve - ERR02'
            
            enddo
            
            deallocate(Species%PopulationOutput%Unit)
            
        end if

        if (Species%BySizeOutput) then
        
            do iIndexOutput = 1, Me%nIndexOutputs
        
                call WriteDataLine(Species%SizeDistributionOutput%Unit(iIndexOutput), '<EndTimeSerie>')
                
                call UnitsManager(Species%SizeDistributionOutput%Unit(iIndexOutput), CLOSE_FILE, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine CloseFiles - ModuleBivalve - ERR03'
            
            enddo
            
            deallocate(Species%SizeDistributionOutput%Unit)

        end if
        
        if (Me%Testing_Parameters) then
        
            do iIndexOutput = 1, Me%nIndexOutputs

                call UnitsManager(Species%TestingParametersOutput%Unit(iIndexOutput), CLOSE_FILE, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine CloseFiles - ModuleBivalve - ERR03'
            
            enddo
            
            deallocate(Species%TestingParametersOutput%Unit)

!                call UnitsManager(Species%MakeRplotsPopulation%Unit, CLOSE_FILE, STAT = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine CloseFiles - ModuleBivalve - ERR04'
!
!                call UnitsManager(Species%MakeRplotsSizeDistribution%Unit, CLOSE_FILE, STAT = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine CloseFiles - ModuleBivalve - ERR05'

        end if


    end subroutine CloseFiles

    !-------------------------------------------------------------------------------

    subroutine DeallocateInstance ()

        !Arguments------------------------------------------------------------------

        !Local----------------------------------------------------------------------
        type (T_Bivalve), pointer          :: AuxObjBivalve
        type (T_Bivalve), pointer          :: PreviousObjBivalve
       
        !Begin----------------------------------------------------------------------

        !Updates pointers
        if (Me%InstanceID == FirstObjBivalve%InstanceID) then
            FirstObjBivalve => FirstObjBivalve%Next
        else
            PreviousObjBivalve => FirstObjBivalve
            AuxObjBivalve      => FirstObjBivalve%Next
            
            do while (AuxObjBivalve%InstanceID /= Me%InstanceID)
                PreviousObjBivalve => AuxObjBivalve
                AuxObjBivalve      => AuxObjBivalve%Next
            enddo

            !Now update linked list
            PreviousObjBivalve%Next => AuxObjBivalve%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 


    end subroutine DeallocateInstance

    !-------------------------------------------------------------------------------
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !-------------------------------------------------------------------------------

    subroutine Ready (ObjBivalve_ID, ready_) 

        !Arguments------------------------------------------------------------------
        integer         :: ObjBivalve_ID
        integer         :: ready_

        !---------------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjBivalve_ID > 0) then
            call LocateObjBivalve (ObjBivalve_ID)
            ready_ = VerifyReadLock (mBivalve_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !---------------------------------------------------------------------------

    end subroutine Ready

    !-------------------------------------------------------------------------------

    subroutine LocateObjBivalve (ObjBivalveID)

        !Arguments------------------------------------------------------------------
        integer         :: ObjBivalveID

        !Local----------------------------------------------------------------------

        Me => FirstObjBivalve
        do while (associated (Me))
            if (Me%InstanceID == ObjBivalveID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleBivalve - LocateObjBivalve - ERR01'

    end subroutine LocateObjBivalve

    !-------------------------------------------------------------------------------

    end module ModuleBivalve