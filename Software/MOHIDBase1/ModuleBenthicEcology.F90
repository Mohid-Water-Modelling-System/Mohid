
!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : BenthicEcology
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2012
! REVISION      : Isabella Ascione Kenov
! DESCRIPTION   : Module to compute simple benthic ecology processes
!
!
!
!
! LOG OF INPUTS
! 24-09-2012 - Generic DEB algorithm for bivalves growth (T Dabrowski code from EASYCO), by M Mateus
!             
!              ******************************************************************
!
!              Shellfish growth model based on DEB theory developed at the Irish 
!              Marine Institute as part of the EASYCO project funded under
!              the Interreg IVB Atlantic Area Programme.
!
!              References:
!	           Dabrowski, T., Lyons, K., Curé, M., Berry, A., Nolan, G., (Accepted)
!              Numerical modelling of spatio-temporal variability of growth of Mytilus edulis (L.)
!              and influence of its cultivation on ecosystem functioning. Journal of Sea Research.
!
!              ******************************************************************
!
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
!DataFile
! ==================MINERALIZATION ====================================
!OUTPUT_TIME              : 0  3600.
!DT_OUTPUT_TIME           : 3600.

!DT                : 120.
!PELAGIC_MODEL     : WaterQuality
!NCRATIO           : 0.18
!PCRATIO           : 0.024
!NITROGEN_USE      : 0
!CARBON_USE        : 1
!PHOSPHORUS_USE    : 0
!PARTICWATERFOOD   : 1
!<end_food>

!<begin_food>
!NAME              : particulate organic
!NCRATIO           : 0.18
!PCRATIO           : 0.024
!NITROGEN_USE      : 1
!CARBON_USE        : 0
!PHOSPHORUS_USE    : 1
!PARTICWATERFOOD   : 1
!<end_food>

!<end_consumer>
!! ========CONSUMERS: Deposit feeders========================================================
!<begin_consumer>
!NAME              : deposit feeders
!BETA20            : 0.002    ! mortality rate in 1/day
!GRMAX             : 0.05     ! maximum ingestion rate 1/day
!RESP20            : 0.0018   ! respiration rate in day^-1
!ASS_EFIC          : 0.25     ! assimilation efficiency (-)
!TFAC              : 1.08     ! dimensionless temperature decay factor
!NCRATIO           : 0.18     ! ratio N:C in benthic feeder
!PCRATIO           : 0.024    ! ratio P:C in benthic feeder
!TOPTMIN           : 13.      ! Minimum temperature of the optimal interval for benthic feeders growth ºC 
!TOPTMAX           : 28.      ! Maximum temperature of the optimal interval for benthic feeders growth ºC 
!TMIN              : 6.       ! Minimum temparature benthic feeders growth ºC 
!TMAX              : 37.      ! Maximum temparature benthic feeders growth ºC 
!K1                : 0.3      ! Constant to control temperature response curve shape 
!K2                : 0.98     ! Constant to control temperature response curve shape 
!K3                : 0.98     ! Constant to control temperature response curve shape 
!K4                : 0.02     ! Constant to control temperature response curve shape 
!COHESIVESED       : 0        ! is the grazing affected by cohesive sediment in water? 
!                               0 (for deposit feeders) , 1 for suspension feeders 
!KFOODC            : 0.001    !  half saturation constant for Carbon uptake by benthic feeder (Kg C/m2)   
!SMIN              : 0.002    ! Minimum consumer's biomass that limits the growth rate (Kg C/m2)
!SMAX              : 0.006    ! Maximum consumer's biomass that limits the growth rate (Kg C/m2)  

! food for deposit feeders:
!<begin_food>
!NAME              : microphytobenthos
!NCRATIO           : 0.18     ! 0.18     ratio N:C in the food
!PCRATIO           : 0.024    ! 0.024    ratio P:C in the food
!NITROGEN_USE      : 1        !          is the food expressed as N ? yes =1 No =0 
!CARBON_USE        : 1        !          is the food expressed as C ? yes =1 No =0 
!PHOSPHORUS_USE    : 1        !          is the food expressed as P ? yes =1 No =0 
!PARTICWATERFOOD   : 0        !          is the food particulate water property ? yes =1 No =0
!MINVAL            : 0.2e-3   ! 0.2e-3   minimum biomass for grazing over food
!<end_food>
!<end_consumer>


! Producers
!<begin_producer>
!NAME             : microphytobenthos
!DESCRIPTION      : microphytobenthos
!MORTALITY_RATE   : 0.02      ! 0.02     producer mortality rate 1/day
!RESPFRAC         : 0.05      ! 0.05     producer respiration fraction [-]
!VMAX             : 2.        ! 2.       producer maximum growth rate (1/day)
!ALPHA            : 0.025     ! 0.025    Slope of P/I curve 1/[day*(W/m2)] 
!KN               : 0.014e-3  ! 0.014e-3 Producer half saturation constant for N uptake g N/l
!KP               : 0.001e-3  ! 0.001e-3 Producer half saturation constant for P uptake g P/l
!NCRATIO          : 0.18      ! 0.18     ratio N:C in producer (gN/gC)
!PCRATIO          : 0.024     ! 0.024    ratio P:C in producer (gN/gC)
!TOPTFMIN         : 25.       ! 25       Minimum temperature of the optimal interval for microphytobenthos photosynthesis (ºC)
!TOPTFMAX         : 26.5      ! 26.5     Maximum temperature of the optimal interval for microphytobenthos photosynthesis (ºC)
!TFMIN            : 10.       ! 4        Minimum tolerable temperature for microphytobenthos photosynthesis (ºC)
!TFMAX            : 30.       ! 37       Maximum tolerable temperature for microphytobenthos photosynthesis (ºC)
!TFCONST1         : 0.05      ! 0.05     Constant to control temperature response curve shape on microphytobenthos [-]
!TFCONST2         : 0.98      ! 0.98     Constant to control temperature response curve shape on microphytobenthos [-]
!TFCONST3         : 0.98      ! 0.98     Constant to control temperature response curve shape on microphytobenthos [-]
!TFCONST4         : 0.02      ! 0.02     Constant to control temperature response curve shape on microphytobenthos [-]
!MINBIOMASS       : 0.0000001 ! 0.000001  Minimum producer biomass (kg C/m2) 
!EROCRIT          : 2.        ! 2.       Critical shear erosion for benthic producers bed (Pa)
!MINP             : 0.001     ! Minimum producer's biomass that limits the growth rate (Kg C/m2)
!MINP             : 0.005     ! Maximum producer's biomass that limits the growth rate (Kg C/m2)
!<end_producer>


! Seagrasses=========================================

! Seagrasses parameters
!GMAX              : 0.23      ! maximum growth rate 1/day
!MORTL0            : 0.064     ! Mortality rate 1/day
!MORTR0            : 0.035     ! Mortality rate 1/day
!KMAX              : 0.23      ! Maximum leaves density (kg DW/m2)
!NMAX              : 31.       ! Maximum internal nitrogen content (mgN/gdw)
!NMIN              : 5.        ! Minimum internal nitrogen content (mgN/gdw)
!NCRIT             : 16        ! Critical internal nitrogen content (mgN/gdw)
!PMAX              : 3.14      ! Maximum internal phosphorus content (mgP/gdw)
!PMIN              : 0.14      ! Minimum internal phosphorus content (mgP/gdw)
!PCRIT             : 0.8       ! Critical internal phosphorus content (mgP/gdw)
!KTR               : 0.28      ! Transfer coefficient from leaves o roots [%]
!MINLEAVES         : 0.000005  ! Minimum Leaves biomass (Kgdw/m2)
!MINROOTS          : 0.000005  ! Minimum Roots biomass (Kgdw/m2)
!EROCRIT           : 2.        ! Critical shear stress erosion (Pa) (tentative)
!GPKGDW            : 1.8       ! Ratio between gP/Kgdw in seagrasses 
!GNKGDW            : 16.       ! Ratio between gN/Kgdw in seagrasses 
!GCKGDW            : 330.      ! Ratio between gC/Kgdw in seagrasses 
!TOPTMIN           : 13.       ! Minimum temperature of the optimal interval for seagrasses growth ºC 
!TOPTMAX           : 28.       ! Maximum temperature of the optimal interval for seagrasses growth ºC 
!TMIN              : 6.        ! Minimum temparature Seagrasses growth ºC 
!TMAX              : 37.       ! Maximum temparature Seagrasses growth ºC 
!K1                : 0.3       ! Constant to control temperature response curve shape 
!K2                : 0.98      ! Constant to control temperature response curve shape 
!K3                : 0.98      ! Constant to control temperature response curve shape 
!K4                : 0.02      ! Constant to control temperature response curve shape 
!LAT               : 40.5      ! Average latitude of the geographic location (used to calculate the daylight)
!MORT_TYPE         : 1/2       ! mortality term as a function of temperature (1) or daylight (2)


! T Dabrowski algorithm (EASYCO) for bivalves
!(values and description provided by the author)
!
!NINDM2         100.                              ! Number of individuals per m2
!MAXEQUE        2190000000.0                      ! Max. equilibrium energy density in J m-3
!EGCC           1900000000.0                      ! Energetic growth cost per unit growth in structural body volume in J m-3
!MAXSAING       22.7778                           ! Max. surface-area-specific ingestion rate in J m-2 s -1
!ASSEFFIC       0.75                              ! Assimilation efficiency (= Pam/Pxm)
!SOMACOST       277.7778                          ! Somatic maintenance cost in J m -3 s -1
!FLUXRESFRAC    0.70                              ! Fraction of flux from reserve spent on somaticmaintenance
!SATCOEF        1.77                              ! Saturation coefficient - food density at which ingestion rate is half the maximum in mg (chl-a) m -3
!REPOREFFIC     0.9                               ! Reproduction efficiency
!ENERCONT       17550.0                           ! Energy content of 1g of reserve
!WWTODW         0.20                              ! WW to DW converter
!SPAWEFFIC      0.90                              ! proportion of the reproductive buffer emptied at each spawning
!TSPAWN         2.0                               ! temperature threshold trigeering spawning in degC
!GSOMINDEX      0.4                               ! Gonado-somatic index triggering spawning (fraction)                          
!JDAYSPAWN      180                               ! Julian day for spawning, NOTE it is an array, so follow namelist formats
!NSPAWND        1                                 ! Number of spawning days in a year (MUST be same as SpawnDay array dimension)
!VOLADULT       0.00000006                        ! volume specifying a change from juvenile to adult in m 3 
!SHAPEPARAM     0.287                             ! shape parameter
!REFTEMP        293.0                             ! Ref. temp for rate constants in K
!TEMPARRH       5800.0                            ! Arrhenius temperature in K
!TEMPLOW        275.0                             ! Lower boundary of tolerance range in K
!THRESP         296.0                             ! Upper boundary of tolerance range for respiration in K
!THING          296.0                             ! Upper boundary of tolerance range for ingestion in K
!TAL            45430.0                           ! Arrhenius temp. for rate of decrease at lower boundary in K
!TAHRESP        31376.0                           ! Arrhenius temp. for rate of decrease at upper boundary for respiration in K
!TAHING         31376.0                           ! Arrhenius temp. for rate of decrease at upper boundary for ingestion in K
!FECALDECAY     2.1222e-4                         ! Corresponds to T90 of 3hrs assumed after document 'MOHID modules', in s-1
!ENERGPHY       47.7546                           ! energetic value of phyt C in J mg-1 C (Platt and Irwin, 1973)
!ENERGO2        14.3                              ! energetic value of oxygen (=14.3 J mg-1 O2)

! Accessory keywords (not in the original algorithm)
! These Keywords are used for enabling teh original code to comply with the MOHID code architecture

!CCHLABIV       60.0                              ! C to Chla ratio used in the bivalve DEB model to converto Phyto to Chla
!NH4EXCFRAC     0                                 ! To define NH4 excretion as a fraction of N ingested (defiend by the name of the bivalve in the original)
!BILENINIC      0.05                              ! Initial bivalve lenght in m





Module ModuleBenthicEcology

    use ModuleGlobalData
    use ModuleEnterData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    
    public  :: ConstructBenthicEcology
    private ::      AllocateInstance
    private ::      ReadData
    private ::          ConstructGlobalVariables
    private ::          ReadOrganicMatterParameters
    private ::          ReadOxygenParameters
    private ::          ReadNitrogenParameters
    private ::          ReadPhosphorusParameters
    private ::          ReadSilicaParameters
    private ::          ReadSeagrassesParameters
    private ::          ConstructProducers
    private ::              AddProducer
    private ::              ConstructProducerParameters
    private ::          ConstructConsumers
    private ::              AddConsumer
    private ::              ConstructConsumerParameters
    private ::          ConstructBivalveDEB
    private ::              AddBivalveDEB
    private ::              ConstructBivalveDEBParameters
   ! private ::          ConstructDecomposers
   ! private ::              AddDecomposer
   ! private ::              ConstructDecomposerParameters
    private ::              ConstructGrazing
    private ::                  AddFood
    private ::                  ConstructFood
    private ::          PropertyIndexNumber
    private ::          ConstructPropertyList
        
        
    !Selector
    public  :: GetDTBenthicEcology
    public  :: GetBenthicEcologyPropertyList
    public  :: GetBenthicEcologySize
    public  :: GetBenthicEcologyPropIndex
    public  :: UnGetBenthicEcology
    public  :: GetBenthicEcologyRateFlux  
    public  :: UnGetBenthicEcologyRateFlux  
                     
    
    !Modifier
    public  :: ModifyBenthicEcology
    private ::      ComputeBenthicProducers
    private ::      ComputeBenthicConsumers
    private ::      ComputeBenthicBivalveDEB   
    private ::      ComputeBenthicPhyto
    private ::      ComputeSeagrasses
    private ::      ComputeBenthicSilica
    private ::      ComputeBenthicPhosphorus
    private ::      ComputeBenthicNitrogen
    !private ::      Decomposers
   
    
    !Destructor
    public  :: KillBenthicEcology                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjBenthicEcology 
    
    
    
  private :: T_PropIndex
    type       T_PropIndex
         
        integer                                 :: Ammonia      = null_int        
        integer                                 :: Nitrate      = null_int
        integer                                 :: POC          = null_int                 
        integer                                 :: PON          = null_int 
        integer                                 :: POP          = null_int
        integer                                 :: Oxygen       = null_int       
        integer                                 :: Phyto        = null_int
        integer                                 :: Phosphate    = null_int
        integer                                 :: Bacteria     = null_int
        integer                                 :: PONr         = null_int
        integer                                 :: DONnr        = null_int
        integer                                 :: Diatoms      = null_int
        integer                                 :: PON1         = null_int
        integer                                 :: PON2         = null_int
        integer                                 :: PON3         = null_int
        integer                                 :: PON4         = null_int
        integer                                 :: PON5         = null_int
        integer                                 :: POP1         = null_int
        integer                                 :: POP2         = null_int
        integer                                 :: POP3         = null_int
        integer                                 :: POP4         = null_int
        integer                                 :: POP5         = null_int
        integer                                 :: DissolvedSilica   = null_int
        integer                                 :: BioSilica         = null_int
        integer                                 :: Leaves            = null_int 
        integer                                 :: Roots             = null_int 
        integer                                 :: Nint              = null_int 
        integer                                 :: Pint              = null_int       
    end type T_PropIndex

     type T_StoredIndex
            integer                             :: NintFactor      = 1
            integer                             :: PintFactor      = 2
            integer                             :: NintFactorR     = 3
            integer                             :: RootsMort       = 4
            integer                             :: PintFactorR     = 5
       end type T_StoredIndex

       type      T_Leaves
            real                                :: NC_Ratio         = null_real 
            real                                :: PC_Ratio         = null_real
            real                                :: MortalityRate    = null_real
            real                                :: MaxGrowthRate    = null_real
            real                                :: Ktr              = null_real
            real                                :: ErosCritShear    = null_real
            real                                :: MinimumBiomass   = null_real
            real                                :: gNKgDW           = null_real
            real                                :: gPKgDW           = null_real
            real                                :: gCKgDW           = null_real
            real                                :: K1               = null_real
            real                                :: K2               = null_real
            real                                :: K3               = null_real
            real                                :: K4               = null_real
            real                                :: TOptMin          = null_real
            real                                :: TOptMax          = null_real
            real                                :: TMin             = null_real
            real                                :: TMax             = null_real
            real                                :: Lat              = null_real
            real                                :: KMAX             = null_real
            integer                             :: MortType         = null_int
       end type T_Leaves

      type      T_Roots
            real                                :: MortalityRate    = null_real
            real                                :: MinimumBiomass   = null_real
       end type T_Roots
       
       type T_Nint
           real                                 :: Nmax              = null_real
           real                                 :: Nmin              = null_real 
           real                                 :: Ncrit             = null_real  
       end Type T_Nint
        
        
           type T_Pint
           real                                 :: Pmax              = null_real
           real                                 :: Pmin              = null_real 
           real                                 :: Pcrit             = null_real  
       end Type T_Pint
    
    
     private :: T_External
    type       T_External
        real, pointer, dimension(:  )       :: Temperature
        real, pointer, dimension(:,:)       :: MassInKgFromWater
        real, pointer, dimension(:  )       :: Sediment
        real(8), pointer, dimension(:  )       :: WaterVolume
        real, pointer, dimension(:  )       :: CellArea
        real, pointer, dimension(:  )       :: ShortWaveTop
        real, pointer, dimension(:  )       :: ShortWaveAverage
        real, pointer, dimension(:  )       :: LightExtCoefField
        real, pointer, dimension(:  )       :: Thickness
        real, pointer, dimension(:  )       :: ShearStress
        real, pointer, dimension(:,:)       :: Mass
        real, pointer, dimension(:  )       :: UptakeNH4NO3w
        real, pointer, dimension(:  )       :: UptakeNH4s
        real, pointer, dimension(:  )       :: UptakePO4w
        real, pointer, dimension(:  )       :: UptakePO4s
        real, pointer, dimension(:  )       :: LightFactor   
    end type T_External
  
  
  type     T_ComputeOptions
        logical                             :: Nitrogen                      = .false.
        logical                             :: Phosphorus                    = .false.
        logical                             :: Producers                     = .false.     
        logical                             :: Consumers                     = .false. 
        logical                             :: Silica                        = .false.
        logical                             :: Diatoms                       = .false.  ! CONTROLLARE
        logical                             :: Phyto                         = .false.  ! ! CONTROLLARE
        logical                             :: PhytoFilt                     = .false.  ! ! CONTROLLARE
        logical                             :: Pompools                      = .false.
        logical                             :: Seagrasses                    = .false.
  end type     T_ComputeOptions
  
    type     T_OrganicMatter
        real                                        :: NC_Ratio         = null_real 
        real                                        :: PC_Ratio         = null_real 
    end type T_OrganicMatter

    type     T_Nitrogen
        real                                        :: PONDecayRate     = null_real 
        real                                        :: PONDecayTFactor  = null_real
    end type T_Nitrogen
    
    type     T_Phosphorus
        real                                        :: POPDecayRate     = null_real 
        real                                        :: POPDecayTFactor  = null_real
    end type T_Phosphorus
    
    type     T_Oxygen
        real                                        :: Minimum          = null_real 
    end type T_Oxygen

    type     T_Silica
        real                                        :: BioSiDecayRate   = null_real 
    end type T_Silica

    type     T_Phyto
        real                                        :: NC_Ratio         = null_real 
        real                                        :: PC_Ratio         = null_real 
        real                                        :: MortalityRate    = null_real 
    end type T_Phyto
    
   private :: T_BenthicEcology
    type       T_BenthicEcology
        integer                                      :: InstanceID
        type (T_Size1D)                              :: Prop
        type (T_Size1D)                              :: Size
        real,    dimension(:,:,:),  pointer          :: Matrix
        real                                         :: DT, DTDay
        character(len=StringLength)                  :: PelagicModel
        integer, dimension(:), pointer               :: PropertyList
        real,   dimension(:,:),    pointer           :: StoredArray
        type(T_Size1D       )                        :: Array
        type(T_PropIndex    )                        :: PropIndex
        type(T_Producer     ), pointer               :: FirstProducer
        type(T_Consumer     ), pointer               :: FirstConsumer
        type(T_BivalveDEB   ), pointer               :: FirstBivalveDEB
        type(T_External     )                        :: ExternalVar
        integer                                      :: ObjEnterData = 0
        type(T_BenthicEcology), pointer              :: Next
        type(T_ComputeOptions)                       :: ComputeOptions
        type(T_Silica        )                       :: Silica
        type(T_Oxygen        )                       :: Oxygen
        type(T_Phyto         )                       :: Phyto
        type(T_OrganicMatter )                       :: OrganicMatter
        type(T_Nitrogen      )                       :: Nitrogen
        type(T_Phosphorus    )                       :: Phosphorus
        type(T_Leaves        )                       :: Leaves
        type(T_Roots         )                       :: Roots
        type(T_Nint          )                       :: Nint
        type(T_Pint          )                       :: Pint
        type(T_StoredIndex   )                       :: StoredIndex
        integer                                      :: JulianDay
    end type  T_BenthicEcology
  
  
  
  
  private :: T_ID
        type       T_ID
            integer                         :: ID
            character(len=StringLength)     :: Name
            character(len=StringLength)     :: Description
        end type   T_ID
        
  
   
  


    private :: T_PoolIndex
        type       T_PoolIndex
            integer                                 :: Carbon      = null_int
            integer                                 :: Nitrogen    = null_int         
            integer                                 :: Phosphorus  = null_int         
            integer                                 :: Silica      = null_int
            integer                                 :: Chlorophyll = null_int        
        end type T_PoolIndex
   
   
   private :: T_Producer
   type     T_Producer
   type(T_ID)                                       :: ID
   type(T_PoolIndex)                                :: PoolIndex
        real                                        :: RespirationRate  = null_real   
        real                                        :: MortalityRate    = null_real 
        real                                        :: Vmax             = null_real  
        real                                        :: RespFrac         = null_real
        real                                        :: alpha            = null_real
        real                                        :: KN               = null_real 
        real                                        :: KP               = null_real 
        real                                        :: KTRANS           = null_real 
        real                                        :: KLIGHT           = null_real
        real                                        :: NCratio          = null_real   
        real                                        :: PCratio          = null_real
        real                                        :: TOptPhytoMin     = null_real
        real                                        :: TOptPhytoMax     = null_real
        real                                        :: TPhytoMin        = null_real
        real                                        :: TPhytoMax        = null_real
        real                                        :: FK1              = null_real
        real                                        :: FK2              = null_real
        real                                        :: FK3              = null_real
        real                                        :: FK4              = null_real
        real                                        :: MINP             = null_real
        real                                        :: MAXP             = null_real
        real                                        :: EroCritShear     = null_real
        real                                        :: MinimumBiomass   = null_real
   type(T_Producer    ), pointer                    :: Next 
    end type T_Producer
 
 private :: T_Food
        type   T_Food
            type(T_ID)                  :: ID
            type(T_Food), pointer       :: Next
            real                        :: NCratio
            real                        :: PCratio
            real                        :: Minval
            logical                     :: Use_Carbon
            logical                     :: Use_Nitrogen
            logical                     :: Use_Phosphorus
            logical                     :: ParticulateWaterFood
        end type T_Food


    private :: T_Grazing
        type   T_Grazing
            real                        :: Ass_Efic         = FillValueReal !assimilation efficiency
            real                        :: NCRatio          = FillValueReal
            logical                     :: Cohesivesed                   ! 1 effect of cohesivesed on grazing
            real                        :: MaxSediment      = null_real 
            real                        :: KFoodC           = FillValueReal
            type(T_Food), pointer       :: FirstFood
        end type T_Grazing


! Bivalve DEB model 


    private :: T_BiSV
        type       T_BiSV
            integer                                 :: Bivol       = null_int
            integer                                 :: Bivenerg    = null_int         
            integer                                 :: Bivenerep   = null_int
            integer                                 :: Bivinddw    = null_int         
            integer                                 :: Bivtotaldw  = null_int
                        
        end type T_BiSV

   private :: T_BivalveDEB 
   type     T_BivalveDEB
        type(T_ID)                                  :: ID
        type(T_BiSV)                                :: BiSv
        type(T_BivalveDEB    ), pointer             :: Next
                                              
        real                                        :: Em               = null_real  ! Max. equilibrium energy density in J m-3
        real                                        :: Eg               = null_real  ! Energetic growth cost per unit growth in structural
                                                                                           ! body volume in J m-3
        real                                        :: Pxm              = null_real  ! Max. surface-area-specific ingestion rate
                                                                                           ! in J m-2 s -1
        real                                        :: ka               = null_real  ! Assimilation efficiency (= Pam/Pxm)
        real                                        :: Pm               = null_real  ! Somatic maintenance cost in J m -3 s -1
        real                                        :: kappa            = null_real  ! Fraction of flux from reserve spent on somatic
                                                                                            ! maintenance
        real                                        :: Xk               = null_real  ! Saturation coefficient - food density at which
                                                                                           ! ingestion rate is half the maximum
                                                                                           ! in mg (chl-a) m -3
        real                                        :: kR               = null_real  ! Reproduction efficiency
        real                                        :: Conv_fac         = null_real  ! Energy content of 1g of reserve
        real                                        :: St_DW_perc       = null_real  ! WW to DW converter
        real                                        :: spawn_eff        = null_real  ! proportion of the reproductive buffer emptied at each spawning
        real                                        :: Tspawn           = null_real  ! temperature threshold triggering spawning
        real                                        :: RGS              = null_real  ! Gonado-somatic index triggering spawning (fraction)
 
        integer                                     :: spawnD                        ! Number of spawning days in a year (MUST be same as SpawnDay array dimension)
        integer, dimension(:), pointer              :: spawnDay                      ! Julian Day(s) to spawn (if no spawning then set spawnD to zero)
        integer                                     :: NInd                          ! Number of individuals per m2
        
        real                                        :: volp             = null_real  ! volume specifying a change from juvenile to adult
                                                                                           ! in m 3
        real                                        :: delm             = null_real  ! shape parameter
        real                                        :: Tref             = null_real  ! Ref. temp for rate constants in K
        real                                        :: Ta               = null_real  ! Arrhenius temperature in K
        real                                        :: TL               = null_real  ! Lower boundary of tolerance range in K
        real                                        :: THresp           = null_real  ! Upper boundary of tolerance range for respiration in K
        real                                        :: THing            = null_real  ! Upper boundary of tolerance range for ingestion in K
        real                                        :: TAL              = null_real  ! Arrhenius temp. for rate of decrease at lower boundary in K
        real                                        :: TAHresp          = null_real  ! Arrhenius temp. for rate of decrease at upper boundary for respiration in K
        real                                        :: TAHing           = null_real  ! Arrhenius temp. for rate of decrease at upper boundary for ingestion in K
        
        real                                        :: kd               = null_real  ! Corresponds to T90 of 3hrs assumed after document 'MOHID modules', in s-1
        
        real                                        :: Ec               = null_real  ! energetic value of phyt C in J mg-1 C (Platt and Irwin, 1973)
        
        real                                        :: etaO2            = null_real  ! energetic value of oxygen (=14.3 J mg-1 O2)
 
        !  Accessory variables (not in the original algorithm)
        real                                        :: CtoChla          = null_real  ! C to Chla ratio used in the bivalve DEB model to converto Phyto to Chla 
        integer                                     :: NH4Frac                       ! To define NH4 excretion as a fraction of N ingested
        real                                        :: BiLen1           = null_real  ! Initial bivalve lenght in m
 
 end type T_BivalveDEB
 
 
 
 private :: T_Consumer 
   type     T_Consumer
        type(T_ID)                                       :: ID
        type(T_Food), pointer                            :: FirstFood 
        type(T_PoolIndex)                                :: PoolIndex
        type(T_Consumer    ), pointer                    :: Next 
        type(T_Grazing     )                             :: Grazing
        real                                             :: RespirationRate   = null_real  
        real                                             :: MortalityRate     = null_real
        real                                             :: BETA20            = null_real
        real                                             :: TemperatureFactor = null_real
        real                                             :: GRMAX             = null_real ! l/day/gC
        real                                             :: RESP20            = null_real  ! [1/day]
        real                                             :: NCratio           = null_real
        real                                             :: PCratio           = null_real    
        real                                             :: TOptMin           = null_real
        real                                             :: TOptMax           = null_real
        real                                             :: TMin              = null_real
        real                                             :: TMax              = null_real
        real                                             :: K1                = null_real
        real                                             :: K2                = null_real
        real                                             :: K3                = null_real
        real                                             :: K4                = null_real
        real                                             :: SMIN             = null_real
        real                                             :: SMAX             = null_real
        real                                             :: KO2              = null_real
 end type T_Consumer

    
    
    
    

!Global Module Variables
    type (T_BenthicEcology), pointer                          :: FirstObjBenthicEcology
    type (T_BenthicEcology), pointer                          :: Me

  !--------------------------------------------------------------------------
    
     contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 
   subroutine ConstructBenthicEcology(ObjBenthicEcologyID, FileName, ILB,IUB, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjBenthicEcologyID  
        character(len=StringLength)                     :: FileName
        integer, optional, intent(OUT)                  :: STAT
        integer          , intent(IN )                  :: ILB, IUB

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mBenthicEcology_)) then
            nullify (FirstObjBenthicEcology)
            call RegisterModule (mBenthicEcology_) 
        endif
        
        call Ready(ObjBenthicEcologyID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%Size%ILB = ILB
            Me%Size%IUB = IUB
            
            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBenthicEcology - ModuleBenthicEcology - ERROR #1'

            call ReadData

            call PropertyIndexNumber
        
            call ConstructPropertyList
            
            call ConstructMatrix
            
            call ConstructRates

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBenthicEcology - ModuleBenthicEcology - ERROR #2'

            !Returns ID
            ObjBenthicEcologyID          = Me%InstanceID

            STAT_ = SUCCESS_
        else 
            
            stop 'ModuleBenthicEcology - ConstructBenthicEcology - ERROR #1' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructBenthicEcology
   
   
   
    !----------------------------------------------------------------------
    
        subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_BenthicEcology), pointer                         :: NewObjBenthicEcology
        type (T_BenthicEcology), pointer                         :: PreviousObjBenthicEcology


        !Allocates new instance
        allocate (NewObjBenthicEcology)
        nullify  (NewObjBenthicEcology%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjBenthicEcology)) then
            FirstObjBenthicEcology         => NewObjBenthicEcology
            Me                   => NewObjBenthicEcology
        else
            PreviousObjBenthicEcology      => FirstObjBenthicEcology
            Me                   => FirstObjBenthicEcology%Next
            do while (associated(Me))
                PreviousObjBenthicEcology  => Me
                Me               => Me%Next
            enddo
            Me                   => NewObjBenthicEcology
            PreviousObjBenthicEcology%Next => NewObjBenthicEcology
        endif

        Me%InstanceID = RegisterNewInstance (mBenthicEcology_)

    end subroutine AllocateInstance


!____________________________________________________________________________________________________
!____________________________________________________________________


subroutine ReadData

        !Arguments-------------------------------------------------------------
                                                           
        !Local-----------------------------------------------------------------
        
        call ConstructGlobalVariables
        
        call ReadOrganicMatterParameters
        call ReadOxygenParameters
        call ReadNitrogenParameters
        call ReadPhosphorusParameters
        call ReadSilicaParameters
        call ReadPhytoParameters
        call ReadSeagrassesParameters
        call ConstructProducers

        call ConstructConsumers
        
        call ConstructBivalveDEB
                
       ! call ConstructDecomposers

    end subroutine ReadData
    
    
    
       subroutine PropertyIndexNumber

        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Consumer),      pointer             :: Consumer
        type(T_BivalveDEB),    pointer             :: BivalveDEB
       ! type(T_Decomposer),    pointer             :: Decomposer
        integer                                    :: Index
        !Local-----------------------------------------------------------------
        
        Me%Prop%ILB = 1
        Me%Prop%IUB = 0

        Index               = 0
        


        !Producer index number      
            Producer => Me%FirstProducer
            do while(associated(Producer))
                
                Index                               = Index + 1
                Producer%PoolIndex%Carbon           = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                Producer%PoolIndex%Nitrogen         = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                Producer%PoolIndex%Phosphorus       = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                Producer => Producer%Next
            end do
        
        !Consumer index number
            Consumer => Me%FirstConsumer
            do while(associated(Consumer))
                
                Index                               = Index + 1
                Consumer%PoolIndex%Carbon           = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                Consumer%PoolIndex%Nitrogen         = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                                
                Index                               = Index + 1
                Consumer%PoolIndex%Phosphorus       = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1

                Consumer => Consumer%Next
            end do
   
   
           !Bivalve (DEB) index number
            BivalveDEB => Me%FirstBivalveDEB
            do while(associated(BivalveDEB))
                
                Index                               = Index + 1
                BivalveDEB%BiSV%Bivol               = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                BivalveDEB%BiSV%Bivenerg            = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                
                Index                               = Index + 1
                BivalveDEB%BiSV%Bivenerep           = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                BivalveDEB%BiSV%Bivinddw            = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                BivalveDEB%BiSV%Bivtotaldw          = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1

                BivalveDEB => BivalveDEB%Next
            end do
   
       
      
   
   
        !Decomposer index number
          !  Decomposer => Me%FirstDecomposer
          !  do while(associated(Decomposer))
                
            !    Index                               = Index + 1
             !    Decomposer%PoolIndex%Carbon         = Index
              !   Me%Prop%IUB                         = Me%Prop%IUB + 1


            !    Decomposer => Decomposer%Next
           ! end do
           
           
           
        
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen             = Me%Prop%IUB
        
        if(Me%ComputeOptions%Nitrogen)then

        !Nitrogen index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia        = Me%Prop%IUB
            
            if(Me%ComputeOptions%Producers)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Nitrate        = Me%Prop%IUB
            endif
        
        if (Me%PelagicModel == LifeModel) then
                !Particulate organic carbon index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POC            = Me%Prop%IUB
        endif

        !OrganicMatter index number

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%PON            = Me%Prop%IUB
      
                  if(Me%ComputeOptions%Pompools)then
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON1           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON2           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON3           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON4           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON5           = Me%Prop%IUB
                       
            end if   
   end if   
      
      !end if

        if(Me%ComputeOptions%Phosphorus)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phosphate      = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POP            = Me%Prop%IUB
            
            if(Me%ComputeOptions%Pompools)then
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP1           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP2           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP3           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP4           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP5           = Me%Prop%IUB
                       
            end if

        end if

        if(Me%ComputeOptions%Silica)then

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DissolvedSilica= Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%BioSilica      = Me%Prop%IUB

        end if


        if(Me%ComputeOptions%PhytoFilt .or. Me%ComputeOptions%Phyto)then
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phyto          = Me%Prop%IUB
        endif
        
        if(Me%ComputeOptions%Seagrasses)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Leaves         = Me%Prop%IUB
           
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Roots          = Me%Prop%IUB
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Nint           = Me%Prop%IUB
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Pint           = Me%Prop%IUB

        end if
        !----------------------------------------------------------------------

    end subroutine PropertyIndexNumber

    !--------------------------------------------------------------------------
    
    
    subroutine ConstructPropertyList

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Consumer),      pointer             :: Consumer
        type(T_BivalveDEB),    pointer             :: BivalveDEB
        !type(T_Decomposer),    pointer             :: Decomposer
        integer                                    :: Index
        !Local-----------------------------------------------------------------
        
        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))
        
        Index = 0
        
         Me%PropertyList(Me%PropIndex%Oxygen)                = Oxygen_

        if(Me%ComputeOptions%Nitrogen)then
            Me%PropertyList(Me%PropIndex%Ammonia)           = Ammonia_
            
                if(Me%ComputeOptions%Producers)then
                Me%PropertyList(Me%PropIndex%Nitrate)           = Nitrate_
                endif
                
            Me%PropertyList(Me%PropIndex%PON)               = PON_
            
            if(Me%ComputeOptions%Pompools)then
                Me%PropertyList(Me%PropIndex%PON1)               = PON1_
                Me%PropertyList(Me%PropIndex%PON2)               = PON2_
                Me%PropertyList(Me%PropIndex%PON3)               = PON3_
                Me%PropertyList(Me%PropIndex%PON4)               = PON4_
                Me%PropertyList(Me%PropIndex%PON5)               = PON5_
            end if
        end if

        if(Me%ComputeOptions%Phosphorus)then
            Me%PropertyList(Me%PropIndex%Phosphate)         = Inorganic_Phosphorus_
            Me%PropertyList(Me%PropIndex%POP)               = POP_
            
            if(Me%ComputeOptions%Pompools)then
                Me%PropertyList(Me%PropIndex%POP1)               = POP1_
                Me%PropertyList(Me%PropIndex%POP2)               = POP2_
                Me%PropertyList(Me%PropIndex%POP3)               = POP3_
                Me%PropertyList(Me%PropIndex%POP4)               = POP4_
                Me%PropertyList(Me%PropIndex%POP5)               = POP5_
            end if
        end if


       if(Me%ComputeOptions%Silica)then

            if(Me%PelagicModel .eq. WaterQualityModel) then
                Me%PropertyList(Me%PropIndex%DissolvedSilica)   = DSilica_
            endif

            if(Me%PelagicModel .eq. LifeModel) then            
                Me%PropertyList(Me%PropIndex%DissolvedSilica)   = Silicate_
            endif

            Me%PropertyList(Me%PropIndex%BioSilica)         = BioSilica_
        
        end if

        if(Me%ComputeOptions%Phyto .or. Me%ComputeOptions%PhytoFilt)then
            Me%PropertyList(Me%PropIndex%Phyto)             = Phytoplankton_
        end if
        
        
       if(Me%ComputeOptions%Seagrasses)then
            Me%PropertyList(Me%PropIndex%Leaves)           = SeagrassesLeaves_
            Me%PropertyList(Me%PropIndex%Roots)            = SeagrassesRoots_
            Me%PropertyList(Me%PropIndex%Nint)             = SeagrassesN_
            Me%PropertyList(Me%PropIndex%Pint)             = SeagrassesP_
        end if

        !Producer index number      
            Producer => Me%FirstProducer
            do while(associated(Producer))
                

               
                Me%PropertyList(Producer%PoolIndex%Carbon)     = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" carbon")
                Me%PropertyList(Producer%PoolIndex%Nitrogen)   = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" nitrogen")
                Me%PropertyList(Producer%PoolIndex%Phosphorus) = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" phosphorus")


                Producer => Producer%Next
            end do
        
        !Consumer index number
            Consumer => Me%FirstConsumer
            do while(associated(Consumer))
                

                Me%PropertyList(Consumer%PoolIndex%Carbon)     = &
                            GetPropertyIDNumber(trim(Consumer%ID%name)//" carbon")
                Me%PropertyList(Consumer%PoolIndex%Nitrogen)   = &
                            GetPropertyIDNumber(trim(Consumer%ID%name)//" nitrogen")
                Me%PropertyList(Consumer%PoolIndex%Phosphorus) = &
                            GetPropertyIDNumber(trim(Consumer%ID%name)//" phosphorus")


                Consumer => Consumer%Next
            end do
            
            
         !Bivalve index number
            BivalveDEB => Me%FirstBivalveDEB
            do while(associated(BivalveDEB))
                

                Me%PropertyList(BivalveDEB%BiSV%Bivol)      = &
                            GetPropertyIDNumber(trim(BivalveDEB%ID%name)//" volume")
                Me%PropertyList(BivalveDEB%BiSV%Bivenerg)   = &
                            GetPropertyIDNumber(trim(BivalveDEB%ID%name)//" energy")
                Me%PropertyList(BivalveDEB%BiSV%Bivenerep)  = &
                            GetPropertyIDNumber(trim(BivalveDEB%ID%name)//" reproduct energy")
                Me%PropertyList(BivalveDEB%BiSV%Bivinddw)  = &
                            GetPropertyIDNumber(trim(BivalveDEB%ID%name)//" individual dw")
                Me%PropertyList(BivalveDEB%BiSV%Bivtotaldw)  = &
                            GetPropertyIDNumber(trim(BivalveDEB%ID%name)//" total dw")

               BivalveDEB => BivalveDEB%Next
            end do
   
   
        !Decomposer index number
          !  Decomposer => Me%FirstDecomposer
          !  do while(associated(Decomposer))
                
             !   Me%PropertyList(Decomposer%PoolIndex%Carbon)     = &
                            !GetPropertyIDNumber(trim(Decomposer%ID%name)//" carbon")
               ! Me%PropertyList(Decomposer%PoolIndex%Nitrogen)   = &
                            !GetPropertyIDNumber(trim(Decomposer%ID%name)//" nitrogen")
               ! Me%PropertyList(Decomposer%PoolIndex%Phosphorus) = &
                            !GetPropertyIDNumber(trim(Decomposer%ID%name)//" phosphorus")

               ! Decomposer => Decomposer%Next
           ! end do

        !Nitrogen index number



    
    end subroutine ConstructPropertyList
    
    !----------------------------------------------------------------------
    subroutine ConstructRates

        !Begin-----------------------------------------------------------------

        allocate(Me%Matrix(Me%Size%ILB:Me%Size%IUB,                           &
                             Me%Prop%ILB:Me%Prop%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB))
                             
                             

      !NintFactor      = 1
      !PintFactor      = 2
      !NintfactorR     = 3
      !RootsMort       = 4
      !PintfactorR     = 5

    
    allocate(Me%StoredArray(Me%Size%ILB:Me%Size%IUB,1:5)) 
    
    end subroutine ConstructRates


    !----------------------------------------------------------------------
        subroutine ConstructMatrix

        !Begin-----------------------------------------------------------------

        allocate(Me%Matrix(Me%Size%ILB:Me%Size%IUB,                           &
                             Me%Prop%ILB:Me%Prop%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB))

    
    end subroutine ConstructMatrix
     !----------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromFile 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)
        
        
         call GetData(Me%DT,                                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DT',                                               &
                     Default      = 3600.,                                              &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #1'

        Me%DTDay = Me%DT / (3600. * 24.)
        
 
        call GetData(Me%PelagicModel,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PELAGIC_MODEL',                                    &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #9'
        if(iflag==0)then
            write(*,*)'Please define the pelagic model to couple with ModuleBenthicEcology'
            stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #9'
        end if

        if((Me%PelagicModel .ne. WaterQualityModel .and. Me%PelagicModel .ne. LifeModel))then
            write(*,*)'Pelagic model to couple with ModuleBenthicEcology must be one of the following:'
            write(*,*)trim(WaterQualityModel)
            write(*,*)trim(LifeModel)
            stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #10'
        endif
                 
        call GetData(Me%ComputeOptions%Nitrogen,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NITROGEN',                                         &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #11'

        call GetData(Me%ComputeOptions%Phosphorus,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHOSPHORUS',                                       &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR50'


        call GetData(Me%ComputeOptions%Silica,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SILICA',                                           &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR60'
        
         call GetData(Me%ComputeOptions%Phyto,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO',                                            &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR90'
        
                 call GetData(Me%ComputeOptions%PhytoFilt,                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTOFILT',                                        &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR100'
        
        if (Me%ComputeOptions%PhytoFilt .AND. Me%ComputeOptions%Phyto) &
        stop ' you can choose only one between PHYTO =1 and PHYTOFILT =1'
        
          call GetData(Me%ComputeOptions%Seagrasses,                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SEAGRASSES',                                       &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR110'            

  call GetData(Me%ComputeOptions%Producers,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PRODUCERS',                                        &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR115'
    end subroutine ConstructGlobalVariables



!----------------------------------------------------------------------

    subroutine ReadOrganicMatterParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%OrganicMatter%NC_Ratio,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NC_RATIO',                                         &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOrganicMatterParameters - ModuleBenthos - ERR01'

        call GetData(Me%OrganicMatter%PC_Ratio,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PC_RATIO',                                         &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOrganicMatterParameters - ModuleBenthos - ERR10'

    end subroutine ReadOrganicMatterParameters

    !--------------------------------------------------------------------------
    
     

    subroutine ReadOxygenParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Oxygen%Minimum,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MIN_OXYGEN',                                       &
                     Default      = 1e-5,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOxygenParameters - ModuleBenthos - ERR01'


    end subroutine ReadOxygenParameters
    
    !--------------------------------------------------------------------------
    
       subroutine ReadSeagrassesParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
        ! Maximum growth rate 1/day
        call GetData(Me%Leaves%MaxGrowthRate,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GMAX',                                             &
                     Default      = 0.23,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR01'
 
       ! Leaves Mortality rate 1/day
        call GetData(Me%Leaves%MortalityRate,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MORTL0',                                           &
                     Default      = 0.064,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR03'
       ! Roots Mortality rate 1/day
        call GetData(Me%Roots%MortalityRate,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MORTR0',                                           &
                     Default      = 0.035,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR04'
        
        
        ! Maximum internal Nitrogen Content gN/KgDW
         call GetData(Me%Nint%Nmax,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NMAX',                                             &
                     Default      = 31.,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR06'
        
        ! Minimum internal Nitrogen Content Nmin = gN/KgDW
         call GetData(Me%Nint%Nmin,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NMIN',                                             &
                     Default      = 5.,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR07'
        ! Critical  internal Nitrogen Content Ncrit = gN/KgDW
        call GetData(Me%Nint%Ncrit,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NCRIT',                                            &
                     Default      = 16.,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR08'
       
        ! Maximum internal Phosphorus Content Pmax = gP/KgDW
         call GetData(Me%Pint%Pmax,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PMAX',                                             &
                     Default      = 3.14,                                                 &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR09'
        ! Minimum internal Phosphorus Content Pmin = gP/KgDW
          call GetData(Me%Pint%Pmin,                                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PMIN',                                             &
                     Default      = 0.14,                                                 &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR10'
        !Critical  internal  Phosphorus Pcrit = gP/KgDW
        call GetData(Me%Pint%Pcrit,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PCRIT',                                            &
                     Default      = 0.8,                                                 &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR11'
       ! Ktr is the traslocation coefficient between roots and leaves (dimensionless)
        call GetData(Me%Leaves%Ktr,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'KTR',                                              &
                     Default      = 0.28,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR12'
    
  
       !Critical shear stress for macroalgae dettachment to occur (in Pa)
        call GetData(Me%Leaves%ErosCritShear,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'EROCRIT',                                          &
                     Default      = 2.,                                              &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR17'
        
            !minimum leaves biomass (kgdw /m2)
        call GetData(Me%Leaves%MinimumBiomass,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MINLEAVES',                                        &
                     Default      = 0.0001,                                             &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR18'   
        
        !minimum roots biomass (kgdw/m2)
        call GetData(Me%Roots%MinimumBiomass,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MINROOTS',                                         &
                     Default      = 0.0001,                                             &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR19'
   
                !minimum roots biomass (kgdw/m2)
        call GetData(Me%Leaves%gNKgDW,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GNKGDW',                                           &
                     Default      = 16.,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR20'
        
        call GetData(Me%Leaves%gPKgDW,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GPKGDW',                                           &
                     Default      = 1.8,                                                 &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR21'
        
         call GetData(Me%Leaves%gCKgDW,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GCKGDW',                                           &
                     Default      = 330.,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR22'
        
        call GetData(Me%Leaves%TOptMin,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TOPTMIN',                                          &
                     Default      = 13.,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR23'
        
        call GetData(Me%Leaves%TOptMax,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TOPTMAX',                                          &
                     Default      = 28.,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR24'
        
        call GetData(Me%Leaves%TMin,                                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TMIN',                                             &
                     Default      = 6.,                                                 &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR25'
        
        call GetData(Me%Leaves%TMax,                                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TMAX',                                             &
                     Default      = 37.,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR26'
        
        call GetData(Me%Leaves%K1,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K1',                                               &
                     Default      = 0.3,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR27'
        
         call GetData(Me%Leaves%K2,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K2',                                               &
                     Default      = 0.98,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR28'
        
         call GetData(Me%Leaves%K3,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K3',                                               &
                     Default      = 0.98,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR29'
        
         call GetData(Me%Leaves%K4,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K4',                                               &
                     Default      = 0.02,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR30'
        
                 call GetData(Me%Leaves%Lat,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'LAT',                                              &
                     Default      = 40.5,                                               &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR31'
        
                call GetData(Me%Leaves%KMAX,                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'KMAX',                                 &
                     Default      = 0.23,                                            &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR32'
   
               call GetData(Me%Leaves%MortType,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MORT_TYPE',                                          &
                     Default      = 1,                                                &
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassesParameters - ModuleBenthicEcology - ERR33'

   
    end subroutine ReadSeagrassesParameters
    
    !--------------------------------------------------------------------------
    
   
    subroutine ReadNitrogenParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Nitrogen%PONDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PON_DECAY_RATE',                                   &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadNitrogenParameters - ModuleBenthos - ERR01'
        
        call GetData(Me%Nitrogen%PONDecayTFactor,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PON_DECAY_TFACTOR',                                &
                     Default      = 1.02,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadNitrogenParameters - ModuleBenthos - ERR02'

    end subroutine ReadNitrogenParameters

    !--------------------------------------------------------------------------

    subroutine ReadPhosphorusParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Phosphorus%POPDecayRate,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POP_DECAY_RATE',                                   &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhosphorusParameters - ModuleBenthos - ERR01'

        call GetData(Me%Phosphorus%POPDecayTFactor,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POP_DECAY_TFACTOR',                                &
                     Default      = 1.08,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhosphorusParameters - ModuleBenthos - ERR02'


    end subroutine ReadPhosphorusParameters

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine ReadPhytoParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetData(Me%Phyto%MortalityRate,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO_MORTALITY',                                  &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhytoParameters - ModuleBenthos - ERR01'

        call GetData(Me%Phyto%NC_Ratio,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO_NC_RATIO',                                   &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhytoParameters - ModuleBenthos - ERR10'


        call GetData(Me%Phyto%PC_Ratio,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO_PC_RATIO',                                   &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhytoParameters - ModuleBenthos - ERR20'


    end subroutine ReadPhytoParameters

    !--------------------------------------------------------------------------

    subroutine ReadSilicaParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Silica%BioSiDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BIOSI_DECAY_RATE',                                 &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSilicaParameters - ModuleBenthos - ERR01'


    end subroutine ReadSilicaParameters
    
    !--------------------------------------------------------------------------
        
    subroutine ConstructProducers

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type (T_Producer),      pointer           :: NewProducer
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_producer>',   &
                                        block_end       = '<end_producer>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    call AddProducer                    (NewProducer)

                    call ConstructProducerParameters    (NewProducer)

                    nullify(NewProducer)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructProducers - ModuleBenthicEcology - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructProducers - ModuleBenthicEcology - ERROR #2'
            else cd1
                    stop       'ConstructProducers - ModuleBenthicEcology - ERROR #3'
            end if cd1
        end do do1

    end subroutine ConstructProducers


    !--------------------------------------------------------------------------
subroutine AddProducer (ObjProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer),      pointer           :: ObjProducer
        !Local-----------------------------------------------------------------
        type (T_Producer),      pointer           :: PreviousProducer
        type (T_Producer),      pointer           :: NewProducer
        integer, save                             :: NextProducerID = 1

        !Allocates new Producer
        allocate (NewProducer)
        nullify  (NewProducer%Next)

        !Insert new Producer into list and makes current algae point to it
        if (.not. associated(Me%FirstProducer)) then
            Me%FirstProducer            => NewProducer
            ObjProducer                 => NewProducer
        else
            PreviousProducer            => Me%FirstProducer
            ObjProducer                 => Me%FirstProducer%Next

            do while (associated(ObjProducer))
                PreviousProducer        => ObjProducer
                ObjProducer             => ObjProducer%Next
            enddo
            ObjProducer                 => NewProducer
            PreviousProducer%Next       => NewProducer
        endif

        !Attributes ID
        ObjProducer%ID%ID               = NextProducerID

        NextProducerID                  = NextProducerID + 1


    end subroutine AddProducer
    
    !--------------------------------------------------------------------------
    
        subroutine ConstructProducerParameters (NewProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer),      pointer           :: NewProducer
        

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                   :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewProducer%ID%Name,                       &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NAME',                     &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #1'

        call GetData(NewProducer%ID%Description,                &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'DESCRIPTION',              &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #2'
        
        !producer mortality rate 1/day
     
        call GetData(NewProducer%MortalityRate,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'MORTALITY_RATE',               &
                     Default      = 0.02,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #4'
     
        ! producer maximum growth rate (1/day)
        call GetData(NewProducer%Vmax,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'VMAX',               &
                     Default      = 2.,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #5'

        !Producer half saturation constant for N uptake g N/l
        call GetData(NewProducer%KN,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'KN',               &
                     Default      = 0.014e-3,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #6'
        
       !Producer half saturation constant for N uptake g P/l
       call GetData(NewProducer%KP,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'KP',               &
                     Default      = 0.001e-3,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #7'
        
        ! Slope of P/I curve 1/[day*(W/m2)] 
        call GetData(NewProducer%alpha,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'ALPHA',               &
                     Default      = 0.025,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #8'
        
         !ratio N:C in producer (gN/gC)
         call GetData(NewProducer%NCratio,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NCRATIO',               &
                     Default      = 0.18,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #9'
        !ratio N:C in producer (gP/gC)
        call GetData(NewProducer%PCratio,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'PCRATIO',               &
                     Default      = 0.024,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #10'
        
         ! minimum temperature of the optimal interval for the microphytobenthos 
            !growth, ºC
        call GetData(NewProducer%TOptPhytoMin,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TOPTFMIN',               &
                     Default      = 25.,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #11'
        
       
        !TOptPhytoMax, maximum temperature of the optimal interval for the microphytobenthos 
        !growth, oC
        call GetData(NewProducer%TOptPhytoMax,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TOPTFMAX',               &
                     Default      = 26.5,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #12' 



        !TPhytoMin, minimum tolerable temperature of the  interval for the microphytobenthos 
        !growth, ºC
            
            call GetData(NewProducer%TPhytoMin,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TFMIN',               &
                     Default      = 10.,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #13' 
        
          
        !TPhytoMax, maximum tolerable temperature of the  interval for the microphytobenthos 
        !growth, oC
             call GetData(NewProducer%TPhytoMax,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TFMAX',               &
                     Default      = 30.,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #14' 
            
        !FK1, constant to control temperature response curve shape
           call GetData(NewProducer%FK1,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TFCONST1',               &
                     Default      = 0.05,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #15' 

       !FK2, constant to control temperature response curve shape
           call GetData(NewProducer%FK2,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TFCONST2',               &
                     Default      = 0.98,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #16' 

       !FK3, constant to control temperature response curve shape
           call GetData(NewProducer%FK3,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TFCONST3',               &
                     Default      = 0.98,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #17' 

       !FK4, constant to control temperature response curve shape
           call GetData(NewProducer%FK4,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'TFCONST4',               &
                     Default      = 0.02,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #18' 
        
               !EROCRIT, shear stress at which erosion of producers occurs (Pa)
           call GetData(NewProducer%EroCritShear,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'EROCRIT',               &
                     Default      = 2.,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #19' 
        
               !MINBIOMASS, minimum producers biomass (kgC/m2)
           call GetData(NewProducer%MinimumBiomass,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'MINBIOMASS',               &
                     Default      = 0.1e-5,                        &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #20'         

       ! Fraction of respiration (-)
        call GetData(NewProducer%RespFrac,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'RESPFRAC',               &
                     Default      = 0.05,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #25'
        
        !Minimum producer's biomass that limits the growth rate (Kg C/m2)
        call GetData(NewProducer%MINP,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'MINP',               &
                     Default      = 0.001,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #27'
        
       !Maximum producer's biomass that limits the growth rate (Kg C/m2)
                call GetData(NewProducer%MAXP,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'MAXP',               &
                     Default      = 0.005,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #29'
        
    end subroutine ConstructProducerParameters


!_______________________________________________________________________________________________________
!__________________________________________________________________________________ 
!_______________________________________________________________________________________________________
!________________________________________________________________________________

    subroutine ConstructConsumers

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type (T_Consumer),      pointer           :: NewConsumer
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_consumer>',   &
                                        block_end       = '<end_consumer>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    call AddConsumer                    (NewConsumer)

                    call ConstructConsumerParameters    (NewConsumer, ClientNumber)

                    nullify(NewConsumer)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructConsumer - ModuleBenthicEcology - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructConsumer - ModuleBenthicEcology - ERROR #2'
            else cd1
                    stop       'ConstructConsumer - ModuleBenthicEcology - ERROR #3'
            end if cd1
        end do do1

    end subroutine ConstructConsumers


    !--------------------------------------------------------------------------

!--------------------------------------------------------------------------


    subroutine AddConsumer (ObjConsumer)

        !Arguments-------------------------------------------------------------
        type (T_Consumer),      pointer           :: ObjConsumer
        !Local-----------------------------------------------------------------
        type (T_Consumer),      pointer           :: PreviousConsumer
        type (T_Consumer),      pointer           :: NewConsumer
        integer, save                             :: NextConsumerID = 1

        !Allocates new Consumer
        allocate (NewConsumer)
        nullify  (NewConsumer%Next)

        !Insert new Consumer into list and makes current ?? point to it
        if (.not. associated(Me%FirstConsumer)) then
            Me%FirstConsumer            => NewConsumer
            ObjConsumer                 => NewConsumer
        else
            PreviousConsumer            => Me%FirstConsumer
            ObjConsumer                 => Me%FirstConsumer%Next

            do while (associated(ObjConsumer))
                PreviousConsumer        => ObjConsumer
                ObjConsumer             => ObjConsumer%Next
            enddo
            ObjConsumer                 => NewConsumer
            PreviousConsumer%Next       => NewConsumer
        endif

        !Attributes ID
        ObjConsumer%ID%ID               = NextConsumerID

        NextConsumerID                  = NextConsumerID + 1


    end subroutine AddConsumer
    
    !--------------------------------------------------------------------------


    subroutine ConstructConsumerParameters (NewConsumer, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Consumer),      pointer           :: NewConsumer
        integer                                   :: ClientNumber
        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                   :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewConsumer%ID%Name,                       &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NAME',                     &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #1'
        
        

        if(.not. CheckPropertyName(trim(NewConsumer%ID%Name)//" carbon"))&
       stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #2'
 
     
 
      ! mortality rate in 1/day
       call GetData(NewConsumer%BETA20,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'BETA20',                      &   ! 
                     Default      = 0.002,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #5'
      
      ! if filter feeder, it is the filtration rate, for example 0.216e3 l/d/g C 
      ! If deposit feeder, it is the ingestion rate for example  0.05 1/day 
       call GetData(NewConsumer%GRMAX,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'GRMAX',                        &
                     Default      = 0.05,                         &   
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #6'
        
      
       ! respiration rate in d-1
       call GetData(NewConsumer%RESP20,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'RESP20',                &   ! 
                     Default      = 0.013,                         & 
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #7'
        
         ! temperature decay factor, dimensionless
             call GetData(NewConsumer%TemperatureFactor,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'TFAC',                &   ! 
                     Default      = 1.08,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #8' 
        
        ! Ratio N:C in consumer gN/gC
        call GetData(NewConsumer%NCRatio,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'NCRATIO',                &   ! 
                     Default      = 0.18,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #9'  

       ! Ratio P:C in consumer
        call GetData(NewConsumer%PCRatio,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'PCRATIO',                &   ! 
                     Default      = 0.024,                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #10' 
        
        
        ! Minimum temperature of the optimal interval for consumer growth ºC
        call GetData(NewConsumer%TOptMin,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'TOPTMIN',                                         &
                     Default      = 13.,                                                &                                
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #11'
        
        ! Maximum temperature of the optimal interval for consumer growth ºC
        call GetData(NewConsumer%TOptMax,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'TOPTMAX',                                           &
                     Default      = 28.,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #12'
       
        ! Minimum temparature for consumer growth ºC 
        call GetData(NewConsumer%TMin,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'TMIN',                                           &
                     Default      = 6.,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #13'
        
        ! Maximum temparature consumers growth ºC 
        call GetData(NewConsumer%TMax,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'TMAX',                                           &
                     Default      = 37.,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #14'
        
        ! Constant to control temperature response curve shape
        call GetData(NewConsumer%K1,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K1',                                           &
                     Default      = 0.3,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #15'
        
          ! Constant to control temperature response curve shape
         call GetData(NewConsumer%K2,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K2',                                           &
                     Default      = 0.98,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #16'
        
           ! Constant to control temperature response curve shape
         call GetData(NewConsumer%K3,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K3',                                           &
                     Default      = 0.98,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #20'
        
          ! Constant to control temperature response curve shape
         call GetData(NewConsumer%K4,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'K4',                                           &
                     Default      = 0.02,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #25'
       
        ! Min consumer density that limits the growth rate (kgC/m2)
        call GetData(NewConsumer%SMIN,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SMIN',                                           &
                     Default      = 0.005,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #30'
        
        ! Max consumer density that limits the growth rate (kgC/m2)
         call GetData(NewConsumer%SMAX,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SMAX',                                           &
                     Default      = 0.020,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #35'
        
        ! Oxygen concentration limitation constant (mg O2/l)
          call GetData(NewConsumer%KO2,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'KO2',                                           &
                     Default      = 0.5,                                                & 
                     ClientModule = 'ModuleBenthicEcology',                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #40'
        
        
        call ConstructGrazing               (NewConsumer%Grazing , ClientNumber)
        
    
    
    
 end subroutine ConstructConsumerParameters
!--------------------------------------------------------------------------


! Bivalve (DEB) model

   subroutine ConstructBivalveDEB

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type (T_BivalveDEB),      pointer         :: NewBivalveDEB
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_BivalveDEB>', &
                                        block_end       = '<end_BivalveDEB>',   &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    call AddBivalveDEB                  (NewBivalveDEB)

                    call ConstructBivalveDEBParameters  (NewBivalveDEB, ClientNumber)

                    nullify(NewBivalveDEB)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructBivalveDEB - ModuleBenthicEcology - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructBivalveDEB - ModuleBenthicEcology - ERROR #2'
            else cd1
                    stop       'ConstructBivalveDEB - ModuleBenthicEcology - ERROR #3'
            end if cd1
        end do do1

    end subroutine ConstructBivalveDEB


   !--------------------------------------------------------------------------


   !--------------------------------------------------------------------------


    subroutine AddBivalveDEB (ObjBivalveDEB)

        !Arguments-------------------------------------------------------------
        type (T_BivalveDEB),      pointer           :: ObjBivalveDEB
        !Local-----------------------------------------------------------------
        type (T_BivalveDEB),      pointer           :: PreviousBivalveDEB
        type (T_BivalveDEB),      pointer           :: NewBivalveDEB
        integer, save                               :: NextBivalveDEBID = 1

        !Allocates new BivalveDEB
        allocate (NewBivalveDEB)
        nullify  (NewBivalveDEB%Next)

        !Insert new BivalveDEB into list and makes current ?? point to it
        if (.not. associated(Me%FirstBivalveDEB)) then
            Me%FirstBivalveDEB            => NewBivalveDEB
            ObjBivalveDEB                 => NewBivalveDEB
        else
            PreviousBivalveDEB            => Me%FirstBivalveDEB
            ObjBivalveDEB                 => Me%FirstBivalveDEB%Next

            do while (associated(ObjBivalveDEB))
                PreviousBivalveDEB        => ObjBivalveDEB
                ObjBivalveDEB             => ObjBivalveDEB%Next
            enddo
            ObjBivalveDEB                 => NewBivalveDEB
            PreviousBivalveDEB%Next       => NewBivalveDEB
        endif

        !Attributes ID
        ObjBivalveDEB%ID%ID               = NextBivalveDEBID

        NextBivalveDEBID                  = NextBivalveDEBID + 1


    end subroutine AddBivalveDEB
    
    !--------------------------------------------------------------------------


subroutine ConstructBivalveDEBParameters (NewBivalveDEB, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_BivalveDEB),      pointer           :: NewBivalveDEB
        integer                                     :: ClientNumber
        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: FromBlock 
        integer                                     :: STATUS, i
        real, dimension(:), allocatable             :: AuxVector

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewBivalveDEB%ID%Name,                     &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NAME',                     &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #1'
        
        

        if(.not. CheckPropertyName(trim(NewBivalveDEB%ID%Name)//" volume"))&
       stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #2'
 
       
       call GetData(NewBivalveDEB%NInd,                                    &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'NINDM2',                              &  
                     Default      = 100,                                   &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #4'
 
       call GetData(NewBivalveDEB%EM,                                      &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'MAXEQUE',                             &  
                     Default      = 2190.0e6,                              &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #5'
      
       call GetData(NewBivalveDEB%Eg,                                      &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'EGCC',                                &  
                     Default      = 1900.0e6,                              &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #6'
      
        call GetData(NewBivalveDEB%Pxm,                                    &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'MAXSAING',                            &  
                     Default      = 22.7778,                               &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #7'

        call GetData(NewBivalveDEB%ka,                                     &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'ASSEFFIC',                            &  
                     Default      = 0.75,                                  &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #8'
         
        call GetData(NewBivalveDEB%Pm,                                     &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'SOMACOST',                            &  
                     Default      = 277.7778,                              &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #9'
        
        call GetData(NewBivalveDEB%kappa,                                  &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'FLUXRESFRAC',                         &  
                     Default      = 0.70,                                  &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #10'
                                      
        call GetData(NewBivalveDEB%Xk,                                      &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'SATCOEF',                                &  
                     Default      = 1.77,                              &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #11'
                                  
        call GetData(NewBivalveDEB%kR,                                     &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'REPOREFFIC',                          &  
                     Default      = 0.9,                                   &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #12'
                                  
        call GetData(NewBivalveDEB%Conv_fac,                               &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'ENERCONT',                            &  
                     Default      = 17550.0,                               &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #13'
                              
        call GetData(NewBivalveDEB%St_DW_perc,                             &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'WWTODW',                              &  
                     Default      = 0.20,                                  &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #14'
                            
        call GetData(NewBivalveDEB%spawn_eff,                              &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'SPAWEFFIC',                           &  
                     Default      = 0.90,                                  &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #15'
                                   
         call GetData(NewBivalveDEB%Tspawn,                                &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'TSPAWN',                              &  
                     Default      = 2.0,                                   &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #16'
                                
         call GetData(NewBivalveDEB%RGS,                                   &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'GSOMINDEX',                           &  
                     Default      = 0.4,                                   &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #17'
                
        call GetData(NewBivalveDEB%spawnD,                                 &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'NSPAWND',                             &  
                     Default      = 1,                                     &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #18'
        
        allocate (NewBivalveDEB%spawnDay(1:NewBivalveDEB%spawnD), STAT = STATUS)
        allocate (AuxVector(1:NewBivalveDEB%spawnD), STAT = STATUS)
        
            call GetData(AuxVector, Me%ObjEnterData, iflag,                    &
                         SearchType   = FromBlock,                             &
                         keyword      = 'JDAYSPAWN',                           &
                         Default      = 180.,                                  & 
                         STAT         = STATUS)
            if (STATUS .NE. SUCCESS_ .or. iflag == 0) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #19'
                         
                         do i = 1, NewBivalveDEB%spawnD
                            NewBivalveDEB%spawnDay(i) = AuxVector(i)
                         enddo

        deallocate (AuxVector)                       
                       
                        
        call GetData(NewBivalveDEB%volp,                                   &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'VOLADULT',                            &  
                     Default      = 0.06e-6,                               &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #20'
                                 
        call GetData(NewBivalveDEB%delm,                                   &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'SHAPEPARAM',                          &  
                     Default      = 0.287,                                 &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #21'
                                 
        call GetData(NewBivalveDEB%Tref,                                   &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'REFTEMP',                             &  
                     Default      = 293.0,                                 &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #22'
                          
        call GetData(NewBivalveDEB%Ta,                                     &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'TEMPARRH',                            &  
                     Default      = 5800.0,                                &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #23'

        call GetData(NewBivalveDEB%TL,                                     &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'TEMPLOW',                             &  
                     Default      = 275.0,                                 &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #24'

        call GetData(NewBivalveDEB%THresp,                                 &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'THRESP',                              &  
                     Default      = 296.0,                                 &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #25'
                                  
        call GetData(NewBivalveDEB%THing,                                  &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'THING',                               &  
                     Default      = 296.0,                                 &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #26'

        call GetData(NewBivalveDEB%TAL,                                    &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'TAL',                                 &  
                     Default      = 45430.0,                               &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #27'
                             
        call GetData(NewBivalveDEB%TAHresp,                                &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'TAHRESP',                             &  
                     Default      = 31376.0,                               &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #28'

        call GetData(NewBivalveDEB%TAHing,                                 &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'TAHING',                              &  
                     Default      = 31376.0,                               &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #29'
                             
        call GetData(NewBivalveDEB%kd,                                     &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'FECALDECAY',                          &  
                     Default      = 2.1222e-4,                             &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #30'

        call GetData(NewBivalveDEB%Ec,                                     &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'ENERGPHY',                            &  
                     Default      = 47.7546,                               &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #31'
                                  
            call GetData(NewBivalveDEB%etaO2,                              &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'ENERGO2',                             &  
                     Default      = 14.3,                                  &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #32'
                                
            call GetData(NewBivalveDEB%CtoChla,                            &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'CCHLABIV',                            &  
                     Default      = 60.0,                                  &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #33'

            call GetData(NewBivalveDEB%NH4Frac,                            &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'NH4EXCFRAC',                          &  
                     Default      = 0,                                     &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #34'

               call GetData(NewBivalveDEB%BiLen1,                          &
                     Me%ObjEnterData, iflag,                               &
                     SearchType   = FromBlock,                             &
                     keyword      = 'BILENINIC',                           &  
                     Default      = 0.05,                                  &  
                     ClientModule = MohidModules(mBenthicEcology_)%Name,   &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBivalveDEBParameters - ModuleBenthicEcology - ERROR #35'
 
    
 end subroutine ConstructBivalveDEBParameters
!--------------------------------------------------------------------------



subroutine ConstructGrazing (Grazing, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Grazing)                                :: Grazing
        integer                                         :: ClientNumber

        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        logical                                         :: BlockInBlockFound
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 
        integer                                         :: FirstLine, LastLine
        type (T_Food), pointer                          :: NewFood

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

       ! 
        call GetData(Grazing%Ass_Efic,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ASS_EFIC',                &
                     Default      = 0.25,                                                &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #3'
        
      call GetData(Grazing%Cohesivesed,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'COHESIVESED',                &
                     Default      = .false.,                                                &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #4'
        
     
     If (Grazing%Cohesivesed) then 
     call GetData(Grazing%MaxSediment,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'SEDMAX',                &
                     Default      = 0.1,                                                &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #5'
     endif
        
         call GetData(Grazing%KFoodC,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'KFOODC',                &
                     Default      = 0.001,                                                &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #6'
       
  
        
        

        do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,         &
                                       '<begin_food>', '<end_food>',          &
                                       BlockInBlockFound,                     &
                                       FirstLine = FirstLine,                 &
                                       LastLine  = LastLine,                  &
                                       STAT      = STAT_CALL)
            if      (STAT_CALL .EQ. SUCCESS_) then    
            
                if (BlockInBlockFound) then
            
                    if     (((LastLine + 1) - (FirstLine - 1)) .GE. 1) then

                        call AddFood        (Grazing, NewFood)

                        call ConstructFood  (NewFood)

                        nullify(NewFood)

                    else
                        write(*,*)  
                        write(*,*) 'Error counting Foods. '
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #4'
                    end if
                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #5'
                    
                    exit 

                end if 

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBlock. '
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #6'
            end if

        end do
        
    end subroutine ConstructGrazing
    
    
    ! -----------------------------------------------------------------------------------------
    subroutine AddFood (Grazing, Food)

        !Arguments-------------------------------------------------------------
        type (T_Grazing)                          :: Grazing
        type (T_Food),          pointer           :: Food

        !Local-----------------------------------------------------------------
        type (T_Food),          pointer           :: PreviousFood
        type (T_Food),          pointer           :: NewFood
        integer, save                             :: NextFoodID = 1

        !Allocates new Producer
        allocate (NewFood)
        nullify  (NewFood%Next)

        !Insert new Food into list and makes current algae point to it
        if (.not. associated(Grazing%FirstFood)) then
            NextFoodID              = 1
            Grazing%FirstFood       => NewFood
            Food                    => NewFood
        else
            PreviousFood            => Grazing%FirstFood
            Food                    => Grazing%FirstFood%Next

            do while (associated(Food))
                PreviousFood        => Food
                Food                => Food%Next
            enddo
            Food                    => NewFood
            PreviousFood%Next       => NewFood
        endif

        !Attributes ID
        Food%ID%ID               = NextFoodID

        NextFoodID               = NextFoodID + 1


    end subroutine AddFood
!--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------


    subroutine ConstructFood (NewFood)

        !Arguments-------------------------------------------------------------
        type (T_Food), pointer                          :: NewFood

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag
        integer                                         :: FromBlockInBlock      

        !Local-----------------------------------------------------------------
        type(T_Producer), pointer                       :: Producer

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlockInBlock = FromBlockInBlock)

        call GetData(NewFood%ID%Name,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'NAME',                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #1'
        
                call GetData(NewFood%NCRatio,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'NCRATIO',                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #2'
                 
                 call GetData(NewFood%PCRatio,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'PCRATIO',                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #3'
        
        call GetData(NewFood%Use_Carbon,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'CARBON_USE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #4'
        
        call GetData(NewFood%Use_Nitrogen,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'NITROGEN_USE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #5'
        
                call GetData(NewFood%Use_Phosphorus,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'PHOSPHORUS_USE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #6'
        
        call GetData(NewFood%ParticulateWaterFood,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'PARTICWATERFOOD',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #7'
        
        
       call GetData(NewFood%Minval,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'MINVAL',                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFood - ModuleBenthicEcology - ERROR #8'
        
        if((Me%PelagicModel == LifeModel))then
        if(.not. CheckPropertyName(trim(NewFood%ID%Name)//" nitrogen"))&
        stop 'ConstructFood - ModuleBenthicEcology - ERROR #9'
        endif        

        Producer => Me%FirstProducer
        do while(associated(Producer))

            !if(NewFood%ID%Name == trim(Producer%ID%Name))then

             !   NewFood%Use_Chl = .true.

              !  if(Producer%Use_Silica)NewFood%Use_Silica = .true.

            !end if

         

            Producer => Producer%Next
        end do

    end subroutine ConstructFood
  
   
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine GetDTBenthicEcology(BenthicEcology_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcology_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcology_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetDTBenthicEcology
    
    !--------------------------------------------------------------------------

    
    subroutine GetBenthicEcologyPropertyList(Life_ID, PropertyList, STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: Life_ID
        integer, dimension(:), pointer                          :: PropertyList
        integer, optional, intent(OUT)                          :: STAT

        !External--------------------------------------------------------------
        integer                                                 :: ready_              

        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Life_ID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBenthicEcology_, Me%InstanceID)

            PropertyList => Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetBenthicEcologyPropertyList
    
    
        !--------------------------------------------------------------------------

    subroutine GetBenthicEcologySize(BenthicEcology_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcology_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcology_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBenthicEcologySize
    
    !--------------------------------------------------------------------------
    
     subroutine GetBenthicEcologyPropIndex (BenthicEcology_ID, PropertyIDNumber, PropertyIndex, STAT)

                                     

        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcology_ID
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

        call Ready(BenthicEcology_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
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

    end subroutine GetBenthicEcologyPropIndex


    !--------------------------------------------------------------------------

   
                
                
    integer function SearchPropIndex (PropIDNumber)

        !Arguments-------------------------------------------------------------
        integer,  intent(IN )                         :: PropIDNumber
        !Local-----------------------------------------------------------------
        integer                                       :: CurrentIndex

        !----------------------------------------------------------------------

        SearchPropIndex = UNKNOWN_

        do CurrentIndex = Me%Prop%ILB, Me%Prop%IUB

            if (PropIDNumber == Me%PropertyList(CurrentIndex))then
                SearchPropIndex = CurrentIndex
                exit
            end if
                    
        end do

    end function SearchPropIndex

    !--------------------------------------------------------------------------

       
       
        subroutine UnGetBenthicEcology(BenthicEcologyID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BenthicEcologyID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcologyID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBenthicEcology_, Me%InstanceID, "UnGetBenthicEcology3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBenthicEcology
    
     !--------------------------------------------------------------------------
        
        
        
        subroutine GetBenthicEcologyRateFlux(BenthicEcologyID, FirstProp, SecondProp, RateFlux, STAT)


        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcologyID
        integer,           intent(IN )      :: FirstProp, SecondProp
        real,    dimension(:), pointer      :: RateFlux
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcologyID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBenthicEcology_, Me%InstanceID)

            select case(FirstProp)
                 case(NintFactor_)
                     RateFlux => Me%StoredArray(:, Me%StoredIndex%NintFactor)
                 
                 case(PintFactor_)
                     RateFlux => Me%StoredArray(:, Me%StoredIndex%PintFactor)       
                
                 case(NintFactorR_)
                     RateFlux => Me%StoredArray(:, Me%StoredIndex%NintFactorR)
                 
                 case(RootsMort_)
                     RateFlux => Me%StoredArray(:, Me%StoredIndex%RootsMort)
                 
                 case(PintFactorR_)
                     RateFlux => Me%StoredArray(:, Me%StoredIndex%PintFactorR)
                             

                case default

                    RateFlux => Me%Matrix    (:, FirstProp, SecondProp      )

            end select

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBenthicEcologyRateFlux

    !--------------------------------------------------------------------------
    
      
      
      subroutine UnGetBenthicEcologyRateFlux(BenthicEcologyID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BenthicEcologyID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcologyID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBenthicEcology_, Me%InstanceID, "UnGetBenthicEcologyRateFlux")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBenthicEcologyRateFlux
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    subroutine ModifyBenthicEcology(ObjBenthicEcologyID, & 
                                    Temperature,         & ! Temperature of the cell (k=1), in ºC
                                    WaterVolume,         & ! Volume of the water cell (k=1) (m3)
                                    CellArea,            & ! Area of the water cell (k=1)(m2)
                                    MassInKgFromWater,   & ! MassInKgFromWater is a matrix containing water state variables (Kg)
                                    Sediment,            & ! Sediment concentration in the cell (k=1) (Kg/m3)
                                    ShortWaveAverage,    & ! Incident shortwave radiation at the bottom surface (W/m2)
                                    ShearStress,         & ! Bottom shear stress in Pa
                                    UptakeNH4s,          & ! Uptake of ammonia from seagrasses roots (gN/day) 
                                    UptakeNH4NO3w,       & ! Uptake of Nh4 and No3 by seagrasses leaves (gN/day)
                                    UptakePO4w,          & ! Uptake of phosphate from seagrasses leaves (gP/day) 
                                    UptakePO4s,          & ! Uptake of phosphate from seagrasses roots (gP/day) 
                                    LightFactor,         & ! Light limitation for seagrasses only , dimensionless
                                    JulianDay,           & ! Day of year (days)
                                    OpenPoints,          & ! Calculation points (0=No calculation, 1 =Calculation)
                                    Mass,                & ! Mass is the matrix of bottom state variables, in Kg
                                    STAT)
    
    
            !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthicEcologyID
        real,    dimension(:  ), pointer            :: Temperature
        real,    dimension(: ,: ), pointer          :: MassInKgFromWater
        real,    dimension(:  ), pointer            :: Sediment
        real,    dimension(:  ), pointer, optional  :: UptakeNH4s
        real,    dimension(:  ), pointer, optional  :: UptakeNH4NO3w
        real,    dimension(:  ), pointer, optional  :: UptakePO4w
        real,    dimension(:  ), pointer, optional  :: UptakePO4s
        real,    dimension(:  ), pointer, optional  :: LightFactor
        real(8), dimension(:  ), pointer            :: WaterVolume
        real,    dimension(:  ), pointer            :: CellArea
        integer, dimension(:  ), pointer, optional  :: OpenPoints
        real,    dimension(:  ), pointer, optional  :: ShortWaveAverage
        real,    dimension(:  ), pointer, optional  :: ShearStress
        real,    dimension(:,:), pointer            :: Mass
        integer, intent(OUT), optional              :: STAT !whatever
        integer, intent(IN)                         :: JulianDay
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: Index
        logical                                     :: Compute
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBenthicEcologyID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
        
            Me%JulianDay = JulianDay

            Me%ExternalVar%Temperature  => Temperature
            if (.not. associated(Me%ExternalVar%Temperature))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR01'
                
             Me%ExternalVar%ShortWaveAverage  => ShortWaveAverage
            if (.not. associated(Me%ExternalVar%ShortWaveAverage))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR02'
                
                
             Me%ExternalVar%Sediment  => Sediment
            if (.not. associated(Me%ExternalVar%Sediment))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR03'
                
             Me%ExternalVar%MassInKgFromWater  => MassInKgFromWater
            if (.not. associated(Me%ExternalVar%MassInKgFromWater))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR04'
                
            Me%ExternalVar%WaterVolume  => WaterVolume
            if (.not. associated(Me%ExternalVar%WaterVolume))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR05'
                
            Me%ExternalVar%CellArea  => CellArea
            if (.not. associated(Me%ExternalVar%CellArea))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR06'
                
           Me%ExternalVar%ShearStress  => ShearStress
           if (.not. associated(Me%ExternalVar%ShearStress))       &
               stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR07'
               
            Me%ExternalVar%UptakeNH4s  => UptakeNH4s
            if (.not. associated(Me%ExternalVar%UptakeNH4s))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR08'  
                
            Me%ExternalVar%UptakeNH4NO3w  => UptakeNH4NO3w
            if (.not. associated(Me%ExternalVar%UptakeNH4NO3w))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR09' 
           
            Me%ExternalVar%UptakePO4w  => UptakePO4w
            if (.not. associated(Me%ExternalVar%UptakePO4w))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR10' 
            
               Me%ExternalVar%UptakePO4s  => UptakePO4s
            if (.not. associated(Me%ExternalVar%UptakePO4s))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR11' 
            
            Me%ExternalVar%LightFactor  => LightFactor
           if (.not. associated(Me%ExternalVar%LightFactor))       &
               stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR12' 
              
              
            Me%ExternalVar%Mass         => Mass
            if (.not. associated(Me%ExternalVar%Mass))              &
                stop 'ModifyBenthos - ModuleBenthicEcology      - ERR013'


            do Index = Me%Size%ILB, Me%Size%IUB

                if (present(OpenPoints)) then
                    if (OpenPoints(Index) == OpenPoint) then
                        Compute = .true.
                    else
                        Compute = .false.
                    endif
                else
                    Compute = .true.
                endif

                Me%Matrix(Index,:,:) = 0.
                
                
                                if(Compute)then
                    
                    call ComputeBenthicProducers      (Index)
                    
                    call ComputeBenthicConsumers      (Index)
                    
                    call ComputeBenthicBivalveDEB     (Index)
                    
                
                    if(Me%ComputeOptions%Nitrogen  ) call ComputeBenthicNitrogen    (Index)
                    
                    if(Me%ComputeOptions%Phosphorus) call ComputeBenthicPhosphorus  (Index)
                    
                    if(Me%ComputeOptions%Silica    ) call ComputeBenthicSilica      (Index)
                    
                    if(Me%ComputeOptions%Seagrasses) call ComputeSeagrasses         (index)
                    
                    if(Me%ComputeOptions%Phyto)      call ComputeBenthicPhyto       (index)
                
                endif

            enddo

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
    
    end subroutine ModifyBenthicEcology
    
    !--------------------------------------------------------------------------------------
    
    
    subroutine ComputeBenthicProducers (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        
       
        integer             :: Producer_N, Producer_C, Producer_P
        integer             :: AM, NA, IP, O2
        integer             :: PON, POP, POC
        real                :: AverageRadiation
        real                :: Lightlim, NLim, PLim, NutLim, GrowthRate
        real                :: UptakeNA, UptakeAM, UptakeP, RespirationRate
        real                :: RespirationC, RespirationN, RespirationP
        real                :: MortalityC, MortalityN, MortalityP
        real                :: FluxToPON, FluxToPOP, FluxToPOC
        real                :: x1,x2,x3,x4, AmmoniaPreferenceFactor
        real                :: s1, s2, xa,xb,ya,yb, TemperatureLim,TemperatureFDecay
        real                :: ConcAM, ConcNA, ConcIP, ConcO2
        real                :: ProducerDensity,DensityLim
        integer             :: Zone
        integer, parameter  :: NoLimitation = 1
        integer, parameter  :: Erosion      = 2

    !------------------------------------------------------------------------
       
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        IP      = Me%PropIndex%Phosphate
        O2      = Me%PropIndex%Oxygen
        
        

        
        
        if(Me%PelagicModel == LifeModel) then
        POC     = Me%PropIndex%POC
        endif


        Producer => Me%FirstProducer

d1:     do while(associated(Producer))

        if(Me%ExternalVar%ShearStress(Index)       > Producer%EroCritShear     )then
                    
                    Zone = Erosion
        else
                    
                    Zone = NoLimitation
        endif

        ! Average radiation
        ! I put the maximum between 0.0 and the value of shortwave radiation because during the first temporal iteration, 
        ! the short wave radiation is not calculated yet and it is retreived from ModuleWaterProperties as negative number.
        ! In the other temporal iterations this would not be necessary.
        AverageRadiation = max(0.0, Me%ExternalVar%ShortWaveAverage(index)) !W/m2
        
        ! TemperatureLim is the function expressing the dependence on the temperature
        ! Dimensionless factor 
        ! temperature limitation function is the same as used in modulewaterquality for phytoplankton
        
        !TemperatureLim, temperature effect on producer assimilation rate
        s1 = (1.0 / (Producer%TOptPhytoMin - Producer%TPhytoMin)) * log((Producer%FK2 * (1.0 - Producer%FK1))              &
                                                          / (Producer%FK1 * (1.0 - Producer%FK2)))

        s2 = (1.0 / (Producer%TPhytoMax - Producer%TOptPhytoMax)) * log((Producer%FK3 * (1.0 - Producer%FK4))              &
                                                          / (Producer%FK4 * (1.0 - Producer%FK3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Producer%TPhytoMin))
        yb = exp(s2 * (Producer%TPhytoMax - Me%ExternalVar%Temperature(index)))

        xa = (Producer%FK1 * ya) / (1.0 + Producer%FK1 * (ya - 1.0))
        xb = (Producer%FK4 * yb) / (1.0 + Producer%FK4 * (yb - 1.0))

        TemperatureLim = xa * xb
        
        !-------------------------------------------------------------------------------------------
        
        Producer_N   = Producer%PoolIndex%Nitrogen
        Producer_C   = Producer%PoolIndex%Carbon
        Producer_P   = Producer%PoolIndex%Phosphorus
        
        ! Calculate concentrations only once to reduce number of calculations
                
        ConcAM=Me%ExternalVar%MassInKgFromWater(AM, index)/Me%ExternalVar%WaterVolume(Index)
        ConcNA=Me%ExternalVar%MassInKgFromWater(NA, index)/Me%ExternalVar%WaterVolume(Index)
        ConcIP=Me%ExternalVar%MassInKgFromWater(IP,Index)/Me%ExternalVar%WaterVolume(Index)
        ConcO2=Me%ExternalVar%MassInKgFromWater(O2,Index)/Me%ExternalVar%WaterVolume(Index)
        
        ProducerDensity=Me%ExternalVar%Mass(Producer_C, Index)/Me%ExternalVar%CellArea(Index)
        
        !Density Limitation Le Pape, 1999
       
        DensityLim= 1. - max(0., min(1., (ProducerDensity - Producer%MINP)/(Producer%MAXP-Producer%MINP) ) )
        ! Evans and Parslow model (1985) model 
        ! Producer%alpha has units of 1/(day*(W/m2))
        ! AverageRadiation has units of W/m2
        ! Producer%Vmax has units of 1/day 
        ! Lightlim is dimensionless
        Lightlim =  Producer%alpha* AverageRadiation / &
                    sqrt(Producer%Vmax*Producer%Vmax + (Producer%alpha**2) * ( AverageRadiation**2))
        
        ! Nutrients Limitation (the same as for Phytoplankton in ModuleWaterQuality)
        ! Nlim is dimensionless ; Producer%KN is in Kg N/m3 (that is the same as g N /L)
        NLim =   (ConcAM + ConcNA )                 /  &
                 (Producer%KN + ConcAM + ConcNA)
       
        ! Plim is dimensionless
        ! Producer%KP is in Kg P/m3 (that is the same as g P /L)
        Plim =  ConcIP                             /  &
               (Producer%KP + ConcIP)
        
        ! NutLim is dimensionless
        
        NutLim = min(NLim,PLim)
        
        ! ammonia preference factor (AmmoniaPreferenceFactor,dimensionless)
        ! following the same calculations made in module WaterQuality for phytoplankton
        
            x1 = ConcAM * ConcNA

            x2 = (Producer%KN + ConcAM)   &
               * (Producer%KN + ConcNA) 

            x3 = Producer%KN * ConcAM

            x4 = (ConcAM + ConcNA)        &
               * (Producer%KN + ConcNA)

            if ((x1 .EQ. 0.0) .AND. (x3 .EQ. 0.0)) then
                AmmoniaPreferenceFactor = 0.0                 
            else 
                AmmoniaPreferenceFactor = (x1 / x2) + (x3 / x4)
            end if 
        
       

       
        select case(Zone)
        
        case(NoLimitation)
        
      
        ! growth rate 1/day
        
        GrowthRate=(1.-Producer%RespFrac)*Producer%Vmax*LightLim*NutLim*TemperatureLim*DensityLim
        
        ! uptake of ammonia 
        ! KgN     = 1/day  * (KgN/KgC) * KgC*day
        UptakeAM = AmmoniaPreferenceFactor * GrowthRate * Producer%NCRatio * &
                  Me%ExternalVar%Mass(Producer_C, Index)* Me%DTDay ! KgN
        
        UptakeNA =(1.- AmmoniaPreferenceFactor) * GrowthRate * Producer%NCRatio * &
                  Me%ExternalVar%Mass(Producer_C, Index)* Me%DTDay ! KgN
        
        UptakeP = GrowthRate * Producer%PCRatio * Me%ExternalVar%Mass(Producer_C, Index)* Me%DTDay ! KgP
        
        
        !KgC = KgC *day * 1/day * [-]   
        
        TemperatureFDecay = 1.08 ** (Me%ExternalVar%Temperature(index)-20.) 
        
        if (ProducerDensity .LT. Producer%MinimumBiomass ) then
        
        MortalityC = 0.
        MortalityN = 0.
        MortalityP = 0.
        
        Me%ExternalVar%Mass(Producer_C,     Index) = Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)
        Me%ExternalVar%Mass(Producer_N,     Index) = Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)* Producer%NCRatio
        Me%ExternalVar%Mass(Producer_P,     Index) = Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)* Producer%PCRatio
        
        else   
                                                                                       
        MortalityC = Me%ExternalVar%Mass(Producer_C, Index) * Me%DTDay * Producer%MortalityRate * TemperatureFDecay
        MortalityN = MortalityC*Producer%NCRatio  ! KgN
        MortalityP = MortalityC*Producer%PCRatio  ! KgP
        
        endif
        ! units of mass
        
        

        RespirationRate=Producer%RespFrac*Producer%Vmax*LightLim*NutLim*TemperatureLim*DensityLim
         
        if (MAX( ConcO2 * 1000., Me%Oxygen%Minimum).eq.Me%Oxygen%Minimum) then
                RespirationRate = -1.0 / null_real
        endif
        
       
        if (ProducerDensity .LT. Producer%MinimumBiomass ) then
        
        ! Respiration is set to 0 if the producer's biomass falls below a user-defined minimum value (Producer%MinimumBiomass)
        ! this steps are necessary to avoid negative concentrations, and also to avoid 
        ! errors due to single precision when using the release single executable 
        ! In general it is better to use the release double executable when running interfacesedimentwater,
        ! because values are very small (since they are expressed in Kg)
         
        RespirationC = 0.
        RespirationN = 0.
        RespirationP = 0.
        
        Me%ExternalVar%Mass(Producer_C,     Index) = Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)
        Me%ExternalVar%Mass(Producer_N,     Index) = Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)* Producer%NCRatio
        Me%ExternalVar%Mass(Producer_P,     Index) = Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)* Producer%PCRatio
        
        else 
        
        RespirationC = Me%ExternalVar%Mass(Producer_C, Index) * RespirationRate * Me%DTDay  !KgC
        RespirationN = RespirationC*Producer%NCRatio  ! KgN
        RespirationP = RespirationC*Producer%PCRatio  ! KgP
        
        endif
                  
        ! mass balance
        
        ! KgC
        Me%ExternalVar%Mass(Producer_C,     Index) = Me%ExternalVar%Mass(Producer_C,     Index) + &  ! KgC
                                                     GrowthRate                                 * &  ! 1/day
                                                     Me%ExternalVar%Mass(Producer_C, Index)     * &  ! KgC
                                                     Me%DTDay                                   - &  ! day
                                                     MortalityC                                 - &  ! KgC
                                                     RespirationC                                    ! KgC
        ! KgN
        Me%ExternalVar%Mass(Producer_N,     Index) = Me%ExternalVar%Mass(Producer_N,     Index) + &  ! KgN 
                                                     GrowthRate                                 * &  ! 1/day
                                                     Me%ExternalVar%Mass(Producer_N, Index)     * &  ! KgN 
                                                     Me%DTDay                                   - &  ! day
                                                     MortalityN                                 - &  ! KgN
                                                     RespirationN                                    ! Kg N
        ! KgP
        Me%ExternalVar%Mass(Producer_P,     Index) = Me%ExternalVar%Mass(Producer_P,     Index) + &  ! KgP
                                                     GrowthRate                                 * &  ! 1/day
                                                     Me%ExternalVar%Mass(Producer_P, Index)     * &  ! KgP
                                                     Me%DTDay                                   - &  ! day
                                                     MortalityP                                 - &  ! KgP
                                                     RespirationP                                    ! KgP
           
        if(Me%ComputeOptions%Nitrogen)then
            
            !what passes from The benthic producer to PON
            Me%ExternalVar%Mass(PON,     Index)             = Me%ExternalVar%Mass(PON,     Index) + MortalityN
            !what passes from The benthic producer to Ammonia
            Me%ExternalVar%MassInKgFromWater(AM,     Index) = Me%ExternalVar%MassInKgFromWater(AM,     Index) + RespirationN
            !what passes from nitrate to the benthic producer
            Me%ExternalVar%MassInKgFromWater(NA,     Index) = Me%ExternalVar%MassInKgFromWater(NA,     Index) - UptakeNA
            !what passes from ammonia to benthic producer 
            Me%ExternalVar%MassInKgFromWater(AM,     Index) = Me%ExternalVar%MassInKgFromWater(AM,     Index) - UptakeAM
        end if
           
       if(Me%ComputeOptions%Phosphorus)then
            
            !what passes from The benthic producer to POP
            Me%ExternalVar%Mass(POP,     Index)             = Me%ExternalVar%Mass(POP,     Index) + MortalityP
            !what passes from The benthic producer to Phosphate
            Me%ExternalVar%MassInKgFromWater(IP,     Index) = Me%ExternalVar%MassInKgFromWater(IP,     Index) + RespirationP
            !what passes from Phosphate to benthic producer 
            Me%ExternalVar%MassInKgFromWater(IP,     Index) = Me%ExternalVar%MassInKgFromWater(IP,     Index) - UptakeP
       end if
       
       Me%ExternalVar%MassInKgFromWater(O2, Index)          = Me%ExternalVar%MassInKgFromWater(O2, Index) + &  ! Kg O2 +
                                                             (GrowthRate                                 * &   ! (1/day *
                                                              Me%ExternalVar%Mass(Producer_C, Index)     * &   !  Kg C  *
                                                              Me%DTDay                                   - &   !  day   -
                                                              RespirationC)*32./12.                            !  Kg C) * KgO2/KgC
       
       ! the ratio 32/12 is the ratio between grams of O2 produced and grams of Carbon in the photosynthesis reaction
       ! the same ratio is applied in the ModuleWaterquality for phytoplankton, and in the ModuleMacroalgae for macroalgae
       ! the same ratio is applied to calculate the oxygen consumed by benthic algae

        if(Me%PelagicModel == LifeModel) then
       ! What passes from The benthic producer to POC (only if ModuleLife is active, 
       ! because ModuleWaterQuality does not account for POC)
        Me%ExternalVar%Mass(POC,     Index) = Me%ExternalVar%Mass(POC,     Index) + MortalityC
        endif
        
        
        case(Erosion)
        
        
         !Not all benthic producers mass is eroded. to make sure 
         !there's always enough to grow back again (minimum concentration)
         ! This is similar to ModuleMacroalgae
                   
                     if(Me%ExternalVar%Mass(Producer_C, Index) > Producer%MinimumBiomass)then
                    
                    
                      ! what passes from benthic producer to POM

                      ! KgC  (MinimumBiomass is expressed as KgC/m2, so it is necessary to multiply by
                      ! the cell area to get the value in KgC)
                      FluxToPOC = Me%ExternalVar%Mass(Producer_C, Index) - &
                                  Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)
                                                  

                     ! KgN  (MinimumBiomass is expressed as KgC/m2, so it is necessary to multiply by the
                     ! ratio  N:C and by cell area to get the value in KgN)
                      FluxToPON = Me%ExternalVar%Mass(Producer_N, Index) - &
                                  Producer%MinimumBiomass*Producer%NCRatio*Me%ExternalVar%CellArea(Index)
                                                  
                     ! KgP  (MinimumBiomass is expressed as KgC/m2, so it is necessary to multiply by the
                     ! ratio  P:C and by cell area to get the value in KgP)
                      FluxToPOP   = Me%ExternalVar%Mass(Producer_P, Index) - &
                                   Producer%MinimumBiomass*Producer%PCRatio*Me%ExternalVar%CellArea(Index)
                                                  
                                             
                     else
                     !what passes from benthic Producer to POM:
                     !( if the benthic producer reaches its minimum value, nothing passes)
                     
                     FluxToPOC=0.
                     FluxToPON=0.
                     FluxToPOP=0.
                     end if
                     

                    
                    !only minimum biomass remains:
                    Me%ExternalVar%Mass(Producer_C, Index)=Producer%MinimumBiomass*Me%ExternalVar%CellArea(Index)
                    Me%ExternalVar%Mass(Producer_N, Index)=Producer%MinimumBiomass*Producer%NCRatio*Me%ExternalVar%CellArea(Index)
                    Me%ExternalVar%Mass(Producer_P, Index)=Producer%MinimumBiomass*Producer%PCRatio*Me%ExternalVar%CellArea(Index)
        
                    if(Me%PelagicModel == LifeModel) then
                     
                    !what passes from The benthic producer to POC (only if ModuleLife is active)
                    Me%ExternalVar%Mass(POC,     Index) = Me%ExternalVar%Mass(POC,     Index) + FluxToPOC
                    endif
                    
                    Me%ExternalVar%Mass(PON,     Index) = Me%ExternalVar%Mass(PON,     Index) + FluxToPON
                    Me%ExternalVar%Mass(POP,     Index) = Me%ExternalVar%Mass(POP,     Index) + FluxToPOP
        
        end select


            Producer => Producer%Next
        end do d1

      
        end subroutine ComputeBenthicProducers
    
    
    
    
    !---------------------------------------------------------------------------------------
    
       subroutine ComputeBenthicConsumers (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Consumer),      pointer             :: Consumer
        type(T_Food),          pointer             :: Food
        integer                                    :: PON
        integer                                    :: POP
        integer                                    :: POC
        integer                                    :: Phyto
        integer                                    :: AM
        integer                                    :: IP
        integer                                    :: O2
        integer                                    :: Consumer_N
        integer                                    :: Consumer_C
        integer                                    :: Consumer_P
        integer                                    :: FoodIndexC
        integer                                    :: FoodIndexN
        integer                                    :: FoodIndexP
        real                                       :: TemperatureDependence
        real                                       :: TemperatureFDecay
        real                                       :: MortalityC
        real                                       :: MortalityN
        real                                       :: MortalityP
        real                                       :: MortalityRate
        real                                       :: RespirationC
        real                                       :: RespirationN
        real                                       :: RespirationP
        real                                       :: OxygenLimitation
        real                                       :: SedimentLimitation
        real                                       :: RespirationRate
        real                                       :: IngestionC
        real                                       :: IngestionN
        real                                       :: IngestionP
        real                                       :: IngestionC_tot
        real                                       :: IngestionN_tot
        real                                       :: IngestionP_tot
        real                                       :: EgestionC
        real                                       :: EgestionN
        real                                       :: EgestionP
        real                                       :: EgestionC_tot
        real                                       :: EgestionN_tot
        real                                       :: EgestionP_tot
        real                                       :: IngestionRate,DensityLimitation
        real                                       :: FoodDens,ConsumerDensity
        real                                       :: s1, s2, xa, xb, ya, yb
        
        

        
        AM      = Me%PropIndex%Ammonia
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        Phyto   = Me%PropIndex%Phyto
        IP      = Me%PropIndex%Phosphate
        O2      = Me%PropIndex%Oxygen
        
        if(Me%PelagicModel == LifeModel) then
        POC     = Me%PropIndex%POC
        endif
        

        
        Consumer => Me%FirstConsumer

d1:     do while(associated(Consumer))

        ! temperature effect on consumer
        s1 = (1. / (Consumer%TOptMin - Consumer%TMin)) * log((Consumer%K2 * (1.0 - Consumer%K1))                    &
                                                     / (Consumer%K1 * (1.0 - Consumer%K2)))

        s2 = (1. / (Consumer%TMax - Consumer%TOptMax)) * log((Consumer%K3 * (1.0 - Consumer%K4))                    &
                                                     / (Consumer%K4 * (1.0 - Consumer%K3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Consumer%TMin))
        yb = exp(s2 * (Consumer%TMax - Me%ExternalVar%Temperature(index)))

        xa = (Consumer%K1 * ya) / (1.0 + Consumer%K1 * (ya - 1.0))
        xb = (Consumer%K4 * yb) / (1.0 + Consumer%K4 * (yb - 1.0))
        ! temperature dependence follows a bell shaped function when multiplied by the maximum growth rate
        TemperatureDependence = xa * xb
        
        !  respiration and mortality increase with temperature
        
        TemperatureFDecay = Consumer%TemperatureFactor ** (Me%ExternalVar%Temperature(index)-20.)
        
        Consumer_N   = Consumer%PoolIndex%Nitrogen
        Consumer_C   = Consumer%PoolIndex%Carbon
        Consumer_P   = Consumer%PoolIndex%Phosphorus 
       
        ConsumerDensity=Me%ExternalVar%Mass(Consumer_C, Index)/Me%ExternalVar%CellArea(Index)
        
        DensityLimitation= 1. - max(0., min(1., (ConsumerDensity - Consumer%SMIN)/(Consumer%SMAX-Consumer%SMIN) ) )
        
        ! cohesive sediment is  g/l
        
        ! 
        
        ! Oxygen Limitation
        !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max((Me%ExternalVar%MassInKgFromWater(O2,Index)/ &
        Me%ExternalVar%WaterVolume(index))*1000., Me%Oxygen%Minimum)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + Consumer%KO2)
        
        MortalityRate=(1.- OxygenLimitation)*Consumer%BETA20*TemperatureFDecay
         ! 1/day
        
        !
        ! 
        if (Consumer%Grazing%Cohesivesed) then
        
        ! if the user provided the maximum growth rate as a filtration rate, the units will be:
        ! IngestionRate = m3/day/KgC =      m3/day /KgC    *  [-]      
        ! if the user provided the maximum growth rate as a grazing rate, the units will be:
        ! 1/day/=      1/day     *  [-] 
        
        SedimentLimitation=max(0., 1. - Me%ExternalVar%Sediment(Index)/Consumer%Grazing%MaxSediment)
                                                                                                                   
        IngestionRate=Consumer%GRMAX *  TemperatureDependence * &
                                          OxygenLimitation     * &
                                          SedimentLimitation * DensityLimitation
                                          
        else 
        
        IngestionRate=Consumer%GRMAX *  TemperatureDependence* &
                                          OxygenLimitation * DensityLimitation
        endif                                  
        
        !  1/day =           l/day     *[-]                                                                     
        RespirationRate=Consumer%RESP20*TemperatureFDecay
        
       
      !     KgC    =    1/day      * KgC                                  * day
        RespirationC=RespirationRate*Me%ExternalVar%Mass(Consumer_C, Index)* Me%DTDay 
        RespirationN=RespirationC*Consumer%NCRatio
        RespirationP=RespirationC*Consumer%PCRatio
        
        !   KgC   =     1/day     * KgC                                   *day
        MortalityC = MortalityRate*Me%ExternalVar%Mass(Consumer_C, Index) * Me%DTDay 
        MortalityN = MortalityC*Consumer%NCRatio 
        MortalityP = MortalityC*Consumer%PCRatio 

        IngestionC_tot=0.
        IngestionN_tot=0.
        IngestionP_tot=0.
       
        EgestionC_tot=0.
        EgestionN_tot=0.
        EgestionP_tot=0.
       
       
       Food => Consumer%Grazing%FirstFood
       
       
  
  
       d2:   do while(associated(Food))

          ! In this part of the code, the properties which are sources of food are identified on the basis of :
          ! 1) the name
          ! 2) the element (N,C,P) : if coupled with module waterquality, the particulate organic matter
          !    is expressed only as N and P, and the phytoplankton is expressed only as C. 
          !    If coupled with module Life, POM and phytoplankton will be expressed in terms of C,N, and P
          !    This could create confusion when expressing the grazing over the food when a different model is used
          !    for water column biogeochemistry
          !    So it is necessary to identify if the Food is expressed as C, N, P , or only one or two of 
          !    these three elements. 
          !    The keywords USE_CARBON, USE_NITROGEN and USE_PHOSPHORUS are meant to enable this identification
          !
          if ((Food%ID%Name=='phytoplankton').AND.(Me%PelagicModel == WaterQualityModel)) then
               
               FoodIndexC  = SearchPropIndex(GetPropertyIDNumber(trim(Food%ID%Name)))
          
          else if ((Food%ID%Name=='particulate organic').AND.(Me%PelagicModel == WaterQualityModel)) then
               
               if(Food%Use_Nitrogen) then
               FoodIndexN  = SearchPropIndex(GetPropertyIDNumber(trim(Food%ID%Name)//" nitrogen"))
               endif
               if(Food%Use_Phosphorus) then
               FoodIndexP  = SearchPropIndex(GetPropertyIDNumber(trim(Food%ID%Name)//" phosphorus"))
               endif
         else     
               
             if(Food%Use_Carbon)  then
             
                  FoodIndexC  = SearchPropIndex(GetPropertyIDNumber(trim(Food%ID%Name)//" carbon"))
             endif
                  
             if(Food%Use_Nitrogen)  then
                  FoodIndexN  = SearchPropIndex(GetPropertyIDNumber(trim(Food%ID%Name)//" nitrogen"))
             endif     
             
             if(Food%Use_Phosphorus)  then
                  FoodIndexP  = SearchPropIndex(GetPropertyIDNumber(trim(Food%ID%Name)//" phosphorus"))
             endif
         
         endif
          
           
        
                         
        if(Food%Use_Nitrogen)  then
        
        ! In this part of the code, more information is used for the identification of the source of food
        ! In MOhid, at the interface sediment-water, the same property (with the same ID and the same name) 
        ! can be on the bottom surface or in the water. 
        ! The benthic feeder needs to know if it will feed on food that is on the bottom surface or in the water
        ! the keyword PARTICULATEWATERFOOD  provides such information 
        ! (1=food is in water, 0 = food is on the bottom surface)
            if((Food%ID%Name=='particulate organic') .and. (Food%ParticulateWaterFood)) then
                ! the food is eaten  from the water (filter feeder) 
                ! In this case Ingestion rate is a filtration rate given as m3/Kg C /day  (or liters/g C/day)
                !                                              (divide by NC ratio to obtain right units)
 
            Me%Matrix(Index, FoodIndexN, Consumer_N) = (IngestionRate/Food%NCRatio)                         *  &  ! m3/Kg N /day 
                                                       Me%ExternalVar%Mass(Consumer_N, Index)               *  &  ! Kg N
                                                      (Me%ExternalVar%MassInKgFromWater(FoodIndexN, Index)  /  &  ! Kg N
                                                         Me%ExternalVar%WaterVolume(index))                 *  &  ! m3
                                                         Me%DTDay                                                 ! day                                            
                                                           
                Me%ExternalVar%MassInKgFromWater(FoodIndexN, Index)= Me%ExternalVar%MassInKgFromWater(FoodIndexN, Index)      - &
                                                                     Me%Matrix(Index, FoodIndexN, Consumer_N)
                else
                ! the food is eaten  from the bottom (deposit feeder)
                ! In this case the ingestion rate is in 1/day
              
               if(Food%Use_Carbon)  then
               
               FoodDens=Me%ExternalVar%Mass(FoodIndexC, Index)/ Me%ExternalVar%CellArea(Index)
               else
               FoodDens=Me%ExternalVar%Mass(FoodIndexN, Index)/ Me%ExternalVar%CellArea(Index)/Food%NCRatio
               endif
               
               
                if(FoodDens < Food%Minval) then
                Me%Matrix(Index, FoodIndexN, Consumer_N) =0.

                else
               
               
                ! Kg N = 1/day * Kg N * day * [-]
               Me%Matrix(Index, FoodIndexN, Consumer_N) = IngestionRate*Me%ExternalVar%Mass(Consumer_C, Index)        *   &  ! 1/day * Kg N
                                                          Consumer%NCRatio                                            *   &
                                                          Me%DTDay                                                    *   &  ! day
                                                          FoodDens / ( FoodDens + Consumer%Grazing%KFoodC)                   ! [-]
                endif
                
                
                Me%ExternalVar%Mass(FoodIndexN, Index)  = Me%ExternalVar%Mass(FoodIndexN, Index)- &
                                                           Me%Matrix(Index, FoodIndexN, Consumer_N)  ! kgN
                                            
                

            endif
        
            IngestionN=Consumer%Grazing%Ass_Efic*Me%Matrix(Index, FoodIndexN, Consumer_N)  ! kgN
            EgestionN=(1.-Consumer%Grazing%Ass_Efic)*Me%Matrix(Index, FoodIndexN, Consumer_N) ! ! kgN
        
        endif
        
        
            if(Food%Use_Phosphorus)  then
        
            if((Food%ID%Name=='particulate organic') .and. (Food%ParticulateWaterFood)) then
                ! food eaten from the water
                 Me%Matrix(Index, FoodIndexP, Consumer_P) =(IngestionRate/Food%PCRatio)                         * &
                                                         Me%ExternalVar%Mass(Consumer_P, Index)                * &
                                                        (Me%ExternalVar%MassInKgFromWater(FoodIndexP, Index)  / &
                                                         Me%ExternalVar%WaterVolume(index))* Me%DTDay       ! kgP
                
                Me%ExternalVar%MassInKgFromWater(FoodIndexP, Index)= Me%ExternalVar%MassInKgFromWater(FoodIndexP, Index)- &
                                                                     Me%Matrix(Index, FoodIndexP, Consumer_P)
                else
               ! the food is eaten  from the bottom (deposit feeder)
               ! Consumer%Grazing%KFoodC has units of kg C /m2
               ! the food limitation is calculated by using food biomass (FoodDens, in kg C/m2) 
                
               if(Food%Use_Carbon)  then
               
               FoodDens=Me%ExternalVar%Mass(FoodIndexC, Index)/ Me%ExternalVar%CellArea(Index)
               else
               FoodDens=Me%ExternalVar%Mass(FoodIndexP, Index)/ Me%ExternalVar%CellArea(Index)/Food%PCRatio
               endif
                
                
                ! for deposit feeders, the Michaelis-Menten kinetic was used
                ! IngestionRate has dimension of 1/day
                
                if(FoodDens < Food%Minval) then
                ! if the 
                Me%Matrix(Index, FoodIndexP, Consumer_P) = 0.
                
                else
                
                ! Kg P
                Me%Matrix(Index, FoodIndexP, Consumer_P) = IngestionRate*Me%ExternalVar%Mass(Consumer_C, Index)     *   &  ! 1/day * Kg C
                                                           Consumer%PCRatio                                         *   &
                                                           Me%DTDay                                                 *   &  ! day
                                                           FoodDens /( FoodDens + Consumer%Grazing%KFoodC)                ! [-]
                endif
                
                Me%ExternalVar%Mass(FoodIndexP, Index)=    Me%ExternalVar%Mass(FoodIndexP, Index)                   -   &
                                                           Me%Matrix(Index, FoodIndexP, Consumer_P)

            
            endif
            
            IngestionP =  Consumer%Grazing%Ass_Efic      * Me%Matrix(Index, FoodIndexP, Consumer_P)
            EgestionP  = (1.- Consumer%Grazing%Ass_Efic) * Me%Matrix(Index, FoodIndexP, Consumer_P)
        endif
        
        if(Food%Use_Carbon)  then
        
            if(Food%ParticulateWaterFood)then
                ! the food is eaten  from the water. 
                
                
                Me%Matrix(Index, FoodIndexC, Consumer_C) = IngestionRate*Me%ExternalVar%Mass(Consumer_C, Index)  * &
                                                          (Me%ExternalVar%MassInKgFromWater(FoodIndexC, Index)    / &
                                                          Me%ExternalVar%WaterVolume(index))* Me%DTDay  
                
                Me%ExternalVar%MassInKgFromWater(FoodIndexC, Index)= Me%ExternalVar%MassInKgFromWater(FoodIndexC, Index)- &
                                                                    Me%Matrix(Index, FoodIndexC, Consumer_C)
                else
                ! the food is eaten  from the bottom (deposit feeder)
                ! Consumer%Grazing%KFoodC has units of kg C /m2
                
                FoodDens=Me%ExternalVar%Mass(FoodIndexC, Index)/ Me%ExternalVar%CellArea(Index)
              
                if(FoodDens < Food%Minval) then
                    Me%Matrix(Index, FoodIndexC, Consumer_C) =0.
                    
                else
                ! for deposit feeders, the Michaelis-Menten kinetic was used
                ! IngestionRate has dimension of 1/day
                ! Kg C
                Me%Matrix(Index, FoodIndexC, Consumer_C) =IngestionRate*Me%ExternalVar%Mass(Consumer_C, Index)      *  & ! 1/day * Kg C
                                                          Me%DTDay                                                  *  & ! day
                                                          FoodDens/ ( FoodDens + Consumer%Grazing%KFoodC)                ! [-]
                                                        
                endif
                
                
                Me%ExternalVar%Mass(FoodIndexC, Index)= Me%ExternalVar%Mass(FoodIndexC, Index)                      - &
                                                        Me%Matrix(Index, FoodIndexC, Consumer_C)

            endif
        
            IngestionC=Consumer%Grazing%Ass_Efic*Me%Matrix(Index, FoodIndexC, Consumer_C)
            EgestionC=(1.-Consumer%Grazing%Ass_Efic)*Me%Matrix(Index, FoodIndexC, Consumer_C)
        
        endif
        
    
       
       if ((Food%ID%Name=='phytoplankton').AND.(Me%PelagicModel == WaterQualityModel)) then
           IngestionN=IngestionC*Consumer%NCRatio 
           IngestionP=IngestionC*Consumer%PCRatio
           EgestionN=EgestionC*Consumer%NCRatio
           EgestionP=EgestionC*Consumer%PCRatio
       endif
       
       
       if ((Food%ID%Name=='particulate organic').AND.(Me%PelagicModel == WaterQualityModel)) then
       
            if(Food%Use_Nitrogen) then
                 ! particulate organic nitrogen
               IngestionC=IngestionN/Consumer%NCRatio
               EgestionC=EgestionN/Consumer%NCRatio
               IngestionP=IngestionN*(Consumer%PCRatio/Consumer%NCRatio)
               EgestionP=EgestionN*(Consumer%PCRatio/Consumer%NCRatio)
           endif
           
           if(Food%Use_Phosphorus) then
               ! particulate organic phosphorus
               IngestionC=IngestionP/Food%PCRatio
               EgestionC=EgestionP/Food%PCRatio
               
               IngestionN=IngestionP*(Food%NCRatio/Food%PCRatio)
               EgestionN=EgestionP*(Food%NCRatio/Food%PCRatio)
           endif
         endif
         

        
           IngestionC_tot=IngestionC_tot+IngestionC   ! Ingestion of each Food is added to the total. 
           IngestionN_tot=IngestionN_tot+IngestionN
           IngestionP_tot=IngestionP_tot+IngestionP
           
           
           EgestionC_tot=EgestionC_tot+EgestionC      ! egestion of each Food is added to the total.
           EgestionN_tot=EgestionN_tot+EgestionN
           EgestionP_tot=EgestionP_tot+EgestionP
           
     
       
       Food => Food%Next
       
        end do d2
            
        Me%ExternalVar%Mass(Consumer_C, Index) = Me%ExternalVar%Mass(Consumer_C, Index)+   IngestionC_tot &
                                                                                       -   MortalityC     &
                                                                                       -   RespirationC

        Me%ExternalVar%Mass(Consumer_N, Index) = Me%ExternalVar%Mass(Consumer_C, Index)* Consumer%NCRatio &
                                                                                       +   IngestionN_tot &
                                                                                       -   MortalityN     &
                                                                                       -   RespirationN
                                                                                       
        Me%ExternalVar%Mass(Consumer_P, Index) = Me%ExternalVar%Mass(Consumer_C, Index)* Consumer%PCRatio &
                                                                                       + IngestionP_tot   &
                                                                                       - MortalityP       &
                                                                                       - RespirationP
                   
      
     
    
      Me%ExternalVar%Mass(PON, Index)=Me%ExternalVar%Mass(PON, Index)+EgestionN_tot + &
                                     MortalityN 
      ! PON e POP are Foods, so the predation was calculated already before
      Me%ExternalVar%Mass(POP, Index)=Me%ExternalVar%Mass(POP, Index)+EgestionP_tot + &
                                      MortalityP
      if (Me%PelagicModel == LifeModel) then
        Me%ExternalVar%Mass(POC, Index)=Me%ExternalVar%Mass(POC, Index)+EgestionC_tot + MortalityC
      endif
      
      Me%ExternalVar%MassInKgFromWater(O2, Index)=Me%ExternalVar%MassInKgFromWater(O2, Index)-RespirationC*32./12.
      Me%ExternalVar%MassInKgFromWater(AM, Index)=Me%ExternalVar%MassInKgFromWater(AM, Index)+RespirationN
      Me%ExternalVar%MassInKgFromWater(IP, Index)=Me%ExternalVar%MassInKgFromWater(IP, Index)+RespirationP
      
   
      
      
      Consumer => Consumer%Next
        end do d1

        
        end subroutine ComputeBenthicConsumers
    
    
    !-----------------------------------------------------------------------------------
    
    
    
    
    !---------------------------------------------------------------------------------------
    
    subroutine ComputeBenthicBivalveDEB (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_BivalveDEB),      pointer           :: BivalveDEB
        integer                                    :: PON
        integer                                    :: Phyto
        integer                                    :: AM
        integer                                    :: O2
        
        ! Add to comply to MOHID structure
        integer                                    :: BivalveDEB_V
        integer                                    :: BivalveDEB_E
        integer                                    :: BivalveDEB_ER
        integer                                    :: BivalveDEB_IDW
        integer                                    :: BivalveDEB_TDW
        integer                                    :: NInd
        
        real                                       :: ConcPhy
        real                                       :: ConcChla                       ! Available food in in mg (chl-a) m -3
        real                                       :: consump_O2                     ! O2 consuption rate in kg O2
        
        ! DEB model parameters
        real                                       :: Pxm_T,Pam_T,Pm_T,Jw_T          ! Temperature adjusted values of Pxm, Pam, Pm and Jw
        real                                       :: fresp                          ! Holling II Organism response function to available food
        real                                       :: Px                             ! Energy ingestion rate in J s-1
        real                                       :: Pa                             ! The assimilation rate in J s-1
        real                                       :: Pc                             ! Energy utilisation rate i.e. energy both fixed and
                                                                                     ! dissipated consumed by body tissues from reserve in J s-1
        real                                       :: Pj                             ! Energy flow to maturity maintanence in J s-1
        real                                       :: Pr                             ! Energy flow to reproductive buffer in J s-1
        real                                       :: shrink_Pm, shrink_Pj           ! Shrinking (maintanence, maturity, 
        real                                       :: shrink_Er, shrink_V, shrink    ! repro, structure, total) in J s-1
        real                                       :: dV                             ! Change in the volume of bivalve in m-3 s-1
        real                                       :: dE                             ! Change in the energy reserve in J s-1
        real                                       :: dEr                            ! Change in reproductive energy store in J s-1
        real                                       :: bi_area_scale                  ! Bivalve area scale in m2 - NOT ACTUAL SHELL AREA
        real                                       :: Pam                            ! Max. surface area specific assimilation rate in J m-2 s-1
        real                                       :: energy_surplus                 ! Energy which goes to structural growth  in J s-1
        real                                       :: e_dens                         ! Energy density in J m-3
        real                                       :: cff1,cff2                      ! temporary variables
        real                                       :: Tresp,Ting                     ! temperature correction factor
        real                                       :: bi_len        
        integer                                    :: seq
 
 
        !output variables in the original code
        real                                       :: tot_DW                ! total DW in g,
        real                                       :: Er_DW_perc            ! repro energy percentage of DW
        real                                       :: bi_consump_chla       ! chlorophyll_a consumption rate in mg chla s-1
        real                                       :: Nexcrt                ! N-NH4 excretion rate mmol N s-1 
        real                                       :: Ceg, Neg              ! egested C and N in faeces in mmolC s-1 and mmolN s-1
        real                                       :: Cresp                 ! respired C in mmol C s-1                               
        real                                       :: Cs                    ! Fecal contamination in bivalves in CFU (or MPN)  
        real                                       :: bi_consump_oxy_mass   ! oxygen consumption rate in mmolO2 s-1
        
      
        
            
            !********* OXYGEN ***********************
        real                                       :: bi_consump_oxy                 ! oxygen consumption rate in mgO2 s-1
        
            !********* FAECES PRODUCTION ***************
        real                                       :: Cing                           ! ingested C in mmol C s-1
        real                                       :: Ning                           ! ingested N in mmol N s-1     
        real                                       :: Cst, Nst, Ce, Ne, Cer, Ner     !all in mmol s-1
        real                                       :: Cflesh, Nflesh, Csh, Nsh       !all in mmol s-1

           !********* PATHOGEN ************************
        !real                                       :: ku, Csdep           ! filtration rate (m3 s-1 g-1), new CFU in bivalves
        !real                                       :: bi_consump_CH2O     ! E coli consumption rate in CFU s-1 (per individual mussel)
        
        
        BivalveDEB => Me%FirstBivalveDEB
        
        PON     = Me%PropIndex%PON
        AM      = Me%PropIndex%Ammonia
        Phyto   = Me%PropIndex%Phyto
        O2      = Me%PropIndex%Oxygen
        
        
        
        !Calculate concentration from mass and volume in each cell
        ConcPhy = Me%ExternalVar%MassInKgFromWater(Phyto, index) / Me%ExternalVar%WaterVolume(Index)  ! gC/l
        
        ConcChla = ConcPhy * (1. / BivalveDEB%CtoChla) * 10.0e6                                        ! mg (chl-a) m -3
        
       
            
!   bi_vol = (NewBivalveDEB%BiLen1 * NewBivalveDEB%delm)**3.0   !volume in m3 (0.05 is L in m)
!   bi_energy = 0.8 * NewBivalveDEB%EM * bi_vol           
               
        
        
        

d1:     do while(associated(BivalveDEB))

        !State-variables for the DEB model
        BivalveDEB_V   = BivalveDEB%BiSV%Bivol
        
        BivalveDEB_E   = BivalveDEB%BiSV%Bivenerg
        
        BivalveDEB_ER  = BivalveDEB%BiSV%Bivenerep
        
        BivalveDEB_IDW = BivalveDEB%BiSV%Bivinddw
        
        BivalveDEB_TDW = BivalveDEB%BiSV%Bivtotaldw  



    !     % Temperature correction factors
    Tresp = exp(BivalveDEB%Ta / BivalveDEB%Tref - BivalveDEB%Ta / (273.15 + Me%ExternalVar%Temperature(index))) /           &
            (1.0 + exp(BivalveDEB%TAL / (273.15 + Me%ExternalVar%Temperature(index)) - BivalveDEB%TAL / BivalveDEB%TL) +    &
            exp(BivalveDEB%TAHresp / BivalveDEB%THresp - BivalveDEB%TAHresp / (273.15 + Me%ExternalVar%Temperature(index))))

    Ting = exp(BivalveDEB%Ta / BivalveDEB%Tref - BivalveDEB%Ta / (273.15 + Me%ExternalVar%Temperature(index))) /            &
            (1.0 + exp(BivalveDEB%TAL / (273.15 + Me%ExternalVar%Temperature(index)) - BivalveDEB%TAL / BivalveDEB%TL) +    &
            exp(BivalveDEB%THing / BivalveDEB%THing - BivalveDEB%THing / (273.15 + Me%ExternalVar%Temperature(index))))




    ! Update rates for temperature:
    
    Pxm_T = BivalveDEB%Pxm * Ting
    
    Pm_T  = BivalveDEB%Pm * Tresp


     fresp = ConcChla/(BivalveDEB%Xk + ConcChla)
          

     bi_area_scale = (Me%ExternalVar%Mass(BivalveDEB_V, Index) / Me%ExternalVar%CellArea(Index)) ** (2.0/3.0)  

     Px = Pxm_T * fresp * bi_area_scale

     Pam_T = BivalveDEB%ka * Pxm_T
     
     Pa = BivalveDEB%ka * Px       ! assimilation flux




    ! Energy density:
    
     e_dens = (Me%ExternalVar%Mass(BivalveDEB_E, Index)/ Me%ExternalVar%CellArea(Index)) / &
               (Me%ExternalVar%Mass(BivalveDEB_V, Index)/ Me%ExternalVar%CellArea(Index))

       if ((Me%ExternalVar%Mass(BivalveDEB_E, Index)/ Me%ExternalVar%CellArea(Index)) .eq. 0.0) then
         
             Pc = 0.0
         
           else
           
             cff1 = (BivalveDEB%Eg * Pam_T * bi_area_scale / BivalveDEB%Em) + Pm_T * &
                     (Me%ExternalVar%Mass(BivalveDEB_V, Index)/ Me%ExternalVar%CellArea(Index))
                       
             cff2 = e_dens/(BivalveDEB%Eg + BivalveDEB%kappa * e_dens)
             
             Pc = cff2 * cff1
       
       endif

   energy_surplus = max(0.0 , BivalveDEB%kappa * Pc - Pm_T * (Me%ExternalVar%Mass(BivalveDEB_V, Index)/ &
                     Me%ExternalVar%CellArea(Index)))  ! (J/s) flow to structural growth



    ! % maturity maintenance and reproduction :
    
        if ((Me%ExternalVar%Mass(BivalveDEB_V, Index)/ Me%ExternalVar%CellArea(Index)) .LT. BivalveDEB%volp) then
        
            Pj = Me%ExternalVar%Mass(BivalveDEB_V, Index) / Me%ExternalVar%CellArea(Index)* &
                 ((1.0-BivalveDEB%kappa) / BivalveDEB%kappa) * Pm_T  !Flow to maturity maintanence
            
        else
            
                Pj = BivalveDEB%volp * ((1.0 - BivalveDEB%kappa) / BivalveDEB%kappa) * Pm_T    !Flow to maturity maintenance
                
                Pr = max(0.0 , (1.0 - BivalveDEB%kappa) * Pc - Pj)       !Allocation to reprod buffer
        
        endif




    ! % Lysis and shrinking :

    shrink_Pm = max(0.0, (Pm_T * (Me%ExternalVar%Mass(BivalveDEB_V, Index) / Me%ExternalVar%CellArea(Index)) - &
                   BivalveDEB%kappa * Pc))  ! to pay maint. struct.
    
    shrink_Pj = max(0.0, (Pj - (1.0 - BivalveDEB%kappa) * Pc))                    ! to pay maint. maturity
    
    shrink = shrink_Pm + shrink_Pj


       if (shrink .gt. 0.0) then                 !Take energy from repro buffer to pay maintenance costs
         
         shrink_Er = shrink / BivalveDEB%kR      !pay from reproductive buffer
         
         shrink_V = 0.0                          !don't pay from structure
       
       else
             
             shrink_Er = 0.0
             
             shrink_V = 0.0
       
       endif

    dE = Pa - Pc
    
    dV = energy_surplus / BivalveDEB%Eg - shrink_V 
    
    dEr = (Pr * BivalveDEB%kR) - shrink_Er



    ! Calculate total DW in g

    tot_DW = (Me%ExternalVar%Mass(BivalveDEB_V, Index) / Me%ExternalVar%CellArea(Index))* BivalveDEB%St_DW_perc * 1.0e6 + &
             (Me%ExternalVar%Mass(BivalveDEB_E, Index)/ Me%ExternalVar%CellArea(Index) + &
              Me%ExternalVar%Mass(BivalveDEB_ER, Index)/ Me%ExternalVar%CellArea(Index)) / BivalveDEB%Conv_fac   ! *1e6 is for converting bi_vol from m3 to cm3

    ! Calculate individusl DW in kg
    Me%ExternalVar%Mass(BivalveDEB_IDW, Index) = (tot_DW/1000) * Me%ExternalVar%CellArea(Index)
    ! Total DW in Kg/m2
    Me%ExternalVar%Mass(BivalveDEB_TDW, Index) = Me%ExternalVar%Mass(BivalveDEB_IDW, Index) * (BivalveDEB%NInd) !LLP
    ! SPAWNING

    ! % Gono-somatic index
    
    Er_DW_perc = ((Me%ExternalVar%Mass(BivalveDEB_ER, Index) / Me%ExternalVar%CellArea(Index)) / BivalveDEB%Conv_fac) / &
                      tot_DW
    


    ! NOTE Jday is a global variable and must be calculated before the call to the subroutine
    ! Jday is the Julian day of the current year (i.e. referred to start of the current year)
   

    !  Integrate the state variables
 
        ! conversion of DT from days (DEB) to s (MOHID)
        Me%ExternalVar%Mass(BivalveDEB_V, Index) = (Me%ExternalVar%Mass(BivalveDEB_V, Index) + &
                                                    ((dV * Me%ExternalVar%CellArea(Index))* Me%DT)) !LLP
       

        Me%ExternalVar%Mass(BivalveDEB_E, Index) = max(0.0 , (Me%ExternalVar%Mass(BivalveDEB_E, Index) + &
                                                    ((dE * Me%ExternalVar%CellArea(Index))* Me%DT))) !LLP max in order to avoid negative values

        Me%ExternalVar%Mass(BivalveDEB_ER, Index) = max(0.0 ,(Me%ExternalVar%Mass(BivalveDEB_ER, Index) + &
                                                     ((dEr * Me%ExternalVar%CellArea(Index))* Me%DT))) !LLP
                                                     
    ! test LLP
        if (Me%ExternalVar%Mass(BivalveDEB_ER, Index) .LT. 0) then
            write(*,*)'ER negativo1', Index, Me%ExternalVar%Mass(BivalveDEB_ER, Index), dEr, Pr, BivalveDEB%kR, shrink_Er
            stop
        endif
            


     !Update of BivalveDEB_ER if spawning occurs
d2:  do seq = 1, BivalveDEB%spawnD
     
         if (BivalveDEB%spawnDay(seq) .GT. 0) then
         
            if ((BivalveDEB%spawnDay(seq) .EQ. Me%JulianDay) .AND. (Er_DW_perc .GT. BivalveDEB%RGS) .AND.   &
                     (Me%ExternalVar%Temperature(index) .GT. BivalveDEB%Tspawn)) then
            
                    Me%ExternalVar%Mass(BivalveDEB_ER, Index) = max(0.0, Me%ExternalVar%Mass(BivalveDEB_ER, Index) - &
                                                                Me%ExternalVar%Mass(BivalveDEB_ER, Index) * &
                                                                BivalveDEB%spawn_eff) !LLP add max(0.0,...) to avoid negative values
                    
                    tot_DW = (Me%ExternalVar%Mass(BivalveDEB_V, Index) / Me%ExternalVar%CellArea(Index))* BivalveDEB%St_DW_perc * 1.0e6 + &
                             (Me%ExternalVar%Mass(BivalveDEB_E, Index)/ Me%ExternalVar%CellArea(Index) + &
                              Me%ExternalVar%Mass(BivalveDEB_ER, Index)/ Me%ExternalVar%CellArea(Index)) / BivalveDEB%Conv_fac!LLP
                              
                    Me%ExternalVar%Mass(BivalveDEB_IDW, Index) = (tot_DW*1000) * Me%ExternalVar%CellArea(Index)
                    
                    Me%ExternalVar%Mass(BivalveDEB_TDW, Index) = Me%ExternalVar%Mass(BivalveDEB_IDW, Index) * (BivalveDEB%NInd) !LLP
            
            if (Me%ExternalVar%Mass(BivalveDEB_ER, Index) .LT. 0.0) then
              write(*,*)'ER negativo2', Index, Me%ExternalVar%Mass(BivalveDEB_ER, Index),dEr, Me%DT, Me%ExternalVar%CellArea(Index)
              stop
            endif
         
            endif
                 
         endif
    
     end do d2

  !_________________________________THIS IS THE END OF GROWTH & BIOENERGETICS 




  ! _________________________________CALCULATE FEEDBACKS 
  
  
   !CALCULATE CHLA CONSUMPTION 
   
   !convert chla consumption rate from mg chla s-1 (energy based approach) to kg chla s-1
   bi_consump_chla = (Px / (BivalveDEB%Ec * BivalveDEB%CtoChla)) / 10.0e6           
  
    
  ! OXYGEN UPTAKE 
  ! Bourles et al. (2009) (IFREMER's  model) relates O2 consumption to the catabolic energy flux, Pc. Therefore:
   
    bi_consump_oxy = (Pc / BivalveDEB%etaO2)     ! in mgO2 s-1
    
  !Convert from mg O2 s-1 to kg O2 s-1 used by MOHID
    
    bi_consump_oxy_mass = (bi_consump_oxy / 10.0e6) 



  ! INGESTION OF FOOD *********************
  ! Ingested C & N:
  
   Cing = Px / (BivalveDEB%Ec * 12.0)   !ingested C in mmol C s-1
  
   Ning = Cing / 6.625                  !ingested N in mmol N s-1



  ! EGESTION OF FEACES ********************
  ! Egested C & N
  
   Ceg = (1.0 - BivalveDEB%ka) * Cing     !egested C in mmol C s-1
  
   Neg = (1.0 - BivalveDEB%ka) * Ning     !egested N in mmol N s-1


    ! Calculate assimilation for M.edulis, composition of other species unknown
    ! For M.edulis NH4 excretion is then calculated from mass balance equation
    ! For other species NH4 excreted is set as a fraction of N ingested
    ! In MOHID this option is activated by using a KEYWORD 
    !     (NH4EXCFRAC = 1 then NH4 is excreted as a fraction of N ingested) 
  
  
  
  
   if(BivalveDEB%NH4Frac .eq. 0)then

    !******* ASSIMILATION OF MATTER ****************
    ! C & N in structural tissue (always >= 0)
   
       Cst = 0.4 * energy_surplus * 1000.0 / (23000.0 * 12.0)        !stuctural tissue in mmol C s-1
       
       Nst = 0.4 * energy_surplus * 1000.0 / (23000.0 * 4.82 * 14.0) !structural tissue in mmol N s-1


    ! C & N in energy reserves
       
       Ce = 0.4 * dE * 1000.0 /(23000.0 * 12.0)          !energy reserves in mmol C s-1
    
       Ne = 0.4 * dE * 1000.0 /(23000.0 * 4.82 * 14.0)   !energy reserves in mmol N s-1


    ! C & N in reproductive energy
   
       Cer = 0.4 * dEr * 1000.0 /(23000.0 * 12.0)         !reproductive energy in mmol C s-1
   
       Ner = 0.4 * dEr * 1000.0 /(23000.0 * 4.82 * 14.0)  !reproductive energy in mmol N s-1


    ! Total C & N in flesh
   
       Cflesh = Cst + Ce + Cer
   
       Nflesh = Nst + Ne + Ner



    ! NOTE: If total energy increased or unchanged (i.e. Cflesh >= 0) then excess C is respired and excess N is excreted
    !       If total energy decreased, then Cresp and Nexcrt are greater than assimilated (i.e. lysis)

    ! C respired and N excreted:
        if (Cflesh.GE.0.0) then
        
             Csh = Cflesh * 0.08 / 0.92
             
             Nsh = Nflesh * 0.12 / 0.88                    !ratios after Rodhouse and Roden (1984) (Figure 7)
             
             Cresp = max(0.0, Cing - Ceg - Cflesh - Csh)
             
             Nexcrt = max(0.0, Ning - Neg - Nflesh - Nsh)  !max here just in case it goes negative
       
            else
         
                 Cresp = Cing - Ceg - Cflesh
         
                 Nexcrt = Ning - Neg - Nflesh
       
         end if

     else
    
    ! if other than M.edulis
    
    Nexcrt = 0.45 * Ning                !factor can be modified by the user

  end if

!_____________________________________________ END OF FEEDBACKS
  
  
  
  
  !___________________________________________ SHELLFISH - PATHOGEN MODEL
  
!   ku = bi_consump_chla / (Xf * bi_vol * 1e6)       ! Volume filtration rate per g wet weight and time (m3 g-1 s-1)

    ! Calculate E coli uptake rate in CFU s-1 (per individua)l
!    bi_consump_CH2O = ku * CH2O * bi_vol * 1e6

!    Csdep = (ku * CH2O - kd * Cs) * dt + Cs          ! New Cs in CFU g-1 (g is g wet weight)
 
!    if(Csdep .lt. 0.0) Csdep=0.0
        
!    Cs = Csdep
!
! ****  END OF SHELLFISH - PATHOGEN MODEL *******
  
  
  
  
 !Upadate of the properties  
    
 consump_O2 = bi_consump_oxy_mass * Me%DT * (BivalveDEB%NInd * Me%ExternalVar%CellArea(Index)) !kg O2
 if (consump_O2 .GE. Me%ExternalVar%MassInKgFromWater(O2, Index)) consump_O2 = 0.0
 
 
 
 Me%ExternalVar%MassInKgFromWater(O2, Index) = Me%ExternalVar%MassInKgFromWater(O2, Index) - &
                                               consump_O2 !LLP                                                    
 
 !converted from mmol N to kg N

 Me%ExternalVar%MassInKgFromWater(AM, Index) = Me%ExternalVar%MassInKgFromWater(AM, Index) + &
                                               (Nexcrt * 14.01 / 10.0e6) * Me%DT * &          
                                               (BivalveDEB%NInd * Me%ExternalVar%CellArea(Index)) !LLP
                                               
 !converted from kg chla s-1 to kg C
 Me%ExternalVar%MassInKgFromWater(Phyto, Index) = Me%ExternalVar%MassInKgFromWater(Phyto, Index)- &
                                                 bi_consump_chla * Me%DT * BivalveDEB%CtoChla * &
                                                 (BivalveDEB%NInd * Me%ExternalVar%CellArea(Index)) !LLP

 !converted from mmol N to kg N
 Me%ExternalVar%Mass(PON, Index) = Me%ExternalVar%Mass(PON, Index) + &
                                   (Neg * 14.01 / 10.0e6) * Me%DT * &
                                   (BivalveDEB%NInd * Me%ExternalVar%CellArea(Index)) !LLP
  
  
    !Only for output purposes in cm
    ! Recalculate length as the output from deb_bivalve_evolution is volume of the animal
    bi_len = ((Me%ExternalVar%Mass(BivalveDEB_V, Index) * Me%ExternalVar%CellArea(Index))**(1.0/3.0) / &
                     BivalveDEB%delm) * 100.0
  
    
      BivalveDEB => BivalveDEB%Next
        end do d1

        
    end subroutine ComputeBenthicBivalveDEB  
    
       
    !-----------------------------------------------------------------------------------
    
    
    
    !-----------------------------------------------------------------------------------
        
    subroutine ComputeBenthicNitrogen(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: AM, PON, O2
        integer                                     :: PON1, PON2, PON3, PON4, PON5
        real                                        :: MineralizationRate
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation

        !Begin-----------------------------------------------------------------

        AM  = Me%PropIndex%Ammonia
        PON = Me%PropIndex%PON
        O2  = Me%PropIndex%Oxygen
        

        
        if(Me%ComputeOptions%Pompools)then
            PON1 = Me%PropIndex%PON1
            PON2 = Me%PropIndex%PON2
            PON3 = Me%PropIndex%PON3
            PON4 = Me%PropIndex%PON4
            PON5 = Me%PropIndex%PON5        
        end if

        !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max((Me%ExternalVar%MassInKgFromWater(O2,Index)/ &
        Me%ExternalVar%WaterVolume(index))*1000., Me%Oxygen%Minimum)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Nitrogen%PONDecayRate       *  &
                             Me%Nitrogen%PONDecayTFactor    ** &
                            (Me%ExternalVar%Temperature(Index) - 20.0)

       
        !kgN * day * day-1 (what passes from PON to ammonia)
        Me%Matrix(Index, PON, AM) = Me%ExternalVar%Mass(PON, Index) * Me%DTDay * &
                                    MineralizationRate * OxygenLimitation

        if(Me%ComputeOptions%Pompools)then
            Me%Matrix(Index, PON1, AM) = Me%ExternalVar%Mass(PON1, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON2, AM) = Me%ExternalVar%Mass(PON2, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON3, AM) = Me%ExternalVar%Mass(PON3, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON4, AM) = Me%ExternalVar%Mass(PON4, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON5, AM) = Me%ExternalVar%Mass(PON5, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
        end if


        if(.NOT. Me%ComputeOptions%Pompools)then

            Me%ExternalVar%MassInKgFromWater(AM,  Index) = Me%ExternalVar%MassInKgFromWater(AM , Index) + Me%Matrix(Index, PON, AM)

            Me%ExternalVar%Mass(PON, Index) = Me%ExternalVar%Mass(PON, Index) - Me%Matrix(Index, PON, AM)

            !what is consumed of oxygen due to mineralization of PON
            OxygenConsumption               = Me%Matrix(Index, PON, AM) * 1. / Me%OrganicMatter%NC_Ratio * &
                                              32. / 12.

        else
         
            Me%ExternalVar%MassInKgFromWater(AM,  Index) = Me%ExternalVar%MassInKgFromWater(AM , Index) + &
                                              Me%Matrix(Index, PON, AM) + &
                                              Me%Matrix(Index, PON1, AM) + Me%Matrix(Index, PON2, AM)     + &
                                              Me%Matrix(Index, PON3, AM) + Me%Matrix(Index, PON4, AM)     + &
                                              Me%Matrix(Index, PON5, AM)  

            Me%ExternalVar%Mass(PON, Index) = Me%ExternalVar%Mass(PON, Index) - Me%Matrix(Index, PON, AM)
            
            Me%ExternalVar%Mass(PON1, Index) = Me%ExternalVar%Mass(PON1, Index) - Me%Matrix(Index, PON1, AM)
            Me%ExternalVar%Mass(PON2, Index) = Me%ExternalVar%Mass(PON2, Index) - Me%Matrix(Index, PON2, AM)
            Me%ExternalVar%Mass(PON3, Index) = Me%ExternalVar%Mass(PON3, Index) - Me%Matrix(Index, PON3, AM)
            Me%ExternalVar%Mass(PON4, Index) = Me%ExternalVar%Mass(PON4, Index) - Me%Matrix(Index, PON4, AM)
            Me%ExternalVar%Mass(PON5, Index) = Me%ExternalVar%Mass(PON5, Index) - Me%Matrix(Index, PON5, AM)

            !what is consumed of oxygen due to mineralization of PON
            OxygenConsumption               = (Me%Matrix(Index, PON, AM) + Me%Matrix(Index, PON1, AM)  + &
                                              Me%Matrix(Index, PON2, AM) + Me%Matrix(Index, PON3, AM)  + &
                                              Me%Matrix(Index, PON4, AM) + Me%Matrix(Index, PON5, AM)) * &
                                              1. / Me%OrganicMatter%NC_Ratio * 32. / 12.
        
        end if
        
        

        Me%ExternalVar%MassInKgFromWater(O2, Index ) = Me%ExternalVar%MassInKgFromWater(O2, Index ) - &
                                                       OxygenConsumption


    end subroutine ComputeBenthicNitrogen
    
    !--------------------------------------------------------------------------


    subroutine ComputeBenthicPhosphorus(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: IP, POP, O2
        integer                                     :: POP1, POP2, POP3, POP4, POP5
        real                                        :: MineralizationRate
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation

        !Begin-----------------------------------------------------------------
        
        IP  = Me%PropIndex%Phosphate
        POP = Me%PropIndex%POP
        O2  = Me%PropIndex%Oxygen
        

        
        if(Me%ComputeOptions%Pompools)then
            POP1 = Me%PropIndex%POP1
            POP2 = Me%PropIndex%POP2
            POP3 = Me%PropIndex%POP3
            POP4 = Me%PropIndex%POP4
            POP5 = Me%PropIndex%POP5        
        end if
        
        
        !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max((Me%ExternalVar%MassInKgFromWater(O2,Index)/ &
        Me%ExternalVar%WaterVolume(index))*1000., Me%Oxygen%Minimum)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Phosphorus%POPDecayRate       *  &
                             Me%Phosphorus%POPDecayTFactor    ** &
                            (Me%ExternalVar%Temperature(Index) - 20.0)


        !kgP * day * day-1 (what passes from POP to inorganic phosphorus)
        Me%Matrix(Index, POP, IP) = Me%ExternalVar%Mass(POP, Index) * Me%DTDay * &
                         MineralizationRate * OxygenLimitation
        
        
            !------------------------------------------POM POOLS                 
            if(Me%ComputeOptions%Pompools)then
            
                Me%Matrix(Index, POP1, IP) = Me%ExternalVar%Mass(POP1, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                
                Me%Matrix(Index, POP2, IP) = Me%ExternalVar%Mass(POP2, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                
                Me%Matrix(Index, POP3, IP) = Me%ExternalVar%Mass(POP3, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                
                Me%Matrix(Index, POP4, IP) = Me%ExternalVar%Mass(POP4, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                                             
                Me%Matrix(Index, POP5, IP) = Me%ExternalVar%Mass(POP5, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation                                    
            end if


        if(.NOT. Me%ComputeOptions%Pompools)then
        
            Me%ExternalVar%MassInKgFromWater(IP,  Index) = Me%ExternalVar%MassInKgFromWater(IP , Index) +  &
                                                           Me%Matrix(Index, POP, IP)

            Me%ExternalVar%Mass(POP, Index) = Me%ExternalVar%Mass(POP, Index) - Me%Matrix(Index, POP, IP)

            OxygenConsumption               = Me%Matrix(Index, POP, IP) * 1. / Me%OrganicMatter%PC_Ratio * &
                                              32. / 12.
        
        else
        
            Me%ExternalVar%MassInKgFromWater(IP,  Index) = Me%ExternalVar%MassInKgFromWater(IP , Index)   + &
                                              Me%Matrix(Index, POP, IP) +  &
                                              Me%Matrix(Index, POP1, IP) + Me%Matrix(Index, POP2, IP)     +  &
                                              Me%Matrix(Index, POP3, IP) + Me%Matrix(Index, POP4, IP)     +  &
                                              Me%Matrix(Index, POP5, IP)

            Me%ExternalVar%Mass(POP, Index) = Me%ExternalVar%Mass(POP, Index) - Me%Matrix(Index, POP, IP)
            
            Me%ExternalVar%Mass(POP1, Index) = Me%ExternalVar%Mass(POP1, Index) - Me%Matrix(Index, POP1, IP)
            Me%ExternalVar%Mass(POP2, Index) = Me%ExternalVar%Mass(POP2, Index) - Me%Matrix(Index, POP2, IP)
            Me%ExternalVar%Mass(POP3, Index) = Me%ExternalVar%Mass(POP3, Index) - Me%Matrix(Index, POP3, IP)
            Me%ExternalVar%Mass(POP4, Index) = Me%ExternalVar%Mass(POP4, Index) - Me%Matrix(Index, POP4, IP)
            Me%ExternalVar%Mass(POP5, Index) = Me%ExternalVar%Mass(POP5, Index) - Me%Matrix(Index, POP5, IP)
            

            OxygenConsumption               = (Me%Matrix(Index, POP, IP) + Me%Matrix(Index, POP1, IP)  +  &
                                              Me%Matrix(Index, POP2, IP) + Me%Matrix(Index, POP3, IP)  +  &
                                              Me%Matrix(Index, POP4, IP) + Me%Matrix(Index, POP5, IP)) *  &
                                              1. / Me%OrganicMatter%PC_Ratio * 32. / 12.

        end if
        
        !if(.NOT.Me%ComputeOptions%Nitrogen)then ! if N is activated, O2 consumption was already calculated 
        ! for PON decomposition
        Me%ExternalVar%MassInKgFromWater(O2, Index ) = Me%ExternalVar%MassInKgFromWater(O2, Index ) - &
                                                       OxygenConsumption 
        !endif                 
                         
    end subroutine ComputeBenthicPhosphorus
    
    !--------------------------------------------------------------------------


    subroutine ComputeBenthicSilica(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: BioSi, Sil

        !Begin-----------------------------------------------------------------


        
        Sil   = Me%PropIndex%DissolvedSilica
        BioSi = Me%PropIndex%BioSilica

        !kg * day * day-1 (what passes from biogenic silica to inorganic dissolved silica)
        Me%Matrix(Index, BioSi, Sil) = Me%ExternalVar%Mass(BioSi, Index) * Me%DTDay * Me%Silica%BioSiDecayRate

        Me%ExternalVar%Mass(Sil, Index)   = Me%ExternalVar%Mass(Sil,   Index) + Me%Matrix(Index, BioSi, Sil)

        Me%ExternalVar%Mass(BioSi, Index) = Me%ExternalVar%Mass(BioSi, Index) - Me%Matrix(Index, BioSi, Sil)

    
    end subroutine ComputeBenthicSilica
    
    !-------------------------------------------------------------------------------------
    
     subroutine ComputeSeagrasses(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: PON, POP, L,R, NINT, PINT
        integer                                     :: PO4,AM,O2
        integer                                     :: Zone
        real                                        :: LeavesMortality
        real                                        :: LeavesBioMass
        real                                        :: RootsMortality
        real                                        :: TemperatureLim
        real                                        :: SpaceLimitation
        real                                        :: InternalNitrogen
        real                                        :: InternalNitrogenFunction
        real                                        :: InternalPhosphorus
        real                                        :: InternalPhosphorusFunction
        real                                        :: NitrogenLimitation
        real                                        :: PhosphorusLimitation
        real                                        :: NutrientsLimitation
        real                                        :: UptakeNH4s
        real                                        :: UptakePO4s
        real                                        :: UptakeNH4NO3w
        real                                        :: UptakePO4w
        real                                        :: UptakeN
        real                                        :: UptakeP
        real                                        :: Growth
        real                                        :: LightFactor
        integer, parameter                          :: NoLimitation = 1
        integer, parameter                          :: Erosion      = 2
        real                                        :: TotalPlantMass
        real                                        :: phi,Daylight,fm
        real, parameter                             :: pi = 3.14159
        ! Description
        
        real    :: s1, s2, xa, xb, ya, yb
        ! This subroutine calculates the Seagrasses sources and sinks terms
        ! Seagrasses are divided into Leaves and Roots. 
        !Begin-----------------------------------------------------------------


        
        
        O2    = Me%PropIndex%Oxygen   ! Index of the property Oxygen
        
        L     = Me%PropIndex%Leaves   ! Index of the property Leaves
        R     = Me%PropIndex%Roots    ! Index of the property Roots
        
        
        
        
        if(Me%ComputeOptions%Nitrogen)then
        PON   = Me%PropIndex%PON     ! Index of the property Particulate Organic Nitrogen
        AM    = Me%PropIndex%Ammonia ! Index of the property Ammonia
        NINT  = Me%PropIndex%Nint    ! Index of the property Seagrasses Internal Nitrogen Content
        endif
        
        if(Me%ComputeOptions%Phosphorus)then 
        POP   = Me%PropIndex%POP     ! Index of the property Particulate Organic Phosphorus
        PO4    = Me%PropIndex%Phosphate ! Index of the property Phosphate
        PINT  = Me%PropIndex%Pint    ! Index of the property Seagrasses Internal Phosphorus Content
        endif
        
        ! Calculate if there is erosion:
        
          if(Me%ExternalVar%ShearStress(Index)       > Me%Leaves%ErosCritShear     )then
                    
                    Zone = Erosion
                else
                    Zone = NoLimitation
                endif
                
                
        ! Initialize
          InternalNitrogenFunction =0.
          InternalPhosphorusFunction=0.
          RootsMortality =0.
        
        select case(Zone)
        
        case(NoLimitation)
        
            select case (Me%Leaves%MortType)
            case(1)
            !
            ! Leaves and roots decay rate expressed as a function of the daylight, 
            ! Calculated as described in Forsythe WC, Rykiel Jr EJ, Stahl RS, Wu H, Schoolfield RM.
            ! A model comparison for daylength as a function of latitude and day of year.
            ! Ecological Modeling 1995; 80:87-95.
            
            phi=asin(.39795*cos(.2163108 + 2*atan(.9671396*tan(.00860*(Me%JulianDay-186)))))
            Daylight=24-(24/pi)*acos((sin(0.8333*pi/180) + sin(Me%Leaves%Lat*pi/180)*sin(phi))/ &
                      (cos(Me%Leaves%Lat*pi/180)*cos(phi)))
            
            
            fm=exp(9.5-Daylight)
            
            case(2)
            !Leaves and roots decay expressed as a function of temperature (tentative)
            
            
            
            fm=min (1., exp(15. - Me%ExternalVar%Temperature(index)))
            end select
        
        ! Leaves uptake nutrients from the water column.
        ! This uptake is not calculated here in the module BenthicEcology,
        ! but in the module SeagrassWaterInteraction, where uptake of nutrients
        ! is stored in the array Me%Rates (See module SeagrassWaterInteraction),
        ! and then it is passed to the module BenthicEcology through the module
        ! InterfaceSedimentWater.
        ! Me%ExternalVar%UptakeNH4NO3w = sum of ammonia and nitrate uptake by leaves
        ! Me%ExternalVar%UptakePO4w    = uptake of phosphate by leaves
        !
        ! gN/day

        UptakeNH4NO3w= Me%ExternalVar%UptakeNH4NO3w(index)
        ! gP/day
        UptakePO4w= Me%ExternalVar%UptakePO4w(index)
        
        
        ! Roots uptake nutrients from the sediment column. 
        ! This uptake is not calculated here in the module BenthicEcology,
        ! but in the module SeagrassSedimInteraction, where the uptake of nutrients
        ! is stored in the array Me%Rates (See module SeagrassSedimInteraction),
        ! and then it is passed to the module BenthicEcology through 
        ! the InterfaceSedimentWater 
        ! Me%ExternalVar%UptakeNH4s = uptake of ammonia by roots
        ! Me%ExternalVar%UptakePO4s = Uptake of PO4 from sediment 
        
        ! gN/day
        UptakeNH4s= Me%ExternalVar%UptakeNH4s(index)
        
        ! gN/day
        UptakePO4s= Me%ExternalVar%UptakePO4s(index)
        
        ! Leaves growth depends on light in the water column
        ! the Light Limitation is calculated in
        ! the module SeagrassWaterInteraction and input here as Me%ExternalVar%LightFactor(index). 
        ! 
       
        LightFactor=max(0., Me%ExternalVar%LightFactor(index))
        
        ! gN/day
        !
        UptakeN   = UptakeNH4s +  UptakeNH4NO3w
        
        ! gP/day
        UptakeP   =  UptakePO4w + UptakePO4s
        

        
       
        !Temperature effect on seagrass
        s1 = (1. / (Me%Leaves%TOptMin - Me%Leaves%TMin)) * log((Me%Leaves%K2 * (1.0 - Me%Leaves%K1))                    &
                                                     / (Me%Leaves%K1 * (1.0 - Me%Leaves%K2)))

        s2 = (1. / (Me%Leaves%TMax - Me%Leaves%TOptMax)) * log((Me%Leaves%K3 * (1.0 - Me%Leaves%K4))                    &
                                                     / (Me%Leaves%K4 * (1.0 - Me%Leaves%K3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Me%Leaves%TMin))
        yb = exp(s2 * (Me%Leaves%TMax - Me%ExternalVar%Temperature(index)))

        xa = (Me%Leaves%K1 * ya) / (1.0 + Me%Leaves%K1 * (ya - 1.0))
        xb = (Me%Leaves%K4 * yb) / (1.0 + Me%Leaves%K4 * (yb - 1.0))

        TemperatureLim = xa * xb
        
        


        ! Total plant mass in KgDW
        TotalPlantMass = max(0.0000001, (Me%ExternalVar%Mass(L, Index)+Me%ExternalVar%Mass(R, Index)))

        
        if(Me%ComputeOptions%Nitrogen) then
        
        !  Me%ExternalVar%Mass(NINT, Index) is The Internal Nitrogen pool
        !  InternalNitrogen is the The Internal Nitrogen Content 
        !  expressed as gN/KgDW 
        InternalNitrogen=(Me%ExternalVar%Mass(NINT, Index)*1000.)/TotalPlantMass
        
        
        ! InternalNitrogenFunction is a dimensionless factor which limits the uptake of nutrients
        ! (ammonia and nitrate)
        ! Me%Nint%Nmax is the maximum internal nitrogen content in gN/KgDW
        ! Me%Nint%Nmin is the minimum internal nitrogen content in gN/KgDW
        InternalNitrogenFunction=max(0., min(1.,(Me%Nint%Nmax-InternalNitrogen)/(Me%Nint%Nmax-Me%Nint%Nmin)))
        
        
        ! NitrogenLimitation is a dimensionless factor which limits the growth of the plant. 
        ! Me%Nint%Nmin is the maximum internal nitrogen content in gN/KgDW
        ! Me%Nint%Ncrit is the critical internal nitrogen content in gN/KgDW
        NitrogenLimitation = max(0., min(1.,(InternalNitrogen-Me%Nint%Nmin)/(Me%Nint%Ncrit-Me%Nint%Nmin)))
        
        else
        InternalNitrogenFunction=1.
        NitrogenLimitation      = 1.
        
        endif


        if(Me%ComputeOptions%Phosphorus)then
        
        !  Me%ExternalVar%Mass(PINT, Index) is The Internal Phosphorus pool of the plant in KgP
        !  InternalPhosphorus is the the Internal Phosphorus Content 
        !  expressed as gP/KgDW 
        InternalPhosphorus=(Me%ExternalVar%Mass(PINT, Index)*1000.)/TotalPlantMass
        
        
        ! InternalPhosphorusFunction is a dimensionless factor which limits the uptake of phosphate
        ! Me%Pint%Pmax is the maximum Phosphorus quota in gP/KgDW 
        ! Me%Pint%Pmin is the minimum Phosphorus quota in gP/KgDW 
        InternalPhosphorusFunction= max(0., min(1.,(Me%Pint%Pmax-InternalPhosphorus)/(Me%Pint%Pmax-Me%Pint%Pmin)))
        
                        
        ! PhosphorusLimitation is a dimensionless factor which limits the growth of the plant. 
        ! Me%Pint%Pmin is the maximum Phosphorus quota in gP/KgDW 
        ! Me%Pint%Pcrit is the critical Phosphorus quota in gP/KgDW 
        PhosphorusLimitation = max(0., min(1.,(InternalPhosphorus-Me%Pint%Pmin)/(Me%Pint%Pcrit-Me%Pint%Pmin)))
        
        else
        InternalPhosphorusFunction=1.
        PhosphorusLimitation = 1.
        endif
        
        LeavesBioMass=max(0., (Me%ExternalVar%Mass(L, Index)/Me%ExternalVar%CellArea(Index)))
        

         ! If there are onlly seagrasses
        
        SpaceLimitation= 1.- max(0., min(1., (LeavesBioMass/Me%Leaves%KMAX)))
        
        

        NutrientsLimitation = min(NitrogenLimitation,PhosphorusLimitation)
        
        
        
        !Leaves Mortality 
        !KgDW * day-1
        LeavesMortality  = Me%ExternalVar%Mass(L, Index)      * &  ! KgDW   *
                           fm                                 * &  ! [-]   *
                           Me%Leaves%MortalityRate               ! 1/day 
           ! 1/day 
                           
       
                           
      ! Roots Mortality
      ! KgDW * day-1
       RootsMortality  =   Me%ExternalVar%Mass(R, Index)      * &  ! KgDW   *
                           fm                                 * &  ! [-]   *
                           Me%Roots%MortalityRate                ! 1/day 
             !day 
                          

       if ((InternalNitrogenFunction .EQ. 0.).OR.(InternalPhosphorusFunction .EQ. 0.)) then
       LeavesMortality=0.
       RootsMortality=0.
       endif

       ! gross growth
       ! KgDW * day-1
       Growth  =          Me%Leaves%MaxGrowthRate             * & ! 1/day *
                          NutrientsLimitation                 * & ! [-]   *
                          LightFactor                         * & ! [-]   *
                          TemperatureLim                      * & ! [-]   *
                          SpaceLimitation                     * & ! [-]   *
                          Me%ExternalVar%Mass(L, Index)           ! KgDW   
                          
       
                                 

      
        
        if(Me%ComputeOptions%Nitrogen)then
            
         !what passes from Leaves to PON (KgN)  =  KgDW/day      * gN/kgDW /1000             * day
            Me%Matrix(Index, L, PON)          = LeavesMortality * (Me%Leaves%gNKgDW/1000)*Me%DTDay

            Me%ExternalVar%MassInKgFromWater(PON,   Index) = Me%ExternalVar%MassInKgFromWater(PON,   Index) + &
                                                Me%Matrix(Index, L, PON)
          ! mortality of roots is not transfered to PON at the interface sediment-water, but at PON
          ! in the sediment column (see below)
            
        end if

        if(Me%ComputeOptions%Phosphorus)then

            !what passes from Leaves to POP
            !                                   kg DW            *gP/gDW              *day
           Me%Matrix(Index, L, POP)          = LeavesMortality * (Me%Leaves%gPKgDW/1000)*Me%DTDay

           Me%ExternalVar%MassInKgFromWater(POP,   Index) = Me%ExternalVar%MassInKgFromWater(POP,   Index) + &
                                                Me%Matrix(Index, L, POP)
                                                
          ! mortality of roots is not transfered to POP at the interface sediment-water, but at POP
          ! in the sediment column (see below)

        end if

          ! divide by 1000. to convert g to Kg
         Me%ExternalVar%Mass(NINT, Index)     = Me%ExternalVar%Mass(NINT, Index)          + &  ! KgN            
                                               (UptakeN/1000.                             - &  !(KgN/day         
                                                Growth*(Me%Leaves%gNKgDW/1000.) )         * &   ! kgdW/day * gN/KgdW/1000  
                                                 Me%DTDay                                 !- &     ! day 
                                               
          
          ! divide by 1000. to convert g to Kg
         Me%ExternalVar%Mass(PINT, Index)     = Me%ExternalVar%Mass(PINT, Index)          + &  ! KgP            
                                               (UptakeP/1000.                             - &  !(KgP/day         
                                                Growth*(Me%Leaves%gPKgDW/1000.))          * &  ! kgdW/day * gP/KgdW/1000  
                                                Me%DTDay                                  !- &   ! day 
                                               
                                                

         Me%ExternalVar%Mass(L, Index)  =     Me%ExternalVar%Mass(L, Index)               + &  ! kgdw     +
                                              ((1.-Me%Leaves%Ktr) *Growth                 - &  ! (Kgdw/day -
                                              LeavesMortality )                           * &  ! kgdw/day)* 
                                              Me%DTDay                                         ! day
         
         Me%ExternalVar%Mass(R, Index)  =     Me%ExternalVar%Mass(R, Index)               + &   ! kgdw     +
                                              (Me%Leaves%Ktr *Growth                      - &   ! (Kgdw/day -
                                              RootsMortality)                             * &   ! kgdw/day)* 
                                              Me%DTDay                                          ! day
         
         
        Me%ExternalVar%MassInKgFromWater(O2, Index)  = Me%ExternalVar%MassInKgFromWater(O2, Index)  +  &  ! KgO2         
                                                       Growth *  &                                        ! KgdW/day      
                                                       (Me%Leaves%gCKgDW/1000) *  &        ! gC/KgdW/1000. = KgC/KgO2 = gO2/gC 
                                                       (32.0 / 12.0)  *  &                 ! gO2/gC        = KgO2/KgC 
                                                       Me%DTDay                                        ! day
         
          case(Erosion)
          
              
       
             !Not all seagrasses mass is detached to make sure 
            !there's always enough to grow back again (minimum concentration)
            ! BEGIN ---This part of the algorithm is the same as in the macroalgae module--------
                   
                     if(Me%ExternalVar%Mass(L, Index) > Me%Leaves%MinimumBiomass)then
                    
                    
                     !what passes from seagrasses to POM
                     Me%Matrix(Index, L, PON)   = (Me%ExternalVar%Mass(L, Index) - &
                                                  Me%Leaves%MinimumBiomass)*(Me%Leaves%gNKgDW/1000.)
                                                  
                     Me%Matrix(Index, L, POP)   = (Me%ExternalVar%Mass(L, Index) - &
                                                  Me%Leaves%MinimumBiomass)*(Me%Leaves%gPKgDW/1000.)
                                               
                     else
                     !what passes from seagrasses to POM
                     Me%Matrix(Index, L, PON)   = 0.
                     Me%Matrix(Index, L, POP)   = 0.
                     
                     end if
                    
                    if(Me%ExternalVar%Mass(R, Index) > Me%Roots%MinimumBiomass) then 
                     
                    ! This is not a loss from the system:
                    ! Roots mortality will be stored in the array Me%StoredArray (see below)
                    ! and transferred to sediment properties
                     
                     RootsMortality   = Me%ExternalVar%Mass(R, Index) - &
                                                  Me%Roots%MinimumBiomass
                                                  

                    else
                     Me%Matrix(Index, R, PON)   = 0.
                     Me%Matrix(Index, R, POP)   = 0.
                    endif
                                
                    Me%ExternalVar%Mass(L, Index)  = Me%Leaves%MinimumBiomass
                    
                    Me%ExternalVar%Mass(R, Index)  = Me%Roots%MinimumBiomass
                         
                    
                    Me%ExternalVar%MassInKgFromWater(PON, Index) = Me%ExternalVar%MassInKgFromWater(PON, Index) + &
                                                      Me%Matrix(Index, L, PON)   
                                                      
                    Me%ExternalVar%MassInKgFromWater(POP, Index) = Me%ExternalVar%MassInKgFromWater(POP, Index) + &
                                                      Me%Matrix(Index, L, POP)   

          end select
          ! END  -----This part of the algorithm is the same as in the macroalgae module--------
         ! 
         ! Store quantities to be used in the other modules


          Me%StoredArray(Index, Me%StoredIndex%NintFactor)  = InternalNitrogenFunction       ! [-]
          Me%StoredArray(Index, Me%StoredIndex%NintFactorR) = InternalNitrogenFunction      ! [-]
          Me%StoredArray(Index, Me%StoredIndex%PintFactor)  = InternalPhosphorusFunction     ! [-]
          Me%StoredArray(Index, Me%StoredIndex%PintFactorR) = InternalPhosphorusFunction      ! [-]
          Me%StoredArray(Index, Me%StoredIndex%RootsMort)   = RootsMortality                  ! [KgDW/day]


    end subroutine ComputeSeagrasses
    
    
    !------------------------------------------------------------------------------
    
    subroutine ComputeBenthicPhyto(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: PON, POP, Phyto
        real                                        :: Mortality
        
        !Begin-----------------------------------------------------------------
        
        PON   = Me%PropIndex%PON
        POP   = Me%PropIndex%POP
        Phyto = Me%PropIndex%Phyto

        !kg * day * day-1
        Mortality = Me%ExternalVar%Mass(Phyto, Index) * Me%DTDay * Me%Phyto%MortalityRate

        Me%ExternalVar%Mass(Phyto, Index) = Me%ExternalVar%Mass(Phyto, Index) - Mortality

        if(Me%ComputeOptions%Nitrogen)then
            
            !what passes from Phyto to PON
            Me%Matrix(Index, Phyto, PON)      = Mortality * Me%Phyto%NC_Ratio

            Me%ExternalVar%Mass(PON,   Index) = Me%ExternalVar%Mass(PON,   Index) + &
                                                Me%Matrix(Index, Phyto, PON)

        end if

        if(Me%ComputeOptions%Phosphorus)then

            !what passes from Phyto to POP
            Me%Matrix(Index, Phyto, POP)      = Mortality * Me%Phyto%PC_Ratio

            Me%ExternalVar%Mass(POP,   Index) = Me%ExternalVar%Mass(POP,   Index) + &
                                                Me%Matrix(Index, Phyto, POP)

        end if


    end subroutine ComputeBenthicPhyto
    
    !--------------------------------------------------------------------------
    
    
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



        subroutine KillBenthicEcology(ObjBenthicEcologyID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjBenthicEcologyID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBenthicEcologyID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mBenthicEcology_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%PropertyList)

                !Deallocates Instance
                call DeallocateInstance ()


                ObjBenthicEcologyID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine KillBenthicEcology
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_BenthicEcology), pointer          :: AuxObjBenthicEcology
        type (T_BenthicEcology), pointer          :: PreviousObjBenthicEcology

        !Updates pointers
        if (Me%InstanceID == FirstObjBenthicEcology%InstanceID) then
            FirstObjBenthicEcology => FirstObjBenthicEcology%Next
        else
            PreviousObjBenthicEcology => FirstObjBenthicEcology
            AuxObjBenthicEcology      => FirstObjBenthicEcology%Next
            do while (AuxObjBenthicEcology%InstanceID /= Me%InstanceID)
                PreviousObjBenthicEcology => AuxObjBenthicEcology
                AuxObjBenthicEcology      => AuxObjBenthicEcology%Next
            enddo

            !Now update linked list
            PreviousObjBenthicEcology%Next => AuxObjBenthicEcology%Next

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

    subroutine Ready (ObjBenthicEcology_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthicEcology_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjBenthicEcology_ID > 0) then
            call LocateObjBenthicEcology (ObjBenthicEcology_ID)
            ready_ = VerifyReadLock (mBenthicEcology_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjBenthicEcology (ObjBenthicEcologyID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthicEcologyID

        !Local-----------------------------------------------------------------

        Me => FirstObjBenthicEcology
        do while (associated (Me))
            if (Me%InstanceID == ObjBenthicEcologyID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleBenthicEcology - LocateObjBenthicEcology - ERR01'

    end subroutine LocateObjBenthicEcology
    
end module ModuleBenthicEcology

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------


